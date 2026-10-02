#!/usr/bin/env python3
"""Controller-authorized, serial RH parameter diagnosis; never evolves binary data."""
import sys
sys.dont_write_bytecode = True
import csv, importlib.util, json, math, subprocess
from fractions import Fraction
from pathlib import Path
import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('work', HERE/'t17-work.py')
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)
ROOT = w.ROOT/'finder-diagnosis'
KINDS = ('unboosted', 'boosted', 'binary')
METHODS = ((32, 1), (128, 1), (512, 1), (128, 16))
AUTHOR = tuple(sorted((w.REPO/'Source/RHFinder').glob('*.hpp')))

def rows(path):
    return list(csv.DictReader(path.open()))

def checked(d, allowed=(0,)):
    r = json.loads((d/(d.name+'.resources.json')).read_text())
    assert r['returncode'] in allowed and r['child_measurement']['returncode'] == r['returncode']
    assert not r['gate_reason'] and r['peak_rss_bytes'] < 6e9 and r['tree_peak_bytes'] < 6.5e9
    assert int((d/'done.exit').read_text()) == r['returncode'] % 256
    return r

def pins():
    for item in json.loads((ROOT/'registration.json').read_text())['pins']:
        assert w.digest(Path(item['path'])) == item['sha256'], item['path']

def register():
    ROOT.mkdir(exist_ok=True)
    pinfiles = [w.PROFILE, w.CTT, w.ROOT/'rh.ex', w.ROOT/'production-t17.ex',
                w.ROOT/'dense-production.ex', *AUTHOR, HERE/'t17-rh.cpp', Path(__file__)]
    registration = dict(methods=METHODS, baseline=dict(chase=1, quota=1, cap=2000, floor=64, seconds=1795),
        cap_only=dict(chase=1, quota=1, cap=100000, floor=100001, seconds=235),
        accelerated=dict(cap=100000, floor=100001, seconds=595, confirmation_seconds=1795),
        angular_points=[48, 96], thresholds=[1e-7, 1e-10, 1e-12],
        compact_control_grid=dict(L=112, N1=64, N2=32, h0=1.75, center=40, translation=2184,
                                  levels=[12, 13], static_t0_only=True),
        seed_sensitivity=dict(radius_factors=[.9, 1.1], center_offsets_Rh=[-.1, .1], seconds=595),
        original_verdict_preserved=True,
        pins=[dict(path=str(f), sha256=w.digest(f), bytes=f.stat().st_size) for f in pinfiles])
    target=ROOT/'registration.json'
    if target.exists(): assert json.loads(target.read_text()) == registration
    else: target.write_text(json.dumps(registration, indent=2)+'\n')
    d=ROOT/'pipeline'; d.mkdir(exist_ok=True)
    (d/'plan.json').write_text(json.dumps(dict(jobs=[dict(name='finder-diagnosis', directory=str(d/'job'),
        command=[w.PY, str(Path(__file__)), 'worker'])]), indent=2)+'\n')
    print(d/'plan.json', flush=True)

def local_faces(cp, center):
    with h5py.File(cp) as f:
        level=int(f.attrs['num_levels'])-1; g=f[f'level_{level}']; h=float(g.attrs['dx'])
        axis=sorted((int(b['lo_i'])*h,(int(b['hi_i'])+1)*h) for b in g['boxes'][:] if int(b['lo_j'])==0)
        merged=[]
        for lo,hi in axis:
            if merged and lo<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],hi)
            else: merged.append([lo,hi])
        lo,hi=next(b for b in merged if b[0]<=center<b[1])
        top=max((int(b['hi_j'])+1)*h for b in g['boxes'][:] if int(b['lo_i'])*h<=center<(int(b['hi_i'])+1)*h)
        return dict(level=level,h_M=h,left_face_M=center-lo,right_face_M=hi-center,y_face_M=top,
                    center_cell_phase=float((Fraction(center)/Fraction(h)-Fraction(1,2))%1))

def initialize(rung, kind):
    level=12 if rung=='coarse' else 13
    text=w.p.replace(w.study_base(level), dict(N1=64,N2=32,L=112,center='56 0',star_centre='40 0',
        mass_extraction_center='40 0',binary='false',ems_ctt_data_path='""',separation=0,
        boosted=str(kind=='boosted').lower(),boost_rapidity=.05238600881798246 if kind=='boosted' else 0,
        ems_binary_refinement='false',ems_track_punctures='false',checkpoint_interval=0,
        plot_interval=-1, max_steps=0,stop_time=0,t14_dense_initial_tags='false'))
    d=ROOT/'initialization'/f'{rung}-{kind}-native'
    rc=w.run_case(d.name,d,w.ROOT/'production-t17.ex',text)
    r=checked(d,(0,-11))
    if rc:
        assert r['returncode']==-11
        d=ROOT/'initialization'/f'{rung}-{kind}-dense'
        text=w.p.replace(text,dict(t14_dense_initial_tags='true'))
        assert w.run_case(d.name,d,w.ROOT/'dense-production.ex',text)==0
        checked(d)
    assert 'GRChombo finished.' in (d/'run.log').read_text()
    cp=d/'chk/EMS_000000.2d.hdf5'; assert cp.exists()
    observed=local_faces(cp,40.)
    binary=Path(json.loads((w.ROOT/'choices.json').read_text())[rung]['geometric_checkpoint'])
    reference=local_faces(binary,2224.)
    # At fixed spacing the dyadic translation must preserve all horizon stencil faces and phase.
    assert observed==reference, (observed,reference)
    return dict(text=text,cp=cp,centres=[40.],faces=observed)

def state(rung, kind, step=0):
    if kind!='binary': return initialize(rung,kind)
    timing=w.ROOT/'timing'/f'{rung}-binary'
    text=(timing/'params.txt').read_text()
    if step:
        cp=timing/'chk/EMS_000001.2d.hdf5'
        records=[line.split() for line in (timing/'punctures.dat').read_text().splitlines()
                 if line.strip() and not line.lstrip().startswith('#')]
        rec=next(x for x in records if abs(float(x[0])-.4375)<1e-8)
        centres=[float(rec[1]),float(rec[3])]
    else:
        cp=Path(json.loads((w.ROOT/'choices.json').read_text())[rung]['geometric_checkpoint'])
        centres=[2224.,2256.]
    assert cp.exists()
    return dict(text=text,cp=cp,centres=centres,faces=local_faces(cp,centres[0]))

def probe(rung, kind, n, method, step=0, radius=1., offset=0., seconds=595):
    speed,quota,cap,floor=method
    s=state(rung,kind,step); count=len(s['centres'])
    name=f'{rung}-{kind}-t{step}-n{n}-s{speed}-q{quota}-c{cap}-f{floor}-r{radius:g}-o{offset:g}-wall{seconds}'
    d=ROOT/'probes'/name
    centres=[c+(1 if k==0 else -1)*offset*w.RH for k,c in enumerate(s['centres'])]
    params=w.p.replace(s['text'],dict(restart_file=s['cp'], RH_activate='false',RH_num_horizons=count,
        RH_initial_centre=' '.join(map(str,centres)),RH_initial_radii=' '.join([format(w.RH*radius,'.17g')]*count),
        RH_num_points=' '.join([str(n)]*count),RH_level=' '.join(['0']*count),
        RH_start_times=' '.join(['0']*count),RH_time_step_freq=' '.join([str(quota)]*count),
        RH_chase_speeds=' '.join([str(speed)]*count),RH_newton_crit=' '.join(['0']*count),
        offline_points=n,offline_seconds=seconds,offline_max_updates=cap,offline_floor_window=floor,
        checkpoint_interval=-1,plot_interval=-1,ems_track_punctures='false',t13_launch_stop_time=0))
    rc=w.run_case(name,d,w.ROOT/'rh.ex',params)
    r=checked(d,(0,1,2))
    if rc==2: assert 'offline surface died' in (d/'run.log').read_text(), d
    raw=rows(d/'surfaces.csv'); finals=[x for x in raw if int(x['stage'])==2]
    passed=(rc==0 and len(finals)==count and all(x['status']=='FOUND' and
            math.isfinite(float(x['expansion_squared'])) and float(x['expansion_squared'])<=1e-12 for x in finals))
    assert (rc==0)==passed, d
    data=[]
    for x in raw:
        shape=np.loadtxt(d/f"shape-{x['search_index']}-{x['stage']}.dat")[2:]
        data.append(dict(rung=rung,kind=kind,step=step,speed=speed,quota=quota,cap=cap,floor=floor,
            radius_factor=radius,center_offset_Rh=offset,wall_limit=seconds,strict_PASS=passed,
            peak_RSS_bytes=r['peak_rss_bytes'],child_returncode=r['returncode'],directory=str(d),
            **x,expansion_RMS=math.sqrt(float(x['expansion_squared'])),radius_min=float(shape.min()),
            radius_max=float(shape.max()),shape_span_fraction=float(np.ptp(shape)/np.mean(shape)),
            centre_offset_from_puncture_M=float(x['centre'])-s['centres'][int(x['search_index'])]))
    return dict(passed=passed, directory=d, raw=raw, data=data,
        worst=max((float(x['expansion_squared']) for x in raw[-count:]),default=math.inf) if rc!=2 else math.inf)

def publish():
    data=[]
    for d in sorted((ROOT/'probes').glob('*')):
        if not (d/'done.exit').exists(): continue
        r=checked(d,(0,1,2))
        p=w.p.parameters((d/'params.txt').read_text()); kind=d.name.split('-')[1]; rung=d.name.split('-')[0]
        step=int(d.name.split('-')[2][1:]); count=int(p['RH_num_horizons'])
        raw=rows(d/'surfaces.csv'); last=[x for x in raw if int(x['stage'])==2]
        passed=r['returncode']==0 and len(last)==count and all(x['status']=='FOUND' and float(x['expansion_squared'])<=1e-12 for x in last)
        for x in raw:
            shape=np.loadtxt(d/f"shape-{x['search_index']}-{x['stage']}.dat")[2:]
            data.append(dict(rung=rung,kind=kind,step=step,speed=p['RH_chase_speeds'].split()[0],quota=p['RH_time_step_freq'].split()[0],
                max_updates=p['offline_max_updates'],floor_window=p['offline_floor_window'],wall_limit=p['offline_seconds'],
                strict_PASS=passed,peak_RSS_bytes=r['peak_rss_bytes'],child_returncode=r['returncode'],directory=str(d),
                **x, expansion_RMS=math.sqrt(float(x['expansion_squared'])),radius_min=float(shape.min()),
                radius_max=float(shape.max()),shape_span_fraction=float(np.ptp(shape)/np.mean(shape))))
    if data: w.save('t17-finder-diagnosis.csv',data)

def admission(confirmation):
    angular=[]; spatial=[]; found={}; verdict={}
    for rung in ('coarse','fine'):
        ok=True
        for step in (0,1):
            pair=[next(x for x in confirmation if x['rung']==rung and x['kind']=='binary' and
                       x['step']==step and x['N_theta']==n) for n in (48,96)]
            ok &= all(x['PASS'] for x in pair)
            if not all(x['PASS'] for x in pair): continue
            surfaces=[{int(x['search_index']):x for x in rows(Path(c['directory'])/'surfaces.csv')
                       if int(x['stage'])==2} for c in pair]
            for hole in (0,1):
                A,Q=float(surfaces[1][hole]['A']),float(surfaces[1][hole]['Q'])
                da=abs(float(surfaces[0][hole]['A'])-A)/abs(A)
                dq=abs(float(surfaces[0][hole]['Q'])-Q)/abs(Q)
                angular.append(dict(rung=rung,step=step,hole=hole,relative_A=da,relative_Q=dq,PASS=da<=1e-3 and dq<=1e-3))
                ok &= da<=1e-3 and dq<=1e-3
                found[rung,step,hole]=dict(A=A,Q=Q,e=Q/math.sqrt(A/(4*math.pi)))
        verdict[rung]=bool(ok)
    spatial_ok=True
    for step in (0,1):
        for hole in (0,1):
            if ('coarse',step,hole) not in found or ('fine',step,hole) not in found:
                spatial_ok=False; continue
            c,f=found['coarse',step,hole],found['fine',step,hole]
            da,dq=abs(c['A']-f['A'])/abs(f['A']),abs(c['Q']-f['Q'])/abs(f['Q'])
            de=abs(c['e']-f['e']); ok=da<=.01 and dq<=.01 and de<=.1
            spatial_ok &= ok
            spatial.append(dict(step=step,hole=hole,relative_A=da,relative_Q=dq,delta_e=de,PASS=ok))
    if angular: w.save('t17-finder-diagnosis-angular.csv',angular)
    if spatial: w.save('t17-finder-diagnosis-spatial.csv',spatial)
    endpoint=rows(HERE/'t17-launch-endpoint.csv')
    for rung in verdict:
        required=[x for x in endpoint if x['rung']==rung and x['field']=='Gamma']
        verdict[rung] &= len(required)==8 and all(x['qualified']=='True' and float(x['upper_5x'])<1 for x in required)
    recommendation='coarse' if verdict['coarse'] and verdict['fine'] and spatial_ok else 'fine' if verdict['fine'] else None
    return dict(binary_horizon_and_launch=verdict,coarse_spatial_PASS=spatial_ok,recommendation=recommendation)

def self_check():
    for level in (12,13):
        h=1.75/2**level
        assert (Fraction(2224)/Fraction(h)-Fraction(40)/Fraction(h)).denominator==1
        assert (Fraction(2224)/Fraction(h)-Fraction(1,2))%1 == (Fraction(40)/Fraction(h)-Fraction(1,2))%1
    pins()
    print('PASS pinned inputs/author source and exact control-grid translation/phase.',flush=True)

def worker():
    pins()
    print('Registered diagnosis: original T17 verdict stays unchanged.',flush=True)
    face=[]
    for rung in ('coarse','fine'):
        for kind in ('unboosted','boosted'):
            s=state(rung,kind);face.append(dict(rung=rung,kind=kind,**s['faces'],checkpoint=str(s['cp']),matched_binary_left=True))
            print('Matched grid',rung,kind,s['faces'],flush=True)
    w.save('t17-finder-control-grids.csv',face)
    for rung in ('coarse','fine'):
        for kind in ('unboosted','boosted'):
            for n in (48,96):
                probe(rung,kind,n,(1,1,2000,64),seconds=1795);publish()
    for kind in KINDS:
        for n in (48,96):
            probe('coarse',kind,n,(1,1,100000,100001),seconds=235);publish()
    best=None
    for speed,quota in METHODS:
        results=[]
        for kind in KINDS:
            for n in (48,96):
                result=probe('coarse',kind,n,(speed,quota,100000,100001));publish();results.append(result)
        score=max(x['worst'] for x in results)
        if best is None or score<best[0]: best=(score,speed,quota,results)
        if all(x['passed'] for x in results): break
    _,speed,quota,results=best
    method=(speed,quota,100000,100001)
    # Fixed-budget seed tests distinguish a basin problem from slow/limited chasing.
    if not all(x['passed'] for x in results):
        for radius,offset in ((.9,0.),(1.1,0.),(1.,-.1),(1.,.1)):
            probe('coarse','binary',96,method,radius=radius,offset=offset);publish()
    confirmation=[]
    for rung in ('coarse','fine'):
        for kind in KINDS:
            for n in (48,96):
                # Reuse a scientifically passing 595-s probe, with identical settings.
                if rung=='coarse':
                    old=next(x for x in results if x['data'][0]['kind']==kind and int(x['raw'][0]['N_theta'])==n)
                    result=old if old['passed'] else probe(rung,kind,n,method,seconds=1795)
                else: result=probe(rung,kind,n,method,seconds=1795)
                confirmation.append(dict(rung=rung,kind=kind,N_theta=n,step=0,PASS=result['passed'],directory=str(result['directory'])))
                publish()
        for n in (48,96):
            result=probe(rung,'binary',n,method,step=1,seconds=1795)
            confirmation.append(dict(rung=rung,kind='binary',N_theta=n,step=1,PASS=result['passed'],directory=str(result['directory'])))
            publish()
    pins()
    w.save('t17-finder-confirmation.csv',confirmation)
    decision=admission(confirmation)
    if decision['recommendation']=='coarse':
        # A newly qualifying cheaper rung needs its own full-domain single-hole rate control.
        w.timed_single(json.loads((w.ROOT/'choices.json').read_text()),'coarse')
        checked(w.ROOT/'timing/coarse-single')
    (HERE/'t17-finder-diagnosis-status.json').write_text(json.dumps(dict(complete=True,speed=speed,quota=quota,
        cap=100000,floor_window=100001,thresholds=[1e-7,1e-10,1e-12],confirmation=confirmation,
        original_verdict_preserved=True,admission=decision),indent=2)+'\n')
    print('Diagnosis receipts complete; exit zero is not horizon admission.',flush=True)

if __name__=='__main__':
    {'register':register,'worker':worker,'publish':publish,'self-check':self_check}[sys.argv[1]]()
