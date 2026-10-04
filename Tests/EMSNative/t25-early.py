#!/usr/bin/env python3
"""Bounded T25 early-checkpoint maps and one bracket-seeded search per hole.

No evolution; one serial checkpoint process at a time. Native maps run at
N96 and N192. Geometric radial spacing is 4%; measured sign intervals are
refined by midpoint bisection (up to four passes). If a sphere has no
pointwise pair, the fixed fallback is offsets {-4,-2,0,2,4} finest cells and
c/a in {0.5,0.75,0.9,1,1.1,1.25,1.5,2}, never using an area target.
One best negative-to-positive pointwise bracket is selected; if absent,
one negative-to-positive mean bracket may seed a search, explicitly marked
weaker. All finder output remains trial-surface data unless all stages pass.
CLI: t25-early.py initialize | maps | finders | summarize | all
"""
import argparse,csv,hashlib,importlib.util,json,math,os,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t25/early')
INPUT=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/early-chk')
TRACK=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall/ev1/runs/exp-0024/merger/punctures.dat')
EXE=Path('/private/tmp/ems-t25/build/t25-expansion-map.ex')
H=1.75/4096
PRODUCTION=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall/ev1/runs/exp-0024/merger/runtime/params-production.txt')
def load_module(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
w=load_module('t25-work');b=load_module('t25-brackets')
def read(p):return list(csv.DictReader(Path(p).open()))
def write(p,rows):b.write(p,rows)
def atomic_json(p,obj):
    p=Path(p);p.with_suffix('.tmp').write_text(json.dumps(obj,indent=2)+'\n');p.with_suffix('.tmp').replace(p)
def verified(d,label,expected,complete,cap=3000000000):
    rc=int((d/'done.exit').read_text());r=json.loads((d/(label+'.resources.json')).read_text())
    assert rc==r['returncode']==r['child_measurement']['returncode']==expected,(d,r)
    assert not r['gate_reason'] and r['peak_rss_bytes']<cap and r['cap_GB']==cap/1e9
    assert complete in (d/'run.log').read_text(),d
    return r
def samples():
    track={float(x[0]):[float(x[1]),float(x[3])] for line in TRACK.read_text().splitlines()
           if (x:=line.split()) and not line.startswith('#')}
    out=[dict(tag='step000000',time=0.,centres=track[0.],
        cp=ROOT/'initialization-d16/chk/EMS_000000.2d.hdf5',
        base=PRODUCTION,configuration='exp0024_d16_regenerated_t0')]
    for step in range(4,33,4):
        t=step*.875;out.append(dict(tag=f'step{step:06}',time=t,centres=track[t],
            cp=INPUT/f'EMS_{step:06}.2d.hdf5',
            base=PRODUCTION,configuration='exp0024_d16'))
    return out
def initialize():
    d=ROOT/'initialization-d16';d.mkdir(parents=True,exist_ok=True)
    trumpet=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo-binary/a0.9-e8.trumpet')
    companion=Path('/Users/auroradysis/Workspace/EMS/.data/echo-binary/tranche3/d16/r32-a10/n32-r6.ctt')
    inputs=[]
    for p,expected in [(trumpet,'4d7cc0b975d4dde0c4da8400ef0bbad8c3e57716c9708c160d19e56288a76f90'),
        (companion,'107370ae00d8b69dc3122dc023e34bbff23b90afc8401bb233743c1182e6083c')]:
        with p.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
        assert actual==expected
        inputs.append(dict(path=str(p),sha256=actual,bytes=p.stat().st_size,role='t0 initialization only'))
    (d/'chk').mkdir(exist_ok=True);(d/'plt').mkdir(exist_ok=True)
    changes=dict(ems_data_path=trumpet,ems_ctt_data_path=companion,max_steps=0,stop_time=0,
        data_subpath='.',data_path='./',output_path='./',hdf5_subpath='""',
        chk_prefix='chk/EMS_',plot_prefix='plt/EMS_Plot_',plot_interval=-1,plot_period=-1,checkpoint_interval=0,
        RH_activate='false',RH_num_horizons=0,ems_track_punctures='false',t2_guard_initial_data_after_t0='true')
    (d/'params.txt').write_text(w.replace(PRODUCTION.read_text(),changes))
    if not (d/'done.exit').exists():assert w.measure('initialization',d,['/private/tmp/ems-t17/production-t17.ex','params.txt'],6000000000)==0
    r=json.loads((d/'initialization.resources.json').read_text())
    assert int((d/'done.exit').read_text())==r['returncode']==r['child_measurement']['returncode']==0
    assert not r['gate_reason'] and r['peak_rss_bytes']<6e9 and 'GRChombo finished.' in (d/'run.log').read_text()
    cp=d/'chk/EMS_000000.2d.hdf5'
    with cp.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
    atomic_json(d/'inputs.json',dict(inputs=inputs,checkpoint_sha256=sha,peak_rss_bytes=r['peak_rss_bytes'],
        initializer='/private/tmp/ems-t17/production-t17.ex',initializer_sha256=hashlib.sha256(Path('/private/tmp/ems-t17/production-t17.ex').read_bytes()).hexdigest(),
        t0_only=True,no_advances=True))
    print('D16_T0_INITIALIZATION_COMPLETE',r['peak_rss_bytes'],flush=True)
def verify_inputs():
    wanted={x[1].lstrip('*'):x[0] for line in (INPUT/'remote.sha256').read_text().splitlines() if (x:=line.split())}
    out=[]
    for s in samples():
        p=s['cp'];expected=wanted[p.name] if s['time'] else json.loads((ROOT/'initialization-d16/inputs.json').read_text())['checkpoint_sha256']
        with p.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
        assert actual==expected,(p,actual,expected)
        out.append(dict(path=str(p),sha256=actual,bytes=p.stat().st_size,time=s['time'],configuration=s['configuration']))
        print('INPUT_VERIFIED',s['tag'],flush=True)
    atomic_json(ROOT/'inputs.json',out)
def radii(spacing=1.04):
    rs=[3*H];r=3*H
    while r*spacing<.05:r*=spacing;rs.append(r)
    return rs+[.05]
def specs(s,fallback=False):
    rows=[]
    for hole,centre in enumerate(s['centres']):
        offsets=[-4,-2,0,2,4] if fallback else [0]
        aspects=[.5,.75,.9,1.,1.1,1.25,1.5,2.] if fallback else [1.]
        for off in offsets:
            for aspect in aspects:
                family=f'h{hole}-o{off:+d}-q{aspect:g}'
                for j,r in enumerate(radii(1.08 if fallback else 1.04)):
                    # Minimum extent is at least three finest cells.
                    if min(r,r*aspect)<3*H:continue
                    rows.append(dict(id=f'{family}-r{j:03}',family=family,centre=centre+off*H,a=r,c=r*aspect))
    return rows
def mapping(s,part,rows,n):
    d=ROOT/s['tag']/part/f'N{n}';d.mkdir(parents=True,exist_ok=True)
    w.surfaces(d/'surfaces.csv',rows)
    changes=dict(restart_file=s['cp'],hdf5_subpath='""',ems_data_path='/T25-FORBIDDEN-static',ems_ctt_data_path='""',
        max_steps=0,stop_time=0,RH_activate='false',RH_num_horizons=0,ems_track_punctures='false',
        map_points=n,map_surfaces=d/'surfaces.csv',map_output='map.csv',map_angles='angles.csv',map_find='false',
        map_punctures=' '.join(map(str,s['centres'])))
    text=w.replace(s['base'].read_text(),changes);(d/'params.txt').write_text(text)
    if not (d/'done.exit').exists():
        assert w.measure('early-map',d,[EXE,'params.txt'],3000000000)==0,d
    verified(d,'early-map',0,'T25_EXPANSION_MAP_COMPLETE; no advances')
    rs=read(d/'map.csv');assert len(rs)==len(rows) and all(float(r['time'])==s['time'] for r in rs)
    assert all(r['id']==x['id'] and all(float(r[k])==float(x[k]) for k in ('centre','a','c')) for r,x in zip(rs,rows)),d
    return rs
def strong(rs):return b.pointwise_barriers(rs)
def refine(s,rs,prefix='refine'):
    # Bisect edges of both uniform-sign sets and each average crossing.
    for iteration in range(4):
        added=[];groups={}
        for r in rs:groups.setdefault(r['family'],[]).append(r)
        for family,xs in groups.items():
            xs.sort(key=lambda x:float(x['a']))
            for x,y in zip(xs,xs[1:]):
                if x['status']!='RESOLVED' or y['status']!='RESOLVED':continue
                def signs(z):return (float(z['theta_max'])<0,float(z['theta_min'])>0,float(z['theta_mean'])<0)
                if signs(x)==signs(y):continue
                if float(y['a'])/float(x['a'])<=1.005:continue
                r=(float(x['a'])+float(y['a']))/2
                added.append(dict(id=f'{family}-{prefix}{iteration}-{len(added):04}',family=family,
                    centre=x['centre'],a=r,c=r*float(x['c'])/float(x['a'])))
        if not added:break
        rs+=mapping(s,f'{prefix}-{iteration}',added,192)
    return rs
def select(rs,hole):
    xs=[r for r in rs if r['family'].startswith(f'h{hole}-')]
    candidates=strong(xs)
    if not candidates:
        candidates=[x for x in b.brackets(xs) if x['status']=='AVERAGE_ONLY_WEAKER' and
                    x['theta_inner_mean']<0 and x['theta_outer_mean']>0]
    if not candidates:return None
    # Angular RMS of the two actual endpoint surfaces, then width; no area target.
    byid={r['id']:r for r in xs}
    return min(candidates,key=lambda x:(max(float(byid[x[k]]['theta_rms']) for k in ('inner','outer')),
                                       float(x['a_outer'])/float(x['a_inner'])))
def maps():
    ROOT.mkdir(parents=True,exist_ok=True);verify_inputs()
    for s in samples():
        d=ROOT/s['tag'];print('MAP_BEGIN',s['tag'],flush=True)
        if (d/'selection.json').exists():
            print('MAP_REUSE',s['tag'],flush=True);continue
        initial=specs(s);r96=mapping(s,'spheres',initial,96);rs=mapping(s,'spheres',initial,192)
        if any(not strong([x for x in rs if x['family'].startswith(f'h{k}-')]) for k in range(2)):
            extra=[x for x in specs(s,True) if x['family'] not in {r['family'] for r in rs}]
            # Thin two-resolution fallback; all angular rows retained externally.
            r96+=mapping(s,'fallback',extra,96);rs+=mapping(s,'fallback',extra,192)
        rs=refine(s,rs)
        selections=[]
        for hole in range(2):
            chosen=select(rs,hole)
            if chosen:
                sd=d/f'hole{hole}-root';sd.mkdir(parents=True,exist_ok=True)
                write(sd/'bracket.csv',[chosen])
                root=[dict(id=f'h{hole}-root',family=f'h{hole}-root',centre=chosen['centre'],
                    a=chosen['mean_linear_a_root'],c=chosen['mean_linear_c_root'])]
                rr=mapping(s,f'hole{hole}-root',root,192)[0]
                mapping(s,f'hole{hole}-root',root,96)
                endpoints=[dict(id=chosen[k],family=chosen['family'],centre=chosen['centre'],
                    a=chosen['a_'+k],c=chosen['c_'+k]) for k in ('inner','outer')]
                checked=mapping(s,f'hole{hole}-endpoints',endpoints,96)
                if chosen['status']=='POINTWISE_NEGATIVE_TO_POSITIVE':
                    assert float(checked[0]['theta_max'])<0 and float(checked[1]['theta_min'])>0,('angular sign change',s,chosen)
                else:assert float(checked[0]['theta_mean'])<0 and float(checked[1]['theta_mean'])>0
                # Permit a weak trial seed only under the explicit Tests-only switch.
                (sd/'finder-params.txt').write_text(w.replace((sd/'N192/params.txt').read_text(),dict(
                    map_surfaces=sd/'N192/surfaces.csv',map_points=48,map_find='true',map_finder_bracket=sd/'bracket.csv',
                    map_find_allow_average='true' if chosen['status']=='AVERAGE_ONLY_WEAKER' else 'false',
                    map_find_seconds=235,map_find_max_updates=100000,map_find_floor_window=100001,
                    map_find_chase=1,map_find_thresholds='1e-7 1e-10 1e-12')))
                selections.append(dict(hole=hole,bracket=chosen,root=rr,seed_dir=str(sd)))
            else:selections.append(dict(hole=hole,bracket=None,root=None))
        if s['time']==28:
            mid=sum(s['centres'])/2;half=abs(s['centres'][1]-s['centres'][0])/2
            common=[]
            for q in [1.,2.,4.,8.]:
                for j,c in enumerate([1.05*half,1.25*half,1.5*half,2*half,4*half]):
                    common.append(dict(id=f'mid-q{q:g}-{j}',family=f'mid-q{q:g}',centre=mid,a=c/q,c=c))
            midpoint=mapping(s,'midpoint',common,192)
            write(HERE/'t25-early-midpoint-t28.csv',midpoint)
            write(HERE/'t25-early-midpoint-t28-brackets.csv',b.brackets(midpoint))
        # Packet CSV stays small; large angular dumps remain outside the worktree.
        selected_families={x['bracket']['family'] for x in selections if x['bracket']}
        packet=[r for r in rs if r['family'] in selected_families or r['family'].endswith('-o+0-q1')]
        write(HERE/f't25-early-{s["tag"]}.csv',packet)
        write(HERE/f't25-early-{s["tag"]}-brackets.csv',b.brackets(rs))
        write(HERE/f't25-early-{s["tag"]}-barriers.csv',strong(rs))
        write(d/'map-N96.csv',r96);write(d/'map-N192.csv',rs)
        atomic_json(d/'selection.json',dict(time=s['time'],configuration=s['configuration'],selections=selections))
        print('MAP_DONE',s['tag'],[(x['hole'],x['bracket']['status'] if x['bracket'] else 'NO_BRACKET') for x in selections],flush=True)
    summarize()
    atomic_json(ROOT/'maps-complete.json',dict(status='PASS',no_advances=True))
def finders():
    improve_wide_barriers()
    for s in samples():
        data=json.loads((ROOT/s['tag']/'selection.json').read_text())
        for item in data['selections']:
            if not item['bracket']:continue
            sd=Path(item['seed_dir']);d=sd/'finder';d.mkdir(exist_ok=True)
            (d/'params.txt').write_text((sd/'finder-params.txt').read_text())
            print('FINDER_BEGIN',s['tag'],item['hole'],item['bracket']['status'],flush=True)
            if not (d/'done.exit').exists():
                rc=w.measure('early-finder',d,[EXE,'params.txt'],6000000000)
                assert rc in (0,1),(d,rc)
            rc=int((d/'done.exit').read_text());verified(d,'early-finder',rc,'T25_FINDER_COMPLETE',6000000000)
            stages=read(d/'finder.csv');assert stages and all(float(r['time'])==s['time'] for r in stages)
            assert all(math.isfinite(float(r[k])) for r in stages for k in
                ('expansion_squared','A','Q','seconds','centre','r_min','r_max')),d
            if rc==0:assert len(stages)==3 and all(r['status']=='FOUND' and float(r['expansion_squared'])<=float(r['threshold']) for r in stages)
            else:assert stages[-1]['status'] in ('TIME_CAP','FLOOR','UPDATE_CAP')
            write(HERE/f't25-early-{s["tag"]}-h{item["hole"]}-finder.csv',stages)
            print('FINDER_DONE',s['tag'],item['hole'],stages[-1]['status'],stages[-1]['expansion_squared'],flush=True)
            summarize()
    atomic_json(ROOT/'finders-complete.json',dict(status='COMPLETE',every_probe_receipt_verified=True))
def improve_wide_barriers():
    """Attempt the requested few-percent barrier before any physical search.

    A radial bisection cannot shrink an angular root band. If the sphere
    pair remains wider than 3%, add the same fixed shifted/spheroid family
    used when no sphere pair exists. Existing maps are never repeated.
    If a broad pair remains, report that fact instead of fabricating a
    narrow pointwise bracket. The trial root uses the narrow adjacent
    mean crossing inside the enclosing pointwise barriers where available.
    """
    for s in samples():
        d=ROOT/s['tag']
        if (d/'improved.json').exists():continue
        data=json.loads((d/'selection.json').read_text())
        # This runs before finders only; do not invalidate a completed probe.
        assert not any((Path(x['seed_dir'])/'finder/done.exit').exists() for x in data['selections'] if x['bracket'])
        rs=read(d/'map-N192.csv');r96=read(d/'map-N96.csv')
        needs=any(not x['bracket'] or x['bracket']['status']!='POINTWISE_NEGATIVE_TO_POSITIVE' or
            float(x['bracket']['a_outer'])/float(x['bracket']['a_inner'])>1.03 for x in data['selections'])
        if needs:
            existing={r['family'] for r in rs}
            extra=[x for x in specs(s,True) if x['family'] not in existing]
            if extra:
                print('WIDE_BARRIER_SHAPE_MAP',s['tag'],len(extra),flush=True)
                r96+=mapping(s,'wide-fallback',extra,96);rs+=mapping(s,'wide-fallback',extra,192)
            rs=refine(s,rs,'wide-refine')
        choices=[];byid={r['id']:r for r in rs}
        for hole in range(2):
            candidates=strong([r for r in rs if r['family'].startswith(f'h{hole}-')])
            narrow=[x for x in candidates if float(x['a_outer'])/float(x['a_inner'])<=1.03]
            if narrow:candidates=narrow
            if candidates:
                chosen=min(candidates,key=lambda x:(float(x['a_outer'])/float(x['a_inner']),
                    max(float(byid[x[k]]['theta_rms']) for k in ('inner','outer'))))
            else:chosen=select(rs,hole)
            if not chosen:choices.append(dict(hole=hole,bracket=None,root=None));continue
            mean_pairs=[x for x in b.brackets([r for r in rs if r['family']==chosen['family']])
                if x['theta_inner_mean']<0 and x['theta_outer_mean']>0 and
                float(x['a_inner'])>=float(chosen['a_inner']) and float(x['a_outer'])<=float(chosen['a_outer'])]
            pair=min(mean_pairs,key=lambda x:float(x['a_outer'])/float(x['a_inner'])) if mean_pairs else chosen
            sd=d/f'hole{hole}-root-final';sd.mkdir(parents=True,exist_ok=True)
            write(sd/'bracket.csv',[chosen]);write(sd/'mean-root-pair.csv',[pair])
            surface=[dict(id=f'h{hole}-root-final',family=f'h{hole}-root-final',centre=chosen['centre'],
                a=pair['mean_linear_a_root'],c=pair['mean_linear_c_root'])]
            root=mapping(s,f'hole{hole}-root-final',surface,192)[0]
            root96=mapping(s,f'hole{hole}-root-final',surface,96)[0]
            endpoints=[dict(id=chosen[k],family=chosen['family'],centre=chosen['centre'],
                a=chosen['a_'+k],c=chosen['c_'+k]) for k in ('inner','outer')]
            checked=mapping(s,f'hole{hole}-final-endpoints',endpoints,96)
            if chosen['status']=='POINTWISE_NEGATIVE_TO_POSITIVE':
                assert float(checked[0]['theta_max'])<0 and float(checked[1]['theta_min'])>0
            else:assert float(checked[0]['theta_mean'])<0 and float(checked[1]['theta_mean'])>0
            assert root['status']==root96['status']=='RESOLVED'
            (sd/'finder-params.txt').write_text(w.replace((sd/'N192/params.txt').read_text(),dict(
                map_surfaces=sd/'N192/surfaces.csv',map_points=48,map_find='true',map_finder_bracket=sd/'bracket.csv',
                map_find_allow_average='true' if chosen['status']=='AVERAGE_ONLY_WEAKER' else 'false',
                map_find_seconds=235,map_find_max_updates=100000,map_find_floor_window=100001,
                map_find_chase=1,map_find_thresholds='1e-7 1e-10 1e-12')))
            choices.append(dict(hole=hole,bracket=chosen,root=root,seed_dir=str(sd),mean_root_pair=pair,
                barrier_fractional_width=float(chosen['a_outer'])/float(chosen['a_inner'])-1,
                angular_root_N96=root96))
        atomic_json(d/'selection.json',dict(time=s['time'],configuration=s['configuration'],selections=choices))
        families={x['bracket']['family'] for x in choices if x['bracket']}
        write(HERE/f't25-early-{s["tag"]}.csv',[r for r in rs if r['family'] in families or r['family'].endswith('-o+0-q1')])
        write(HERE/f't25-early-{s["tag"]}-brackets.csv',b.brackets(rs));write(HERE/f't25-early-{s["tag"]}-barriers.csv',strong(rs))
        write(d/'map-N192.csv',rs);write(d/'map-N96.csv',r96)
        atomic_json(d/'improved.json',dict(status='COMPLETE',physical_finder_probes=0))
        print('FINAL_BARRIERS',s['tag'],[(x['hole'],x['barrier_fractional_width'] if x['bracket'] else None) for x in choices],flush=True)
        summarize()
def summarize():
    rows=[];receipts=[]
    for s in samples():
        p=ROOT/s['tag']/'selection.json'
        if not p.exists():continue
        d=json.loads(p.read_text())
        for x in d['selections']:
            r=x['root'];br=x['bracket']
            out=dict(time=s['time'],configuration=s['configuration'],hole=x['hole'],checkpoint=str(s['cp']),
                bracket_status=br['status'] if br else 'NO_RESOLVED_BRACKET',
                r_cyl=r['a'] if r else '',r_axial=r['c'] if r else '',cells=r['min_extent_over_dx'] if r else '',
                r_cyl_finest_cells=float(r['a'])/H if r else '',r_axial_finest_cells=float(r['c'])/H if r else '',
                local_dx_min=r['dx_min'] if r else '',local_dx_max=r['dx_max'] if r else '',
                A_trial=r['A'] if r else '',Q_trial=r['Q'] if r else '',lapse_mean=r['lapse_mean'] if r else '',
                chi_mean=r['chi_mean'] if r else '',K_mean=r['K_mean'] if r else '',K_minus_2Theta_mean=r['K_minus_2Theta_mean'] if r else '',
                theta_min=r['theta_min'] if r else '',theta_max=r['theta_max'] if r else '',theta_rms=r['theta_rms'] if r else '',
                centre=r['centre'] if r else '',a_inner=br['a_inner'] if br else '',a_outer=br['a_outer'] if br else '',
                finder_status='PENDING' if br else 'NOT_RUN_NO_BRACKET',finder_squared='',finder_A='',finder_Q='',finder_seconds='')
            if br:
                fp=Path(x['seed_dir'])/'finder/finder.csv'
                if fp.exists() and (fp.parent/'done.exit').exists():
                    stage=read(fp)[-1];out.update(finder_status=stage['status'],finder_squared=stage['expansion_squared'],
                        finder_A=stage['A'],finder_Q=stage['Q'],finder_seconds=stage['seconds'])
            rows.append(out)
    if rows:write(HERE/'t25-early-time-table.csv',rows)
    for p in ROOT.rglob('*.resources.json'):
        if 'build' in p.parts:continue
        r=json.loads(p.read_text());receipts.append(dict(directory=r['directory'],process=r['process'],returncode=r['returncode'],
            peak_rss_bytes=r['peak_rss_bytes'],wall_seconds=r['wall_seconds'],gate_reason=r['gate_reason']))
    if receipts:write(HERE/'t25-early-receipts.csv',receipts)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('operation',choices=['initialize','maps','finders','summarize','all']);a=p.parse_args()
    if a.operation in ('initialize','all'):initialize()
    if a.operation in ('maps','all'):maps()
    if a.operation in ('finders','all'):finders()
    if a.operation=='summarize':summarize()
