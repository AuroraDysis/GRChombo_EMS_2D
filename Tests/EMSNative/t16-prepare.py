#!/usr/bin/env python3
"""T16 registration, measured builds, t=0 controls and native audit."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, json, math, os, subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
ROOT=Path('/private/tmp/ems-t16');PY='/Users/auroradysis/miniconda3/bin/python'
PROFILE=Path('/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-echobin/artifacts/echo-binary/a0.9-e8.trumpet')
V=.0523381404947;ETA=math.atanh(V);RH=.0060404520035922523
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
p=module('t13prepare',HERE/'t13-prepare.py')
def save(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def timed(name,directory,command):
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/(name+'.log')).open('wb') as log:
        subprocess.run([PY,str(HERE/'t16-run.py'),'--measure',name,'--directory',str(directory),'--',*map(str,command)],
            stdout=log,stderr=subprocess.STDOUT,check=True)
def build():
    subprocess.run(['make','-j1','all'],cwd=REPO/'Examples/EMS',check=True)
    b=module('t14prepare',HERE/'t14-prepare.py');b.ROOT=ROOT
    (HERE/'t16-native.hpp').write_text((HERE/'T11RHS.cpp').read_text().split('\nint main(')[0]+'\n')
    b.build('t16-launch.cpp','launch.ex')
    small=module('t7ebuild',HERE/'t7e-localize.py');small.TMP=ROOT
    small.build('t16-audit.cpp','audit.ex');small.build('t16-ranges.cpp','ranges.ex')
    # Reuse the current-state replay with the actual progenitor coupling.
    text=(HERE/'T13Replay.cpp').read_text().replace('-.8','-.9')
    (HERE/'t16-replay.cpp').write_text(text)
    small.build('t16-replay.cpp','replay.ex')
    census=(HERE/'T14Census.cpp').read_text().replace('    if(!p.emsbh_params.use_maximal_initial_lapse) return 3;\n','')
    (HERE/'t16-census.cpp').write_text(census);b.build('t16-census.cpp','census.ex')
def build_audit():
    b=module('t7ebuild',HERE/'t7e-localize.py');b.TMP=ROOT;b.build('t16-audit.cpp','audit.ex')
def prepare():
    ROOT.mkdir(exist_ok=True);rows=[];jobs=[]
    inputs=[PROFILE,ROOT/'baseline/launch.ex',Path('/Users/auroradysis/Workspace/EMS-deps/worktrees/consult-8-reply.md')]
    inputs+=[Path('/Users/auroradysis/Workspace/EMS/docs')/s for s in ('binary-ctt.md','binary-ctt-format.md')]
    (HERE/'t16-inputs.json').write_text(json.dumps([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),bytes=f.stat().st_size) for f in inputs],indent=2)+'\n')
    base=(HERE/'params/t13-E-mid-maximal.txt').read_text()
    for level in (13,14):
        h=.875/2**level
        for mode in ('sqrt','geometric'):
            name=f'L{level}-{mode}';d=ROOT/'evolution'/name
            for sub in ('plt','chk'):(d/sub).mkdir(parents=True,exist_ok=True)
            changes=dict(N1=768,N2=384,max_level=level,regrid_interval=' '.join(['0']*level),
                ems_data_path=PROFILE,ems_f2=-.9,bh_charge=.21329649676624957,boosted='true',boost_rapidity=ETA,
                ems_use_maximal_initial_lapse='false',ems_use_geometric_initial_lapse=str(mode=='geometric').lower(),
                RH_activate='false',RH_num_horizons=0,activate_mq_extraction=0,stop_time=10.5,max_steps=2,
                t14_dense_initial_tags='true',num_mass_extraction_radii=level,
                mass_extraction_levels=' '.join(map(str,range(1,level+1))),
                mass_extraction_radii=' '.join(str(112/2**l*.8) for l in range(1,level+1)))
            param=HERE/f't16-{name}.txt';param.write_text(p.replace(base,changes))
            jobs.append(dict(name=name,directory=str(d),command=[str(ROOT/'launch.ex'),str(param)]))
            rows.append(dict(run=name,h_M=h,Rh_over_h=RH/h,face_M=112/2**level,
                dt_M=h/4,steps=int(.002/(h/4)),common_steps=74 if level==13 else 148,
                native_t_stop_M=int(.002/(h/4))*(h/4),
                common_t_stop_M=74*.875/8192/4,estimated_peak_RSS_GB=1.3 if level==13 else 1.45,
                planned_wall_upper_s=120 if level==13 else 240,parameters=str(param)))
    save('t16-registration.csv',rows)
    (ROOT/'evolution/plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    print('Registered four isolated launch legs; upper estimate 720 s; launch separately after admission checks.',flush=True)
def audit():
    inp=ROOT/'audit-points.txt'
    inp.write_text(''.join(f'{.875/2**l:.17g} {lo+(hi-lo)*k/32:.17g} {ray} {domain}\n'
        for l in (13,14) for ray in (0,1) for domain,(lo,hi) in enumerate(((.00075,.0025),(.01,.02))) for k in range(33)))
    for mode in ('sqrt','geometric'):
        timed('audit-'+mode,ROOT/'audit',[ROOT/'audit.ex',PROFILE,mode,inp,HERE/f't16-audit-{mode}.csv'])
    names='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split()
    summary=[]
    for mode in ('sqrt','geometric'):
        allrows=list(csv.DictReader((HERE/f't16-audit-{mode}.csv').open()))
        for domain in (0,1):
            for ray in (0,1):
                for n,field in enumerate(names):
                    norms=[]
                    for l in (13,14):
                        a=[r for r in allrows if int(r['domain'])==domain and int(r['ray'])==ray and int(r['component'])==n and float(r['h'])==.875/2**l]
                        metrics={}
                        for key in ('physical_residual','total_residual','physical_rhs','translation','ko'):
                            v=[float(r[key]) for r in a];metrics[key+'_peak']=max(map(abs,v));metrics[key+'_RMS']=math.sqrt(sum(z*z for z in v)/len(v))
                        norms.append(metrics)
                    for norm in ('peak','RMS'):
                        e0=norms[0]['physical_residual_'+norm];e1=norms[1]['physical_residual_'+norm]
                        summary.append(dict(mode=mode,domain='W' if domain==0 else 'smooth',ray='axis' if ray==0 else 'diagonal',field=field,
                            identity=n not in range(13,18),norm=norm,residual_h=e0,residual_half=e1,
                            p=math.log2(e0/e1) if e0>0 and e1>0 else '',
                            total_residual_h=norms[0]['total_residual_'+norm],total_residual_half=norms[1]['total_residual_'+norm],
                            translation_half=norms[1]['translation_'+norm],rhs_half=norms[1]['physical_rhs_'+norm],KO_half=norms[1]['ko_'+norm]))
    save('t16-audit-summary.csv',summary)
    for r in summary:
        if r['mode']=='geometric' and r['domain']=='smooth' and r['identity'] and r['norm']=='RMS':print(r,flush=True)
def controls():
    rows=[];base=(Path('/private/tmp/ems-t13/controls/E-point-off/params.txt')).read_text()
    for transfer,boosted in [(t,b) for t in ('point','legacy') for b in (False,True)]:
        runs={}
        for variant in ('baseline','omitted','false'):
            name=transfer+'-'+('boosted' if boosted else 'unboosted')+'-'+variant;d=ROOT/'controls'/name
            for sub in ('plt','chk'):(d/sub).mkdir(parents=True,exist_ok=True)
            text=p.replace(base,dict(ems_data_path=PROFILE,ems_f2=-.9,bh_charge=.21329649676624957,
                boosted=str(boosted).lower(),boost_rapidity=ETA if boosted else 0.,RH_activate='false',RH_num_horizons=0,amr_transfer=transfer))
            if variant=='false':text=p.replace(text,dict(ems_use_geometric_initial_lapse='false'))
            param=d/'params.txt';param.write_text(text)
            exe=ROOT/'baseline/launch.ex' if variant=='baseline' else ROOT/'launch.ex'
            timed(name,d,[exe,param]);runs[variant]=d
        for variant in ('omitted','false'):
            for rel in ('plt/EMS_Plot_000000.2d.hdf5','chk/EMS_000002.2d.hdf5'):
                count,bits=p.compare(runs['baseline']/rel,runs[variant]/rel);assert bits==0
                rows.append(dict(transfer=transfer,boosted=boosted,variant=variant,file=rel,Float64_values=count,bit_mismatches=bits))
    save('t16-controls.csv',rows);print('Default identity: 16 comparisons PASS',flush=True)
def census():
    for l in (13,14):
        d=ROOT/'census'/f'L{l}'
        for sub in ('plt','chk'):(d/sub).mkdir(parents=True,exist_ok=True)
        timed(f'census-L{l}',d,[ROOT/'census.ex',HERE/f't16-L{l}-geometric.txt',HERE/f't16-census-L{l}.csv'])
    margins()
    print('Full native initialization census complete',flush=True)
def margins():
    margins=[]
    for l in (13,14):
        rows=list(csv.DictReader((HERE/f't16-census-L{l}.csv').open()));h=.875/2**l
        speed=max(float(r['max_shift_speed']) for r in rows if int(r['level'])==l)
        for ray,factor in (('axis',1.),('diagonal',math.sqrt(2))):
            gap=112/2**l-.0025/factor-8*h
            margins.append(dict(level=l,ray=ray,face_M=112/2**l,buffer_M=8*h,
                initial_max_family_speed=speed,provisional_bound=1.1,arrival_bound_M=gap/1.1,
                margin_beyond_window_M=gap/1.1-.002,condition='conditional characteristic estimate; FD/KO tails are not compact'))
            assert gap/1.1>.002 and speed<1.1
    save('t16-characteristic-margins.csv',margins)
def ranges():
    for l in (13,14):
        rows=list(csv.DictReader((HERE/f't16-census-L{l}.csv').open()))
        inp=ROOT/f'ranges-L{l}.txt';inp.write_text(''.join(' '.join(r[k] for k in ('level','box','x0','y0','x1','y1','h_M'))+'\n' for r in rows))
        timed(f'ranges-L{l}',ROOT/'ranges',[ROOT/'ranges.ex',PROFILE,inp,HERE/f't16-ranges-L{l}.csv'])
    print('Per-hole grid and near-companion ranges complete',flush=True)
def guards():
    base=(ROOT/'controls/point-boosted-omitted/params.txt').read_text();rows=[]
    tests={'geometric-boost':(dict(ems_use_geometric_initial_lapse='true'),True),
        'geometric-binary':(dict(ems_use_geometric_initial_lapse='true',binary='true',separation=32),True),
        'maximal-boost':(dict(ems_use_maximal_initial_lapse='true'),False),
        'maximal-binary':(dict(ems_use_maximal_initial_lapse='true',boosted='false',boost_rapidity=0,binary='true',separation=32),False),
        'both':(dict(ems_use_geometric_initial_lapse='true',ems_use_maximal_initial_lapse='true'),False),
        'legacy':(dict(ems_use_geometric_initial_lapse='true',ems_data_format='legacy_dat'),False),
        'RN':(dict(ems_use_geometric_initial_lapse='true',ems_not_rn='false'),False)}
    for name,(changes,accept) in tests.items():
        changes['check_params']='true';d=ROOT/'guards'/name;d.mkdir(parents=True,exist_ok=True)
        param=d/'params.txt';param.write_text(p.replace(base,changes))
        with (d/'run.log').open('wb') as f:
            q=subprocess.run([PY,str(HERE/'t16-run.py'),'--measure',name,'--directory',str(d),'--',str(ROOT/'launch.ex'),str(param)],stdout=f,stderr=subprocess.STDOUT)
        assert (q.returncode==0)==accept,(name,q.returncode)
        rows.append(dict(case=name,expected_accept=accept,returncode=q.returncode,result='PASS'))
    save('t16-guards.csv',rows)
def filter_check():
    import lzma,struct
    import numpy as np
    decoder=module('t7decoder',HERE/'t7-check.py')
    fr=next(decoder.frames(ROOT/'census/L13/t13-t7-stage-L13.xz'))
    d=ROOT/'filter-check';d.mkdir(exist_ok=True);raw=d/'input.raw';out=d/'output.xz'
    previous=None;expected=[]
    with raw.open('wb') as f:
        f.write(b'T7OP0002')
        for index,(phase,t) in enumerate(((50,0.),(50,.875/8192/8),(50,.875/8192/4),(3,.875/8192/8))):
            values=fr['values']+index*.125;bits=values.reshape(-1).view('u8')
            payload=bits if previous is None else bits^previous;previous=bits.copy()
            f.write(struct.pack('<11i',phase,13,fr['source'],-1,28,len(values),int(index>0),*fr['valid']))
            meta=fr['meta'].copy();meta[0]=t;f.write(meta.tobytes());f.write(fr['cells'].astype('<i4').tobytes())
            f.write(payload.view('u1').reshape(-1,8).T.copy().tobytes())
            if index!=1:
                keep=np.ones(len(values),bool) if phase!=50 else ((abs((fr['cells'][:,0]+.5)*meta[1]-336)<=.0035)&(abs((fr['cells'][:,1]+.5)*meta[1])<=.0035))
                expected.append((phase,t,fr['cells'][keep],values[keep]))
    with raw.open('rb') as src,out.open('wb') as dst:
        subprocess.run([PY,str(HERE/'t16-xz.py')],stdin=src,stdout=dst,check=True)
    count=0
    for got,want in zip(decoder.frames(out),expected,strict=True):
        phase,t,cells,values=want
        assert got['phase']==phase and got['meta'][0]==t and np.array_equal(got['cells'],cells)
        assert np.array_equal(got['values'].view('u8'),values.view('u8'));count+=values.size
    save('t16-filter-check.csv',[dict(kept_frames=3,dropped_half_clock_frames=1,Float64_values=count,bit_mismatches=0,result='PASS')])
    bindir=ROOT/'bin';bindir.mkdir(exist_ok=True);xz=bindir/'xz'
    xz.write_bytes((HERE/'t16-xz.py').read_bytes());xz.chmod(0o755)
    print('Output filter: bitwise round trip and skipped-clock delta chain PASS',flush=True)
def launch_plan():
    assert all(int(r['bit_mismatches'])==0 for r in csv.DictReader((HERE/'t16-controls.csv').open()))
    assert (HERE/'t16-filter-check.csv').exists()
    plan=ROOT/'evolution/plan.json';obj=json.loads(plan.read_text())
    for job in obj['jobs']:
        job['command']=[PY,str(HERE/'t16-launch-job.py'),str(ROOT/'launch.ex'),str(HERE/f"t16-{job['name']}.txt")]
    plan.write_text(json.dumps(obj,indent=2)+'\n')
    print(plan,flush=True)
if __name__=='__main__':
    {'build':build,'build-audit':build_audit,'prepare':prepare,'audit':audit,'controls':controls,'census':census,'margins':margins,'ranges':ranges,
     'guards':guards,'filter-check':filter_check,'launch-plan':launch_plan}[sys.argv[1]]()
