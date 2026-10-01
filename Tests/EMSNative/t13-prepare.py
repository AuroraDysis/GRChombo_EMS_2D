#!/usr/bin/env python3
"""Frozen-parameter T13 launch registration and small bit-identity control."""
import sys
sys.dont_write_bytecode=True
import csv,hashlib,json,os,re,subprocess
from pathlib import Path
import h5py
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
TMP=Path('/private/tmp/ems-t13')
PY='/Users/auroradysis/miniconda3/bin/python'
EXE=REPO/'Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex'
BASE=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0020/submissions/exp-0020/params-E-mid.txt')
PROFILE=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet')

def parameters(text):
    return {m[1]:m[2].strip() for s in text.splitlines()
            if (m:=re.match(r'^\s*([\w]+)\s*=\s*([^#]*)',s))}

def replace(text,changes):
    for key,value in changes.items():
        pattern=r'(?m)^\s*'+re.escape(key)+r'\s*=.*$'
        line=f'{key} = {value}'
        text=re.sub(pattern,line,text) if re.search(pattern,text) else text+'\n'+line+'\n'
    return text

def prepare():
    assert hashlib.sha256(PROFILE.read_bytes()).hexdigest().startswith('2a8de074')
    exe=TMP/'launch.ex';exe.write_bytes(EXE.read_bytes());exe.chmod(0o755)
    original=BASE.read_text();jobs=[];rows=[]
    common=dict(ems_data_path=str(PROFILE),checkpoint_interval=-1,plot_interval=-1,
                data_subpath='.',t13_launch_stop_time=.002)
    for name,maximal,dt in [('original',False,.25),('maximal',True,.25),
                             ('original-half',False,.125),('maximal-half',True,.125)]:
        d=TMP/'evolution'/name;d.mkdir(parents=True,exist_ok=True)
        for sub in ('chk','plt'): (d/sub).mkdir(exist_ok=True)
        changes=dict(common,dt_multiplier=dt,ems_use_maximal_initial_lapse=str(maximal).lower())
        p=HERE/'params'/f't13-E-mid-{name}.txt';p.write_text(replace(original,changes))
        diff={k:(parameters(original).get(k),v) for k,v in parameters(p.read_text()).items()
              if parameters(original).get(k)!=v}
        assert set(diff)<=set(changes),diff
        assert parameters(p.read_text())['sigma']=='1'
        assert parameters(p.read_text())['t2_guard_initial_data_after_t0']=='true'
        command=[str(exe),str(p)]
        jobs.append(dict(name=name,directory=str(d),command=command))
        h=7/12/4096;fine_dt=h*dt;steps=int(.002/fine_dt)
        rows.append(dict(run=name,h12_M=h,dt12_M=fine_dt,steps=steps,
                         requested_stop_M=.002,actual_stop_M=steps*fine_dt,
                         maximal=maximal,parameters=str(p),differences=json.dumps(diff)))
    plan=TMP/'evolution/plan.json';plan.write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    with (HERE/'t13-registration.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print(plan,flush=True)

def compare(a,b):
    count=0;mismatch=0
    with h5py.File(a) as x,h5py.File(b) as y:
        def visit(name,obj):
            nonlocal count,mismatch
            if not isinstance(obj,h5py.Dataset):return
            assert name in y and obj.shape==y[name].shape and obj.dtype==y[name].dtype,name
            aa=obj[...];bb=y[name][...]
            if aa.dtype.kind=='f':
                count+=aa.size;mismatch+=np.count_nonzero(aa.view('u8')!=bb.view('u8'))
            else:assert np.array_equal(aa,bb),name
        x.visititems(visit)
    return int(count),int(mismatch)

def preflight():
    rows=[];n=[];storage=0;vmax=0.
    path=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0020/E-mid/plt/EMS_Plot_000000.2d.hdf5')
    with h5py.File(path) as f:
        assert float(f.attrs['time'])==0.
        for level in range(13):
            g=f[f'level_{level}'];h=float(g.attrs['dx'])
            boxes=np.array([[b[k] for k in ('lo_i','lo_j','hi_i','hi_j')] for b in g['boxes'][:]])
            widths=boxes[:,2:]-boxes[:,:2]+1
            cells=int(np.prod(widths,axis=1).sum());halo=int(np.prod(widths+6,axis=1).sum())
            parent=int(np.prod(widths//2+8,axis=1).sum()) if level else 0
            # Native states/diagnostics, the recursive RK rhs+temporary, point
            # coarsened work plus RK Taylor coefficients and diff12. The final
            # budget adds allocator/stencil/diagnostic overhead explicitly.
            payload=halo*(4*28+18)*8+parent*(8*28)*8
            storage+=payload;n.append(cells)
            face=float((boxes[:,3].max()+1)*h)
            rows.append(dict(level=level,h_M=h,boxes=len(boxes),valid_cells=cells,
                cells_with_halo=halo,face_M=face,estimated_array_bytes=payload))
            if level==12:
                offsets=g['data:offsets=0'][:];data=g['data:datatype=0']
                for i,width in enumerate(widths):
                    size=int(np.prod(width+6));start=int(offsets[i])
                    bx=data[start+14*size:start+15*size];by=data[start+15*size:start+16*size]
                    vmax=max(vmax,float((1+np.hypot(bx,by)).max()))
        assert abs(rows[12]['face_M']-.02734375)<1e-12
    with (HERE/'t13-box-census.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    full=sum(c*2**l for l,c in enumerate(n));cost=[]
    for name,steps in [('dt',56),('half',112)]:
        count=sum(c*max(1,int(np.ceil(steps/2**(12-l)))) for l,c in enumerate(n))
        # Cluster observed E-mid: 162 s, 128 cores. Serial 2-thread scaling is
        # an estimate, not a claimed local performance measurement.
        estimate=162*64*count/full
        cost.append(dict(policy=name,finest_steps=steps,cell_steps=count,
            full_coarse_cell_steps=full,cost_fraction=count/full,
            core_scaling_wall_s=estimate,planned_wall_low_s=estimate*.5,
            planned_wall_high_s=estimate*2+180,estimated_peak_RSS_GB=storage*1.35/1e9,
            fixed_array_bytes=storage,initial_max_shift_gauge_speed=vmax))
    with (HERE/'t13-resource-plan.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=cost[0]);w.writeheader();w.writerows(cost)
    assert storage*1.35<6_000_000_000,cost
    margins=[]
    for level in (12,13,14):
        h=7/12/2**level;face=112/2**level
        for ray in ('axis','diagonal'):
            support=7*h;gap=face-.0025/(1 if ray=='axis' else np.sqrt(2))-support
            margins.append(dict(level=level,ray=ray,face_M=face,combined_stencil_buffer_M=support,
                face_normal_gap_M=gap,speed_bound=1.1,travel_time_bound_M=gap/1.1,
                time_margin_M=gap/1.1-.002,condition='initial shift-speed max measured; 1.1 bound must be checked on evolved recorder states'))
    with (HERE/'t13-characteristic-margins.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=margins[0]);w.writeheader();w.writerows(margins)
    print(json.dumps(dict(cost=cost,margins=margins),indent=2),flush=True)

def controls():
    rows=[]
    for leg in ('E','ref'):
        base=(HERE/'params'/f't4-control-{leg}.txt').read_text()
        for mode in ('legacy','point'):
            for variant in ('baseline','off','capture'):
                if mode=='legacy' and variant=='capture':continue
                d=TMP/'controls'/f'{leg}-{mode}-{variant}';d.mkdir(parents=True,exist_ok=True)
                for sub in ('chk','plt'): (d/sub).mkdir(exist_ok=True)
                text=replace(base,dict(amr_transfer=mode if mode=='point' else 'legacy',
                        num_plot_vars=28,plot_vars='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'))
                # Existing mode spelling is "legacy"? The parser accepts "legacy".
                if variant=='capture':text=replace(text,dict(t13_launch_stop_time=1.))
                p=d/'params.txt';p.write_text(text)
                exe=TMP/'baseline-430470e.ex' if variant=='baseline' else TMP/'launch.ex'
                with (d/'run.log').open('wb') as log:
                    rc=subprocess.run([PY,str(HERE/'t13-run.py'),'--measure',d.name,
                        '--directory',str(d),'--',str(exe),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
                assert rc==0,(d,rc)
            b=TMP/'controls'/f'{leg}-{mode}-baseline'
            for variant in ('off','capture'):
                if mode=='legacy' and variant=='capture':continue
                d=TMP/'controls'/f'{leg}-{mode}-{variant}'
                for filename in ('plt/EMS_Plot_000000.2d.hdf5','chk/EMS_000002.2d.hdf5'):
                    n,m=compare(b/filename,d/filename)
                    rows.append(dict(leg=leg,mode=mode,variant=variant,file=filename,
                                     Float64_values=n,bit_mismatches=m,pass_control=m==0))
                    assert m==0,rows[-1]
    with (HERE/'t13-controls.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print('bit-identity controls:',rows,flush=True)

def early_control():
    d=TMP/'controls/early-stop';d.mkdir(parents=True,exist_ok=True)
    text=replace((HERE/'params/t4-control-E.txt').read_text(),
        dict(amr_transfer='point',t13_launch_stop_time=.022,plot_interval=-1,checkpoint_interval=-1,
             t2_guard_initial_data_after_t0='true'))
    p=d/'params.txt';p.write_text(text)
    with (d/'run.log').open('wb') as log:
        subprocess.run([PY,str(HERE/'t13-run.py'),'--measure','early-stop',
            '--directory',str(d),'--',str(TMP/'launch.ex'),str(p)],
            stdout=log,stderr=subprocess.STDOUT,check=True)
    actual=float((d/'t13-stop.csv').read_text().splitlines()[1])
    assert actual==.015625 and not list(d.glob('**/*.hdf5'))
    assert 'T13 clean native stop' in (d/'run.log').read_text()
    print('clean unsynchronized early stop:',actual,'M, no HDF5 output',flush=True)

def manifest():
    resources=[]
    for p in sorted(TMP.rglob('*.resources.json')):
        q=json.loads(p.read_text())
        resources.append(dict(process=q['process'],wall_seconds=q['wall_seconds'],
            peak_rss_bytes=q['peak_rss_bytes'],peak_rss_GB=q['peak_rss_bytes']/1e9,
            total_sampled_tree_GB=q['tree_peak_bytes']/1e9,returncode=q['returncode'],
            time_returncode=q.get('time_returncode',''),gate_reason=q['gate_reason'],record=str(p)))
    with (HERE/'t13-resources.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=resources[0]);w.writeheader();w.writerows(resources)
    targets=[HERE/'README.md',HERE/'T13Replay.cpp',REPO/'Source/GRChomboCore/T13LaunchRecorder.hpp',
             REPO/'Source/GRChomboCore/GRAMRLevel.hpp',REPO/'Source/GRChomboCore/GRAMRLevel.cpp',
             REPO/'Examples/EMS/EMSBH2DLevel.hpp',REPO/'Examples/EMS/EMSBH2DLevel.cpp',
             REPO/'Examples/EMS/Main_EMSBH2DBH.cpp',BASE,PROFILE,TMP/'baseline-430470e.ex',TMP/'launch.ex',TMP/'replay.ex',
             TMP/'evolution/plan.json',Path('/Users/auroradysis/Workspace/EMS-deps/worktrees/consult-7-reply.md')]
    targets+=list(HERE.glob('t13-*'))+list((HERE/'params').glob('t13-*'))
    targets+=list(TMP.rglob('*.resources.json'))+list(TMP.rglob('*.child.json'))
    targets+=list(TMP.rglob('*.time'))
    targets+=list((TMP/'controls').rglob('*.xz'))
    targets+=list((TMP/'controls').rglob('*.csv'))
    targets+=list((HERE/'figures').glob('t13-*'))
    targets+=[p for p in (TMP/'evolution').rglob('*') if p.is_file()]
    targets+=[p for p in (TMP/'analysis').rglob('*') if p.is_file() and p.suffix in
              ('.bin','.npz','.sha256','.log')]
    targets+=[REPO/'Source/Cartoon/CCZ4Cartoon.impl.hpp',REPO/'Source/CCZ4/ExperimentalGauge.hpp',
              HERE/'T11RHS.cpp',HERE/'t7-check.py',REPO/'Source/GRChomboCore/T7OperationRecorder.hpp',
              REPO/'Source/TaggingCriteria/EMSExtractionTaggingCriterion.hpp']
    text=['T13 uncommitted diagnostics; baseline HEAD 430470e; no production equations/gauge/KO/transfer change.',
          'Resource limits: 6,000,000,000 RSS/footprint bytes per child; 8,000,000,000 total; OMP=2 and two single-thread xz streams.',
          '/usr/bin/time -l invoked for each command; sandbox denies kern.clockrate and returns 1 even after child success.',
          'Actual exit status and peak RSS use explicit child metadata and wait4 RUSAGE_CHILDREN; libproc monitors live RSS/footprint and tree total.',
          'Initial build time-l RSS was unavailable; the complete subsequent rebuild is measured in build.resources.json and the resource CSV.',
          'All four Stage B evolutions completed with done.exit=0; final runtime frames/logs and bounded analysis caches are retained and hashed.',
          'Stage C parameter files and plan are prepared ONLY; no Stage C evolution has run. New boxes remain a block-aligned prediction to census before advancing.',
          'Screen endpoint is the last native step inside 0.002 M: 0.001993815104166667 M (56 nominal / 112 control steps).',
          'Earlier measure peaks retained in tool output before repeat archiving was added: 265584640 and 252657664 bytes; final/repeated process JSONs are hashed below.',
          'The currently active manifest writer timer/log/metadata is excluded from hashes; its final measured metadata remains under analysis, and the prior writer is pinned.',
          'SHA256  bytes  path']
    for p in sorted(set(targets)):
        if not p.is_file():continue
        if p.parent==TMP/'analysis' and p.name.startswith(os.environ.get('T13_CURRENT_MEASURE','NO_ACTIVE_MEASUREMENT')+'.'):continue
        h=hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
        text.append(f'{h.hexdigest()}  {p.stat().st_size}  {p}')
    (HERE/'COMMIT-MANIFEST-T13.txt').write_text('\n'.join(text)+'\n')
    print('T13 manifest/resources refreshed:',len(resources),'measured commands',flush=True)

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='controls':controls()
    elif sys.argv[1]=='preflight':preflight()
    elif sys.argv[1]=='early-control':early_control()
    elif sys.argv[1]=='manifest':manifest()
    elif sys.argv[1]=='build-replay':
        import importlib.util
        s=importlib.util.spec_from_file_location('builder',HERE/'t7e-localize.py')
        m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
        m.TMP=TMP;m.build('T13Replay.cpp','replay.ex')
    else:raise ValueError(sys.argv[1])
