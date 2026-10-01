#!/usr/bin/env python3
"""T14 real production initialization census and serial launch registration."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, json, os, subprocess
from pathlib import Path
import h5py
import numpy as np

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
ROOT=Path('/private/tmp/ems-t14');OLD=Path('/private/tmp/ems-t13')
PY='/Users/auroradysis/miniconda3/bin/python'
spec=importlib.util.spec_from_file_location('prepare',HERE/'t13-prepare.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)

def save(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def timed(label,directory,command):
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/(label+'.log')).open('wb') as log:
        subprocess.run([PY,str(HERE/'t14-run.py'),'--measure',label,
            '--directory',str(directory),'--',*command],stdout=log,stderr=subprocess.STDOUT,check=True)

def build(source='T14Census.cpp',binary='census.ex'):
    ch=Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
    hdf=Path('/Users/auroradysis/Workspace/EMS-deps/hdf5-1.14')
    obj=REPO/'Examples/EMS/o/2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
    cmd=['/opt/homebrew/bin/g++-16','-O3','-std=c++17','-fopenmp','-fno-access-control']
    cmd+=['-D'+s for s in ('CH_SPACEDIM=2','CH_Darwin','NDEBUG','CH_USE_COMPLEX',
        'CH_USE_64','CH_USE_DOUBLE','CH_USE_HDF5','H5_USE_16_API','CH_USE_LAPACK',
        'CH_FORT_UNDERSCORE','CH_LANG_CC')]
    cmd+=['-I'+str(REPO/'Examples/EMS'),'-I'+str(hdf/'include')]
    cmd+=['-I'+str(d) for d in (REPO/'Source').rglob('*') if d.is_dir() and not d.name.startswith('.')]
    cmd+=['-I'+str(ch/'src'/d) for d in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
    cmd+=[str(HERE/source)]+[str(q) for q in sorted(obj.glob('*.o')) if q.stem!='Main_EMSBH2DBH']
    cmd+=['-L'+str(ch)]+['-l'+s+'2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
        for s in ('amrtimedependent','amrtools','boxtools','basetools')]
    cmd+=['-L'+str(hdf/'lib'),'-Wl,-rpath,'+str(hdf/'lib'),'-lhdf5','-lz',
        '-L/opt/homebrew/lib/gcc/16','-lgfortran','-lm','-lgomp','-framework','Accelerate',
        '-o',str(ROOT/binary)]
    (ROOT/(binary+'-build-command.json')).write_text(json.dumps(cmd,indent=2)+'\n')
    subprocess.run(cmd,check=True)

def prepare():
    jobs=[];rows=[]
    # Reuse the tested, frozen T13 executable, not a new production build.
    exe=ROOT/'launch.ex'
    production=REPO/'Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex'
    assert hashlib.sha256((OLD/'launch.ex').read_bytes()).digest()==hashlib.sha256(production.read_bytes()).digest()
    for name,level,dt in [('max13',13,.25),('max14',14,.25),('max14-half',14,.125)]:
        text=(HERE/'params'/f't13-stageC-E-mid-max{level}.txt').read_text()
        text=p.replace(text,dict(dt_multiplier=dt,t14_dense_initial_tags='true'))
        param=HERE/'params'/f't14-E-mid-{name}.txt';param.write_text(text)
        d=ROOT/'evolution'/name;d.mkdir(parents=True,exist_ok=True)
        for sub in ('plt','chk'):(d/sub).mkdir(exist_ok=True)
        dd=p.parameters(text);base=p.parameters((HERE/'params/t13-E-mid-maximal.txt').read_text())
        changes={k:[base.get(k),v] for k,v in dd.items() if base.get(k)!=v}
        assert set(changes)<= {'dt_multiplier','max_level','regrid_interval',
            'num_mass_extraction_radii','mass_extraction_levels','mass_extraction_radii','t14_dense_initial_tags'}
        assert dd['sigma']=='1' and dd['t2_guard_initial_data_after_t0']=='true'
        assert dd['amr_transfer']=='point' and dd['ems_use_maximal_initial_lapse']=='true'
        h=7/12/2**level;fine_dt=h*dt;steps=int(.002/fine_dt)
        jobs.append(dict(name=name,directory=str(d),command=[str(exe),str(param)]))
        rows.append(dict(run=name,level=level,h_M=h,dt0_M=7/12*dt,dt_finest_M=fine_dt,
            steps=steps,planned_native_endpoint_M=steps*fine_dt,
            common_endpoint_M=56*(7/12/4096)/4,
            wall_estimate_s=80 if level==13 else 156 if dt==.25 else 300,
            planned_wall_upper_s=120 if level==13 else 233 if dt==.25 else 450,
            estimated_peak_RSS_GB=2.15 if level==13 else 2.33,
            cap_RSS_GB=3,parameters=str(param),differences=json.dumps(changes)))
    plan=ROOT/'evolution/plan.json';plan.write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    save('t14-registration.csv',rows)
    print(plan,flush=True)

def census():
    for level in (13,14):
        d=ROOT/'census'/f'max{level}';d.mkdir(parents=True,exist_ok=True)
        for sub in ('plt','chk'):(d/sub).mkdir(exist_ok=True)
        timed(f'census-max{level}',d,[str(ROOT/'census.ex'),
            str(HERE/'params'/f't14-E-mid-max{level}.txt'),str(d/'boxes.csv')])
    verify()

def controls():
    rows=[]
    for member in ('E','ref'):
        for mode in ('legacy','point'):
            baseline=OLD/'controls'/f'{member}-{mode}-off'
            for dense in (False,True):
                name=f'{member}-{mode}-'+('dense' if dense else 'off')
                d=ROOT/'controls'/name;d.mkdir(parents=True,exist_ok=True)
                for sub in ('plt','chk'):(d/sub).mkdir(exist_ok=True)
                text=p.replace((baseline/'params.txt').read_text(),
                    dict(t14_dense_initial_tags=str(dense).lower()))
                param=d/'params.txt';param.write_text(text)
                timed('control-'+name,d,[str(ROOT/'launch.ex'),str(param)])
                for relative in ('plt/EMS_Plot_000000.2d.hdf5','chk/EMS_000002.2d.hdf5'):
                    count,bits=p.compare(baseline/relative,d/relative)
                    assert bits==0
                    rows.append(dict(member=member,mode=mode,dense=dense,
                        file=relative,Float64_values=count,bit_mismatches=bits,
                        condition='Tests factory against frozen T13 production; control finest level 1'))
    save('t14-controls.csv',rows)
    print('T14 Tests-only factory identity:',len(rows),'files, zero mismatches',flush=True)

def verify():
    reference={}
    with h5py.File('/Users/auroradysis/Workspace/EMS/.data/exp-0020/E-mid/plt/EMS_Plot_000000.2d.hdf5') as f:
        assert float(f.attrs['time'])==0
        for l in range(13):
            g=f[f'level_{l}'];reference[l]=(float(g.attrs['dx']),
                sorted(tuple(int(b[k]) for k in ('lo_i','lo_j','hi_i','hi_j')) for b in g['boxes'][:]))
    summary=[];margins=[];allboxes=[]
    for top in (13,14):
        rows=list(csv.DictReader((ROOT/'census'/f'max{top}'/'boxes.csv').open()))
        for l in range(top+1):
            group=[a for a in rows if int(a['level'])==l]
            assert group
            boxes=sorted(tuple(int(a[k]) for k in ('x0','y0','x1','y1')) for a in group)
            h=float(group[0]['h_M']);lo=np.min(np.array(boxes)[:,:2],axis=0)
            hi=np.max(np.array(boxes)[:,2:],axis=0)+1
            assert all(all(b[k]%8==0 for k in (0,1)) and all((b[k]+1)%8==0 for k in (2,3)) for b in boxes)
            count=sum(int(a['valid_cells']) for a in group)
            assert count==int(np.prod(hi-lo)),(top,l,'union has holes/overlap')
            old=l<=12
            if old:assert h==reference[l][0] and boxes==reference[l][1],(top,l,'old boxes changed')
            if top==14 and l==13:
                prior=[a for a in allboxes if a['run']=='max13' and int(a['level'])==13]
                assert boxes==sorted(tuple(int(a[k]) for k in ('x0','y0','x1','y1')) for a in prior)
            face=(hi[1])*h
            assert abs(lo[0]*h-(336-face))<1e-12 and abs(hi[0]*h-(336+face))<1e-12 and lo[1]==0
            if l:assert abs(face-112/2**l)<1e-12
            assert not any(int(a[k]) for a in group for k in ('nonfinite','chi_activations','lapse_activations'))
            summary.append(dict(run=f'max{top}',level=l,h_M=h,boxes=len(boxes),valid_cells=count,
                face_M=face,x_lo_M=lo[0]*h,x_hi_M=hi[0]*h,y_lo_M=0.,y_hi_M=face,
                old_boxes_bit_identical=old,block_aligned=True,
                chi_min=min(float(a['chi_min']) for a in group),lapse_min=min(float(a['lapse_min']) for a in group),
                max_shift_speed=max(float(a['max_shift_speed']) for a in group),
                max_light_speed=max(float(a['max_light_speed']) for a in group),
                max_lapse_speed=max(float(a['max_lapse_speed']) for a in group)))
            allboxes.extend(dict(run=f'max{top}',**a) for a in group)
            if l==top:
                for ray,v in [('axis',1.),('diagonal',np.sqrt(2))]:
                    buffer=7*h;gap=face-.0025/v-buffer
                    margins.append(dict(run=f'max{top}',ray=ray,face_M=face,
                        ghost_support_M=3*h,interpolation_plus_RHS_buffer_M=buffer,
                        gap_M=gap,speed_bound=1.1,arrival_bound_M=gap/1.1,
                        margin_after_stop_M=gap/1.1-.002,
                        actual_initial_max_speed=summary[-1]['max_shift_speed'],
                        condition='physical characteristic bound; verify on evolved states; FD/KO tails not compact'))
                    assert gap/1.1>.002 and summary[-1]['max_shift_speed']<1.1
    save('t14-boxes.csv',allboxes);save('t14-census.csv',summary);save('t14-characteristic-margins.csv',margins)
    # Real max13 initialization reached the same finest level through the
    # original factory. Compare its entire retained parent/finest native t=0
    # recorder windows, including ghosts, against the opt-in Tests factory.
    s=importlib.util.spec_from_file_location('decoder',HERE/'t7-check.py')
    decoder=importlib.util.module_from_spec(s);s.loader.exec_module(decoder)
    identity=[]
    for l in (12,13):
        native=ROOT/'census/native-max13'/f't13-t7-stage-L{l}.xz'
        test=ROOT/'census/max13'/f't13-t7-stage-L{l}.xz'
        count=bits=frames=0
        for a,b in zip(decoder.frames(native),decoder.frames(test),strict=True):
            assert a['phase']==b['phase']==50 and a['meta'][0]==b['meta'][0]==0
            assert np.array_equal(a['cells'],b['cells'])
            count+=a['values'].size
            bits+=int(np.count_nonzero(a['values'].view('u8')!=b['values'].view('u8')));frames+=1
        assert count>0 and bits==0
        identity.append(dict(level=l,frames=frames,Float64_values=count,bit_mismatches=bits,
            condition='real max13 t=0 retained full recorder windows, all 28 components and ghosts'))
    save('t14-census-identity.csv',identity)
    (ROOT/'census/verified.json').write_text(json.dumps(dict(old_levels=[0,12],
        actual_census=True,old_boxes_identical=True,block_alignment=8,floor_activations=0,
        faces=[.013671875,.0068359375],margins=margins),indent=2)+'\n')
    print(json.dumps(dict(margins=margins,identity=identity),indent=2),flush=True)

def manifest():
    resources=[]
    for q in sorted(ROOT.rglob('*.resources.json')):
        a=json.loads(q.read_text())
        sampled=max((v[0] for v in a.get('per_pid_peaks',{}).values()),default=0)
        peak=max(a['peak_rss_bytes'],sampled)
        resources.append(dict(process=a['process'],
            peak_rss_bytes=peak,peak_RSS_GB=peak/1e9,
            reported_peak_rss_bytes=a['peak_rss_bytes'],sampled_peak_rss_bytes=sampled,
            sampled_tree_GB=a['tree_peak_bytes']/1e9,wall_s=a['wall_seconds'],
            returncode=a['returncode'],gate_reason=a['gate_reason'],record=str(q)))
    if resources:save('t14-resources.csv',resources)
    targets=list(HERE.glob('t14-*'))+list(HERE.glob('T14*'))+[HERE/'t13-run.py',HERE/'README.md',ROOT/'launch.ex',ROOT/'census.ex']
    targets+=list((HERE/'params').glob('t14-*'))+list((HERE/'figures').glob('t14-*'))
    targets+=[OLD/'launch.ex',OLD/'replay.ex',HERE/'T13Replay.cpp',HERE/'t13-analyze.py',HERE/'t13-check.py',HERE/'t13-report.py',
        HERE/'t13-controls.csv',HERE/'params/t13-stageC-reading.md',
        Path('/Users/auroradysis/Workspace/EMS-deps/worktrees/consult-7-reply.md'),p.BASE,p.PROFILE]
    targets+=list((ROOT/'census').rglob('*')) if (ROOT/'census').exists() else []
    targets+=list(ROOT.glob('*.json'))+list(ROOT.glob('*.time'))+list(ROOT.glob('*.log'))
    targets+=list((ROOT/'controls').rglob('*')) if (ROOT/'controls').exists() else []
    # Running output streams are pinned after completion, never hash moving files.
    for plan in (ROOT/'evolution').rglob('plan.json'):
        targets.append(plan)
        targets+=list(plan.parent.glob('*.json'))+list(plan.parent.glob('*.exit'))
        targets+=list(plan.parent.glob('*.log'))+list(plan.parent.glob('*.time'))
        for job in json.loads(plan.read_text())['jobs']:
            d=Path(job['directory'])
            if (d/'done.exit').exists():targets+=list(d.rglob('*'))
    targets+=list((ROOT/'evolution/max14-half.stopped-ceiling').rglob('*'))
    targets+=list((ROOT/'runner-check').rglob('*'))
    targets+=list((ROOT/'analysis').rglob('*'))
    targets += [HERE/'COMMIT-MANIFEST-T13.txt',HERE/'params/t13-E-mid-maximal.txt']
    for name in ('maximal','maximal-half'):
        targets+=list((OLD/'analysis').glob(name+'-*profiles.npz'))
    targets+=list((OLD/'analysis/maximal-stages').glob('*'))
    targets += [OLD/'analysis/maximal-meta.npz',OLD/'analysis/maximal-q.bin',OLD/'analysis/maximal-initial.npz']
    header=['T14 Stage C; baseline HEAD 949a744; no production C++ change or commit.',
        'Per-process RSS/footprint hard gate 3,000,000,000 bytes; sampled active child tree gate also 3,000,000,000 bytes; serial OMP=2, BLAS=1, xz T1.',
        'All numerical processes invoke /usr/bin/time -l; sandbox kern.clockrate failure uses child wait4 RUSAGE_CHILDREN plus live libproc sampling, as T13.',
        'Actual real initialization census before evolution; current-field recorder only after t=0. Input static data read-only.',
        'Original queue output gate 5,700,000,000 bytes across /private/tmp/ems-t14. Half-step restart: 6,000,000,000 bytes scoped to evolution/max14-half; partial output retained as max14-half.stopped-ceiling. T13 untouched.',
        'Completed T14: endpoint contraction passes; full registered history is inconclusive (53 sampling-unqualified nonzero entries); no Stage D draft or production admission. Full analysis in README T14.',
        'All completed recorder streams, replay caches, removed-payload hashes, scientific tables/figures and every measured process peak are pinned. T13 baseline caches are read-only; attribution now streams one stage at a time.',
        'Gate-killed wrappers can underreport wait4 RSS: resource CSV retains raw measurement and sampled RSS, using their maximum for peak.']
    header += [f"RESOURCE {a['process']} peak_RSS_bytes={a['peak_rss_bytes']} returncode={a['returncode']} record={a['record']}" for a in resources]
    header.append('SHA256  bytes  path')
    for q in sorted(set(targets)):
        if not q.is_file():continue
        if q.name.startswith(os.environ.get('T13_CURRENT_MEASURE','NO_ACTIVE_MEASUREMENT')+'.'):continue
        digest=hashlib.sha256()
        with q.open('rb') as f:
            for chunk in iter(lambda:f.read(1<<20),b''):digest.update(chunk)
        header.append(f'{digest.hexdigest()}  {q.stat().st_size}  {q}')
    (HERE/'COMMIT-MANIFEST-T14.txt').write_text('\n'.join(header)+'\n')
    print('T14 measured processes',len(resources),flush=True)

if __name__=='__main__':
    ROOT.mkdir(exist_ok=True)
    {'build':build,'build-launch':lambda:build('T14Launch.cpp','launch.ex'),
        'prepare':prepare,'census':census,'controls':controls,'verify':verify,'manifest':manifest}[sys.argv[1]]()
