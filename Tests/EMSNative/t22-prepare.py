#!/usr/bin/env python3
"""Bounded local T22 source snapshots/builds and serial evidence plans. No SSH."""
import argparse, hashlib, importlib.util, json, os, re, shutil, subprocess, sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];ROOT=Path('/private/tmp/ems-t22')
PY=sys.executable
CONTRACT=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0023/submissions/exp-0023')
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def replace(text,changes):
    for k,v in changes.items():
        pat=r'^'+re.escape(k)+r'\s*=.*$';line=f'{k} = {v}'
        text=re.sub(pat,line,text,flags=re.M) if re.search(pat,text,re.M) else text+'\n'+line+'\n'
    return text
def snapshot():
    ROOT.mkdir(exist_ok=True)
    for name,source in [('baseline',Path('/private/tmp/ems-t21/src')),('src',REPO)]:
        dst=ROOT/name
        if dst.exists():continue
        dst.mkdir()
        # T21's baseline object/source archive is already qualified. Keep both
        # snapshots separate; do not change source under a running probe.
        for directory in ('Source','Examples/EMS','Tests/EMSNative'):
            for p in (source/directory).rglob('*'):
                if p.is_file() and p.suffix in ('.hpp','.cpp','.H','.inc'):
                    out=dst/p.relative_to(source);out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
        for p in HERE.glob('t22-*'):
            if p.is_file() and p.suffix in ('.hpp','.cpp'):
                shutil.copy2(p,dst/'Tests/EMSNative'/p.name)
    obj=REPO/'Examples/EMS/o/2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
    (ROOT/'reused-objects.json').write_text(json.dumps({str(p):sha(p) for p in obj.glob('*.o')},indent=2)+'\n')
    records={}
    for directory in ('src','baseline'):
        records[directory]={str(p.relative_to(ROOT/directory)):sha(p) for p in (ROOT/directory).rglob('*') if p.is_file()}
    (ROOT/'source-snapshots.json').write_text(json.dumps(records,indent=2)+'\n')
def build():
    snapshot();ch=Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib');hdf=Path('/Users/auroradysis/Workspace/EMS-deps/hdf5-1.14')
    base=['/opt/homebrew/bin/g++-16','-O3','-std=c++17','-fopenmp','-fno-access-control']
    base+=['-D'+s for s in ('CH_SPACEDIM=2','CH_Darwin','NDEBUG','CH_USE_COMPLEX','CH_USE_64','CH_USE_DOUBLE','CH_USE_HDF5','H5_USE_16_API','CH_USE_LAPACK','CH_FORT_UNDERSCORE','CH_LANG_CC')]
    obj=REPO/'Examples/EMS/o/2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
    objects=[str(p) for p in sorted(obj.glob('*.o')) if p.stem not in ('Main_EMSBH2DBH','EMSBH2DLevel')]
    libs=['-L'+str(ch)]+['-l'+s+'2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC' for s in ('amrtimedependent','amrtools','boxtools','basetools')]
    libs+=['-L'+str(hdf/'lib'),'-Wl,-rpath,'+str(hdf/'lib'),'-lhdf5','-lz','-L/opt/homebrew/lib/gcc/16','-lgfortran','-lm','-lgomp','-framework','Accelerate']
    commands=[]
    def call(cmd):
        commands.append(cmd);(ROOT/'build-commands.json').write_text(json.dumps(commands,indent=2)+'\n');subprocess.run(cmd,check=True)
    for name in ('baseline','src'):
        src=ROOT/name
        flags=base+['-I'+str(src/'Examples/EMS'),'-I'+str(hdf/'include')]
        flags+=['-I'+str(d) for d in (src/'Source').rglob('*') if d.is_dir() and not d.name.startswith('.')]
        flags+=['-I'+str(ch/'src'/d) for d in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
        call(flags+[str(src/'Tests/EMSNative/t22-audit.cpp'),*libs,'-o',str(ROOT/(name+'-generic.ex'))])
        if name=='baseline':
            # Qualified pre-selector native level, exact recorded object.
            level=Path('/private/tmp/ems-t21/level-recording.o')
            call(flags+[str(src/'Tests/EMSNative/T14Launch.cpp'),str(level),*objects,*libs,'-o',str(ROOT/'baseline-dense.ex')])
        else:
            level=ROOT/'level.o'
            call(flags+['-c',str(src/'Examples/EMS/EMSBH2DLevel.cpp'),'-o',str(level)])
            moving=ROOT/'moving.o'
            call(flags+['-c',str(src/'Examples/EMS/EMSBH2DMovingGauge.cpp'),'-o',str(moving)])
            for source,binary in [('T14Launch.cpp','dense.ex'),('t22-initial-audit.cpp','initial.ex'),('t21-native-replay.cpp','replay.ex')]:
                call(flags+[str(src/'Tests/EMSNative'/source),str(level),str(moving),*objects,*libs,'-o',str(ROOT/binary)])
    print('T22_BUILD_COMPLETE',flush=True)
def plans():
    seed=Path('/private/tmp/ems-t13/controls/E-point-off/params.txt').read_text()
    common=dict(ems_use_maximal_initial_lapse='true',max_steps=1,checkpoint_interval=1,plot_interval=1,t13_launch_stop_time=0,t6_interface_diagnostics='false',t7_diagnostics='false',t2_guard_initial_data_after_t0='true',amr_transfer='point',verbosity=0)
    recorder=(HERE/'t21-recording-params.txt').read_text()
    recorder=replace(recorder,dict(t21_levels='0 1',t21_puncture_radius=.2,t21_profile_radius=1))
    jobs=[]
    def job(name,directory,cmd):
        directory.mkdir(parents=True,exist_ok=True);jobs.append(dict(name=name,directory=str(directory),command=list(map(str,cmd))))
    for case in ('baseline','omitted','explicit','capture','moving','moving-capture'):
        d=ROOT/'controls'/case
        for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
        text=replace(seed,common)
        if case=='explicit':text=replace(text,dict(ems_gauge='experimental'))
        if case.startswith('moving'):text=replace(text,dict(ems_gauge='moving_puncture',lapse_advec_coeff=0,shift_advec_coeff=0,shift_Gamma_coeff=1))
        if case in ('capture','moving-capture'):
            cap=recorder
            if case=='moving-capture':cap=cap.replace('experimental','moving_puncture').replace('integrated_offset_decaying_0p1','differential_gamma_driver')
            text+='\n'+replace(cap,dict(t21_rhs_prefix=d/'native'))
        p=d/'params.txt';p.write_text(text)
        exe=ROOT/('baseline-dense.ex' if case=='baseline' else 'dense.ex')
        job(case,d,[exe,p])
    for source in ('baseline','src'):
        for gauge in ('experimental','moving_puncture'):
            d=ROOT/'generic'/(source+'-'+gauge)
            job(source+'-'+gauge,d,[ROOT/(source+'-generic.ex'),gauge,d/'native.bin'])
    for mode,checkpoint,gauge,cap in [('restart-same','moving','moving_puncture',None),('restart-legacy','baseline','experimental',None),('restart-cross','baseline','moving_puncture',None),('restart-cross-tagged','moving','experimental',None),('restart-same-experimental','explicit','experimental',None)]:
        d=ROOT/'controls'/mode;d.mkdir(parents=True,exist_ok=True)
        for sub in ('chk','plt'):(d/sub).mkdir(exist_ok=True)
        text=replace(seed,dict(common,restart_file=ROOT/f'controls/{checkpoint}/chk/EMS_000001.2d.hdf5',ems_gauge=gauge,ems_data_path='/T22_STATIC_FILE_MUST_NOT_BE_READ',max_steps=1))
        if gauge=='moving_puncture':text=replace(text,dict(lapse_advec_coeff=0,shift_advec_coeff=0,shift_Gamma_coeff=1))
        p=d/'params.txt';p.write_text(text)
        job(mode,d,[PY,HERE/'t22-verify.py','restart-probe',ROOT/'dense.ex',p,mode])
    job('quick-qualification',ROOT/'quick-qualification',[PY,HERE/'t22-verify.py','quick'])
    d=ROOT/'quick-pipeline';d.mkdir(exist_ok=True);(d/'plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    # High hierarchy initialization + sampled native audit: no time advancement.
    jobs=[]
    for gauge in ('experimental','moving_puncture'):
        d=ROOT/'equilibrium'/gauge
        for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
        seed=(CONTRACT/'production/params-E-high.txt').read_text()
        changes=dict(ems_gauge=gauge,ems_data_path='/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet',max_steps=0,plot_interval=-1,checkpoint_interval=1,verbosity=0,RH_activate='false',activate_mq_extraction=0,t22_audit_levels='12 13 14',t22_audit_puncture_cells=4,t22_audit_radii='.0004 .001 .003 .0055 .0060286129392158804 .0075 .012 .02')
        if gauge=='moving_puncture':changes.update(lapse_advec_coeff=0,shift_advec_coeff=0,shift_Gamma_coeff=1)
        p=d/'params.txt';p.write_text(replace(seed,changes));job('equilibrium-'+gauge,d,[ROOT/'initial.ex',p])
    job('equilibrium-qualification',ROOT/'equilibrium-qualification',[PY,HERE/'t22-verify.py','equilibrium'])
    # Native T13 window on E-mid: sqrt(chi) is NOT changed; both use alpha_K.
    # Two gauges plus dt/2, clocks matched at 56 and 112 finest steps.
    seed=(HERE/'params/t13-E-mid-maximal.txt').read_text()
    for gauge in ('experimental','moving_puncture'):
        for half in (False,True):
            name=gauge+('-half' if half else '')
            d=ROOT/'launch'/name
            for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
            changes=dict(ems_gauge=gauge,ems_use_maximal_initial_lapse='true',dt_multiplier=.125 if half else .25,ems_data_path='/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet',checkpoint_interval=-1,plot_interval=-1,verbosity=0,RH_activate='false',activate_mq_extraction=0)
            if gauge=='moving_puncture':changes.update(lapse_advec_coeff=0,shift_advec_coeff=0,shift_Gamma_coeff=1)
            p=d/'params.txt';p.write_text(replace(seed,changes));job(name,d,[ROOT/'dense.ex',p])
    job('launch-qualification',ROOT/'launch-qualification',[PY,HERE/'t22-verify.py','launch'])
    job('seal',ROOT/'seal',[PY,HERE/'t22-seal.py'])
    d=ROOT/'evidence-pipeline';d.mkdir(exist_ok=True);(d/'plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    print(ROOT/'quick-pipeline/plan.json');print(ROOT/'evidence-pipeline/plan.json')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['build','plans','snapshot']);a=ap.parse_args();globals()[a.mode]()
