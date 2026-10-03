#!/usr/bin/env python3
"""Isolated commit/build matrix and bounded native diagnostic-ghost probes."""
import argparse,hashlib,json,os,re,shutil,subprocess,sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];ROOT=Path('/private/tmp/ems-t23');PY=sys.executable
STAGES={'old':'5da576b','setter':'31176be','tracking':'d507438','recording':'09ef48f','candidate':'758f0d2'}
COMMON=['GRAMRLevel','GRAMR','GRLevelData','BoundaryConditions','SmallDataIO','PunctureTracker','PETScCommunicator']
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def replace(text,changes):
    for k,v in changes.items():
        pat=r'^'+re.escape(k)+r'\s*=.*$';line=f'{k} = {v}'
        text=re.sub(pat,line,text,flags=re.M) if re.search(pat,text,re.M) else text+'\n'+line+'\n'
    return text
def prepare():
    ROOT.mkdir(exist_ok=True);inputs={};sources=ROOT/'sources';sources.mkdir(exist_ok=True)
    for name,commit in STAGES.items():
        inputs[name]=dict(commit=subprocess.check_output(['git','rev-parse',commit],cwd=REPO,text=True).strip(),overlay='qualified T21 recorder only' if name=='recording' else None)
        src=sources/name
        if src.exists():continue
        src.mkdir();tar=ROOT/(name+'.tar')
        with tar.open('wb') as f:subprocess.run(['git','archive',commit],cwd=REPO,stdout=f,check=True)
        subprocess.run(['tar','-xf',str(tar),'-C',str(src)],check=True)
        if name=='recording':
            subprocess.run(['patch','-p1','-i',str(HERE/'t21-recorder.patch')],cwd=src,check=True)
            shutil.copy2(HERE/'t21-recorder.hpp',src/'Tests/EMSNative/t21-recorder.hpp')
    # Diagnostic-only private instrumentation. No workspace/production source
    # or identity check is edited, and the seed never enters evolved storage.
    src=sources/'poison'
    if not src.exists():
        shutil.copytree(sources/'old',src)
        p=src/'Source/GRChomboCore/GRAMRLevel.cpp';text=p.read_text();needle='        LevelData<FArrayBox> plot_data(levelGrids, num_states, iv_ghosts);'
        assert text.count(needle)==1
        addition='''
        // T23 Tests-only: witness which saved values are never written.
        if (const char *seed = std::getenv("T23_PLOT_BUFFER_SEED"))
        {
            const double value = std::strtod(seed, nullptr);
            for (DataIterator it=plot_data.dataIterator();it.ok();++it)
                plot_data[it()].setVal(value);
        }
'''
        p.write_text('#include <cstdlib>\n'+text.replace(needle,needle+addition))
        import difflib
        (HERE/'t23-plot-buffer-probe.patch').write_text(''.join(difflib.unified_diff(text.splitlines(True),p.read_text().splitlines(True),fromfile='a/Source/GRChomboCore/GRAMRLevel.cpp',tofile='b/Source/GRChomboCore/GRAMRLevel.cpp')))
    for name in (*STAGES,'poison'):
        src=sources/name
        inputs.setdefault(name,dict(commit=inputs['old']['commit'],overlay='private parametric plot-buffer witness'))
        inputs[name]['source_sha256']={str(p.relative_to(src)):sha(p) for d in ('Source','Examples/EMS') for p in (src/d).rglob('*') if p.is_file() and p.suffix in ('.cpp','.hpp','.H','.inc')}
    (HERE/'t23-source-inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    # Build every native unit against its own parameter layout: cross-commit
    # GRAMRLevel.o reuse would violate SimulationParameters/level ABI.
    variants=[('old','old',False),('setter','setter',False),('tracking','tracking',False),('recording','recording',False),('candidate','candidate',True),('candidate-native-flags','candidate',False),('old-all-access','old',True),('poison','poison',False)]
    jobs=[]
    def job(name,directory,command):
        d=ROOT/directory;d.mkdir(parents=True,exist_ok=True);jobs.append(dict(name=name,directory=str(d),command=list(map(str,command))))
    for name,stage,access in variants:job('build-'+name,'builds/'+name,[PY,HERE/'t23-prepare.py','build',name,stage,'all' if access else 'main-only'])
    job('link-old-extra','builds/old-extra',[PY,HERE/'t23-prepare.py','extra','old-extra'])
    seed=Path('/private/tmp/ems-t22/controls/omitted/params.txt').read_text()
    common=dict(num_plot_vars=34,plot_vars='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi Ham Mom1 Mom2 Mom GaussE GaussB',write_plot_ghosts='true',max_steps=1,checkpoint_interval=1,plot_interval=1,ems_data_path='/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet',RH_activate='false',t14_dense_initial_tags='true',t21_rhs_capture='false',t2_guard_initial_data_after_t0='true')
    for name in ('old','setter','tracking','recording','candidate','candidate-native-flags','old-all-access','old-extra'):
        for factory in ('production','dense'):
            run=name+'-'+factory;d=ROOT/'runs'/run
            for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
            p=d/'params.txt';p.write_text(replace(seed,common));exe=ROOT/'builds'/name/(factory+'.ex')
            job(run,'runs/'+run,[exe,p])
    for name,value in [('poison-a','1.25'),('poison-b','2.5')]:
        d=ROOT/'runs'/name
        for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
        p=d/'params.txt';p.write_text(replace(seed,common))
        job(name,'runs/'+name,['/usr/bin/env','T23_PLOT_BUFFER_SEED='+value,ROOT/'builds/poison/production.ex',p])
    # The high hierarchy is initialization-only, never an expensive coarse
    # advance. It is used only if ordinary small-grid diagnostic differences
    # were absent; the small cases always cover the one-step check.
    job('high-fallback','high-fallback',[PY,HERE/'t23-compare.py','fallback'])
    job('qualification','qualification',[PY,HERE/'t23-compare.py','analyze'])
    d=ROOT/'pipeline';d.mkdir(exist_ok=True);(d/'plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    print('T23_PLAN_READY',len(jobs),'serial jobs',flush=True)
def build(name,stage,mode):
    src=ROOT/'sources'/stage;d=ROOT/'builds'/name;d.mkdir(parents=True,exist_ok=True)
    ch=Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib');hdf=Path('/Users/auroradysis/Workspace/EMS-deps/hdf5-1.14')
    base=['/opt/homebrew/bin/g++-16','-O3','-fno-lto','-std=c++17','-fopenmp']
    base+=['-D'+s for s in ('CH_SPACEDIM=2','CH_Darwin','NDEBUG','CH_USE_COMPLEX','CH_USE_64','CH_USE_DOUBLE','CH_USE_HDF5','H5_USE_16_API','CH_USE_LAPACK','CH_FORT_UNDERSCORE','CH_LANG_CC')]
    includes=['-I'+str(src/'Examples/EMS'),'-I'+str(hdf/'include')]+['-I'+str(p) for p in (src/'Source').rglob('*') if p.is_dir()]
    includes+=['-I'+str(ch/'src'/s) for s in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
    libs=['-L'+str(ch)]+['-l'+s+'2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC' for s in ('amrtimedependent','amrtools','boxtools','basetools')]
    libs+=['-L'+str(hdf/'lib'),'-Wl,-rpath,'+str(hdf/'lib'),'-lhdf5','-lz','-L/opt/homebrew/lib/gcc/16','-lgfortran','-lm','-lgomp','-framework','Accelerate']
    flags=base+(['-fno-access-control'] if mode=='all' else [])+includes;commands=[];objects=[]
    # The candidate checkpoint overrides call private GRAMRLevel header
    # methods. Main-only access suppression cannot build that pinned source.
    # Keep all other native units unsuppressed for the narrow flag control.
    access_units=[]
    def call(cmd):
        commands.append(cmd);(d/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');subprocess.run(cmd,check=True)
    for stem in COMMON+['EMSBH2DLevel']+(['EMSBH2DMovingGauge'] if stage=='candidate' else []):
        matches=list(src.rglob(stem+'.cpp'));assert len(matches)==1,(stem,matches)
        extra=['-fno-access-control'] if stage=='candidate' and mode!='all' and stem=='EMSBH2DLevel' else []
        if extra:access_units.append(stem)
        obj=d/(stem+'.o');call(flags+extra+['-c',str(matches[0]),'-o',str(obj)]);objects.append(str(obj))
    for factory,source in [('production',src/'Examples/EMS/Main_EMSBH2DBH.cpp'),('dense',src/'Tests/EMSNative/T14Launch.cpp')]:
        extra=['-fno-access-control'] if factory=='dense' and mode!='all' else []
        call(flags+extra+[str(source),*objects,*libs,'-o',str(d/(factory+'.ex'))])
    (d/'build-spec.json').write_text(json.dumps(dict(stage=stage,mode=mode,additional_access_units=access_units,scope_note='candidate requires EMSBH2DLevel access flag in addition to dense main' if access_units else '',flags=flags,libs=libs,objects=objects,executables={f:sha(d/(f+'.ex')) for f in ('production','dense')},compiler=subprocess.check_output([base[0],'--version'],text=True).splitlines()[0]),indent=2)+'\n')
    print('T23_BUILD_COMPLETE',name,flush=True)
def extra(name):
    original=ROOT/'builds/old';d=ROOT/'builds'/name;d.mkdir(exist_ok=True)
    cmds=json.loads((original/'commands.json').read_text());added=ROOT/'builds/candidate/EMSBH2DMovingGauge.o';assert added.exists()
    new=[]
    for factory,cmd in zip(('production','dense'),cmds[-2:]):
        cmd=cmd.copy();cmd[cmd.index('-o')+1]=str(d/(factory+'.ex'));cmd.insert(cmd.index('-o'),str(added));subprocess.run(cmd,check=True);new.append(cmd)
    (d/'commands.json').write_text(json.dumps(new,indent=2)+'\n');(d/'build-spec.json').write_text(json.dumps(dict(stage='old',mode='old build + actual candidate unused moving object',added_object=str(added),added_sha256=sha(added)),indent=2)+'\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','build','extra']);ap.add_argument('args',nargs='*');a=ap.parse_args();globals()[a.mode](*a.args)
