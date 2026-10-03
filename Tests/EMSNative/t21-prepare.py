#!/usr/bin/env python3
"""Prepare a recording-only overlay and bounded local qualification; no level-15 run."""
import argparse, csv, difflib, hashlib, json, os, resource, shutil, subprocess, sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];ROOT=Path('/private/tmp/ems-t21')
PY=sys.executable
FIELDS='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split()
CONTRACT=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0023/submissions/exp-0023')

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def replace(text,changes):
    import re
    for k,v in changes.items():
        line=f'{k} = {v}'
        if re.search(r'^'+re.escape(k)+r'\s*=',text,re.M):text=re.sub(r'^'+re.escape(k)+r'\s*=.*$',line,text,flags=re.M)
        else:text+='\n'+line+'\n'
    return text

def overlay():
    path=REPO/'Examples/EMS/EMSBH2DLevel.hpp';old=path.read_text()
    new=old.replace('#include "EMSCouplingFunction.hpp"','#include "EMSCouplingFunction.hpp"\n#include "../../Tests/EMSNative/t21-recorder.hpp"')
    new=new.replace('    bool m_t6_interface_diagnostics = [] {','    T21RHSRecorder m_t21;\n    bool m_t6_interface_diagnostics = [] {')
    assert old!=new
    diffs=list(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/Examples/EMS/EMSBH2DLevel.hpp',tofile='b/Examples/EMS/EMSBH2DLevel.hpp'))
    path=REPO/'Examples/EMS/EMSBH2DLevel.cpp';old=path.read_text()
    needle='    if (t13 || (m_t7.enabled && m_level>=4 && m_level<=6))\n'
    assert old.count(needle)==1
    insertion='''    const bool t21=m_t21.selected(m_level,m_time,m_dt);
    if (t21 && (t13 || m_t7.enabled))
        MayDay::Error("T21 production recorder cannot be combined with T7/T13 capture");
    if (t21)
    {
        m_t21.next_call();
        m_t21.set_kernel_type(typeid(my_ccz4_cartoon).name());
        for (DataIterator it=a_soln.dataIterator();it.ok();++it)
        {
            const Box valid=a_soln.disjointBoxLayout()[it()];
            FArrayBox capture(a_rhs[it()].box(),2*NUM_VARS);
            std::vector<Box> regions{valid};
            my_ccz4_cartoon.t7_capture(&capture,&regions);
            BoxLoops::loop(make_compute_pack(my_ccz4_cartoon,set_analysis_vars_zero),a_soln[it()],a_rhs[it()],valid);
            const auto &record_state=m_t21.checkpoint_source?(*m_t21.checkpoint_source)[it()]:a_soln[it()];
            m_t21.write(m_level,m_rk_stage,a_time,m_time,m_dx,m_dt,
                        {m_p.center[0],m_p.center[1]},valid,record_state,capture,a_rhs[it()]);
        }
    }
    else if (t13 || (m_t7.enabled && m_level>=4 && m_level<=6))
'''
    new=old.replace(needle,insertion).replace('#include <cstdint>','#include <cstdint>\n#include <typeinfo>')
    before='    if (m_t7.enabled) t7_record(2,a_soln);'
    assert new.count(before)==1
    new=new.replace(before,'''    if (m_t21.selected(m_level,m_time,m_dt))
        m_t21.before_projection(m_level,m_rk_stage,a_time,m_dx,
                               {m_p.center[0],m_p.center[1]},a_soln,m_p.min_chi,m_p.min_lapse);
'''+before)
    diffs+=list(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/Examples/EMS/EMSBH2DLevel.cpp',tofile='b/Examples/EMS/EMSBH2DLevel.cpp'))
    (HERE/'t21-recorder.patch').write_text(''.join(diffs))

def register():
    ROOT.mkdir(exist_ok=True);overlay()
    cfg=dict(status='recording_only; old-gauge level15 plan superseded by operator/consult9b',
        dt_multiplier=.25,ray_radius_min_M=1e-5,ray_radius_max_M=.02,ray_points=2048,
        rays=[dict(name='axis',direction=[1.,0.]),dict(name='equator',direction=[0.,1.]),dict(name='diagonal',direction=[2**-.5,2**-.5])],
        fixed_radii_M=[.00005,.0001,.0002,.0005,.001,.002,.00341796875,.004833737762017416,.0055,.0060286129392158804,.0075,.01,.013671875,.02],
        fields=FIELDS,levels=[12,13,14],stage_steps_per_level=4,puncture_cells=4,puncture_radius_M=.0004,
        strip_half_cells=12,diagnostic_clocks_M=[.875*k for k in range(13)],
        interpolation='I8 central; abs(I8-I10) spread; 128*Float64 epsilon*max(1,abs(value)) floor; fivefold significance',
        B_comparison='B1/B2 are gauge-labelled driver/offset variables and cannot be compared across gauge packages',
        future_gauge_control='not implemented, run or admitted by T21')
    (HERE/'t21-recording.json').write_text(json.dumps(cfg,indent=2)+'\n')
    params='''# Opt-in recording only. Merge with the selected run's frozen parameters.
t21_rhs_capture = true
t21_gauge_label = experimental
t21_driver_semantics = integrated_offset_decaying_0p1
t21_rhs_prefix = stage/native
t21_initial_rk_steps = 4
t21_puncture_cells = 4
t21_puncture_radius = 0.0004
t21_profile_radius = 0.02
t21_strip_half_cells = 12
t21_levels = 12 13 14
t21_ray_directions = 1 0 0 1 0.7071067811865475244 0.7071067811865475244
t21_fields = '''+' '.join(FIELDS)+'\n'
    (HERE/'t21-recording-params.txt').write_text(params)
    records=[]
    import numpy as np
    data=list(csv.DictReader((HERE/'t20-fixed-radius.csv').open()))
    for field in ('A11','K','lapse'):
        for ray in ('axis','equator','diagonal'):
            for radius in (.0055,.00608,.0075):
                a=sorted([r for r in data if r['field']==field and r['ray']==ray and float(r['target_radius_M'])==radius and float(r['time_M'])>0],key=lambda r:float(r['time_M']))
                assert len(a)==12
                d1=np.array([float(r['low_mid']) for r in a]);d2=np.array([float(r['mid_high']) for r in a])
                correlation=float(d1@d2/(np.linalg.norm(d1)*np.linalg.norm(d2)))
                records.append(dict(field=field,ray=ray,target_radius_M=radius,actual_radius_M=a[0]['radius_M'],positive_clocks=12,
                    signed_history_inner_product=correlation,
                    descriptive_norm_ratio_order=float(np.log(np.linalg.norm(d1)/np.linalg.norm(d2))/np.log(1.5)),
                    interpretation='signed differences anti-correlate; norm ratio is not a common-leading-error witness' if correlation<0 else 'signed correlation reported; norm ratio alone does not establish a common leading error'))
    with (HERE/'t21-signed-history-check.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    inputs=[REPO/'Examples/EMS/EMSBH2DLevel.hpp',REPO/'Examples/EMS/EMSBH2DLevel.cpp',REPO/'Source/Cartoon/CCZ4Cartoon.impl.hpp',REPO/'Source/CCZ4/ExperimentalGauge.hpp',REPO/'Source/CCZ4/MovingPunctureGauge.hpp',
        CONTRACT/'production/params-E-high.txt',CONTRACT/'submit-contract.md',HERE/'t20-fixed-radius.csv',
        Path('/Users/auroradysis/Workspace/EMS-deps/worktrees/consult-9-reply.md'),Path('/Users/auroradysis/Workspace/EMS-deps/worktrees/consult-9b-reply.md')]
    inputs.extend([REPO/'Tests/EMSNative/T14DenseTags.hpp',REPO/'Tests/EMSNative/T14Launch.cpp',Path('/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet')])
    receipt=dict(status='SUPERSEDED_OLD_GAUGE_LEVEL15_NO_CENSUS_OR_P_PARAMS_COMPLETED',
        fork_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        overlay_targets_uncommitted_preexisting=True,inputs=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in inputs],
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1)
    (HERE/'t21-registration.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

def snapshot():
    # This local build overlay is not a fork commit and is not a submission authorization.
    src=ROOT/'src'
    if src.exists():return
    src.mkdir()
    archive=ROOT/'HEAD.tar'
    with archive.open('wb') as f:subprocess.run(['git','archive','HEAD'],cwd=REPO,stdout=f,check=True)
    subprocess.run(['tar','-xf',str(archive),'-C',str(src)],check=True)
    for p in ('Examples/EMS/EMSBH2DLevel.cpp','Examples/EMS/Main_EMSBH2DBH.cpp','Examples/EMS/SimulationParameters.hpp','Source/BlackHoles/PunctureTracker.cpp'):
        shutil.copy2(REPO/p,src/p)
    # Source-only dependencies of the pre-existing T17 diagnostics, without
    # importing compiled objects, outputs or data into the snapshot.
    for directory in ('Source','Examples/EMS'):
        for p in (REPO/directory).rglob('*'):
            if p.is_file() and p.suffix in ('.hpp','.cpp','.H','.inc'):
                target=src/p.relative_to(REPO);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    for p in HERE.glob('t21-*'):
        if p.is_file():shutil.copy2(p,src/'Tests/EMSNative'/p.name)
    # Existing frozen objects are reused except the changed level TU and Main TU.
    manifest={}
    obj=REPO/'Examples/EMS/o/2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
    for p in obj.glob('*.o'):manifest[str(p)]=sha(p)
    (ROOT/'reused-objects.json').write_text(json.dumps(manifest,indent=2)+'\n')

def build():
    snapshot();src=ROOT/'src';ch=Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib');hdf=Path('/Users/auroradysis/Workspace/EMS-deps/hdf5-1.14')
    base=['/opt/homebrew/bin/g++-16','-O3','-std=c++17','-fopenmp','-fno-access-control']
    base+=['-D'+s for s in ('CH_SPACEDIM=2','CH_Darwin','NDEBUG','CH_USE_COMPLEX','CH_USE_64','CH_USE_DOUBLE','CH_USE_HDF5','H5_USE_16_API','CH_USE_LAPACK','CH_FORT_UNDERSCORE','CH_LANG_CC')]
    base+=['-I'+str(src/'Examples/EMS'),'-I'+str(hdf/'include')]
    base+=['-I'+str(d) for d in (src/'Source').rglob('*') if d.is_dir() and not d.name.startswith('.')]
    base+=['-I'+str(ch/'src'/d) for d in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
    obj=REPO/'Examples/EMS/o/2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
    objects=[str(p) for p in sorted(obj.glob('*.o')) if p.stem not in ('Main_EMSBH2DBH','EMSBH2DLevel')]
    libs=['-L'+str(ch)]+['-l'+s+'2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC' for s in ('amrtimedependent','amrtools','boxtools','basetools')]
    libs+=['-L'+str(hdf/'lib'),'-Wl,-rpath,'+str(hdf/'lib'),'-lhdf5','-lz','-L/opt/homebrew/lib/gcc/16','-lgfortran','-lm','-lgomp','-framework','Accelerate']
    commands=[]
    def call(cmd):
        commands.append(cmd);(ROOT/'build-commands.json').write_text(json.dumps(commands,indent=2)+'\n');subprocess.run(cmd,check=True)
    call(base+['-c',str(src/'Examples/EMS/EMSBH2DLevel.cpp'),'-o',str(ROOT/'level-original.o')])
    call(base+[str(src/'Examples/EMS/Main_EMSBH2DBH.cpp'),str(ROOT/'level-original.o'),*objects,*libs,'-o',str(ROOT/'original.ex')])
    subprocess.run(['patch','-p1','--batch','--forward','-i',str(HERE/'t21-recorder.patch')],cwd=src,check=True)
    call(base+['-c',str(src/'Examples/EMS/EMSBH2DLevel.cpp'),'-o',str(ROOT/'level-recording.o')])
    call(base+[str(src/'Examples/EMS/Main_EMSBH2DBH.cpp'),str(ROOT/'level-recording.o'),*objects,*libs,'-o',str(ROOT/'recording.ex')])
    call(base+[str(src/'Tests/EMSNative/t21-native-replay.cpp'),str(ROOT/'level-recording.o'),*objects,*libs,'-o',str(ROOT/'replay.ex')])
    print('T21_RECORDING_BUILD_PASS',flush=True)

def plans():
    jobs=[dict(name='build',directory=str(ROOT/'build'),command=[PY,str(HERE/'t21-prepare.py'),'build'])]
    seed=(Path('/private/tmp/ems-t13/controls/E-point-off/params.txt')).read_text()
    changes=dict(ems_use_maximal_initial_lapse='true',max_steps=1,checkpoint_interval=1,plot_interval=1,
        t13_launch_stop_time=0,t6_interface_diagnostics='false',t7_diagnostics='false',t2_guard_initial_data_after_t0='true',amr_transfer='point',verbosity=0)
    recording=(HERE/'t21-recording-params.txt').read_text()
    recording=replace(recording,dict(t21_levels='0 1',t21_puncture_radius=.2,t21_profile_radius=1,t21_strip_half_cells=12))
    for mode in ('original','omitted','off','on'):
        d=ROOT/'controls'/mode
        for s in ('plt','chk'):(d/s).mkdir(parents=True,exist_ok=True)
        text=replace(seed,changes)
        if mode=='off':text+='\nt21_rhs_capture = false\n'
        if mode=='on':text+='\n'+replace(recording,dict(t21_rhs_prefix=str(d/'native')))
        param=d/'params.txt';param.write_text(text)
        exe=ROOT/('original.ex' if mode=='original' else 'recording.ex')
        jobs.append(dict(name=mode,directory=str(d),command=[str(exe),str(param)]))
    for step in (0,1):
        d=ROOT/'replay'/f'step{step:06}';d.mkdir(parents=True,exist_ok=True)
        text=replace(seed,changes)+'\n'+replace(recording,dict(t21_rhs_prefix=str(d/'native')))
        text=replace(text,dict(restart_file=ROOT/f'controls/original/chk/EMS_{step:06}.2d.hdf5',ems_data_path='/T21_STATIC_FILE_MUST_NOT_BE_READ'))
        param=d/'params.txt';param.write_text(text)
        jobs.append(dict(name=f'replay-{step}',directory=str(d),command=[str(ROOT/'replay.ex'),str(param)]))
    jobs.append(dict(name='qualification',directory=str(ROOT/'qualification'),command=[PY,str(HERE/'t21-qualify.py')]))
    jobs.append(dict(name='seal',directory=str(ROOT/'seal'),command=[PY,str(HERE/'t21-manifest.py')]))
    pipe=ROOT/'pipeline';pipe.mkdir(exist_ok=True)
    (pipe/'plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n');print(pipe/'plan.json')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['register','build','plans']);a=ap.parse_args();globals()[a.mode]()
