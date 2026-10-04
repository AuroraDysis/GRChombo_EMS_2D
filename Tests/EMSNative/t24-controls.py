#!/usr/bin/env python3
"""Strict pinned-source safe-abort build, native bits and restart dt regressions."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
ROOT=Path('/private/tmp/ems-t24/controls')
PIN=Path('/private/tmp/ems-t24/source')
PY=sys.executable
os.environ.update(TMPDIR='/private/tmp',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',
    T13_OUTPUT_ROOT='/private/tmp/ems-t24',T13_PROCESS_CAP_BYTES='6000000000',T13_TREE_CAP_BYTES='7500000000',T13_DISK_CAP_BYTES='40000000000')


def module(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def replace(text,changes):
    for k,v in changes.items():
        line=f'{k} = {v}'
        pat=r'^'+re.escape(k)+r'\s*=.*$'
        text=re.sub(pat,lambda _:line,text,flags=re.M) if re.search(pat,text,re.M) else text+'\n'+line+'\n'
    return text


def measured(name,cmd,d,expected=0):
    d.mkdir(parents=True,exist_ok=True)
    with (d/'run.log').open('wb') as log:
        rc=subprocess.run([PY,str(HERE/'t13-run.py'),'--measure',name,'--directory',str(d),'--',*map(str,cmd)],stdout=log,stderr=subprocess.STDOUT).returncode
    (d/'done.exit').write_text(str(rc)+'\n')
    q=json.loads((d/(name+'.resources.json')).read_text())
    assert q['returncode']==q['child_measurement']['returncode']==expected and not q['gate_reason'],(name,q)
    assert q['peak_rss_bytes']<6e9
    return dict(case=name,native_returncode=expected,peak_RSS_bytes=q['peak_rss_bytes'],wall_s=q['wall_seconds'],receipt=str(d/(name+'.resources.json')))


def build():
    ROOT.mkdir(parents=True,exist_ok=True); src=ROOT/'source'
    if src.exists(): shutil.rmtree(src)
    shutil.copytree(PIN,src)
    shutil.copyfile(REPO/'Source/BoxUtils/SafeNanAbort.hpp',src/'Source/BoxUtils/SafeNanAbort.hpp')
    p=src/'Examples/EMS/SimulationParameters.hpp';text=p.read_text()
    old='        pp.load("ems_puncture_tracking_level", ems_puncture_tracking_level, 0);'
    current=(REPO/'Examples/EMS/SimulationParameters.hpp').read_text()
    addition=current[current.index(old)+len(old):current.index('\n\n        pp.load("G_Newton"')]
    assert 'ems_safe_nan_abort' in addition
    text=text.replace(old,old+addition)
    text=text.replace('    // tagging\n','    // tagging\n    bool ems_safe_nan_abort;\n    std::string ems_nan_abort_prefix;\n    double ems_nan_max_abs;\n    int ems_nan_abort_exit_code;\n')
    p.write_text(text)
    p=src/'Examples/EMS/EMSBH2DLevel.cpp';text=p.read_text()
    current=(REPO/'Examples/EMS/EMSBH2DLevel.cpp').read_text()
    start=current.index('    // Check for nan\'s\n');end=current.index('\n}\n',start)
    b=text.index('    // Check for nan\'s\n');e=text.index('\n}\n',b)
    text=text[:b]+current[start:end]+text[e:]
    text=text.replace('#include "NanCheck.hpp"','#include "NanCheck.hpp"\n#include "SafeNanAbort.hpp"')
    p.write_text(text)
    spec=json.loads(Path('/private/tmp/ems-t24/build/build-spec.json').read_text())
    flags=[s.replace(str(PIN),str(src)) for s in spec['flags']]
    flags.insert(1,'-I'+str(REPO/'Source/BoxUtils'))
    libs=spec['libraries'];bd=ROOT/'build';bd.mkdir(exist_ok=True)
    res=[]
    # GRAMRLevel embeds SimulationParameters by value: rebuild every native TU,
    # never mix objects compiled against two different parameter layouts.
    objs=[]
    for name in ('GRAMRLevel','GRAMR','GRLevelData','BoundaryConditions','SmallDataIO','PunctureTracker','PETScCommunicator','EMSBH2DLevel'):
        matches=list(src.rglob(name+'.cpp'));assert len(matches)==1,matches
        obj=bd/(name+'.o');objs.append(obj)
        cmd=flags+['-c',str(matches[0]),'-o',str(obj)]
        res.append(measured('build-'+name,cmd,bd/name))
    res.append(measured('build-main',flags+[str(src/'Examples/EMS/Main_EMSBH2DBH.cpp'),*objs,*libs,'-o',str(bd/'safe.ex')],bd/'main'))
    res.append(measured('build-probe',flags+[str(HERE/'t24-safe-nan-probe.cpp'),*libs,'-o',str(bd/'nan.ex')],bd/'probe'))
    (ROOT/'build-spec.json').write_text(json.dumps(dict(flags=flags,libraries=libs,pin='d507438',source=str(src),resources=res),indent=2)+'\n')
    print('T24_SAFE_BUILD_COMPLETE',flush=True)


def regression():
    cmp=module('t23-compare');res=[];comparisons=[]
    base=Path('/private/tmp/ems-t24/regression/plain-regrid/params.txt').read_text()
    for name,flag in [('abort-omitted',None),('abort-off','false'),('abort-on','true')]:
        d=ROOT/name
        for sub in ('chk','plt'): (d/sub).mkdir(parents=True,exist_ok=True)
        text=base if flag is None else replace(base,dict(ems_safe_nan_abort=flag,ems_nan_abort_prefix=str(d/'abort')))
        (d/'params.txt').write_text(text)
        res.append(measured(name,[ROOT/'build/safe.ex',d/'params.txt'],d))
        assert 'GRChombo finished.' in (d/'run.log').read_text()
        assert not list(d.glob('abort.rank*.json'))
        for step in range(3):
            for prefix in ('chk/EMS_','plt/EMS_Plot_'):
                file=f'{prefix}{step:06d}.2d.hdf5'
                rr=cmp.comparisons(Path('/private/tmp/ems-t24/regression/plain-regrid')/file,d/file,name)
                defined=[r for r in rr if r['field'] in 'chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split() or r['region']=='valid']
                assert not any(r['bit_mismatches'] for r in defined),(name,file)
                comparisons.append(dict(case=name,file=file,defined_values=sum(r['values'] for r in defined),bit_mismatches=0))
    for name,mode,expect in [('finite','finite',0),('nan-worker','nan',86)]:
        d=ROOT/name;d.mkdir(exist_ok=True)
        res.append(measured(name,[ROOT/'build/nan.ex',mode,d/'abort','0'],d,expect))
        b=json.loads((d/'abort.bits.rank0.json').read_text());assert b['bit_identical'] and b['values']==112
        if expect:
            q=json.loads((d/'abort.rank0.json').read_text());assert q['workers_joined'] and q['exit_code']==86 and q['component']=='K'
            assert not list(d.glob('*.tmp.*'))
    checkpoint=Path('/private/tmp/ems-t24/regression/plain-regrid/chk/EMS_000001.2d.hdf5')
    dtrows=[]
    for label,exe,mult,steps in [('old-same',Path('/private/tmp/ems-t24/build/production.ex'),.25,2),('new-omitted',ROOT/'build/safe.ex',.25,2),('new-off',ROOT/'build/safe.ex',.25,2),('new-on',ROOT/'build/safe.ex',.25,2),('half',ROOT/'build/safe.ex',.125,3)]:
        d=ROOT/('restart-'+label)
        for sub in ('chk','plt'): (d/sub).mkdir(parents=True,exist_ok=True)
        text=replace(base,dict(restart_file=checkpoint,ems_data_path='/private/tmp/T24_STATIC_PATH_ABSENT',dt_multiplier=mult,max_steps=steps,verbosity=1,ems_nan_abort_prefix=str(d/'abort')))
        if label!='new-omitted':text=replace(text,dict(ems_safe_nan_abort='true' if label=='new-on' else 'false'))
        (d/'params.txt').write_text(text)
        res.append(measured('restart-'+label,[exe,d/'params.txt'],d))
        import h5py
        for step in range(2,steps+1):
            with h5py.File(d/f'chk/EMS_{step:06d}.2d.hdf5') as f:
                for level in range(int(f.attrs['num_levels'])):
                    g=f[f'level_{level}']; dx=float(g.attrs['dx']);dt=float(g.attrs['dt']);t=float(g.attrs['time'])
                    assert dt==mult*dx and t==.03125+(step-1)*mult*.125
                    dtrows.append(dict(case=label,step=step,level=level,multiplier=mult,dx=dx,dt=dt,time=t))
        if label in ('new-omitted','new-off','new-on'):
            for prefix in ('chk/EMS_','plt/EMS_Plot_'):
                file=f'{prefix}000002.2d.hdf5'
                rr=cmp.comparisons(ROOT/'restart-old-same'/file,d/file,label)
                defined=[r for r in rr if r['field'] in 'chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split() or r['region']=='valid']
                assert not any(r['bit_mismatches'] for r in defined)
                comparisons.append(dict(case='restart-'+label,file=file,defined_values=sum(r['values'] for r in defined),bit_mismatches=0))
    # Half steps reach the same physical endpoint as the nominal restart.
    assert dtrows[-1]['time']==next(r['time'] for r in dtrows if r['case']=='old-same' and r['level']==1)
    for name,rr in [('t24-control-resources.csv',res),('t24-control-bits.csv',comparisons),('t24-restart-dt.csv',dtrows)]:
        with (HERE/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rr[0]);w.writeheader();w.writerows(rr)
    out=dict(status='SERIAL_NATIVE_PASS_MPI_PENDING',resources=res,comparisons=comparisons,restart_dt=dtrows,
        note='restart already recomputes dt: no new restart parameter or evolution change',
        probe='2 OpenMP threads, 112 Float64 read-only values including a worker NaN; closed abort record before exit 86')
    (HERE/'t24-controls-qualification.json').write_text(json.dumps(out,indent=2)+'\n')
    print('T24_NATIVE_CONTROL_REGRESSIONS_PASS',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=('build','regression','all'));a=ap.parse_args()
    ROOT.mkdir(parents=True,exist_ok=True)
    (ROOT/'done.exit').unlink(missing_ok=True)
    try:
        if a.mode in ('build','all'):build()
        if a.mode in ('regression','all'):regression()
    except Exception:
        (ROOT/'done.exit').write_text('1\n');raise
    (ROOT/'done.exit').write_text('0\n')
if __name__=='__main__':main()
