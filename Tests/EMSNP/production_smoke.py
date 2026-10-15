"""Build two real production executables and smoke a supplied non-echo hole."""
import hashlib
import csv
import math
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
from prepare_t2 import replace
from config import EMS, LOG_ROOT, CHOMBO_HOME, CXX

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE=Path(os.environ.get('EMS_NP_BASE_ROOT',str(EMS/'.ariadne/worktrees/fork-np-base')))
LOG=LOG_ROOT/'t5'
OUT=HERE/'results-t5'


def call(args,cwd,log):
    with (LOG/log).open('w') as stream:
        subprocess.run([sys.executable,str(HERE/'measure_t3.py'),*map(str,args)],
            cwd=cwd,stdout=stream,stderr=subprocess.STDOUT,check=True)


def prepare():
    profile=Path(os.environ.get('EMS_NP_SMOKE_TRUMPET',str(EMS/'artifacts/reference-alpha20-q0p7.trumpet')))
    metadata={}
    for line in profile.read_text().splitlines():
        if line.startswith('# ') and '=' in line:
            key,value=line[2:].split('=',1);metadata[key]=value
    assert float(metadata['M'])==1 and abs(float(metadata['Q'])-.7)<1e-15
    assert float(metadata['ems_f2'])==-20 and abs(float(metadata['ems_alpha'])-4*__import__('math').pi)<1e-14
    record={'path':str(profile),'bytes':profile.stat().st_size,
        'sha256':hashlib.sha256(profile.read_bytes()).hexdigest(),'header':metadata}
    expected=json.loads((OUT/'reference-input.json').read_text())
    assert record['sha256']==expected['sha256'] and record['bytes']==expected['bytes'],profile
    (OUT/'reference-input.json').write_text(json.dumps(record,indent=2)+'\n')
    source=(HERE/'results-t3/rn-n256/params.txt').read_text()
    # Three short coarse steps of the non-echo reference; KO stays enabled.
    common={'N1':64,'N2':32,'L':16,'center':'8 0','star_centre':'8 0',
        'ems_data_path':profile,'ems_f2':-20,'bh_charge':.7,
        'max_steps':3,'stop_time':.2,'dt_multiplier':.1,'sigma':1,
        'checkpoint_interval':3,'plot_interval':-1,'calculate_constraint_violations':'true',
        'num_extraction_radii':2,'extraction_radii':'2 4','extraction_levels':'0 0',
        'num_points_theta':33,'ems_np_activate':'false','ems_weyl_matter':'false'}
    for name in ('baseline-off','current-off','current-on'):
        folder=OUT/name;folder.mkdir(exist_ok=True)
        (folder/'chk').mkdir(exist_ok=True)
        values=dict(common)
        if name=='current-on':values.update(ems_np_activate='true',ems_weyl_matter='true')
        (folder/'params.txt').write_text(replace(source,values))


def authenticate_baseline():
    entries=subprocess.check_output(['git','ls-tree','-rz','00d7171'],cwd=ROOT).split(b'\0')
    count=0
    for entry in filter(None,entries):
        meta,name=entry.split(b'\t',1);mode,kind,expected=meta.decode().split()
        if kind!='blob':continue
        relative=name.decode();path=BASE/relative
        data=os.readlink(path).encode() if mode=='120000' else path.read_bytes()
        if relative=='Source/GRChomboCore/GRAMRLevel.hpp':
            original=subprocess.check_output(['git','show','00d7171:'+relative],cwd=ROOT).decode()
            for declaration in ('virtual void writeCheckpointHeader(HDF5Handle &a_handle) const;',
                                'virtual void readCheckpointHeader(HDF5Handle &a_handle);'):
                original=original.replace('    '+declaration,'  protected:\n    '+declaration+'\n  private:')
            assert data==original.encode(),relative
        else:
            assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==expected,relative
        count+=1
    record={'commit':subprocess.check_output(['git','rev-parse','00d7171'],cwd=ROOT,text=True).strip(),
        'verified_blobs':count,'exception':'two checkpoint header access declarations only'}
    (OUT/'baseline-source.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def main():
    LOG.mkdir(parents=True,exist_ok=True);OUT.mkdir(exist_ok=True)
    if '--prepare' in sys.argv:prepare();return
    if '--verify' in sys.argv:verify();return
    if not BASE.exists():
        BASE.mkdir(parents=True)
        archive=subprocess.Popen(['git','archive','00d7171'],cwd=ROOT,stdout=subprocess.PIPE)
        subprocess.run(['tar','-x','-C',str(BASE)],stdin=archive.stdout,check=True)
        assert archive.wait()==0
    header=BASE/'Source/GRChomboCore/GRAMRLevel.hpp'
    before=header.read_text()
    for declaration in ('virtual void writeCheckpointHeader(HDF5Handle &a_handle) const;',
                        'virtual void readCheckpointHeader(HDF5Handle &a_handle);'):
        assert before.count(declaration)==1
        if '  protected:\n    '+declaration not in before:
            before=before.replace('    '+declaration,'  protected:\n    '+declaration+'\n  private:')
    if header.read_text()!=before:header.write_text(before)
    source_record=authenticate_baseline()
    env=os.environ
    env.setdefault('CHOMBO_HOME',CHOMBO_HOME)
    env['OMP_NUM_THREADS']='1'
    for name,tree in [('baseline',BASE),('current',ROOT)]:
        call(['make','-C',tree/'Examples/EMS','all','DIM=2','CXX='+CXX,'-j','2'],ROOT,name+'-production-build.log')
        exes=list((tree/'Examples/EMS').glob('Main_EMSBH2DBH*.ex'))
        assert len(exes)==1,exes
        for case in (['baseline-off'] if name=='baseline' else ['current-off','current-on']):
            call([exes[0],'params.txt'],OUT/case,case+'.log')
    (OUT/'production-runs.json').write_text(json.dumps({'baseline_commit':'00d7171',
        'baseline_root':str(BASE),'baseline_portability_change':'only the two checkpoint headers protected',
        'baseline_checkout':source_record,
        'valid_smoke_profile':json.loads((OUT/'reference-input.json').read_text())['path'],
        'runs':{'baseline-off':0,'current-off':0,'current-on':0}},indent=2)+'\n')


def verify():
    # 1e-12 covers Float64 arithmetic/derivative conditioning over three RK4 steps.
    # It is numerical agreement, not byte identity of HDF5 containers.
    tolerance=1e-12
    base=OUT/'baseline-off';current=OUT/'current-off';enabled=OUT/'current-on'
    result={'absolute_tolerance':tolerance,'basis':'Float64 conditioning at dx=.25, three RK4 steps',
        'checkpoints':[],'text_outputs':[]}
    for step in (0,3):
        name=f'EMS_{step:06d}.2d.hdf5'
        checkpoint=lambda folder:folder/'chk'/name if (folder/'chk'/name).exists() else LOG/'smoke-checkpoints'/folder.name/name
        command=[os.environ.get('H5DIFF','h5diff'),'-d',str(tolerance),checkpoint(base),checkpoint(current)]
        with (LOG/f'h5diff-{step}.log').open('w') as stream:
            status=subprocess.run(list(map(str,command)),stdout=stream,stderr=subprocess.STDOUT).returncode
        assert status==0,name
        result['checkpoints'].append({'step':step,'exit':status})
    for p in sorted(base.iterdir()):
        if p.suffix not in ('.dat','.csv') or p.name=='params.txt':continue
        other=current/p.name
        assert other.exists(),p
        if p.suffix=='.csv':
            a=list(csv.reader(p.open()));b=list(csv.reader(other.open()));assert a.pop(0)==b.pop(0)
        else:
            read=lambda f:[line.split() for line in f.read_text().splitlines() if line.strip() and not line.startswith('#')]
            a=read(p);b=read(other)
        assert len(a)==len(b),p
        maximum=0
        for x,y in zip(a,b):
            assert len(x)==len(y)
            for s,t in zip(x,y):
                x0=float(s);x1=float(t)
                assert math.isfinite(x0) and math.isfinite(x1)
                maximum=max(maximum,abs(x0-x1))
        assert maximum<=tolerance,p
        result['text_outputs'].append({'file':p.name,'rows':len(a),'max_absolute_difference':maximum})
    assert (base/'ems-gauge.txt').read_text()==(current/'ems-gauge.txt').read_text()
    assert not (base/'ems_np.csv').exists() and not (current/'ems_np.csv').exists()
    with (enabled/'ems_np.csv').open() as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==56 and len({float(d['t']) for d in rows})==4
    assert all(math.isfinite(float(v)) for d in rows for v in d.values())
    assert abs(max(float(d['t']) for d in rows)-.075)<=64*sys.float_info.epsilon*.075
    result.update(status='PASS',finite_NP_rows=len(rows),times=4,coarse_steps=3,
        max_time=max(float(d['t']) for d in rows),default_off_new_files=False)
    (OUT/'production-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    artifacts=[]
    for folder in (base,current,enabled):
        target=LOG/'smoke-checkpoints'/folder.name;target.mkdir(parents=True,exist_ok=True)
        for p in (folder/'chk').glob('*.hdf5'):
            destination=target/p.name
            shutil.copy2(p,destination)
            assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256(destination.read_bytes()).digest()
            p.unlink()
        for p in target.glob('*.hdf5'):
            artifacts.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
    (OUT/'smoke-checkpoints.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
