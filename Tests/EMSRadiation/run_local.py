"""Bounded serial builds/evolutions; all regenerable bulk stays outside the tree."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(sys.argv[1]).resolve()
ENV=dict(os.environ,OMP_NUM_THREADS='1',CHOMBO_HOME='/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
OUT.mkdir(parents=True,exist_ok=True)
records=[]


def run(name,cmd,cwd=ROOT,cap=240):
    path=OUT/name; path.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    with (path/'run.log').open('w') as log:
        try:
            code=subprocess.run(list(map(str,cmd)),cwd=cwd,env=ENV,stdout=log,
                                stderr=subprocess.STDOUT,timeout=cap).returncode
        except subprocess.TimeoutExpired:
            code=124
    record=dict(name=name,exit=code,seconds=time.monotonic()-start)
    records.append(record)
    (OUT/'status.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    if code: raise RuntimeError(f'{name} failed: {code}')


def executable(directory,name):
    return next((ROOT/directory).glob(name+'2d.*.ex'))


def params(path,**changes):
    p={}
    for line in (ROOT/'Examples/EMS/params-trumpet-evolve-n32.txt').read_text().splitlines():
        if '=' in line and not line.startswith('#'):
            k,v=line.split('=',1); p[k.strip()]=v.strip()
    p.update(verbosity=0,plot_interval=-1,checkpoint_interval=-1,num_plot_vars=0,
             plot_vars='',N1=64,N2=32,L=16,center='8 0',star_centre='8 0',
             max_box_size=64,stop_time=.25,max_steps=4)
    p.update(changes)
    if not p.get('plot_vars'): p.pop('plot_vars',None)
    path.mkdir(parents=True,exist_ok=True)
    (path/'chk').mkdir(exist_ok=True); (path/'plt').mkdir(exist_ok=True)
    (path/'params.txt').write_text('\n'.join(f'{k} = {v}' for k,v in p.items())+'\n')


def main():
    mode=sys.argv[2]
    if mode=='build':
        for directory in ('Examples/EMS','Tests/EMSRadiation',
                          'Tests/EMSRHFinder','Tests/EMSTrumpet','Tests/EMSCTT'):
            run('build-'+directory.replace('/','-'),['make','-C',directory,'all','DIM=2','-j','2'])
    elif mode=='pulse-grid':
        exe=executable('Tests/EMSRadiation','EMSRadiationGridTest')
        for n in (64,128,256):
            name=f'pulse-grid-{n}'; path=OUT/name
            params(path,N1=n,N2=n//2,L=24,center='12 0',star_centre='12 0',
                   activate_extraction=1,extraction_center='12 0',num_extraction_radii=3,
                   extraction_radii='4 6 8',extraction_levels='0 0 0',
                   num_points_theta=129,num_points_phi=2,max_steps=0,stop_time=0)
            run(name,[exe,'params.txt'],path)
    elif mode in ('evolve','stationary'):
        exe=executable('Examples/EMS','Main_EMSBH2DBH') if mode=='evolve' else \
            executable('Tests/EMSRadiation','EMSStationaryRadiationTest')
        for n,steps,stop in ((128,4,.25),(128,48,3),(256,96,3),(512,192,3)):
            name=f'stationary-{n}-{steps}'; path=OUT/name
            params(path,N1=n,N2=n//2,L=32,center='16 0',star_centre='16 0',
                   ems_radiation_activate=1,activate_extraction=1,
                   extraction_center='16 0',num_extraction_radii=3,
                   extraction_radii='4 6 8',extraction_levels='0 0 0',
                   num_points_theta=129,num_points_phi=2,
                   num_modes=5,modes='2 0 3 0 4 0 5 0 6 0',
                   max_steps=steps,stop_time=stop)
            run(name,[exe,'params.txt'],path)
    elif mode=='off':
        # Build an immutable fork-main snapshot, never compare EMS to itself.
        baseline=OUT/'fork-main'
        baseline.mkdir(exist_ok=True)
        revision='4bf79b8e829159f5edbfb2a769e1fcbe1ad888f7'
        archive=OUT/'fork-main.tar'
        with archive.open('wb') as f:
            subprocess.run(['git','archive',revision],cwd=ROOT,stdout=f,check=True)
        with tarfile.open(archive) as f:
            f.extractall(baseline,filter='data')
        run('build-fork-main',['make','-C',baseline/'Examples/EMS','all','DIM=2','-j','2'])
        (OUT/'baseline.json').write_text(json.dumps(dict(revision=revision))+'\n')
        # HDF5 stores whole-second object timestamps. Run these tiny cases in
        # the same wall-clock second so raw-file identity includes metadata.
        time.sleep(1-time.time()%1)
        for name,directory,base,explicit in (
                ('legacy',baseline/'Examples/EMS','Main_EMSBH2DBH',{}),
                ('default','Examples/EMS','Main_EMSBH2DBH',{}),
                ('explicit-off','Examples/EMS','Main_EMSBH2DBH',
                 {'ems_radiation_activate':0,'ems_rh_expansion_threshold':'1e-7'})):
            path=OUT/name
            params(path,checkpoint_interval=4,**explicit)
            run(name,[executable(directory,base),'params.txt'],path)
        checks={}
        for file in sorted((OUT/'legacy').rglob('*')):
            if file.suffix not in ('.dat','.hdf5'): continue
            rel=file.relative_to(OUT/'legacy')
            expected=file.read_bytes()
            for name in ('default','explicit-off'):
                assert (OUT/name/rel).read_bytes()==expected, (name,rel)
            checks[str(rel)]=hashlib.sha256(expected).hexdigest()
        assert len(checks)==5
        assert not any((OUT/name/'ems_radiation.csv').exists()
                       for name in ('legacy','default','explicit-off'))
        (OUT/'off-identity.json').write_text(json.dumps(checks,indent=2)+'\n')
    elif mode=='amr':
        path=OUT/'amr'
        params(path,N1=64,N2=32,L=32,center='16 0',star_centre='16 0',
               ems_radiation_activate=1,activate_extraction=1,
               extraction_center='16 0',num_extraction_radii=3,
               extraction_radii='4 6 8',extraction_levels='1 1 1',
               num_points_theta=65,num_points_phi=2,max_level=1,regrid_interval=2,
               ems_radiation_wave_level=1,ems_radiation_wave_radius=9,
               max_steps=4,stop_time=.5)
        run('amr',[executable('Examples/EMS','Main_EMSBH2DBH'),'params.txt'],path)
    elif mode=='parameters':
        exe=executable('Examples/EMS','Main_EMSBH2DBH')
        pilot=ROOT/'Examples/EMS/params-radiation-pilot.txt'
        cases=[('valid',[]),('threshold-zero',['ems_rh_expansion_threshold=0']),
               ('threshold-loose',['ems_rh_expansion_threshold=1e-6']),
               ('no-spheres',['activate_extraction=0']),
               ('wrong-axis',['extraction_center=768 1']),
               ('bad-wave-level',['ems_radiation_wave_level=7']),
               ('too-few-angles',['num_points_theta=9'])]
        results=[]
        for name,extra in cases:
            result=subprocess.run([str(exe),str(pilot),'check_params=1',*extra],
                                  env=ENV,capture_output=True,text=True,timeout=10)
            assert (result.returncode==0)==(name=='valid'),name
            results.append(dict(case=name,exit=result.returncode,passed=True))
        p=dict(line.split('=',1) for line in pilot.read_text().splitlines()
               if '=' in line and not line.startswith('#'))
        p={k.strip():v.strip() for k,v in p.items()}
        for key,value in {'RH_num_horizons':'3','RH_level':'3 3 3',
                'RH_start_times':'0 0 70','RH_initial_centre':'752 784 768',
                'RH_initial_radii':'0.63593977642346233 0.63593977642346233 4',
                'RH_num_points':'96 96 96','RH_chase_speeds':'0.125 0.125 0.125',
                'RH_time_step_freq':'400 400 400','ems_rh_expansion_threshold':'1e-12',
                'max_box_size':'16','block_factor':'16',
                'activate_mq_extraction':'1'}.items():
            assert p[key]==value,(key,p[key])
        assert 'mod_F' in p['plot_vars'].split()
        assert len(p['plot_vars'].split())==int(p['num_plot_vars'])
        assert all(p[k].startswith('/HPC_STAGING/') for k in ('ems_data_path','ems_ctt_data_path'))
        (OUT/'parameter-checks.json').write_text(json.dumps(results,indent=2)+'\n')
    elif mode=='fixtures':
        ems=Path('/Users/auroradysis/Workspace/EMS')
        run('trumpet',[executable('Tests/EMSTrumpet','EMSTrumpetFixtureTest'),
                      ROOT/'Tests/EMSTrumpet/fixtures',ROOT/'Tests/EMSTrumpet/fixtures-echo/B',
                      ROOT/'Tests/EMSTrumpet/fixtures-echo/E'])
        run('ctt',[executable('Tests/EMSCTT','EMSCTTFixtureTest'),
                   ems/'artifacts/reference-alpha20-qfile.trumpet',
                   ems/'.data/binary-ctt-t6/n28-r6.ctt',
                   ems/'test/fixtures/emsctt1',OUT/'ctt'])
        # Reuse the existing 26-job driver without changing its archived outputs.
        import importlib.util
        spec=importlib.util.spec_from_file_location('regressions',ROOT/'Tests/EMSRHFinder/regressions.py')
        module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        module.OUT=OUT/'grid'; module.main()
    elif mode=='threshold':
        import importlib.util
        spec=importlib.util.spec_from_file_location('rh',ROOT/'Tests/EMSRHFinder/run_cases.py')
        rh=importlib.util.module_from_spec(spec); spec.loader.exec_module(rh)
        for threshold in ('1e-7','1e-10','1e-12'):
            status=rh.run(OUT/threshold,'static',48,32,
                          ems_rh_expansion_threshold=threshold,RH_chase_speeds='0.125')
            if status: raise RuntimeError('threshold run failed')
    else: raise ValueError(mode)
    (OUT/'complete.json').write_text(json.dumps({'mode':mode,'passed':True})+'\n')


if __name__=='__main__': main()
