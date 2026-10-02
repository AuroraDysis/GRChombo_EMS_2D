#!/usr/bin/env python3
"""Build and register short current-state captures; reuse T18 resource gates."""
import csv,hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t18')
PY='/Users/auroradysis/miniconda3/bin/python'
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
w=module('t18work',HERE/'t18-work.py')
replace=w.module('params',HERE/'t13-prepare.py').replace
def build():
    original=json.loads((ROOT/'native-window.ex-build-command.json').read_text())
    for source,binary in [('t18-launch.cpp','launch.ex'),('t16-replay.cpp','launch-replay.ex')]:
        command=[]
        for arg in original:
            if source=='t16-replay.cpp' and arg.endswith('.o'):continue
            if arg.endswith('/t18-native.cpp'):arg=str(HERE/source)
            elif '/Examples/EMS/o/' in arg:arg=str(ROOT/'production-o'/arg.split('/Examples/EMS/o/')[1])
            elif arg==str(ROOT/'native-window.ex'):arg=str(ROOT/binary)
            command.append(arg)
        (ROOT/(binary+'-build-command.json')).write_text(json.dumps(command,indent=2)+'\n')
        w.measured('launch-build-'+binary,command)
def prepare():
    census=list(csv.DictReader((HERE/'t18-census-boxes.csv').open()))
    base=(w.T17/'t17-coarse-geometric.txt').read_text()
    jobs=[]
    for top,dt in [(12,.25),(12,.375),(12,.5),(13,.25)]:
        name=('coarse' if top==12 else 'fine')+'-dt'+str(dt)
        directory=ROOT/'launch'/name;directory.mkdir(parents=True,exist_ok=True)
        boxes=HERE/('t18-launch-'+name+'-boxes.txt')
        lines=[]
        for r in census:
            l=int(r['level'])
            if int(r['max_box_size'])!=128 or not 1<=l<=top:continue
            shift=1024*2**l
            coords=[int(r[k]) for k in ('x0','y0','x1','y1')];coords[0]-=shift;coords[2]-=shift
            assert 0<=coords[0]<=coords[2]<512*2**l and 0<=coords[1]<=coords[3]<256*2**l
            lines.append(' '.join(map(str,[l,*coords])))
        boxes.write_text('\n'.join(lines)+'\n')
        radii=' '.join(format(x,'.17g') for x in ([224/1.2,144/1.2]+[224/2**i/1.2 for i in range(3,top+1)]))
        param=HERE/('t18-launch-'+name+'.txt')
        param.write_text(replace(base,dict(N1=512,N2=256,L=896,center='448 0',star_centre='448 0',
            mass_extraction_center='448 0',max_level=top,max_steps=1,stop_time=2,
            regrid_interval=' '.join(['0']*(top+1)),checkpoint_interval=-1,plot_interval=-1,
            ems_track_punctures='false',ems_binary_refinement='false',dt_multiplier=dt,
            t18_monitor='false',t18_capture_stop=.001922607421875,t18_capture_boxes=str(boxes),
            num_mass_extraction_radii=top,mass_extraction_levels=' '.join(map(str,range(1,top+1))),
            mass_extraction_radii=radii,t13_launch_stop_time=0)))
        jobs.append(dict(name=name,directory=str(directory),command=['env','OMP_NUM_THREADS=2',str(ROOT/'launch.ex'),str(param)]))
    plan=ROOT/'launch/plan.json';plan.write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    record=dict(purpose='Complete constraint and temporal-versus-spatial launch evidence omitted from the window histories',
        binary=True,h0_M=1.75,centre_M=448,translation_M=1792,base_shape=[512,256],
        refined_boxes='exact T17 census, integer translation on every refined level',
        unchanged='native physics, Float64, KO=1, geometric lapse, point transfers, initial-data guard',
        endpoint_M=.001922607421875,common_positive_clocks_M=[.000640869140625,.00128173828125,.001922607421875],
        window_W_M=[.00075,.0025],rays=['axis','diagonal'],holes=['left','right'],norms=['peak','RMS'],
        temporal_measure='norm of the difference between native t0-subtracted profiles at coarse dt and coarse dt0.25',
        spatial_measure='norm of the difference between native t0-subtracted profiles at coarse and fine dt0.25',
        subdominant_rule='temporal error plus 5x P8/P10 spread of the matched temporal difference < spatial difference minus 5x P8/P10 spread of the matched spatial difference; roundoff/unresolved cases labelled',
        production_control='Compare saved profiles with full T17 geometric captures at exactly the same clocks',
        process_RSS_cap_bytes=4000000000,threads=2,jobs=jobs)
    (HERE/'t18-launch-registration.json').write_text(json.dumps(record,indent=2)+'\n')
    print(plan)
if __name__=='__main__':
    if sys.argv[1]=='build':build()
    else:prepare()
