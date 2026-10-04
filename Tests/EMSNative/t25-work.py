#!/usr/bin/env python3
"""Local T25 initial-data regeneration and offline native maps; no evolution.

Maps consume only checkpoint geometry plus current numerical centres/coupling
parameters. Initializers alone may read trumpet/CTT, at t=0, with the guard on.
Individual operations have wait4/RSS receipts and an atomic done marker.
Native maps are capped at 3GB; t=0 initialization at 6GB. This is NOT a finder
sweep. Surface parameters are public CSV inputs; change them on the CLI.
"""
import argparse,csv,hashlib,json,os,resource,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t25')
PY='/Users/auroradysis/miniconda3/bin/python'
def replace(text,changes):
    seen=set();out=[]
    for line in text.splitlines():
        key=line.split('=',1)[0].strip() if '=' in line and not line.lstrip().startswith('#') else None
        if key in changes:out.append(f'{key} = {changes[key]}');seen.add(key)
        else:out.append(line)
    out.extend(f'{k} = {v}' for k,v in changes.items() if k not in seen)
    return '\n'.join(out)+'\n'
def measure(name,directory,command,cap):
    directory.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,T13_OUTPUT_ROOT=str(ROOT),T13_PROCESS_CAP_BYTES=str(cap),T13_TREE_CAP_BYTES='7000000000',T13_DISK_CAP_BYTES='12000000000')
    with (directory/'run.log').open('wb') as log:
        rc=subprocess.run([PY,str(HERE/'t13-run.py'),'--measure',name,'--directory',str(directory),'--',*map(str,command)],stdout=log,stderr=subprocess.STDOUT,env=env).returncode
    (directory/'done.exit.tmp').write_text(str(rc)+'\n');(directory/'done.exit.tmp').replace(directory/'done.exit')
    return rc
def surfaces(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=['id','family','centre','a','c'],lineterminator='\n');w.writeheader();w.writerows(rows)
def prepare():
    ROOT.mkdir(parents=True,exist_ok=True)
    for kind,src in [('single',Path('/private/tmp/ems-t17/finder-diagnosis/initialization/coarse-unboosted-native/params.txt')),
                     ('binary',Path('/private/tmp/ems-t17/census/coarse-native/params.txt'))]:
        d=ROOT/'initialization'/kind;d.mkdir(parents=True,exist_ok=True)
        (d/'params.txt').write_text(replace(src.read_text(),dict(max_steps=0,stop_time=0,plot_interval=-1,checkpoint_interval=0,
            RH_activate='false',RH_num_horizons=0,ems_track_punctures='false',t2_guard_initial_data_after_t0='true')))
    radii=[.0013,.002,.003,.004,.005,.0055,.0058,.0059,.006,.00602,.00604,.00606,.00608,.0061,.0062,.0064,.0068,.0075,.009,.012,.02,.04,.08,.1]
    rows=[dict(id=f'single-{j:03}',family='single-sphere',centre=40,a=r,c=r) for j,r in enumerate(radii)]
    surfaces(HERE/'t25-surfaces-single.csv',rows)
    rows=[]
    for hole,centre in enumerate((2224.,2256.)):
        for offset in [-.001,-.0005,0,.0005,.001]:
            for aspect in [.85,.95,1.,1.05,1.15]:
                for j,r in enumerate(radii):
                    fam=f'initial-h{hole}-o{offset:+.4f}-c{aspect:.2f}'
                    rows.append(dict(id=f'{fam}-r{j:03}',family=fam,centre=centre+offset,a=r,c=r*aspect))
    surfaces(HERE/'t25-surfaces-binary.csv',rows)
    track=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall/ev1/runs/exp-0024/merger/punctures.dat')
    trow=next(line.split() for line in track.read_text().splitlines() if not line.startswith('#') and line.split() and float(line.split()[0])==133.)
    # Cartoon tracker prints t,x1,y1,x2,y2. Fail rather than silently change schema.
    assert len(trow)==5,trow
    centres=[float(trow[1]),float(trow[3])];rows=[]
    for hole,centre in enumerate(centres):
        for aspect in [.75,1.,1.5,2.]:
            for j,r in enumerate([3*1.75/4096,.002,.003,.004,.006,.008,.01,.015,.02,.03,.04,.047,.058,.07,.08,.09,.1]):
                fam=f'late-h{hole}-c{aspect:.2f}'
                rows.append(dict(id=f'{fam}-r{j:03}',family=fam,centre=centre,a=r,c=r*aspect))
    mid=sum(centres)/2
    for aspect in [1.,1.5,2.,3.]:
        for j,r in enumerate([.25,.3,.4,.5,.6,.8,1.,1.5,2.,3.,4.,6.,8.]):
            fam=f'late-common-c{aspect:.2f}'
            rows.append(dict(id=f'{fam}-r{j:03}',family=fam,centre=mid,a=r,c=r*aspect))
    surfaces(HERE/'t25-surfaces-late.csv',rows)
    (HERE/'t25-centres.json').write_text(json.dumps(dict(single=[40.],binary=[2224.,2256.],late=centres,late_time=133.,track=str(track)),indent=2)+'\n')
def initialize(kind):
    d=ROOT/'initialization'/kind
    (d/'chk').mkdir(exist_ok=True);(d/'plt').mkdir(exist_ok=True)
    rc=measure('initialization',d,['/private/tmp/ems-t17/production-t17.ex','params.txt'],6000000000)
    assert rc==0 and 'GRChombo finished.' in (d/'run.log').read_text()
    assert (d/'chk/EMS_000000.2d.hdf5').is_file()
def mapping(kind,points):
    d=ROOT/'maps'/f'{kind}-N{points}';d.mkdir(parents=True,exist_ok=True)
    if kind=='late':
        base=Path('/private/tmp/ems-t24/horizon/params.txt').read_text()
        cp='/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall/chk/EMS_000152.2d.hdf5'
    else:
        base=(ROOT/'initialization'/kind/'params.txt').read_text()
        cp=str(ROOT/'initialization'/kind/'chk/EMS_000000.2d.hdf5')
    # Deliberately absent initial inputs: proves maps cannot read static/CTT files.
    centres=json.loads((HERE/'t25-centres.json').read_text())[kind]
    changes=dict(restart_file=cp,hdf5_subpath='""',ems_data_path='/T25-FORBIDDEN-static',ems_ctt_data_path='""',
        max_steps=0,stop_time=0,RH_activate='false',RH_num_horizons=0,ems_track_punctures='false',
        map_points=points,map_surfaces=HERE/f't25-surfaces-{kind}.csv',map_output='map.csv',map_angles='angles.csv',
        map_punctures=' '.join(map(str,centres)))
    (d/'params.txt').write_text(replace(base,changes))
    rc=measure('native-map',d,[ROOT/'build/t25-expansion-map.ex','params.txt'],3000000000)
    assert rc==0 and 'T25_EXPANSION_MAP_COMPLETE; no advances' in (d/'run.log').read_text()
    rows=list(csv.DictReader((d/'map.csv').open()));expected=list(csv.DictReader((HERE/f't25-surfaces-{kind}.csv').open()))
    assert len(rows)==len(expected) and len({r['id'] for r in rows})==len(rows)
    (HERE/f't25-map-{kind}-N{points}.csv').write_text((d/'map.csv').read_text())
def seed(kind,family):
    import importlib.util
    spec=importlib.util.spec_from_file_location('brackets',HERE/'t25-brackets.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
    rs=list(csv.DictReader((HERE/f't25-map-{kind}-N192.csv').open()))
    bs=b.pointwise_barriers(rs)
    if family:bs=[r for r in bs if r['family']==family]
    else:bs=[r for r in bs if r['status']=='POINTWISE_NEGATIVE_TO_POSITIVE']
    assert len(bs)==1, 'select exactly one measured barrier family'
    r=bs[0];tag=family or kind;d=ROOT/'seeds'/tag;d.mkdir(parents=True,exist_ok=True)
    b.write(d/'bracket.csv',[r]);surfaces(d/'surface.csv',[dict(id='root',family='bracket-root',centre=r['centre'],a=r['mean_linear_a_root'],c=r['mean_linear_c_root'])])
    base=(ROOT/'maps'/f'{kind}-N192/params.txt').read_text()
    (d/'params.txt').write_text(replace(base,dict(map_surfaces=d/'surface.csv',map_points=192,map_find='false')))
    rc=measure('root-map',d,[ROOT/'build/t25-expansion-map.ex','params.txt'],3000000000)
    assert rc==0 and 'T25_EXPANSION_MAP_COMPLETE; no advances' in (d/'run.log').read_text()
    (HERE/f't25-root-{tag}.csv').write_text((d/'map.csv').read_text())
def finder(kind,family,points):
    tag=family or kind;sd=ROOT/'seeds'/tag;d=ROOT/'finders'/f'{tag}-N{points}';d.mkdir(parents=True,exist_ok=True)
    (d/'params.txt').write_text(replace((sd/'params.txt').read_text(),dict(map_points=points,map_find='true',map_finder_bracket=sd/'bracket.csv',
        map_find_seconds=235,map_find_max_updates=100000,map_find_floor_window=100001,map_find_chase=1)))
    rc=measure('finder',d,[ROOT/'build/t25-expansion-map.ex','params.txt'],3000000000)
    assert rc in (0,1) and 'T25_FINDER_COMPLETE' in (d/'run.log').read_text()
    (HERE/f't25-finder-{tag}-N{points}.csv').write_text((d/'finder.csv').read_text())
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('operation',choices=['prepare','initialize','map','seed','finder']);p.add_argument('kind',nargs='?',choices=['single','binary','late']);p.add_argument('--points',type=int,default=96);p.add_argument('--family')
    a=p.parse_args()
    if a.operation=='prepare':prepare()
    elif a.operation=='initialize':initialize(a.kind)
    elif a.operation=='map':mapping(a.kind,a.points)
    elif a.operation=='seed':seed(a.kind,a.family)
    else:finder(a.kind,a.family,a.points)
