#!/usr/bin/env python3
"""Bounded author-finder searches on numerical step152; never evolve/initialize."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t24/horizon')
EVIDENCE=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--seconds',type=float,default=90)
    ap.add_argument('--updates',type=int,default=100)
    ap.add_argument('--radii',type=float,nargs='+',default=[.002,.004,.008,.016,.032,.064,.096])
    ap.add_argument('--common',action='store_true',help='numerical midpoint seeds')
    ap.add_argument('--points',type=int,default=48)
    ap.add_argument('--label',default='individual')
    ap.add_argument('--floor-window',type=int,default=32)
    ap.add_argument('--chase',type=float,default=1.)
    a=ap.parse_args()
    root=ROOT if a.label=='individual' else ROOT.parent/('horizon-'+a.label)
    root.mkdir(parents=True,exist_ok=True)
    s=importlib.util.spec_from_file_location('controls',HERE/'t24-controls.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    tracks=list(csv.DictReader((HERE/'t24-last-punctures.csv').open()));t=next(r for r in tracks if float(r['time_M'])==133.)
    radii=a.radii*2;centres=[float(t['x1'])]*len(a.radii)+[float(t['x2'])]*len(a.radii);n=len(radii)
    if a.common:
        radii=a.radii;centres=[.5*(float(t['x1'])+float(t['x2']))]*len(radii);n=len(radii)
    text=(EVIDENCE/'ev1/runs/exp-0024/merger/runtime/params-production.txt').read_text()
    text=m.replace(text,dict(restart_file=EVIDENCE/'chk/EMS_000152.2d.hdf5',ems_data_path='/private/tmp/T24_NO_STATIC_DATA',
        ems_ctt_data_path='/private/tmp/T24_NO_CTT_DATA',RH_activate='true',RH_num_horizons=n,RH_initial_radii=' '.join(map(str,radii)),
        RH_initial_centre=' '.join(map(str,centres)),RH_num_points=' '.join([str(a.points)]*n),RH_level=' '.join(['0']*n),
        RH_time_step_freq=' '.join(['1']*n),RH_chase_speeds=' '.join([str(a.chase)]*n),RH_start_times=' '.join(['0']*n),
        RH_newton_crit=' '.join(['0']*n),offline_points=a.points,offline_max_updates=a.updates,offline_seconds=a.seconds,
        offline_floor_window=a.floor_window,t6_checkpoint_diagnostics='false',t7_diagnostics='false',checkpoint_interval=0,plot_interval=0,
        max_steps=0,verbosity=0,data_path='./',hdf5_path='./',output_path='./'))
    (root/'params.txt').write_text(text)
    env=dict(os.environ,T13_OUTPUT_ROOT='/private/tmp/ems-t24',T13_PROCESS_CAP_BYTES='6000000000',T13_TREE_CAP_BYTES='7500000000',
        T13_DISK_CAP_BYTES='40000000000',OMP_NUM_THREADS='2',TMPDIR='/private/tmp')
    with (root/'run.log').open('wb') as log:
        rc=subprocess.run([sys.executable,str(HERE/'t13-run.py'),'--measure','horizon','--directory',str(root),'--',
            '/private/tmp/ems-t17/rh.ex',str(root/'params.txt')],stdout=log,stderr=subprocess.STDOUT,env=env).returncode
    (root/'done.exit').write_text(str(rc)+'\n')
    q=json.loads((root/'horizon.resources.json').read_text())
    rows=list(csv.DictReader((root/'surfaces.csv').open())) if (root/'surfaces.csv').exists() else []
    for r in rows:
        idx=int(r['search_index']);r['hole']='common' if a.common else 1+idx//len(a.radii);r['seed_radius']=radii[idx]
        shape=root/f'shape-{idx}-{r["stage"]}.dat'
        if shape.exists():
            import numpy as np
            v=np.loadtxt(shape);f=v[2:];theta=(np.arange(len(f))+.5)*np.pi/len(f)
            r.update(rho_max=float(np.max(f*np.sin(theta))),axial_min=float(np.min(v[1]+f*np.cos(theta))),axial_max=float(np.max(v[1]+f*np.cos(theta))))
    if rows:
        with (HERE/('t24-horizon-search.csv' if a.label=='individual' else f't24-horizon-{a.label}.csv')).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    receipt=dict(native=q['child_measurement'],resources=q,rows=rows,radii=a.radii,centres=centres,
        advances=0,static_data_reads=0,scope=f'fresh numerical-radius spheres, unchanged RHUnion/RHSurf, chase={a.chase}/quota=1/N{a.points}',
        status='NO_GROUP_QUALIFICATION_IN_BOUNDED_SEARCH' if rc else f'ALL_THREE_STAGES_N{a.points}_FOUND',
        limit='N96 agreement not tested; failure to find does not prove absence of an apparent horizon')
    (HERE/('t24-horizon-search.json' if a.label=='individual' else f't24-horizon-{a.label}.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    assert not q['gate_reason'] and q['peak_rss_bytes']<6e9 and q['returncode'] in (0,1,2),(rc,q)
    print(receipt['status'],q['peak_rss_bytes'],len(rows),flush=True)
if __name__=='__main__':main()
