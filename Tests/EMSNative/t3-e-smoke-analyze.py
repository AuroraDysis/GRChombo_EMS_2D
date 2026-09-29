#!/usr/bin/env python3
"""Compare E's restricted two-step class residuals with T2 controls."""
import csv
import importlib.util
import sys
from pathlib import Path

import h5py

here=Path(__file__).parent
spec=importlib.util.spec_from_file_location('t2_localize',here/'t2-localize.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
classes=('all','convex_corner','patch_edge','axis_patch_junction','cartoon_axis','box_seam','outer_boundary','interior')
runroot=Path(sys.argv[1]);rows=[]
for scale in ('low','mid'):
    run=f't3-E-{scale}-restrict-2step'
    for step in (0,1,2):
        path=runroot/run/'plt'/f'EMS_Plot_{step:06d}.2d.hdf5'
        if not path.exists():
            continue
        with h5py.File(path) as f: time=float(f['level_0'].attrs['time'])
        for _,_,mask,cl,field,n,rms,maximum in module.measure(run,path,'E',224.,classes):
            rows.append((scale,step,time,mask,cl,field,n,rms,maximum))
with (here/'t3-e-smoke-localization.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(('scale','coarse_step','time_M','mask','cell_class','constraint','cells','cylindrical_rms','maximum'));w.writerows(rows)
t2=list(csv.DictReader((here/'t2-smoke-comparison.csv').open()))
lookup={(x['mask'],x['cell_class'],x['constraint']):x for x in t2}
with (here/'t3-e-class-comparison.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(('scale','coarse_step','time_M','mask','cell_class','constraint','restricted_rms','t2_unfixed_low_rms','t2_pointwise_low_rms','restricted_over_unfixed_low','restricted_over_pointwise_low'))
    for scale,step,time,mask,cl,field,_,rms,_ in rows:
        control=lookup.get((mask,cl,field)) if scale=='low' else None
        old=float(control[f'unfixed_{"t0" if step==0 else f"step{step}"}']) if control else None
        point=float(control[f'pointwise_{"t0" if step==0 else f"step{step}"}']) if control else None
        w.writerow((scale,step,time,mask,cl,field,rms,old,point,rms/old if old else '',rms/point if point else ''))
print('analyzed',len(rows),'E class norms')
