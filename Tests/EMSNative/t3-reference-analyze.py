#!/usr/bin/env python3
"""Analyze the registered reference evolution using the T2 mask/class rules."""
import csv
import importlib.util
import math
import sys
from pathlib import Path

import h5py
import numpy as np

here = Path(__file__).parent
spec = importlib.util.spec_from_file_location('t2_localize', here / 't2-localize.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
classes = ('all', 'convex_corner', 'patch_edge', 'axis_patch_junction',
           'cartoon_axis', 'box_seam', 'outer_boundary', 'interior')
runroot = Path(sys.argv[1])
rows, layouts, gaussb = [], [], []
for scale, count, interval in (('low', 40, 4), ('mid', 60, 6), ('high', 90, 9)):
    for flag in ('off', 'on'):
        run = f't3-ref-match-{scale}-{flag}'
        if scale == 'high':
            marker = runroot / run / 'done.exit'
            if not marker.exists() or marker.read_text().strip() != '0':
                continue
        for step in range(0, count + 1, interval):
            path = runroot / run / 'plt' / f'EMS_Plot_{step:06d}.2d.hdf5'
            with h5py.File(path) as f:
                time = float(f['level_0'].attrs['time'])
                bmax = 0.
                for level in range(int(f.attrs['num_levels'])):
                    g = f[f'level_{level}']
                    offsets = g['data:offsets=0'][:]
                    for bi, box in enumerate(g['boxes'][:]):
                        cells = (int(box['hi_i'])-int(box['lo_i'])+1) * (int(box['hi_j'])-int(box['lo_j'])+1)
                        start = int(offsets[bi]) + 8*cells
                        bmax = max(bmax, float(np.max(np.abs(g['data:datatype=0'][start:start+cells]))))
                if step == 0:
                    for level in range(int(f.attrs['num_levels'])):
                        g = f[f'level_{level}']
                        boxes = g['boxes'][:]
                        a = np.array([[int(b[k]) for k in ('lo_i','lo_j','hi_i','hi_j')] for b in boxes])
                        lo, hi = a[:,:2].min(axis=0), a[:,2:].max(axis=0)
                        h = float(g.attrs['dx'])
                        layouts.append((run,level,h,len(boxes),lo[0]*h-256,(hi[0]+1)*h-256,lo[1]*h,(hi[1]+1)*h))
            gaussb.append((run,scale,flag,step,time,bmax))
            for _, _, mask, cell_class, field, n, rms, maximum in module.measure(run, path, 'R', 256., classes):
                rows.append((run,scale,flag,step,time,mask,cell_class,field,n,rms,maximum))
with (here/'t3-reference-timeseries.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(('run','scale','restriction','coarse_step','time_M','mask','cell_class','constraint','cells','cylindrical_rms','maximum'));w.writerows(rows)
with (here/'t3-reference-layout.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(('run','level','dx','boxes','xlo_M','xhi_M','ylo_M','yhi_M'));w.writerows(layouts)
with (here/'t3-reference-gaussb.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(('run','scale','restriction','coarse_step','time_M','global_max_abs_GaussB'));w.writerows(gaussb)
lookup={(flag,scale,round(t,8),mask,cl,field):(n,rms) for _,scale,flag,_,t,mask,cl,field,n,rms,_ in rows}
with (here/'t3-reference-orders.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(('restriction','time_M','mask','cell_class','constraint','low_rms','mid_rms','high_rms','order_low_mid','order_mid_high','low_cells','mid_cells','high_cells'))
    for flag in ('off','on'):
        for t in (1.,2.5,5.):
            keys=sorted({(mask,cl,field) for _,_,fl,_,tt,mask,cl,field,_,_,_ in rows if fl==flag and abs(tt-t)<1e-8})
            for mask,cl,field in keys:
                v=[lookup.get((flag,scale,t,mask,cl,field)) for scale in ('low','mid','high')]
                if not all(v[:2]):continue
                r=[x[1] if x else '' for x in v]
                o=[math.log(r[0]/r[1])/math.log(1.5) if r[0]>0 and r[1]>0 else float('nan'),
                   math.log(r[1]/r[2])/math.log(1.5) if v[2] and r[1]>0 and r[2]>0 else '']
                w.writerow((flag,t,mask,cl,field,*r,*o,*(x[0] if x else '' for x in v)))
print('analyzed',len(rows),'time-series norms and',len(layouts),'levels')
