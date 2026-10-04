#!/usr/bin/env python3
"""Exact native-valid annular support and uncompressed stage recording sizes."""
import argparse
import csv
import resource
from pathlib import Path
import h5py
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint',type=Path)
    p.add_argument('--boxes-csv',type=Path)
    p.add_argument('--radii',type=float,nargs=2,default=(.04,.09))
    p.add_argument('--levels',type=int,nargs='+',default=(11,12))
    p.add_argument('--centres',type=float,nargs='+',default=(2239.7740195,2240.2259805))
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'t24-annulus-size.csv')
    a=p.parse_args();rr=[]
    saved=list(csv.DictReader(a.boxes_csv.open())) if a.boxes_csv else None
    with h5py.File(a.checkpoint,'r') if saved is None else __import__('contextlib').nullcontext(None) as f:
        for level in a.levels:
            boxes=[r for r in saved if int(r['level'])==level] if saved is not None else f[f'level_{level}/boxes'][:]
            dx=float(boxes[0]['h_M']) if saved is not None else float(f[f'level_{level}'].attrs['dx'])
            n=0;per_hole=np.zeros(len(a.centres),dtype=int)
            for b in boxes:
                x0,y0,x1,y1=[int(b[k]) for k in ('lo_i','lo_j','hi_i','hi_j')]
                xx,yy=np.meshgrid((np.arange(x0,x1+1)+.5)*dx,(np.arange(y0,y1+1)+.5)*dx)
                union=np.zeros(xx.shape,bool)
                for i,c in enumerate(a.centres):
                    r2=(xx-c)**2+yy*yy;take=(r2>=a.radii[0]**2)&(r2<=a.radii[1]**2)
                    union|=take;per_hole[i]+=take.sum()
                n+=union.sum()
            # 15 Float64 fields; two int32 cell coordinates + uint32 hole mask.
            per_stage=int(n)*132
            rr.append(dict(level=level,dx=dx,valid_union_cells=int(n),cells_per_hole=';'.join(map(str,per_hole)),
                payload_bytes_per_stage=per_stage,RK_stages_per_coarse_step=4*2**level,
                payload_bytes_per_coarse_step=per_stage*4*2**level,
                window_duration=.041015625,window_fine_steps=int(.041015625/(.5*dx)),
                payload_bytes_per_window=per_stage*4*int(.041015625/(.5*dx)),
                peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    with a.output.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rr[0]);w.writeheader();w.writerows(rr)
    print(rr)
if __name__=='__main__':main()
