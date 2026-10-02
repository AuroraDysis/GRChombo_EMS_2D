#!/usr/bin/env python3
"""Compare saved native t0 streams across the full and cropped launch domains."""
import sys
sys.dont_write_bytecode=True
import csv,json,re
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t18c')
def layout(scope):
    records=list(csv.DictReader((ROOT/'evolution'/scope/'t18c-ranges.csv').open()));offset=0;result={}
    for r in records:
        l=int(r['level']);box=tuple(int(r[k]) for k in ('x0','y0','x1','y1'));n=int(r['nonlapse_values'])
        result[(l,*box)]=(offset,n,float(r['h']));offset+=n
    p=ROOT/'evolution'/scope/'t18c-nonlapse.bin';assert p.stat().st_size==offset*8
    return result,np.memmap(p,dtype='f8',mode='r',shape=(offset,))
def main():
    full,u=layout('full-sqrt');local,v=layout('launch-sqrt');values=bits=boxes=0
    for (l,x0,y0,x1,y1),(off,n,h) in local.items():
        shift=1024*2**l;key=(l,x0+shift,y0,x1+shift,y1)
        if l==0 and key not in full:continue # cropped physical boundary boxes
        assert key in full
        oldoff,oldn,oldh=full[key];assert (oldn,oldh)==(n,h)
        nx=x1-x0+7;ny=y1-y0+7
        x=u[oldoff:oldoff+n].reshape(ny,nx,27);y=v[off:off+n].reshape(ny,nx,27)
        if l==0:
            X=(np.arange(x0-3,x1+4)+.5)*h-448
            Y=(np.arange(y0-3,y1+4)+.5)*h
            mask=(abs(X)[None,:]<=200)&(abs(Y)[:,None]<=200)
            x=x[mask];y=y[mask]
        values+=x.size;bits+=int(np.count_nonzero(x.view('u8')!=y.view('u8')));boxes+=1
    assert values>0 and bits==0
    counts={}
    for mode in ('sqrt','geometric'):
        log=(ROOT/'evolution'/('launch-'+mode)/'run.log').read_text()
        z=[len(re.findall(r'GRAMRLevel::advance level '+str(l)+r' at time ',log)) for l in range(13)]
        assert min(z)>0;counts[mode]=z
    # A deliberately generous footprint estimate, not a full evolved-domain
    # bit comparison: per level allow 32 cells per started RK update plus 16
    # for initial derivatives, point interpolation and diagnostic stencils.
    bound=max(sum((32*n+16)*1.75/2**l for l,n in enumerate(z)) for z in counts.values())
    assert bound+8+.01<448
    result=dict(full_vs_cropped_initial_nonlapse_values=values,bit_mismatches=bits,compared_boxes=boxes,
        scope='all ghosts/valid cells on every refined box; base cells in abs(x-midpoint)<=200, abs(y)<=200',
        native_advance_counts=counts,conservative_coordinate_footprint_estimate_M=bound,
        cropped_outer_boundary_distance_from_nearer_puncture_M=440,
        limit='No full-production evolved-state comparison; the captured finest t0 comparison separately checks all 28 fields, including lapse.')
    (HERE/'t18c-domain-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
