#!/usr/bin/env python3
"""Recheck T6 default controls and the bounded, completed strip probes."""
import csv
import struct
from pathlib import Path

import h5py
import numpy as np

ROOT = Path('/private/tmp/ems-t6')
HERE = Path(__file__).resolve().parent


def equal_value(a, b):
    a,b = np.asarray(a),np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape
    assert np.array_equal(a,b) if a.dtype.hasobject else a.tobytes() == b.tobytes()


def equal_hdf(a, b):
    with h5py.File(a,'r') as f,h5py.File(b,'r') as g:
        def check(name):
            u,v = f[name],g[name]
            assert type(u)==type(v) and set(u.attrs)==set(v.attrs)
            for k in u.attrs: equal_value(u.attrs[k],v.attrs[k])
            if isinstance(u,h5py.Dataset): equal_value(u[:],v[:])
            else: assert set(u)==set(v)
        check('/');f.visit(check)


def strip(path):
    with path.open('rb') as f:
        assert f.read(8)==b'T6STRIP1'
        time,h,face,child = struct.unpack('<4d',f.read(32))
        level,n,ncol = struct.unpack('<3I',f.read(12))
        assert n==28 and ncol==36 and level in (4,5,6) and h>0
        values = np.fromfile(f,dtype='<f8').reshape(-1,ncol)
    assert len(values) and np.isfinite(values).all()
    assert set(values[:,2]) <= {0.,1.} and np.all(values[:,-1]==0.)
    iv = np.rint(np.column_stack(((values[:,0]+256)/h-.5,values[:,1]/h-.5))).astype(int)
    assert len(np.unique(iv,axis=0))==len(iv)
    return time,h,level,iv,values


def main():
    rows = []
    controls = [('reference_mid_7_levels_3_steps','control-ref',('default','strip-on')),
                ('E_legacy_regrid_2_steps','control-legacy-E',('default',)),
                ('ref_legacy_regrid_2_steps','control-legacy-ref',('default',))]
    for label,directory,variants in controls:
        for variant in variants:
            for kind in ('plt','chk'):
                a = sorted((ROOT/directory/'baseline'/kind).glob('*.hdf5'))
                b = sorted((ROOT/directory/variant/kind).glob('*.hdf5'))
                assert len(a)==len(b)>0
                for u,v in zip(a,b): assert u.name==v.name;equal_hdf(u,v)
                rows.append(dict(control=label,variant=variant,kind=kind,files=len(a),result='BIT_IDENTICAL'))
    a,b = ROOT/'rh-pilot',ROOT/'rh-control-default-final'
    def physical(d):
        return [{k:v for k,v in r.items() if k!='seconds'} for r in csv.DictReader((d/'surfaces.csv').open())]
    assert physical(a)==physical(b)
    for f in a.glob('shape-*.dat'): assert f.read_bytes()==(b/f.name).read_bytes()
    rows.append(dict(control='frozen_E_low_terminal_400_chase_iterations',variant='default',
                     kind='physical_CSV_and_shapes',files=1+len(list(a.glob('shape-*.dat'))),result='BIT_IDENTICAL'))
    for d in (ROOT/'control-ref'/'strip-on',ROOT/'probes'/'face3-smoke'):
        paths = sorted((d/'strips').glob('*.bin'))
        assert len(paths)==9
        for path in paths: strip(path)
        plot = sorted((d/'plt').glob('*.hdf5'))[-1]
        with h5py.File(plot,'r') as f:
            assert int(f.attrs['num_components'])==33
            for path in paths:
                time,h,lev,iv,v = strip(path)
                if abs(time-.25)>1e-12: continue
                g = f[f'level_{lev}'];assert h==float(g.attrs['dx'])
                offsets = g['data:offsets=0'][:]
                seen = np.zeros(len(v),bool)
                for i,box in enumerate(g['boxes'][:]):
                    x0,y0,x1,y1 = (int(box[k]) for k in ('lo_i','lo_j','hi_i','hi_j'))
                    use = (iv[:,0]>=x0)&(iv[:,0]<=x1)&(iv[:,1]>=y0)&(iv[:,1]<=y1)&(v[:,2]==0)
                    data = g['data:datatype=0'][int(offsets[i]):int(offsets[i+1])].reshape(33,y1-y0+1,x1-x0+1)
                    equal_value(v[use,3:31],data[:28,iv[use,1]-y0,iv[use,0]-x0].T)
                    seen |= use
                assert np.all(seen[v[:,2]==0])
        print(f'PASS {d.name}: 9 strip frames, finite Float64, GaussB zero; terminal uncovered evolved values bit-identical to plot')
    with (HERE/'t6-controls.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    print('PASS: all HDF5 datasets/attributes and default finder physical outputs bit-identical')


if __name__=='__main__': main()
