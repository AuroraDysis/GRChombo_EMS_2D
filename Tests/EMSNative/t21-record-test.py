#!/usr/bin/env python3
"""Codec/sampling failure checks, independent of any EMS trajectory."""
import importlib.util, json, resource, struct, tempfile
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('record',HERE/'t21-record.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
for n in (8,10):
    for x in (-.5,.0,.125,.5,1.25):
        nodes,w=m.weights(x,n)
        for degree in range(n):
            actual=w@((nodes/n)**degree);expected=(x/n)**degree
            assert abs(actual-expected)<1e-12,(n,x,degree,actual,expected)
with tempfile.TemporaryDirectory(prefix='t21-codec-',dir='/private/tmp') as d:
    directory=Path(d);p=directory/'native.bin'
    ij=np.array([(i,j) for j in range(21) for i in range(-20,21)],dtype='i4')
    a=np.zeros((len(ij),19));x=ij[:,0]+.5;y=ij[:,1]+.5
    a[:,0]=x*x+y*y;a[:,1]=y;a[:,-1]=1
    meta=dict(fields=['phi','B2'],parities=[0,2],gauge='fixture',driver_semantics='test_driver',native_kernel_type='fixture_only',geometry=m.GEOMETRY)
    p.with_suffix('.json').write_text(json.dumps(meta))
    dtype=np.dtype([('ij','i4',(2,)),('values','f8',(19,))]);records=np.empty(len(ij),dtype=dtype);records['ij']=ij;records['values']=a
    p.write_bytes(b'T21RHS01'+struct.pack('=7I',0,0,0,1,2,len(ij),1)+struct.pack('=6d',0,0,1,.25,0,0)+records.tobytes())
    assert len(list(m.frames(p)))==1
    clock,meta,grids=m.load_grids(directory)
    for n in (8,10):
        for x,y in ((.75,0.),(.25,1.25),(1.3,2.3)):
            value=m.sample(grids[0],x,y,n)
            assert abs(value[0]-(x*x+y*y))<1e-12,(value,x,y)
            assert abs(value[1]-y)<1e-12,(value,x,y)
        assert m.sample(grids[0],100,100,n) is None
    assert m.labels(meta)[1]=='state:B2[fixture:test_driver]'
    bad=directory/'bad.bin';bad.write_bytes(p.read_bytes()[:-1]);bad.with_suffix('.json').write_text(json.dumps(meta))
    try:list(m.frames(bad));raise AssertionError('truncation was accepted')
    except ValueError:pass
    bad.unlink();bad.with_suffix('.json').unlink()
    original=p.read_bytes();corrupt=bytearray(original);corrupt[92:100]=struct.pack('=d',float('nan'));p.write_bytes(corrupt)
    try:list(m.frames(p));raise AssertionError('nonfinite native input was accepted')
    except ValueError:pass
    p.write_bytes(original)
    other=directory/'other.bin';other.write_bytes(p.read_bytes());meta['gauge']='other_gauge';other.with_suffix('.json').write_text(json.dumps(meta))
    try:m.load_grids(directory);raise AssertionError('mixed gauges were accepted')
    except ValueError:pass
receipt=dict(status='PASS',checks=['degree<n I8/I10 polynomials','axis parity on native-cell data','missing support rejected','truncated payload rejected','nonfinite native state/RHS rejected','mixed gauge metadata rejected','driver labels include gauge semantics'],
    peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='codec/sampler synthetic fixture only; no EMS or evolution claim')
(HERE/'t21-record-test.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
