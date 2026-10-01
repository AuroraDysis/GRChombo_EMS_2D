#!/Users/auroradysis/miniconda3/bin/python
"""Output-only lossless T7 filter: common clocks and sufficient W stencils."""
import lzma, struct, sys
import numpy as np
src=sys.stdin.buffer;dst=sys.stdout.buffer
compressor=lzma.LZMACompressor(format=lzma.FORMAT_XZ,preset=3)
def emit(a):dst.write(compressor.compress(a))
def read(n):
    a=src.read(n);assert len(a)==n,(len(a),n);return a
assert read(8)==b'T7OP0002';emit(b'T7OP0002')
before={};after={};clock=.875/8192/4
while head:=src.read(44):
    assert len(head)==44
    p,lev,source,stage,n,count,delta,x0,y0,x1,y1=struct.unpack('<11i',head)
    meta=read(80);t,h=struct.unpack('<2d',meta[:16])
    coords=read(count*8);iv=np.frombuffer(coords,dtype='<i4').reshape(count,2)
    key=(lev,source,n,p if p==40 or 20<=p<40 else 0,coords)
    payload=np.frombuffer(read(count*n*8),dtype='u1').copy()
    if delta:payload^=before[key]
    before[key]=payload.copy()
    if p==50 and abs(t/clock-round(t/clock))>1e-7:continue
    if p==50:
        take=(abs((iv[:,0]+.5)*h-336)<=.0035)&(abs((iv[:,1]+.5)*h)<=.0035)
        if not np.any(take):continue
        coords=iv[take].tobytes();payload=payload.reshape(8,count,n)[:,take,:].reshape(-1).copy();count=int(take.sum())
    key=(lev,source,n,p if p==40 or 20<=p<40 else 0,coords)
    previous=after.get(key);encoded=payload.copy();delta=int(previous is not None)
    if delta:encoded^=previous
    after[key]=payload.copy()
    emit(struct.pack('<11i',p,lev,source,stage,n,count,delta,x0,y0,x1,y1));emit(meta);emit(coords);emit(encoded.tobytes())
dst.write(compressor.flush());dst.flush()
