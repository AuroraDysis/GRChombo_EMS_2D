#!/usr/bin/env python3
"""Lossless recorder decoder and runnable production/control checks."""
import csv,gzip,lzma,struct,importlib.util,sys,os,subprocess
from pathlib import Path
from collections import Counter
import numpy as np
import h5py
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t7')
spec=importlib.util.spec_from_file_location('t6check',HERE/'t6-check.py')
t6=importlib.util.module_from_spec(spec);spec.loader.exec_module(t6)

def frames(path):
    previous={}
    with (lzma.open(path,'rb') if path.suffix=='.xz' else gzip.open(path,'rb')) as f:
        assert f.read(8)==b'T7OP0002'
        while raw:=f.read(44):
            assert len(raw)==44
            p,lev,src,stage,n,count,delta,x0,y0,x1,y1=struct.unpack('<11i',raw)
            meta=np.frombuffer(f.read(80),dtype='<f8').copy()
            cells=np.frombuffer(f.read(count*8),dtype='<i4').reshape(count,2).copy()
            key=(lev,src,n,p if p==40 or 20<=p<40 else 0,cells.tobytes())
            bits=np.frombuffer(f.read(count*n*8),dtype='u1').reshape(8,-1).T.copy().view('<u8').ravel()
            if delta:bits^=previous[key]
            previous[key]=bits.copy();values=bits.view('<f8').reshape(count,n)
            assert np.isfinite(values).all() and np.isfinite(meta).all()
            yield dict(phase=p,level=lev,source=src,stage=stage,meta=meta,cells=cells,values=values,valid=(x0,y0,x1,y1))

def save(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def recorder(d):
    rows=[]
    for path in sorted(list(d.glob('t7-stage*.gz'))+list(d.glob('t7-stage*.xz'))):
        counts=Counter();states={};filled={};parts={};parent={};worst=0.;prolong_error=0.;prolong_count=0;deps=0;parentdeps=0
        for fr in frames(path):
            p,src,m=fr['phase'],fr['source'],fr['meta'];counts[p]+=1
            if p in (2,3,40):
                lookup={tuple(iv):v for iv,v in zip(fr['cells'],fr['values'])}
                (states if p==3 else filled if p==2 else parent)[src]=(fr,lookup)
            if p in (20,21):parts[(src,p)]=fr
            if p==20:
                st,lookup=states[src];assert st['stage']==fr['stage'] and st['meta'][3]==m[3]
                for iv in fr['cells']:
                    for dx in range(-3,4):
                        for dy in range(-3,4):
                            if abs(dx)>2 or abs(dy)>2:
                                if not ((abs(dx)==3 and dy==0) or (abs(dy)==3 and dx==0)):continue
                            assert (iv[0]+dx,iv[1]+dy) in lookup;deps+=1
                # Own-face windows use dense parent support. Child-face windows use same-level state.
                if src%16<3 and src in parent:
                    pf,pl=parent[src];assert pf['stage']==fr['stage'] and pf['meta'][3]==m[3]
                    x0,y0,x1,y1=st['valid'];own=m[8];h=m[1];centre=256/h
                    for iv in st['cells']:
                        x=(iv[0]+.5-centre)*h;y=(iv[1]+.5)*h
                        if y<0 or (abs(x)<own and y<own):continue
                        start=np.floor((iv+.5)/2-.5).astype(int)-2
                        q=(iv+.5)/2-.5
                        w=[np.array([np.prod([(q[z]-start[z]-j)/(i-j) for j in range(6) if j!=i]) for i in range(6)]) for z in (0,1)]
                        base=pl[tuple(start+2)];correction=np.zeros(28)
                        for b in range(6):
                            for a in range(6):
                                assert (start[0]+a,start[1]+b) in pl;parentdeps+=1
                                correction+=w[0][a]*w[1][b]*(pl[(start[0]+a,start[1]+b)]-base)
                        actual=filled[src][1][tuple(iv)]
                        prolong_error=max(prolong_error,float(np.max(abs(base+correction-actual)/np.maximum(1.,abs(actual)))))
                        prolong_count+=1
            if p==22:
                a,b=parts[(src,20)],parts[(src,21)]
                assert a['cells'].tobytes()==b['cells'].tobytes()==fr['cells'].tobytes()
                err=abs(fr['values']-(a['values']+b['values']))
                scale=abs(a['values'])+abs(b['values'])+np.finfo(float).tiny
                worst=max(worst,float(np.max(err/np.maximum(1.,scale))))
        assert worst<16*np.finfo(float).eps and prolong_error<32*np.finfo(float).eps and deps>0
        rows.append(dict(case=d.name,file=path.name,frames=sum(counts.values()),phase_counts=str(dict(sorted(counts.items()))),
                         finite_Float64='PASS',RHS_sum_relative_max=worst,fine_stencil_cells_checked=deps,coarse_support_cells_checked=parentdeps,
                         compressed_bytes=path.stat().st_size,prolongation_ghost_cells_checked=prolong_count,prolongation_replay_max_error=prolong_error))
    for path in sorted(list(d.glob('t7-snapshot*.gz'))+list(d.glob('t7-snapshot*.xz'))):
        n=0
        for fr in frames(path):assert fr['phase']==50 and fr['values'].shape[1]==33;n+=1
        rows.append(dict(case=d.name,file=path.name,frames=n,phase_counts='{50: '+str(n)+'}',finite_Float64='PASS',
                         RHS_sum_relative_max=0,fine_stencil_cells_checked=0,coarse_support_cells_checked=0,compressed_bytes=path.stat().st_size,prolongation_ghost_cells_checked=0,prolongation_replay_max_error=0))
    return rows

def controls():
    rows=[]
    for case in ('legacy-E','legacy-ref','point-ref'):
        base=ROOT/'controls'/case/'baseline';new=ROOT/'controls'/case/'default-final'
        for kind in ('plt','chk'):
            aa,bb=sorted((base/kind).glob('*.hdf5')),sorted((new/kind).glob('*.hdf5'))
            assert len(aa)==len(bb)>0
            for a,b in zip(aa,bb):assert a.name==b.name;t6.equal_hdf(a,b)
            rows.append(dict(control=case,variant='default-off',kind=kind,files=len(aa),result='BIT_IDENTICAL'))
    # Recording alone must not alter either the initial plot or final two-step checkpoint.
    for case in ('clock-final-xz','space-final-xz'):
        for kind in ('plt','chk'):
            aa=sorted((ROOT/'probes'/case/kind).glob('*.hdf5'))
            bb=sorted((ROOT/'controls'/case.replace('-final-xz','-final')/'baseline'/kind).glob('*.hdf5'))
            assert len(aa)==len(bb)>0
            for a,b in zip(aa,bb):assert a.name==b.name;t6.equal_hdf(a,b)
            rows.append(dict(control=case,variant='recording-on',kind=kind,files=len(aa),result='BIT_IDENTICAL'))
    def physical(d):
        return [{k:v for k,v in r.items() if k!='seconds'} for r in csv.DictReader((d/'surfaces.csv').open())]
    a=Path('/private/tmp/ems-t6/rh-control-default-final');b=ROOT/'controls/finder-default'
    assert physical(a)==physical(b)
    for f in a.glob('shape-*.dat'):assert f.read_bytes()==(b/f.name).read_bytes()
    rows.append(dict(control='frozen-E-terminal-finder',variant='default-off',kind='physical_CSV_and_shapes',files=1+len(list(a.glob('shape-*.dat'))),result='BIT_IDENTICAL'))
    save('t7-controls.csv',rows)


def replay(d):
    # First four synchronization times, all three regions and covered boundary cells.
    rows=[]
    for level in (5,6):
        fine={};times=[]
        path=d/f't7-stage-L{level}.xz'
        for fr in frames(path):
            if fr['phase']!=10:continue
            t=float(fr['meta'][0])
            if t not in times and len(times)<4:times.append(t)
            if t in times:fine.setdefault(t,[]).append(fr)
        count=0;worst=0.
        weights=np.array([np.prod([(.5-j)/(i-j) for j in range(-2,4) if j!=i]) for i in range(-2,4)])
        for fr in frames(d/f't7-stage-L{level-1}.xz'):
            if fr['phase']!=11:continue
            t=float(fr['meta'][0]);match=next((q for q in times if abs(q-t)<1e-13),None)
            if match is None:continue
            for iv,value in zip(fr['cells'],fr['values']):
                owner=next((q for q in fine[match] if q['valid'][0]<=2*iv[0]<=q['valid'][2] and q['valid'][1]<=2*iv[1]<=q['valid'][3]),None)
                if owner is None:continue
                lookup={tuple(x):v for x,v in zip(owner['cells'],owner['values'])}
                # Only compare the parent's child-face core captured by fine support.
                anchor=tuple(2*iv)
                if anchor not in lookup:continue
                base=lookup[anchor];correction=np.zeros(28)
                for y in range(6):
                    for x in range(6):
                        cell=tuple(2*iv+np.array([x-2,y-2]));assert cell in lookup
                        correction+=weights[x]*weights[y]*(lookup[cell]-base)
                worst=max(worst,float(np.max(abs(base+correction-value)/np.maximum(1.,abs(value)))));count+=1
        assert count>0 and worst<32*np.finfo(float).eps
        rows.append(dict(case=d.name,operation='restriction_outermost_covered_row_and_corner',fine_level=level,
                         times_checked=len(times),cells_checked=count,max_scaled_error=worst,result='PASS'))
    return rows

def reproduce_controls():
    # Existing parameter files are the concrete bounded test fixtures.
    jobs=[]
    for case in ('legacy-E','legacy-ref','point-ref'):
        for variant in ('baseline','default-final'):
            jobs.append((ROOT/'controls'/case/variant,ROOT.parent/'ems-t6/evolution.ex' if variant=='baseline' else ROOT/'evolution.ex'))
    for case in ('clock-final-xz','space-final-xz'):
        jobs.append((ROOT/'controls'/case.replace('-final-xz','-final')/'baseline',ROOT.parent/'ems-t6/evolution.ex'))
        jobs.append((ROOT/'probes'/case,ROOT/'evolution.ex'))
    for d,exe in jobs:
        for kind in ('plt','chk'):(d/kind).mkdir(exist_ok=True)
        with (d/'run.log').open('w') as log:
            result=subprocess.run([str(exe),'params.txt'],cwd=d,env={**os.environ,'OMP_NUM_THREADS':'4'},stdout=log,stderr=subprocess.STDOUT,timeout=175)
        assert result.returncode==0

if __name__=='__main__':
    if '--reproduce-controls' in sys.argv:reproduce_controls()
    if '--recorder-only' not in sys.argv:controls()
    rows=[]
    for case in ('clock-final-xz','space-final-xz'):rows+=recorder(ROOT/'probes'/case)
    save('t7-recorder-checks.csv',rows)
    rr=[]
    for case in ('clock-final-xz','space-final-xz'):rr+=replay(ROOT/'probes'/case)
    save('t7-operation-replay.csv',rr)
    print('PASS: bit identity, lossless finite frames, actual pre-KO/KO/total RHS and complete spatial support')
