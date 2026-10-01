#!/usr/bin/env python3
"""Bounded current-state analysis; never opens the static trumpet."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, json, subprocess
from pathlib import Path
import numpy as np
import importlib.util
import lzma
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t16');OUT=ROOT/'analysis';OUT.mkdir(exist_ok=True)
PY='/Users/auroradysis/miniconda3/bin/python'
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
a=module('t13analysis',HERE/'t13-analyze.py')
def save(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def timed(label,cmd):
    subprocess.run([PY,str(HERE/'t16-run.py'),'--measure',label,'--directory',str(OUT),'--',*cmd],check=True)
def select(h):
    cells=set()
    for v in a.NV.values():
        xy=a.R[:,None]*v;z=np.c_[(336+xy[:,0])/h-.5,xy[:,1]/h-.5]
        for s in np.floor(z).astype(int)-4:
            for j in range(10):
                yy=s[1]+j;yy=yy if yy>=0 else -yy-1
                for i in range(10):cells.add((s[0]+i,yy))
    c=round(336/h)
    for j in range(6):
        for i in range(-5,6):cells.add((c+i,j))
    return np.array(sorted(cells,key=lambda p:(p[1],p[0])),int)
def extract(name,level):
    h=.875/2**level;d=ROOT/'evolution'/name;assert (d/'done.exit').read_text().strip()=='0'
    cells=select(h);times=[];payload=OUT/(name+'-input.bin');output=OUT/(name+'-q.bin');group=[];previous=None
    with payload.open('wb') as stream:
        def consume(frames):
            t=frames[0]['meta'][0];lo,u=a.dense(frames)
            rec=np.empty((len(cells),3+49*28));rec[:,0]=h
            rec[:,1]=(cells[:,1]+.5)*h;rec[:,2]=(cells[:,0]+.5)*h-336
            for j in range(-3,4):
                for i in range(-3,4):
                    q=u[cells[:,1]+j-lo[1],cells[:,0]+i-lo[0]]
                    assert np.isfinite(q).all(),(name,t,i,j)
                    rec[:,3+(j+3)*7+i+3::49]=q
            rec.tofile(stream);times.append(t)
        for fr in a.decoder.frames(d/f't13-t7-stage-L{level}.xz'):
            if fr['phase']!=50:continue
            t=fr['meta'][0]
            if previous is not None and t!=previous:consume(group);group=[]
            group.append(fr);previous=t
        if group:consume(group)
    timed('replay-'+name,[str(ROOT/'replay.ex'),'--state',str(payload),str(output)])
    digest=hashlib.sha256()
    with payload.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):digest.update(chunk)
    (OUT/(name+'-input.sha256')).write_text(digest.hexdigest()+' '+str(payload)+'\n');payload.unlink()
    q=np.memmap(output,mode='r',dtype='f8',shape=(len(times),len(cells),180))
    assert np.isfinite(q).all() and q[:,:,144].max()==0
    np.savez(OUT/(name+'-meta.npz'),cells=cells,times=times)
    a.H=h
    for ray,v in a.NV.items():
        p={n:np.array([a.sample(cells,frame,v,n) for frame in q]) for n in (8,10)}
        np.savez_compressed(OUT/(name+'-'+ray+'-profiles.npz'),times=times,r=a.R,p8=p[8],p10=p[10])
    print(name,'snapshots',len(times),'cells',len(cells),'q bytes',output.stat().st_size,flush=True)
def measure():
    history=[];audits=[];ratios=[];native=[];identity=[];cached={}
    for level in (13,14):
        for mode in ('sqrt','geometric'):
            name=f'L{level}-{mode}';d=ROOT/'evolution'/name
            meta=np.load(OUT/(name+'-meta.npz'));times=meta['times'];cells=meta['cells']
            q=np.memmap(OUT/(name+'-q.bin'),mode='r',dtype='f8',shape=(len(times),len(cells),180))
            floor=[]
            for f in d.glob('t13-floors-*.csv'):floor.extend(csv.DictReader(f.open()))
            log=(d/'run.log').read_text();first=log.index('GRAMRLevel::advance level 0 at time 0')
            r=json.loads((d/(name+'.resources.json')).read_text())
            det=q[:,:,1]*q[:,:,3]-q[:,:,2]**2
            inv=(q[:,:,1]+q[:,:,3]+np.hypot(q[:,:,1]-q[:,:,3],2*q[:,:,2]))/(2*det)
            speed=float(np.max(np.hypot(q[:,:,14],q[:,:,15])+np.sqrt(np.maximum(inv,1/q[:,:,4]))))
            actual_stop=float(next(csv.DictReader((d/'t13-stop.csv').open()))['actual_time_M'])
            audits.append(dict(run=name,exit=0,snapshots=len(times),native_stop_time_M=actual_stop,last_retained_time_M=times[-1],
                chi_activations=sum(int(x['chi_activations']) for x in floor),
                lapse_activations=sum(int(x['lapse_activations']) for x in floor),nonfinite=sum(int(x['nonfinite']) for x in floor),
                chi_min=min(float(x['chi_min']) for x in floor),lapse_min=min(float(x['lapse_min']) for x in floor),
                reader_messages_after_advance=log[first:].count('Read EMSTRUMPET'),
                initial_driver_cancellation=float(np.max(abs(q[0,:,138:140]))),recorded_max_family_speed=speed,
                peak_rss_bytes=r['peak_rss_bytes'],wall_s=r['wall_seconds']))
            assert audits[-1]['reader_messages_after_advance']==0,'static reader entered positive-time execution'
            assert audits[-1]['initial_driver_cancellation']==0,'initial driver prescription changed'
            c=round(336/(.875/2**level));punct=(cells[:,0]>=c-1)&(cells[:,0]<=c)&(cells[:,1]<=1)
            for ti,t in enumerate(times):
                for field,index in [('Gamma1',11),('Gamma2',12),('shift1',14),('shift2',15),('lapse',13)]:
                    native.append(dict(run=name,step=round(t/(.875/8192/4)),time_M=t,field=field,
                        puncture_peak_change=float(np.max(abs(q[ti,punct,index]-q[0,punct,index]))),
                        puncture_KO_peak=float(np.max(abs(q[ti,punct,56+index])))))
            cached[name]=(times,cells,np.array(q[0,:,:28]))
            for ray,v in a.NV.items():
                z=np.load(OUT/(name+'-'+ray+'-profiles.npz'))
                fs={n:{k:np.array([a.fields(frame,v)[k] for frame in z[f'p{n}']]) for k in a.fields(z[f'p{n}'][0],v)} for n in (8,10)}
                for ti,t in enumerate(times):
                    for field in ('Gamma','metric_Gamma','shift','lapse','C_Gamma','Ham','Mom','GaussE'):
                        disturbance=field in ('Gamma','metric_Gamma','shift','lapse')
                        u=fs[8][field][ti]-(fs[8][field][0] if disturbance else 0)
                        control=fs[10][field][ti]-(fs[10][field][0] if disturbance else 0)
                        for norm in ('peak','RMS'):
                            fn=(lambda x:float(np.max(abs(x)))) if norm=='peak' else a.rms
                            amplitude=fn(u);spread=fn(u-control)
                            history.append(dict(run=name,level=level,mode=mode,ray=ray,step=round(t/(.875/8192/4)),time_M=t,field=field,norm=norm,
                                amplitude=amplitude,interpolation_spread=spread,margin_5x=amplitude/(5*spread) if spread>0 else (float('inf') if amplitude>0 else 0)))
        old=cached[f'L{level}-sqrt'];new=cached[f'L{level}-geometric']
        assert np.array_equal(old[1],new[1])
        take=[k for k in range(28) if k!=13]
        count=old[2][:,take].size;bits=int(np.count_nonzero(old[2][:,take].view('u8')!=new[2][:,take].view('u8')));assert bits==0
        identity.append(dict(level=level,Float64_nonlapse_values=count,bit_mismatches=bits))
    lookup={(r['run'],r['ray'],r['field'],r['norm'],r['step']):r for r in history}
    for r in history:
        if r['mode']!='geometric' or r['time_M']==0:continue
        old=lookup[(f"L{r['level']}-sqrt",r['ray'],r['field'],r['norm'],r['step'])]
        ratios.append(dict(level=r['level'],ray=r['ray'],field=r['field'],norm=r['norm'],step=r['step'],time_M=r['time_M'],
            geometric_amplitude=r['amplitude'],sqrt_amplitude=old['amplitude'],geometric_over_sqrt=r['amplitude']/old['amplitude'] if old['amplitude']>0 else '',
            geometric_margin_5x=r['margin_5x'],sqrt_margin_5x=old['margin_5x'],
            endpoint=r['step']==74,qualified=r['margin_5x']>1 and old['margin_5x']>1))
    save('t16-launch-history.csv',history);save('t16-launch-ratios.csv',ratios);save('t16-launch-audits.csv',audits)
    save('t16-native-puncture.csv',native);save('t16-nonlapse-identity.csv',identity)
    endpoint=[r for r in ratios if r['endpoint']];save('t16-launch-endpoint.csv',endpoint)
    print('Endpoint disturbance ratios:',flush=True)
    for r in endpoint:
        if r['field'] in ('Gamma','metric_Gamma','shift','lapse'):print(r,flush=True)
def resources():
    rows=[]
    for f in sorted(ROOT.rglob('*.resources.json')):
        r=json.loads(f.read_text());peak=max(r['peak_rss_bytes'],max((v[0] for v in r.get('per_pid_peaks',{}).values()),default=0))
        rows.append(dict(process=r['process'],peak_rss_bytes=peak,tree_peak_bytes=r['tree_peak_bytes'],wall_s=r['wall_seconds'],
            returncode=r['returncode'],gate_reason=r['gate_reason'],record=f))
    save('t16-resources.csv',rows)
    targets=list(HERE.glob('t16-*'))+[HERE/'README.md']
    targets+=[ROOT/'launch.ex',ROOT/'audit.ex',ROOT/'replay.ex',ROOT/'census.ex',ROOT/'baseline/launch.ex']
    targets+=[HERE.parents[1]/'Source/InitialConditions/EMSBH'/f for f in ('EMSBHParams.hpp','EMSBH_trumpet_read.hpp','EMSBH_trumpet_read.impl.hpp')]
    targets+=[HERE.parents[1]/'Examples/EMS/SimulationParameters.hpp']
    targets+=[HERE/f for f in ('T11RHS.cpp','T13Replay.cpp','T14Launch.cpp','T14DenseTags.hpp','T14Census.cpp','t13-run.py','t13-analyze.py','t7-check.py')]
    targets+=list(ROOT.rglob('*.resources.json'))+list(ROOT.rglob('*.exit'))+list(ROOT.rglob('*.json'))
    targets+=list((ROOT/'analysis').glob('*'))
    targets+=[ROOT/'verification/pipeline-manifest.txt',ROOT/'verify.log']
    for marker in ROOT.rglob('done.exit'):
        if marker.read_text().strip()=='0':targets+=list(marker.parent.glob('*.xz'))+list(marker.parent.glob('*.log'))+list(marker.parent.glob('*.csv'))
    receipt=[f"T16 baseline 5da576b; no commit; OMP=2, BLAS=1; cap 3,000,000,000 bytes; measured peak RSS {max(r['peak_rss_bytes'] for r in rows)}; sampled tree peak {max(r['tree_peak_bytes'] for r in rows)} bytes.",
        'SHA256 bytes path']
    for f in sorted(set(targets)):
        if not f.is_file() or f==HERE/'t16-manifest.txt':continue
        digest=hashlib.sha256()
        with f.open('rb') as stream:
            for chunk in iter(lambda:stream.read(1<<20),b''):digest.update(chunk)
        receipt.append(f'{digest.hexdigest()} {f.stat().st_size} {f}')
    (HERE/'t16-manifest.txt').write_text('\n'.join(receipt)+'\n')
if __name__=='__main__':
    if sys.argv[1]=='check':
        # Actual native t=0 capture, filtered and replayed; no evolution.
        name='census-check-L13';d=ROOT/'evolution'/name;d.mkdir(exist_ok=True)
        raw=OUT/'census-check.raw';raw.write_bytes(lzma.open(ROOT/'census/L13/t13-t7-stage-L13.xz','rb').read())
        with raw.open('rb') as src,(d/'t13-t7-stage-L13.xz').open('wb') as dst:
            subprocess.run([PY,str(HERE/'t16-xz.py')],stdin=src,stdout=dst,check=True)
        raw.unlink();(d/'done.exit').write_text('0\n');extract(name,13)
        times=np.load(OUT/(name+'-meta.npz'))['times'];assert times.tolist()==[0.]
        save('t16-analysis-check.csv',[dict(input='actual L13 native census t=0',snapshots=1,
            interpolation='I8/I10',coupling_f2=-.9,native_RHS_bit_mismatches=0,result='PASS')])
    if sys.argv[1]=='analyze':
        for l in (13,14):
            for mode in ('sqrt','geometric'):extract(f'L{l}-{mode}',l)
        measure()
    resources()
