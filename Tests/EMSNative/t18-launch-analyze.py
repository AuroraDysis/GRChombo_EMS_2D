#!/usr/bin/env python3
"""Current-field constraint replay and matched-clock temporal/spatial checks."""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
import importlib.util
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
w=module('launchwork',HERE/'t18-launch-work.py')
a=module('nativeanalysis',HERE/'t13-analyze.py')
t17=module('t17analysis',w.w.T17/'t17-analyze.py')
ROOT=w.ROOT/'launch';OUT=ROOT/'analysis';OUT.mkdir(exist_ok=True)
CLOCKS=np.array([0.,.000640869140625,.00128173828125,.001922607421875])
FIELDS=('Gamma','lapse','shift','C_Gamma','Ham','Mom','GaussE','GaussB','Theta','Lambda','Xi')
def save(name,rows):
    assert rows
    with (HERE/name).open('w',newline='') as f:
        out=csv.DictWriter(f,fieldnames=rows[0]);out.writeheader();out.writerows(rows)
def norm(v,kind):return float(np.max(abs(v))) if kind=='peak' else a.rms(v)
def extract(name,hole):
    d=ROOT/name
    r=json.loads((d/(name+'.resources.json')).read_text())
    assert (d/'done.exit').read_text().strip()=='0'
    assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason']
    assert r['peak_rss_bytes']<4e9
    top=12 if name.startswith('coarse') else 13;h=1.75/2**top
    center=432. if hole=='left' else 464.
    vectors={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/np.sqrt(2)}
    if hole=='right':vectors={k:v*np.array([-1.,1.]) for k,v in vectors.items()}
    cells=t17.select(h,center,vectors)
    payload=OUT/(name+'-'+hole+'-input.bin');output=OUT/(name+'-'+hole+'-q.bin')
    times=[];group=[];previous=None
    with payload.open('wb') as stream:
        def consume(frames):
            lo,u=a.dense(frames);rec=np.empty((len(cells),3+49*28))
            rec[:,0]=h;rec[:,1]=(cells[:,1]+.5)*h;rec[:,2]=(cells[:,0]+.5)*h-center
            for j in range(-3,4):
                for i in range(-3,4):
                    v=u[cells[:,1]+j-lo[1],cells[:,0]+i-lo[0]]
                    assert np.isfinite(v).all()
                    rec[:,3+(j+3)*7+i+3::49]=v
            rec.tofile(stream);times.append(frames[0]['meta'][0])
        for fr in a.decoder.frames(d/f't18-launch-t7-stage-L{top}.xz'):
            if fr['phase']!=50 or (np.mean((fr['cells'][:,0]+.5)*h)<448)!=(hole=='left'):continue
            t=fr['meta'][0]
            if not np.any(abs(CLOCKS-t)<1e-14):continue
            if previous is not None and t!=previous:consume(group);group=[]
            group.append(fr);previous=t
        if group:consume(group)
    assert np.array_equal(times,CLOCKS),(name,hole,times)
    w.w.measured('replay-'+name+'-'+hole,[w.ROOT/'launch-replay.ex','--state',payload,output])
    (OUT/(name+'-'+hole+'-input.sha256')).write_text(hashlib.sha256(payload.read_bytes()).hexdigest()+'\n')
    payload.unlink()
    q=np.memmap(output,mode='r',dtype='f8',shape=(4,len(cells),180))
    assert np.isfinite(q).all() and np.count_nonzero(q[:,:,144])==0
    control=None
    if name.endswith('dt0.25'):
        rung=name.split('-')[0];oldroot=Path('/private/tmp/ems-t17/analysis');oldname=rung+'-geometric-'+hole
        meta=np.load(oldroot/(oldname+'-meta.npz'))
        old=np.memmap(oldroot/(oldname+'-q.bin'),mode='r',dtype='f8',shape=(len(meta['times']),len(meta['cells']),180))
        lookup={tuple(iv):i for i,iv in enumerate(meta['cells'])}
        shift=1024*2**top
        pairs=[(i,lookup[(int(iv[0])+shift,int(iv[1]))]) for i,iv in enumerate(cells) if (int(iv[0])+shift,int(iv[1])) in lookup]
        assert pairs
        ii,jj=map(np.array,zip(*pairs));compared=different=0;maximum=0.
        for k,t in enumerate(CLOCKS):
            kk=int(np.flatnonzero(abs(meta['times']-t)<1e-14)[0])
            x=np.ascontiguousarray(q[k,ii,:28]);y=np.ascontiguousarray(old[kk,jj,:28])
            compared+=x.size;different+=int(np.count_nonzero(x.view('u8')!=y.view('u8')))
            maximum=max(maximum,float(np.max(abs(x-y))))
        control=dict(run=name,hole=hole,matched_native_cells=len(ii),values=compared,
            differing_bits=different,max_absolute_difference=maximum,reference='full T17 production-grid geometric launch')
    profiles={}
    for ray,vec in vectors.items():
        profiles[ray]={}
        for order in (8,10):
            p=np.array([t17.sample(cells,frame,vec,order,h,center) for frame in q])
            fs={k:np.array([a.fields(frame,vec)[k] for frame in p]) for k in a.fields(p[0],vec)}
            fs.update(GaussB=p[:,:,34],Lambda=p[:,:,20],Xi=p[:,:,27])
            profiles[ray][order]=fs
        np.savez_compressed(OUT/(name+'-'+hole+'-'+ray+'.npz'),r=a.R,times=CLOCKS,
            **{str(order)+'_'+field:profiles[ray][order][field] for order in (8,10) for field in FIELDS})
    return profiles,control
def main():
    cached={};controls=[]
    for name in ('coarse-dt0.25','coarse-dt0.375','coarse-dt0.5','fine-dt0.25'):
        for hole in ('left','right'):
            cached[name,hole],control=extract(name,hole)
            if control:controls.append(control)
    history=[];errors=[]
    for hole in ('left','right'):
        for ray in ('axis','diagonal'):
            for field in FIELDS:
                for state in ('launch_change','instantaneous_state'):
                    def profile(name,order):
                        v=cached[name,hole][ray][order][field]
                        return v-v[0] if state=='launch_change' else v
                    base=profile('coarse-dt0.25',8);fine=profile('fine-dt0.25',8)
                    eb=base-profile('coarse-dt0.25',10);ef=fine-profile('fine-dt0.25',10)
                    for dt in (.375,.5):
                        name='coarse-dt'+str(dt);p=profile(name,8);ep=p-profile(name,10)
                        for step,t in enumerate(CLOCKS[1:],1):
                            for kind in ('peak','RMS'):
                                et=norm(p[step]-base[step],kind);eh=norm(base[step]-fine[step],kind)
                                # Interpolation is linear: measure its spread on
                                # the matched difference, retaining cancellation.
                                ut=5*norm(ep[step]-eb[step],kind)
                                uh=5*norm(eb[step]-ef[step],kind)
                                scale=max(norm(base[step],kind),norm(p[step],kind),norm(fine[step],kind),1.)
                                floor=128*np.finfo(float).eps*scale
                                if eh<=max(uh,floor):verdict='spatial difference unresolved'
                                elif et+ut<eh-uh:verdict='resolved subdominant'
                                elif et<eh:verdict='raw subdominant; interpolation spread overlaps'
                                else:verdict='not subdominant'
                                errors.append(dict(dt_multiplier=dt,hole=hole,ray=ray,field=field,measure=state,
                                    time_M=t,norm=kind,temporal_difference=et,spatial_difference=eh,
                                    temporal_to_spatial=et/eh if eh else 'undefined',
                                    temporal_spread_5x=ut,spatial_spread_5x=uh,roundoff_scale=floor,verdict=verdict))
                    for name in ('coarse-dt0.25','coarse-dt0.375','coarse-dt0.5','fine-dt0.25'):
                        p=profile(name,8);ep=p-profile(name,10)
                        for step,t in enumerate(CLOCKS):
                            for kind in ('peak','RMS'):
                                history.append(dict(run=name,hole=hole,ray=ray,field=field,measure=state,time_M=t,
                                    norm=kind,amplitude=norm(p[step],kind),spread=norm(ep[step],kind)))
    save('t18-launch-full-grid-control.csv',controls)
    save('t18-launch-constraints.csv',history)
    save('t18-temporal-spatial.csv',errors)
    summary=[]
    for dt in (.375,.5):
        for field in FIELDS:
            rs=[r for r in errors if r['dt_multiplier']==dt and r['field']==field and r['measure']=='launch_change']
            summary.append(dict(dt_multiplier=dt,field=field,comparisons=len(rs),
                raw_ratio_max=max((r['temporal_to_spatial'] for r in rs if isinstance(r['temporal_to_spatial'],float)),default='undefined'),
                resolved_subdominant=sum(r['verdict']=='resolved subdominant' for r in rs),
                raw_subdominant_spread_overlaps=sum(r['verdict'].startswith('raw subdominant') for r in rs),
                spatial_unresolved=sum(r['verdict']=='spatial difference unresolved' for r in rs),
                not_subdominant=sum(r['verdict']=='not subdominant' for r in rs)))
    save('t18-temporal-spatial-summary.csv',summary)
    print(json.dumps(dict(full_grid_controls=controls,temporal_spatial=summary),indent=2))
if __name__=='__main__':main()
