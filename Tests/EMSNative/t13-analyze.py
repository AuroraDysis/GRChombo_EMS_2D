#!/usr/bin/env python3
"""Bounded current-state T13 analysis; no static solution is opened."""
import sys
sys.dont_write_bytecode=True
import csv,hashlib,importlib.util,json,math,subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;TMP=Path('/private/tmp/ems-t13')
OUT=TMP/'analysis';OUT.mkdir(exist_ok=True)
PY='/Users/auroradysis/miniconda3/bin/python';LEGS=('original','maximal','original-half','maximal-half')
N=28;COL=180;H=7/12/4096;DT=H/4;CENTRE=round(336/H)
R=np.unique(np.r_[np.linspace(.00075,.0025,257),.0015]);NV={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/np.sqrt(2)}
ODD_FIELDS=[2,7,12,15,17,22,23,25,26]
ODD=[o+c for o in (0,28,56,84,152) for c in ODD_FIELDS]+[113+2*g for g in range(13)]+[139,142,146,149]
SELECT=list(range(28))+list(range(145,152));ODD_C=ODD_FIELDS+[29,32]
spec=importlib.util.spec_from_file_location('decoder',HERE/'t7-check.py')
decoder=importlib.util.module_from_spec(spec);spec.loader.exec_module(decoder)

def save(name,rows):
    assert rows,name
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def timed(label,cmd):
    subprocess.run([PY,str(HERE/'t13-run.py'),'--measure',label,'--directory',str(OUT),'--',*cmd],check=True)

def dense(group):
    iv=np.concatenate([a['cells'] for a in group]);lo=iv.min(0);hi=iv.max(0)
    a=np.full((hi[1]-lo[1]+1,hi[0]-lo[0]+1,N),np.nan)
    # Ghosts first; native valid cells have priority at same-level box seams.
    for valid_only in (False,True):
        for fr in group:
            p=fr['cells'];x0,y0,x1,y1=fr['valid']
            take=(p[:,0]>=x0)&(p[:,0]<=x1)&(p[:,1]>=y0)&(p[:,1]<=y1) if valid_only else np.ones(len(p),bool)
            a[p[take,1]-lo[1],p[take,0]-lo[0]]=fr['values'][take]
    return lo,a

def selected():
    # All tensor P8 nodes on both W rays, native corroboration, and puncture cells.
    cells=set()
    for vec in NV.values():
        xy=R[:,None]*vec
        z=np.c_[(336+xy[:,0])/H-.5,xy[:,1]/H-.5]
        start=np.floor(z).astype(int)-3
        for s in start:
            for j in range(8):
                for i in range(8):cells.add((s[0]+i,abs(s[1]+j) if s[1]+j>=0 else -s[1]-j-1))
    for j in range(6):
        for i in range(-5,6):cells.add((CENTRE+i,j))
    return np.array(sorted(cells,key=lambda p:(p[1],p[0])),int)

def extract(leg):
    directory=TMP/'evolution'/leg;assert (directory/'done.exit').read_text().strip()=='0'
    cells=selected();times=[];initial=None;group=[];previous=None
    payload=OUT/(leg+'-input.bin');output=OUT/(leg+'-q.bin')
    with payload.open('wb') as stream:
        def consume(group):
            nonlocal initial
            t=group[0]['meta'][0];lo,a=dense(group)
            if t==0:initial=(lo,a.copy());np.savez(OUT/(leg+'-initial.npz'),lo=lo,a=a)
            record=np.empty((len(cells),3+49*N));record[:,0]=H
            record[:,1]=(cells[:,1]+.5)*H;record[:,2]=(cells[:,0]+.5)*H-336
            for j in range(-3,4):
                for i in range(-3,4):
                    v=a[cells[:,1]+j-lo[1],cells[:,0]+i-lo[0]]
                    assert np.isfinite(v).all(),(leg,t,i,j)
                    record[:,3+(j+3)*7+i+3::49]=v
            record.tofile(stream);times.append(t)
        for fr in decoder.frames(directory/'t13-t7-stage-L12.xz'):
            if fr['phase']!=50:continue
            t=fr['meta'][0]
            if previous is not None and t!=previous:consume(group);group=[]
            group.append(fr);previous=t
        if group:consume(group)
    timed('native-snapshots-'+leg,[str(TMP/'replay.ex'),'--state',str(payload),str(output)])
    q=np.memmap(output,mode='r',dtype='f8',shape=(len(times),len(cells),COL))
    assert np.isfinite(q).all() and q[:,:,144].max()==0
    np.savez(OUT/(leg+'-meta.npz'),cells=cells,times=np.array(times))
    digest=hashlib.sha256()
    with payload.open('rb') as f:
        for a in iter(lambda:f.read(1<<20),b''):digest.update(a)
    (OUT/(leg+'-input.sha256')).write_text(digest.hexdigest()+' '+str(payload)+'\n')
    payload.unlink() # Own regenerable intermediate, hash retained.
    print(leg,'snapshots',len(times),'native cells',len(cells),'q bytes',output.stat().st_size,flush=True)

def weights(z,start,n):
    return np.array([np.prod([(z-start-j)/(i-j) for j in range(n) if i!=j],axis=0) for i in range(n)])

def sample(cells,q,vec,n):
    lookup={tuple(k):i for i,k in enumerate(cells)}
    xy=R[:,None]*vec;z=np.c_[(336+xy[:,0])/H-.5,xy[:,1]/H-.5]
    start=np.floor(z).astype(int)-n//2+1
    wx=weights(z[:,0],start[:,0],n);wy=weights(z[:,1],start[:,1],n)
    ans=np.zeros((len(R),len(SELECT)));anchor=None
    for j in range(n):
        yy=start[:,1]+j;iy=np.where(yy<0,-yy-1,yy)
        for i in range(n):
            indices=np.array([lookup[(x,y)] for x,y in zip(start[:,0]+i,iy)])
            v=np.array(q[np.ix_(indices,SELECT)])
            v[np.ix_(np.flatnonzero(yy<0),ODD_C)]*=-1
            if anchor is None:anchor=v.copy();ans=anchor.copy()
            ans+=(v-anchor)*(wx[i]*wy[j])[:,None]
    return ans

def fields(a,vec):
    return dict(Gamma=a[:,11:13]@vec,metric_Gamma=a[:,28:30]@vec,
        shift=a[:,14:16]@vec,lapse=a[:,13],C_Gamma=(a[:,11:13]-a[:,28:30])@vec,
        Ham=a[:,30],Mom=np.hypot(a[:,31],a[:,32]),Theta=a[:,10],K=a[:,5],chi=a[:,0],
        B=a[:,16:18]@vec,GaussE=a[:,33])

def rms(a):return float(np.sqrt(np.trapezoid(a*a,R)/(R[-1]-R[0])))

def measure():
    history=[];profiles=[];audits=[];cached={};initials={}
    for leg in LEGS:
        d=TMP/'evolution'/leg;meta=np.load(OUT/(leg+'-meta.npz'));cells=meta['cells'];times=meta['times']
        q=np.memmap(OUT/(leg+'-q.bin'),mode='r',dtype='f8',shape=(len(times),len(cells),COL))
        initial=np.load(OUT/(leg+'-initial.npz'));initials[leg]=(initial['lo'],initial['a'])
        pp=(HERE/'params'/('t13-E-mid-'+leg+'.txt')).read_text()
        assert 'sigma = 1\n' in pp and 'amr_transfer = point' in pp and 't2_guard_initial_data_after_t0 = true' in pp
        floor_rows=[]
        for f in d.glob('t13-floors-*.csv'):
            floor_rows.extend(csv.DictReader(f.open()))
        log=(d/'run.log').read_text();first=log.index('GRAMRLevel::advance level 0 at time 0')
        det=q[:,:,1]*q[:,:,3]-q[:,:,2]**2
        inv_max=(q[:,:,1]+q[:,:,3]+np.sqrt((q[:,:,1]-q[:,:,3])**2+4*q[:,:,2]**2))/(2*det)
        gauge_speed=float(np.max(np.sqrt(np.maximum(inv_max,1/q[:,:,4]))+np.hypot(q[:,:,14],q[:,:,15])))
        audits.append(dict(run=leg,exit=0,snapshots=len(times),end_time_M=times[-1],
            chi_activations=sum(int(a['chi_activations']) for a in floor_rows),
            lapse_activations=sum(int(a['lapse_activations']) for a in floor_rows),
            nonfinite_count=sum(int(a['nonfinite']) for a in floor_rows),
            chi_min=min(float(a['chi_min']) for a in floor_rows),lapse_min=min(float(a['lapse_min']) for a in floor_rows),
            reader_messages_after_first_advance=log[first:].count('Read EMSTRUMPET'),
            current_shift_speed_proxy_max=gauge_speed,driver_initial_absmax=float(np.max(abs(q[0,:,138:140]))),
            peak_rss_GB=json.loads((d/(leg+'.resources.json')).read_text())['peak_rss_bytes']/1e9,
            wall_s=json.loads((d/(leg+'.resources.json')).read_text())['wall_seconds']))
        for ray,vec in NV.items():
            samples={n:np.array([sample(cells,frame,vec,n) for frame in q]) for n in (6,8)}
            np.savez_compressed(OUT/(leg+'-'+ray+'-profiles.npz'),times=times,r=R,p6=samples[6],p8=samples[8])
            arrays={n:{k:np.array([fields(a,vec)[k] for a in samples[n]]) for k in fields(samples[n][0],vec)} for n in (6,8)}
            cached[leg,ray]=(times,arrays)
            for ti,t in enumerate(times):
                f=arrays[6];f8=arrays[8];gamma=f['Gamma'][ti]-f['Gamma'][0];gamma8=f8['Gamma'][ti]-f8['Gamma'][0]
                row=dict(run=leg,ray=ray,time_M=t,step=ti,Gamma_peak=float(max(abs(gamma))),Gamma_RMS=rms(gamma),
                    interpolation_peak=float(max(abs(gamma-gamma8))),interpolation_RMS=rms(gamma-gamma8),
                    probe_Gamma=float(np.interp(.0015,R,gamma)),probe_C_Gamma=float(np.interp(.0015,R,f['C_Gamma'][ti])),
                    metric_Gamma_peak=float(max(abs(f['metric_Gamma'][ti]-f['metric_Gamma'][0]))),
                    C_Gamma_RMS=rms(f['C_Gamma'][ti]),C_Gamma_peak=float(max(abs(f['C_Gamma'][ti]))),
                    Ham_RMS=rms(f['Ham'][ti]),Ham_peak=float(max(abs(f['Ham'][ti]))),
                    Mom_RMS=rms(f['Mom'][ti]),Mom_peak=float(max(abs(f['Mom'][ti]))),
                    **{'probe_'+k:float(np.interp(.0015,R,f[k][ti])) for k in
                       ('metric_Gamma','shift','lapse','Ham','Mom','Theta','K','chi','B','GaussE')})
                history.append(row)
                if ti in {0,round(len(times)*.25),round((len(times)-1)*.5),round((len(times)-1)*.75),len(times)-1}:
                    for j,r in enumerate(R):
                        profiles.append(dict(run=leg,ray=ray,time_M=t,R_M=r,
                            **{k:float(v[ti,j]) for k,v in f.items()},Gamma_disturbance=float(gamma[j]),
                            Gamma_interpolation_spread=float(gamma[j]-gamma8[j])))
    identity=[]
    for a,b in [('original','maximal'),('original-half','maximal-half')]:
        lo,x=initials[a];lb,y=initials[b];assert np.array_equal(lo,lb) and x.shape==y.shape
        use=np.isfinite(x).all(-1)&np.isfinite(y).all(-1);keep=[c for c in range(N) if c!=13]
        xx=x[use][:,keep].copy();yy=y[use][:,keep].copy()
        count=int(np.count_nonzero(xx.view('u8')!=yy.view('u8')))
        identity.append(dict(pair=a+'/'+b,non_lapse_values=xx.size,bit_mismatches=count));assert count==0
    save('t13-run-audit.csv',audits);save('t13-initial-identity.csv',identity)
    save('t13-history.csv',history);save('t13-profiles.csv',profiles)
    screen=[];temporal=[]
    for ray in NV:
        for policy,aa,bb in [('nominal','original','maximal'),('half','original-half','maximal-half')]:
            for norm,col in [('peak','Gamma_peak'),('RMS','Gamma_RMS')]:
                old=[a for a in history if a['run']==aa and a['ray']==ray];new=[a for a in history if a['run']==bb and a['ray']==ray]
                ratios=[b[col]/a[col] for a,b in zip(old[1:],new[1:])]
                spread=sum(a['interpolation_'+norm] for a in (old[-1],new[-1]))
                tempor=0.;temporal_history=np.zeros(len(old))
                for run in (aa,bb):
                    run=run.removesuffix('-half')
                    control=[a for a in history if a['run']==run+'-half' and a['ray']==ray]
                    base=[a for a in history if a['run']==run and a['ray']==ray]
                    deltas=[abs(a[col]-b[col]) for a,b in zip(base,control[::2])]
                    native=cached[run,ray][1][6]['Gamma'];half=cached[run+'-half',ray][1][6]['Gamma'][::2]
                    differences=(native-native[0])-(half-half[0])
                    profile_deltas=[float(max(abs(z))) if norm=='peak' else rms(z) for z in differences]
                    tempor+=profile_deltas[-1]
                    temporal_history+=np.array(profile_deltas) if policy=='nominal' else np.interp(cached[aa,ray][0],cached[run,ray][0],profile_deltas)
                    if policy=='nominal':
                        temporal.append(dict(run=run,ray=ray,norm=norm,endpoint_nominal=base[-1][col],
                            endpoint_half=control[-1][col],relative_change=control[-1][col]/base[-1][col]-1,
                            max_absolute_history_change=max(deltas),endpoint_profile_difference=profile_deltas[-1],
                            max_history_profile_difference=max(profile_deltas)))
                reduction=old[-1][col]-new[-1][col]
                reductions=np.array([x[col]-y[col] for x,y in zip(old[1:],new[1:])])
                spreads=np.array([x['interpolation_'+norm]+y['interpolation_'+norm] for x,y in zip(old[1:],new[1:])])
                common=slice(None) if policy=='nominal' else slice(1,None,2)
                screen.append(dict(policy=policy,ray=ray,norm=norm,original=old[-1][col],maximal=new[-1][col],
                    endpoint_ratio=new[-1][col]/old[-1][col],history_sup_ratio=max(b[col] for b in new)/max(a[col] for a in old),
                    worst_nonzero_time_ratio=max(ratios),interpolation_bound=spread,temporal_bound=tempor,
                    reduction_over_interpolation=reduction/spread if spread else float('inf'),
                    reduction_over_temporal=reduction/tempor if tempor else float('inf'),
                    minimum_history_reduction_over_interpolation=float(np.min(reductions/np.maximum(spreads,1e-300))),
                    minimum_history_reduction_over_temporal=float(np.min(reductions[common]/np.maximum(temporal_history[1:][common],1e-300)))))
    save('t13-screen.csv',screen);save('t13-temporal.csv',temporal)
    print('screen',json.dumps(screen,indent=2),flush=True)

def qualify():
    histories=list(csv.DictReader((HERE/'t13-history.csv').open()));rows=[];native=[]
    for policy,a,b in [('nominal','original','maximal'),('half','original-half','maximal-half')]:
        for ray,vec in NV.items():
            x=np.load(OUT/(a+'-'+ray+'-profiles.npz'));y=np.load(OUT/(b+'-'+ray+'-profiles.npz'))
            xx={n:{k:np.array([fields(v,vec)[k] for v in x[f'p{n}']]) for k in fields(x['p6'][0],vec)} for n in (6,8)}
            yy={n:{k:np.array([fields(v,vec)[k] for v in y[f'p{n}']]) for k in fields(y['p6'][0],vec)} for n in (6,8)}
            for constraint in ('Ham','Mom','C_Gamma','GaussE'):
                for norm in ('peak','RMS'):
                    calc=lambda v:float(np.max(abs(v))) if norm=='peak' else rms(v)
                    va=np.array([calc(v) for v in xx[6][constraint]]);vb=np.array([calc(v) for v in yy[6][constraint]])
                    ea=np.array([calc(z-w) for z,w in zip(xx[6][constraint],xx[8][constraint])])
                    eb=np.array([calc(z-w) for z,w in zip(yy[6][constraint],yy[8][constraint])])
                    rows.append(dict(policy=policy,ray=ray,constraint=constraint,norm=norm,
                        endpoint_original=va[-1],endpoint_maximal=vb[-1],endpoint_ratio=vb[-1]/va[-1],
                        history_max_original=va.max(),history_max_maximal=vb.max(),history_ratio=vb.max()/va.max(),
                        endpoint_interpolation_bound=ea[-1]+eb[-1],
                        worst_excess_over_joint_interpolation=float(np.max((vb-va)/(ea+eb+1e-300))),
                        times_significantly_worse=int(np.count_nonzero(vb>va+ea+eb))))
    for leg in LEGS:
        meta=np.load(OUT/(leg+'-meta.npz'));cells=meta['cells'];times=meta['times']
        q=np.memmap(OUT/(leg+'-q.bin'),mode='r',dtype='f8',shape=(len(times),len(cells),COL))
        x=(cells[:,0]+.5)*H-336;y=(cells[:,1]+.5)*H;r=np.hypot(x,y)
        for ray,vec in NV.items():
            mask=((cells[:,1]==0)&(x>0)) if ray=='axis' else (cells[:,0]-CENTRE==cells[:,1])
            mask&=(r>=R[0])&(r<=R[-1]);delta=(q[-1,mask,11:13]-q[0,mask,11:13])@vec
            native.append(dict(run=leg,ray=ray,native_cells=int(mask.sum()),
                Gamma_peak=float(abs(delta).max()),Gamma_RMS=float(np.sqrt(np.mean(delta*delta))),
                condition='axis first-row proxy y=h/2; diagonal native centres; unweighted cell RMS',
                GaussB_nonzero=int(np.count_nonzero(q[:,:,151]))))
    save('t13-constraint-screen.csv',rows);save('t13-native-screen.csv',native)
    print(json.dumps(rows,indent=2),flush=True)

def stages():
    s=importlib.util.spec_from_file_location('closure',HERE/'t13-check.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    rows=[]
    for leg in LEGS:rows.append(m.closure(TMP/'evolution'/leg,12))
    save('t13-stage-closure.csv',rows)

if __name__=='__main__':
    if sys.argv[1]=='extract':
        for leg in LEGS:extract(leg)
    elif sys.argv[1]=='measure':measure()
    elif sys.argv[1]=='stages':stages()
    elif sys.argv[1]=='qualify':qualify()
    else:raise ValueError(sys.argv[1])
