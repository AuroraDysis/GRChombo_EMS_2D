#!/usr/bin/env python3
"""Read-only exp-0020 admission diagnostics; no static reader or evolution path."""
import csv, gzip, hashlib, importlib.util, json, math, re, resource, sys, time
from collections import Counter, defaultdict
from pathlib import Path
sys.dont_write_bytecode=True
import h5py
import numpy as np
from scipy.signal import find_peaks, peak_widths

HERE=Path(__file__).resolve().parent
DATA=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0020')
SUB=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0020/submissions/exp-0020')
TMP=Path('/private/tmp/ems-t8');TMP.mkdir(exist_ok=True)
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
audit=module('t8audit',SUB/'audit.py')
t6=module('t8t6',HERE/'t6c-analyze.py')
LEGS=('E-low','E-mid','E-high');TIMES=np.arange(13)*.875
F=tuple(audit.EVOLVED)+('Ham','Mom1','Mom2','Mom','GaussE','GaussB')
INDEX={f:i for i,f in enumerate(F)}
ODD=[INDEX[f] for f in ('h12','A12','Gamma2','shift2','B2','By','Bz','Ey','Ez','Mom2')]
VECS=np.array(((1,0),(0,1),(2**-.5,2**-.5)));RAYS=('axis','equator','corner')
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def save(name,rows,directory=HERE):
    assert rows,name
    with (directory/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def params(leg):
    return {k.strip():v.split('#',1)[0].strip() for k,v in (l.split('=',1) for l in (SUB/f'params-{leg}.txt').read_text().splitlines() if '=' in l and not l.startswith('#'))}
def common_time(t):
    q=round(float(t)/.875)*.875;assert abs(float(t)-q)<1e-10;return q

def tables():
    inputs=[];runrows=[];rates=[];plots=[];horizons=[];native=[];mq=[]
    for leg in LEGS:
        d=DATA/leg;p=params(leg);steps=[];failures=[];read_before=read_after=0;firststep=False
        assert p['sigma']=='1' and p['amr_transfer']=='point' and p['nan_check']=='1' and p['t2_guard_initial_data_after_t0']=='true'
        assert set(p['regrid_interval'].split())=={'0'}
        for name in ('run.exit','evolution.exit','reduction.exit'):assert (d/name).read_text().strip()=='0'
        for name in ('run.done','chain.done','qualification.done'):assert (d/name).is_file()
        runlog=(d/'run.log').read_text();assert 'number_procs = 32' in runlog and 'threads = 4' in runlog
        with gzip.open(d/'pout.0.gz','rt') as f:
            for line in f:
                m=re.search(r'coarse time step\s+(\d+)\s+old time = (\S+)\s+old dt = (\S+)\s+wallclocktime = (\S+)',line)
                if m:
                    firststep=True;steps.append(dict(leg=leg,step=int(m[1]),old_time_M=float(m[2]),dt_M=float(m[3]),seconds=float(m[4])))
                if 'Read EMSTRUMPET' in line:
                    if firststep:read_after+=1
                    else:read_before+=1
                if any(s in line.lower() for s in ('nancheck in','mayday::error','not found','segmentation fault')):failures.append(line.strip())
        assert not failures and len(steps)==int(p['max_steps']) and read_after==0 and read_before==91
        rates+=steps;steady=np.array([r['seconds'] for r in steps[1:]])
        paths=sorted((d/'plt').glob('*.hdf5'));assert len(paths)==13
        for path in paths:
            with h5py.File(path) as f:
                t=common_time(f['level_0'].attrs['time']);assert audit.names(f)==list(F) and int(f.attrs['num_levels'])==13
                for lev in range(13):
                    g=f[f'level_{lev}'];assert g['data:datatype=0'].dtype==np.dtype('f8')
                    assert abs(float(g.attrs['time'])-t)<1e-10
                    bs=g['boxes'][:];h=float(g.attrs['dx']);lo=np.array([[int(b[k]) for k in ('lo_i','lo_j')] for b in bs]).min(0);hi=np.array([[int(b[k]) for k in ('hi_i','hi_j')] for b in bs]).max(0)
                    r=336 if lev==0 else 112/2**lev
                    assert np.max(abs(np.array([lo[0]*h-336,(hi[0]+1)*h-336,lo[1]*h,(hi[1]+1)*h])-[-r,r,0,r]))<1e-10
                plots.append(dict(leg=leg,time_M=t,file=str(path),bytes=path.stat().st_size,components=34,levels=13,geometry='REGISTERED_COMMON_FACES',checksum_provenance='CONTROLLER_RCLONE_CHECK_ZERO_DIFFERENCES'))
        assert {q['time_M'] for q in plots if q['leg']==leg}==set(TIMES)
        row=json.loads((d/'run.done').read_text());assert row['checkpoint_step']==int(p['max_steps']) and abs(row['checkpoint_time_M']-10.5)<256*np.finfo(float).eps*10.5
        runrows.append(dict(leg=leg,exit=0,endpoint_M=row['checkpoint_time_M'],steps=len(steps),sigma=1,point=True,fixed_hierarchy=True,
            missing_or_failure_lines=0,plots=13,qualified_rows=26,steady_median_s=float(np.median(steady)),steady_min_s=float(steady.min()),
            steady_max_s=float(steady.max()),sum_coarse_step_s=sum(r['seconds'] for r in steps),params_sha256=digest(SUB/f'params-{leg}.txt'),
            MPI_ranks=32,OMP_threads=4,startup_reader_messages=read_before,reader_messages_after_first_step=read_after,
            sigma_evidence='SEALED_PARAMETER_FILE;VALUE_NOT_ECHOED_BY_VERBOSITY_1',post_t0_reader_guard=True))
        hh=read(d/'qualified-horizons.csv');assert len(hh)==26
        for r in hh:
            r['time_M']=common_time(r['time_M']);r['N_theta']=int(r['N_theta']);assert r['status']=='FOUND' and int(r['stage'])==2 and float(r['expansion_squared'])<=1e-12 and float(r['theta_minus'])<0
        for n in (48,96):assert {r['time_M'] for r in hh if r['N_theta']==n}==set(TIMES)
        h96={r['time_M']:r for r in hh if r['N_theta']==96};h48={r['time_M']:r for r in hh if r['N_theta']==48}
        for t,r in h96.items():
            out=dict(leg=leg,time_M=t,A=float(r['A']),Q=float(r['Q']),expansion_RMS=math.sqrt(float(r['expansion_squared'])))
            for field in ('A','Q'):
                z=float(r[field]);base=float(h96[0.][field]);drift=z/base-1
                other=float(h48[t][field])/float(h48[0.][field])-1
                out.update({f'relative_{field}':drift,f'angular_absolute_{field}':float(r[f'angular_delta_{field}']),
                    f'angular_drift_{field}':abs(drift-other),f'stopping_absolute_{field}':float(r[f'stopping_delta_{field}']),
                    f'stopping_drift_bound_{field}':(float(r[f'stopping_delta_{field}'])+abs(z/base)*float(h96[0.][f'stopping_delta_{field}']))/base})
            horizons.append(out)
        for r in read(d/'finder-values.csv'):
            t=float(r['t']);a=float(r['Area']);q=float(r['Q_charge']);native.append(dict(leg=leg,time_M=t,A=a,Q=q,
                relative_A=a/float(h96[0.]['A'])-1,relative_Q=q/float(h96[0.]['Q'])-1,expansion_RMS=math.sqrt(float(r['err'])),mode=r['mode']))
        for r in read(d/'mq-values.csv'):mq.append(dict(leg=leg,**{k:float(v) for k,v in r.items()}))
    save('t8-run-audit.csv',runrows);save('t8-coarse-step-times.csv',rates);save('t8-plot-census.csv',plots)
    save('t8-qualified-horizons.csv',horizons);save('t8-native-horizon-history.csv',native);save('t8-mq-charge.csv',mq)
    s=read(DATA/'self-differences.csv');assert len(s)==37*11*13
    fields=[]
    for r in s:
        p=float(r['measured_order']) if r['measured_order'] else math.nan
        d=min(float(r['difference_low_mid']),float(r['difference_mid_high']));spread=float(r['interpolation_spread_4_vs_6'])
        status='UNDEFINED_ZERO_OR_FLOOR' if not math.isfinite(p) else 'NEGATIVE' if p<0 else 'LOW_LT_3' if p<3 else 'GE_3'
        fields.append({**r,'flag':status,'interpolation_spread_over_smaller_difference':spread/d if d else math.nan})
    save('t8-field-orders.csv',fields)
    groups=defaultdict(list)
    for r in fields:groups[(r['mask'],r['field'])].append(r)
    summaries=[]
    for (mask,field),rr in sorted(groups.items()):
        measured=[float(r['measured_order']) for r in rr if r['measured_order']]
        summaries.append(dict(mask=mask,field=field,measured_times=len(measured),undefined_times=13-len(measured),
             minimum_order=min(measured) if measured else '',maximum_order=max(measured) if measured else '',
             negative_times=sum(p<0 for p in measured),low_positive_times=sum(0<=p<3 for p in measured),
             maximum_interpolation_fraction=max((float(r['interpolation_spread_over_smaller_difference']) for r in rr if r['measured_order']),default=math.nan)))
    save('t8-field-order-summary.csv',summaries)
    norms={leg:{(common_time(r['time_M']),r['mask'],r['constraint']):r for r in read(DATA/leg/'constraint-rms.csv')} for leg in LEGS}
    orders=[];spots=[]
    for r in read(DATA/'constraint-orders.csv'):
        key=(common_time(r['time_M']),r['mask'],r['constraint']);v=[float(norms[leg][key]['rms']) for leg in LEGS]
        assert v==[float(r[k]) for k in ('low','mid','high')]
        for i,k in enumerate(('order_low_mid','order_mid_high')):
            q=math.log(v[i]/v[i+1])/math.log(1.5) if min(v[i:i+2])>1e-13 else None
            assert (q is None and r[k]=='') or (q is not None and abs(q-float(r[k]))<1e-12)
        out={**r,'flag':'NEGATIVE' if any(float(r[k])<0 for k in ('order_low_mid','order_mid_high') if r[k]) else 'LOW_LT_3' if any(float(r[k])<3 for k in ('order_low_mid','order_mid_high') if r[k]) else 'UNDEFINED' if not all(r[k] for k in ('order_low_mid','order_mid_high')) else 'GE_3'}
        orders.append(out)
        if key in ((10.5,'far','Ham'),(10.5,'outer_boundary_shell','Theta')):
            spots.append(dict(time_M=key[0],mask=key[1],constraint=key[2],low=v[0],mid=v[1],high=v[2],
                order_low_mid=float(r['order_low_mid']),order_mid_high=float(r['order_mid_high']),
                low_cells=norms['E-low'][key]['cells'],mid_cells=norms['E-mid'][key]['cells'],high_cells=norms['E-high'][key]['cells'],check='RECOMPUTED_FROM_PER_LEG_RMS'))
    assert len(orders)==13*11*12 and len(spots)==2
    save('t8-constraint-orders.csv',orders);save('t8-reduction-spots.csv',spots)
    summary=[]
    for leg in LEGS:
        hh=[r for r in horizons if r['leg']==leg];nn=[r for r in native if r['leg']==leg];tt=np.array([r['time_M'] for r in hh]);rate=np.polyfit(tt,[r['relative_Q'] for r in hh],1)
        late=[r for r in hh if r['time_M']>=1.75];pfit=np.polyfit([r['time_M'] for r in late],[r['relative_Q'] for r in late],1)
        summary.append(dict(leg=leg,A0=hh[0]['A'],Q0=hh[0]['Q'],last_A=hh[-1]['A'],last_Q=hh[-1]['Q'],
            max_qualified_A_drift=max(abs(r['relative_A']) for r in hh),max_qualified_Q_drift=max(abs(r['relative_Q']) for r in hh),
            max_angular_A_drift=max(r['angular_drift_A'] for r in hh),max_angular_Q_drift=max(r['angular_drift_Q'] for r in hh),
            max_stopping_A_bound=max(r['stopping_drift_bound_A'] for r in hh),max_stopping_Q_bound=max(r['stopping_drift_bound_Q'] for r in hh),
            max_native_A_drift=max(abs(r['relative_A']) for r in nn),max_native_Q_drift=max(abs(r['relative_Q']) for r in nn),
            native_largest_A_time_M=max(nn,key=lambda r:abs(r['relative_A']))['time_M'],
            Q_relative_rate_all_per_M=rate[0],Q_relative_rate_late_per_M=pfit[0],Q_late_fit_rms=float(np.sqrt(np.mean((np.array([r['relative_Q'] for r in late])-np.polyval(pfit,[r['time_M'] for r in late]))**2))),
            projected_Q_drift_100M=pfit[0]*100,condition='LINEAR_RATE_PROJECTION_NOT_A_PREDICTION_OR_CONTINUUM_LIMIT'))
    save('t8-horizon-summary.csv',summary)
    projections=[]
    for r in summary:
        delta=r['last_Q']/r['Q0']-1;rate=r['Q_relative_rate_late_per_M']
        for endpoint in (25.375,49.875,100.,100.625):
            projections.append(dict(leg=r['leg'],endpoint_M=endpoint,measured_t10p5_relative_Q=delta,endpoint_average_rate_per_M=delta/10.5,
                all_time_fit_rate_per_M=r['Q_relative_rate_all_per_M'],late_fit_rate_per_M=rate,late_fit_RMS_relative_Q=r['Q_late_fit_rms'],
                anchored_late_linear_projection=delta+rate*(endpoint-10.5),endpoint_average_linear_projection=delta/10.5*endpoint,
                budget_crossing_time_late_linear_M=10.5+(2e-4-delta)/rate if r['leg']=='E-high' else '',
                condition='CONDITIONAL_LINEAR_PROJECTION;NO_CONTINUUM_OR_100M_PREDICTION'))
    save('t8-charge-rate-projections.csv',projections)
    rates=np.array([r['Q_relative_rate_late_per_M'] for r in summary]);p=math.log(abs((rates[0]-rates[1])/(rates[1]-rates[2])))/math.log(1.5)
    save('t8-charge-rate-model.csv',[dict(apparent_difference_order=p,fitted_limit_rate_per_M=rates[2]+(rates[2]-rates[1])/(1.5**p-1),
        condition='THREE_RATES_SINGLE_POWER_FIT_ONLY;NOT_A_CONTINUUM_LIMIT_MEASUREMENT')])
    for r in runrows:print('RUN',r)
    for r in summary:print('HORIZON',r)
    for mask in audit.MASKS:
        rr=[r for r in fields if r['mask']==mask and float(r['time_M'])>0 and r['field'] in audit.EVOLVED]
        print('SELF',mask,Counter(r['flag'] for r in rr))
    for r in spots:print('SPOT',r)

# Fixed diagnostic points; neither targets nor run-time mesh placement.
R=np.concatenate([np.linspace(.003,.02734375,2048,endpoint=False)]+[
    np.linspace(112/2**l/2,112/2**l,2048,endpoint=False) for l in range(11,2,-1)])
R=R[R<=24]
def load_plot(path):
    levels=[]
    with h5py.File(path) as f:
        cs=audit.names(f);assert cs==list(F)
        for l in range(13):
            g=f[f'level_{l}'];h=float(g.attrs['dx']);bs=g['boxes'][:]
            xy=np.array([[int(b[k]) for k in audit.BOX_KEYS] for b in bs]);lo=xy[:,:2].min(0);hi=xy[:,2:].max(0)
            v=np.full((len(F),hi[1]-lo[1]+1,hi[0]-lo[0]+1),np.nan)
            for key,raw in audit.blocks(f,l):
                x0,y0,x1,y1=key;assert np.isfinite(raw).all()
                v[:,y0-lo[1]:y1-lo[1]+1,x0-lo[0]:x1-lo[0]+1]=raw
            assert np.isfinite(v).all()
            levels.append(dict(h=h,lo=lo,hi=hi,a=v,bounds=(lo[0]*h-336,(hi[0]+1)*h-336,0,(hi[1]+1)*h)))
        return common_time(f['level_0'].attrs['time']),levels

def sample(levels,xy,n):
    # T6 compensated tensor point sampling, generalized only to centre=336/34 components.
    chosen=np.zeros(len(xy),int)
    for l,g in enumerate(levels):
        x0,x1,y0,y1=g['bounds'];chosen[(xy[:,0]>=x0)&(xy[:,0]<=x1)&(xy[:,1]>=y0)&(xy[:,1]<=y1)]=l
    out=np.empty((len(F),len(xy)))
    for l,g in enumerate(levels):
        use=chosen==l
        if not use.any():continue
        h,lo,v=g['h'],g['lo'],g['a'];assert lo[1]==0
        x=(xy[use,0]+336)/h-.5-lo[0];y=xy[use,1]/h-.5
        sx=np.clip(np.floor(x).astype(int)-(n//2-1),0,v.shape[2]-n);sy=np.minimum(np.floor(y).astype(int)-(n//2-1),v.shape[1]-n)
        ix=sx[:,None]+np.arange(n);raw=sy[:,None]+np.arange(n);iy=np.where(raw<0,-raw-1,raw)
        wx,wy=t6.t5.weights(x,sx,n),t6.t5.weights(y,sy,n)
        anchor=v[:,iy[:,n//2],ix[:,n//2]];z=anchor.copy()
        for j in range(n):
            for i in range(n):
                q=v[:,iy[:,j],ix[:,i]].copy();q[ODD]*=np.where(raw[:,j]<0,-1.,1.)
                z+=(q-anchor)*(wy[:,j]*wx[:,i])[None,:]
        out[:,use]=z
    return out,chosen

def profiles(leg):
    begin=time.monotonic();values=[];v8=[];meta=[];boundary=[];native=[];tails=[];outer=[];outer8=[]
    xy=(R[None,:,None]*VECS[:,None,:]).reshape(-1,2)
    depth=np.arange(.125,24.001,.125)
    outerxy=np.array([np.column_stack((336-depth,np.full_like(depth,y))) for y in (84,168,252)]+[
        np.column_stack((np.full_like(depth,x),336-depth)) for x in (0,168,-168)])
    for path in sorted((DATA/leg/'plt').glob('*.hdf5')):
        t,levels=load_plot(path);print('PROFILE',leg,t,flush=True)
        q,l=sample(levels,xy,6);qq,_=sample(levels,xy,8)
        values.append(q.reshape(len(F),3,-1).transpose(1,0,2));v8.append(qq.reshape(len(F),3,-1).transpose(1,0,2))
        ob,_=sample(levels,outerxy.reshape(-1,2),6);ob8,_=sample(levels,outerxy.reshape(-1,2),8)
        outer.append(ob.reshape(len(F),6,-1).transpose(1,0,2));outer8.append(ob8.reshape(len(F),6,-1).transpose(1,0,2))
        meta.append(dict(leg=leg,time_M=t,chi_min=min(g['a'][0].min() for g in levels),lapse_min=min(g['a'][13].min() for g in levels),
            floor_cells=sum(int(np.count_nonzero(g['a'][[0,13]]<=1e-12)) for g in levels),finite=True))
        # Native fourth-order radial Gamma curvature on each level; covered fields retained.
        for lev,g in enumerate(levels[1:],1):
            h=g['h'];x=(np.arange(g['lo'][0],g['hi'][0]+1)+.5)*h-336;y=(np.arange(g['hi'][1]+1)+.5)*h
            centre=int(np.argmin(abs(x-h/2)))
            for ri,ray in enumerate(RAYS):
                if ri==0:
                    use=x>0;rad=x[use];z=g['a'][11,0,use];step=h
                elif ri==1:
                    rad=y;z=g['a'][12,:,centre];step=h
                else:
                    ii=centre+np.arange(len(y));rad=np.sqrt(2)*y;z=(g['a'][11,np.arange(len(y)),ii]+g['a'][12,np.arange(len(y)),ii])/np.sqrt(2);step=np.sqrt(2)*h
                curvature=(-z[:-4]+16*z[1:-3]-30*z[2:-2]+16*z[3:-1]-z[4:])/(12*step*step);rad=rad[2:-2]
                annulus=rad>=.5*112/2**lev*(np.sqrt(2) if ri==2 else 1)
                for sign in (-1,1):
                    peaks,pr=find_peaks(sign*curvature,prominence=0)
                    for j,k in enumerate(peaks):
                        if rad[k]<.5*112/2**lev*(np.sqrt(2) if ri==2 else 1):continue
                        if pr['prominences'][j]<max(abs(curvature[annulus]))*1e-3:continue
                        w=peak_widths(sign*curvature,[k],prominence_data=(pr['prominences'][j:j+1],pr['left_bases'][j:j+1],pr['right_bases'][j:j+1]))
                        native.append(dict(leg=leg,time_M=t,level=lev,ray=ray,sign=sign,radius_M=rad[k],curvature=curvature[k],
                            prominence=pr['prominences'][j],width_M=w[0][0]*step,width_fine_cells=w[0][0]*step/h,width_receiving_cells=w[0][0]*step/(2*h),
                            basin_left_M=rad[pr['left_bases'][j]],basin_right_M=rad[pr['right_bases'][j]],
                            basin_clipped=bool(pr['left_bases'][j]==0 or pr['right_bases'][j]==len(rad)-1),
                            axis_condition='EXACT_DIAGONAL' if ri==2 else 'NATIVE_FIRST_ROW_OFFSET_H_OVER_2'))
        g=levels[0];h=g['h'];v=g['a'];x,y=audit.coordinates((*g['lo'],*g['hi']),h);d=np.minimum.reduce((x+336,336-x,336-y))
        use=(d>=0)&(d<24)&(y>1);bins=np.floor(d[use]/.25).astype(int);w=2*np.pi*y[use]*h*h
        ws=np.bincount(bins,weights=w,minlength=96);cnt=np.bincount(bins,minlength=96)
        for field in ('Ham','Mom','GaussE','Theta','lapse','chi','K','Gamma1','Gamma2','shift1','shift2','B1','B2','phi','Xi'):
            z=v[INDEX[field]][use];ss=np.bincount(bins,weights=w*z*z,minlength=96);mean=np.bincount(bins,weights=w*z,minlength=96)
            for bi in range(96):
                if ws[bi]:boundary.append(dict(leg=leg,time_M=t,field=field,depth_M=(bi+.5)*.25,cells=int(cnt[bi]),
                    volume_weight=ws[bi],rms=np.sqrt(ss[bi]/ws[bi]),mean=mean[bi]/ws[bi]))
        # Source BC uses unit radial speed and a 1/r tail for every component.
        if t==0:
            sl=(slice(4,-4),slice(4,-4));xx=x[sl];yy=y[sl];r=np.hypot(xx,yy);sel=(np.minimum.reduce((xx+336,336-xx,336-yy))<4)&(yy>1)
            asym={c:1. for c in ('chi','h11','h22','hww','lapse')}
            for field in ('chi','lapse','K','shift1','shift2','Gamma1','Gamma2','B1','B2','Ex','Ey','phi','Xi','Theta'):
                a=v[INDEX[field]];dx=(a[4:-4,5:-3]-a[4:-4,3:-5])/(2*h);dy=(a[5:-3,4:-4]-a[3:-5,4:-4])/(2*h)
                rhs=-(xx*dx+yy*dy+a[sl]-asym.get(field,0))/r
                ww=yy[sel];tails.append(dict(leg=leg,field=field,near_boundary_valid_cells=int(sel.sum()),
                    current_t0_sommerfeld_rhs_rms=float(np.sqrt(np.sum(ww*rhs[sel]**2)/sum(ww))),maximum=float(np.max(abs(rhs[sel]))),
                    condition='VALID_INTERIOR_TAIL_BAND_4_TO_5_CELLS_FROM_BOUNDARY;NATIVE_SECOND_ORDER_GRADIENT;NO_STATIC_READER'))
        del levels
    np.savez_compressed(TMP/f'{leg}-rays.npz',time=TIMES,radius=R,fields=F,values=values,values8=v8)
    np.savez_compressed(TMP/f'{leg}-outer.npz',time=TIMES,depth=depth,fields=F,values=outer,values8=outer8)
    save(f'{leg}-boundary.csv',boundary,TMP);save(f'{leg}-native-widths.csv',native,TMP);save(f'{leg}-snapshot-audit.csv',meta,TMP);save(f'{leg}-bc-tail.csv',tails,TMP)
    save(f'{leg}-performance.csv',[dict(leg=leg,seconds=time.monotonic()-begin,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,
        cache=str(TMP/f'{leg}-rays.npz'),cache_bytes=(TMP/f'{leg}-rays.npz').stat().st_size)],TMP)
    print('COMPLETE',leg,time.monotonic()-begin,'RSS',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)

def features():
    rows=[];tracks=[];native=[];boundary=[];boundaryfits=[];health=[];tails=[];amplitudes=[]
    # Fixed physical windows selected from the observed early outgoing ridge in the map.
    # They do not use a predicted characteristic speed or resolution-dependent positions.
    events=[(6,'axis',.875,.875,1.25,-1),(6,'equator',.875,.875,1.25,-1),
            (6,'corner',1.75,1.5,2.1,1),(5,'axis',2.625,1.9,3.3,-1),
            (5,'equator',2.625,1.9,3.3,-1),(5,'corner',4.375,3.8,4.8,-1),
            (4,'axis',6.125,4.5,6.9,-1),(4,'equator',6.125,4.5,6.9,-1),
            (4,'corner',8.75,7.5,9.6,-1),
            (5,'axis',8.75,2.5,3.4,-1),(5,'equator',8.75,2.5,3.4,-1),
            (5,'corner',10.5,3.9,4.8,-1)]
    for leg in LEGS:
        p=np.load(TMP/f'{leg}-rays.npz');v=p['values'];v8=p['values8'];rad=p['radius'];h0=float(params(leg)['L'])/int(params(leg)['N1'])
        nr=read(TMP/f'{leg}-native-widths.csv');native+=nr
        for lev,ray,t,left,right,sign in events:
            ri=RAYS.index(ray);ti=round(t/.875);sel=(rad>left)&(rad<right);r=rad[sel];out=[]
            for data in (v,v8):
                z=data[ti,ri];gn=VECS[ri,0]*z[11]+VECS[ri,1]*z[12]
                curv=np.gradient(np.gradient(gn,rad),rad)[sel];pk,pr=find_peaks(sign*curv,prominence=0)
                assert len(pk)>0
                j=int(np.argmax(pr['prominences']));k=pk[j]
                w=peak_widths(sign*curv,[k],prominence_data=(pr['prominences'][j:j+1],pr['left_bases'][j:j+1],pr['right_bases'][j:j+1]))
                wl=np.interp(w[2][0],np.arange(len(r)),r);wr=np.interp(w[3][0],np.arange(len(r)),r)
                out.append(dict(r=r[k],width=wr-wl,peak=curv[k],prom=pr['prominences'][j],z=z[:,np.where(sel)[0][k]],
                    basin_clipped=bool(pr['left_bases'][j]==0 or pr['right_bases'][j]==len(r)-1)))
            a,b=out;z=a['z'];vec=VECS[ri];det=z[1]*z[3]-z[2]*z[2]
            branch='LATE_BROAD' if t>=8.75 and lev==5 else 'EARLY_GAMMA_SHIFT'
            g=(vec[0]**2*z[3]-2*vec[0]*vec[1]*z[2]+vec[1]**2*z[1])/det;bn=vec[0]*z[14]+vec[1]*z[15]
            candidates=[q for q in nr if int(q['level'])==lev and q['ray']==ray and float(q['time_M'])==t and int(q['sign'])==sign and left<float(q['radius_M'])<right]
            current=min(candidates,key=lambda q:abs(float(q['radius_M'])-a['r'])) if candidates else None
            rows.append(dict(leg=leg,branch=branch,face_M=112/2**lev,level=lev,ray=ray,time_M=t,sign=sign,radius_M=a['r'],
                 width_M=a['width'],width_P8_M=b['width'],width_relative_sampling_spread=abs(b['width']/a['width']-1),
                 phase_sampling_spread_M=abs(a['r']-b['r']),fine_h_M=h0/2**lev,receiving_h_M=h0/2**(lev-1),
                 fine_cells=a['width']/(h0/2**lev),receiving_cells=a['width']/(h0/2**(lev-1)),
                 native_width_M=float(current['width_M']) if current else '',native_condition=current['axis_condition'] if current else '',
                 curvature_peak=a['peak'],curvature_prominence=a['prom'],basin_clipped=a['basin_clipped'],
                 lapse_speed=-bn+np.sqrt(1.8*z[13]*z[0]*g),shift_longitudinal_speed=-bn+np.sqrt(g),
                 shift_transverse_speed=-bn+np.sqrt(.75*g),light_speed=-bn+z[13]*np.sqrt(z[0]*g),
                 Gamma_n=vec[0]*z[11]+vec[1]*z[12],shift_n=bn,driver_B_n=vec[0]*z[16]+vec[1]*z[17],
                 lapse=z[13],K=z[5],Theta=z[10],chi=z[0],phi=z[18],Pi=z[19],Xi=z[27],
                 ROI_left_M=left,ROI_right_M=right,condition='UNQUALIFIED_UPSTREAM_FACE_P6_P8' if t==.875 else 'PRE_FACE;HALF_PROMINENCE_CURVATURE_WIDTH'))
            local=(rad>a['r']-a['width'])&(rad<a['r']+a['width'])
            for field in ('lapse','K','Theta','chi','phi','Pi','Xi','Ham','Gamma_n','shift_n','driver_B_n'):
                if field.endswith('_n'):
                    i={'Gamma_n':11,'shift_n':14,'driver_B_n':16}[field];q=vec[0]*v[ti,ri,i]+vec[1]*v[ti,ri,i+1];q8=vec[0]*v8[ti,ri,i]+vec[1]*v8[ti,ri,i+1]
                else:q=v[ti,ri,INDEX[field]];q8=v8[ti,ri,INDEX[field]]
                rr=rad[local];yy=q[local];line=yy[0]+(yy[-1]-yy[0])*(rr-rr[0])/(rr[-1]-rr[0]);defect=yy-line
                amplitudes.append(dict(leg=leg,branch=branch,face_M=112/2**lev,ray=ray,time_M=t,field=field,
                    current_field_range=float(np.ptp(yy)),local_detrended_amplitude=float(np.max(abs(defect))),
                    signed_peak_deviation=float(defect[np.argmax(abs(defect))]),P6_P8_absmax=float(np.max(abs(q[local]-q8[local]))),
                    condition='CURRENT_LOCAL_LINEAR_BACKGROUND;NO_STATIC_OR_TARGET_SUBTRACTION'))
        for ri,ray in enumerate(RAYS):
            a=[r for r in rows if r['leg']==leg and r['ray']==ray and r['branch']=='EARLY_GAMMA_SHIFT']
            good=[r for r in a if r['time_M']>=1.75]
            if len(good)>=2:
                x=[r['time_M'] for r in good];y=[r['radius_M'] for r in good];q=np.polyfit(x,y,1)
                tracks.append(dict(leg=leg,ray=ray,branch='EARLY_GAMMA_SHIFT',samples=len(good),speed=q[0],intercept=q[1],
                    phase_fit_RMS_M=float(np.sqrt(np.mean((np.array(y)-np.polyval(q,x))**2))),
                    characteristic_speed_mean=float(np.mean([r['shift_longitudinal_speed'] for r in good])),
                    lapse_speed_mean=float(np.mean([r['lapse_speed'] for r in good])),light_speed_mean=float(np.mean([r['light_speed'] for r in good]))))
        po=np.load(TMP/f'{leg}-outer.npz');d=po['depth'];a=po['values'];b=po['values8']
        for ray in range(6):
            for field in ('Ham','Mom','GaussE','Theta','lapse','chi','Gamma1','Gamma2'):
                for ti,t in enumerate(TIMES):
                    if field not in ('Ham','Mom','GaussE','Theta') and ti<2:continue
                    q=a[ti,ray,INDEX[field]];q8=b[ti,ray,INDEX[field]]
                    z=abs(q if field in ('Ham','Mom','GaussE','Theta') else q-a[ti-1,ray,INDEX[field]])
                    z8=abs(q8 if field in ('Ham','Mom','GaussE','Theta') else q8-b[ti-1,ray,INDEX[field]])
                    m=z.max()
                    for fraction in (.01,.1):
                        ii=np.where(z>m*fraction)[0];jj=np.where(z8>z8.max()*fraction)[0]
                        boundary.append(dict(leg=leg,time_M=t,ray=ray,field=field,threshold_fraction=fraction,
                            peak=float(m),peak_depth_M=float(d[z.argmax()]),inward_extent_M=float(d[ii[-1]]) if len(ii) else '',
                            extent_P8_M=float(d[jj[-1]]) if len(jj) else '',current_P6_P8_absmax=float(np.max(abs(q-q8))),
                            reference='SUCCESSIVE_POSITIVE_TIME_SNAPSHOTS' if field not in ('Ham','Mom','GaussE','Theta') else 'ABSOLUTE_CURRENT_CONSTRAINT'))
        for ray in range(6):
            for field in ('Ham','Theta','lapse'):
                rr=[r for r in boundary if r['leg']==leg and r['ray']==ray and r['field']==field and r['threshold_fraction']==.1 and r['time_M']>=1.75]
                for quantity in ('peak_depth_M','inward_extent_M'):
                    xx=np.array([r['time_M'] for r in rr]);yy=np.array([r[quantity] for r in rr]);q=np.polyfit(xx,yy,1)
                    boundaryfits.append(dict(leg=leg,ray=ray,field=field,feature=quantity,speed_normal_M_per_M=q[0],
                        intercept_M=q[1],fit_RMS_M=float(np.sqrt(np.mean((yy-np.polyval(q,xx))**2))),threshold_fraction=.1))
        health+=read(TMP/f'{leg}-snapshot-audit.csv');tails+=read(TMP/f'{leg}-bc-tail.csv')
    save('t8-pulse-widths.csv',rows);save('t8-pulse-speed-fits.csv',tracks);save('t8-native-curvature-candidates.csv',native)
    save('t8-pulse-field-amplitudes.csv',amplitudes)
    save('t8-boundary-features.csv',boundary);save('t8-boundary-speeds.csv',boundaryfits)
    save('t8-current-snapshot-audit.csv',health);save('t8-boundary-tail.csv',tails)
    for leg in LEGS:
        a=[r for r in rows if r['leg']==leg];print('WIDTHS',leg,[(r['face_M'],r['ray'],r['width_M'],r['receiving_cells'],r['width_relative_sampling_spread']) for r in a])
    for r in tracks:print('TRACK',r)
    for r in boundaryfits:
        if r['ray']==0:print('BOUNDARY',r)

def shell():
    m=module('t8maxwell',HERE/'t7-maxwell.py');from scipy.interpolate import CubicSpline
    m.selfcheck();rows=[]
    for leg in LEGS:
        for ti in (0,12):
            path=sorted((DATA/leg/'plt').glob('*.hdf5'))[ti]
            with h5py.File(path) as f:
                lev=f['level_12'];boxes=lev['boxes'][:];xy=np.array([[int(b[k]) for k in audit.BOX_KEYS] for b in boxes]);lo=xy[:,:2].min(0);hi=xy[:,2:].max(0)
                a=np.full((len(F),hi[1]-lo[1]+1,hi[0]-lo[0]+1),np.nan)
                for (x0,y0,x1,y1),raw in audit.blocks(f,12):a[:,y0-lo[1]:y1-lo[1]+1,x0-lo[0]:x1-lo[0]+1]=raw
                g=dict(a=a,lo=lo,hi=hi,h=float(lev.attrs['dx']),centre=336.,time=common_time(lev.attrs['time']))
            shape_path=DATA/leg/'qualified'/f't{ti:03d}'/'n96'/'shape-0-2.dat'
            s=np.fromstring(shape_path.read_text(),sep=' ');assert len(s)==98 and abs(s[1]-336)<1e-10
            theta=(np.arange(96)+.5)*np.pi/96
            shape=CubicSpline(np.r_[-theta[:4][::-1],theta,2*np.pi-theta[-4:][::-1]],np.r_[s[2:6][::-1],s[2:],s[-4:][::-1]])
            inv,ff,vol=m.geometry(a);D=m.displacement(a);density=vol*a[INDEX['GaussE']];sphere=CubicSpline([0,np.pi],[.02,.02])
            for nt,nr in ((96,64),(192,128)):
                th,rh,wt,r,w=m.quadrature(g,shape,nt,nr);x,y=r*np.cos(th[:,None]),r*np.sin(th[:,None])
                for n in (6,8):
                    integral=m.NORM*np.sum(w*m.sample(density[None],g,x,y,n,odd=[])[0]);qout=m.flux(D,g,sphere,nt,n);qh=m.flux(D,g,shape,nt,n)
                    rows.append(dict(leg=leg,time_M=g['time'],N_theta=nt,N_r=nr,point_nodes=n,dx_M=g['h'],
                        horizon_min_M=float(s[2:].min()),horizon_max_M=float(s[2:].max()),sphere_M=.02,faces_in_shell=0,
                        Q_H_same_quadrature=qh,Q_out_same_quadrature=qout,flux_gap=qout-qh,signed_Gauss_volume=integral,
                        volume_minus_gap=integral-(qout-qh),condition='CURRENT_PLOT_DENSITY;NUMERICAL_QUALIFIED_SHAPE;SAME_GAUSS_LEGENDRE_FLUX_QUADRATURE'))
            print('SHELL',leg,g['time'],rows[-1],flush=True)
    save('t8-gauss-shell.csv',rows)

def finish():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scienceplots
    plt.style.use(['science','no-latex'])
    plt.rcParams.update({'font.size':9,'figure.dpi':130,'savefig.dpi':180})
    figdir=HERE/'figures';figdir.mkdir(exist_ok=True)
    def figure(fig,name):
        fig.savefig(figdir/f't8-{name}.png',bbox_inches='tight');fig.savefig(figdir/f't8-{name}.pdf',bbox_inches='tight');plt.close(fig)
    cc=('tab:blue','tab:orange','tab:green')
    hh=read(HERE/'t8-qualified-horizons.csv');nn=read(HERE/'t8-native-horizon-history.csv');ww=read(HERE/'t8-pulse-widths.csv')
    oo=read(HERE/'t8-constraint-orders.csv');ss=read(HERE/'t8-field-orders.csv')
    summary=[];stats=[];mq=[];cost=[];arrivals=[]
    for mask in audit.MASKS:
        for kind,rr in (('28_evolved_fields',[r for r in ss if r['mask']==mask and r['field'] in audit.EVOLVED and float(r['time_M'])>0]),
                        ('12_constraints',[r for r in oo if r['mask']==mask and float(r['time_M'])>0])):
            flags=Counter(r['flag'] for r in rr)
            summary.append(dict(mask=mask,kind=kind,post_t0_rows=len(rr),negative=flags['NEGATIVE'],low_positive_lt3=flags['LOW_LT_3'],
                undefined=flags['UNDEFINED']+flags['UNDEFINED_ZERO_OR_FLOOR'],GE3=flags['GE_3'],
                interpretation='P_LT3_IS_A_REPORTED_DIAGNOSTIC_NOT_AN_ADDED_CONTRACT_THRESHOLD'))
        for con in audit.CONSTRAINTS:
            rr=[r for r in oo if r['mask']==mask and r['constraint']==con]
            v1=[float(r['order_low_mid']) for r in rr if r['order_low_mid']];v2=[float(r['order_mid_high']) for r in rr if r['order_mid_high']]
            stats.append(dict(mask=mask,constraint=con,min_low_mid=min(v1) if v1 else '',max_low_mid=max(v1) if v1 else '',
                min_mid_high=min(v2) if v2 else '',max_mid_high=max(v2) if v2 else '',negative_low_mid=sum(p<0 for p in v1),
                negative_mid_high=sum(p<0 for p in v2),undefined_low_mid=13-len(v1),undefined_mid_high=13-len(v2),
                endpoint_low=rr[-1]['low'],endpoint_mid=rr[-1]['mid'],endpoint_high=rr[-1]['high'],
                endpoint_order_low_mid=rr[-1]['order_low_mid'],endpoint_order_mid_high=rr[-1]['order_mid_high']))
    for leg in LEGS:
        rr=[r for r in read(HERE/'t8-mq-charge.csv') if r['leg']==leg]
        for radius in (20,50,100):
            key=f'Q{radius}';q=np.array([float(r[key]) for r in rr]);mq.append(dict(leg=leg,radius_M=radius,native_times=len(q),
                first_time_M=float(rr[0]['time_M']),last_time_M=float(rr[-1]['time_M']),first_Q=q[0],last_Q=q[-1],minimum_Q=q.min(),maximum_Q=q.max(),
                relative_change_first_to_last=q[-1]/q[0]-1,printed_charge_precision=1e-10,
                qualification='AUTHOR_MQ_NATIVE_QUADRATURE;NO_ANGULAR_EXTRAPOLATION;NO_T0_ROW_WRITTEN'))
        times={round(float(r['time_M']),9) for r in rr}
        for t in TIMES[1:]:assert any(abs(float(r['time_M'])-t)<1e-8 for r in rr)
        run=next(r for r in read(HERE/'t8-run-audit.csv') if r['leg']==leg);dt={'E-low':7/32,'E-mid':7/48,'E-high':7/72}[leg]
        for endpoint in (25.375,49.875,100.625):
            extra=round(endpoint/dt)-int(run['steps']);cost.append(dict(leg=leg,endpoint_M=endpoint,additional_coarse_steps=extra,
                median_step_seconds=float(run['steady_median_s']),maximum_step_seconds=float(run['steady_max_s']),
                additional_wall_hours_median=extra*float(run['steady_median_s'])/3600,
                additional_wall_hours_observed_max=extra*float(run['steady_max_s'])/3600,
                allocated_cores=128,condition='UNCHANGED_CHAIN;RATE_EXTRAPOLATION_NO_EXTRA_IO_OR_QUEUE_MARGIN'))
        fit=[r for r in read(HERE/'t8-pulse-speed-fits.csv') if r['leg']==leg]
        for ray in RAYS:
            f=next(r for r in fit if r['ray']==ray);speed=float(f['speed']);intercept=float(f['intercept'])
            for lev in range(1,13):
                face=112/2**lev;rr=face*(np.sqrt(2) if ray=='corner' else 1);arrival=(rr-intercept)/speed
                arrivals.append(dict(leg=leg,ray=ray,level=lev,face_M=face,ray_intersection_radius_M=rr,arrival_phase_fit_M=arrival if arrival>=.875 else '',
                     phase_time_sampling_uncertainty_M=.4375,observed_before_endpoint=bool(rr<speed*10.5+intercept),
                     qualification='NOT_CROSSED' if arrival>10.5 else 'BEFORE_FIRST_POST_T0_PLOT_UNRESOLVED' if arrival<.875 else 'PHASE_FIT_NOT_EXACT_CROSSING_CAPTURE'))
    save('t8-order-flag-counts.csv',summary);save('t8-constraint-summary.csv',stats);save('t8-mq-summary.csv',mq)
    save('t8-continuation-cost.csv',cost);save('t8-face-arrivals.csv',arrivals)
    fig,ax=plt.subplots(2,1,figsize=(8,6.5),sharex=True)
    for leg,color in zip(LEGS,cc):
        h=[r for r in hh if r['leg']==leg];n=[r for r in nn if r['leg']==leg]
        for j,c in enumerate(('A','Q')):
            t=np.array([float(r['time_M']) for r in h]);v=np.array([float(r[f'relative_{c}']) for r in h]);u=np.array([float(r[f'angular_drift_{c}'])+float(r[f'stopping_drift_bound_{c}']) for r in h])
            ax[j].plot(t,v,'o-',color=color,label=leg+' qualified N96');ax[j].fill_between(t,v-u,v+u,color=color,alpha=.2)
            ax[j].plot([float(r['time_M']) for r in n],[float(r[f'relative_{c}']) for r in n],':',color=color,alpha=.5,label=leg+' native, unqualified' if j==0 else None)
        ax[0].set_ylabel(r'$A/A_0-1$');ax[1].set_ylabel(r'$Q/Q_0-1$')
    for a,budget in zip(ax,(1e-3,2e-4)):
        a.axhline(budget,color='k',ls='--',lw=.8);a.axhline(-budget,color='k',ls='--',lw=.8);a.axhline(0,color='.5',lw=.5)
    zoom=ax[1].inset_axes([.08,.12,.43,.43]);high=[r for r in hh if r['leg']=='E-high']
    zoom.plot([float(r['time_M']) for r in high],[float(r['relative_Q']) for r in high],'o-',color=cc[2],ms=3)
    zoom.axhline(2e-4,color='k',ls='--',lw=.7);zoom.set_ylim(-2e-5,2.2e-4);zoom.set_title('E-high charge and budget',fontsize=8);zoom.tick_params(labelsize=7)
    ax[0].legend(ncol=2,fontsize=7);ax[1].set_xlabel('time / M');figure(fig,'horizons')
    fig,axs=plt.subplots(4,3,figsize=(12,11),sharex=True)
    for a,mask in zip(axs.flat,audit.MASKS):
        for color,leg,key in zip(cc,LEGS,('low','mid','high')):
            for con,style in zip(('Ham','Mom','GaussE'),('-','--',':')):
                r=[r for r in oo if r['mask']==mask and r['constraint']==con];y=[max(float(z[key]),1e-30) for z in r]
                a.semilogy([float(z['time_M']) for z in r],y,color=color,ls=style,label=f'{leg} {con}' if mask=='horizon' else None)
        a.set_title(mask.replace('_',' '));a.set_ylim(bottom=1e-19);a.set_xlabel('time / M');a.set_ylabel('coordinate-volume RMS')
    axs.flat[-1].axis('off');h,l=axs.flat[0].get_legend_handles_labels();axs.flat[-1].legend(h,l,loc='center',fontsize=8)
    figure(fig,'constraint-rms')
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for leg,color in zip(LEGS,cc):
        for ray,marker in zip(('axis','corner'),('o','s')):
            r=[r for r in ww if r['leg']==leg and r['ray']==ray and r['branch']=='EARLY_GAMMA_SHIFT' and r['condition']!='UNQUALIFIED_UPSTREAM_FACE_P6_P8']
            ax[0].plot([float(z['face_M']) for z in r],[float(z['width_M']) for z in r],marker+'-',color=color,label=f'{leg} {ray}')
            ax[1].plot([float(z['face_M']) for z in r],[float(z['receiving_cells']) for z in r],marker+'-',color=color)
    ax[1].axhline(4,color='k',ls='--');ax[0].set_ylabel('early curvature width / M');ax[1].set_ylabel('width / receiving h');
    for a in ax:
        a.set_xlabel('square face / M');a.set_xscale('log',base=2);a.set_xticks([1.75,3.5,7]);a.set_xticklabels(['1.75','3.5','7'])
    ax[0].legend(fontsize=7);figure(fig,'face-widths')
    fig,axs=plt.subplots(2,3,figsize=(12,7),sharex=True)
    for col,(leg,color) in enumerate(zip(LEGS,cc)):
        p=np.load(TMP/f'{leg}-outer.npz');d=p['depth'];v=p['values']
        for row,field in enumerate(('Ham','Theta')):
            ax=axs[row,col];z=np.log10(np.maximum(abs(v[:,0,INDEX[field]]),1e-19))
            im=ax.pcolormesh(d,TIMES,z,vmin=-14,vmax=-7,shading='auto',cmap='magma');ax.set_title(leg+' '+field);ax.set_xlabel('inward depth at x=336, y=84 / M');ax.set_ylabel('time / M');ax.set_xlim(0,16)
    fig.colorbar(im,ax=list(axs.flat),label='log10 |current constraint|',shrink=.7);figure(fig,'outer-boundary')
    fig,axs=plt.subplots(2,3,figsize=(12,7),sharex=True,sharey=True)
    p=np.load(TMP/'E-high-rays.npz');r=p['radius'];v=p['values'];sel=(r>=.4)&(r<=11)
    for ax,field in zip(axs.flat,('Gamma1','shift1','lapse','K','Theta','Ham')):
        z=v[:,0,INDEX[field]][:,sel]
        # Current successive-snapshot differences identify transport; no static targets.
        if field in ('Gamma1','shift1','lapse','K'):z=np.r_[np.zeros_like(z[:2]),np.diff(z[1:],axis=0)]
        im=ax.pcolormesh(r[sel],TIMES,np.log10(np.maximum(abs(z),1e-22)),shading='auto',cmap='magma');ax.set_title(field+(' current time increment' if field in ('Gamma1','shift1','lapse','K') else ' current'))
        for face in (1.75,3.5,7.):ax.axvline(face,color='white',lw=.5,ls='--')
        ax.set_xlabel('axis radius / M');ax.set_ylabel('time / M');fig.colorbar(im,ax=ax,shrink=.7)
    fig.tight_layout();figure(fig,'pulse-map')
    fig,axs=plt.subplots(3,2,figsize=(10,8),sharex=True)
    for color,leg in zip(cc,LEGS):
        p=np.load(TMP/f'{leg}-rays.npz');r=p['radius'];v=p['values'];sel=(r>2.3)&(r<2.95)
        for ax,field in zip(axs.flat,('Gamma1','shift1','lapse','K','Theta','B1')):
            z=v[3,0,INDEX[field]][sel];ax.plot(r[sel],z,color=color,label=leg);ax.set_title(field)
    axs[0,0].legend(fontsize=8)
    for ax in axs.flat:ax.set_xlabel('axis radius at t=2.625 M');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    fig.tight_layout();figure(fig,'incident-profiles')
    print('FINISH',len(summary),'flag summaries',len(stats),'constraint summaries',len(cost),'cost rows')

def check():
    assert len(read(HERE/'t8-field-orders.csv'))==5291 and len(read(HERE/'t8-constraint-orders.csv'))==1716
    assert len(read(HERE/'t8-qualified-horizons.csv'))==39 and len(read(HERE/'t8-plot-census.csv'))==39
    assert set(r['floor_cells'] for r in read(HERE/'t8-current-snapshot-audit.csv'))=={'0'}
    assert all(float(r['rms'])==0 for leg in LEGS for r in read(DATA/leg/'constraint-rms.csv') if r['constraint'] in ('GaussB','Lambda'))
    for leg in LEGS:
        # Sample the same polynomial at the exact axis and inside a rectilinear patch.
        h=1/64;lo=np.array([0,0]);hi=np.array([79,39]);x=(np.arange(80)+.5)*h-336;y=(np.arange(40)+.5)*h
        X,Y=np.meshgrid(x,y);a=np.zeros((len(F),40,80));a[0]=(X+336)**4+Y**4;a[2]=((X+336)**2)*Y**3
        g=dict(h=h,lo=lo,hi=hi,a=a,bounds=(-336,-334.75,0,.625));xy=np.array([[-335.625,0],[-335.5,.1875]])
        for n in (6,8):
            v,l=sample([g],xy,n);assert np.max(abs(v[0]-((xy[:,0]+336)**4+xy[:,1]**4)))<2e-15
            assert np.max(abs(v[2]-(xy[:,0]+336)**2*xy[:,1]**3))<2e-15
    print('CHECK PASS: all rows, finite/floor census, zero GaussB/Lambda, tensor point/parity polynomial sampling')

def localize():
    rows=[]
    for leg in LEGS:
        paths=sorted((DATA/leg/'plt').glob('*.hdf5'))
        for ti in (5,8,12):
            t,gs=load_plot(paths[ti])
            for mask in ('far','far_core','exterior_wake'):
                # Composite coordinate-volume norm partition, identical physical collars.
                parts=[]
                for l,g in enumerate(gs):
                    h=g['h'];x,y=audit.coordinates((*g['lo'],*g['hi']),h);valid=audit.mask_cells(mask,x,y)
                    ext=np.maximum(abs(x),y)
                    if l<12:valid &= ext>=112/2**(l+1)
                    if not valid.any():continue
                    weight=2*np.pi*y*h*h
                    for region in ('whole','face3.5_collar','face7_collar','face14_collar'):
                        selection=valid.copy()
                        if region!='whole':
                            face=float(region[4:-7]);low_receiving=(7/8)/(112/face/2);selection &= abs(ext-face)<=2*low_receiving
                        if not selection.any():continue
                        for field in ('Ham','Mom','GaussE','Theta'):
                            z=g['a'][INDEX[field]];w=weight[selection];q=z[selection];energy=float(np.sum(w*q*q));idx=np.flatnonzero(selection)[np.argmax(abs(q))];ij=np.unravel_index(idx,selection.shape)
                            parts.append(dict(leg=leg,time_M=t,mask=mask,level=l,region=region,field=field,cells=int(selection.sum()),
                                coordinate_volume=float(w.sum()),RMS=float(np.sqrt(energy/w.sum())),squared_error_integral=energy,
                                maximum=float(np.max(abs(q))),peak_x_M=float(x[ij]),peak_y_M=float(y[ij]),dx_M=h))
                for r in parts:
                    total=sum(a['squared_error_integral'] for a in parts if a['region']=='whole' and a['field']==r['field']);r['mask_squared_error_fraction']=r['squared_error_integral']/total if total else 0
                rows+=parts
            print('LOCALIZE',leg,t,flush=True);del gs
    save('t8-interface-localization.csv',rows)

if __name__=='__main__':globals()[sys.argv[1]](*sys.argv[2:])
