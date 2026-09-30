#!/usr/bin/env python3
"""Read-only fixed-window exp-0020 front diagnostics. No static reader."""
import argparse, itertools, json, math, resource, sys, time
from pathlib import Path
sys.dont_write_bytecode=True
import h5py
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import CubicSpline
import importlib.util
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('t9t8',HERE/'t8-analyze.py')
t8=importlib.util.module_from_spec(spec);spec.loader.exec_module(t8)
audit=t8.audit
TMP=Path('/private/tmp/ems-t9');TMP.mkdir(exist_ok=True)
NAMES=('beta','Gamma','Bdriver','Gamma_metric','CGamma','Gamma_metric_2D','Gamma_metric_ww','chi','lapse','K','Theta','g_normal')
FIELDS=('beta','Gamma','Wout','Win','Gamma_metric','CGamma','Gamma_metric_2D','Gamma_metric_ww')
BASE=('shift1','shift2','Gamma1','Gamma2','B1','B2','chi','lapse','K','Theta','h11','h12','h22')
ALL=BASE+('Gmetric1','Gmetric2','G2D1','G2D2','Gww1','Gww2')
ODD=[ALL.index(k) for k in ('shift2','Gamma2','B2','h12','Gmetric2','G2D2','Gww2')]
TIMES=(0.,.875,1.75,2.625,4.375)
def window(t):
    return (.0001, .04) if t==0 else (t-.65,t+.65)

def load(path):
    """Native fourth-order h derivatives with file ghosts and ww cartoon terms."""
    levels=[]
    with h5py.File(path) as f:
        cs=audit.names(f);idx={c:i for i,c in enumerate(cs)}
        for l in range(13):
            g=f[f'level_{l}'];h=float(g.attrs['dx']);gx,gy=map(int,g['data_attributes'].attrs['outputGhost'])
            box=np.array([[int(b[k]) for k in audit.BOX_KEYS] for b in g['boxes'][:]])
            lo=box[:,:2].min(0);hi=box[:,2:].max(0)
            out=np.full((len(ALL),hi[1]-lo[1]+1,hi[0]-lo[0]+1),np.nan)
            for key,a in audit.blocks(f,l,True):
                v,inv,dh,dw=audit.z_fields(a,cs,h,gx,gy);ny,nx=v.shape[1:]
                yy=(np.arange(key[1],key[3]+1)+.5)[:,None]*h
                geom=[];direct=[];cartoon=[]
                for i in range(2):
                    direct_i=np.zeros((ny,nx))
                    for j,k,m in itertools.product(range(2),repeat=3):
                        direct_i+=.5*inv[j,k]*inv[i,m]*(dh[j,m,k]+dh[k,m,j]-dh[m,j,k])
                    ww=((1. if i==1 else 0.)-inv[i,1]*v[idx['hww']])/yy
                    for j in range(2):ww-=.5*inv[i,j]*dw[j]
                    ww/=v[idx['hww']]
                    direct.append(direct_i);cartoon.append(ww);geom.append(direct_i+ww)
                q=np.stack([v[idx[k]] for k in BASE]+geom+direct+cartoon)
                assert np.isfinite(q).all()
                x0,y0,x1,y1=key;out[:,y0-lo[1]:y1-lo[1]+1,x0-lo[0]:x1-lo[0]+1]=q
            assert np.isfinite(out).all()
            levels.append(dict(a=out,h=h,lo=lo,hi=hi,bounds=(lo[0]*h-336,(hi[0]+1)*h-336,0,(hi[1]+1)*h)))
        return t8.common_time(f['level_0'].attrs['time']),levels

def basis(q,first,n,maxder=3):
    """Product-rule Lagrange weights, no generated arithmetic artifact."""
    z=q-first;ans=np.zeros((maxder+1,len(z),n))
    for i in range(n):
        other=[j for j in range(n) if j!=i];den=math.prod(i-j for j in other)
        for d in range(maxder+1):
            for removed in itertools.combinations(other,d):
                term=np.ones(len(z))
                for j in other:
                    if j not in removed:term*=z-j
                ans[d,:,i]+=math.factorial(d)*term/den
    return ans

def sample(levels,r,vec,n,force=None):
    xy=r[:,None]*vec[None,:];chosen=np.zeros(len(r),int)
    for l,g in enumerate(levels):
        x0,x1,y0,y1=g['bounds'];chosen[(xy[:,0]>=x0)&(xy[:,0]<=x1)&(xy[:,1]>=y0)&(xy[:,1]<=y1)]=l
    if force is not None:chosen[:]=force
    out=np.zeros((4,len(ALL),len(r)))
    for l in np.unique(chosen):
        use=chosen==l;g=levels[l];v=g['a'];h=g['h'];lo=g['lo']
        x=(xy[use,0]+336)/h-.5-lo[0];y=xy[use,1]/h-.5
        sx=np.clip(np.floor(x).astype(int)-n//2+1,0,v.shape[2]-n);sy=np.minimum(np.floor(y).astype(int)-n//2+1,v.shape[1]-n)
        ix=sx[:,None]+np.arange(n);raw=sy[:,None]+np.arange(n);iy=np.where(raw<0,-raw-1,raw)
        assert ix.max()<v.shape[2] and iy.max()<v.shape[1]
        wx=basis(x,sx,n);wy=basis(y,sy,n)
        anchor=v[:,iy[:,n//2],ix[:,n//2]].copy()
        if raw[:,n//2].min()<0:anchor[ODD]*=np.where(raw[:,n//2]<0,-1.,1.)
        result=np.zeros((4,len(ALL),len(x)));result[0]=anchor
        for j,i in itertools.product(range(n),repeat=2):
            q=v[:,iy[:,j],ix[:,i]].copy();q[ODD]*=np.where(raw[:,j]<0,-1.,1.)
            q-=anchor
            for d in range(4):
                weight=sum(math.comb(d,k)*vec[0]**k*vec[1]**(d-k)*wx[k,:,i]*wy[d-k,:,j] for k in range(d+1))/h**d
                result[d]+=q*weight
        out[:,:,use]=result
    # Longitudinal components and local inverse normal conformal metric.
    def dot(a,i):return vec[0]*a[:,i]+vec[1]*a[:,i+1]
    projected=np.zeros((4,len(NAMES),len(r)))
    for name,i in (('beta',0),('Gamma',2),('Bdriver',4),('Gamma_metric',13),('Gamma_metric_2D',15),('Gamma_metric_ww',17)):
        projected[:,NAMES.index(name)]=dot(out,i)
    projected[:,NAMES.index('CGamma')]=projected[:,1]-projected[:,3]
    for name,i in (('chi',6),('lapse',7),('K',8),('Theta',9)):projected[:,NAMES.index(name)]=out[:,i]
    h11,h12,h22=out[0,10:13];det=h11*h22-h12*h12
    projected[0,-1]=(vec[0]**2*h22-2*vec[0]*vec[1]*h12+vec[1]**2*h11)/det
    projected[1:,-1]=np.nan  # Only the current value of g is used; g derivatives are not computed.
    return projected,chosen

def projection(q,frontindex):
    b=q[0,0,frontindex];g=q[0,-1,frontindex];c=math.sqrt(g);mu=.75
    assert abs(b)<c
    pcoef=c/mu;bout=c/(mu*(c-b));bin_=c/(mu*(c+b))
    value=np.stack((q[:,0],q[:,1],q[:,1],q[:,1],q[:,3],q[:,4],q[:,5],q[:,6]),axis=1)
    value[3,2:4]=np.nan  # W''' needs beta''''; it is not requested or computed.
    for d in range(3):
        value[d,2]=q[d,1]-pcoef*q[d+1,0]-bout*q[d,2]
        value[d,3]=q[d,1]+pcoef*q[d+1,0]-bin_*q[d,2]
    return value,dict(beta_f=b,g_f=g,c_f=c,outgoing_speed=c-b,incoming_speed=-c-b,p_coefficient=pcoef,Bout_coefficient=bout,Bin_coefficient=bin_)

def background(r,u,t,derivative=0):
    use=(abs(r-t)>=.4)&(abs(r-t)<=.6)
    fit=np.polynomial.Polynomial.fit(r[use]-t,u[use],3)
    return fit.deriv(derivative)(r-t)

def signed_integrals(spline,left,right):
    """Exact positive/negative integrals of the C2 cubic's linear second derivative."""
    x=np.r_[left,spline.x[(spline.x>left)&(spline.x<right)],right]
    y=spline.derivative(2)(x);pos=neg=0.
    for a,b,ya,yb in zip(x[:-1],x[1:],y[:-1],y[1:]):
        if ya*yb<0:
            root=a+(b-a)*(-ya)/(yb-ya)
            areas=((root-a)*ya/2,(b-root)*yb/2)
        else:areas=((b-a)*(ya+yb)/2,)
        for area in areas:
            if area>=0:pos+=area
            else:neg+=area
    return pos,neg

def front(r,q,t):
    curvature=q[2,1]-background(r,q[0,1],t,2)
    sel=np.flatnonzero(abs(r-t)<=.3)
    k=sel[np.argmax(abs(curvature[sel]))]
    peaks,props=find_peaks(abs(curvature),prominence=0)
    assert k in peaks,(t,r[k])
    j=int(np.flatnonzero(peaks==k)[0]);w=peak_widths(abs(curvature),[k],prominence_data=(props['prominences'][j:j+1],props['left_bases'][j:j+1],props['right_bases'][j:j+1]))
    return int(k),float(w[0][0]*(r[1]-r[0])),float(curvature[k])

def profiles(leg):
    start=time.monotonic();paths={}
    for p in sorted((t8.DATA/leg/'plt').glob('*.hdf5')):
        with h5py.File(p) as f:paths[t8.common_time(f['level_0'].attrs['time'])]=p
    assert len(paths)==13
    metadata=[]
    for t in TIMES:
        tt,levels=load(paths[t]);assert tt==t;lo,hi=window(t)
        r=np.linspace(lo,hi,2601) if t else np.geomspace(lo,hi,1601)
        for ri,(ray,vec) in enumerate(zip(t8.RAYS,t8.VECS)):
            q6,lev=sample(levels,r,vec,6);q8,_=sample(levels,r,vec,8)
            if t:
                k,w,amp=front(r,q6,t);k8,w8,amp8=front(r,q8,t)
            else:k=int(np.argmin(abs(r-.006)));w=w8=amp=amp8=math.nan;k8=k
            v6,coeff=projection(q6,k);v8,_=projection(q8,k)
            # Same frozen coefficients for P6/P8: sensitivity is interpolation only.
            for d in range(3):
                v8[d,2]=q8[d,1]-coeff['p_coefficient']*q8[d+1,0]-coeff['Bout_coefficient']*q8[d,2]
                v8[d,3]=q8[d,1]+coeff['p_coefficient']*q8[d+1,0]-coeff['Bin_coefficient']*q8[d,2]
            parent=max(l for l,g in enumerate(levels) if (hi*vec[0]<=g['bounds'][1] and hi*vec[1]<=g['bounds'][3]))
            parent6,parentlev=sample(levels,r,vec,6,force=parent)
            parentv,_=projection(parent6,k)
            np.savez_compressed(TMP/f'{leg}-{t:g}-{ray}.npz',r=r,values=v6,values8=v8,q=q6,q8=q8,levels=lev,parent=parentv,parent_level=parent,
                front_index=k,width=w,front_index8=k8,width8=w8,coefficients=json.dumps(coeff),fields=FIELDS,source_fields=NAMES)
            metadata.append(dict(leg=leg,time_M=t,ray=ray,front_M=r[k],front8_M=r[k8],width_M=w,width8_M=w8,
                Gamma_curvature_signed=amp,Gamma_curvature_signed8=amp8,level=int(lev[k]),h_M=levels[int(lev[k])]['h'],parent_level=parent,
                levels_in_delta0p4=','.join(map(str,np.unique(lev[abs(r-r[k])<=.4]))),native_ghosts='FILE_GHOSTS_GE_2',
                **coeff,path=str(paths[t])))
        print('PROFILE',leg,t,'elapsed',round(time.monotonic()-start,2),flush=True)
        del levels
    t8.save(f'{leg}-metadata.csv',metadata,TMP)
    print('DONE',leg,'seconds',time.monotonic()-start,'peak_RSS_GB',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9,flush=True)

def measure():
    rows=[];jumps=[];scaling=[];profile_rows=[];launch=[]
    for leg in t8.LEGS:
        for meta in t8.read(TMP/f'{leg}-metadata.csv'):
            t=float(meta['time_M']);ray=meta['ray'];d=np.load(TMP/f'{leg}-{t:g}-{ray}.npz')
            r=d['r'];v=d['values'];vv=d['values8'];levels=d['levels'];k=int(d['front_index']);rf=r[k]
            base=dict(leg=leg,time_M=t,ray=ray,front_M=rf,h_M=float(meta['h_M']),level=int(meta['level']))
            # Signed profiles, moderately thinned; exact dense cache retained/hash recorded.
            keep=np.unique(np.r_[np.arange(0,len(r),20),np.flatnonzero(abs(r-t)<=.2)[::4]]) if t else np.arange(0,len(r),5)
            for j in keep:
                row=dict(**base,r_M=r[j],native_level=int(levels[j]))
                for fi,f in enumerate(FIELDS):
                    row[f]=v[0,fi,j];row[f+'_P8_minus_P6']=vv[0,fi,j]-v[0,fi,j]
                    if fi<3:row[f+'_d1']=v[1,fi,j];row[f+'_d2']=v[2,fi,j]
                profile_rows.append(row)
            if t==0:
                for fi,f in enumerate(FIELDS[:6]):
                    launch.append(dict(**base,field=f,maximum_abs=float(max(abs(v[0,fi]))),minimum=float(min(v[0,fi])),maximum=float(max(v[0,fi])),
                        max_abs_first_derivative=float(max(abs(v[1,fi]))),P6_P8_spread=float(max(abs(vv[0,fi]-v[0,fi]))),condition='PUNCTURE_NOT_SAMPLED_R_GE_0P0001'))
                continue
            use=abs(r-rf)<=.1
            for fi,f in enumerate(FIELDS):
                bg=background(r,v[0,fi],t)
                feature=v[0,fi]-bg;curv=v[2,fi]-background(r,v[0,fi],t,2)
                peak=np.flatnonzero(abs(r-rf)<=.1)[np.argmax(abs(curv[use]))]
                peaks,pr=find_peaks(abs(curv),prominence=0)
                if peak in peaks:
                    j=int(np.flatnonzero(peaks==peak)[0]);wd=peak_widths(abs(curv),[peak],prominence_data=(pr['prominences'][j:j+1],pr['left_bases'][j:j+1],pr['right_bases'][j:j+1]))[0][0]*(r[1]-r[0])
                else:wd=math.nan
                rows.append(dict(**base,field=f,front_width_Gamma_M=float(d['width']),field_curvature_peak_M=r[peak],field_curvature_halfprom_width_M=wd,
                    current_value_at_front=v[0,fi,k],maximum_abs_current_fixed_window=float(max(abs(v[0,fi]))),
                    detrended_peak_fixed_delta0p1=float(max(abs(feature[use]))),signed_detrended_at_front=feature[k],
                    max_abs_P6_P8_fixed_delta0p1=float(max(abs(vv[0,fi,use]-v[0,fi,use]))),
                    curvature_peak_signed=curv[peak],fine_cells=wd/float(meta['h_M']),receiving_cells=wd/(2*float(meta['h_M'])),
                    profile_bound_condition='THREE_GRIDS_ONLY_NOT_CONTINUUM_CERTIFICATE'))
            for delta in (.1,.2,.4):
                left=rf-delta;right=rf+delta
                use=(r>=left)&(r<=right);rr=r[use];lv=levels[use]
                for fi,f in enumerate(FIELDS[:6]):
                    val=v[:,fi];val8=vv[:,fi]
                    spline=CubicSpline(r,val[0]);spline8=CubicSpline(r,val8[0])
                    endpoint=spline.derivative()(np.array([left,right]));ep8=spline8.derivative()(np.array([left,right]))
                    analytic=np.interp(right,r,val[1])-np.interp(left,r,val[1])
                    # C2 reconstruction joins the P6/P8 samples. Its J is exactly the
                    # endpoint difference; lobe integrals include interpolation switch points.
                    J=float(endpoint[1]-endpoint[0]);J8=float(ep8[1]-ep8[0])
                    positive,negative=signed_integrals(spline,left,right)
                    assert abs((positive+negative)-J)<1e-12*max(1.,abs(positive),abs(negative))
                    half=CubicSpline(r[::2],val[0,::2]);Jhalf=float(half.derivative()(right)-half.derivative()(left))
                    fitvals=[];fitders=[]
                    for side in (-1,1):
                        which=((r-rf)*side>=delta/2)&((r-rf)*side<=delta)
                        fit=np.polynomial.Polynomial.fit(r[which]-rf,val[0,which],3)
                        fitvals.append(float(fit(0)));fitders.append(float(fit.deriv()(0)))
                    parentspline=CubicSpline(r,d['parent'][0,fi]);parentJ=float(parentspline.derivative()(right)-parentspline.derivative()(left))
                    jumps.append(dict(**base,field=f,delta_M=delta,J_endpoint=J,J_P8=J8,J_P8_minus_P6=J8-J,
                        integral_positive=positive,integral_negative=negative,J_integrated=positive+negative,
                        quadrature_minus_endpoint=positive+negative-J,analytic_piecewise_endpoint=analytic,C2_minus_analytic_endpoint=J-analytic,
                        J_half_density=Jhalf,J_half_minus_full=Jhalf-J,
                        one_sided_left_value=np.interp(left,r,val[0]),
                        one_sided_right_value=np.interp(right,r,val[0]),one_sided_left_d1=endpoint[0],one_sided_right_d1=endpoint[1],
                        extrapolated_left_value=fitvals[0],extrapolated_right_value=fitvals[1],extrapolated_value_jump=fitvals[1]-fitvals[0],
                        extrapolated_left_d1=fitders[0],extrapolated_right_d1=fitders[1],extrapolated_d1_jump=fitders[1]-fitders[0],
                        parent_J_endpoint=parentJ,parent_minus_composite_J=parentJ-J,
                        levels_in_integral=','.join(map(str,np.unique(lv))),same_valid_level=len(np.unique(lv))==1,
                        condition='RAW_J_INCLUDES_SMOOTH_BACKGROUND;DELTA_SENSITIVITY_REQUIRED'))
    t8.save('t9-front-metrics.csv',rows);t8.save('t9-jump-integrals.csv',jumps);t8.save('t9-signed-profiles.csv',profile_rows);t8.save('t9-launch.csv',launch)
    t8.save('t9-front-census.csv',[r for leg in t8.LEGS for r in t8.read(TMP/f'{leg}-metadata.csv')])
    for t,ray,f in itertools.product(TIMES[1:],t8.RAYS,FIELDS):
        subset=[next(r for r in rows if r['leg']==leg and r['time_M']==t and r['ray']==ray and r['field']==f) for leg in t8.LEGS]
        z=[r['detrended_peak_fixed_delta0p1'] for r in subset];w=[r['field_curvature_halfprom_width_M'] for r in subset]
        scaling.append(dict(time_M=t,ray=ray,field=f,low_peak=z[0],mid_peak=z[1],high_peak=z[2],
            amplitude_order_low_mid=math.log(z[0]/z[1])/math.log(1.5),amplitude_order_mid_high=math.log(z[1]/z[2])/math.log(1.5),
            low_width=w[0],mid_width=w[1],high_width=w[2],width_order_low_mid=math.log(w[0]/w[1])/math.log(1.5) if min(w[:2])>0 else math.nan,
            width_order_mid_high=math.log(w[1]/w[2])/math.log(1.5) if min(w[1:])>0 else math.nan))
    t8.save('t9-fixed-window-scaling.csv',scaling)
    gamma=[];differences=[];sensitivity=[];probes=[]
    for leg,ray in itertools.product(t8.LEGS,t8.RAYS):
        d=np.load(TMP/f'{leg}-0-{ray}.npz');r=d['r'];v=d['values'];q=d['q'];q8=d['q8']
        for rr in (.0001,.0003,.001,.003,.006,.02):
            probes.append(dict(leg=leg,ray=ray,r_M=rr,beta=float(np.interp(rr,r,v[0,0])),beta_r=float(np.interp(rr,r,v[1,0])),
                Gamma=float(np.interp(rr,r,v[0,1])),Gamma_metric=float(np.interp(rr,r,v[0,4])),CGamma=float(np.interp(rr,r,v[0,5])),
                Bdriver=float(np.interp(rr,r,q[0,2])),beta_P6_P8=float(np.interp(rr,r,q8[0,0]-q[0,0])),
                beta_r_P6_P8=float(np.interp(rr,r,q8[1,0]-q[1,0])),condition='T0_CURRENT_VALUES;R0_NOT_AVAILABLE'))
    t8.save('t9-initial-probes.csv',probes)
    for t,ray in itertools.product(TIMES[1:],t8.RAYS):
        ds=[np.load(TMP/f'{leg}-{t:g}-{ray}.npz') for leg in t8.LEGS];r=ds[0]['r'];common=abs(r-t)<=.2
        for leg,d in zip(t8.LEGS,ds):
            vv=d['values'];rf=r[int(d['front_index'])];use=abs(r-rf)<=.1
            features=[vv[0,i]-background(r,vv[0,i],t) for i in (1,4,5)]
            norms=[float(np.sqrt(np.mean(f[use]**2))) for f in features]
            gamma.append(dict(leg=leg,time_M=t,ray=ray,Gamma_feature_RMS=norms[0],metric_Gamma_feature_RMS=norms[1],CGamma_feature_RMS=norms[2],
                CGamma_over_Gamma_RMS=norms[2]/norms[0],metric_over_Gamma_RMS=norms[1]/norms[0],
                Gamma_metric_profile_correlation=float(np.corrcoef(features[0][use],features[1][use])[0,1]),
                CGamma_current_RMS=float(np.sqrt(np.mean(vv[0,5,use]**2))),CGamma_current_peak=float(max(abs(vv[0,5,use])))))
            for fi,f in enumerate(FIELDS[:6]):
                for degree in (3,5,7):
                    side=(abs(r-t)>=.4)&(abs(r-t)<=.6)
                    fit=np.polynomial.Polynomial.fit(r[side]-t,vv[0,fi,side],degree)
                    feat=vv[0,fi]-fit(r-t)
                    sensitivity.append(dict(leg=leg,time_M=t,ray=ray,field=f,side_fit_degree=degree,feature_peak_delta0p1=float(max(abs(feat[use]))),
                        feature_RMS_delta0p1=float(np.sqrt(np.mean(feat[use]**2))),condition='SAME_CURRENT_PROFILE_FIXED_SIDEBANDS_NO_STATIC_TARGET'))
        for fi,f in enumerate(FIELDS[:6]):
            dd=[float(np.sqrt(np.mean((ds[i]['values'][0,fi,common]-ds[i+1]['values'][0,fi,common])**2))) for i in (0,1)]
            ss=[float(np.sqrt(np.mean((d['values8'][0,fi,common]-d['values'][0,fi,common])**2))) for d in ds]
            differences.append(dict(time_M=t,ray=ray,field=f,low_mid_RMS=dd[0],mid_high_RMS=dd[1],
                self_order=math.log(dd[0]/dd[1])/math.log(1.5) if min(dd)>0 else math.nan,
                maximum_P6_P8_RMS=max(ss),interpolation_over_smaller_difference=max(ss)/min(dd) if min(dd)>0 else math.nan,
                condition='UNALIGNED_COMMON_PHYSICAL_WINDOW_TIME_PLUS_MINUS_0P2'))
    t8.save('t9-gamma-content.csv',gamma);t8.save('t9-profile-differences.csv',differences);t8.save('t9-background-sensitivity.csv',sensitivity)
    summary=[]
    for t,ray,delta,f in itertools.product(TIMES[1:],t8.RAYS,(.1,.2,.4),('beta','Gamma','Wout')):
        rr=[next(j for j in jumps if j['leg']==leg and j['time_M']==t and j['ray']==ray and j['delta_M']==delta and j['field']==f) for leg in t8.LEGS]
        summary.append(dict(time_M=t,ray=ray,delta_M=delta,field=f,low_J=rr[0]['J_endpoint'],mid_J=rr[1]['J_endpoint'],high_J=rr[2]['J_endpoint'],
            max_abs_P6_P8=max(abs(j['J_P8_minus_P6']) for j in rr),same_valid_level_all=all(j['same_valid_level'] for j in rr),
            low_positive=rr[0]['integral_positive'],mid_positive=rr[1]['integral_positive'],high_positive=rr[2]['integral_positive'],
            low_negative=rr[0]['integral_negative'],mid_negative=rr[1]['integral_negative'],high_negative=rr[2]['integral_negative']))
    t8.save('t9-jump-summary.csv',summary)
    qualification=[]
    for leg in t8.LEGS:
        for meta in t8.read(TMP/f'{leg}-metadata.csv'):
            t=float(meta['time_M']);ray=meta['ray']
            if not t:continue
            spread=abs(float(meta['width8_M'])/float(meta['width_M'])-1)
            # No retrospective pass threshold: expose measured spreads and geometry.
            qualification.append(dict(leg=leg,time_M=t,ray=ray,width_P6_P8_relative_spread=spread,
                phase_P6_P8_M=abs(float(meta['front8_M'])-float(meta['front_M'])),
                axial_face_overlap=(ray!='corner' and t in (.875,1.75)),
                full_delta0p4_single_valid_level=len(meta['levels_in_delta0p4'].split(','))==1,
                condition='FACE_OVERLAP_UNQUALIFIED_FRONT_WIDTH' if ray!='corner' and t in (.875,1.75) else 'MEASURED_SPREAD_RECORDED_NOT_CONTINUUM_ORDER_CERTIFICATE'))
    t8.save('t9-sampling-qualification.csv',qualification)

def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    import scienceplots
    plt.style.use(['science','no-latex'])
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'legend.fontsize':6,'figure.dpi':120})
    directory=HERE/'figures';directory.mkdir(exist_ok=True)
    colors=('#2166ac','#d97706','#20864a')
    def emit(fig,name):
        fig.savefig(directory/(name+'.png'),dpi=190);fig.savefig(directory/(name+'.pdf'));plt.close(fig)
    for t in (1.75,2.625,4.375):
        fig,axs=plt.subplots(3,3,figsize=(10,7.3),layout='constrained')
        for ri,ray in enumerate(t8.RAYS):
            for leg,color in zip(t8.LEGS,colors):
                d=np.load(TMP/f'{leg}-{t:g}-{ray}.npz');r=d['r'];v=d['values'];vv=d['values8']
                for fi in range(3):
                    ax=axs[ri,fi];ax.plot(r,v[0,fi],color=color,label=leg[2:]);ax.fill_between(r,v[0,fi],vv[0,fi],color=color,alpha=.15)
                    if t==1.75 and ri<2:ax.plot(r,d['parent'][0,fi],color=color,ls=':',lw=.6)
            for fi in range(3):
                ax=axs[ri,fi];ax.set(xlim=(t-.3,t+.3),title=f'{ray}: {FIELDS[fi]}',xlabel='r / M');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
                if ri<2 and t==1.75:ax.axvline(1.75,color='.5',ls='--',lw=.6)
        axs[0,0].legend(ncol=3);fig.suptitle(f'Signed current profiles, t={t:g} M; P6/P8 shading'+('; dotted: current single-parent level 5' if t==1.75 else ''))
        emit(fig,f't9-signed-profiles-{str(t).replace(".","p")}')
    fig,axs=plt.subplots(3,3,figsize=(10,7.3),layout='constrained')
    figc,axc=plt.subplots(3,3,figsize=(10,7.3),layout='constrained')
    for ti,t in enumerate((1.75,2.625,4.375)):
        for ri,ray in enumerate(t8.RAYS):
            for leg,color in zip(t8.LEGS,colors):
                d=np.load(TMP/f'{leg}-{t:g}-{ray}.npz');r=d['r'];v=d['values'];vv=d['values8']
                for fi,style,label in ((1,'-','Gamma'),(4,'--','metric Gamma')):
                    z=v[0,fi]-background(r,v[0,fi],t)
                    axs[ri,ti].plot(r,z,color=color,ls=style,label=leg[2:]+' '+label)
                axc[ri,ti].plot(r,v[0,5],color=color,label=leg[2:]);axc[ri,ti].fill_between(r,v[0,5],vv[0,5],color=color,alpha=.2)
                if t==1.75 and ri<2:axc[ri,ti].plot(r,d['parent'][0,5],color=color,ls=':',lw=.6)
            for ax in (axs[ri,ti],axc[ri,ti]):
                ax.set(xlim=(t-.2,t+.2),title=f'{ray}, t={t:g}',xlabel='r / M');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
                for l in range(1,13):
                    face=112/2**l*(math.sqrt(2) if ri==2 else 1)
                    if t-.2<face<t+.2:ax.axvline(face,color='.6',lw=.6,ls=':')
    axs[0,0].legend(ncol=2);fig.suptitle('Signed front features: evolved and native metric Gamma; same fixed sidebands')
    axc[0,0].legend(ncol=3);figc.suptitle('Native conformal Gamma constraint; P6/P8 shading, fixed physical intervals')
    emit(fig,'t9-metric-gamma-features');emit(figc,'t9-gamma-constraint')
    jj=t8.read(HERE/'t9-jump-integrals.csv')
    for t in (1.75,2.625,4.375):
        fig,axs=plt.subplots(3,3,figsize=(10,7),layout='constrained')
        for fi,f in enumerate(('beta','Gamma','Wout')):
            for di,delta in enumerate((.1,.2,.4)):
                ax=axs[fi,di]
                for ray,col,style in zip(t8.RAYS,colors,('-','--',':')):
                    rr=[next(r for r in jj if float(r['time_M'])==t and r['ray']==ray and r['field']==f and float(r['delta_M'])==delta and r['leg']==leg) for leg in t8.LEGS]
                    h=np.array([float(r['h_M']) for r in rr]);j=np.array([float(r['J_endpoint']) for r in rr]);j8=np.array([float(r['J_P8']) for r in rr])
                    ax.plot(h,j,style,color=col,marker='o',ms=3,label=ray);ax.fill_between(h,j,j8,color=col,alpha=.2)
                ax.set(xscale='log',xlabel='native h / M',title=f'{f}, delta={delta:g} M');ax.axhline(0,color='.6',lw=.5);ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
                ticks=sorted({float(r['h_M']) for r in jj if float(r['time_M'])==t})
                ax.xaxis.set_major_locator(FixedLocator(ticks));ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:.4f}'))
                ax.xaxis.set_minor_formatter(NullFormatter());ax.tick_params(axis='x',labelsize=6)
        axs[0,0].legend(ncol=3);fig.suptitle(f'Signed jump integrals, t={t:g} M; C2 reconstruction of current P6/P8 samples')
        emit(fig,f't9-jumps-{str(t).replace(".","p")}')
    fig,axs=plt.subplots(1,3,figsize=(10,3),layout='constrained')
    for leg,col in zip(t8.LEGS,colors):
        d=np.load(TMP/f'{leg}-0-axis.npz');r=d['r'];v=d['values']
        axs[0].plot(r,v[0,0],color=col,label=leg[2:]);axs[1].plot(r,v[1,0],color=col)
        d=np.load(TMP/f'{leg}-0.875-corner.npz');axs[2].plot(d['r'],d['values'][0,1],color=col)
    axs[0].set(xscale='log',xlabel='r / M',title='t=0, axis: beta');axs[1].set(xscale='log',xlabel='r / M',title='t=0, axis: beta derivative')
    axs[2].set(xlim=(.72,1.02),xlabel='r / M',title='t=0.875, diagonal: Gamma');axs[0].legend(ncol=3)
    fig.suptitle('Launch is unsampled between t=0 and 0.875 M; r=0 is not a native cell centre')
    emit(fig,'t9-launch')

def check():
    worst=[0.,0.,0.,0.]
    for n in (6,8):
        q=np.array([-.5,.25,2.25,n-.75]);first=np.zeros(len(q));w=basis(q,first,n)
        for deg in range(n):
            for d in range(min(3,deg)+1):
                got=w[d]@np.arange(n,dtype=float)**deg
                target=math.factorial(deg)/math.factorial(deg-d)*q**(deg-d)
                err=max(abs(got-target))/max(1,max(abs(target)));worst[d]=max(worst[d],err)
        assert np.max(abs(w[0]-t8.t6.t5.weights(q,first,n)))<1e-13
    assert max(worst)<1e-10,worst
    h=.125;lo=np.array([2672,0]);xx,yy=np.meshgrid((np.arange(32)+lo[0]+.5)*h-336,(np.arange(16)+.5)*h)
    a=np.zeros((len(ALL),16,32));radius2=xx*xx+yy*yy
    a[0]=xx*radius2;a[1]=yy*radius2;a[2]=xx*radius2**2;a[3]=yy*radius2**2;a[4]=xx;a[5]=yy
    a[6]=a[7]=a[10]=a[12]=1
    grid=[dict(a=a,h=h,lo=lo,bounds=(-2,2,0,2))];rr=np.array([.25,.33,.52,.75]);end_to_end=0.
    for n,vec in itertools.product((6,8),t8.VECS):
        got,_=sample(grid,rr,vec,n)
        for fi,degree in ((0,3),(1,5),(2,1)):
            for derivative in range(4):
                expected=math.factorial(degree)/math.factorial(degree-derivative)*rr**(degree-derivative) if derivative<=degree else np.zeros_like(rr)
                error=float(max(abs(got[derivative,fi]-expected)))/max(1.,max(abs(expected)))
                end_to_end=max(end_to_end,error)
                assert error<1e-10,(n,vec,fi,derivative,error)
    cubic=CubicSpline(np.linspace(-1,1,17),np.linspace(-1,1,17)**3)
    positive,negative=signed_integrals(cubic,-1,1)
    assert abs(positive-3)<1e-13 and abs(negative+3)<1e-13
    # Native cartoon contraction check against the sealed independent reducer.
    metric_error=0.
    for leg in t8.LEGS:
        path=sorted((t8.DATA/leg/'plt').glob('*.hdf5'))[3]
        with h5py.File(path) as f:
            l=5;g=f[f'level_{l}'];h=float(g.attrs['dx']);gx,gy=map(int,g['data_attributes'].attrs['outputGhost']);cs=audit.names(f)
            key,a=next(audit.blocks(f,l,True));v,q=audit.quantities(a,cs,h,key,gx,gy)
            # The loaded derived metric must reproduce Z=chi*C_Gamma/2 exactly to roundoff.
            _,levels=load(path);arr=levels[l]['a'];lo=levels[l]['lo'];x0,y0,x1,y1=key
            for i in range(2):
                cg=arr[2+i,y0-lo[1]:y1-lo[1]+1,x0-lo[0]:x1-lo[0]+1]-arr[13+i,y0-lo[1]:y1-lo[1]+1,x0-lo[0]:x1-lo[0]+1]
                err=np.max(abs(cg-2*q[f'Z{i+1}']/v[cs.index('chi')]))/max(1.,np.max(abs(cg)))
                metric_error=max(metric_error,float(err))
                assert err<16*np.finfo(float).eps,err
            del levels
    t8.save('t9-sampling-checks.csv',[dict(check='POLYNOMIAL_DERIVATIVES_DEGREE_N_MINUS_1',derivative=d,max_relative_error=e,status='PASS') for d,e in enumerate(worst)]+[
        dict(check='DIRECTIONAL_TENSOR_POLYNOMIALS_PARITY_ALL_RAYS',derivative=3,max_relative_error=end_to_end,status='PASS'),
        dict(check='NATIVE_CARTOON_CGAMMA_VS_SEALED_Z',derivative=0,max_relative_error=metric_error,status='PASS')])
    print('CHECK PASS',worst)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('profiles','measure','check','figures'));p.add_argument('leg',nargs='?');a=p.parse_args()
    if a.action=='profiles':profiles(a.leg)
    else:globals()[a.action]()
