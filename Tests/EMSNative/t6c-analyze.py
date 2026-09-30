#!/usr/bin/env python3
"""Frozen T6 C states: independent events, complete ray fields, native strips."""
import csv
import hashlib
import importlib.util
import json
import os
import resource
import time
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import find_peaks, peak_widths

HERE = Path(__file__).resolve().parent
ROOT = Path('/private/tmp/ems-t6/evolution')
CACHE = Path('/private/tmp/ems-t6/t6c-rays.npz')
spec = importlib.util.spec_from_file_location('t5', HERE/'t5-analyze.py')
t5 = importlib.util.module_from_spec(spec);spec.loader.exec_module(t5)
spec = importlib.util.spec_from_file_location('check', HERE/'t6-check.py')
check = importlib.util.module_from_spec(spec);spec.loader.exec_module(check)
VARS = t5.EVOLVED_COMPONENTS+('Ham','Mom1','Mom2','GaussE','GaussB')
t5.COMPONENTS = VARS
ODD = [VARS.index(v) for v in ('h12','A12','Gamma2','shift2','B2','By','Bz','Ey','Ez','Mom2')]
RAYS = ('axis_plus','axis_minus','equator','diagonal')
VECTORS = np.array(((1,0),(-1,0),(0,1),(2**-.5,2**-.5)))
R = np.arange(.25,8.,1/192)
T = np.arange(25)*.25


def save(name, rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)


def sample(levels,xy,n):
    chosen=np.zeros(len(xy),int)
    for l,g in enumerate(levels):
        a,b,c,d=g['bounds'];chosen[(xy[:,0]>=a)&(xy[:,0]<=b)&(xy[:,1]>=c)&(xy[:,1]<=d)]=l
    out=np.empty((len(VARS),len(xy)))
    for l,g in enumerate(levels):
        use=chosen==l
        if not use.any():continue
        h,lo,a=g['h'],g['lo'],g['a'];assert lo[1]==0
        x=(xy[use,0]+256)/h-.5-lo[0];y=xy[use,1]/h-.5
        sx=np.clip(np.floor(x).astype(int)-(n//2-1),0,a.shape[2]-n)
        sy=np.minimum(np.floor(y).astype(int)-(n//2-1),a.shape[1]-n)
        ix=sx[:,None]+np.arange(n);raw=sy[:,None]+np.arange(n);iy=np.where(raw<0,-raw-1,raw)
        wx,wy=t5.weights(x,sx,n),t5.weights(y,sy,n)
        anchor=a[:,iy[:,n//2],ix[:,n//2]];v=anchor.copy()
        for j in range(n):
            for i in range(n):
                q=a[:,iy[:,j],ix[:,i]].copy();q[ODD]*=np.where(raw[:,j]<0,-1.,1.)
                v+=(q-anchor)*(wy[:,j]*wx[:,i])[None,:]
        out[:,use]=v
    return out,chosen


def selfcheck():
    x=np.arange(250,262)+.5-256;y=np.arange(12)+.5
    xy=np.array(((.125,0.),(-.625,.125),(1.375,.375)))
    g=dict(h=1.,lo=np.array((250,0)),bounds=(-6,6,0,12))
    even=[i for i in range(len(VARS)) if i not in ODD]
    for n in (6,8):
        for p in range(n):
            for q in range(n):
                data=(x[None,:]/8)**p*(y[:,None]/16)**q
                g['a']=np.zeros((len(VARS),12,12))
                chosen=ODD if q%2 else even;g['a'][chosen]=data
                got,_=sample([g],xy,n)
                expected=np.zeros_like(got);expected[chosen]=(xy[:,0]/8)**p*(xy[:,1]/16)**q
                assert np.max(abs(got-expected))<3e-13,(n,p,q)
    # A later prominence basin must index its growth locally, not globally.
    log=np.array((2,1,0,1,2,1,0,0,2,3,0),float)
    pk,pr=find_peaks(log,prominence=.5)
    lo=int(pr['left_bases'][-1]);k=pk[-1];d=np.diff(log[lo:k+1]);j=np.argmax(d)
    assert lo>0 and d[j]==2 and lo+j+1==8


def cache():
    start=time.monotonic();audit=[];norms=[];out=[];out8=[];local=[]
    for case in ('face2','face3'):
        plots=sorted((ROOT/case/'plt').glob('*.hdf5'));assert len(plots)==25
        lines=[];lines8=[];levels_at=[];initial=None
        for ti,p in enumerate(plots):
            t,levels=t5.read(p);assert abs(t-T[ti])<1e-10
            geometry=[g['boxes'] for g in levels]
            if initial is None:initial=geometry
            assert all(np.array_equal(a,b) for a,b in zip(initial,geometry))
            xy=(R[None,:,None]*VECTORS[:,None,:]).reshape(-1,2)
            q,l=sample(levels,xy,6);q8,_=sample(levels,xy,8)
            lines.append(q.reshape(len(VARS),4,-1).transpose(1,0,2))
            lines8.append(q8.reshape(len(VARS),4,-1).transpose(1,0,2));levels_at.append(l.reshape(4,-1))
            nodes,v,lev,_=t5.nodes(levels,maximum=8.)
            radius=np.hypot(nodes[:,0],nodes[:,1]);h=np.array([g['h'] for g in levels])[lev]
            w=2*np.pi*nodes[:,1]*h*h
            for name,lo,hi in (('horizon',.63593977642346233,1.2718795528469247),('near_hole',.1,2),('ring',2,4),('far',4,8)):
                use=(radius>=lo)&(radius<=hi)
                for namec,c in (('Ham',v[28]),('Mom',np.hypot(v[29],v[30])),('GaussE',v[31])):
                    norms.append(dict(case=case,time_M=t,mask=name,constraint=namec,cells=int(use.sum()),
                                      weighting='coordinate_volume',rms=np.sqrt(np.sum(w[use]*c[use]**2)/w[use].sum())))
            assert all(np.count_nonzero(g['a'][-1])==0 for g in levels)
            audit.append(dict(case=case,time_M=t,path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                              components=33,levels=7,GaussB_nonzero=0,
                              chi_min=min(g['a'][0].min() for g in levels),lapse_min=min(g['a'][13].min() for g in levels),
                              chi_floor_cells=sum(np.count_nonzero(g['a'][0]<=1e-12) for g in levels),
                              lapse_floor_cells=sum(np.count_nonzero(g['a'][13]<=1e-12) for g in levels),
                              driver_B_absmax=max(np.abs(g['a'][16:18]).max() for g in levels),
                              magnetic_B_absmax=max(np.abs(g['a'][21:24]).max() for g in levels)))
        out.append(lines);out8.append(lines8);local.append(levels_at)
    np.savez(CACHE,time=T,radius=R,values=out,values8=out8,level=local)
    save('t6c-input-audit.csv',audit);save('t6c-mask-norms.csv',norms)
    print('CACHE',time.monotonic()-start,'s',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'RSS bytes')


def strips():
    rows=[];events=[];audit=[]
    for case,face in (('face2',2.),('face3',3.)):
        paths=sorted((ROOT/case/'strips').glob('*.bin'));assert len(paths)==216
        for path in paths:
            t,h,lev,iv,v=check.strip(path)
            fields=v[:,3:];x,y=v[:,0],v[:,1]
            audit.append(dict(case=case,time_M=t,level=lev,cells=len(v),covered_cells=int(v[:,2].sum()),
                              file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            # Geometry-only ray collars; times searched over the entire run.
            for loc,a in (('old_2M',2.),('new_3M',3.),('parent_4M',4.)):
                for ray in RAYS:
                    if ray=='axis_plus':use=(x>0)&(y<2*h)&(abs(x-a)<2*h+1e-10)
                    elif ray=='axis_minus':use=(x<0)&(y<2*h)&(abs(x+a)<2*h+1e-10)
                    elif ray=='equator':use=(abs(x)<2*h)&(abs(y-a)<2*h+1e-10)
                    else:use=(abs(x-a)<2*h+1e-10)&(abs(y-a)<2*h+1e-10)
                    for covered in (0,1):
                        sel=use&(v[:,2]==covered)
                        if not sel.any():continue
                        for c in ('Ham','Theta'):
                            k=VARS.index(c);q=fields[sel,k];imax=np.argmax(abs(q));rad=np.hypot(x[sel],y[sel])
                            rows.append(dict(case=case,location=loc,face_M=a,ray=ray,level=lev,dx_M=h,covered=covered,
                                             time_M=t,cells=int(sel.sum()),field=c,abs_peak=abs(q[imax]),signed_peak=q[imax],
                                             peak_radius_M=rad[imax],rms=np.sqrt(np.mean(q*q))))
        groups={(r['location'],r['ray'],r['level'],r['covered'],r['field']) for r in rows if r['case']==case}
        for key in sorted(groups):
            q=sorted([r for r in rows if r['case']==case and (r['location'],r['ray'],r['level'],r['covered'],r['field'])==key],key=lambda r:r['time_M'])
            times=np.array([r['time_M'] for r in q]);amp=np.array([r['abs_peak'] for r in q])
            log=np.log10(np.maximum(amp,1e-30))
            peaks,prop=find_peaks(log,prominence=.5)
            for i,k in enumerate(peaks):
                # Global log-growth time for this independently found peak's basin.
                lo=int(prop['left_bases'][i]);growth=np.diff(log[lo:k+1]);gj=int(np.argmax(growth));j=lo+gj
                events.append(dict(case=case,location=key[0],ray=key[1],level=key[2],covered=key[3],field=key[4],
                                   peak_time_M=times[k],peak_radius_M=q[k]['peak_radius_M'],abs_peak=amp[k],
                                   prominence_decades=prop['prominences'][i],
                                   largest_growth_time_M=times[j+1],growth_factor=10**growth[gj],
                                   selection='ALL_LOCAL_LOG_PEAKS_PROMINENCE_GE_0.5_FULL_INTERVAL'))
    save('t6c-strip-audit.csv',audit);save('t6c-strip-profiles.csv',rows);save('t6c-events.csv',events)
    print('STRIPS',len(rows),'profiles',len(events),'independent peaks')


def dynamics():
    f=np.load(CACHE);tracks=[];fits=[];state=[];decay=[]
    for ci,case in enumerate(('face2','face3')):
        for method,values in (('P6',f['values']),('P8',f['values8'])):
            for ray,(nx,ny) in zip(RAYS,VECTORS):
                di=RAYS.index(ray);a=values[ci,:,di]
                gamma=nx*a[:,11]+ny*a[:,12];beta=nx*a[:,14]+ny*a[:,15]
                det=a[:,1]*a[:,3]-a[:,2]**2
                g=(nx*nx*a[:,3]-2*nx*ny*a[:,2]+ny*ny*a[:,1])/det
                assert np.all(det>0) and np.all(g>0) and np.all(a[:,0]>0) and np.all(a[:,13]>0)
                y=np.gradient(np.gradient(gamma,R,axis=-1),R,axis=-1)
                dt=np.gradient(a,T,axis=0)
                for sign in (1,-1):
                    line=[]
                    for ti,t in enumerate(T):
                        if t<.75:continue
                        peaks,pr=find_peaks(sign*y[ti],prominence=0)
                        good=np.where((R[peaks]>.5)&(R[peaks]<7.5))[0]
                        j=good[np.argmax(pr['prominences'][good])];k=peaks[j]
                        width=peak_widths(sign*y[ti],[k])[0][0]*(R[1]-R[0])
                        face=float(case[-1]);extent=R[k]*max(abs(nx),abs(ny))
                        ahead=extent<face-4/48
                        row=dict(case=case,method=method,ray=ray,lobe='positive' if sign==1 else 'negative',
                                 time_M=t,radius_M=R[k],curvature_peak=y[ti,k],prominence=pr['prominences'][j],width_M=width,
                                 width_fine_cells=48*width,width_receiving_cells=24*width,
                                 width_crosses_face=(R[k]+width/2)*max(abs(nx),abs(ny))>=face,
                                 ahead_of_face=ahead,level=int(f['level'][ci,ti,di,k]),
                                 lapse=a[ti,13,k],chi=a[ti,0,k],beta_n=beta[ti,k],h_inverse_nn=g[ti,k],
                                 lapse_speed=-beta[ti,k]+np.sqrt(1.8*a[ti,13,k]*a[ti,0,k]*g[ti,k]),
                                 shift_transverse_speed=-beta[ti,k]+np.sqrt(.75*g[ti,k]),
                                 shift_longitudinal_speed=-beta[ti,k]+np.sqrt(g[ti,k]),
                                 light_speed=-beta[ti,k]+a[ti,13,k]*np.sqrt(a[ti,0,k]*g[ti,k]),
                                 Gamma_n=gamma[ti,k],Gamma_t=-ny*a[ti,11,k]+nx*a[ti,12,k],
                                 beta_t=-ny*a[ti,14,k]+nx*a[ti,15,k])
                        line.append(row)
                        if method=='P6' and sign==-1:
                            state.append(dict(case=case,ray=ray,time_M=t,radius_M=R[k],
                                              **{name:a[ti,i,k] for i,name in enumerate(VARS)},
                                              **{'dt_'+name:dt[ti,i,k] for i,name in enumerate(VARS)}))
                    for i,row in enumerate(line):
                        before,after=line[max(i-1,0)],line[min(i+1,len(line)-1)]
                        row['measured_local_speed']=(after['radius_M']-before['radius_M'])/(after['time_M']-before['time_M'])
                    tracks.extend(line)
                    good=[q for q in line if q['ahead_of_face'] and not q['width_crosses_face']]
                    tt=np.array([q['time_M'] for q in good]);rr=np.array([q['radius_M'] for q in good])
                    if len(tt)>=3:
                        (v,b),cov=np.polyfit(tt,rr,1,cov=True)
                        errors={name:np.sqrt(np.mean([(q['measured_local_speed']-q[name])**2 for q in good[1:-1]])) for name in ('lapse_speed','shift_transverse_speed','shift_longitudinal_speed','light_speed')}
                        fits.append(dict(case=case,method=method,ray=ray,lobe='positive' if sign==1 else 'negative',
                                         samples=len(tt),time_start_M=tt[0],time_end_M=tt[-1],speed=v,
                                         fit_standard_error=np.sqrt(cov[0,0]),radius_fit_rms_M=np.sqrt(np.mean((rr-v*tt-b)**2)),
                                         **{'local_speed_mismatch_rms_'+k:v for k,v in errors.items()}))
                for ti,t in enumerate(T[1:-1],1):
                    use=(R>1)&(R<2.5)
                    bn=nx*a[ti,16]+ny*a[ti,17];dbn=nx*dt[ti,16]+ny*dt[ti,17]
                    ratio=dbn[use]/bn[use]
                    decay.append(dict(case=case,method=method,ray=ray,time_M=t,
                                      median_dt_B_over_B=float(np.median(ratio)),max_deviation_from_minus_p1=float(np.max(abs(ratio+.1)))))
    save('t6c-front-tracks.csv',tracks);save('t6c-front-fits.csv',fits)
    save('t6c-field-at-front.csv',state);save('t6c-driver-decay.csv',decay)
    print('DYNAMICS',len(tracks),'tracks',len(fits),'pre-face fits')


def native():
    bands=[];widths=[];replay=[]
    w=t5.weights(np.array([.5]),-2,6)[0]
    assert np.array_equal(w,np.array([3,-25,150,150,-25,3])/256)
    for case in ('face2','face3'):
        for p in sorted((ROOT/case/'plt').glob('*.hdf5')):
            t,levels=t5.read(p)
            xy,v,lev,_=t5.nodes(levels,maximum=8.)
            x,y=xy.T;collar=1/24
            for location,a in (('old_2M',2.),('new_3M',3.),('parent_4M',4.),('new_seam_1p5',1.5)):
                for ray in RAYS:
                    if ray=='axis_plus':use=(x>0)&(y<collar)&(abs(x-a)<collar)
                    elif ray=='axis_minus':use=(x<0)&(y<collar)&(abs(x+a)<collar)
                    elif ray=='equator':use=(abs(x)<collar)&(abs(y-a)<collar)
                    else:use=(abs(x-a)<collar)&(abs(y-a)<collar)
                    if not use.any():continue
                    for field,c in (('Ham',28),('Theta',10)):
                        k=np.argmax(abs(v[c,use]));r=np.hypot(x[use],y[use]);q=v[c,use]
                        bands.append(dict(case=case,location=location,ray=ray,time_M=t,field=field,cells=int(use.sum()),
                                          peak_abs=abs(q[k]),peak_radius_M=r[k],rms=np.sqrt(np.mean(q*q)),
                                          finest_level=int(lev[use].max()),coarsest_level=int(lev[use].min())))
            if t>0:
                coarse,fine=levels[5],levels[6];h=coarse['h'];a=fine['bounds'][1]
                # Inner second covered coarse row: all six fine nodes are valid.
                for ray in ('axis_plus','equator','diagonal'):
                    target=np.array((a-1.5*h,4.5*h)) if ray=='axis_plus' else np.array((4.5*h,a-1.5*h)) if ray=='equator' else np.array((a-1.5*h,a-1.5*h))
                    iv=np.rint((target+np.array((256,0)))/h-.5).astype(int)
                    fi=2*iv-fine['lo'];co=iv-coarse['lo']
                    assert np.all(fi-2>=0) and np.all(fi+3<=fine['hi']-fine['lo'])
                    base=fine['a'][:28,fi[1],fi[0]];corr=np.zeros(28)
                    for j in range(6):
                        for i in range(6):corr+=w[i]*w[j]*(fine['a'][:28,fi[1]+j-2,fi[0]+i-2]-base)
                    actual=coarse['a'][:28,co[1],co[0]]
                    delta=actual-(base+corr)
                    replay.append(dict(case=case,time_M=t,ray=ray,coarse_x_M=(iv[0]+.5)*h-256,
                                       coarse_y_M=(iv[1]+.5)*h,max_abs_replay_difference=np.max(abs(delta)),
                                       max_scaled_eps_difference=np.max(abs(delta)/(np.finfo(float).eps*np.maximum(1,abs(actual)))),
                                       bit_differences=np.count_nonzero(actual.view('u8')!=(base+corr).view('u8')),
                                       condition='CURRENT_FULLY_VALID_FINE_STENCIL_SECOND_COVERED_COARSE_ROW'))
            for lev in (5,6):
                g=levels[lev];h=g['h'];r=(np.arange(g['lo'][0],g['hi'][0]+1)+.5)*h-256
                gamma=g['a'][11,0];d2=(-gamma[4:]+16*gamma[3:-1]-30*gamma[2:-2]+16*gamma[1:-3]-gamma[:-4])/(12*h*h)
                rr=r[2:-2]
                for sign in (1,-1):
                    pk,pr=find_peaks(sign*d2,prominence=0)
                    good=np.where((rr[pk]>.5)&(rr[pk]<float(case[-1])-.1))[0]
                    if not len(good):continue
                    j=good[np.argmax(pr['prominences'][good])];k=pk[j]
                    width=peak_widths(sign*d2,[k])[0][0]*h
                    widths.append(dict(case=case,time_M=t,level=lev,dx_M=h,lobe='positive' if sign==1 else 'negative',
                                       radius_M=rr[k],curvature_peak=d2[k],width_M=width,width_cells=width/h,
                                       prominence_baseline_at_domain_end=pr['right_bases'][j]==len(d2)-1 or pr['left_bases'][j]==0,
                                       condition='NATIVE_FIRST_Y_ROW_FOURTH_ORDER_D2_LINEAR_HALF_PROMINENCE_CROSSINGS'))
    save('t6c-native-bands.csv',bands);save('t6c-native-widths.csv',widths);save('t6c-restriction-replay.csv',replay)
    assert max(r['max_scaled_eps_difference'] for r in replay)<16
    print('NATIVE',len(bands),'bands,',len(replay),'restriction replays, max scaled eps',max(r['max_scaled_eps_difference'] for r in replay))


def carriers():
    f=np.load(CACHE);rows=[]
    for ci,case in enumerate(('face2','face3')):
        for di,ray in enumerate(RAYS):
            a=f['values'][ci,:,di]
            curvature=np.gradient(np.gradient(a,R,axis=-1),R,axis=-1)
            for vi,name in enumerate(VARS):
                for ti,t in enumerate(T):
                    for sign in (-1,1):
                        q=sign*curvature[ti,vi]
                        pk,pr=find_peaks(q,prominence=0)
                        good=np.where((R[pk]>.5)&(R[pk]<7.5)&(pr['prominences']>1e-18))[0]
                        if not len(good):continue
                        j=good[np.argmax(pr['prominences'][good])];k=pk[j]
                        width=peak_widths(q,[k])[0][0]*(R[1]-R[0])
                        rows.append(dict(case=case,ray=ray,field=name,time_M=t,sign=sign,radius_M=R[k],
                                         curvature_peak=curvature[ti,vi,k],prominence=pr['prominences'][j],
                                         width_M=width,level=int(f['level'][ci,ti,di,k]),
                                         condition='GLOBAL_STRONGEST_SPATIAL_CURVATURE_PROMINENCE_P6_NO_SPEED_WINDOW'))
    save('t6c-field-curvature-tracks.csv',rows)
    print('CARRIERS',len(rows),'independent component curvature peaks')


def report():
    os.environ.setdefault('MPLCONFIGDIR','/Users/auroradysis/.cache/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scienceplots
    plt.style.use(['science','no-latex'])
    plt.rcParams.update({'font.size':9,'axes.labelsize':9})
    f=np.load(CACHE);a=f['values'];events=list(csv.DictReader((HERE/'t6c-events.csv').open()))
    tracks=list(csv.DictReader((HERE/'t6c-front-tracks.csv').open()))
    figures=HERE/'figures';figures.mkdir(exist_ok=True)
    def finish(fig,name):
        for suffix in ('png','pdf'):fig.savefig(figures/(name+'.'+suffix),dpi=240)
        plt.close(fig)
    pairs=[]
    for ray in ('axis_plus','axis_minus','equator','diagonal'):
        for field in ('Ham','Theta'):
            pair=[]
            for case,loc in (('face2','old_2M'),('face3','new_3M')):
                q=[e for e in events if e['case']==case and e['location']==loc and e['ray']==ray
                   and e['field']==field and e['level']=='6' and e['covered']=='0']
                e=max(q,key=lambda z:float(z['abs_peak']));pair.append(e)
            pairs.append(dict(ray=ray,field=field,face2_peak_time_M=pair[0]['peak_time_M'],
                              face3_peak_time_M=pair[1]['peak_time_M'],
                              delay_M=float(pair[1]['peak_time_M'])-float(pair[0]['peak_time_M']),
                              face2_abs_peak=pair[0]['abs_peak'],face3_abs_peak=pair[1]['abs_peak'],
                              condition='LARGEST_AMPLITUDE_AMONG_ALL_INDEPENDENT_FULL_INTERVAL_LOCAL_PEAKS'))
    save('t6c-event-pairs.csv',pairs)
    fig,axs=plt.subplots(4,2,figsize=(8,9),layout='constrained')
    for row,(field,ray,low,high) in enumerate((('Ham','axis_plus',-9,-3),('Ham','diagonal',-9,-3),
                                             ('Theta','axis_plus',-10,-5),('Theta','diagonal',-10,-5))):
        for ci,case in enumerate(('face2','face3')):
            ax=axs[row,ci];q=a[ci,:,RAYS.index(ray),VARS.index(field)]
            im=ax.pcolormesh(R,T,np.log10(np.maximum(abs(q),1e-30)),shading='auto',
                             vmin=low,vmax=high,rasterized=True,cmap='magma')
            factor=2**.5 if ray=='diagonal' else 1
            ax.axvline(float(case[-1])*factor,color='cyan',ls='--',lw=.8)
            ax.axvline(4*factor,color='lime',ls=':',lw=.8)
            e=[e for e in events if e['case']==case and e['field']==field and e['ray']==ray
               and e['location']==('old_2M' if ci==0 else 'new_3M') and e['level']=='6' and e['covered']=='0']
            ax.scatter([float(z['peak_radius_M']) for z in e],[float(z['peak_time_M']) for z in e],
                       facecolors='none',edgecolors='white',s=25,lw=.8)
            ax.set(xlim=(.5,6.5),ylim=(0,6),ylabel='t / M',title=f'{case}: {field}, {ray.replace("_plus", "")}')
            if row==3:ax.set_xlabel('coordinate radius / M')
            fig.colorbar(im,ax=ax,label=f'log10 |{field}|',fraction=.04)
    finish(fig,'t6c-interface-events')
    fig,axs=plt.subplots(3,4,figsize=(11,8),layout='constrained')
    fields=a[1,:,0];dt=np.gradient(fields,T,axis=0)
    panels=[('Gamma1',np.gradient(np.gradient(fields[:,11],R,axis=-1),R,axis=-1)),
            ('shift1',np.gradient(np.gradient(fields[:,14],R,axis=-1),R,axis=-1))]
    panels += [(v,dt[:,VARS.index(v)]) for v in ('chi','lapse','K','B1','phi','Pi','Ex','Bz')]
    panels += [(v,fields[:,VARS.index(v)]) for v in ('Theta','GaussE')]
    assert all(q.shape==(25,len(R)) for _,q in panels)
    for ax,(name,q) in zip(axs.flat,panels):
        # Explicit single scale per panel, fixed over all radii and times.
        vmax=max(np.max(abs(q[:,(R>.5)&(R<6.5)])),1e-18);hi=np.ceil(np.log10(vmax));lo=hi-5
        im=ax.pcolormesh(R,T,np.log10(np.maximum(abs(q),1e-30)),shading='auto',vmin=lo,vmax=hi,
                         rasterized=True,cmap='magma')
        label='d2/dr2' if name in ('Gamma1','shift1') else ('value' if name in ('Theta','GaussE') else 'd/dt')
        ax.set(xlim=(.5,6.5),ylim=(0,6),title=f'{label} {name}',xlabel='r / M',ylabel='t / M')
        ax.axvline(3,color='cyan',ls='--',lw=.7);ax.axvline(4,color='lime',ls=':',lw=.7)
        fig.colorbar(im,ax=ax,label='log10 absolute',fraction=.04)
    fig.suptitle('Face3 axis: current fields and finite-cadence derivatives; individual panel scales')
    finish(fig,'t6c-incoming-fields')
    fig,axs=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
    q=[z for z in tracks if z['case']=='face3' and z['ray']=='axis_plus' and z['method']=='P6'
       and z['lobe']=='negative' and z['ahead_of_face']=='True' and z['width_crosses_face']=='False']
    tt=np.array([float(z['time_M']) for z in q]);rr=np.array([float(z['radius_M']) for z in q])
    assert len(tt)==10
    axs[0].plot(tt,rr,'o',label='independent Gamma curvature crest')
    axs[0].plot(tt,np.polyval(np.polyfit(tt,rr,1),tt),label='linear fit')
    axs[0].set(xlabel='t / M',ylabel='r / M',ylim=(.5,3.2));axs[0].legend(fontsize=7)
    for key,label in (('measured_local_speed','crest'),('shift_longitudinal_speed','longitudinal shift'),
                      ('shift_transverse_speed','transverse shift'),('lapse_speed','1+log lapse'),('light_speed','light')):
        axs[1].plot(tt[1:-1],[float(z[key]) for z in q[1:-1]],label=label)
    axs[1].set(xlabel='t / M',ylabel='outward coordinate speed',ylim=(0,1.1));axs[1].legend(fontsize=7)
    for method,style in (('P6','-'),('P8','--')):
        w=[z for z in tracks if z['case']=='face3' and z['ray']=='axis_plus' and z['method']==method
           and z['lobe']=='negative' and z['ahead_of_face']=='True' and z['width_crosses_face']=='False']
        for key,label in (('width_fine_cells','fine'),('width_receiving_cells','receiving coarse')):
            axs[2].plot([float(z['time_M']) for z in w],[float(z[key]) for z in w],style,label=f'{method} {label}')
    axs[2].set(xlabel='t / M',ylabel='half-prominence curvature width / h',ylim=(0,8));axs[2].legend(fontsize=7)
    finish(fig,'t6c-speed-width')
    print('REPORT',len(pairs),'event pairs; three PNG/PDF figures')


if __name__=='__main__':
    import sys
    selfcheck()
    if sys.argv[-1]=='cache':cache()
    elif sys.argv[-1]=='strips':strips()
    elif sys.argv[-1]=='dynamics':dynamics()
    elif sys.argv[-1]=='native':native()
    elif sys.argv[-1]=='carriers':carriers()
    elif sys.argv[-1]=='report':report()
    else:raise SystemExit('choose cache, strips, dynamics, native, carriers, or report')
