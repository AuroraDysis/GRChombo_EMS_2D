#!/usr/bin/env python3
"""Completed T7 A: current-field snapshots, independent events, operation replay."""
import csv,hashlib,importlib.util,json,os,re,resource,sys,time
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from scipy.signal import find_peaks,peak_widths
HERE=Path(__file__).resolve().parent; TMP=Path('/private/tmp/ems-t7a');TMP.mkdir(exist_ok=True)
ROOT=Path('/private/tmp/ems-t7/evolution')
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
check=module('t7check',HERE/'t7-check.py');t6=module('t6c',HERE/'t6c-analyze.py');t5=t6.t5
V=t6.VARS;R=np.arange(.5,6.5+1/768,1/384);TIMES=np.arange(73)/12
RAYS=('axis_plus','equator','diagonal');NV=np.array(((1,0),(0,1),(2**-.5,2**-.5)))
ODD=t6.ODD
FACES={4:8.,5:4.,6:3.}
def snapshot_time(t):
 q=round(t*12)/12
 assert abs(t-q)<1e-9, 'time outside the registered 1/12 M ladder'
 return q
def artifact(name):
 p=HERE/name
 return p if p.exists() else TMP/'tables'/name
def save(name,rows):
 assert rows,name
 directory=TMP/'tables' if re.search(r'-(clock|space)-L[456]\.csv$',name) else HERE
 directory.mkdir(exist_ok=True)
 with (directory/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for q in iter(lambda:f.read(1<<20),b''):h.update(q)
 return h.hexdigest()
def line_sample(iv,v,h,lev,n):
 out=np.full((3,33,len(R)),np.nan)
 for ri,ray in enumerate(RAYS):
  x=(iv[:,0]+.5)*h-256;y=(iv[:,1]+.5)*h
  use=(iv[:,1]==0)&(x>0) if ri==0 else (abs(x-h/2)<1e-10) if ri==1 else abs(x-y)<1e-10
  xx=(x if ri==0 else y if ri==1 else x*np.sqrt(2))[use];yy=v[use].T
  order=np.argsort(xx);xx=xx[order];yy=yy[:,order]
  if len(xx)<n:continue
  useq=(R>=xx[0]-(xx[1]-xx[0])/2)&(R<=xx[-1]+(xx[1]-xx[0])/2)
  step=h*(np.sqrt(2) if ri==2 else 1);q=(R[useq]-xx[0])/step
  st=np.clip(np.floor(q).astype(int)-n//2+1,0,len(xx)-n)
  weights=t5.weights(q,st,n);base=yy[:,st+n//2];ans=base.copy()
  for j in range(n):ans+=(yy[:,st+j]-base)*weights[:,j]
  out[ri,:,useq]=ans.T
 return out

def tensor_sample(iv,v,h,lev,n):
 # Qualified n x n current-field sampling only within complete saved face collars.
 lo=iv.min(0);hi=iv.max(0);a=np.full((33,hi[1]-lo[1]+1,hi[0]-lo[0]+1),np.nan)
 a[:,iv[:,1]-lo[1],iv[:,0]-lo[0]]=v.T
 out=np.full((3,33,len(R)),np.nan)
 for face in (3.,4.):
  for ri,vec in enumerate(NV):
   xy=R[:,None]*vec; ext=xy.max(1)
   use=(abs(ext-face)<.175)&(ext<FACES[lev])&(ext>0)
   if not use.any():continue
   # Restrict radial coordinates to the filled collar; valid union bounds on this level.
   low=np.ceil((face-.1875+256)/h-.5).astype(int)-lo[0];high=int(np.floor((min(face+.1875,FACES[lev])+256)/h-.5))-lo[0]
   lowy=max(0,int(np.ceil((face-.1875)/h-.5))-lo[1]);highy=int(np.floor(min(face+.1875,FACES[lev])/h-.5))-lo[1]
   if high-low+1<n:continue
   x=(xy[use,0]+256)/h-.5-lo[0];y=xy[use,1]/h-.5-lo[1]
   sx=np.floor(x).astype(int)-n//2+1;sy=np.floor(y).astype(int)-n//2+1
   if ri==0:sx=np.clip(sx,low,high-n+1)
   if ri==1:
    if highy-lowy+1<n:continue
    sy=np.clip(sy,lowy,highy-n+1)
   if ri==2:
    if highy-lowy+1<n:continue
    sx=np.clip(sx,low,high-n+1);sy=np.clip(sy,lowy,highy-n+1)
   sx=np.clip(sx,0,a.shape[2]-n)
   sy=np.minimum(sy,a.shape[1]-n)
   ix=sx[:,None]+np.arange(n);raw=sy[:,None]+np.arange(n);iy=np.where(raw<0,-raw-1,raw)
   if iy.max()>=a.shape[1] or iy.min()<0:continue
   wx=t5.weights(x,sx,n);wy=t5.weights(y,sy,n)
   base=a[:,iy[:,n//2],ix[:,n//2]];z=base.copy()
   for j in range(n):
    for i in range(n):
     q=a[:,iy[:,j],ix[:,i]].copy();q[ODD]*=np.where(raw[:,j]<0,-1.,1.)
     z+=(q-base)*(wx[:,i]*wy[:,j])[None,:]
   out[ri,:,use]=z.T
 return out

def metrics(case,t,h,lev,iv,v):
 x=(iv[:,0]+.5)*h-256;y=(iv[:,1]+.5)*h;extent=np.maximum(abs(x),y)
 valid=extent>FACES.get(lev+1,-1.)+1e-10
 rows=[]
 for face in (3.,4.):
  for ray in ('axis_plus','equator','diagonal','all_face'):
   raymask=(x>0)&(y<1/12) if ray=='axis_plus' else (abs(x)<1/12) if ray=='equator' else (abs(x-face)<1/12)&(abs(y-face)<1/12) if ray=='diagonal' else np.ones(len(x),bool)
   for band,mask in [('packet',abs(extent-face)<1/12),('wake_inner',(extent>face-.1875)&(extent<face-1/12)),('wake_outer',(extent>face+1/12)&(extent<face+.1875))]:
    sel=valid&raymask&mask
    if not sel.any():continue
    w=2*np.pi*y[sel]*h*h
    for name,q in [('Ham',v[:,28]),('Theta',v[:,10]),('Mom',np.hypot(v[:,29],v[:,30])),('GaussE',v[:,31])]:
     z=q[sel];j=np.argmax(abs(z))
     rows.append(dict(case=case,time_M=t,level=lev,face_M=face,ray=ray,band=band,field=name,cells=int(sel.sum()),weight=float(w.sum()),square_sum=float(np.dot(w,z*z)),rms=float(np.sqrt(np.dot(w,z*z)/w.sum())),peak_abs=float(abs(z[j])),peak_radius_M=float(np.hypot(x[sel][j],y[sel][j]))))
 return rows

def snapshots(case,lev):
 start=time.monotonic();p=ROOT/case/f't7-snapshot-L{lev}.xz';allrows=[];aud=[];vals=[];v8=[];ts=[];tensor=[];tensor8=[]
 group=[];current=None
 def consume(group):
  t=group[0]['meta'][0];h=group[0]['meta'][1];iv=np.concatenate([f['cells'] for f in group]);v=np.concatenate([f['values'] for f in group]);assert len(np.unique(iv,axis=0))==len(iv)
  aud.append(dict(case=case,level=lev,time_M=t,cells=len(iv),chi_min=v[:,0].min(),lapse_min=v[:,13].min(),chi_floor_cells=int(np.count_nonzero(v[:,0]<=1e-12)),lapse_floor_cells=int(np.count_nonzero(v[:,13]<=1e-12)),GaussB_nonzero=int(np.count_nonzero(v[:,32])),magnetic_B_absmax=np.max(abs(v[:,21:24])),finite=True))
  allrows.extend(metrics(case,t,h,lev,iv,v));vals.append(line_sample(iv,v,h,lev,6));v8.append(line_sample(iv,v,h,lev,8));tensor.append(tensor_sample(iv,v,h,lev,6));tensor8.append(tensor_sample(iv,v,h,lev,8));ts.append(t)
 for fr in check.frames(p):
  assert fr['phase']==50
  t=fr['meta'][0]
  if current is not None and abs(t-current)>1e-10:consume(group);group=[]
  current=t;group.append(fr)
 consume(group);assert len(ts)==73 and np.max(abs(np.array(ts)-TIMES))<1e-9
 np.savez(TMP/f'{case}-L{lev}.npz',time=ts,radius=R,values=vals,values8=v8,tensor=tensor,tensor8=tensor8)
 save(f't7a-snapshot-metrics-{case}-L{lev}.csv',allrows);save(f't7a-snapshot-audit-{case}-L{lev}.csv',aud)
 save(f't7a-performance-{case}-L{lev}.csv',[dict(case=case,level=lev,seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,file=str(p),sha256=digest(p),bytes=p.stat().st_size)])
 print('SNAPSHOTS',case,lev,len(allrows),time.monotonic()-start,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)

def baseline():
 rows=[];cache={lev:defaultdict(list) for lev in (4,5,6)}
 for p in sorted(Path('/private/tmp/ems-t6/evolution/face3/plt').glob('*.hdf5')):
  t,levels=t5.read(p)
  for lev in (4,5,6):
   g=levels[lev];xx,yy=np.meshgrid(np.arange(g['lo'][0],g['hi'][0]+1),np.arange(g['lo'][1],g['hi'][1]+1));iv=np.column_stack((xx.ravel(),yy.ravel()));v=g['a'].reshape(33,-1).T
   rows+=metrics('T6_face3',t,g['h'],lev,iv,v)
   cache[lev]['time'].append(t)
   for n,key in ((6,'values'),(8,'values8')):cache[lev][key].append(line_sample(iv,v,g['h'],lev,n))
   for n,key in ((6,'tensor'),(8,'tensor8')):cache[lev][key].append(tensor_sample(iv,v,g['h'],lev,n))
 for lev,q in cache.items():np.savez(TMP/f'T6_face3-L{lev}.npz',radius=R,**q)
 save('t7a-baseline-metrics.csv',rows);print('BASELINE',len(rows))


def combined():
 rows=[];events=[]
 for p in sorted(HERE.glob('t7a-snapshot-metrics-*.csv'))+sorted((TMP/'tables').glob('t7a-snapshot-metrics-*.csv'))+[HERE/'t7a-baseline-metrics.csv']:
  for r in csv.DictReader(p.open()):
   for k in ('time_M','face_M','weight','square_sum','rms','peak_abs','peak_radius_M'):r[k]=float(r[k])
   r['level']=int(r['level']);r['cells']=int(r['cells']);rows.append(r)
 save('t7a-native-metrics.csv',rows)
 groups=defaultdict(list)
 for r in rows:
  if r['band']=='packet' and r['ray']!='all_face' and r['field'] in ('Ham','Theta'):groups[(r['case'],r['face_M'],r['ray'],r['level'],r['field'])].append(r)
 for key,rr in sorted(groups.items()):
  rr.sort(key=lambda q:q['time_M']);ts=np.array([r['time_M'] for r in rr]);a=np.array([r['peak_abs'] for r in rr]);log=np.log10(np.maximum(a,1e-30));pk,pr=find_peaks(log,prominence=.5)
  candidates=[(int(k),int(pr['left_bases'][i]),float(pr['prominences'][i]),False) for i,k in enumerate(pk)]
  minima=find_peaks(-log)[0];lo=int(minima[-1]) if len(minima) else 0
  if log[-1]>log[-2] and log[-1]-log[lo]>=.5:candidates.append((len(a)-1,lo,float(log[-1]-log[lo]),True))
  for k,lo,prom,censored in candidates:
   eligible=np.where((a[lo:k]>=.01*a[k])&(a[lo+1:k+1]>=.01*a[k]))[0]
   if not len(eligible):continue
   j=lo+eligible[np.argmax(np.diff(log[lo:k+1])[eligible])]+1
   events.append(dict(case=key[0],face_M=key[1],ray=key[2],level=key[3],field=key[4],peak_time_M=ts[k],amplitude=a[k],growth_time_M=ts[j],growth_factor=10**(log[j]-log[j-1]),prominence_decades=prom,dominant=bool(a[k]==max(a[c[0]] for c in candidates)),censored=censored,selection='ALL_LOG_PEAKS_AND_RISING_ENDPOINTS_PROMINENCE_GE_0.5_FULL_INTERVAL'))
 save('t7a-events.csv',events)
 for r in events:
  if r['dominant'] and r['ray'] in ('axis_plus','diagonal'):print(r)
 # Composite volume norms in a fixed physical collar, independently of box decomposition.
 groups=defaultdict(list)
 for r in rows:groups[(r['case'],snapshot_time(r['time_M']),r['face_M'],r['ray'],r['band'],r['field'])].append(r)
 composite=[]
 for key,rr in sorted(groups.items()):
  w=sum(q['weight'] for q in rr);ss=sum(q['square_sum'] for q in rr)
  composite.append(dict(case=key[0],time_M=key[1],face_M=key[2],ray=key[3],band=key[4],field=key[5],cells=sum(q['cells'] for q in rr),volume_weight=w,rms=np.sqrt(ss/w),peak_abs=max(q['peak_abs'] for q in rr)))
 save('t7a-composite-metrics.csv',composite)
 ratios=[]
 for face in (3.,4.):
  for field in ('Ham','Theta','Mom','GaussE'):
   for band in ('packet','wake_inner','wake_outer'):
    for tt in (4.,6.):
     q={r['case']:r['rms'] for r in composite if r['face_M']==face and r['ray']=='all_face' and r['band']==band and r['field']==field and abs(r['time_M']-tt)<1e-10}
     if len(q)!=3:continue
     ratios.append(dict(face_M=face,field=field,band=band,time_M=tt,T6_face3=q['T6_face3'],clock=q['clock'],space=q['space'],clock_over_T6=q['clock']/q['T6_face3'],space_over_clock=q['space']/q['clock'],temporal_apparent_order=np.log2(q['T6_face3']/q['clock']),spatial_apparent_order=np.log2(q['clock']/q['space'])))
 save('t7a-ratios.csv',ratios)

def profiles():
 tracks=[];sensitivity=[];profilemetrics=[];seams=[]
 for case in ('T6_face3','clock','space'):
  fs={lev:np.load(TMP/f'{case}-L{lev}.npz') for lev in (4,5,6)};ts=fs[6]['time'];nr=len(R)
  combined_fields={}
  for mode in ('values','values8','tensor','tensor8'):
   q=np.full((len(ts),3,33,nr),np.nan)
   for lev in (4,5,6):
    for ri,vec in enumerate(NV):
     mask=(R*max(vec)<FACES[lev])&(R*max(vec)>=FACES.get(lev+1,0))
     q[:,ri,:,mask]=fs[lev][mode][:,ri,:,mask]
   combined_fields[mode]=q
  np.savez(TMP/f'{case}-rays.npz',time=ts,radius=R,**combined_fields)
  for ri,(ray,vec) in enumerate(zip(RAYS,NV)):
   for mode in ('values','values8'):
    a=combined_fields[mode][:,ri];gn=vec[0]*a[:,11]+vec[1]*a[:,12];yn=-np.gradient(np.gradient(gn,R,axis=-1),R,axis=-1)
    for ti,t in enumerate(ts):
     if t<.75:continue
     peaks,pr=find_peaks(yn[ti],prominence=0);good=np.where((R[peaks]>.75)&(R[peaks]*max(vec)<3.-1/12))[0]
     if not len(good):continue
     j=good[np.argmax(pr['prominences'][good])];k=peaks[j];width=peak_widths(yn[ti],[k])[0][0]/384
     # Reject endpoint/patch-clipped prominence basins from incoming-width claims.
     h6=(4/3 if case!='space' else 2/3)/64;h5=2*h6
     tracks.append(dict(case=case,method=mode,ray=ray,time_M=float(t),radius_M=R[k],Gamma_curvature_peak=-yn[ti,k],width_M=width,fine_cells=width/h6,receiving_cells=width/h5,prominence=pr['prominences'][j],basin_crosses_face=bool(R[int(pr['right_bases'][j])]*max(vec)>=3.),lapse=a[ti,13,k],chi=a[ti,0,k],shift_n=vec[0]*a[ti,14,k]+vec[1]*a[ti,15,k]))
   for ti,t in enumerate(ts):
    # All requested fields; pulse amplitudes here are variations across a common pre-face interval.
    if t not in (2.75,3.,4.,6.) and min(abs(t-q) for q in (2.75,3.,4.,6.))>1e-10:continue
    for mode in ('values','tensor'):
     a=combined_fields[mode][ti,ri];b=combined_fields[mode+'8'][ti,ri];use=(R*max(vec)>2.825)&(R*max(vec)<2.98)&np.all(np.isfinite(a),axis=0)&np.all(np.isfinite(b),axis=0)
     if not use.any():continue
     for k,name in enumerate(V):
      profilemetrics.append(dict(case=case,method=mode,ray=ray,time_M=float(t),field=name,samples=int(use.sum()),minimum=float(a[k,use].min()),maximum=float(a[k,use].max()),profile_range=float(np.ptp(a[k,use])),P6_P8_difference_rms=float(np.sqrt(np.mean((a[k,use]-b[k,use])**2))),P6_P8_difference_max=float(np.max(abs(a[k,use]-b[k,use])))))
   for seam in (1.5,2.25):
    use=(abs(R*max(vec)-seam)<1/12)
    for ti,t in enumerate(ts):
     for k,name in ((28,'Ham'),(10,'Theta')):
      z=combined_fields['values'][ti,ri,k,use]
      seams.append(dict(case=case,time_M=float(t),ray=ray,seam_coordinate_M=seam,field=name,rms=float(np.sqrt(np.mean(z*z))),peak_abs=float(np.max(abs(z))),region='COMMON_PHYSICAL_RAY_BAND_HALF_WIDTH_1_12'))
  for face in (3.,4.):
   for ri,ray in enumerate(RAYS):
    a=combined_fields['tensor'][:,ri];b=combined_fields['tensor8'][:,ri];use=(abs(R*max(NV[ri])-face)<.16)
    ok=np.all(np.isfinite(a),axis=1)&np.all(np.isfinite(b),axis=1)&use
    for k,name in enumerate(V):
     z=abs(a[:,k]-b[:,k]);good=ok
     if good.any():sensitivity.append(dict(case=case,face_M=face,ray=ray,field=name,qualified_samples=int(good.sum()),P6_P8_max=float(z[good].max()),P6_P8_RMS=float(np.sqrt(np.mean(z[good]**2))),missing_samples=int(use.sum()*len(ts)-good.sum())))
 save('t7a-incident-tracks.csv',tracks);save('t7a-profile-metrics.csv',profilemetrics);save('t7a-sampling-sensitivity.csv',sensitivity);save('t7a-common-seams.csv',seams)
 print('PROFILES',len(tracks),len(profilemetrics),len(sensitivity))

def window_map(case):
 ev=list(csv.DictReader((HERE/'t7a-events.csv').open()));out={}
 for face in (3.,4.):
  for ray in RAYS:
   rr=[r for r in ev if r['case']==case and float(r['face_M'])==face and r['ray']==ray and r['field']=='Ham' and r['dominant']=='True' and float(r['amplitude'])>1e-7]
   if not rr:continue
   out[(face,ray)]=(min(float(r['growth_time_M']) for r in rr)-.25,max(float(r['peak_time_M']) for r in rr)+1/12)
 return out

def source_region(lev,src):
 region=src%16
 return (4. if lev==4 or lev==5 and region<3 else 3.),RAYS[region%3]

def operations(case,lev):
 import subprocess,pickle,gzip
 start=time.monotonic();p=ROOT/case/f't7-stage-L{lev}.xz';windows=window_map(case);counts=Counter();previous={};states={};parts={};dense={};agg={};selected=[];native_meta=[];expected=[];native_index={};payload=TMP/f'{case}-L{lev}-stencils.bin';ghost_count=0;ghost_error=0.;dependency_count=0;stage_error=0.;fractions=set();projections=defaultdict(float)
 fields=(0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,24,25,27)
 def accumulate(fr,op,c,q):
  face,ray=source_region(lev,fr['source']);tt=fr['meta'][3] if op.startswith('RHS') else fr['meta'][0];key=(int(round(tt*24))/24,face,ray,op,V[c]);ab=abs(q);j=np.argmax(ab)
  if key not in agg:agg[key]=[0,0.,0.,0.,0.]
  z=agg[key];z[0]+=len(q);z[3]+=float(np.dot(q,q));z[4]+=float(q.sum())
  if ab[j]>z[1]:z[1]=float(ab[j]);z[2]=float(q[j])
 with payload.open('wb') as output:
  for fr in check.frames(p):
   counts[fr['phase']]+=1;m=fr['meta'];phase=fr['phase'];src=fr['source'];face,ray=source_region(lev,src);win=windows.get((face,ray));keep=win is not None and win[0]-1e-9<=m[3]<=win[1]+1e-9
   if phase in (2,3,20,21,22,40):
    assert fr['stage'] in range(4)
    stage_error=max(stage_error,abs(m[3]-m[0]-[0,.5,.5,1][fr['stage']]*m[2]))
    fractions.add((round(m[6],8),round((m[3]-m[4])/(m[5]-m[4]),8)))
   if keep:selected.append(fr)
   if phase in (2,4,6,8,9,12):previous[(src,phase)]=fr
   pairs={3:2,5:4,7:6,11:9,13:12}
   if phase in pairs:
    before=previous[(src,pairs[phase])];assert np.array_equal(before['cells'],fr['cells']);delta=fr['values']-before['values'];op={3:'stage_projection',5:'update_projection',7:'end_projection',11:'restriction',13:'ODE_increment'}[phase]
    for c in fields:
     projections[(op,V[c])]=max(projections[(op,V[c])],float(np.max(abs(delta[:,c]))))
     if keep:accumulate(fr,op,c,delta[:,c])
   if phase==40:dense[src]=(fr,{tuple(iv):v for iv,v in zip(fr['cells'],fr['values'])})
   if phase in (2,3):
    lookup={tuple(iv):v for iv,v in zip(fr['cells'],fr['values'])};states[(src,phase)]=(fr,lookup)
    if phase==2 and keep and ((lev==6) or (lev==5 and src%16<3)) and src in dense:
     parent,pl=dense[src];assert parent['stage']==fr['stage'] and parent['meta'][3]==m[3]
     for iv,val in zip(fr['cells'],fr['values']):
      x=(iv[0]+.5)*m[1]-256;y=(iv[1]+.5)*m[1]
      if y<0 or abs(x)<FACES[lev] and y<FACES[lev]:continue
      q=(iv+.5)/2-.5;first=np.floor(q).astype(int)-2;wx=t5.weights(np.array([q[0]]),first[0],6)[0];wy=t5.weights(np.array([q[1]]),first[1],6)[0];base=pl[tuple(first+2)];z=base.copy()
      for j in range(6):
       for i in range(6):z+=(pl[tuple(first+np.array((i,j)))]-base)*wx[i]*wy[j]
      ghost_error=max(ghost_error,float(np.max(abs(z-val)/np.maximum(1,abs(val)))));ghost_count+=1
    if keep:
     x0,y0,x1,y1=fr['valid'];h=m[1];core=[]
     # Core is the two-cell window immediately at this physical face/corner, intersected with this FAB.
     c=int(round(256/h));k=int(round(face/h))
     coords=[(c+k-1,0),(c+k,0)] if ray=='axis_plus' else [(i,j) for j in (k-1,k) for i in (c-1,c)] if ray=='equator' else [(i,j) for j in (k-1,k) for i in (c+k-1,c+k)]
     for iv in coords:
      if not(x0<=iv[0]<=x1 and y0<=iv[1]<=y1):continue
      assert iv in lookup
      a=np.zeros((28,7,7))
      for dy in range(-3,4):
       for dx in range(-3,4):
        needed=abs(dx)<=2 and abs(dy)<=2 or dy==0 and abs(dx)==3 or dx==0 and abs(dy)==3
        if needed:
         assert (iv[0]+dx,iv[1]+dy) in lookup;dependency_count+=1;a[:,dy+3,dx+3]=lookup[(iv[0]+dx,iv[1]+dy)]
      np.r_[h,(iv[1]+.5)*h,a.ravel()].tofile(output)
      idx=len(native_meta);native_meta.append(dict(case=case,level=lev,face_M=face,ray=ray,phase=phase,time_M=m[3],start_time_M=m[0],stage=fr['stage'],coarse_fraction=(m[3]-m[4])/(m[5]-m[4]),start_fraction=m[6],i=iv[0],j=iv[1],lapse=lookup[iv][13],Theta=lookup[iv][10]))
      expected.append(np.full(56,np.nan))
      if phase==3:native_index[(src,m[0],fr['stage'],iv)]=idx
   if phase in (20,21,22):
    if keep:
     opname={20:'RHS_pre_KO',21:'RHS_KO',22:'RHS_total'}[phase]
     for c in fields:accumulate(fr,opname,c,fr['values'][:,c])
     if phase in (20,21):
      for iv,val in zip(fr['cells'],fr['values']):
       idx=native_index[(src,m[0],fr['stage'],tuple(iv))];expected[idx][0 if phase==20 else 28:28 if phase==20 else 56]=val
    if phase in (20,21):parts[(src,phase)]=fr
    else:
     a,b=parts[(src,20)],parts[(src,21)];assert np.array_equal(a['cells'],fr['cells']);err=float(np.max(abs(a['values']+b['values']-fr['values'])));projections[('RHS_sum','ALL')]=max(projections[('RHS_sum','ALL')],err)
 rows=[]
 for (t,face,ray,op,field),z in sorted(agg.items()):rows.append(dict(case=case,level=lev,time_bin_M=t,face_M=face,ray=ray,operation=op,field=field,operands=z[0],abs_peak=z[1],signed_at_peak=z[2],RMS=np.sqrt(z[3]/z[0]),mean=z[4]/z[0],condition='INDEPENDENT_H_EVENT_WINDOW;RHS_UNITS_PER_M;STATE_DELTAS_PER_ACTUAL_OPERATION'))
 save(f't7a-operations-{case}-L{lev}.csv',rows)
 with gzip.open(TMP/f'{case}-L{lev}-events.pkl.gz','wb',compresslevel=3) as f:pickle.dump(selected,f,protocol=5)
 out=TMP/f'{case}-L{lev}-replay.bin';subprocess.run([str(TMP/'replay.ex'),str(payload),str(out)],check=True,timeout=100)
 actual=np.fromfile(out,dtype='f8').reshape(-1,57);assert len(actual)==len(native_meta) and np.isfinite(actual).all();expected=np.array(expected);selected3=np.array([r['phase']==3 for r in native_meta]);assert np.isfinite(expected[selected3]).all()
 error=float(np.max(abs(actual[selected3,:56]-expected[selected3])/np.maximum(1,abs(expected[selected3]))));assert error<1e-11,error
 nr=[]
 for meta,vals,ex in zip(native_meta,actual,expected):
  nr.append({**meta,'Ham_input':vals[56],'Theta_pre_KO':vals[10],'Theta_KO':vals[38],'Theta_total':vals[10]+vals[38],'half_alpha_Ham':.5*meta['lapse']*vals[56],'RHS_replay_error':float(np.nanmax(abs(vals[:56]-ex))) if meta['phase']==3 else 0.})
 save(f't7a-stage-native-{case}-L{lev}.csv',nr)
 save(f't7a-stage-audit-{case}-L{lev}.csv',[dict(case=case,level=lev,frames=sum(counts.values()),phase_counts=str(dict(counts)),stage_time_error=stage_error,coarse_fractions=str(sorted(fractions)),selected_frames=len(selected),replayed_core_states=len(actual),ghosts_replayed=ghost_count,ghost_replay_max_error=ghost_error,native_rhs_replay_scaled_error=error,stencil_dependencies=dependency_count,chi_floor_hits=int(projections[('stage_projection','chi')]>0 or projections[('end_projection','chi')]>0),lapse_floor_hits=int(projections[('stage_projection','lapse')]>0 or projections[('end_projection','lapse')]>0),finite=True,seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,sha256=digest(p),bytes=p.stat().st_size)])
 save(f't7a-global-operation-deltas-{case}-L{lev}.csv',[dict(case=case,level=lev,operation=op,field=field,abs_max=val) for (op,field),val in sorted(projections.items())])
 # Hash own bulky scratch before removal; retained event frames and CSVs remain.
 save(f't7a-scratch-hashes-{case}-L{lev}.csv',[dict(path=str(q),sha256=digest(q),bytes=q.stat().st_size,action='REMOVED_AFTER_NATIVE_REPLAY') for q in (payload,out)])
 payload.unlink();out.unlink();print('OPERATIONS',case,lev,len(nr),'seconds',time.monotonic()-start,'rss',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)

def native_widths(case,lev=6):
 rows=[]
 def consume(group):
  t=group[0]['meta'][0];h=group[0]['meta'][1];iv=np.concatenate([f['cells'] for f in group]);v=np.concatenate([f['values'] for f in group]);x=(iv[:,0]+.5)*h-256;y=(iv[:,1]+.5)*h
  for ri,(ray,vec) in enumerate(zip(RAYS,NV)):
   mask=(iv[:,1]==0)&(x>0) if ri==0 else (abs(x-h/2)<1e-10) if ri==1 else abs(x-y)<1e-10
   rr=(x if ri==0 else y if ri==1 else x*np.sqrt(2))[mask];vv=v[mask];order=np.argsort(rr);rr=rr[order];vv=vv[order];step=h*(np.sqrt(2) if ri==2 else 1.)
   gn=vec[0]*vv[:,11]+vec[1]*vv[:,12];beta=vec[0]*vv[:,14]+vec[1]*vv[:,15]
   curv=(-gn[4:]+16*gn[3:-1]-30*gn[2:-2]+16*gn[1:-3]-gn[:-4])/(12*step*step);r=rr[2:-2]
   for sign in (-1,1):
    pk,pr=find_peaks(sign*curv,prominence=0);good=np.where((r[pk]*max(vec)>(.75 if lev==6 else 3.2))&(r[pk]*max(vec)<FACES[lev]-.1))[0]
    if not len(good):continue
    j=good[np.argmax(pr['prominences'][good])];k=pk[j];width=peak_widths(sign*curv,[k])[0][0]*step;z=vv[k+2];det=z[1]*z[3]-z[2]**2;g=(vec[0]**2*z[3]-2*vec[0]*vec[1]*z[2]+vec[1]**2*z[1])/det;bn=beta[k+2]
    yy=sign*curv;half=yy[k]/2;left=np.where(yy[:k]<=half)[0];right=np.where(yy[k+1:]<=half)[0];closed=yy[k]>0 and len(left)>0 and len(right)>0;ww=np.nan
    if closed:
     l=left[-1];rgt=k+1+right[0];xl=r[l]+(half-yy[l])/(yy[l+1]-yy[l])*step;xr=r[rgt-1]+(half-yy[rgt-1])/(yy[rgt]-yy[rgt-1])*step;ww=xr-xl
    rows.append(dict(case=case,ray=ray,time_M=t,lobe='negative' if sign==-1 else 'positive',radius_M=r[k],width_M=width,fine_cells=width/h,receiving_cells=width/(2*h),halfheight_width_M=ww,halfheight_closed=bool(closed),halfheight_receiving_cells=ww/(2*h),curvature_peak=curv[k],prominence=pr['prominences'][j],domain_end_baseline=bool(pr['left_bases'][j]==0 or pr['right_bases'][j]==len(curv)-1),lapse_speed=-bn+np.sqrt(1.8*z[13]*z[0]*g),shift_transverse_speed=-bn+np.sqrt(.75*g),shift_longitudinal_speed=-bn+np.sqrt(g),light_speed=-bn+z[13]*np.sqrt(z[0]*g),**{name:z[V.index(name)] for name in ('lapse','K','Theta','chi','phi','Pi','Xi')},Gamma_n=gn[k+2],shift_n=bn,Bdriver_n=vec[0]*z[16]+vec[1]*z[17],electric_n=vec[0]*z[24]+vec[1]*z[25]))
 group=[];current=None
 for fr in check.frames(ROOT/case/f't7-snapshot-L{lev}.xz'):
  t=fr['meta'][0]
  if current is not None and abs(t-current)>1e-10:consume(group);group=[]
  current=t;group.append(fr)
 consume(group);save(f't7a-native-widths-{case}-L{lev}.csv',rows)
 for r in rows:
  if r['lobe']=='negative' and r['ray'] in ('axis_plus','diagonal') and abs(r['time_M']-2.5)<1e-9:print(r)

def rk_budget(case,lev):
 import gzip,pickle
 with gzip.open(TMP/f'{case}-L{lev}-events.pkl.gz','rb') as f:frames=pickle.load(f)
 steps={};rows=[]
 for fr in frames:
  p=fr['phase'];m=fr['meta'];src=fr['source'];idx=int(round(m[0]/m[2]))-(1 if p in (6,7) else 0);key=(src,idx)
  if p==3 and fr['stage']==0:
   steps.setdefault(key,{})['input']=fr
  if p in (20,21):steps.setdefault(key,{})[(p,fr['stage'])]=fr
  if p==7:steps.setdefault(key,{})['end']=fr
 for key,q in steps.items():
  if not all(k in q for k in ['input','end']+[(p,s) for p in (20,21) for s in range(4)]):continue
  end=q['end'];dt=end['meta'][2];lookup={tuple(iv):v for iv,v in zip(q['input']['cells'],q['input']['values'])};init=np.array([lookup[tuple(iv)][10] for iv in end['cells']]);actual=end['values'][:,10]-init
  phys=dt/6*sum(w*q[(20,s)]['values'][:,10] for s,w in enumerate((1,2,2,1)));ko=dt/6*sum(w*q[(21,s)]['values'][:,10] for s,w in enumerate((1,2,2,1)));error=float(np.max(abs(actual-phys-ko)));assert error<1e-18,error
  face,ray=source_region(lev,end['source']);k=np.argmax(abs(actual))
  rows.append(dict(case=case,level=lev,face_M=face,ray=ray,end_time_M=end['meta'][0],dt_M=dt,cells=len(init),Theta_before_at_peak=init[k],Theta_after_at_peak=end['values'][k,10],physical_increment_at_peak=phys[k],KO_increment_at_peak=ko[k],total_increment_at_peak=actual[k],physical_increment_absmax=float(np.max(abs(phys))),KO_increment_absmax=float(np.max(abs(ko))),total_increment_absmax=float(np.max(abs(actual))),RK4_replay_error=error))
 save(f't7a-rk-budget-{case}-L{lev}.csv',rows)
 print('RK',case,lev,len(rows),max(r['RK4_replay_error'] for r in rows))

def restriction(case,lev):
 import gzip,pickle
 def read(l):
  with gzip.open(TMP/f'{case}-L{l}-events.pkl.gz','rb') as f:return pickle.load(f)
 parent=read(lev);child=read(lev+1);supports=defaultdict(list);before={};rows=[]
 for fr in child:
  if fr['phase']==10:supports[int(round(fr['meta'][0]*1e9))].append(fr)
 w=np.array((3,-25,150,150,-25,3))/256
 for fr in parent:
  if fr['phase']==9:before[fr['source']]=fr
  if fr['phase']!=11:continue
  old=before.get(fr['source']);fine=supports.get(int(round(fr['meta'][0]*1e9)),[])
  if old is None or not fine:continue
  for k,(iv,actual) in enumerate(zip(fr['cells'],fr['values'])):
   h=fr['meta'][1];x=(iv[0]+.5)*h-256;y=(iv[1]+.5)*h
   if max(abs(x),y)>=FACES[lev+1]:continue
   anchor=2*iv;owner=next((q for q in fine if q['valid'][0]<=anchor[0]<=q['valid'][2] and q['valid'][1]<=anchor[1]<=q['valid'][3]),None)
   if owner is None:continue
   lookup={tuple(a):v for a,v in zip(owner['cells'],owner['values'])}
   if tuple(anchor) not in lookup:continue
   base=lookup[tuple(anchor)];z=base.copy()
   for j in range(6):
    for i in range(6):
     node=tuple(anchor+np.array((i-2,j-2)));assert node in lookup;z+=(lookup[node]-base)*w[i]*w[j]
   err=float(np.max(abs(z-actual)/np.maximum(1,abs(actual))));assert err<32*np.finfo(float).eps
   face,ray=source_region(lev,fr['source']);rows.append(dict(case=case,coarse_level=lev,face_M=face,ray=ray,time_M=fr['meta'][0],i=iv[0],j=iv[1],Theta_before=old['values'][k,10],Theta_after=actual[10],Theta_jump=actual[10]-old['values'][k,10],fine_Theta_support_absmax=float(max(abs(v[10]) for v in lookup.values())),replay_max_error=err,condition='COVERED_ONLY_CURRENT_FINE_VALID_AND_ACTUAL_GHOST_SUPPORT'))
 save(f't7a-restriction-{case}-L{lev}.csv',rows);print('RESTRICTION',case,lev,len(rows),max(r['replay_max_error'] for r in rows))

def common_packet_peaks():
 rows=[]
 def consume(case,lev,t,h,iv,v):
  x=(iv[:,0]+.5)*h-256;y=(iv[:,1]+.5)*h;ext=np.maximum(abs(x),y);valid=ext>FACES.get(lev+1,-1)+1e-10
  for face in (3.,4.):
   width=1/24 if face==3 else 1/12
   for ray in RAYS:
    mask=(x>0)&(y<width) if ray=='axis_plus' else abs(x)<width if ray=='equator' else (abs(x-face)<width)&(abs(y-face)<width)
    sel=valid&mask&(abs(ext-face)<width)
    if not sel.any():continue
    for name,k in (('Ham',28),('Theta',10)):rows.append(dict(case=case,time_M=t,face_M=face,ray=ray,level=lev,field=name,cells=int(sel.sum()),peak_abs=float(np.max(abs(v[sel,k]))),ROI_half_width_M=width))
 for p in sorted(Path('/private/tmp/ems-t6/evolution/face3/strips').glob('*.bin')):
  t,h,lev,iv,v=t6.check.strip(p);valid=v[:,2]==0;consume('T6_face3',lev,t,h,iv[valid],v[valid,3:])
 for case in ('clock','space'):
  for lev in (4,5,6):
   group=[];current=None
   def flush(group):consume(case,lev,group[0]['meta'][0],group[0]['meta'][1],np.concatenate([f['cells'] for f in group]),np.concatenate([f['values'] for f in group]))
   for fr in check.frames(ROOT/case/f't7-snapshot-L{lev}.xz'):
    t=fr['meta'][0]
    if current is not None and abs(t-current)>1e-10:flush(group);group=[]
    current=t;group.append(fr)
   flush(group)
 save('t7a-common-packet-history.csv',rows)
 ratios=[]
 for face in (3.,4.):
  for ray in RAYS:
   for lev in (int(8-face),int(9-face)):
    for field in ('Ham','Theta'):
     q={case:[r for r in rows if r['case']==case and r['face_M']==face and r['ray']==ray and r['level']==lev and r['field']==field] for case in ('T6_face3','clock','space')}
     if not all(q.values()):continue
     peaks={case:max(rr,key=lambda r:r['peak_abs']) for case,rr in q.items()};a,b,c=[peaks[k]['peak_abs'] for k in ('T6_face3','clock','space')]
     ratios.append(dict(face_M=face,ray=ray,level=lev,field=field,T6_peak=a,clock_peak=b,space_peak=c,T6_peak_time_M=peaks['T6_face3']['time_M'],clock_peak_time_M=peaks['clock']['time_M'],space_peak_time_M=peaks['space']['time_M'],clock_over_T6=b/a,space_over_clock=c/b,temporal_apparent_order=np.log2(a/b),spatial_apparent_order=np.log2(b/c),condition='SAME_PHYSICAL_ROI_AND_PRE_PARENT_RESTRICTION_SAMPLING_CADENCE;PEAK_SELECTED_FULL_INTERVAL'))
 save('t7a-common-packet-ratios.csv',ratios)
 for r in ratios:
  if r['level']==int(9-r['face_M']) and r['ray'] in ('axis_plus','diagonal'):print(r)

def summaries():
 audit=[];defaults=[]
 for case in ('clock','space'):
  d=ROOT/case;txt=(d/'run.log').read_text();params={}
  for line in (d/'params.txt').read_text().splitlines():
   if '=' in line:k,v=line.split('=',1);params[k.strip()]=v.strip()
  assert (d/'done.exit').read_text().strip()=='0' and params['sigma']=='1' and params['amr_transfer']=='point' and params['t2_guard_initial_data_after_t0']=='true'
  assert params['nan_check']=='1' and set(params['regrid_interval'].split())=={'0'}
  assert not any(s in txt.lower() for s in ('nancheck in','mayday::error','segmentation fault','nonfinite'))
  for line in txt.splitlines():
   if 'not found' in line:defaults.append(dict(case=case,message=line))
  ss=[r for lev in (4,5,6) for r in csv.DictReader(artifact(f't7a-snapshot-audit-{case}-L{lev}.csv').open())]
  st=[next(csv.DictReader(artifact(f't7a-stage-audit-{case}-L{lev}.csv').open())) for lev in (4,5,6)]
  audit.append(dict(case=case,exit=0,sigma=1,amr_transfer='point',fixed_hierarchy=True,static_guard=True,stage_frames=sum(int(r['frames']) for r in st),snapshots_per_level=73,missing_evolved_components=0,nonfinite='NONE_IN_ALL_DECODED_FRAMES;GLOBAL_RUNTIME_NAN_CHECK_ON',chi_min_captured=min(float(r['chi_min']) for r in ss),lapse_min_captured=min(float(r['lapse_min']) for r in ss),floor_cells_captured=sum(int(r['chi_floor_cells'])+int(r['lapse_floor_cells']) for r in ss),GaussB_nonzero=sum(int(r['GaussB_nonzero']) for r in ss),wall_seconds=float((d/'wall_seconds').read_text()),evolution_peak_rss_bytes=int((d/'peak_rss_bytes').read_text()),stage_time_error_max=max(float(r['stage_time_error']) for r in st),ghost_replay_error_max=max(float(r['ghost_replay_max_error']) for r in st),native_rhs_replay_error_max=max(float(r['native_rhs_replay_scaled_error']) for r in st),params_sha256=digest(d/'params.txt'),log_sha256=digest(d/'run.log')))
 save('t7a-run-audit.csv',audit);save('t7a-default-parameters.csv',defaults)
 source=[];deltas=[];native=[]
 for case in ('clock','space'):
  for lev in (4,5,6):
   rr=list(csv.DictReader(artifact(f't7a-stage-native-{case}-L{lev}.csv').open()));pp=defaultdict(dict)
   for r in rr:
    key=(r['face_M'],r['ray'],r['stage'],r['time_M'],r['i'],r['j']);pp[key][int(r['phase'])]=r
   for key,qr in pp.items():
    if len(qr)==2:deltas.append(dict(case=case,level=lev,face_M=float(key[0]),ray=key[1],projection_Ham_delta=float(qr[3]['Ham_input'])-float(qr[2]['Ham_input'])))
   for face in (3.,4.):
    for ray in RAYS:
     q=[r for r in rr if int(r['phase'])==3 and float(r['face_M'])==face and r['ray']==ray]
     if not q:continue
     peak=max(q,key=lambda r:abs(float(r['Theta_total'])));peakabs=abs(float(peak['Theta_total']));first=min((r for r in q if abs(float(r['Theta_total']))>=.01*peakabs),key=lambda r:float(r['time_M']))
     rk=[r for r in csv.DictReader(artifact(f't7a-rk-budget-{case}-L{lev}.csv').open()) if float(r['face_M'])==face and r['ray']==ray]
     net=max(rk,key=lambda r:float(r['total_increment_absmax']))
     source.append(dict(case=case,level=lev,face_M=face,ray=ray,first_1pct_source_time_M=float(first['time_M']),peak_source_time_M=float(peak['time_M']),peak_stage=int(peak['stage']),peak_stage_coarse_fraction=float(peak['coarse_fraction']),peak_start_fraction=float(peak['start_fraction']),Ham_at_peak=float(peak['Ham_input']),Theta_pre_KO_at_peak=float(peak['Theta_pre_KO']),Theta_KO_at_peak=float(peak['Theta_KO']),Theta_total_at_peak=float(peak['Theta_total']),half_alpha_Ham_at_peak=float(peak['half_alpha_Ham']),KO_over_physical_at_peak=float(peak['Theta_KO'])/float(peak['Theta_pre_KO']),net_RK_step_time_M=float(net['end_time_M']),net_physical_increment_at_peak=float(net['physical_increment_at_peak']),net_KO_increment_at_peak=float(net['KO_increment_at_peak']),net_total_increment_at_peak=float(net['total_increment_at_peak']),net_step_error=float(net['RK4_replay_error'])))
     for stage in range(4):
      z=[r for r in q if int(r['stage'])==stage]
      native.append(dict(case=case,level=lev,face_M=face,ray=ray,stage=stage,physical_Theta_RHS_absmax=max(abs(float(r['Theta_pre_KO'])) for r in z),KO_Theta_RHS_absmax=max(abs(float(r['Theta_KO'])) for r in z),Ham_input_absmax=max(abs(float(r['Ham_input'])) for r in z)))
 save('t7a-source-summary.csv',source);save('t7a-stage-summary.csv',native)
 pd=[]
 for case in ('clock','space'):
  for lev in (4,5,6):
   rr=[r for r in deltas if r['case']==case and r['level']==lev]
   z=list(csv.DictReader(artifact(f't7a-global-operation-deltas-{case}-L{lev}.csv').open()))
   pd.append(dict(case=case,level=lev,stage_projection_Ham_absmax=max(abs(r['projection_Ham_delta']) for r in rr),stage_projection_Theta_absmax=max(float(r['abs_max']) for r in z if r['operation']=='stage_projection' and r['field']=='Theta'),end_projection_Theta_absmax=max(float(r['abs_max']) for r in z if r['operation']=='end_projection' and r['field']=='Theta'),stage_projection_A_absmax=max(float(r['abs_max']) for r in z if r['operation']=='stage_projection' and r['field'] in ('A11','A12','A22','Aww'))))
 save('t7a-projection-summary.csv',pd)
 rs=[]
 for case in ('clock','space'):
  for lev in (4,5):
   raw=list(csv.DictReader(artifact(f't7a-restriction-{case}-L{lev}.csv').open()))
   for ray in RAYS:
    q=[r for r in raw if r['ray']==ray];peak=max(q,key=lambda r:abs(float(r['Theta_jump'])))
    rs.append(dict(case=case,coarse_level=lev,face_M=float(peak['face_M']),ray=ray,operations=len(q),peak_jump_time_M=float(peak['time_M']),Theta_jump_at_peak=float(peak['Theta_jump']),Theta_before_at_peak=float(peak['Theta_before']),Theta_after_at_peak=float(peak['Theta_after']),fine_Theta_support_absmax_at_peak=float(peak['fine_Theta_support_absmax']),replay_max_error=max(float(r['replay_max_error']) for r in q),condition='COVERED_COARSE_ONLY;CURRENT_FINE_PACKET_EXISTS_BEFORE_COPY'))
 save('t7a-restriction-summary.csv',rs)
 # Phase uncertainty is radial sampling sensitivity, not a transverse-error bound.
 tracks=list(csv.DictReader((HERE/'t7a-incident-tracks.csv').open()));phase=[]
 for case in ('T6_face3','clock','space'):
  for ray in RAYS:
   q={r['method']:r for r in tracks if r['case']==case and r['ray']==ray and abs(float(r['time_M'])-2.5)<1e-9}
   a,b=q['values'],q['values8']
   phase.append(dict(case=case,ray=ray,time_M=2.5,radius_P6_M=float(a['radius_M']),radius_P8_M=float(b['radius_M']),P6_P8_phase_difference_M=abs(float(a['radius_M'])-float(b['radius_M'])),sampling_step_M=1/384,Gamma_curvature_peak_P6=float(a['Gamma_curvature_peak']),Gamma_curvature_peak_P8=float(b['Gamma_curvature_peak']),width_P6_M=float(a['width_M']),width_P8_M=float(b['width_M']),condition='EXACT_DIAGONAL' if ray=='diagonal' else 'NATIVE_FIRST_ROW_H_OVER_2_OFFSET;TRANSVERSE_ERROR_NOT_BOUNDED'))
 save('t7a-phase-summary.csv',phase)
 # Measured widths and characteristic speed fits; no static-profile speed estimates.
 fits=[];widths=[]
 for case in ('clock','space'):
  for lev in (5,6):
   p=artifact(f't7a-native-widths-{case}-L{lev}.csv')
   if not p.exists():p=artifact(f't7a-native-widths-{case}.csv')
   rr=list(csv.DictReader(p.open()))
   for ray in RAYS:
    q=[r for r in rr if r['ray']==ray and r['lobe']=='negative' and (lev==5 or r['domain_end_baseline']=='False')]
    if lev==6:q=[r for r in q if .75<=float(r['time_M'])<=2.5 and float(r['radius_M'])>1.]
    else:
     q=[r for r in q if float(r['time_M'])>3.];branch=[]
     for row in q:
      if branch and float(row['radius_M'])<float(branch[-1]['radius_M'])-.01:break
      branch.append(row)
     q=branch
    if len(q)>=3:
     ts=np.array([float(r['time_M']) for r in q]);rads=np.array([float(r['radius_M']) for r in q]);(vel,intercept),cov=np.polyfit(ts,rads,1,cov=True)
     fits.append(dict(case=case,level=lev,ray=ray,samples=len(q),t_start_M=ts.min(),t_end_M=ts.max(),speed=vel,fit_standard_error=np.sqrt(cov[0,0]),fit_radius_RMS_M=np.sqrt(np.mean((rads-vel*ts-intercept)**2)),**{k:np.mean([float(r[k]) for r in q]) for k in ('lapse_speed','shift_transverse_speed','shift_longitudinal_speed','light_speed')}))
    for tt in ((2.5,) if lev==6 else (3.75,5.5)):
     use=[r for r in rr if r['ray']==ray and r['lobe']=='negative' and abs(float(r['time_M'])-tt)<1e-9]
     if use:widths.append(dict(level=lev,face_M=FACES[lev],**use[0]))
 save('t7a-speed-fits.csv',fits);save('t7a-width-summary.csv',widths)
 # Relative profile contrast on a common fixed ray; tensor support only at the face collar.
 contrasts=[]
 data={c:np.load(TMP/f'{c}-rays.npz') for c in ('T6_face3','clock','space')}
 for tt in (2.5,3.,4.,6.):
  for ri,ray in enumerate(RAYS):
   for method in ('values','tensor'):
    a={c:d[method][np.argmin(abs(d['time']-tt)),ri] for c,d in data.items()};a8={c:d[method+'8'][np.argmin(abs(d['time']-tt)),ri] for c,d in data.items()}
    use=np.all(np.isfinite(np.array(list(a.values()))),axis=(0,1))&np.all(np.isfinite(np.array(list(a8.values()))),axis=(0,1))&(R>.75)&(R*max(NV[ri])<2.98)
    if not use.any():continue
    for k,name in enumerate(V):
     scale=np.ptp(a['clock'][k,use]);dc=a['clock'][k,use]-a['T6_face3'][k,use];ds=a['space'][k,use]-a['clock'][k,use]
     contrasts.append(dict(time_M=tt,ray=ray,method=method,field=name,samples=int(use.sum()),clock_range=scale,clock_vs_T6_RMS=float(np.sqrt(np.mean(dc*dc))),space_vs_clock_RMS=float(np.sqrt(np.mean(ds*ds))),space_vs_clock_normalized_RMS=float(np.sqrt(np.mean(ds*ds))/scale) if scale else 0.,clock_vs_T6_max=float(np.max(abs(dc))),space_vs_clock_max=float(np.max(abs(ds))),P6_P8_max=max(float(np.max(abs(a[c][k,use]-a8[c][k,use]))) for c in a),ray_condition='EXACT_DIAGONAL' if ray=='diagonal' and method=='values' else 'QUALIFIED_FIXED_XY_TENSOR_COLLAR' if method=='tensor' else 'NATIVE_FIRST_ROW_OFFSET_H_OVER_2;TRANSVERSE_ERROR_NOT_IN_P6_P8'))
 save('t7a-profile-contrasts.csv',contrasts)
 print('SUMMARY',len(source),'sources',len(contrasts),'contrasts')

def ray_wakes():
 data={c:np.load(TMP/f'{c}-rays.npz') for c in ('T6_face3','clock','space')};rows=[]
 for tt,lo,hi in ((4.,3.125,3.5),(6.,3.125,5.)):
  for ri,ray in enumerate(RAYS):
   for field,k in (('Ham',28),('Theta',10),('GaussE',31)):
    for mode in ('values','values8'):
     rms={};use=(R>lo)&(R<hi)
     for case,f in data.items():
      q=f[mode][np.argmin(abs(f['time']-tt)),ri,k,use];rms[case]=np.sqrt(np.mean(q*q))
     rows.append(dict(time_M=tt,ray=ray,field=field,method=mode,r_min_M=lo,r_max_M=hi,T6_rms=rms['T6_face3'],clock_rms=rms['clock'],space_rms=rms['space'],clock_over_T6=rms['clock']/rms['T6_face3'],space_over_clock=rms['space']/rms['clock'],spatial_apparent_order=np.log2(rms['clock']/rms['space']),sampling_condition='FIXED_DIAGONAL_EXACT_RAY' if ray=='diagonal' else 'NATIVE_FIRST_ROW_OFFSET;P6_P8_RADIAL_ONLY'))
 save('t7a-ray-wake-ratios.csv',rows)


def figures():
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 import scienceplots
 plt.style.use(['science','no-latex']);plt.rcParams.update({'font.size':9})
 figdir=HERE/'figures';data={}
 for c in ('T6_face3','clock','space'):
  with np.load(TMP/f'{c}-rays.npz') as f:data[c]={k:f[k] for k in ('time','values','values8')}
 colors={'T6_face3':'#666666','clock':'#2166ac','space':'#b2182b'}
 labels={'T6_face3':'T6 face3','clock':'clock','space':'space'}
 fig,axs=plt.subplots(4,3,figsize=(10,10),layout='constrained')
 choices=[('lapse',13),('chi',0),('K',5),('Gamma_n',None),('shift_n',None),('driver B_n',None),('Theta',10),('phi',18),('Pi',19),('E_n',None),('Xi',27),('Gamma curvature',None)]
 for ax,(name,k) in zip(axs.flat,choices):
  for case,f in data.items():
   ti=np.argmin(abs(f['time']-2.5));q=f['values'][ti,2];q8=f['values8'][ti,2]
   def value(q):
    if k is not None:return q[k]
    if name=='Gamma curvature':return np.gradient(np.gradient((q[11]+q[12])/np.sqrt(2),R),R)
    i={'Gamma_n':11,'shift_n':14,'driver B_n':16,'E_n':24}[name];return (q[i]+q[i+1])/np.sqrt(2)
   z=value(q);z8=value(q8);shown=(R>=1.75)&(R<=2.75);ax.plot(R[shown],z[shown],color=colors[case],label=labels[case],lw=1.3);ax.fill_between(R[shown],z[shown],z8[shown],color=colors[case],alpha=.16)
  ax.set_xlim(1.75,2.75);ax.set_title(name);ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,3));ax.set_xlabel('r/M, fixed diagonal ray')
 axs.flat[0].legend(loc='best');fig.suptitle('Incident fields at t=2.5 M; shade = six/eight-node sampling difference')
 for ext in ('png','pdf'):fig.savefig(figdir/f't7a-incident-profiles.{ext}',dpi=220)
 plt.close(fig)
 hist=list(csv.DictReader((HERE/'t7a-common-packet-history.csv').open()));wake=list(csv.DictReader((HERE/'t7a-composite-metrics.csv').open()))
 fig,axs=plt.subplots(3,2,figsize=(10,9),layout='constrained')
 for col,face in enumerate((3.,4.)):
  fine=int(9-face)
  for row,field in enumerate(('Ham','Theta')):
   ax=axs[row,col]
   for case in data:
    q=sorted([r for r in hist if r['case']==case and float(r['face_M'])==face and int(r['level'])==fine and r['field']==field and r['ray']=='axis_plus'],key=lambda r:float(r['time_M']))
    ax.semilogy([float(r['time_M']) for r in q],np.maximum([float(r['peak_abs']) for r in q],1e-30),color=colors[case],label=labels[case])
   ax.set_ylim((1e-10,1e-3) if field=='Ham' else (1e-12,1e-5));ax.set_title(f'{face:g} M face: axial {field} packet peak');ax.set_xlim(0,6);ax.set_xlabel('t/M')
  ax=axs[2,col]
  for case in data:
   for band,style in (('wake_inner','-'),('wake_outer','--')):
    q=sorted([r for r in wake if r['case']==case and float(r['face_M'])==face and r['band']==band and r['field']=='Ham' and r['ray']=='all_face'],key=lambda r:float(r['time_M']))
    ax.semilogy([float(r['time_M']) for r in q],[float(r['rms']) for r in q],style,color=colors[case],label=labels[case]+(' inner' if style=='-' else ' outer'))
  ax.set_ylim(1e-10,1e-3);ax.set_xlim(0,6);ax.set_title(f'{face:g} M: volume Ham in fixed wake collars');ax.set_xlabel('t/M')
 axs[0,0].legend();axs[2,0].legend(ncol=2,fontsize=7)
 for ext in ('png','pdf'):fig.savefig(figdir/f't7a-packet-wake.{ext}',dpi=220)
 plt.close(fig)
 fig,axs=plt.subplots(2,4,figsize=(14,7),layout='constrained')
 for col,(face,ray,lev) in enumerate(((3.,'axis_plus',6),(3.,'diagonal',6),(4.,'axis_plus',5),(4.,'diagonal',5))):
  for case in ('clock','space'):
   rr=[r for r in csv.DictReader(artifact(f't7a-stage-native-{case}-L{lev}.csv').open()) if r['phase']=='3' and float(r['face_M'])==face and r['ray']==ray];groups=defaultdict(list)
   for r in rr:groups[round(float(r['time_M'])*24)/24].append(r)
   for field,style in (('Theta_pre_KO','-'),('Theta_KO',':'),('Theta_total','--')):
    ts=sorted(groups);vals=[float(max(groups[t],key=lambda r:abs(float(r['Theta_total'])))[field]) for t in ts]
    axs[0,col].plot(ts,vals,style,color=colors[case],label=case+' '+field.replace('Theta_',''))
   rk=[r for r in csv.DictReader(artifact(f't7a-rk-budget-{case}-L{lev}.csv').open()) if float(r['face_M'])==face and r['ray']==ray]
   groups=defaultdict(list)
   for r in rk:groups[round(float(r['end_time_M'])*24)/24].append(r)
   for field,style in (('physical_increment_at_peak','-'),('KO_increment_at_peak',':'),('total_increment_at_peak','--')):
    ts=sorted(groups);vals=[float(max(groups[t],key=lambda r:abs(float(r['total_increment_at_peak'])))[field]) for t in ts]
    axs[1,col].plot(ts,vals,style,color=colors[case])
  axs[0,col].set_title(f'{face:g} M {ray}: stage Theta RHS'+ ('\nendpoint censored' if col==3 else ''));axs[1,col].set_title('Actual RK4 step increments')
  for ax in axs[:,col]:ax.set_xlabel('t/M');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.axhline(0,color='k',lw=.4)
 axs[0,0].legend(fontsize=6,ncol=2);fig.suptitle('Native replay: physical RHS, KO and total; lower row uses actual RK4 weights')
 for ext in ('png','pdf'):fig.savefig(figdir/f't7a-operation-replay.{ext}',dpi=220)
 plt.close(fig)
 print('FIGURES PASS')

if __name__=='__main__':
 tasks={name:globals()[name] for name in ('snapshots','baseline','combined','profiles','operations','native_widths','rk_budget','restriction','summaries','common_packet_peaks','ray_wakes','figures')}
 if len(sys.argv)<2 or sys.argv[1] not in tasks:raise SystemExit('Choose one bounded analysis operation: '+', '.join(tasks))
 if sys.argv[1] in ('snapshots','operations','native_widths','rk_budget','restriction'):tasks[sys.argv[1]](sys.argv[2],int(sys.argv[3]))
 else:tasks[sys.argv[1]]()
