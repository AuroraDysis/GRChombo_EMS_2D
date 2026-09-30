#!/usr/bin/env python3
"""Stream sampled native stencils, current ghosts only; radial data strictly t=0."""
import csv,importlib.util,json,os,subprocess,time,hashlib,resource,math
from pathlib import Path
import numpy as np,h5py
HERE=Path(__file__).resolve().parent;TMP=Path('/private/tmp/ems-t7e');INPUT=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0019')
spec=importlib.util.spec_from_file_location('t4b',HERE/'t4b-diagnose.py');t4b=importlib.util.module_from_spec(spec);spec.loader.exec_module(t4b)
MASKS=t4b.MASKS
ODD=[2,7,12,15,17,22,23,25,26]
PROFILE=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet')
def save(name,rows):
 directory=TMP/'tables' if any(name.startswith('t7e-'+prefix) for prefix in ('local-','profile-','coordinate-control-','gamma-control-','performance-')) else HERE
 directory.mkdir(exist_ok=True)
 with (directory/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def build(source="T7EInitial.cpp",binary="initial.ex"):
 repo=HERE.parents[1];ch=Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
 cmd=['/opt/homebrew/bin/g++-16','-O3','-std=c++17','-fopenmp','-DCH_SPACEDIM=2','-DCH_Darwin','-DCH_LANG_CC','-DNDEBUG','-DCH_USE_64','-DCH_USE_DOUBLE','-I'+str(repo/'Examples/EMS')]
 cmd+=['-I'+str(d) for d in (repo/'Source').rglob('*') if d.is_dir() and not d.name.startswith('.')]
 cmd+=['-I'+str(ch/'src'/d) for d in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
 cmd+=[str(HERE/source),'-L'+str(ch),'-lboxtools2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC','-lbasetools2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC','-o',str(TMP/binary)]
 with (TMP/(binary+'-build.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
def interp6(ix,iy,a,origin):
 q=[ix/2.-.25,iy/2.-.25];start=[np.floor(v).astype(int)-2 for v in q];weights=[]
 for v,lo in zip(q,start):
  ww=[]
  for k in range(6):
   w=np.ones(v.shape)
   for j in range(6):
    if j!=k:w*=(v-lo-j)/(k-j)
   ww.append(w)
  weights.append(ww)
 base=a[:,start[1]+2-origin[1],start[0]+2-origin[0]];out=base.copy()
 for j in range(6):
  for i in range(6):out+=(a[:,start[1]+j-origin[1],start[0]+i-origin[0]]-base)*(weights[0][i]*weights[1][j])
 return out
def sample_checkpoint(scale,cpname):
 begin=time.monotonic();cp=INPUT/('E-'+scale)/'chk'/cpname;parts=[];gamma_parts=[];payload=TMP/'stencils.bin';out=TMP/'budgets.bin';initial=TMP/'initial.bin';xyfile=TMP/'xy.bin'
 # A bounded deterministic census: stride 4 in each direction; all interface/axis strips retained.
 with h5py.File(cp) as f,payload.open('wb') as stream:
  t=float(f.attrs['time']);coarse=origin=None
  for lev in range(int(f.attrs['num_levels'])):
   g=f[f'level_{lev}'];assert abs(float(g.attrs['time'])-t)<=64*np.finfo(float).eps*max(1,abs(t));b=t4b.boxes(g);lo=b[:,:2].min(axis=0);hi=b[:,2:].max(axis=0);h=float(g.attrs['dx']);assert np.sum((b[:,2]-b[:,0]+1)*(b[:,3]-b[:,1]+1))==(hi[0]-lo[0]+1)*(hi[1]-lo[1]+1);off=g['data:offsets=0'][:]
   a=np.empty((28,hi[1]-lo[1]+7,hi[0]-lo[0]+7));blocks=[]
   for k,(x0,y0,x1,y1) in enumerate(b):
    v=g['data:datatype=0'][int(off[k]):int(off[k+1])].reshape(28,y1-y0+7,x1-x0+7);blocks.append(v)
    a[:,y0-lo[1]:y1-lo[1]+7,x0-lo[0]:x1-lo[0]+7]=v
   for v,(x0,y0,x1,y1) in zip(blocks,b):a[:,y0-lo[1]+3:y1-lo[1]+4,x0-lo[0]+3:x1-lo[0]+4]=v[:,3:-3,3:-3]
   X,Y=np.meshgrid(np.arange(lo[0]-3,hi[0]+4),np.arange(lo[1]-3,hi[1]+4))
   if lev:
    cf=((X<lo[0])|(X>hi[0])|(Y>hi[1]))&(Y>=0)
    a[:,cf]=interp6(X[cf],Y[cf],coarse,origin)
   if lo[1]==0:
    a[:,:3]=a[:,5:2:-1];a[ODD,:3]*=-1
   fine=t4b.boxes(f[f'level_{lev+1}'])//2 if lev+1<int(f.attrs['num_levels']) else np.empty((0,4),int)
   xx,yy=X[3:-3,3:-3],Y[3:-3,3:-3];r=np.hypot((xx+.5)*h-224,(yy+.5)*h);valid=np.ones(xx.shape,bool);near=np.zeros(xx.shape,bool)
   for x0,y0,x1,y1 in fine:
    cov=(xx>=x0)&(xx<=x1)&(yy>=y0)&(yy<=y1);valid&=~cov
    near|=(xx>=x0-2)&(xx<=x1+2)&(yy>=y0-2)&(yy<=y1+2)&~cov
   edge=((xx-lo[0]<2)|(hi[0]-xx<2)|(hi[1]-yy<2))&(lev>0)
   axis=yy<2;interior=~(edge|near|axis)
   region=np.zeros(xx.shape,bool)
   names=['inside_inner_ring','cavity','between_rings','cavity_core','far','far_core']
   for name in names:region|=(r>=MASKS[name][0])&(r<=MASKS[name][1])
   take=valid&region&(((xx%4==0)&(yy%4==0))|edge|near|axis|((r>=4)&(r<=8)))
   jj,ii=np.where(take)
   if len(ii):
    records=np.empty((len(ii),703));records[:,0]=h;records[:,1]=(yy[take]+.5)*h;records[:,2]=(xx[take]+.5)*h-224
    for z in range(-2,3):
     for x in range(-2,3):records[:,3+(z+2)*5+x+2::25]=a[:,jj+3+z,ii+3+x].T
    records.tofile(stream)
    p=dict(level=np.full(len(ii),lev),h=np.full(len(ii),h),x=records[:,2],y=records[:,1],r=r[take],edge=edge[take],near=near[take],axis=axis[take],interior=interior[take],chi=a[0,jj+3,ii+3])
    p['fields']=a[:,jj+3,ii+3].T.copy();p['weight']=np.where(edge[take]|near[take]|axis[take]|((r[take]>=4)&(r[take]<=8)),1.,16.)
    for name in names:p[name]=(r[take]>=MASKS[name][0])&(r[take]<=MASKS[name][1])
    parts.append(p)
   if t>0 and lev in (3,4):
    take9=valid&(r>=4)&(r<=8)&(xx%16==0)&(yy%16==0)&(xx>=lo[0]+4)&(xx<=hi[0]-4)&(yy>=lo[1]+4)&(yy<=hi[1]-4)
    jj9,ii9=np.where(take9)
    if len(ii9):
     rec=np.empty((len(ii9),2271));rec[:,0]=h;rec[:,1]=(yy[take9]+.5)*h;rec[:,2]=(xx[take9]+.5)*h-224
     for z in range(-4,5):
      for x in range(-4,5):rec[:,3+(z+4)*9+x+4::81]=a[:,jj9+3+z,ii9+3+x].T
     gamma_parts.append((lev,rec))
   coarse=a;origin=lo-3
 data={k:np.concatenate([p[k] for p in parts]) for k in parts[0]}
 print(scale,cpname,'sampled',len(data['r']),'stencil_MB',payload.stat().st_size/1e6,flush=True)
 start=time.monotonic();subprocess.run(['/private/tmp/ems-t4/t4b/roundoff.ex',str(payload),str(out)],check=True,timeout=115)
 budget=np.fromfile(out,dtype='f8').reshape(-1,13);assert len(budget)==len(data['r']) and np.isfinite(budget).all()
 subprocess.run([str(TMP/'native.ex'),str(payload),str(TMP/'native.bin')],check=True,timeout=115)
 native=np.fromfile(TMP/'native.bin',dtype='f8').reshape(-1,3);(TMP/'native.bin').unlink()
 if t==0:
  sampled=np.memmap(payload,dtype='f8',mode='r',shape=(len(native),703))
  ix=np.rint((data['x']+224)/data['h']-.5).astype(int);iy=np.rint(data['y']/data['h']-.5).astype(int)
  pick=data['interior']&(ix%32==0)&(iy%32==0);indexes=np.where(pick)[0];probe=sampled[pick].copy();del sampled
  positions=[]
  for z in range(-2,3):
   for x in range(-2,3):positions.append(np.column_stack((ix[pick]+x,iy[pick]+z,data['h'][pick])))
  positions=np.stack(positions,axis=1);positions.tofile(TMP/'fused-coords.bin')
  with (TMP/'fused-coords.bin').open('rb') as inp,(TMP/'fused-values.bin').open('wb') as dest,(TMP/'fused.log').open('w') as err:
   subprocess.run([str(TMP/'initial.ex'),str(PROFILE),'--t0-fused'],stdin=inp,stdout=dest,stderr=err,check=True,timeout=115)
  values=np.fromfile(TMP/'fused-values.bin',dtype='f8').reshape(len(probe),25,35)
  probe[:,3:28]=values[:,:,0];probe.tofile(TMP/'fused-stencils.bin')
  subprocess.run([str(TMP/'native.ex'),str(TMP/'fused-stencils.bin'),str(TMP/'fused-native.bin')],check=True,timeout=115)
  fused=np.fromfile(TMP/'fused-native.bin',dtype='f8').reshape(-1,3);control=[]
  for name in names:
   ch=data[name][pick]
   if not ch.any():continue
   w=data['y'][pick][ch];rms=lambda a:np.sqrt(np.sum(w*a*a)/np.sum(w))
   control.append(dict(scale=scale,mask=name,time_M=0.,bulk_cells=int(ch.sum()),native_Ham=rms(native[pick,0][ch]),chi_fused_Ham=rms(fused[ch,0]),native_Mom=rms(native[pick,1][ch]),chi_fused_Mom=rms(fused[ch,1]),max_chi_change=np.max(abs(values[ch,:,0]-np.memmap(payload,dtype='f8',mode='r',shape=(len(native),703))[indexes[ch],3:28]))))
  save(f't7e-coordinate-control-{scale}.csv',control)
  for filename in ('fused-coords.bin','fused-values.bin','fused-stencils.bin','fused-native.bin'):(TMP/filename).unlink()
 payload.unlink();out.unlink()
 cont=None
 if t==0:
  # No static-file evaluation is reachable for an evolved checkpoint.
  np.column_stack((data['x'],data['y'])).tofile(xyfile)
  with xyfile.open('rb') as inp,initial.open('wb') as dest,(TMP/'initial.log').open('w') as err:
   subprocess.run([str(TMP/'initial.ex'),str(PROFILE),'--t0'],stdin=inp,stdout=dest,stderr=err,check=True,timeout=115)
  cont=np.fromfile(initial,dtype='f8').reshape(-1,35);assert len(cont)==len(budget) and np.isfinite(cont).all();xyfile.unlink();initial.unlink()
 if gamma_parts:
  np.concatenate([v for lev,v in gamma_parts]).tofile(TMP/'gamma.bin')
  subprocess.run([str(TMP/'native.ex'),str(TMP/'gamma.bin'),str(TMP/'gamma-out.bin'),'--gamma'],check=True,timeout=115)
  value=np.fromfile(TMP/'gamma-out.bin',dtype='f8').reshape(-1,4);offset=0;rr=[]
  for lev,v in gamma_parts:
   q=value[offset:offset+len(v)];offset+=len(v);w=v[:,1];rr.append(dict(scale=scale,time_M=t,level=lev,cells=len(v),stored_Gamma_Ham=np.sqrt(np.sum(w*q[:,0]**2)/sum(w)),metric_Gamma_Ham=np.sqrt(np.sum(w*q[:,3]**2)/sum(w)),Ham_change=np.sqrt(np.sum(w*(q[:,0]-q[:,3])**2)/sum(w))))
  save(f't7e-gamma-control-{scale}-{cpname[4:10]}.csv',rr)
  (TMP/'gamma.bin').unlink();(TMP/'gamma-out.bin').unlink()
 rows=[]
 def rms(v,w):return np.sqrt(np.sum(w*v*v)/np.sum(w))
 for name in names:
  for subset,flag in [('all',np.ones(len(budget),bool)),('interface',data['edge']),('near_finer',data['near']),('axis',data['axis']),('bulk',data['interior'])]+[(f'level_{lev}',data['level']==lev) for lev in np.unique(data['level'])]:
   choose=data[name]&flag
   if not choose.any():continue
   w=data['y'][choose]*data['weight'][choose];a=budget[choose]
   wv=w*data['h'][choose]**2
   row=dict(scale=scale,checkpoint=cpname,time_M=t,mask=name,subset=subset,sampled_cells=int(choose.sum()),estimated_cells=int(data['weight'][choose].sum()),levels=';'.join(map(str,np.unique(data['level'][choose]))),min_h_M=float(data['h'][choose].min()),max_h_M=float(data['h'][choose].max()),volume_Ham_rms=rms(native[choose,0],wv),volume_Mom_rms=rms(native[choose,1],wv),volume_GaussE_rms=rms(native[choose,2],wv),Ham_rms=rms(native[choose,0],w),Mom_rms=rms(native[choose,1],w),GaussE_rms=rms(native[choose,2],w),budget_replay_Ham_rms=rms(a[:,0],w),native_budget_difference=rms(native[choose,0]-a[:,0],w),Ham_roundoff=rms(a[:,3],w),Ham_initial_evaluation_budget=rms(a[:,10],w) if t==0 else np.nan,eps_chi_h2=rms(np.finfo(float).eps*data['chi'][choose]/data['h'][choose]**2,w),metric_Ricci=rms(a[:,6],w),Ham_without_metric_Ricci=rms(a[:,7],w),metric_span_max=a[:,8].max(),gamma_identity_defect=rms(a[:,9],w))
   if cont is not None:
    c=cont[choose];field=data['fields'][choose]
    row.update(continuum_Ham_rms=rms(c[:,28],w),continuum_terms_roundoff=rms(np.finfo(float).eps*np.sum(abs(c[:,29:32]),axis=1),w),setter_chi_relative_rms=rms((c[:,0]-field[:,0])/field[:,0],w),setter_phi_absolute_rms=rms(c[:,18]-field[:,18],w),coordinate_Ham_sensitivity=rms((32/3)*np.finfo(float).eps*224*abs(c[:,32]*data['x'][choose]/data['r'][choose])/data['h'][choose]**2,w),setter_h_absolute_rms=rms(np.max(abs(c[:,1:5]-field[:,1:5]),axis=1),w))
   rows.append(row)
 # Radial bins localize which level/face carries the far floor and its change in time.
 profiles=[]
 for lev in np.unique(data['level']):
  for lo in np.arange(4,8,.125):
   choose=(data['level']==lev)&(data['r']>=lo)&(data['r']<lo+.125)
   if not choose.any():continue
   w=data['y'][choose]*data['weight'][choose];a=budget[choose]
   wv=w*data['h'][choose]**2
   profiles.append(dict(scale=scale,checkpoint=cpname,time_M=t,level=int(lev),r_mid_M=lo+.0625,sampled_cells=int(choose.sum()),volume_Ham_rms=rms(native[choose,0],wv),volume_Mom_rms=rms(native[choose,1],wv),volume_GaussE_rms=rms(native[choose,2],wv),Ham_rms=rms(native[choose,0],w),Mom_rms=rms(native[choose,1],w),metric_Ricci_rms=rms(a[:,6],w),bulk_Ham_rms=rms(native[choose&data['interior'],0],data['y'][choose&data['interior']]*data['weight'][choose&data['interior']]) if (choose&data['interior']).any() else np.nan))
 save(f't7e-local-{scale}-{cpname[4:10]}.csv',rows);save(f't7e-profile-{scale}-{cpname[4:10]}.csv',profiles)
 save(f't7e-performance-{scale}-{cpname[4:10]}.csv',[dict(scale=scale,checkpoint=cpname,wall_seconds=time.monotonic()-begin,parent_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,child_peak_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,threads=1)])
 print('native replay_seconds',time.monotonic()-start,flush=True)
def layout_audit():
 rows=[];layout=[]
 for scale in ('low','mid','high'):
  with h5py.File(INPUT/('E-'+scale)/'chk/EMS_000000.2d.hdf5') as f:
   acc={k:[0.,0.,0,0.] for k in ('inside_inner_ring','cavity','between_rings','cavity_core')}
   for lev in range(13):
    g=f[f'level_{lev}'];h=float(g.attrs['dx']);b=t4b.boxes(g);lo=b[:,:2].min(0);hi=b[:,2:].max(0)
    assert np.sum((b[:,2]-b[:,0]+1)*(b[:,3]-b[:,1]+1))==(hi[0]-lo[0]+1)*(hi[1]-lo[1]+1)
    layout.append(dict(scale=scale,level=lev,h_M=h,x_min_M=float(lo[0]*h-224),x_max_M=float((hi[0]+1)*h-224),y_max_M=float((hi[1]+1)*h),boxes=len(b),x_seams_M=';'.join(map(str,sorted(set(float(z[0]*h-224) for z in b if z[0]!=lo[0])))),y_seams_M=';'.join(map(str,sorted(set(float(z[1]*h) for z in b if z[1]!=lo[1]))))))
    fine=t4b.boxes(f[f'level_{lev+1}'])//2 if lev<12 else []
    for x0,y0,x1,y1 in b:
     xx,yy=np.meshgrid(np.arange(x0,x1+1),np.arange(y0,y1+1));x=(xx+.5)*h-224;y=(yy+.5)*h;r=np.hypot(x,y);valid=np.ones(xx.shape,bool)
     for a,c,d,e in fine:valid&=~((xx>=a)&(xx<=d)&(yy>=c)&(yy<=e))
     if not any(np.any(valid&(r>=MASKS[k][0])&(r<=MASKS[k][1])) for k in acc):continue
     exact=np.fromiter((math.fma(float(i)+.5,h,-224.) for i in range(x0,x1+1)),float);error=x-exact[None,:]
     for mask,v in acc.items():
      a,c=MASKS[mask];ch=valid&(r>=a)&(r<=c);w=y[ch];e=error[ch];v[0]+=float(np.sum(w*e*e));v[1]+=float(sum(w));v[2]+=int(ch.sum());v[3]=max(v[3],float(np.max(abs(e),initial=0)))
   for mask,v in acc.items():rows.append(dict(scale=scale,mask=mask,cells=v[2],coordinate_error_rms_M=math.sqrt(v[0]/v[1]),coordinate_error_max_M=v[3],reference='CORRECTLY_ROUNDED_FMA_GLOBAL_COORDINATE'))
 save('t7e-coordinate-errors.csv',rows);save('t7e-layout.csv',layout)

if __name__=='__main__':
 import sys
 if sys.argv[1]=='layout':layout_audit()
 elif sys.argv[1]=='build':build();build('T7ENative.cpp','native.ex')
 else:sample_checkpoint(sys.argv[1],sys.argv[2])
