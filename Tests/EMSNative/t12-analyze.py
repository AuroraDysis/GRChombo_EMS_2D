#!/usr/bin/env python3
"""T12 Stage A only. Static data enter only this t=0 diagnostic; no evolution."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, math, os, subprocess, time, resource
from pathlib import Path
import numpy as np
import h5py
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];TMP=Path('/private/tmp/ems-t12');TMP.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('t11',HERE/'t11-analyze.py');t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t)
t.TMP=TMP;N=t.N;COL=t.COL;LEGS=('E-low','E-mid','E-high');PROFILE=t.PROFILE
def save(name,rows):
 with (HERE/('t12-'+name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def command(args,log=None,timeout=480):
 with (TMP/(log or 'command.log')).open('w') as f:subprocess.run(list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=timeout,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
def build():
 b=t.module('buildt12',HERE/'t7e-localize.py');b.TMP=TMP;b.build('T12Initial.cpp','initial.ex');b.build('T11RHS.cpp','rhs.ex')
 # Frozen headers precede the production include path; identical flags and source, only baseline opt-in code excluded.
 cmd=['/opt/homebrew/bin/g++-16','-O3','-std=c++17','-fopenmp','-DCH_SPACEDIM=2','-DCH_Darwin','-DCH_LANG_CC','-DNDEBUG','-DCH_USE_64','-DCH_USE_DOUBLE','-DT12_BASELINE','-I'+str(TMP/'baseline'),'-I'+str(ROOT/'Examples/EMS')]
 cmd+=['-I'+str(d) for d in (ROOT/'Source').rglob('*') if d.is_dir() and not d.name.startswith('.')]
 ch=Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib');cmd+=['-I'+str(ch/'src'/d) for d in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
 cmd+=[str(HERE/'T12Initial.cpp'),'-L'+str(ch),'-lboxtools2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC','-lbasetools2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC','-o',str(TMP/'baseline.ex')]
 command(cmd,'baseline-build.log',90)
def points(xy,mode):
 xy=np.asarray(xy,dtype='f8');xy.tofile(TMP/'points.bin')
 command([TMP/'initial.ex','--points',PROFILE,mode,TMP/'points.bin',TMP/'values.bin'],'points.log')
 a=np.fromfile(TMP/'values.bin').reshape(-1,32);assert len(a)==len(xy) and np.isfinite(a).all();return a
def census(leg):
 rows=[]
 with h5py.File(t.DATA/'exp-0020'/leg/'plt/EMS_Plot_000000.2d.hdf5') as f:
  assert float(f.attrs['time'])==0
  for l in range(13):
   g=f[f'level_{l}'];h=float(g.attrs['dx'])
   for k,b in enumerate(g['boxes'][:]):rows.append([l,k,*[int(b[a]) for a in t.KEYS],h])
 inp=TMP/(leg+'-boxes.txt');inp.write_text(''.join(' '.join(map(str,r))+'\n' for r in rows));start=time.monotonic()
 tables={}
 for label,binary,mode in (('baseline','baseline.ex','original'),('original','initial.ex','original'),('maximal','initial.ex','maximal')):
  out=HERE/f't12-census-{leg}-{label}.csv';command([TMP/binary,'--census',PROFILE,mode,inp,out],leg+'-'+label+'.log');tables[label]=list(csv.DictReader(out.open()))
 a,b,c=[tables[k] for k in ('baseline','original','maximal')]
 assert len(a)==len(b)==len(c)
 assert all(x['state_hash']==y['state_hash'] for x,y in zip(a,b)), 'default setter changed'
 assert all(x['nonlapse_hash']==y['nonlapse_hash'] for x,y in zip(b,c)), 'non-lapse setter changed'
 assert all(float(x['driver_cancel_max'])==0 and float(x['alpha_floor_margin'])>0 and float(x['chi_floor_margin'])>0 for x in c)
 print(leg,'census points',sum(int(x['points']) for x in a),'seconds',time.monotonic()-start,'default/non-lapse hashes PASS',flush=True)
def guards():
 rows=[]
 for mode in ('boost','boost_flag','binary','separation','ctt','helper_boost','helper_binary'):
  p=subprocess.run([str(TMP/'initial.ex'),'--guards',str(PROFILE),mode,'unused','unused'],capture_output=True,text=True,env={**os.environ,'OMP_NUM_THREADS':'1'})
  assert p.returncode!=0 and 'ems_use_maximal_initial_lapse' in p.stderr
  rows.append(dict(case=mode,exit=p.returncode,rejected=True,message=p.stderr.strip().splitlines()[-1]))
 save('rejections',rows)
def interp(ix,iy,a,origin):
 # Same anchor, accumulation order and six-point weights as PointAMRTransfer::apply.
 q=[(ix+.5)/2.-.5,(iy+.5)/2.-.5];start=[np.floor(v).astype(int)-2 for v in q];weights=[]
 for z,lo in zip(q,start):
  w=[]
  for i in range(6):
   u=np.ones(z.shape)
   for j in range(6):
    if i!=j:u*=(z-lo-j)/(i-j)
   w.append(u)
  weights.append(w)
 base=a[start[1]+2-origin[1],start[0]+2-origin[0]];correction=np.zeros(base.shape)
 for j in range(6):
  for i in range(6):correction+=weights[0][i]*weights[1][j]*(a[start[1]+j-origin[1],start[0]+i-origin[0]]-base)
 return base+correction
def native(leg):
 old=np.load('/private/tmp/ems-t11/'+leg+'-native.npz');m=old['meta'];coarse=None;grids={};transfer=[];records=TMP/'stencils.bin';allstate_equal=True
 with h5py.File(t.DATA/'exp-0020'/leg/'plt/EMS_Plot_000000.2d.hdf5') as f:
  for l in range(6,13):
   g=f[f'level_{l}'];h=float(g.attrs['dx']);bs=np.array([[b[k] for k in t.KEYS] for b in g['boxes'][:]],int);lo=bs[:,:2].min(0);hi=bs[:,2:].max(0)
   assert sum((b[2]-b[0]+1)*(b[3]-b[1]+1) for b in bs)==np.prod(hi-lo+1)
   xx,yy=np.meshgrid(np.arange(lo[0]-3,hi[0]+4),np.arange(lo[1]-3,hi[1]+4));xy=np.column_stack((((xx+.5)*h-336).ravel(),((yy+.5)*h).ravel()))
   a=points(xy,'maximal')[:,13].reshape(xx.shape)
   if l>6:
    cf=((xx<lo[0])|(xx>hi[0])|(yy>hi[1]))&(yy>=0);a[cf]=interp(xx[cf],yy[cf],coarse,origin)
   a[:3]=a[5:2:-1];assert a.min()>1e-12
   grids[l]=(a,lo-3);coarse=a;origin=lo-3
  for mode in ('original','maximal'):
   with records.open('wb') as stream:
    for l in range(7,13):
     g=f[f'level_{l}'];off=g['data:offsets=0'][:];bs=np.array([[b[k] for k in t.KEYS] for b in g['boxes'][:]],int)
     for bi,(x0,y0,x1,y1) in enumerate(bs):
      mm=m[(m[:,0]==l)&(m[:,6]==bi)]
      if not len(mm):continue
      ii=mm[:,2].astype(int)-x0;jj=mm[:,3].astype(int)-y0
      a=g['data:datatype=0'][int(off[bi]):int(off[bi+1])].reshape(34,y1-y0+7,x1-x0+7)[:N]
      rec=np.empty((len(mm),3+N*49));rec[:,:3]=mm[:,[1,5,4]]
      for dy in range(-3,4):
       for dx in range(-3,4):rec[:,3+(dy+3)*7+dx+3::49]=a[:,jj+3+dy,ii+3+dx].T
      if mode=='maximal':
       alpha,org=grids[l]
       for dy in range(-3,4):
        for dx in range(-3,4):rec[:,3+13*49+(dy+3)*7+dx+3]=alpha[mm[:,3].astype(int)+dy-org[1],mm[:,2].astype(int)+dx-org[0]]
      rec.tofile(stream)
   command([TMP/'rhs.ex','--t0-native',records,TMP/'rhs.bin'],'native-'+leg+'-'+mode+'.log')
   q=np.fromfile(TMP/'rhs.bin').reshape(-1,COL);assert len(q)==len(m) and np.isfinite(q).all() and q[:,-1].sum()==0
   np.savez_compressed(TMP/f'{leg}-{mode}-native.npz',meta=m,q=q)
   if mode=='original':assert np.array_equal(q.view('u8'),old['q'].view('u8')),'T11 default RHS changed'
   else:
    non=[i for i in range(N) if i!=13];assert np.array_equal(q[:,non].copy().view('u8'),old['q'][:,non].copy().view('u8'))
    assert np.array_equal(q[:,2*N+6:2*N+10].copy().view('u8'),old['q'][:,2*N+6:2*N+10].copy().view('u8'))
    assert abs(q[:,138:140]).max()==0 and ((q[:,0]>1e-12)&(q[:,13]>1e-12)).all()
 records.unlink();(TMP/'rhs.bin').unlink()
 print(leg,'native cells',len(m),'default all columns bit-identical; non-lapse, KO(A), driver PASS',flush=True)
def continuum():
 for leg in LEGS:
  m=np.load(TMP/f'{leg}-original-native.npz')['meta']
  for mode in ('original','maximal'):chunked_continuum(m[:,[4,5]],leg+'-'+mode,mode)
 r=np.unique(np.r_[np.geomspace(1e-6,1e3,5001),np.linspace(1e-4,.5,3000),[.0001,.0002,.001,.005,.01,.02,.05,.1,.2,.5]])
 for mode in ('original','maximal'):chunked_continuum(np.column_stack((r/np.sqrt(2),r/np.sqrt(2))),'radial-'+mode,mode)
 a=points(np.column_stack((r,np.zeros(len(r)))),'maximal');z=np.load(TMP/'radial-maximal-continuum.npz');j=z['radial_jets']
 assert (a[:,30]>0).all() and (a[:,31]>0).all() and (a[:,13]>1e-12).all() and (a[:,0]>1e-12).all()
 rows=[dict(r_M=R,Y=a[i,30],G=a[i,31],alpha_K=a[i,13],chi=a[i,0],alpha_over_Rnu=a[i,13]/R**1.3372155112113842,alpha_floor_margin=a[i,13]-1e-12,chi_floor_margin=a[i,0]-1e-12,reader_analytic_alpha_difference=a[i,13]-j[i,1,0]) for i,R in enumerate(r)]
 save('reader-line',rows)
def chunked_continuum(xy,name,mode):
 # ponytail: 1024-point chunks bound every Jet operation; disk-backed arrays avoid accumulating chunks.
 xy=np.asarray(xy);shapes={'xy':xy.shape,'q':(len(xy),COL),'Killing_lapse':(len(xy),),'radial_jets':(len(xy),9,3)}
 arrays={key:np.lib.format.open_memmap(TMP/(name+'-'+key+'.npy'),mode='w+',dtype='f8',shape=shape) for key,shape in shapes.items()}
 memory=[]
 for first in range(0,len(xy),1024):
  last=min(len(xy),first+1024);chunk=name+'-chunk';t.continuum(xy[first:last],chunk,mode)
  with np.load(TMP/(chunk+'-continuum.npz')) as z:
   for key,a in arrays.items():a[first:last]=z[key]
  memory.append(dict(comparison=name,first=first,points=last-first,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
  assert memory[-1]['peak_RSS_bytes']<2*1024**3,'2 GiB diagnostic memory stop'
 for a in arrays.values():a.flush()
 np.savez_compressed(TMP/(name+'-continuum.npz'),**arrays)
 with (TMP/'continuum-memory.csv').open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=memory[0]);
  if f.tell()==0:w.writeheader()
  w.writerows(memory)
 return np.array(arrays['q'])
def measure():
 inner=[];profiles=[];control=[];orders=[];fixed=[];inputs=[]
 for leg in LEGS:
  oo=np.load(TMP/f'{leg}-original-native.npz');m=oo['meta'];o=oo['q'];k=np.load(TMP/f'{leg}-maximal-native.npz')['q']
  control.append(dict(leg=leg,cells=len(m),default_rhs_bit_mismatches=int(np.sum(o.copy().view('u8')!=np.load('/private/tmp/ems-t11/'+leg+'-native.npz')['q'].copy().view('u8'))),nonlapse_bits_changed=int(np.sum(o[:,[i for i in range(N) if i!=13]].copy().view('u8')!=k[:,[i for i in range(N) if i!=13]].copy().view('u8'))),driver_cancel_original=float(abs(o[:,138:140]).max()),driver_cancel_maximal=float(abs(k[:,138:140]).max()),maximal_term_defect_eps=float(k[:,-2].max()),min_lapse=float(k[:,13].min()),min_chi=float(k[:,0].min()),finite=True))
  for ray in ('axis','diagonal'):
   for mode,a in (('original',o),('maximal',k)):
    cc=np.load(TMP/f'{leg}-{mode}-continuum.npz')['q']
    for l in range(7,13):
     use=(m[:,0]==l)&(m[:,4]>0)&(m[:,4]<=.5)&((m[:,3]==0) if ray=='axis' else ((abs(m[:,4]-m[:,5])<m[:,1]/100)&(np.hypot(m[:,4],m[:,5])<=.5)))
     ids=np.where(use)[0];ids=ids[np.argsort(m[ids,4])];rr=m[ids,4] if ray=='axis' else np.hypot(m[ids,4],m[ids,5]);g=t.projected(a[ids],ray,N,11);gc=t.projected(cc[ids],ray,N,11)
     for x,i,R in zip(g,ids,rr):
      row=dict(leg=leg,mode=mode,ray=ray,level=l,h_M=m[i,1],r_M=R,x_M=m[i,4],y_M=m[i,5],physical_R_M=np.hypot(m[i,4],m[i,5]),covered=bool(l<12 and max(m[i,4],m[i,5])<112/2**(l+1)),lapse=a[i,13],Gamma_rhs=x,continuum_Gamma=float(t.projected(cc[i:i+1],ray,N,11)[0]))
      row.update({name:float(t.projected(a[i:i+1],ray,4*N+2*ig,0)[0]) for ig,name in enumerate(t.GROUPS)})
      row.update(alpha_advection=a[i,3*N+13],beta_advection=float(t.projected(a[i:i+1],ray,3*N,14)[0]),B_decay=float(t.projected(a[i:i+1],ray,141,0)[0]),A12_KO=a[i,2*N+7],lapse_KO=a[i,2*N+13],h11_RHS=a[i,N+1],h12_RHS=a[i,N+2]);profiles.append(row)
     if l==12:
      row=dict(profiles[-len(ids)]);row.update(continuum_lapse=cc[ids[0],13],Gamma_error=g[0]-gc[0],h_times_A12_KO=m[ids[0],1]*a[ids[0],2*N+7]);inner.append(row)
 # Reuse T11's fixed physical supports and parity; native file alias is temporary under T12 only.
 r=np.unique(np.r_[np.geomspace(.0002,.025,180),np.linspace(.026,.49,350),[.0002,.001,.005,.01,.02,.05,.1,.2,.4]])
 for ray in ('axis','diagonal'):
  keep=np.ones(len(r),bool);fac=1 if ray=='axis' else np.sqrt(2)
  for l in range(7,13):keep &= abs(r-112/2**l*fac)>5*(.875/2**l)*fac
  for mode in ('original','maximal'):
   samples=[];samples8=[]
   for leg in LEGS:
    alias=TMP/f'{leg}-native.npz';alias.unlink(missing_ok=True);alias.symlink_to(TMP/f'{leg}-{mode}-native.npz')
    samples.append(t.sample(leg,ray,r,6)[0]);samples8.append(t.sample(leg,ray,r,8)[0])
   xy=np.column_stack((r,r*0+1e-14)) if ray=='axis' else np.column_stack((r/np.sqrt(2),r/np.sqrt(2)))
   c=chunked_continuum(xy,'fixed-'+ray+'-'+mode,mode)
   for lo,hi in ((.0002,.001),(.001,.005),(.005,.02),(.02,.05),(.05,.1),(.1,.49)):
    use=keep&(r>=lo)&(r<=hi);v=[t.projected(a,ray,N,11) for a in samples];v8=[t.projected(a,ray,N,11) for a in samples8]
    d01=t.rms((v[0]-v[1])[use]);d12=t.rms((v[1]-v[2])[use]);spread=max(t.rms((a-b)[use]) for a,b in zip(v,v8))
    orders.append(dict(mode=mode,ray=ray,r_min_M=lo,r_max_M=hi,points=int(use.sum()),low_mid_RMS=d01,mid_high_RMS=d12,order=math.log(d01/d12)/math.log(1.5),P6_P8_RMS=spread,interpolation_fraction=spread/d12))
   for R in (.0002,.001,.005,.01,.02,.05,.1,.2,.4):
    i=int(abs(r-R).argmin())
    for leg,a,a8 in zip(LEGS,samples,samples8):fixed.append(dict(leg=leg,mode=mode,ray=ray,r_M=r[i],requested_r_M=R,lapse=a[i,13],Gamma_rhs=float(t.projected(a[i:i+1],ray,N,11)[0]),Gamma_continuum=float(t.projected(c[i:i+1],ray,N,11)[0]),Gamma_P6_P8_spread=float(abs(t.projected(a[i:i+1]-a8[i:i+1],ray,N,11)[0]))))
 save('controls',control);save('puncture-terms',inner);save('native-profiles',profiles);save('rhs-orders',orders);save('fixed-radii',fixed)
def checks():
 assert all(int(x['default_rhs_bit_mismatches'])==int(x['nonlapse_bits_changed'])==0 and float(x['driver_cancel_maximal'])==float(x['driver_cancel_original'])==0 for x in csv.DictReader((HERE/'t12-controls.csv').open()))
 for leg in LEGS:
  a=list(csv.DictReader((HERE/f't12-census-{leg}-baseline.csv').open()));b=list(csv.DictReader((HERE/f't12-census-{leg}-original.csv').open()));c=list(csv.DictReader((HERE/f't12-census-{leg}-maximal.csv').open()))
  assert all(x['state_hash']==y['state_hash'] for x,y in zip(a,b)) and all(x['nonlapse_hash']==y['nonlapse_hash'] for x,y in zip(b,c))
 print('T12 default bit identity, non-lapse identity, reader/floors, native split and driver checks PASS')
def identity_points():
 rows=[]
 for leg in LEGS:
  for mode in ('original','maximal'):
   with np.load(TMP/f'{leg}-{mode}-continuum.npz') as z:
    xy=z['xy'];r=np.hypot(*xy.T);j=z['radial_jets'];q=z['q']
   x,ak,b,k,ch,p,ph=[j[:,i] for i in range(7)];a=x if mode=='original' else ak
   mb=b[:,1]-b[:,0]/r-3*ak[:,0]*k[:,0];mbp=b[:,2]-b[:,1]/r+b[:,0]/r**2-3*(ak[:,1]*k[:,0]+ak[:,0]*k[:,1])
   cm=2*k[:,1]+6*k[:,0]/r-3*k[:,0]*ch[:,1]/ch[:,0]-16*np.pi*p[:,0]*ph[:,1]
   mismatch=4*(k[:,1]*(ak[:,0]-a[:,0])+k[:,0]*(ak[:,1]-a[:,1]))+12*k[:,0]*(ak[:,0]-a[:,0])/r
   residual=4/3*(mbp+3*mb/r)+2*a[:,0]*cm;source=np.sum(q[:,39:41]*xy/r[:,None],axis=1)
   scale=4/3*(abs(b[:,2])+2*abs(b[:,1])/r+2*abs(b[:,0])/r**2)+4*abs(k[:,0]*a[:,1])+6*abs(a[:,0]*k[:,0]*ch[:,1]/ch[:,0])+32*np.pi*abs(a[:,0]*p[:,0]*ph[:,1])
   budget=128*np.finfo(float).eps*np.maximum(1,scale);defect=source-mismatch-residual
   assert (abs(defect)<=budget).all()
   rows.append(dict(leg=leg,mode=mode,points=len(r),r_min_M=r.min(),M_beta_max_abs=abs(mb).max(),M_beta_prime_max_abs=abs(mbp).max(),C_M_max_abs=abs(cm).max(),source_max_abs=abs(source).max(),residual_source_max_abs=abs(residual).max(),identity3_defect_max_abs=abs(defect).max(),defect_budget_ratio_max=float(np.max(abs(defect)/budget))))
 save('native-identity',rows);print(rows)
if __name__=='__main__':
 start=time.monotonic();action=sys.argv[1]
 if action=='build':build()
 elif action=='census':census(sys.argv[2])
 elif action=='native':native(sys.argv[2])
 elif action=='guards':guards()
 elif action=='continuum':continuum()
 elif action=='measure':measure()
 elif action=='check':checks()
 elif action=='identity-points':identity_points()
 else:raise SystemExit('unknown action')
 print(action,'wall_s',time.monotonic()-start,'RSS_GB',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9,flush=True)
