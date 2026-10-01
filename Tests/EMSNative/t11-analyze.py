#!/usr/bin/env python3
"""T11 offline RHS audit. Every static-file operation is restricted to t=0."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, json, math, resource, subprocess, time
from pathlib import Path
import h5py
import numpy as np

HERE=Path(__file__).resolve().parent; TMP=Path('/private/tmp/ems-t11');TMP.mkdir(exist_ok=True)
DATA=Path('/Users/auroradysis/Workspace/EMS/.data'); PROFILE=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet')
N=28; G=13; COL=4*N+2*G+7
GROUPS=('advection','shift_laplacian_2D','shift_graddiv_2D','lapse_A','chi_A','K_gradient','Theta_gradient','geometry_2D','cartoon_geometry','cartoon_shift','reduction','matter','KO')
FIELDS=('chi','h11','h12','h22','hww','K','A11','A12','A22','Aww','Theta','Gamma1','Gamma2','lapse','shift1','shift2','B1','B2','phi','Pi','Lambda','Bx','By','Bz','Ex','Ey','Ez','Xi')
ODD=[2,7,12,15,17,22,23,25,26]
ODDOUT=[offset+c for offset in (0,N,2*N,3*N) for c in ODD]+[4*N+2*g+1 for g in range(G)]+[4*N+2*G+1,4*N+2*G+4]
KEYS=('lo_i','lo_j','hi_i','hi_j')
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def save(name,rows):
 with (HERE/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def build():
 b=module('t11build',HERE/'t7e-localize.py');b.TMP=TMP;b.build('T11RHS.cpp','rhs.ex')
def native(leg):
 start=time.monotonic();path=DATA/'exp-0020'/leg/'plt/EMS_Plot_000000.2d.hdf5';meta=[];performance=[]
 payload=TMP/'stencils.bin';out=TMP/'rhs.bin'
 with h5py.File(path) as f,payload.open('wb') as stream:
  assert float(f.attrs['time'])==0
  assert [f.attrs[f'component_{i}'].decode() for i in range(N)]==list(FIELDS)
  for l in range(7,13):
   g=f[f'level_{l}'];h=float(g.attrs['dx']);boxes=np.array([[b[k] for k in KEYS] for b in g['boxes'][:]],int)
   low=boxes[:,:2].min(0);high=boxes[:,2:].max(0);offset=g['data:offsets=0'][:]
   gx,gy=map(int,g['data_attributes'].attrs['outputGhost']);assert gx==gy==3
   count=0;read=0
   for bi,(x0,y0,x1,y1) in enumerate(boxes):
    xx,yy=np.meshgrid(np.arange(x0,x1+1),np.arange(y0,y1+1));x=(xx+.5)*h-336;y=(yy+.5)*h
    use=(x>=-5*h)&(x<=.5+5*h)&((y<=8*h)|(abs(x-y)<=8.1*h))
    jj,ii=np.where(use)
    if not len(ii):continue
    a=g['data:datatype=0'][int(offset[bi]):int(offset[bi+1])].reshape(34,y1-y0+7,x1-x0+7)[:N]
    read+=a.nbytes;count+=1
    # One box at a time; halo values are copied verbatim, never regenerated from static data.
    records=np.empty((len(ii),3+N*49));records[:,0]=h;records[:,1]=y[use];records[:,2]=x[use]
    for j in range(-3,4):
     for i in range(-3,4):records[:,3+(j+3)*7+i+3::49]=a[:,jj+3+j,ii+3+i].T
    records.tofile(stream)
    for a0,b0,xa,ya in zip(xx[use],yy[use],x[use],y[use]):
     meta.append((l,h,int(a0),int(b0),xa,ya,bi,min(a0-low[0],high[0]-a0,high[1]-b0),min(a0-x0,x1-a0,b0-y0,y1-b0)))
   performance.append(dict(leg=leg,level=l,boxes_read=count,boxes_available=len(boxes),bytes_read=read,h_M=h))
 subprocess.run([str(TMP/'rhs.ex'),'--t0-native',str(payload),str(out)],check=True,timeout=180)
 a=np.fromfile(out).reshape(-1,COL);m=np.array(meta);assert len(a)==len(m) and np.isfinite(a).all()
 np.savez_compressed(TMP/f'{leg}-native.npz',meta=m,q=a)
 payload.unlink();out.unlink()
 save(f't11-io-{leg}.csv',performance)
 print(leg,'cells',len(a),'seconds',time.monotonic()-start,'max term replay eps',a[:,-2].max(),'Base compute bit mismatches',a[:,-1].sum(),'peak RSS',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)

def basis(z,first,n):
 return np.array([math.prod((z-first-j)/(i-j) for j in range(n) if i!=j) for i in range(n)])
def sample(leg,ray,r,n=6):
 z=np.load(TMP/f'{leg}-native.npz');m=z['meta'];q=z['q'];vec=np.array([1.,0.]) if ray=='axis' else np.ones(2)/np.sqrt(2)
 result=np.empty((len(r),COL));levs=np.empty(len(r),int)
 maps={l:{(int(a[2]),int(a[3])):i for i,a in enumerate(m) if int(a[0])==l} for l in range(7,13)}
 hs={l:float(m[m[:,0]==l,1][0]) for l in range(7,13)}
 for k,R in enumerate(r):
  x,y=R*vec; l=min(12,max(7,int(math.floor(math.log2(112/max(x,y))))));h=hs[l]
  zx=(x+336)/h-.5;zy=y/h-.5;sx=math.floor(zx)-n//2+1;sy=math.floor(zy)-n//2+1
  # Use the finest level that contains the full point-value support; class rows separately retain face cells.
  while True:
   nodes=[(sx+i,sy+j if sy+j>=0 else -sy-j-1) for j in range(n) for i in range(n)]
   if all(key in maps[l] for key in nodes):break
   l-=1;assert l>=7,(leg,ray,R,n)
   h=hs[l];zx=(x+336)/h-.5;zy=y/h-.5;sx=math.floor(zx)-n//2+1;sy=math.floor(zy)-n//2+1
  wx=basis(zx,sx,n);wy=basis(zy,sy,n);anchor=q[maps[l][nodes[(n//2)*n+n//2]]].copy()
  if sy+n//2<0:anchor[ODDOUT]*=-1
  v=anchor.copy()
  for j in range(n):
   for i in range(n):
    a=q[maps[l][nodes[j*n+i]]].copy()
    if sy+j<0:a[ODDOUT]*=-1
    v+=(a-anchor)*wx[i]*wy[j]
  result[k]=v;levs[k]=l
 return result,levs

class Jet:
 """Analytic radial value/first/second jets, vectorized; no finite differences."""
 def __init__(self,v,d=0.,dd=0.):self.v=np.asarray(v);self.d=d;self.dd=dd
 def __add__(a,b):
  b=b if isinstance(b,Jet) else Jet(b);return Jet(a.v+b.v,a.d+b.d,a.dd+b.dd)
 __radd__=__add__
 def __neg__(a):return Jet(-a.v,-a.d,-a.dd)
 def __sub__(a,b):return a+-b if isinstance(b,Jet) else a+(-b)
 def __rsub__(a,b):return -a+b
 def __mul__(a,b):
  b=b if isinstance(b,Jet) else Jet(b);return Jet(a.v*b.v,a.d*b.v+a.v*b.d,a.dd*b.v+2*a.d*b.d+a.v*b.dd)
 __rmul__=__mul__
 def __pow__(a,p):return Jet(a.v**p,p*a.v**(p-1)*a.d,p*(p-1)*a.v**(p-2)*a.d*a.d+p*a.v**(p-1)*a.dd)
 def __truediv__(a,b):return a*(b**-1 if isinstance(b,Jet) else 1/b)
 def __rtruediv__(a,b):return a**-1*b
 def exp(a):
  e=np.exp(a.v);return Jet(e,e*a.d,e*(a.dd+a.d*a.d))

def continuum(xy,name,initial_lapse='original'):
 # This function is intentionally separate from every positive-time cache path.
 R=np.hypot(xy[:,0],xy[:,1]);assert (R>0).all();R.tofile(TMP/'r.bin')
 with (TMP/'cheb.bin').open('wb') as out,(TMP/'cheb.log').open('w') as log:
  subprocess.run([str(TMP/'rhs.ex'),'--t0-cheb',str(PROFILE),str(TMP/'r.bin')],stdout=out,stderr=log,check=True,timeout=60)
 a=np.fromfile(TMP/'cheb.bin').reshape(-1,20);s=a[:,1];nu=a[:,2];ks=1+(nu-1)*s;t=1-s
 sr=nu*s*t/(R*ks);srr=sr*sr*(1/s-1/t-(nu-1)/ks)-sr/R
 S=Jet(s,sr,srr);T=1-S;Y=Jet(a[:,8],a[:,9]*sr,a[:,10]*sr**2+a[:,9]*srr)
 D=Jet(a[:,12],a[:,13]*sr,a[:,14]*sr**2+a[:,13]*srr)
 P=Jet(a[:,16],a[:,17]*sr,a[:,18]*sr**2+a[:,17]*srr)
 metadata={k.strip():float(v) for line in PROFILE.read_text().splitlines() if line.startswith('# ') for k,sep,v in [line[2:].partition('=')] if sep and k in ('ell','C','phi_inf','q_native')}
 X=metadata['ell']*S**(1/nu)/Y;chi=X*X
 b=-metadata['C']*metadata['ell']*S**(1/nu)*T*T/(Y**3)
 curvature=metadata['C']*(T*T*D).exp()*T**3/Y**3
 phi=metadata['phi_inf']+T*P
 ps=Jet(-a[:,16]+t*a[:,17],(-2*a[:,17]+t*a[:,18])*sr,(-3*a[:,18]+t*a[:,19])*sr**2+(-2*a[:,17]+t*a[:,18])*srr)
 gg=Jet(a[:,8]+t*a[:,9],t*a[:,10]*sr,(-a[:,10]+t*a[:,11])*sr**2+t*a[:,10]*srr)
 # Keep Jet on the left: ndarray * Jet dispatches elementwise object operations.
 alpha_K=(-(T*T*D)).exp()*nu*S*gg/(Y*(1+S*(nu-1)))
 assert alpha_K.v.shape==R.shape and alpha_K.v.dtype!=object
 assert initial_lapse in ('original','maximal')
 alpha=X if initial_lapse=='original' else alpha_K
 pi=-metadata['C']*(T*T*D).exp()*T**4*ps/(Y**2*gg)
 coupling=(8*math.pi*.8*phi*phi).exp();electric=metadata['q_native']/(coupling*(Y/T)*Jet(R,np.ones_like(R)))
 n=xy/R[:,None];v=np.zeros((len(R),N));d1=np.zeros((len(R),N,2));d2=np.zeros((len(R),N,2,2))
 def scalar(c,z):
  v[:,c]=z.v;d1[:,c]=np.asarray(z.d)[...,None]*n
  for j in range(2):
   for k in range(2):d2[:,c,j,k]=z.dd*n[:,j]*n[:,k]+z.d/R*((1 if j==k else 0)-n[:,j]*n[:,k])
 def vector(c,z):
  f=z/Jet(R,np.ones_like(R))
  for i in range(2):
   v[:,c+i]=z.v*n[:,i]
   for j in range(2):
    d1[:,c+i,j]=z.d*n[:,i]*n[:,j]+z.v/R*((1 if i==j else 0)-n[:,i]*n[:,j])
    for k in range(2):d2[:,c+i,j,k]=xy[:,i]*((f.dd-f.d/R)*n[:,j]*n[:,k]+f.d/R*(1 if j==k else 0))+f.d*((1 if i==j else 0)*n[:,k]+(1 if i==k else 0)*n[:,j])
 for c,z in ((0,chi),(13,alpha),(18,phi),(19,pi)):scalar(c,z)
 for c,z in ((14,b),(16,-b),(24,electric)):vector(c,z)
 v[:,[1,3,4]]=1
 for c,i,j in ((6,0,0),(7,0,1),(8,1,1)):
  v[:,c]=curvature.v*(3*n[:,i]*n[:,j]-(1 if i==j else 0))
  for k in range(2):d1[:,c,k]=curvature.d*n[:,k]*(3*n[:,i]*n[:,j]-(1 if i==j else 0))+3*curvature.v/R*((1 if i==k else 0)*n[:,j]+(1 if j==k else 0)*n[:,i]-2*n[:,i]*n[:,j]*n[:,k])
 scalar(9,-curvature)
 adv=np.einsum('ni,nci->nc',v[:,14:16],d1)
 record=np.column_stack((np.zeros(len(R)),xy[:,1],xy[:,0],v,d1.reshape(len(R),-1),d2.reshape(len(R),-1),adv))
 record.tofile(TMP/'analytic.bin');subprocess.run([str(TMP/'rhs.ex'),'--t0-analytic',str(TMP/'analytic.bin'),str(TMP/'analytic-out.bin')],check=True,timeout=60)
 q=np.fromfile(TMP/'analytic-out.bin').reshape(-1,COL);assert np.isfinite(q).all()
 jets=np.stack([np.column_stack((z.v,np.broadcast_to(z.d,R.shape),np.broadcast_to(z.dd,R.shape))) for z in (X,alpha_K,b,curvature,chi,pi,phi,Y,gg)],axis=1)
 np.savez_compressed(TMP/f'{name}-continuum.npz',xy=xy,q=q,Killing_lapse=a[:,4],radial_jets=jets)
 return q

def rms(x):return float(np.sqrt(np.mean(np.asarray(x)**2)))
def projected(q,ray,offset,c):
 return q[:,offset+c] if ray=='axis' else (q[:,offset+c]+q[:,offset+c+1])/np.sqrt(2)
def measurements():
 from scipy.signal import find_peaks
 controls=[];termrows=[];variable=[];profile=[];peaks=[];convergence=[];continuumrows=[];locations=[];inner=[]
 for leg in ('E-low','E-mid','E-high'):
  z=np.load(TMP/f'{leg}-native.npz');m=z['meta'];q=z['q'];co=np.load(TMP/f'{leg}-continuum.npz')['q']
  controls.append(dict(leg=leg,cells=len(m),term_sum_defect_eps_max=q[:,-2].max(),Base_compute_bit_mismatches=int(q[:,-1].sum()),driver_cancellation_max=abs(q[:,138:140]).max(),lapse_sqrt_chi_max=abs(q[:,13]-np.sqrt(q[:,0])).max(),Gamma_t0_max=abs(q[:,11:13]).max(),floor_cells=int(((q[:,0]<=1e-12)|(q[:,13]<=1e-12)).sum()),finite=bool(np.isfinite(q).all())))
  for ray in ('axis','diagonal'):
   for l in range(7,13):
    take=(m[:,0]==l)&(m[:,4]>0)&(m[:,4]<=.5)&((m[:,3]==0) if ray=='axis' else ((abs(m[:,4]-m[:,5])<m[:,1]/100)&(np.hypot(m[:,4],m[:,5])<=.5)))
    ii=np.where(take)[0];ii=ii[np.argsort(m[ii,4])];r=m[ii,4] if ray=='axis' else np.hypot(m[ii,4],m[ii,5]);h=m[ii,1][0];dr=h*(1 if ray=='axis' else np.sqrt(2));a=q[ii];c=co[ii]
    gamma=projected(a,ray,N,11);gamma_c=projected(c,ray,N,11)
    uncovered=(np.maximum(m[ii,4],m[ii,5])>=112/2**(l+1)) if l<12 else np.ones(len(ii),bool)
    choose=np.where(uncovered)[0]
    if len(choose):
     k=int(choose[np.argmax(abs(gamma[choose]))]);locations.append(dict(leg=leg,ray=ray,level=l,uncovered_cells=int(uncovered.sum()),Gamma_rhs_peak=float(gamma[k]),peak_r_M=float(r[k]),outer_face_ray_M=112/2**l*(1 if ray=='axis' else np.sqrt(2)),axis_spacing_M=h,at_outer_face_Gamma_max=float(np.max(abs(gamma[abs(r-112/2**l*(1 if ray=='axis' else np.sqrt(2)))<3*dr]),initial=0))))
    if l==12:
     row=dict(leg=leg,ray=ray,h_M=h,r_M=r[0],r_over_h=r[0]/h,Gamma_rhs=float(gamma[0]),continuum_Gamma_rhs=float(gamma_c[0]),Gamma_error=float(gamma[0]-gamma_c[0]))
     for g,name in enumerate(GROUPS):
      row[name]=float(projected(a[:1],ray,4*N+2*g,0)[0]);row[name+'_error']=float(projected(a[:1],ray,4*N+2*g,0)[0]-projected(c[:1],ray,4*N+2*g,0)[0])
     row['A12_KO']=float(a[0,2*N+7]);row['h_times_A12_KO']=float(h*a[0,2*N+7]);row['lapse_KO']=float(a[0,2*N+13]);inner.append(row)
    for g,name in enumerate(GROUPS):
     u=projected(a,ray,4*N+2*g,0);k=int(abs(u).argmax());cu=(-u[:-4]+16*u[1:-3]-30*u[2:-2]+16*u[3:-1]-u[4:])/(12*dr*dr)
     kc=int(abs(cu).argmax());aa=bb=kc
     while aa>0 and abs(cu[aa-1])>=abs(cu[kc])/2:aa-=1
     while bb<len(cu)-1 and abs(cu[bb+1])>=abs(cu[kc])/2:bb+=1
     termrows.append(dict(leg=leg,ray=ray,level=l,term=name,peak_signed=u[k],peak_abs=abs(u[k]),peak_r_M=r[k],peak_r_over_h=r[k]/h,curvature_peak=float(abs(cu[kc])),curvature_peak_r_M=float(r[kc+2]),curvature_width_cells=bb-aa+1,curvature_endpoint_lobe=bool(aa==0 or bb==len(cu)-1),inner_3cell_max=float(np.max(abs(u[r<3*dr]),initial=0)),outer_face_M=112/2**l*(1 if ray=='axis' else np.sqrt(2)),at_face_max=float(np.max(abs(u[abs(r-112/2**l*(1 if ray=='axis' else np.sqrt(2)))<3*dr]),initial=0)),RMS=rms(u),continuum_peak_abs=float(np.max(abs(projected(c,ray,4*N+2*g,0))))))
    for j,R in enumerate(r):
     row=dict(leg=leg,ray=ray,level=l,h_M=h,r_M=R,covered_by_finer=bool(l<12 and max(m[ii[j],4],m[ii[j],5])<112/2**(l+1)),near_face=bool(m[ii[j],7]<3),near_box_seam=bool(m[ii[j],8]<3),Gamma_rhs=gamma[j],continuum_Gamma_rhs=gamma_c[j])
     row.update({name:float(projected(a[j:j+1],ray,4*N+2*g,0)[0]) for g,name in enumerate(GROUPS)})
     row.update(beta_advection=float(projected(a[j:j+1],ray,3*N,14)[0]),beta_driver=float(projected(a[j:j+1],ray,138,0)[0]),beta_KO=float(projected(a[j:j+1],ray,2*N,14)[0]),B_decay=float(projected(a[j:j+1],ray,141,0)[0]),B_KO=float(projected(a[j:j+1],ray,2*N,16)[0]),lapse_advection=float(a[j,3*N+13]),lapse_slicing=float(a[j,140]),lapse_KO=float(a[j,2*N+13]))
     profile.append(row)
    # Native radial curvature: boundary maxima are explicitly endpoint lobes, not a fitted FWHM.
    d=np.full(len(r),np.nan);d[2:-2]=(-gamma[:-4]+16*gamma[1:-3]-30*gamma[2:-2]+16*gamma[3:-1]-gamma[4:])/(12*dr*dr)
    valid=np.where(np.isfinite(d))[0];k=int(valid[np.argmax(abs(d[valid]))]);height=abs(d[k]);a0=b0=k
    while a0>valid[0] and abs(d[a0-1])>=height/2:a0-=1
    while b0<valid[-1] and abs(d[b0+1])>=height/2:b0+=1
    peaks.append(dict(leg=leg,ray=ray,level=l,h_M=h,curvature_peak=abs(d[k]),peak_r_M=r[k],peak_r_over_h=r[k]/h,width_native_cells=b0-a0+1,width_M=(b0-a0+1)*dr,endpoint_lobe=bool(k==valid[0] or a0==valid[0]),Gamma_rhs_peak=float(np.max(abs(gamma))),Gamma_rhs_peak_r_M=float(r[np.argmax(abs(gamma))]),native_continuum_RMS=rms(gamma-gamma_c)))
    for j,name in enumerate(FIELDS):
     p=a[:,N+j];ko=a[:,2*N+j];adv=a[:,3*N+j];variable.append(dict(leg=leg,ray=ray,level=l,variable=name,physical_max=float(np.max(abs(p))),KO_max=float(np.max(abs(ko))),advection_max=float(np.max(abs(adv))),physical_peak_r_M=float(r[np.argmax(abs(p))]),KO_peak_r_M=float(r[np.argmax(abs(ko))]),h_times_KO_max=float(h*np.max(abs(ko))),continuum_max=float(np.max(abs(c[:,N+j])))))
 save('t11-harness-controls.csv',controls);save('t11-term-peaks.csv',termrows);save('t11-variable-peaks.csv',variable);save('t11-native-profiles.csv',profile);save('t11-source-widths.csv',peaks);save('t11-uncovered-sources.csv',locations);save('t11-inner-cell-terms.csv',inner)
 # Fixed physical windows; omit interpolation supports that touch any union face.
 for ray in ('axis','diagonal'):
  r=np.unique(np.r_[np.geomspace(.0002,.025,180),np.linspace(.026,.49,350)])
  samples=[sample(leg,ray,r,6)[0] for leg in ('E-low','E-mid','E-high')]
  samples8=[sample(leg,ray,r,8)[0] for leg in ('E-low','E-mid','E-high')]
  keep=np.ones(len(r),bool)
  for l in range(7,13):keep &= abs(r-112/2**l*(1 if ray=='axis' else np.sqrt(2)))>5*(.875/2**l)*(1 if ray=='axis' else np.sqrt(2))
  for lo,hi in ((.0002,.001),(.001,.005),(.005,.02),(.02,.05),(.05,.1),(.1,.49)):
   use=keep&(r>=lo)&(r<=hi)
   for label,c in [('Gamma1',11),('shift1',14),('lapse',13),('chi',0),('K',5),('Theta',10),('h11',1),('A11',6),('B1',16)]:
    vals=[projected(a,ray,N,c) if c in (11,14,16) else a[:,N+c] for a in samples]
    vals8=[projected(a,ray,N,c) if c in (11,14,16) else a[:,N+c] for a in samples8]
    d01=rms((vals[0]-vals[1])[use]);d12=rms((vals[1]-vals[2])[use]);spread=max(rms((v-v8)[use]) for v,v8 in zip(vals,vals8))
    convergence.append(dict(ray=ray,variable=label,r_min_M=lo,r_max_M=hi,points=int(use.sum()),low_mid_RMS=d01,mid_high_RMS=d12,order=math.log(d01/d12)/math.log(1.5) if d01*d12>0 else np.nan,P6_P8_RMS_max=spread,interpolation_fraction=spread/d12 if d12 else np.nan))
  # Continuum at exactly the common physical points; axis is approached analytically through y>0.
  xy=np.column_stack((r,r*0+1e-14)) if ray=='axis' else np.column_stack((r/np.sqrt(2),r/np.sqrt(2)))
  c=continuum(xy,'fixed-'+ray)
  for leg,a,a8 in zip(('E-low','E-mid','E-high'),samples,samples8):
   for lo,hi in ((.0002,.001),(.001,.005),(.005,.02),(.02,.05),(.05,.1),(.1,.49)):
    use=keep&(r>=lo)&(r<=hi);g=projected(a,ray,N,11);gc=projected(c,ray,N,11);g8=projected(a8,ray,N,11)
    continuumrows.append(dict(leg=leg,ray=ray,r_min_M=lo,r_max_M=hi,points=int(use.sum()),native_RMS=rms(g[use]),continuum_RMS=rms(gc[use]),error_RMS=rms((g-gc)[use]),P6_P8_RMS=rms((g-g8)[use]),state_alpha_error_max=float(np.max(abs(a[use,13]-c[use,13]))),state_A_error_max=float(np.max(abs(a[use,6:10]-c[use,6:10])))))
 save('t11-rhs-orders.csv',convergence);save('t11-continuum-errors.csv',continuumrows)

def validation():
 from scipy.interpolate import CubicHermiteSpline
 rows=[];profiles=[];p=DATA/'exp-0022/E-L/rays';dt=7/384
 for ray in ('axis','diagonal'):
  states=[np.load(p/f'E-L-step{i:06d}-{ray}.npz') for i in range(3)]
  r=np.geomspace(.0002,.49,500);a,l=sample('E-mid',ray,r,6);a8,_=sample('E-mid',ray,r,8)
  source_names=list(states[0]['native_fields'])
  for name,c in [('shift1',14),('Gamma1',11),('B1',16),('lapse',13),('chi',0),('K',5),('Theta',10),('h11',1),('h12',2),('h22',3),('hww',4)]:
   j=source_names.index(name);u=[]
   for z in states:
    values=z['native'][0,j];derivative=z['native'][1,j]
    if c in (11,14,16) and ray=='diagonal':
     values=(values+z['native'][0,j+1])/np.sqrt(2);derivative=(derivative+z['native'][1,j+1])/np.sqrt(2)
    u.append(CubicHermiteSpline(z['r'],values,derivative)(r))
   rhs=projected(a,ray,N,c)+projected(a,ray,2*N,c) if c in (11,14,16) else a[:,N+c]+a[:,2*N+c]
   rhs8=projected(a8,ray,N,c)+projected(a8,ray,2*N,c) if c in (11,14,16) else a8[:,N+c]+a8[:,2*N+c]
   v0=projected(a,ray,0,c) if c in (11,14,16) else a[:,c]
   delta=u[1]-u[0];second=u[2]-2*u[1]+u[0];residual=delta-dt*rhs;corrected=residual-.5*second
   keep=np.ones(len(r),bool)
   for level in range(7,13):keep &= abs(r-112/2**level*(1 if ray=='axis' else np.sqrt(2)))>5*(7/12/2**level)*(1 if ray=='axis' else np.sqrt(2))
   for lo,hi in ((.0002,.005),(.005,.02),(.02,.05),(.05,.1),(.1,.49)):
    use=keep&(r>=lo)&(r<=hi);leading=rms((dt*rhs)[use]);curvature=rms((.5*second)[use]);err=rms(residual[use]);corr=rms(corrected[use]);spread=rms((dt*(rhs-rhs8))[use])
    rows.append(dict(ray=ray,variable=name,r_min_M=lo,r_max_M=hi,cells=int(use.sum()),dt_M=dt,dt_RHS_RMS=leading,first_change_RMS=rms(delta[use]),first_minus_dt_RHS_RMS=err,Odt2_from_steps12_RMS=curvature,corrected_residual_RMS=corr,residual_over_Odt2=err/curvature if curvature else np.nan,corrected_over_first= corr/rms(delta[use]) if rms(delta[use]) else np.nan,dt_RHS_interpolation_RMS=spread,initial_state_cache_error_RMS=rms((v0-u[0])[use]),Taylor_condition='OUTSIDE_2DT' if lo>2*dt else 'CADENCE_TOO_LARGE_FOR_SOURCE'))
   if name in ('Gamma1','shift1','lapse'):
    for i,R in enumerate(r):profiles.append(dict(ray=ray,variable=name,r_M=R,delta1=delta[i],dt_RHS=dt*rhs[i],Odt2=.5*second[i],corrected_residual=corrected[i]))
 save('t11-step-validation.csv',rows);save('t11-step-profiles.csv',profiles)
 launch=[]
 raw=list(csv.DictReader((DATA/'exp-0022/E-L/launch-native-widths.csv').open()))
 values=list(csv.DictReader((DATA/'exp-0022/E-L/launch-native.csv').open()))
 for row in raw:
  if row['role']!='front' or not 1<=int(row['coarse_step'])<=6:continue
  step=int(row['coarse_step']);ray=row['ray'];t=float(row['time_M']);v=[x for x in values if int(x['coarse_step'])==step and x['ray']==ray and x['role']=='front']
  k=min(v,key=lambda a:abs(float(a['r_M'])-float(row['dominant_curvature_peak_M'])));rr=float(k['r_M'])
  d=dict(row);d.update(Gamma_at_curvature_peak=float(k['Gamma_n']),curvature_peak=float(k['Gamma_n_d2_native']),C_Gamma_at_peak=float(k['C_Gamma']),backtracked_radius_speed1_M=rr-t)
  launch.append(d)
 save('t11-launch.csv',launch)
 independent=[]
 for ray in ('axis','diagonal'):
  for step in range(1,7):
   z=np.load(p/f'E-L-step{step:06d}-{ray}.npz');r=z['r'];take=np.where((r>=.003)&(r<=.16))[0];k=int(take[np.argmax(abs(z['q'][2,1,take]))]);k8=int(take[np.argmax(abs(z['q8'][2,1,take]))])
   l=int(z['levels'][k]);h=7/12/2**l;fac=1 if ray=='axis' else np.sqrt(2);distance=min(abs(r[k]-112/2**a*fac) for a in range(7,13))
   independent.append(dict(ray=ray,step=step,time_M=float(z['time_M']),level=l,r_peak_M=r[k],r_peak_P8_M=r[k8],P6_P8_position_spread_M=abs(r[k]-r[k8]),Gamma_curvature=z['q'][2,1,k],Gamma_amplitude=z['q'][0,1,k],C_Gamma=z['q'][0,4,k],C_Gamma_share=abs(z['q'][0,4,k]/z['q'][0,1,k]),nearest_face_M=distance,clean=bool(distance>4*h*fac),selection='GLOBAL_ABS_GAMMA_CURVATURE_FIXED_0.003_0.16_M'))
 save('t11-independent-launch.csv',independent)

def figures():
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 try:
  import scienceplots
  plt.style.use(['science','no-latex'])
 except ImportError:plt.rcParams.update({'font.family':'serif','font.size':10})
 out=HERE/'figures';out.mkdir(exist_ok=True)
 def finish(fig,name):
  fig.savefig(out/f'{name}.png',dpi=220);fig.savefig(out/f'{name}.pdf');plt.close(fig)
 z=np.load(TMP/'radial-continuum.npz');r=np.hypot(*z['xy'].T);a=z['q']
 data=list(csv.DictReader((HERE/'t11-native-profiles.csv').open()))
 fig,axs=plt.subplots(2,2,figsize=(9,7),layout='constrained')
 for ai,ray in enumerate(('axis','diagonal')):
  for l in range(7,13):
   v=[x for x in data if x['leg']=='E-mid' and x['ray']==ray and int(x['level'])==l]
   rr=np.array([float(x['r_M']) for x in v]);u=np.array([float(x['Gamma_rhs']) for x in v]);axs[ai,0].plot(rr,u,label=f'L{l}')
  axs[ai,0].plot(r,a[:,39:41].sum(1)/np.sqrt(2),'k--',label='continuum, t=0');axs[ai,0].set_xscale('log');axs[ai,0].set_xlim(1e-4,.5);axs[ai,0].set_ylim(-15,440);axs[ai,0].set_title(f'{ray}: native Γ RHS by level');axs[ai,0].set_xlabel('r/M');axs[ai,0].set_ylabel('∂t Γ normal')
  v=[x for x in data if x['leg']=='E-mid' and x['ray']==ray and x['level']=='12'];rr=np.array([float(x['r_M']) for x in v])
  for name in ('lapse_A','chi_A','shift_laplacian_2D','shift_graddiv_2D','cartoon_shift','matter'):axs[ai,1].plot(rr,[float(x[name]) for x in v],label=name.replace('_',' '))
  axs[ai,1].set_xlim(0,.006);axs[ai,1].set_ylim(-160,330);axs[ai,1].set_title(f'{ray}: L12 term split');axs[ai,1].set_xlabel('r/M')
 axs[0,0].legend(ncol=3,fontsize=8);axs[0,1].legend(ncol=2,fontsize=8);finish(fig,'t11-term-profiles')
 fig,axs=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
 for leg,color in zip(('E-low','E-mid','E-high'),('C0','C1','C2')):
  v=[x for x in data if x['leg']==leg and x['ray']=='diagonal' and x['level']=='12'];rr=np.array([float(x['r_M']) for x in v]);u=np.array([float(x['Gamma_rhs']) for x in v]);co=np.array([float(x['continuum_Gamma_rhs']) for x in v]);h=float(v[0]['h_M'])
  axs[0].plot(rr,u,color=color,label=leg);axs[1].plot(rr/h,u-co,color=color,label=leg)
  raw=np.load(TMP/f'{leg}-native.npz');m=raw['meta'];q=raw['q'];use=(m[:,0]==12)&(m[:,4]>0)&(abs(m[:,4]-m[:,5])<m[:,1]/100);rr=np.hypot(m[use,4],m[use,5]);order=np.argsort(rr);axs[2].plot(rr[order]/h,h*q[use,2*N+7][order],color=color,label=leg)
 axs[0].plot(r,a[:,39:41].sum(1)/np.sqrt(2),'k--',label='continuum');axs[0].set_xlim(0,.003);axs[0].set_ylim(0,430);axs[0].set_xlabel('r/M');axs[0].set_ylabel('Γ RHS');axs[0].legend(fontsize=8)
 axs[1].set_xlim(0,15);axs[1].set_ylim(-35,15);axs[1].set_xlabel('r/h12');axs[1].set_ylabel('native − continuum Γ RHS')
 axs[2].set_xlim(0,12);axs[2].set_ylim(-.5,4.5);axs[2].set_xlabel('r/h12');axs[2].set_ylabel('h12 × KO(A12)');finish(fig,'t11-continuum-puncture')
 conv=list(csv.DictReader((HERE/'t11-rhs-orders.csv').open()));fig,axs=plt.subplots(1,2,figsize=(8,3.5),layout='constrained')
 for ai,ray in enumerate(('axis','diagonal')):
  for name in ('Gamma1','lapse','K','h11'):
   v=[x for x in conv if x['ray']==ray and x['variable']==name];xx=np.arange(len(v));axs[ai].plot(xx,[float(x['order']) for x in v],'.-',label=name)
  axs[ai].axhline(4,color='k',ls='--');axs[ai].set_xticks(range(6),['.0002–.001','.001–.005','.005–.02','.02–.05','.05–.1','.1–.49'],rotation=45,ha='right');axs[ai].set_ylim(-4,7);axs[ai].set_title(ray);axs[ai].set_ylabel('RHS self-difference order');axs[ai].set_xlabel('fixed radial window / M')
 axs[0].legend(fontsize=8);finish(fig,'t11-rhs-convergence')
 val=list(csv.DictReader((HERE/'t11-step-validation.csv').open()));fig,axs=plt.subplots(1,2,figsize=(8,3.6),layout='constrained')
 for ai,ray in enumerate(('axis','diagonal')):
  for name in ('B1','lapse','K','Gamma1','shift1','h11'):
   v=[x for x in val if x['ray']==ray and x['variable']==name];axs[ai].plot(np.arange(len(v)),[float(x['residual_over_Odt2']) for x in v],'.-',label=name)
  axs[ai].set_yscale('log');axs[ai].set_ylim(.01,1000);axs[ai].axhline(1,color='k',ls='--');axs[ai].set_title(ray);axs[ai].set_xticks(range(5),['.0002–.005','.005–.02','.02–.05','.05–.1','.1–.49'],rotation=45,ha='right');axs[ai].set_ylabel('|Δu − Δt RHS| / |½ Δ²u|');axs[ai].set_xlabel('fixed radial window / M')
 axs[0].legend(ncol=2,fontsize=8);finish(fig,'t11-step-validation')

def check():
 controls=list(csv.DictReader((HERE/'t11-harness-controls.csv').open()))
 assert all(float(x['term_sum_defect_eps_max'])<4 and int(x['Base_compute_bit_mismatches'])==0 and float(x['driver_cancellation_max'])==0 and x['finite']=='True' and int(x['floor_cells'])==0 for x in controls)
 for leg in ('E-low','E-mid','E-high'):
  z=np.load(TMP/f'{leg}-native.npz');co=np.load(TMP/f'{leg}-continuum.npz')['q'];m=z['meta'];q=z['q'];assert np.isfinite(q).all() and np.isfinite(co).all();assert np.all(q[:,13]>1e-12)
  # t=0 reconstruction must reproduce the numerical setter, including all A/EMS fields.
  use=m[:,4]>.001;keep=[0,1,2,3,4,6,7,8,9,13,14,15,16,17,18,19,24,25]
  assert np.max(abs(q[use][:,keep]-co[use][:,keep]))<1e-9
 assert Jet(np.array([2.]),1.,0.).__pow__(3).dd==12
 print('T11 native kernel/term/driver/analytic-setter checks PASS')

def ghost_check():
 rows=[]
 for leg in ('E-low','E-mid','E-high'):
  with h5py.File(DATA/'exp-0020'/leg/'plt/EMS_Plot_000000.2d.hdf5') as f:
   assert float(f.attrs['time'])==0
   for l in range(7,13):
    g=f[f'level_{l}'];h=float(g.attrs['dx']);boxes=np.array([[b[k] for k in KEYS] for b in g['boxes'][:]],int);offset=g['data:offsets=0'][:];cache={};counts={s:[0,0,0,0.] for s in ('axis_parity','same_level_copy')}
    def block(bi):
     if bi not in cache:
      x0,y0,x1,y1=boxes[bi];cache[bi]=g['data:datatype=0'][int(offset[bi]):int(offset[bi+1])].reshape(34,y1-y0+7,x1-x0+7)[:N]
     return cache[bi]
    def owner(i,j):
     choose=np.where((boxes[:,0]<=i)&(boxes[:,2]>=i)&(boxes[:,1]<=j)&(boxes[:,3]>=j))[0];assert len(choose)==1;return int(choose[0])
    center=int(round(336/h))
    for i0 in range(center,center+3):
     for j0 in range(3):
      bi=owner(i0,j0);b=boxes[bi];a=block(bi)
      for dj in range(-3,4):
       for di in range(-3,4):
        i,j=i0+di,j0+dj
        if b[0]<=i<=b[2] and b[1]<=j<=b[3]:continue
        jf=-j-1 if j<0 else j;oi=owner(i,jf);ob=boxes[oi];expected=block(oi)[:,jf-ob[1]+3,i-ob[0]+3].copy()
        if j<0:expected[ODD]*=-1
        actual=a[:,j-b[1]+3,i-b[0]+3];values=counts['axis_parity' if j<0 else 'same_level_copy'];values[0]+=N;values[1]+=int(np.sum(actual.view('u8')!=expected.view('u8')));values[2]+=int(np.sum(actual!=expected));values[3]=max(values[3],float(abs(actual-expected).max()))
    for kind,(values,bits,numeric,maximum) in counts.items():rows.append(dict(leg=leg,level=l,kind=kind,compared_values=values,bit_mismatches=bits,numeric_mismatches=numeric,max_difference=maximum))
 save('t11-ghost-controls.csv',rows)
 assert all(x['numeric_mismatches']==0 for x in rows),rows
 print('T11 puncture parity/same-level saved ghost controls PASS')

def provenance():
 rows=[];sub=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0020/submissions/exp-0020')
 for leg in ('E-low','E-mid','E-high'):
  params={line.split('=',1)[0].strip():line.split('=',1)[1].split('#',1)[0].strip() for line in (sub/f'params-{leg}.txt').read_text().splitlines() if '=' in line and not line.strip().startswith('#')}
  for key,value in dict(amr_transfer='point',sigma='1',eta='1',kappa1='0.1',kappa2='0',kappa3='1',formulation='0',num_ghosts='3').items():assert params[key]==value,(leg,key,params[key])
  path=DATA/'exp-0020'/leg/'plt/EMS_Plot_000000.2d.hdf5';m=np.load(TMP/f'{leg}-native.npz')['meta']
  with h5py.File(path) as f:
   assert float(f.attrs['time'])==0
   for l in range(7,13):
    g=f[f'level_{l}'];off=g['data:offsets=0'][:];boxes=g['boxes'][:];digest=hashlib.sha256();ids=np.unique(m[m[:,0]==l,6]).astype(int)
    for bi in ids:
     x0,y0,x1,y1=[int(boxes[bi][k]) for k in KEYS];a=g['data:datatype=0'][int(off[bi]):int(off[bi+1])].reshape(34,y1-y0+7,x1-x0+7)[:N];digest.update(a.tobytes())
    rows.append(dict(leg=leg,level=l,time_M=0.,input_path=str(path),input_file_bytes=path.stat().st_size,boxes=';'.join(map(str,ids)),selected_evolved_fields_sha256=digest.hexdigest(),params_sha256=hashlib.sha256((sub/f'params-{leg}.txt').read_bytes()).hexdigest(),sigma=params['sigma'],amr_transfer=params['amr_transfer']))
 save('t11-input-provenance.csv',rows)
 print('T11 params and selected input-box provenance PASS')

def report():
 def rows(name):return list(csv.DictReader((HERE/name).open()))
 def table(headers,data):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in data])+'\n\n'
 def sci(x):return f'{float(x):.6e}'
 controls=rows('t11-harness-controls.csv');terms=rows('t11-term-peaks.csv');inner=rows('t11-inner-cell-terms.csv');width=rows('t11-source-widths.csv');orders=rows('t11-rhs-orders.csv');val=rows('t11-step-validation.csv');launch=rows('t11-launch.csv');events=rows('t11-independent-launch.csv');cross=json.loads((HERE.parents[1]/'scripts/cas/t11-continuum-crosscheck.json').read_text())
 out=['''## T11 — initialized longitudinal RHS and puncture-end source

**READY-EXCEPT:** the t=0 source is localized to the puncture-end cells, and its leading term is identified. It is a nonregular geometric source already present in the continuum RHS of the **actual initialized state**, with an additional O(1) lapse-gradient discretization error in the first cells and a 1/h KO source on the angular extrinsic-curvature components. The shift-driver combination cancels exactly. The native harness is bit-identical to a direct call of the unmodified production kernel. The supplied first coarse-step snapshots cannot certify the inner-cell Taylor expansion or determine how much of the emitted pulse is seeded by geometric Gamma versus KO on A: they skip 4096 finest-level steps before their first sample. This qualification is not deferred as an unreported validation pass.

This is an offline, read-only diagnostic on HEAD bb2d9f3. No production C++, initialization, evolution, gauge, Chombo, finder, parameters or input file is changed, and no commit or evolution run is made. All static-file operations are explicitly t=0. Positive-time analysis reads only E-L's saved numerical fields. Four standalone scientific figures use the sciplot workflow; this report follows the sciwrite evidence/limitation convention.

### Native kernel, saved ghosts and exact controls

[T11RHS.cpp](T11RHS.cpp) is a Tests-only subclass of `CCZ4Cartoon<ExperimentalGauge,FourthOrderDerivatives,CouplingFunction>`. It calls the production protected RHS on native derivatives, separates KO using the production dissipation method, and independently groups the Gamma equation's contributions. It applies the same trace-A removal and 1e-12 chi/lapse floors to the entire input stencil, including ghosts, as `specificEvalRHS`. For a second control it calls the public production `compute` and compares all 28 outputs bit for bit with the native RHS plus sequential native dissipation additions. Arithmetic order matters when two KO directions cancel; the bit control uses the production addition order, whereas separately reported KO is the sum of the two directions.

All three exp-0020 t=0 plots carry every one of the 28 evolved variables, including A and EMS fields, plus six constraints, and three saved ghost cells in both directions. The supplied exp-0019 E-high checkpoint also carries the 28 evolved variables and three ghosts, but has centre 224 and a different hierarchy; it is not substituted for the E-mid source state with centre 336. Plot time is checked to equal zero before any native or continuum preparation. Only selected ray-support boxes on levels 7–12 are read; no whole level is loaded. The streamed seven-by-seven stencils include current-field parity and coarse/fine ghosts exactly as stored. Float64, point transfers, sigma=1, G=1, eta=1, shift Gamma coefficient 0.75 and kappa=(0.1,0,1) match the supplied parameters. The native covariant damping actually uses 0.1*alpha/(0.005+alpha), as in the frozen source.
''']
 out.append(table(['rung','native cells','production output bit mismatches','Gamma split defect / ε max(1,sum abs terms)','driver cancellation max','floor cells'],[[r['leg'],r['cells'],r['Base_compute_bit_mismatches'],f"{float(r['term_sum_defect_eps_max']):.6f}",r['driver_cancellation_max'],r['floor_cells']] for r in controls]))
 out.append('''There are 82,807 sampled cells and 2,318,596 production-output components. All are finite. Alpha equals sqrt(chi) bit for bit in the stored numerical fields; the largest initial Gamma component is 1.315519e-12. The split error is at most 2.494 machine epsilons times the stated absolute-term scale. [t11-harness-controls.csv](t11-harness-controls.csv) retains the controls. A separate **108,864-value** comparison of the first three puncture-cell stencil halos with valid same-level owners and exact axis parity has **zero numeric and zero bit differences**, on every rung and level 7–12; see [t11-ghost-controls.csv](t11-ghost-controls.csv). Thus the first-cell source does not result from a wrong parity fill or a same-level seam copy. Its stencil is far from a coarse/fine face.

### Actual initialization and term decomposition

The stored trumpet Killing lapse is not the initialized lapse. `EMSBH_trumpet_read::conformal_decomposition` sets alpha(0)=sqrt(chi), even though the radial reader also reconstructs alpha_K. The geometric data have h_ij=delta_ij, chi=X², K=0 and A_ij=k(3n_i n_j−delta_ij) in the unboosted continuum reconstruction. Numerical metric/K deviations are at floating-point precision. This is a correction to the assumed stationary-lapse picture, not a proposed initialization or gauge change.

`ExperimentalGauge::compute` sets B_driver=0.75*Gamma−eta*beta. Consequently 0.75*Gamma−eta*beta−B_driver is **exactly zero** on every native sampled cell at t=0. The actual gauge equations are alpha_t=advec(alpha)−1.8*alpha*(K−2Theta), beta_t=advec(beta)+0.75*Gamma−eta*beta−B_driver, and B_driver,t=−0.1*B_driver, with native KO added to each. There is no advection in the driver equation, although the generic derivative object computes an unused driver-advection value. The lapse coefficient in the parameter file does not replace the hard-coded 1.8 in ExperimentalGauge.

The table gives signed extrema of Gamma-normal contributions over 0<r<=0.5 M on **E-mid L12**. Native axis profiles use the closest row y=h/2 and label its horizontal coordinate x−336 as r, matching E-L's native provenance; the same-cell continuum comparison uses the actual (x−336,y). Diagonal profiles use x−336=y and r=sqrt(2)*y. Fixed physical P6/P8 axis sampling later uses y=0. Cartesian 2D Laplacian and gradient/divergence terms plus cartoon terms together form the 3D shift-second-derivative block; the longitudinal coefficient is 4/3, not an extra 4/3 multiplying the full Laplacian. The reported curvature width counts connected native centres above half the absolute curvature maximum, using the native fourth-order radial difference of each RHS contribution. `E` denotes an endpoint lobe whose inner side is not bracketed; it is not a complete FWHM. Widths of pure roundoff/zero terms are not physical packets.
''')
 data=[]
 for name in GROUPS:
  a=next(x for x in terms if x['leg']=='E-mid' and x['ray']=='axis' and x['level']=='12' and x['term']==name);d=next(x for x in terms if x['leg']=='E-mid' and x['ray']=='diagonal' and x['level']=='12' and x['term']==name)
  w=lambda x:'—' if float(x['peak_abs'])<1e-8 else x['curvature_width_cells']+(' E' if x['curvature_endpoint_lobe']=='True' else '')
  data.append([name,sci(a['peak_signed']),sci(a['peak_r_M']),sci(d['peak_signed']),sci(d['peak_r_M']),w(a)+' / '+w(d)])
 out.append(table(['Gamma term','axis extremum','axis r/M','diagonal extremum','diagonal r/M','curvature centres axis / diagonal'],data))
 out.append('''The dominant positive Gamma terms are −2*A^ij*partial_j alpha and −3*alpha*A^ij*partial_j chi/chi. The shift-second-derivative block contributes the opposing fractional-power correction. K/Theta gradients, Christoffels of h and hww, reduction damping and Gamma advection have no material initial narrow source. The scalar-momentum matter term is nonzero, but its peak is a broader feature near 0.005 M and is much smaller than the puncture source. All signed profiles, including beta advection/driver/KO, driver decay/KO and lapse advection/slicing/KO, are in [t11-native-profiles.csv](t11-native-profiles.csv); every term/level/ray extremum and width is in [t11-term-peaks.csv](t11-term-peaks.csv). The frozen kernel evaluates all 28 RHS variables; [t11-variable-peaks.csv](t11-variable-peaks.csv) retains physical/advection/KO extrema, including every h/A and EMS component.

The following are **same-cell comparisons** at the closest positive diagonal cell. They use the continuum source at that cell's actual radius, not a fit or a target for evolution.
''')
 out.append(table(['rung','h12/M','r/M','Gamma RHS','continuum Gamma RHS','native − continuum','lapse-gradient error','KO(A12)','h12*KO(A12)','KO(alpha)'],[[x['leg'],sci(x['h_M']),sci(x['r_M']),f"{float(x['Gamma_rhs']):.6f}",f"{float(x['continuum_Gamma_rhs']):.6f}",f"{float(x['Gamma_error']):.6f}",f"{float(x['lapse_A_error']):.6f}",sci(x['A12_KO']),f"{float(x['h_times_A12_KO']):.6f}",f"{float(x['lapse_KO']):.6f}"] for x in inner if x['ray']=='diagonal']))
 out.append('''The approximately −30.3 lapse-gradient error stays O(1) as the first cell moves inward with h. This is the fourth-order centred stencil differentiating the conical alpha~r/Y0, not an error at a fixed nonzero physical radius. It occupies the first few cells. Native Gamma KO is only about 1e-9–1e-8 and cannot explain a Gamma RHS of hundreds. **KO on A12 is different:** h*KO(A12) approaches approximately 4.14, and its peak grows as 1/h. A_ij has a bounded but direction-dependent puncture limit; applying Cartesian KO to those angular components creates a singular grid-scale RHS. The corresponding source on lapse is O(1), approximately 0.342. E-mid's diagonal L12 maxima of KO(phi), KO(Pi) and KO(Ex) are 0.0837899, 0.943445 and 119608.4 respectively; the covariant electric field also has a singular puncture representation. These are real current-field stencil effects, not new continuum damping or targets. They demonstrate that a Gamma-only KO audit would miss large coupled sources.

### Where the source is, and how it scales

The first Gamma RHS curvature lobe is at 2.5h on the axis and 3.5355h on the diagonal, with one/two counted native centres; both are endpoint lobes. The first two radial stencil points cannot supply a centred five-point second difference, so these are **first evaluable** curvature locations, not proven locations of the exact maximum. The source itself is present at the innermost cell before any update. Its extremum and its curvature are different measurements.
''')
 out.append(table(['rung','axis max Gamma RHS','diagonal max Gamma RHS','axis curvature peak','diagonal curvature peak','axis / diagonal lobe centres'],[[leg,*(sci(next(x for x in width if x['leg']==leg and x['ray']==ray and x['level']=='12')[key]) for key in ('Gamma_rhs_peak','curvature_peak') for ray in ('axis','diagonal')),'1 / 2 E'] for leg in ('E-low','E-mid','E-high')]))
 loc=rows('t11-uncovered-sources.csv');data=[]
 for l in range(7,13):
  pair=[next((x for x in loc if x['leg']=='E-mid' and x['ray']==ray and int(x['level'])==l),None) for ray in ('axis','diagonal')]
  data.append([l,f'{112/2**l:.9f}',*[sci(x['Gamma_rhs_peak']) if x else 'no uncovered cells at r≤0.5' for x in pair]])
 out.append(table(['level','positive axial face/M','E-mid uncovered axis Gamma extremum','E-mid uncovered diagonal Gamma extremum'],data))
 out.append('''Actual inner faces are R12=0.02734375, R11=0.0546875, R10=0.109375, R9=0.21875 and R8=0.4375 M; R7=0.875 M is outside the requested interval. Their diagonal corners are sqrt(2) times these values. E-mid's largest Gamma RHS in the three-cell collar of R12 is only 0.00607142 on the axis and 0.00309711 at its diagonal corner, against inner maxima 249.705 and 276.898. Parent levels 7–11 can show a puncture-shaped lobe on **covered** diagnostic cells; those are not added to the source on the active AMR hierarchy. [t11-uncovered-sources.csv](t11-uncovered-sources.csv) distinguishes the active cells, and the profile CSV labels covered cells. A seam may geometrically pass through the puncture, but the ghost-copy control and analytic source show that this does not create the leading source.

RHS self-convergence compares native low−mid and mid−high differences on the same physical P6 sample points, with p=log(D_low_mid/D_mid_high)/log(1.5). All supports within five low-grid cell spacings of a face are excluded; a P8 replay gives the interpolation spread. The first two windows are poorly qualified by interpolation and must not be used as asymptotic orders. Low/negative outer orders at tiny cancelling sources are retained.
''')
 data=[]
 for ray in ('axis','diagonal'):
  for x in orders:
   if x['ray']==ray and x['variable']=='Gamma1':data.append([ray,x['r_min_M']+'–'+x['r_max_M'],sci(x['low_mid_RMS']),sci(x['mid_high_RMS']),f"{float(x['order']):.6f}",f"{float(x['interpolation_fraction']):.6f}"])
 out.append(table(['ray','fixed r/M window','low−mid RHS RMS','mid−high RHS RMS','p','P6/P8 spread / mid−high'],data))
 out.append('''On the resolved 0.005–0.02 M window the Gamma RHS orders are 3.99285/4.00758 on axis/diagonal, with interpolation spread below 0.008 of the finer difference. The spatial source is therefore not globally nonconvergent at fixed r. Its conical first-cell error and increasingly singular curvature track the grid **as r~h→0**. The diagonal outer source is very small: at 0.05–0.1 M its high RMS is 4.463516e-4 but the native/continuum discrepancy is 1.680124e-9. The nonconvergent few-e-9 differences in the cancelling Gamma/K/Theta sources there are consistent with initial field/absolute-coordinate and derivative roundoff, as in T7-E; their precise share is inferred, not established by a new fused-coordinate intervention. They are vastly smaller than the inner O(10–100) source/error. [t11-rhs-orders.csv](t11-rhs-orders.csv) carries all variables and negative orders, and [t11-continuum-errors.csv](t11-continuum-errors.csv) carries the continuum and interpolation comparison.

![Native Gamma RHS and term profiles by level](figures/t11-term-profiles.png)

![RHS fixed-window self-convergence, including small-source floors](figures/t11-rhs-convergence.png)

### Continuum reconstruction and regularity

The t=0 continuum diagnostic loads E.trumpet through the actual C++ reader, retaining its inverse compactification, Clenshaw arithmetic and derivative coefficients. Analytic radial jets, followed by analytic Cartesian scalar/vector/tensor derivatives and cartoon contractions, are passed into the **same unmodified CCZ4Cartoon RHS equation**. There is no finite differencing of a radial table and no use of the constraint ODE to force a residual/source to zero. The initialized lapse remains X=sqrt(chi), not alpha_K. Native setter values and analytic values agree to less than 1e-9 on all sampled r>0.001 cells, including A, scalar momentum and electric fields; the readable check asserts this. The profile hash is 2a8de074ae17c0b11d323d4b0933a6bdb7430a5305473c8cc7d4ce37f39fa793.

For radial beta=b(r)n and the unboosted continuum state, the longitudinal Gamma RHS is

`(4/3)*(b''+2*b'/r−2*b/r²) −4*k*alpha' −6*alpha*k*chi'/chi −32*pi*alpha*Pi*phi'`.

Y0=0.13570329494190452, D0=2.2778090227037476 and k0=C*exp(D0)/Y0³=−3.559422033019796. Near the puncture alpha~r/Y0, chi~r²/Y0², beta is linear plus an r^(1+nu) correction, and A has an angular limit k0*(3n_i*n_j−delta_ij). Hence the Gamma-normal source tends to **−16*k0/Y0=419.6711108060989**, with its leading fractional correction scaling as r^(nu−1), nu=1.3372155112113842. The observed correction slope over 1e-8–1e-6 M is 0.337260, consistent with nu−1=0.337216. The source is bounded and locally volume-integrable, but its radial derivative diverges as r^(nu−2) and the vector source has no direction-independent puncture value. It is not a smooth compact continuum pulse.

The radial identities and leading balance have exact CAS witnesses and a separate 60-digit Cartesian differentiation cross-check; [CAS evidence](../../scripts/cas/t11-evidence.md) gives domain, exclusions and status. An independently reconstructed, forward-recurrence 60-digit Chebyshev source agrees with the Float64 C++-reader/CCZ4 result at six radii, with maximum normalized discrepancy 1.626175e-12. That numeric transfer check is CORROBORATED, not a continuum-physics proof. CAS scope: algebraic identity / transfer layer — production physics not certified. [t11-continuum-crosscheck.csv](t11-continuum-crosscheck.csv) retains the samples.

![Continuum source, native first-cell error and the scaled A12 KO source](figures/t11-continuum-puncture.png)

The principal answer is therefore **both**: the continuum initialized RHS itself has a puncture-end nonregular source, and the native puncture-end stencils add a nonuniform first-cell error and singular KO contributions. This is more specific than attributing everything to a smooth fractional-power pulse that merely needs better transport resolution.

### First-step validation and connection to launch

E-L's t0 state is the E-mid source state. Its saved step interval is Δt0=7/384=0.0182291667 M, whereas the finest native step is Δt12=4.45048e-6 M. Its native ray caches carry shift/Gamma/driver, lapse, chi, K, Theta and h; they do **not** carry A or the EMS fields. Those inputs are complete in the t=0 plots, so all native RHS variables are kernel-controlled, but no nonexistent positive-time A/EMS validation is claimed. Current-ray values and first derivatives are remapped to common points by cubic Hermite interpolation, then compared with P6/P8 sampling of the native t=0 RHS. This avoids a linear-profile interpolation error larger than some first-step changes.

For successive saved fields u0,u1,u2, the estimated quadratic term is 0.5*(u2−2*u1+u0). The table uses diagonal supports clean of faces on 0.1≤r≤0.49 M, with 95 common points. All norms are RMS in that fixed physical interval, without cylindrical/volume weighting. The final column subtracts the quadratic estimate and divides by the first change.
''')
 data=[]
 for name in ('B1','lapse','K','Gamma1','shift1','chi','h11','Theta'):
  x=next(x for x in val if x['ray']=='diagonal' and x['variable']==name and x['r_min_M']=='0.1')
  data.append([name,sci(x['dt_RHS_RMS']),sci(x['first_change_RMS']),sci(x['first_minus_dt_RHS_RMS']),sci(x['Odt2_from_steps12_RMS']),f"{float(x['residual_over_Odt2']):.6f}",f"{float(x['corrected_over_first']):.6f}"])
 out.append(table(['variable','Δt*RHS','first change','first minus Δt*RHS','quadratic estimate','remainder / quadratic','corrected / first'],data))
 out.append('''For B/lapse/K the remainder matches the quadratic estimate to approximately 0.1–0.5%, and the corrected fraction is 1.2e-6, 2.2e-5 and 1.4e-4. Gamma, beta, chi and h have material temporal corrections; their remainder/quadratic ratios are 1.18, 1.40, 0.95 and 0.95. The inner time scale and coupling are not assumed small. Theta is at the finite-difference/initial-rounding source floor: a snapshot RHS does not predict the tiny floating-point step change, and its 93.2 ratio is explicitly a failed Taylor/floor comparison. On the inner 0.0002–0.005 M window, Gamma's Δt*RHS RMS is 2.99558 but the first saved change is only 1.28157e-4; the estimated quadratic term is 1.41728e-4. This **does not validate** the inner Taylor expansion. The source has evolved through thousands of native stages by that snapshot; steps 1–2 cannot reconstruct its initial fine-step second derivative. [t11-step-validation.csv](t11-step-validation.csv) retains both rays/all fields/windows and interpolation/state-replay errors.

![Native-RHS versus saved-step remainder and quadratic estimate](figures/t11-step-validation.png)

Launch locations are checked in two ways. [t11-launch.csv](t11-launch.csv) retains the given direct native-cell provenance at steps 1–6. Independently, [t11-independent-launch.csv](t11-independent-launch.csv) selects the largest absolute Gamma curvature on the **same fixed 0.003–0.16 M interval at every time**, without using r=t or a predicted face arrival. Its clean-axis fitted speed is 0.975714 with intercept 0.000521661 M (four samples, fit RMS 0.000241374 M); the six clean diagonal samples give 0.984244 with intercept 0.000352845 M and fit RMS 0.000140814 M. P6/P8 position spread is at most 7.12043e-6 M on the clean axis and zero on the diagonal. These fits are compatible with the initial longitudinal shift characteristic speed approximately 1−beta_n, and backtrack to the sub-0.001 M source region. They do not determine the precise emission time or a continuum signal speed: dominant lobe signs change and finite-resolution dispersion remains.
''')
 out.append(table(['step','t/M','ray','source level','native curvature-lobe r/M','counted centres','Gamma at lobe','C_Gamma at lobe','clean'],[[x['coarse_step'],f"{float(x['time_M']):.9f}",x['ray'],x['level'],f"{float(x['dominant_curvature_peak_M']):.9f}",x['width_native_cells'] or 'unbracketed',sci(x['Gamma_at_curvature_peak']),sci(x['C_Gamma_at_peak']),x['clean']] for x in launch]))
 out.append('''At the first sample the lobe is inside L12, before R12 or its corner, and is mostly metric-consistent Gamma. This supports a puncture launch rather than birth at an inner AMR face. The t=0 source audit establishes that the nonregular source is already present; the cadence cannot show which earliest projection/KO/geometric update supplies the outgoing packet or exclude nonlinear reorganization before 0.018229 M.

### Branch decision and single next experiment

**Branch (ii), puncture-end discretization/regularity first, is supported.** The actual alpha~r cusp and direction-dependent A limit are essential, alongside the fractional correction r^(nu−1); it is not solely the error of differentiating an otherwise smooth source. The geometric Gamma RHS is already hundreds at t=0, its first-cell lapse-gradient error stays O(1), and coupled A KO grows as 1/h. The source and first outgoing packet lie well inside the finest face. Against a claim of complete causal closure, the data do not determine the relative emitted amplitude from these coupled sources.

Branch (i) is not implicated by a cancellation error: the native driver combination is exactly zero, and driver decay/advection/KO have their separately reported sizes. This does not prove that every possible B(0) intervention is ineffective, but the audit supplies no reason to prioritize B(0)=0. Branch (iii)'s **smooth initial RHS** premise is contradicted by the continuum and native puncture-end source, although early-stage instrumentation remains useful to measure its coupled evolution. A new transport enlargement or a 100 M receiving-resolution envelope is not justified by these results.

The next experiment should be **one very short source-refinement ladder with first-stage capture**, preserving the frozen gauge/equations/initialization, point transfers and sigma=1. Keep the existing levels through L12 and their physical faces; compare the baseline with nested L13 and L14 (faces 0.013671875 and 0.0068359375 M), so the new source spacings are h12/2 and h12/4. Use one common Δt0=0.005 M and stop at 0.005 M; the outgoing feature should remain inside every finest face. Compare the emitted current-field profiles on the fixed 0.001–0.004 M interval and at r=0.002 M; retain the innermost cells separately for mechanism attribution. Record the first native RK stages' Gamma split, A/alpha/EMS KO and projections. No positive-time static target or replenished static ghost data is allowed. Diagnostics belong in a Tests-only harness, with no production change.

**Pre-registered reading:** a Cauchy launch requires decreasing successive fixed-window Gamma/metric/shift profile differences (D_13_14/D_12_13≤0.8), differences exceeding five times P6/P8 interpolation spread, and packet amplitudes stabilizing (|A14/A13−1|≤0.1), while clean of all faces. A bounded continuum pulse need not decrease to zero. If the difference ratio is ≥1 while the feature remains ≤3 local centres or its amplitude grows, source refinement has not closed the mechanism; use the stage replay to decide whether the first material change comes from A/alpha KO/projection or the geometric Gamma source. Between 0.8 and 1, or with insufficient interpolation significance, the result is inconclusive. Track both bounded continuum source and h*KO(A): a pointwise 1/h KO peak alone is not a proof of nonconvergence in evolved observables. Do not change B(0), alpha(0), beta(0), gauge, equations or sigma in this test. No experiment was launched in T11.

### Reproduction and provenance

Run `t11-analyze.py build`, `native E-low`, `native E-mid`, `native E-high`, `continuum`, `measure`, `validate`, `ghost-check`, `check`, `figures`, `report` with `/Users/auroradysis/miniconda3/bin/python`, `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`, and `MPLCONFIGDIR=/private/tmp/ems-t11/mpl` for figures. The build reuses T7-E's standalone serial Chombo compiler/link helper; its source is entirely in Tests. Run the two CAS scripts listed in the evidence card after preparing the continuum cache. The readable checks assert native output identity, term reconstruction, exact cancellation, finite/floor controls, saved ghost identity and analytic-setter agreement. Each command completed in seconds, below the detached threshold; no long simulation or analysis process was started. Native replay peak RSS was 0.336 GB, one thread, with selected native field reads of 48.3 / 69.9 / 100.8 MB on low/mid/high. Dense caches are retained under /private/tmp/ems-t11 with hashes; public CSVs total approximately 4 MB. No files from other tranches were removed. [COMMIT-MANIFEST-T11.txt](COMMIT-MANIFEST-T11.txt) records the artifacts and frozen source/input provenance.
''')
 p=HERE/'README.md';old=p.read_text();prefix=old.split('\n## T11 —',1)[0].rstrip();p.write_text(prefix+'\n\n'+''.join(out))

if __name__=='__main__':
 if sys.argv[1]=='build':build()
 elif sys.argv[1]=='native':native(sys.argv[2])
 elif sys.argv[1]=='continuum':
  for leg in ('E-low','E-mid','E-high'):
   z=np.load(TMP/f'{leg}-native.npz');continuum(z['meta'][:,[4,5]],leg)
  r=np.geomspace(1e-8,.5,1800);continuum(np.column_stack((r/np.sqrt(2),r/np.sqrt(2))),'radial')
 elif sys.argv[1]=='measure':measurements()
 elif sys.argv[1]=='validate':validation()
 elif sys.argv[1]=='figures':figures()
 elif sys.argv[1]=='check':check()
 elif sys.argv[1]=='ghost-check':ghost_check()
 elif sys.argv[1]=='report':report()
 elif sys.argv[1]=='provenance':provenance()
