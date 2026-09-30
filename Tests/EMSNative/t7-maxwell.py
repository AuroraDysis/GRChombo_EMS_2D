#!/usr/bin/env python3
"""Current checkpoint Maxwell source/flux budget; no static-profile evaluator."""
import csv, struct, importlib.util, hashlib, re, gzip, os, subprocess, time, sys
from pathlib import Path
import numpy as np
import h5py
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline

HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t7/maxwell')
spec=importlib.util.spec_from_file_location('t5',HERE/'t5-analyze.py')
t5=importlib.util.module_from_spec(spec);spec.loader.exec_module(t5)
VARS=t5.EVOLVED_COMPONENTS
ODD=[VARS.index(v) for v in ('h12','A12','Gamma2','shift2','B2','By','Bz','Ey','Ez')]
ODD_ALL=ODD+[28+i for i in ODD]+[56+i for i in ODD]+[84+i for i in ODD]
NORM=1/np.sqrt(2*np.pi)
EMS_ALPHA=4*np.pi
EMS_F2=-.8  # supplied E parameters, distinct from reference f2=-20

def save(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def read(path):
    chunks=[]
    with (path.open('rb') if path.exists() else gzip.open(path.with_suffix('.bin.gz'),'rb')) as f:
        assert f.read(8)==b'T7MAX001'
        time,h,centre=struct.unpack('<3d',f.read(24));level,n,parts=struct.unpack('<3I',f.read(12))
        assert n==28 and parts==88 and centre==224 and level==12
        while b:=f.read(16):
            x0,y0,x1,y1=struct.unpack('<4i',b);size=(x1-x0+1)*(y1-y0+1)
            data=np.frombuffer(f.read(size*118*8),dtype='<f8').reshape(y1-y0+1,x1-x0+1,118)
            chunks.append(((x0,y0,x1,y1),data))
    lo=np.min([b[:2] for b,a in chunks],axis=0);hi=np.max([b[2:] for b,a in chunks],axis=0)
    a=np.empty((118,hi[1]-lo[1]+1,hi[0]-lo[0]+1));seen=np.zeros(a.shape[1:],int)
    for (x0,y0,x1,y1),v in chunks:
        ys,xs=slice(y0-lo[1],y1-lo[1]+1),slice(x0-lo[0],x1-lo[0]+1)
        a[:,ys,xs]=v.transpose(2,0,1);seen[ys,xs]+=1
    assert np.all(seen==1) and np.isfinite(a).all() and lo[1]==0
    return dict(a=a,lo=lo,hi=hi,h=h,time=time,centre=centre,boxes=[b for b,v in chunks])

def sample(a,g,x,y,n=6,odd=ODD_ALL):
    shape=np.shape(x);x=np.ravel(x);y=np.ravel(y)
    xx=(x+g['centre'])/g['h']-.5-g['lo'][0];yy=y/g['h']-.5
    sx=np.clip(np.floor(xx).astype(int)-(n//2-1),0,a.shape[2]-n)
    sy=np.minimum(np.floor(yy).astype(int)-(n//2-1),a.shape[1]-n)
    ix=sx[:,None]+np.arange(n);raw=sy[:,None]+np.arange(n);iy=np.where(raw<0,-raw-1,raw)
    wx,wy=t5.weights(xx,sx,n),t5.weights(yy,sy,n)
    anchor=a[:,iy[:,n//2],ix[:,n//2]];out=anchor.copy()
    for j in range(n):
        for i in range(n):
            q=a[:,iy[:,j],ix[:,i]].copy();q[odd]*=np.where(raw[:,j]<0,-1,1)
            out+=(q-anchor)*(wy[:,j]*wx[:,i])[None,:]
    return out.reshape((a.shape[0],)+shape)

def geometry(u):
    chi,h11,h12,h22,hww=u[:5];det=h11*h22-h12*h12
    assert np.all(chi>0) and np.all(det>0) and np.all(hww>0)
    inv=np.array(((h22,-h12),(-h12,h11)))/det
    F=np.exp(-2*EMS_ALPHA*EMS_F2*u[18]**2)
    volume=np.sqrt(det*hww)/chi**1.5
    return inv,F,volume

def displacement(u,q=None):
    inv,F,volume=geometry(u);up=np.einsum('ij...,j...->i...',inv,u[24:26])
    D=volume*F*u[0]*up
    if q is None:return D
    qh=np.array(((q[1],q[2]),(q[2],q[3])))
    trace=np.einsum('ij...,ji...->...',inv,qh)
    factor=.5*trace+.5*q[4]/u[4]-.5*q[0]/u[0]-4*EMS_ALPHA*EMS_F2*u[18]*q[18]
    return volume*F*u[0]*(np.einsum('ij...,j...->i...',inv,q[24:26])-
        np.einsum('ij...,jk...,k...->i...',inv,qh,up)+up*factor)

def derivative(a,h,axis,odd_y=False):
    if axis==0:
        pad=np.concatenate((a[1::-1]*(-1 if odd_y else 1),a,a[-1:],a[-1:]),axis=0)
        return (pad[:-4]/12-2*pad[1:-3]/3+2*pad[3:-1]/3-pad[4:]/12)/h
    pad=np.pad(a,((0,0),(2,2)),mode='edge')
    return (pad[:,:-4]/12-2*pad[:,1:-3]/3+2*pad[:,3:-1]/3-pad[:,4:]/12)/h

def divergence(q,g):
    y=(np.arange(q.shape[1])+.5)*g['h']
    return derivative(q[0],g['h'],1)+derivative(q[1],g['h'],0,True)+q[1]/y[:,None]

def quadrature(g,shape,nt,nr):
    z,w=leggauss(nt);theta=(z+1)*np.pi/2;wt=w*np.pi/2
    rh=shape(theta);zr,wr=leggauss(nr)
    r=rh[:,None]+(.02-rh[:,None])*(zr+1)/2
    weight=2*np.pi*np.sin(theta[:,None])*r*r*wt[:,None]*wr*(.02-rh[:,None])/2
    return theta,rh,wt,r,weight

def flux(D,g,shape,nt,n=6):
    z,w=leggauss(nt);theta=(z+1)*np.pi/2;wt=w*np.pi/2
    r=shape(theta);dr=shape(theta,1);x,y=r*np.cos(theta),r*np.sin(theta)
    d=sample(D,g,x,y,n,odd=[1])
    radial=d[0]*np.cos(theta)+d[1]*np.sin(theta)
    tangent=-d[0]*np.sin(theta)+d[1]*np.cos(theta)
    return NORM*2*np.pi*np.sum(wt*np.sin(theta)*(r*r*radial-r*dr*tangent))

def shape_for(scale,checkpoint):
    p=Path('/private/tmp/ems-t6/horizons')/scale/Path(checkpoint).stem/'n96/shape-0-2.dat'
    s=np.fromstring(p.read_text(),sep=' ');assert len(s)==98 and abs(s[1]-224)<2e-13
    theta=(np.arange(96)+.5)*np.pi/96
    # Current numerical shape only; even continuation at both poles.
    x=np.r_[-theta[:4][::-1],theta,2*np.pi-theta[-4:][::-1]]
    y=np.r_[s[2:6][::-1],s[2:],s[-4:][::-1]]
    return CubicSpline(x,y),p

def selfcheck():
    h=.01;x=(np.arange(80)-40+.5)*h;y=(np.arange(40)+.5)*h
    X,Y=np.meshgrid(x,y);q=np.array((2*X,3*Y));g=dict(h=h)
    assert np.max(abs(divergence(q,g)[2:-2,2:-2]-8))<1e-12
    # Independent finite differences of the nonlinear displacement chain rule.
    rng=np.random.default_rng(7007);u=np.zeros((28,8,8));u[0]=.2+rng.random((8,8));u[1]=u[3]=u[4]=1;u[18]=.05
    u[24:26]=rng.normal(size=(2,8,8));q=rng.normal(size=u.shape)*.001
    exact=displacement(u,q)
    eps=1e-3;fd=(displacement(u+eps*q)-displacement(u-eps*q))/(2*eps)
    assert np.max(abs(fd-exact))/np.max(abs(exact))<1e-8

def analyze():
    finals=list(csv.DictReader((HERE/'t6-finder-final.csv').open()))
    spheres=list(csv.DictReader((HERE/'t6-fixed-spheres.csv').open()))
    budgets=[];sources=[];cleaners=[];layout=[];checks=[];boxes=[]
    constant=CubicSpline([0,np.pi],[.02,.02])
    for row in finals:
        scale,cp=row['scale'],row['checkpoint'];g=read(ROOT/scale/Path(cp).stem/'maxwell.bin')
        checkpoint=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0019')/scale/'chk'/cp
        parameters=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0019/submissions/exp-0019')/f'params-{scale}.txt'
        p=dict(re.findall(r'^\s*(\w+)\s*=\s*([^#\n]+)',parameters.read_text(),re.M))
        assert float(p['ems_alpha'])==EMS_ALPHA and float(p['ems_f2'])==EMS_F2 and float(p['ems_f0'])==float(p['ems_f1'])==0
        with h5py.File(checkpoint) as f:
            lev=f['level_12'];off=lev['data:offsets=0'][:]
            for i,b in enumerate(lev['boxes'][:]):
                x0,y0,x1,y1=(int(b[k]) for k in ('lo_i','lo_j','hi_i','hi_j'))
                u0=lev['data:datatype=0'][int(off[i]):int(off[i+1])].reshape(28,y1-y0+7,x1-x0+7)[:,3:-3,3:-3]
                native=g['a'][:28,y0-g['lo'][1]:y1-g['lo'][1]+1,x0-g['lo'][0]:x1-g['lo'][0]+1]
                assert u0.tobytes()==native.tobytes()
            for level in range(int(f.attrs['num_levels'])):
                z=f[f'level_{level}'];h=float(z.attrs['dx'])
                for i,b in enumerate(z['boxes'][:]):
                    x0,y0,x1,y1=(int(b[k]) for k in ('lo_i','lo_j','hi_i','hi_j'))
                    boxes.append(dict(scale=scale,checkpoint=cp,level=level,box=i,dx_M=h,xlo_M=x0*h-224,xhi_M=(x1+1)*h-224,ylo_M=y0*h,yhi_M=(y1+1)*h))
        checks.append(dict(scale=scale,checkpoint=cp,valid_state='BIT_IDENTICAL',GaussB_max=float(np.max(abs(g['a'][117]))),
                           current_input_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest()))
        u=g['a'][:28];phys=g['a'][28:56];ko=g['a'][56:84];clean=g['a'][84:112]
        shape,shape_path=shape_for(scale,cp);inv,F,vol=geometry(u);D=displacement(u)
        up=np.einsum('ij...,j...->i...',inv,u[24:26])*u[0]
        omitted=np.zeros_like(F)
        for axis in (0,1):
            dh=np.array(((derivative(u[1],g['h'],1-axis),derivative(u[2],g['h'],1-axis,True)),
                         (derivative(u[2],g['h'],1-axis,True),derivative(u[3],g['h'],1-axis))))
            trace=np.einsum('ij...,ji...->...',inv,dh)+derivative(u[4],g['h'],1-axis)/u[4]
            omitted-=.5*up[axis]*trace
        pieces={'native_total':displacement(u,phys-clean),'cleaning':displacement(u,clean),'KO_all':displacement(u,ko)}
        ne=np.zeros_like(phys);ne[24:26]=(phys-clean)[24:26]
        ns=np.zeros_like(phys);ns[18]=phys[18]
        pieces['native_E_only']=displacement(u,ne);pieces['native_scalar_only']=displacement(u,ns)
        pieces['native_metric_only']=pieces['native_total']-pieces['native_E_only']-pieces['native_scalar_only']
        qe=np.zeros_like(ko);qe[24:26]=ko[24:26]
        qs=np.zeros_like(ko);qs[18]=ko[18]
        pieces['KO_E_only']=displacement(u,qe);pieces['KO_scalar_only']=displacement(u,qs)
        pieces['KO_geometry_only']=pieces['KO_all']-pieces['KO_E_only']-pieces['KO_scalar_only']
        tendencies={name:divergence(q,g) for name,q in pieces.items()}
        source=np.array(list(tendencies.values()))
        for nt,nr in ((96,64),(192,128)):
            theta,rh,wt,r,weights=quadrature(g,shape,nt,nr)
            x,y=r*np.cos(theta[:,None]),r*np.sin(theta[:,None])
            # Frozen numerical Gauss diagnostic, already multiplied by F.
            density=vol*g['a'][116]
            sampled=sample(density[None],g,x,y,odd=[])[0]
            integral=NORM*np.sum(weights*sampled)
            qout=float(next(z['Q'] for z in spheres if z['scale']==scale and z['checkpoint']==cp and z['N_theta']=='96' and float(z['radius'])==.02))
            gap=qout-float(row['Q96'])
            fluxgap=flux(D,g,constant,nt)-flux(D,g,shape,nt)
            budgets.append(dict(scale=scale,checkpoint=cp,time_M=g['time'],N_theta=nt,N_r=nr,
                                finder_flux_gap=gap,signed_Gauss_volume=integral,volume_minus_finder_gap=integral-gap,
                                independent_density_flux_gap=fluxgap,volume_minus_density_flux_gap=integral-fluxgap,
                                shape_path=str(shape_path),dx_M=g['h']))
            if nt==192:
                data=sample(source,g,x,y,odd=[])
                for j,name in enumerate(pieces):
                    for region,selection in (('shell',np.ones(r.shape,bool)),('near_H',r<.008),('middle', (r>=.008)&(r<.012)),('outer',r>=.012)):
                        w=weights*selection;v=data[j]
                        sources.append(dict(scale=scale,checkpoint=cp,time_M=g['time'],term=name,region=region,
                                            signed_Gauss_density_rate=NORM*np.sum(w*v),L1_Gauss_density_rate=NORM*np.sum(w*abs(v)),
                                            coordinate_RMS_density_rate=np.sqrt(np.sum(w*v*v)/w.sum()),
                                            fixed_surface_flux_gap_rate=flux(pieces[name],g,constant,nt)-flux(pieces[name],g,shape,nt) if region=='shell' else np.nan))
                ce=g['a'][116]/F;delta=g['a'][112]-ce
                a=sample(np.array((ce,delta,u[13],u[27],F,vol,delta-omitted)),g,x,y,odd=[])
                cleaners.append(dict(scale=scale,checkpoint=cp,time_M=g['time'],RMS_CE=np.sqrt(np.sum(weights*a[0]**2)/weights.sum()),
                                     RMS_rhs_minus_diagnostic_CE=np.sqrt(np.sum(weights*a[1]**2)/weights.sum()),
                                     max_rhs_minus_diagnostic_CE=np.max(abs(a[1])),lapse_min=np.min(a[2]),lapse_max=np.max(a[2]),
                                     Xi_absmax=np.max(abs(a[3])),F_min=np.min(a[4]),F_max=np.max(a[4]),
                                     RMS_det_omission_replay_error=np.sqrt(np.sum(weights*a[6]**2)/weights.sum()),
                                     max_det_omission_replay_error=np.max(abs(a[6]))))
        layout.append(dict(scale=scale,checkpoint=cp,level=12,dx_M=g['h'],face_M=(g['hi'][1]+1)*g['h'],
                           horizon_min_M=float(shape((np.arange(96)+.5)*np.pi/96).min()),horizon_max_M=float(shape((np.arange(96)+.5)*np.pi/96).max()),
                           sphere_M=.02,coarse_fine_faces_in_shell=0,boxes=len(g['boxes'])))
        print(scale,cp,'gap',gap,'volume',integral,'KO/native/clean',[(z['term'],z['signed_Gauss_density_rate']) for z in sources[-24:] if z['region']=='shell'],flush=True)
    save('t7-maxwell-budget.csv',budgets);save('t7-maxwell-sources.csv',sources)
    save('t7-cleaner-audit.csv',cleaners);save('t7-E-shell-layout.csv',layout)
    save('t7-maxwell-checks.csv',checks);save('t7-E-boxes.csv',boxes)

def dump():
    """Recreate small frozen-state outputs from read-only checkpoints."""
    exe=Path('/private/tmp/ems-t7/maxwell.ex')
    if not exe.exists():exe=next((HERE.parent/'EMSRHFinder').glob('EMSRHCheckpoint2d.*.ex')).resolve()
    rows=[]
    for scale in ('E-low','E-mid'):
        inputs=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0019')/scale/'chk'
        template=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0019/submissions/exp-0019')/f'params-{scale}.txt'
        p=dict(re.findall(r'^\s*(\w+)\s*=\s*([^#\n]+)',template.read_text(),re.M))
        for k in ('output_path','hdf5_subpath','data_subpath','pout_subpath'):p.pop(k,None)
        for cp in sorted(inputs.glob('*.hdf5')):
            d=ROOT/scale/cp.stem;d.mkdir(parents=True,exist_ok=True)
            q={**p,'restart_file':str(cp),'ems_data_path':'/private/tmp/ems-t7/ABSENT-STATIC-INPUT','checkpoint_interval':-1,
               'plot_interval':-1,'verbosity':0,'t6_checkpoint_diagnostics':'true','t7_diagnostics':'true'}
            (d/'params.txt').write_text(''.join(f'{k} = {v}\n' for k,v in q.items()))
            start=time.monotonic()
            with (d/'run.log').open('w') as log:r=subprocess.run([str(exe),'params.txt'],cwd=d,env={**os.environ,'OMP_NUM_THREADS':'4'},stdout=log,stderr=subprocess.STDOUT,timeout=20)
            assert r.returncode==0
            rows.append(dict(scale=scale,checkpoint=cp.name,wall_seconds=time.monotonic()-start,exit=r.returncode,bytes=(d/'maxwell.bin').stat().st_size,
                             checkpoint_sha256=hashlib.sha256(cp.read_bytes()).hexdigest(),binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()))
    save('t7-maxwell-input-audit.csv',rows)

if __name__=='__main__':
    if '--dump' in sys.argv:dump()
    selfcheck();analyze()
