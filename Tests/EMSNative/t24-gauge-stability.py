#!/usr/bin/env python3
"""Frozen d507438 alpha=0 gauge/metric operator, discrete Fourier/RK4 audit.

The homogeneous background has zero Cartesian derivatives; its cylindrical
source terms are retained. This is a local constant-coefficient diagnostic,
not a variable-coefficient AMR stability proof. Gamma and B perturbations are
scaled by dx, and chi by its background value. No static fields are read.
"""
import argparse
import csv
import json
import resource
import time
from pathlib import Path

import h5py
import mpmath as mp
import numpy as np
import sympy as sy

FIELDS = 'beta1 beta2 B1 B2 Gamma1 Gamma2 h11 h12 h22 hww chi_relative'.split()


def checkpoint_cell(path):
    """Independent blockwise active-valid maximum of H22+.75 H11, L12."""
    with h5py.File(path, 'r') as f:
        names = [f.attrs[f'component_{i}'].decode() for i in range(int(f.attrs['num_components']))]
        l = int(f.attrs['num_levels'])-1
        g = f[f'level_{l}']; dx = float(g.attrs['dx'])
        offsets = g['data:offsets=0'][:]
        gx, gy = map(int, g['data_attributes'].attrs['outputGhost'])
        best = None
        for bi, box in enumerate(g['boxes'][:]):
            x0,y0,x1,y1 = [int(box[k]) for k in ('lo_i','lo_j','hi_i','hi_j')]
            nx,ny = x1-x0+1,y1-y0+1
            a = g['data:datatype=0'][int(offsets[bi]):int(offsets[bi+1])].reshape(len(names),ny+2*gy,nx+2*gx)[:,gy:gy+ny,gx:gx+nx]
            h11,h12,h22 = [a[names.index(k)] for k in ('h11','h12','h22')]
            det = h11*h22-h12*h12
            s = (h11+.75*h22)/det
            j,i = np.unravel_index(s.argmax(),s.shape)
            if best is None or float(s[j,i])>best['S']:
                best = dict(S=float(s[j,i]),level=l,cell=[x0+int(i),y0+int(j)],dx=dx,
                    rho=(y0+int(j)+.5)*dx,values={k:float(a[names.index(k),j,i]) for k in names})
    return best


def symbols(k, beta):
    d = 1j*(8*np.sin(k)-np.sin(2*k))/6
    d2 = (-30+32*np.cos(k)-2*np.cos(2*k))/12
    ko = -(np.sin(k/2)**6)
    pos = sum(w*np.exp(1j*k*j) for j,w in zip((-1,0,1,2,3),(-1/4,-5/6,3/2,-1/2,1/12)))
    neg = sum(w*np.exp(1j*k*j) for j,w in zip((-3,-2,-1,0,1),(-1/12,1/2,-3/2,5/6,1/4)))
    adv = beta*np.where(beta>=0,pos,neg)
    return d,d2,ko,adv


def operator(kx, ky, p, full=True):
    """All entries analytic in perturbations; no numerical differencing."""
    h = np.array([[p['h11'],p['h12']],[p['h12'],p['h22']]])
    H = np.linalg.inv(h); hw=p['hww']; Hw=1/hw
    beta=np.array(p['beta']); G=np.array(p['Gamma'])*p['dx']
    r=p['rho']/p['dx']; invr=1/r if full else 0.; invr2=invr*invr
    ds=[symbols(kx,beta[0]),symbols(ky,beta[1])]
    d=np.stack([s[0] for s in ds],-1)
    dd=np.empty(kx.shape+(2,2),complex)
    dd[...,0,0]=ds[0][1]; dd[...,1,1]=ds[1][1]
    dd[...,0,1]=dd[...,1,0]=d[...,0]*d[...,1]
    adv=ds[0][3]+ds[1][3]; ko=p['sigma']*(ds[0][2]+ds[1][2])
    shape=kx.shape; L=np.zeros(shape+(11,11),complex)
    for col in range(11):
        v=np.zeros(11); v[col]=1
        b=v[:2]; driver=v[2:4]; gam=v[4:6]
        dh=np.array([[v[6],v[7]],[v[7],v[8]]]); dhw=v[9]
        dH=-H@dh@H; dHw=-dhw/(hw*hw)
        div=np.sum(d*b,-1)+b[1]*invr
        divw=np.sum(d*b,-1)-2*b[1]*invr
        out=np.zeros(shape+(11,),complex)
        out[...,:2]=(adv+ko-p['eta']*p['dx'])[...,None]*b+p['cGamma']*gam-driver
        out[...,2:4]=(ko-.1*p['dx'])[...,None]*driver
        out[...,4:6]=(adv+ko+2*beta[1]*invr/3)[...,None]*gam
        for i in range(2):
            # CCZ4 kappa3=1: contracted Christoffel + 2Z/chi = evolved Gamma.
            out[...,4+i] += 2*G[i]*div/3-np.sum(G*d*b[i],-1)
            out[...,4+i] += Hw*(invr*d[...,1]*b[i]-invr2*(i==1)*b[1])
            out[...,4+i] += dHw*(-invr2*(i==1)*beta[1])
            out[...,4+i] += sum(H[j,k]*dd[...,j,k]*b[i] for j in range(2) for k in range(2))
            out[...,4+i] += sum(H[i,j]*dd[...,j,k]*b[k]/3 for j in range(2) for k in range(2))
            out[...,4+i] += sum(H[i,j]*(invr*d[...,j]*b[1]-invr2*(j==1)*b[1])/3 for j in range(2))
            out[...,4+i] += dH[i,1]*(-invr2*beta[1])/3
        for row,(i,j) in zip((6,7,8),((0,0),(0,1),(1,1))):
            out[...,row]=(adv+ko-2*beta[1]*invr/3)*dh[i,j]-2*h[i,j]*div/3
            out[...,row]+=sum(h[k,i]*d[...,j]*b[k]+h[k,j]*d[...,i]*b[k] for k in range(2))
        out[...,9]=(adv+ko+4*beta[1]*invr/3)*dhw-2*hw*divw/3
        out[...,10]=(adv+ko-2*beta[1]*invr/3)*v[10]-2*div/3
        L[...,col]=out
    return L


def scan(p,n,full):
    k=np.linspace(-np.pi,np.pi,n)
    kx,ky=np.meshgrid(k,k,indexing='ij')
    L=operator(kx,ky,p,full)
    eig=np.linalg.eigvals(L)
    z=p['courant']*eig
    amp=1+z+z*z/2+z*z*z/6+z*z*z*z/24
    absamp=np.abs(amp)
    # This homogeneous cylindrical background can have physical zero-order
    # growth. Keep that separate from unstable amplification of decaying modes.
    decay=eig.real<=1e-11
    flat=np.argmax(np.where(decay,absamp,-np.inf)); i,j,m=np.unravel_index(flat,absamp.shape)
    z0=z[i,j,m]; a0=absamp[i,j,m]
    allmax=np.unravel_index(np.argmax(absamp),absamp.shape)
    spectrum=np.linalg.eig(operator(np.array(kx[i,j]),np.array(ky[i,j]),p,full))
    idx=np.argmin(abs(spectrum[0]-eig[i,j,m]))
    vv=abs(spectrum[1][:,idx]); vv/=vv.max()
    nyq=np.linalg.eigvals(operator(np.array(np.pi),np.array(np.pi),p,full))
    zz=p['courant']*nyq; rr=1+zz+zz**2/2+zz**3/6+zz**4/24
    wave=np.argmax(np.where(abs(nyq.imag)>1e-4,abs(rr),-np.inf))
    return dict(max_amplification=float(absamp.max()),max_decaying_amplification=float(a0),
        nyquist_wave_amplification=float(abs(rr[wave])),nyquist_wave_z_real=float(zz[wave].real),nyquist_wave_z_imag=float(zz[wave].imag),
        kx_over_pi=float(kx[i,j]/np.pi),ky_over_pi=float(ky[i,j]/np.pi),
        z_real=float(z0.real),z_imag=float(z0.imag),
        all_max_kx_over_pi=float(kx[allmax[:2]]/np.pi),all_max_ky_over_pi=float(ky[allmax[:2]]/np.pi),
        positive_real_eigenvalues=int((~decay).sum()),
        largest_positive_real_lambda_dx=float(eig.real.max()),
        mode=';'.join(f'{name}:{x:.4g}' for name,x in zip(FIELDS,vv) if x>.01))


def exact_witness():
    y=sy.symbols('y',real=True); z=sy.symbols('z'); q=sy.symbols('q',nonzero=True)
    R=1+z+z*z/2+z**3/6+z**4/24
    d2=sum(sy.Rational(w)*q**j for j,w in zip((-2,-1,0,1,2),(' -1/12','4/3','-5/2','4/3','-1/12')))
    ko=sum(sy.Rational(w)*q**j for j,w in zip(range(-3,4),('1/64','-6/64','15/64','-20/64','15/64','-6/64','1/64')))
    assert d2.subs(q,-1)==-sy.Rational(16,3) and ko.subs(q,-1)==-1
    d1=(q**-2-8*q**-1+8*q-q**2)/12
    up=(-3*q**-1-10+18*q-6*q**2+q**3)/12
    assert d1.subs(q,-1)==0 and up.subs(q,-1)==-sy.Rational(8,3)
    assert sy.expand(ko-(q-1)**6/(64*q**3))==0
    aa,bb,cc,c,l=sy.symbols('H11 H12 H22 c lambda',real=True)
    HH=sy.Matrix([[aa,bb],[bb,cc]])
    QQ=-sy.Rational(16,3)*c*((aa+cc)*sy.eye(2)+HH/3)
    big=sy.zeros(4);big[:2,2:]=c*sy.eye(2);big[2:,:2]=QQ/c
    expected=(l*l+sy.Rational(16,3)*c*(aa+cc+aa/3))*(l*l+sy.Rational(16,3)*c*(aa+cc+cc/3))-(sy.Rational(16,9)*c*bb)**2
    characteristic=big.charpoly()
    assert sy.expand(characteristic.as_expr().subs(characteristic.gen,l)-expected)==0
    eq=sy.expand(R.subs(z,-1+sy.I*y)*R.subs(z,-1-sy.I*y)-1)
    poly=sy.Poly(eq,y)
    positive=[float(sy.re(r)) for r in sy.nroots(poly,maxsteps=200) if abs(float(sy.im(r)))<1e-12 and float(sy.re(r))>0]
    bound=max(positive); S=3*bound**2/4
    mp.mp.dps=70
    def stages(z):
        k1=z; k2=z*(1+k1/2); k3=z*(1+k2/2); k4=z*(1+k3)
        return 1+(k1+2*k2+2*k3+k4)/6
    root=mp.findroot(lambda y:abs(stages(-1+1j*y))-1,(mp.mpf('2.5'),mp.mpf('2.6')))
    assert abs(root-bound)<mp.mpf('2e-14')
    checks=[]
    for xx,yy in [('-1','2.5'),('-1','2.7'),('-.5','2.9'),('0','2.8')]:
        a=mp.mpc(xx,yy)
        assert abs(stages(a)-(1+a+a*a/2+a**3/6+a**4/24))<mp.mpf('1e-65')
        checks.append(dict(z=[xx,yy],stage_amplification=str(abs(stages(a)))))
    root1=mp.findroot(lambda y:abs(stages(mp.mpf('-.5')+1j*y))-1,3)
    return dict(status='exact stencil and RK4 identities verified; numerical matrix scan corroborated',
        rk4_minus_one_polynomial=str(eq),nyquist_d2=str(d2.subs(q,-1)),nyquist_KO_per_direction=-1,
        nyquist_d1=0,nyquist_upwind_positive=-8/3,
        principal_beta_Gamma_characteristic_polynomial=str(sy.expand(expected)),
        nyquist_scalar_limit_S=str(3*root*root/4),nyquist_imag_limit=str(root),
        one_direction_limit_H22=str(3*root1*root1/4),without_KO_limit_S=6,
        independent_70_digit_RK_stage_checks=checks)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--checkpoint',type=Path)
    ap.add_argument('--cell-json',type=Path)
    ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--grid',type=int,default=129)
    ap.add_argument('--h',type=float,nargs=4,metavar=('h11','h12','h22','hww'))
    ap.add_argument('--beta',type=float,nargs=2)
    ap.add_argument('--Gamma',type=float,nargs=2)
    ap.add_argument('--dx',type=float)
    ap.add_argument('--rho',type=float)
    ap.add_argument('--cGamma',type=float,default=.75)
    ap.add_argument('--eta',type=float,default=1)
    ap.add_argument('--sigma',type=float,default=1)
    ap.add_argument('--courant',type=float,default=.5)
    a=ap.parse_args(); start=time.monotonic(); a.output.mkdir(exist_ok=True,parents=True)
    cell=checkpoint_cell(a.checkpoint) if a.checkpoint else json.loads(a.cell_json.read_text()) if a.cell_json else None
    if cell:
        v=cell['values']; h=[v[k] for k in ('h11','h12','h22','hww')]
        beta=[v['shift1'],v['shift2']]; Gamma=[v['Gamma1'],v['Gamma2']]; dx=cell['dx']; rho=cell['rho']
        (a.output/'t24-gauge-cell.json').write_text(json.dumps(cell,indent=2)+'\n')
    else:
        h=a.h; beta=a.beta; Gamma=a.Gamma or [0,0]; dx=a.dx; rho=a.rho
        if h is None or beta is None or dx is None or rho is None: ap.error('provide checkpoint/cell-json or all h,beta,dx,rho')
    p=dict(zip(('h11','h12','h22','hww'),a.h or h),beta=a.beta or beta,Gamma=a.Gamma or Gamma,
        dx=a.dx or dx,rho=a.rho or rho,cGamma=a.cGamma,eta=a.eta,sigma=a.sigma,courant=a.courant)
    def S(p):
        H=np.linalg.inv([[p['h11'],p['h12']],[p['h12'],p['h22']]])
        return float(H[1,1]+.75*H[0,0])
    cases=[('step152',p)]
    for s in (4.9,4.96,5.03,5.5,6.4):
        b=p.copy(); factor=S(p)/s
        for name in ('h11','h12','h22','hww'): b[name]*=factor
        cases.append((f'S_{s}',b))
    for label,key,val in [('half_dt','courant',.25),('sigma_half','sigma',.5),('sigma_zero','sigma',0),('cGamma_half','cGamma',.375)]:
        cases.append((label,dict(p,**{key:val})))
    for s in (25.,26.):
        b=dict(p,courant=.25);factor=S(p)/s
        for name in ('h11','h12','h22','hww'):b[name]*=factor
        cases.append((f'half_dt_S_{s}',b))
    # Controller's isolated scalar criterion: no shift/advection/cartoon terms.
    cases.append(('no_shift_scalar',dict(p,beta=[0,0],Gamma=[0,0])))
    result=[]
    for label,b in cases:
        for full in (False,True):
            row=dict(case=label,operator='full_cartoon' if full else 'principal',S=S(b),metric_U=S(b),weighted_S=S(b)*b['cGamma']/.75,courant=b['courant'],sigma=b['sigma'],cGamma=b['cGamma'])
            row.update(scan(b,a.grid,full)); result.append(row)
    with (a.output/'t24-gauge-stability.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=result[0]); w.writeheader(); w.writerows(result)
    lo,hi=4.6,5.0
    for _ in range(22):
        mid=(lo+hi)/2;b=p.copy();factor=S(p)/mid
        for name in ('h11','h12','h22','hww'):b[name]*=factor
        if scan(b,129,False)['max_decaying_amplification']>1+1e-10:hi=mid
        else:lo=mid
    receipt=dict(parameters=p,scan_grid=a.grid,wavenumber_domain='[-pi,pi]^2',
        sampled_S_transition_interval=[lo,hi],transition_grid=129,
        exact_witness=exact_witness(),results=result,elapsed_s=time.monotonic()-start,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        limitation='zero Cartesian background gradients; cylindrical terms retained; no AMR/regrid or finite-alpha coupling proof')
    (a.output/'t24-gauge-stability.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(elapsed_s=receipt['elapsed_s'],peak_RSS_bytes=receipt['peak_RSS_bytes'],results=result),indent=2))


if __name__=='__main__': main()
