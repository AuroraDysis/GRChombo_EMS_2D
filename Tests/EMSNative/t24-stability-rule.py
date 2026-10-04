#!/usr/bin/env python3
"""RK4 Nyquist rule and sampled principal limits, with explicit scope.

S is weighted: cGamma*(H11+4*H22/3), U=H22+.75*H11.
The Nyquist rule has zero shift and diagonal H22>=H11. A shape/shift-dependent
whole-wave-number scan is also supplied; neither is an AMR stability proof.
"""
import csv
import importlib.util
import json
from pathlib import Path
import resource
import time
import mpmath as mp
import numpy as np
import sympy as sy

HERE=Path(__file__).resolve().parent


def R(z): return 1+z+z*z/2+z**3/6+z**4/24


def boundary(a):
    # Exact polynomial in Y=y², for real a and conjugate RK arguments.
    z,y=sy.symbols('z y',real=True)
    f=1+z+z*z/2+z**3/6+z**4/24
    p=sy.Poly(sy.expand(f.subs(z,a+sy.I*y)*f.subs(z,a-sy.I*y)-1),y)
    assert all(p.nth(i)==0 for i in (1,3,5,7))
    co=[float(p.nth(i)) for i in (8,6,4,2,0)]
    roots=np.roots(co)
    yy=min(float(q.real) for q in roots if abs(q.imag)<1e-9 and q.real>1e-10)
    mp.mp.dps=70
    def stages(z):
        k1=z;k2=z*(1+k1/2);k3=z*(1+k2/2);k4=z*(1+k3)
        return 1+(k1+2*k2+2*k3+k4)/6
    aa=mp.mpf(str(sy.numer(a)))/mp.mpf(str(sy.denom(a)));root=mp.findroot(lambda q:abs(stages(aa+1j*q))-1,mp.sqrt(yy))
    assert abs(root**2-yy)<mp.mpf('1e-10')
    for frac in ('0.9','1.1'):
        q=root*mp.mpf(frac)
        assert abs(stages(aa+1j*q)-(1+(aa+1j*q)+(aa+1j*q)**2/2+(aa+1j*q)**3/6+(aa+1j*q)**4/24))<mp.mpf('1e-65')
    return float(root),dict(real_z=str(a),polynomial_in_y=str(p.as_expr()),positive_root=str(root),
        identity_status='PROVED',root_status='CORROBORATED',precision_digits=70,
        independent_definition='four RK stages',root_residual=str(abs(stages(aa+1j*root))-1))


def principal_scan(p,courant,sigma,cGamma,grid=257):
    k=np.linspace(-np.pi,np.pi,grid);kx,ky=np.meshgrid(k,k,indexing='ij')
    H=np.linalg.inv([[p['h11'],p['h12']],[p['h12'],p['h22']]])
    U=H[1,1]+.75*H[0,0]
    modspec=importlib.util.spec_from_file_location('stability',HERE/'t24-gauge-stability.py')
    mod=importlib.util.module_from_spec(modspec);modspec.loader.exec_module(mod)
    sx,syy=mod.symbols(kx,p['beta'][0]),mod.symbols(ky,p['beta'][1])
    d=[sx[0],syy[0]];D=np.empty(kx.shape+(2,2),complex)
    D[...,0,0]=sx[1];D[...,1,1]=syy[1];D[...,0,1]=D[...,1,0]=d[0]*d[1]
    lap=sum(H[j,l]*D[...,j,l] for j in range(2) for l in range(2))
    QH=np.empty(kx.shape+(2,2),complex)
    for i in range(2):
        for j in range(2):QH[...,i,j]=(lap if i==j else 0)+sum(H[i,l]*D[...,l,j] for l in range(2))/3
    G=np.array(p['Gamma'])*p['dx'];QG=np.empty_like(QH)
    for i in range(2):
        for j in range(2):QG[...,i,j]=2*G[i]*d[j]/3-(sum(G[l]*d[l] for l in range(2)) if i==j else 0)
    adv=sx[3]+syy[3];ko=sigma*(sx[2]+syy[2]);a=adv+ko
    eta=p['eta']*p['dx']
    def maximum(S):
        Q=QH*(S/(U*cGamma/.75))+QG
        tr=(Q[...,0,0]+Q[...,1,1])/2
        detsqrt=np.sqrt(((Q[...,0,0]-Q[...,1,1])/2)**2+Q[...,0,1]*Q[...,1,0])
        e=np.stack([tr+detsqrt,tr-detsqrt],-1)
        rr=np.sqrt(eta*eta/4+cGamma*e)
        eig=np.concatenate([(a-eta/2)[...,None]+rr,(a-eta/2)[...,None]-rr],-1)
        amplitudes=np.where(eig.real<=1e-11,abs(R(courant*eig)),0);idx=np.unravel_index(np.argmax(amplitudes),amplitudes.shape)
        return float(amplitudes[idx]),[float(kx[idx[:2]]/np.pi),float(ky[idx[:2]]/np.pi)]
    lo,hi=0.,2000.
    assert maximum(lo)[0]<=1+1e-10
    for _ in range(42):
        mid=(lo+hi)/2
        if maximum(mid)[0]>1+1e-10:hi=mid
        else:lo=mid
    return lo,hi,maximum(hi)[1],float(max(abs(R(courant*a)).max(),abs(R(courant*(ko-.1*p['dx']))).max()))


def main():
    start=time.monotonic();p=json.loads((HERE/'t24-gauge-stability.json').read_text())['parameters']
    rows=[];witnesses=[]
    for nu in (.5,.375,.25,.125):
        for sigma in (0.,.3,.5,1.):
            y,w=boundary(-sy.Rational(str(2*sigma*nu)));witnesses.append(w)
            limit=3*y*y/(16*nu*nu)
            for cg in (.75,.375,1.):
                lo,hi,mode,other=principal_scan(p,nu,sigma,cg)
                rows.append(dict(dt_dx=nu,sigma=sigma,cGamma=cg,weighted_S_Nyquist_limit=limit,
                    metric_U_Nyquist_limit=limit*.75/cg,imaginary_RK_limit=y,
                    measured_shape_shift_sampled_S_lo=lo,measured_shape_shift_sampled_S_hi=hi,
                    kx_pi=mode[0],ky_pi=mode[1],advected_metric_chi_and_B_max_amplification=other,
                    scope='Nyquist zero shift / sampled decaying principal branches, measured shape and shift, zero Cartesian gradients'))
    growth=[]
    for S in (4.96,5.03):
        amp=float(abs(R(-1+1j*np.sqrt(4*S/3))))
        for level in (10,11,12):
            dx=1.75/2**level;dt=.5*dx
            growth.append(dict(S=S,level=level,dt_M=dt,scalar_gain_per_step=amp,efold_steps=1/np.log(amp),
                doubling_steps=np.log(2)/np.log(amp),efold_time_M=dt/np.log(amp),
                log10_gain_if_held_for_0p77_M=.77*np.log10(amp)/dt))
    # Quasi-frozen Nyquist ramp from actual S153/S154; maxima can migrate.
    slope=(4.962229183-4.818769342)/.875;limit=next(r['weighted_S_Nyquist_limit'] for r in rows if r['dt_dx']==.5 and r['sigma']==1 and r['cGamma']==.75)
    crossing=134.75+(limit-4.962229183)/slope;failure=135.173828125
    t=np.linspace(crossing,failure,20001);ss=4.962229183+(t-134.75)*slope
    logamp=np.log(abs(R(-1+1j*np.sqrt(4*ss/3))))
    ramp=[dict(level=l,log10_amplitude_gain=float(np.trapezoid(logamp,t)/(.5*1.75/2**l)/np.log(10))) for l in (10,11,12)]
    for name,data in [('t24-stability-limits.csv',rows),('t24-stability-growth.csv',growth)]:
        with (HERE/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
    receipt=dict(status='CORROBORATED sampled boundaries; exact polynomial identities PROVED',rows=rows,
        witnesses=witnesses,growth=growth,linear_ramp=dict(crossing_M=crossing,end_M=failure,end_S=float(ss[-1]),
        gains=ramp,scope='constant-coefficient quasi-frozen maximal-S ramp; not measured amplitude or unique-mechanism proof'),
        grid=257,elapsed_s=time.monotonic()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        table_scope='S limits exclude finite alpha, Cartesian background gradients, AMR/subcycling/boundaries; metric U limits require diagonal rho-fast branch')
    (HERE/'t24-stability-rule.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(seconds=receipt['elapsed_s'],RSS=receipt['peak_RSS_bytes'],ramp=receipt['linear_ramp'],limits=[r for r in rows if r['cGamma']==.75]),indent=2))
if __name__=='__main__':main()
