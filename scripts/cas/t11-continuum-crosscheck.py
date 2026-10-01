#!/usr/bin/env python3
"""Independent 60-digit differentiation of the stored Chebyshev representation."""
import csv,json,hashlib
from pathlib import Path
import numpy as np
import mpmath as mp
mp.mp.dps=60
profile=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet');lines=profile.read_text().splitlines()
meta={a[2:].split('=',1)[0]:a.split('=',1)[1] for a in lines if a.startswith('# ')}
cs=[]
for field in 'YDS':
 i=lines.index(field+' 97');cs.append([mp.mpf(float(a)) for a in lines[i+1:i+98]])
def cheb(a,ss):
 x=2*ss-1;t0=mp.mpf(1);t1=x;z=a[0]+a[1]*x
 for c in a[2:]:t0,t1=t1,2*x*t1-t0;z+=c*t1
 return z
ell=mp.mpf(float(meta['ell']));nu=mp.mpf(float(meta['nu']));C=mp.mpf(float(meta['C']));pi=mp.pi
def radial(ss):return ell*ss**(1/nu)/(1-ss)
def alpha(ss):return ell*ss**(1/nu)/cheb(cs[0],ss)
def beta(ss):return -C*ell*ss**(1/nu)*(1-ss)**2/cheb(cs[0],ss)**3
def curv(ss):return C*mp.exp((1-ss)**2*cheb(cs[1],ss))*(1-ss)**3/cheb(cs[0],ss)**3
def phi(ss):return mp.mpf(float(meta['phi_inf']))+(1-ss)*cheb(cs[2],ss)
def Pi(ss):
 Y=cheb(cs[0],ss);G=Y+(1-ss)*mp.diff(lambda z:cheb(cs[0],z),ss)
 return -C*mp.exp((1-ss)**2*cheb(cs[1],ss))*(1-ss)**4*mp.diff(phi,ss)/(Y**2*G)
z=np.load('/private/tmp/ems-t11/radial-continuum.npz');r=np.hypot(*z['xy'].T);q=z['q'];rows=[]
for target in [1e-7,1e-5,1e-4,.001,.01,.1]:
 i=int(abs(r-target).argmin());R=mp.mpf(float(r[i]));ss=mp.findroot(lambda s:mp.log(radial(s)/R),((R/ell)**nu/2,(R/ell)**nu),solver='secant') if R<mp.mpf('.01') else mp.findroot(lambda s:radial(s)-R,(mp.mpf('.1'),mp.mpf('.8')))
 sr=1/mp.diff(radial,ss);srr=-mp.diff(radial,ss,2)*sr**3
 b=beta(ss);bp=mp.diff(beta,ss)*sr;bpp=mp.diff(beta,ss,2)*sr**2+mp.diff(beta,ss)*srr
 al=alpha(ss);ap=mp.diff(alpha,ss)*sr;k=curv(ss);phip=mp.diff(phi,ss)*sr
 g=mp.mpf(4)/3*(bpp+2*bp/R-2*b/R**2)-4*k*ap-12*k*ap-32*pi*al*Pi(ss)*phip
 actual=float(q[i,39:41].sum()/np.sqrt(2));error=abs(mp.mpf(actual)-g)
 rows.append(dict(r_M=float(R),native_continuum_Gamma_rhs=actual,independent_60digit_Gamma_rhs=str(g),absolute_error=float(error),relative_error=float(error/max(1,abs(g)))))
Y0=sum(c*(-1)**i for i,c in enumerate(cs[0]));D0=sum(c*(-1)**i for i,c in enumerate(cs[1]));k0=C*mp.exp(D0)/Y0**3;limit=-16*k0/Y0
assert Y0>0 and nu>1
assert max(x['relative_error'] for x in rows)<1e-9
with Path('Tests/EMSNative/t11-continuum-crosscheck.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
result=dict(status='CORROBORATED',scope='transfer layer — production physics not certified',points=6,precision_digits=60,profile_sha256=hashlib.sha256(profile.read_bytes()).hexdigest(),worst_relative_error=max(x['relative_error'] for x in rows),Y0=str(Y0),D0=str(D0),k0=str(k0),Gamma_rhs_limit=str(limit),nu=str(nu),definition='Independent Chebyshev recurrence + mpmath differentiation, compared to C++ reader jets + unmodified CCZ4Cartoon RHS; real R>0')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
