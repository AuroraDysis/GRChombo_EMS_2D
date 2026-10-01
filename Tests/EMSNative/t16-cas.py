#!/usr/bin/env python3
"""Exact scalar/matrix witnesses, then independent 70-digit inversions."""
import sys
sys.dont_write_bytecode=True
import json, random
from pathlib import Path
import sympy as s
import mpmath as mp

# Contract: real alpha,X>0, real b_i, |q|<1 (q=tanh(eta/2)), J>0,
# w>0. Equality over Q(alpha,X,bx,by,bz,q). Positive square-root branch.
# The 16 residuals g*g^{-1}-I prove the inverse; gtt then fixes the lapse.
a,X=s.symbols('alpha X',positive=True)
bx,by,bz,q=s.symbols('bx by bz q',real=True)
b=s.Matrix([bx,by,bz]); c=(1+q*q)/(1-q*q); sh=2*q/(1-q*q)
g=s.eye(4)/X**2; g[0,0]=(b.dot(b)/X**2-a*a)
for i in range(3): g[0,i+1]=g[i+1,0]=b[i]/X**2
gi=s.zeros(4); gi[0,0]=-1/a**2
for i in range(3):
    gi[0,i+1]=gi[i+1,0]=b[i]/a**2
    for j in range(3):gi[i+1,j+1]=(X**2 if i==j else 0)-b[i]*b[j]/a**2
L=s.eye(4); L[0,0]=L[1,1]=c; L[0,1]=L[1,0]=-sh
Li=L.copy();Li[0,1]=Li[1,0]=sh
inv=Li*gi*Li.T
def closes(v):
    num,den=s.fraction(s.cancel(s.together(v)));assert num==0,'NOT-CLOSED';return str(den)
assert all(closes(z) for z in g*gi-s.eye(4))
assert all(closes(z) for z in L*Li-s.eye(4))
J=(c-sh*bx)**2-sh**2*a**2*X**2
closes(inv[0,0]+J/a**2)
lo,hi=s.symbols('lo hi',positive=True)
closes(lo**2/(1+(lo/hi)**2-lo**2)-lo**2*hi**2/(lo**2+hi**2-lo**2*hi**2))
# Axis consistency witness for native cartoon terms containing D1(f)/y.
# For f=y^5, fourth-order D1 has error -4h^4; on y=h/2 this is -8h^3.
y,h=s.symbols('y h',positive=True)
fd=(-((y+2*h)**5)+8*(y+h)**5-8*(y-h)**5+(y-2*h)**5)/(12*h)
closes(fd-5*y**4+4*h**4)
closes(((fd-5*y**4)/y).subs(y,h/2)+8*h**3)
# No generated numerical kernel: independently assemble and invert the metric.
mp.mp.dps=70;rng=random.Random(1609);worst=mp.mpf(0);count=0
for _ in range(64):
    aa=mp.mpf(str(rng.uniform(.05,1)));xx=mp.mpf(str(rng.uniform(.02,1)))
    bb=[mp.mpf(str(rng.uniform(-.08,.08))) for i in range(3)]
    eta=mp.mpf(str(rng.uniform(-.3,.3)));cc=mp.cosh(eta);ss=mp.sinh(eta)
    gg=mp.eye(4)/xx**2;gg[0,0]=sum(z*z for z in bb)/xx**2-aa**2
    for i in range(3):gg[0,i+1]=gg[i+1,0]=bb[i]/xx**2
    ll=mp.eye(4);ll[0,0]=ll[1,1]=cc;ll[0,1]=ll[1,0]=-ss
    brute=(ll.T*gg*ll)**-1
    lapse_inverse=1/mp.sqrt(-brute[0,0])
    closed=aa/mp.sqrt((cc-ss*bb[0])**2-ss**2*aa**2*xx**2)
    worst=max(worst,abs(lapse_inverse/closed-1));count+=1
assert worst<mp.mpf('1e-60')
record=dict(status='PROVED',witness='g gi=I, L Li=I; transformed inverse gtt=-J/alpha^2; binary squared identity',
    assumptions='alpha,X>0; real b_i,q; |q|<1; J,w>0; binary 0<lo<=hi<=1',
    excluded_locus='alpha=0, X=0, q^2=1, J=0; binary hi=0 or denominator=0',
    numeric_status='CORROBORATED',seed=1609,count=count,decimal_precision=70,
    region='alpha [.05,1], X [.02,1], b_i [-.08,.08], eta [-.3,.3]',worst_relative=str(worst))
record['axis_witness']='PROVED: (D1_h(y^5)-5y^4)/y=-8h^3 at y=h/2; consistency example, not a dominance attribution'
Path(__file__).with_name('t16-cas.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
