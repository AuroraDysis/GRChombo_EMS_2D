#!/usr/bin/env python3
"""Exact radial-tensor witnesses and independent 60-digit numerical checks."""
import json, random
from pathlib import Path
import sympy as s
import mpmath as mp
r=s.Symbol('r',positive=True);nu=s.Symbol('nu',positive=True);q=s.Symbol('s',positive=True)
n=s.symbols('n0:3');b0,b1,b2,k0,k1=s.symbols('b0 b1 b2 k0 k1')
# Algebra Q(r,n,b0,b1,b2,k0,k1), modulo n.n=1; excludes r=0.
def D(f,j):
 return n[j]*(s.diff(f,r)+s.diff(f,b0)*b1+s.diff(f,b1)*b2+s.diff(f,k0)*k1)+sum(((1 if a==j else 0)-n[a]*n[j])*s.diff(f,n[a])/r for a in range(3))
ideal=sum(z*z for z in n)-1
def close(expr):
 num,den=s.fraction(s.cancel(expr));rem=s.rem(s.Poly(num,n[2]),s.Poly(ideal,n[2])).as_expr()
 assert s.cancel(rem)==0,(expr,rem);return str(den)
witness=[]
for i in range(3):
 lap=sum(D(D(b0*n[i],j),j) for j in range(3));gd=D(sum(D(b0*n[j],j) for j in range(3)),i)
 target=(b2+2*b1/r-2*b0/r**2)*n[i]
 witness.extend([close(lap-target),close(gd-target)])
 for j in range(3):
  for a in range(3):
   tij=k0*(3*n[i]*n[j]-(1 if i==j else 0))
   rhs=k1*n[a]*(3*n[i]*n[j]-(1 if i==j else 0))+3*k0/r*((1 if i==a else 0)*n[j]+(1 if j==a else 0)*n[i]-2*n[i]*n[j]*n[a])
   witness.append(close(D(tij,a)-rhs))
sr=nu*q*(1-q)/(r*(1+(nu-1)*q))
srr=s.diff(sr,q)*sr+s.diff(sr,r)
assert s.cancel(srr-(sr**2*(1/q-1/(1-q)-(nu-1)/(1+(nu-1)*q))-sr/r))==0
# Leading balance only: alpha=r/Y0, chi=r^2/Y0^2, A_nn=2*k0, beta=b0*r.
Y0=s.Symbol('Y0',positive=True)
leading=-2*(2*k0)*s.diff(r/Y0,r)-3*(r/Y0)*(2*k0)*s.diff(r*r/(Y0*Y0),r)/(r*r/(Y0*Y0))
assert s.cancel(leading+16*k0/Y0)==0
aa,bb=s.symbols('aa bb');trial=aa*r+bb*r**(1+nu)
shift=(s.diff(trial,r,2)+2*s.diff(trial,r)/r-2*trial/r**2)/r**(nu-1)
assert s.cancel(s.expand_power_base(shift)-bb*nu*(nu+3))==0
# Independent Cartesian derivatives versus radial formula; not generated code.
mp.mp.dps=60; rng=random.Random(11011);errors=[]
def beta(x,y,z,i):
 R=mp.sqrt(x*x+y*y+z*z);return (R+mp.mpf('.3')*R**mp.mpf('2.3372155112113842'))*(x,y,z)[i]/R
for _ in range(12):
 xyz=[mp.mpf('.2')+mp.mpf(str(rng.random())) for j in range(3)];R=mp.sqrt(sum(x*x for x in xyz));p=mp.mpf('2.3372155112113842')
 bv=R+mp.mpf('.3')*R**p;bp=1+mp.mpf('.3')*p*R**(p-1);bpp=mp.mpf('.3')*p*(p-1)*R**(p-2)
 for i in range(3):
  lap=mp.mpf('0')
  for j in range(3):
   def f(t):
    xx=xyz.copy();xx[j]=t;return beta(*xx,i)
   lap+=mp.diff(f,xyz[j],2)
  errors.append(abs(lap-(bpp+2*bp/R-2*bv/R**2)*xyz[i]/R))
assert max(errors)<mp.mpf('1e-50')
record=dict(status='PROVED',scope='algebraic identity — production physics not certified',
 claim='Radial vector Laplacian/grad-div, radial A derivatives, inverse compactification second derivative, leading lapse/chi balance',
 domain='r>0, nu>0, 0<s<1, Y0>0, n.n=1; real positive powers; excludes r=0 and s=0,1',
 witness='Exact rational numerators reduced modulo n.n-1; coefficient cancellation for srr and leading balance',
 exact_components=len(witness)+3,simplifier='cancel/rem only, assertions consumed',
 numerical_points=12,numerical_components=36,precision_digits=60,seed=11011,worst_residual=str(max(errors)))
dest=Path(__file__).with_suffix('.json');dest.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
