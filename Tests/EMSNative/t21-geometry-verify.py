#!/usr/bin/env python3
"""Exact rational witnesses, followed by independent 50-digit four-metric inversion."""
import json, random, resource
from pathlib import Path
import sympy as s
import mpmath as mp

# Claim domain: real ADM data, alpha>0, chi>0, H positive definite, hww>0,
# r>0; principal positive area-density radius. No spherical symmetry assumed
# for the inverse-metric identity; its areal/Misner-Sharp interpretation needs it.
a,c,h11,h12,h22,hw,b1,b2,rt,rx,ry=s.symbols('alpha chi h11 h12 h22 hww beta1 beta2 Rdot Rx Ry',real=True)
H=s.Matrix([[h11,h12,0],[h12,h22,0],[0,0,hw]]);beta=s.Matrix([b1,b2,0]);gamma=H/c
g=s.zeros(4);g[0,0]=-a*a+(beta.T*gamma*beta)[0]
for i in range(3):
    g[0,i+1]=g[i+1,0]=(gamma*beta)[i]
    for j in range(3):g[i+1,j+1]=gamma[i,j]
inv=s.zeros(4);inv[0,0]=-1/a**2
for i in range(3):
    inv[0,i+1]=inv[i+1,0]=beta[i]/a**2
    for j in range(3):inv[i+1,j+1]=c*H.inv()[i,j]-beta[i]*beta[j]/a**2
matrix=g*inv-s.eye(4)
assert all(s.fraction(s.cancel(v))[0]==0 for v in matrix), 'inverse witness NOT-CLOSED'
d=s.Matrix([rt,rx,ry,0]);claim=c*(h22*rx**2-2*h12*rx*ry+h11*ry**2)/(h11*h22-h12**2)-(rt-b1*rx-b2*ry)**2/a**2
assert s.fraction(s.cancel((d.T*inv*d)[0]-claim))[0]==0, 'contraction NOT-CLOSED'
R,r,ht=s.symbols('R r ht',positive=True)
htdot,hwdot,cdot=s.symbols('htdot hwdot chidot',real=True)
# Differentiate the defining polynomial R^4 chi^2 = r^4 ht hww at fixed
# coordinate point and fixed angular basis; use its positive root to divide.
rate=R*(s.Rational(1,4)*(htdot/ht+hwdot/hw)-s.Rational(1,2)*cdot/c)
residual=4*R**3*rate*c**2+2*R**4*c*cdot-r**4*(htdot*hw+ht*hwdot)
num=s.fraction(s.cancel(residual.subs(r**4,R**4*c**2/(ht*hw))))[0]
assert num==0, 'area time-rate NOT-CLOSED'

mp.mp.dps=50;rng=random.Random(21009);worst=mp.mpf(0)
for _ in range(40):
    def v():return mp.mpf(str(rng.uniform(-1,1)))
    A=mp.matrix([[v() for j in range(3)] for i in range(3)])
    spatial=A.T*A+mp.eye(3);beta=mp.matrix([v(),v(),v()]);alpha=mp.mpf('0.05')+abs(v())
    metric=mp.matrix(4);metric[0,0]=-alpha**2+(beta.T*spatial*beta)[0]
    for i in range(3):
        metric[0,i+1]=metric[i+1,0]=(spatial*beta)[i]
        for j in range(3):metric[i+1,j+1]=spatial[i,j]
    derivative=mp.matrix([v(),v(),v(),v()]);gradient=derivative[1:4,0]
    brute=(derivative.T*(metric**-1)*derivative)[0]
    separate=(gradient.T*(spatial**-1)*gradient)[0]-(derivative[0]-(beta.T*gradient)[0])**2/alpha**2
    worst=max(worst,abs(brute-separate)/max(1,abs(brute),abs(separate)))
assert worst<mp.mpf('1e-44'), 'independent numerical disagreement'
receipt=dict(exact_status='PROVED',witnesses=['ADM inverse: every rational numerator zero','4-metric contraction numerator zero','positive area-density-radius time derivative numerator zero'],
    excluded_locus='alpha=0, chi=0, det(H)=0, r=0; positive SPD/area branch required',
    numerical_status='CORROBORATED',samples=40,seed=21009,precision_digits=50,worst_scaled_residual=str(worst),
    physical_scope='R_area_density is a local angular area-density proxy; areal radius and Misner-Sharp mass require spherical symmetry',
    peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
Path(__file__).with_name('t21-geometry-witness.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
