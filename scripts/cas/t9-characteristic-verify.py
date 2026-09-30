#!/usr/bin/env python3
"""Exact conditional shift-block witness; independent 60-digit check."""
import json, random
from pathlib import Path
import sympy as s
import mpmath as mp

# Claim: left eigenvectors of the frozen (beta_r,Gamma_n,B_driver) block.
# Ring QQ[b,c,mu], c>0, mu>0, b real, b != +/-c; g=3c^2/(4mu)>0.
# K/Theta/lapse/metric perturbations and lower-order/cartoon forcing held fixed.
# Witness l P - (b +/- c) l=0, polynomial numerators after rational cancellation.
b,c,mu=s.symbols('b c mu', real=True)
P=s.Matrix([[b,mu,-1],[c*c/mu,b,0],[0,0,0]])
# EXPLORE: nondegenerate special values select the small left-eigenvector route.
assert sorted(map(float,P.subs({b:s.Rational(1,10),c:1,mu:s.Rational(3,4)}).eigenvals()))==[-.9,0.,1.1]
for sign in (-1,1):
    lam=b+sign*c
    left=s.Matrix([[sign*c/mu,1,-sign*c/(mu*lam)]])
    assert all(s.fraction(s.cancel(v))[0]==0 for v in left*P-lam*left)
    assert s.cancel(left[2]+left[0]/lam)==0
mp.mp.dps=60; rng=random.Random(9009); worst=mp.mpf('0')
for _ in range(32):
    cc=mp.mpf(str(rng.uniform(.6,1.4))); mm=mp.mpf(str(rng.uniform(.5,1.))); bb=mp.mpf(str(rng.uniform(-.2,.2)))
    pp=mp.matrix([[bb,mm,-1],[cc*cc/mm,bb,0],[0,0,0]])
    for sign in (-1,1):
        ll=mp.matrix([[sign*cc/mm,1,-sign*cc/(mm*(bb+sign*cc))]])
        ans=ll*pp; target=(bb+sign*cc)*ll
        worst=max(worst,max(abs(ans[i]-target[i]) for i in range(3)))
assert worst<mp.mpf('1e-55')
x=s.Symbol('x');polynomial_checks=0
for n in (6,8):
    basis=[s.prod((x-j)/s.Integer(i-j) for j in range(n) if j!=i) for i in range(n)]
    for deg in range(n):
        residual=sum(basis[i]*s.Integer(i)**deg for i in range(n))-x**deg
        assert s.Poly(s.expand(residual),x).is_zero
        for derivative in range(4):
            assert s.Poly(s.expand(s.diff(residual,x,derivative)),x).is_zero
            polynomial_checks+=1
nx,ny,r=s.symbols('nx ny r',real=True)
for power in (1,2):
    vector_projection=nx*(nx*r)*((nx*r)**2+(ny*r)**2)**power+ny*(ny*r)*((nx*r)**2+(ny*r)**2)**power
    target=(nx*nx+ny*ny)**(power+1)*r**(2*power+1)
    for derivative in range(4):
        assert s.Poly(s.expand(s.diff(vector_projection-target,r,derivative)),r,nx,ny).is_zero
out=dict(status='PROVED',claim='conditional left eigenvectors of 3x3 frozen longitudinal shift block',
    witness='lP-lambda*l rational numerator identically zero, both signs',ring='QQ[b,c,mu]',
    domain='c>0,mu>0,b real,b!=+/-c; g=3c^2/(4mu)>0',branch='c is positive speed relative to shift',
    assumptions='frozen current metric and shift; K/Theta/lapse/metric principal perturbations excluded; cartoon radial terms, damping, KO and matter are forcing',
    outgoing='Gamma_n - c/mu*d_r beta_n - c/(mu*(c-beta_n))*B_driver_n',
    incoming='Gamma_n + c/mu*d_r beta_n - c/(mu*(c+beta_n))*B_driver_n',
    coordinate_speeds='outgoing c-beta_n; incoming -c-beta_n; B has zero principal speed',
    numeric_status='CORROBORATED',seed=9009,points=32,precision_digits=60,
    region='c in [0.6,1.4],mu in [0.5,1],b in [-0.2,0.2]',max_numeric_residual=str(worst),
    sampling_status='PROVED',sampling_witness='sum L_i(x) i^k - x^k polynomial zero in QQ[x], all 0<=k<n,n=6,8; derivatives 0..3',
    sampling_exact_checks=polynomial_checks,directional_exact_checks=8,
    directional_witness='radial cubic/quintic synthetic vectors project to (nx^2+ny^2)^(p+1) r^(2p+1), p=1,2; derivatives 0..3')
print(json.dumps(out,indent=2))
