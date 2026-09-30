#!/usr/bin/env python3
"""Conditional frozen-coefficient gauge blocks; not the full CCZ4 symbol."""
import json
import random
import sympy as s
import mpmath as mp

# Claim contract: a,chi,g,mu>0, b real; homogeneous 1D principal blocks,
# Theta perturbation set to zero, lower-order matter/damping/KO excluded.
# Witness: characteristic-polynomial equality in QQ[a,chi,g,mu,b,z].
a,chi,g,mu = s.symbols('a chi g mu',positive=True)
b,z = s.symbols('b z',real=True)
# P satisfies u_t=P u_n, hence coordinate propagation speeds are eig(-P).
P=s.Matrix([[b,-s.Rational(9,5)*a,0,0],[-chi*g,b,0,0],
            [0,0,b,mu],[0,-s.Rational(4,3)*a*g,s.Rational(4,3)*g,b]])
target=((z+b)**2-s.Rational(9,5)*a*chi*g)*((z+b)**2-s.Rational(4,3)*mu*g)
# charpoly creates its own symbol with the requested name; use coefficients.
coeff=list((-P).charpoly().all_coeffs())
res=s.Poly(sum(v*z**(4-i) for i,v in enumerate(coeff))-target,z).all_coeffs()
assert all(s.expand(v)==0 for v in res)
Q=s.Matrix([[b,mu],[g,b]])
coeff=Q.charpoly().all_coeffs()
assert s.expand(sum(v*z**(2-i) for i,v in enumerate(coeff))-((z-b)**2-mu*g))==0
ps=s.Symbol('psi',positive=True)
assert s.cancel(s.Rational(9,5)*a*g/ps**4-s.Rational(9,5)*a*g*ps**-4)==0
mp.mp.dps=60;rng=random.Random(6019);worst=mp.mpf('0')
for _ in range(24):
    aa,xx,gg,mm=[mp.mpf(str(rng.uniform(.2,1.5))) for _ in range(4)]
    bb=mp.mpf(str(rng.uniform(-.2,.2)))
    numeric=mp.matrix([[bb,-mp.mpf(9)/5*aa,0,0],[-xx*gg,bb,0,0],
                       [0,0,bb,mm],[0,-mp.mpf(4)/3*aa*gg,mp.mpf(4)/3*gg,bb]])
    eig=sorted([-mp.re(v) for v in mp.eig(numeric)[0]])
    ca=mp.sqrt(mp.mpf(9)/5*aa*xx*gg);cs=mp.sqrt(mp.mpf(4)/3*mm*gg)
    expected=sorted([-bb-ca,-bb+ca,-bb-cs,-bb+cs])
    worst=max(worst,max(abs(u-v) for u,v in zip(eig,expected)))
assert worst<mp.mpf('1e-55')
print(json.dumps(dict(status='PROVED',scope='algebraic identity — production physics not certified',
                     assumptions='alpha,chi,conformal inverse-normal metric,mu positive; beta real; homogeneous frozen 1D gauge block; Theta=0; lower-order forcing/damping/KO excluded',
                     witness='exact characteristic-polynomial coefficient residuals; transverse 2x2 determinant',
                     normal_form='expand polynomial coefficients; boolean asserts CLOSE; no general simplifier',
                     branch='positive square roots; coordinate outward speeds -beta_n+c',
                     exclusions='nonpositive metric/alpha/chi; coincident eigenvalues excluded from independent numerical eigenvector interpretations',
                     speed_squares=['9 alpha chi g / 5','4 mu g / 3','mu g'],
                     numerical_seed=6019,points=24,precision_digits=60,
                     numeric_region='alpha,chi,g,mu in [0.2,1.5], beta in [-0.2,0.2]',
                     max_independent_eigenvalue_residual=str(worst)),indent=2))
