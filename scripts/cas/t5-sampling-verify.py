#!/usr/bin/env python3
"""Exact witnesses for T5 diagnostic sampling and speed conversion only."""
import json
import random
import sympy as s
import mpmath as mp

q = s.Symbol('q', real=True)
moments = {}
for n in (6, 8):
    w = [s.prod((q-j)/s.Integer(i-j) for j in range(n) if i != j) for i in range(n)]
    residuals = [s.expand(sum(w[i]*i**k for i in range(n))-q**k) for k in range(n)]
    assert residuals == [0]*n
    moments[n] = [str(v) for v in residuals]
# Equality of squares with nonnegative branches is equivalent to equality.
psi, alpha = s.symbols('psi alpha', positive=True)
chi = psi**-4
for c in (s.Rational(9, 5), s.Integer(2)):
    assert s.cancel(c*alpha/psi**4-c*alpha*chi) == 0
r, h, a = s.Rational(3, 2), s.Symbol('h', positive=True), s.Symbol('a', real=True)
assert s.cancel(a*h**4*(1-r**-4)-r**4*a*h**4*(r**-4-r**-8)) == 0
mp.mp.dps = 60
rng = random.Random(1729)
errors = []
for _ in range(32):
    ps = mp.mpf(str(rng.uniform(.5, 10)))
    al = mp.mpf(str(rng.uniform(.0001, 1)))
    errors.append(abs(mp.sqrt(mp.mpf('1.8')*al)/ps**2-mp.sqrt(mp.mpf('1.8')*al*ps**-4)))
assert max(errors) < mp.mpf('1e-55')
print(json.dumps(dict(status='PROVED', scope='algebraic identity — production physics not certified',
                     assumptions='real q; distinct integer stencil nodes; psi>0; alpha>=0 (zero by continuity); h>0; coefficient a real; ratio=3/2',
                     witnesses='Lagrange moment residuals in QQ[q]; equality of speed squares with positive branches; scaled fourth-order difference residual',
                     excluded_loci='psi=0; stencil denominators i-j=0 excluded by distinct nodes',
                     normal_form='expand polynomial coefficients and cancel rational residuals; asserts CLOSE; no general simplifier',
                     moments=moments, scale='81/16', numeric_seed=1729, numeric_points=32,
                     numeric_precision_digits=60, numeric_region='psi [0.5,10], alpha [0.0001,1]',
                     max_independent_speed_residual=str(max(errors))), indent=2))
