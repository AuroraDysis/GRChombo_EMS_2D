#!/usr/bin/env python3
"""Q-moment witness for the offline P8 probe, plus log-coordinate sensitivity."""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import sympy as s

moments=[]
for q,first in ((s.Rational(-1,4),-4),(s.Rational(1,4),-3)):
    weights=[s.prod((q-first-j)/s.Integer(k-j) for j in range(8) if j!=k)
             for k in range(8)]
    residuals=[s.expand(sum(weights[k]*(first+k)**p for k in range(8))-q**p)
               for p in range(8)]
    assert residuals==[0]*8
    moments.append(dict(q=str(q),first=first,weights=list(map(str,weights)),residuals=list(map(int,residuals))))
u,nu,x,delta=s.symbols('u nu x delta')
f=s.log(u)/nu-s.log(1-u)-x
log_coordinate_witness=s.expand(f.subs(x,x+delta)-f+delta)
assert log_coordinate_witness==0
ff=lambda xx:math.log(.37)/1.61-math.log1p(-.37)-xx
xx=math.log(.37)/1.61-math.log1p(-.37)
numeric_log_error=abs(ff(xx+1e-8)-ff(xx)+1e-8)
assert numeric_log_error<np.finfo(float).eps
here=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('t4b',here/'Tests/EMSNative/t4b-diagnose.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
yy,xx=np.mgrid[:64,:64];ix=np.array([20,31,42,53]);iy=np.array([21,32,43,54])
maximum=0.
for a in range(8):
    for b in range(8):
        coarse=((xx+.5)/64.)**a*((yy+.5)/64.)**b
        actual=module.prolong8(ix,iy,coarse,(0,0))
        expected=((ix+.5)/128.)**a*((iy+.5)/128.)**b
        maximum=max(maximum,float(np.max(abs(actual-expected))))
assert maximum<5e-15
print(json.dumps(dict(status='PROVED',scope='algebraic identities and offline probe — production physics and floor saturation not certified',
 assumptions='P8: eight distinct consecutive nodes, degree<=7 per coordinate; log inversion: 0<u<1, nu>0, x=log(R/ell), R,ell>0',
 witness='sum(w_k*(first+k)^p)-q^p=0 for p=0..7; f(x+delta)-f(x)+delta=0',
 simplifier_roles='deterministic expand CLOSE; no general simplifier',excluded_loci='coincident stencil nodes; u=0,1; nu=0',
 moments=moments,log_coordinate_residual=int(log_coordinate_witness),
 float64_tensor_polynomial_max_error=maximum,float64_log_crosscheck_error=numeric_log_error),indent=2))
