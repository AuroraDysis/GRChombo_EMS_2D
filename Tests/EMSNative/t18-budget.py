#!/usr/bin/env python3
"""Frozen scalar-wave CFL proxy and explicitly assumed separation scaling."""
import json
import math
from pathlib import Path
import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent
y,q=s.symbols('y q',real=True)
z=s.symbols('z')
R=lambda x:1+x+x**2/2+x**3/6+x**4/24
imag=s.expand(R(s.I*y)*R(-s.I*y))
expected=1+y**6*(y**2-8)/576
assert s.cancel(imag-expected)==0
symbol=(q*q-8*q+7)/3
assert symbol.subs(q,-1)==s.Rational(16,3)
assert s.diff(symbol,q).subs(q,1)<0 # derivative increases on [-1,1]
cfl=s.sqrt(8)/s.sqrt(2*symbol.subs(q,-1))
assert s.cancel(cfl*cfl-s.Rational(3,4))==0
mp.mp.dps=80;errors=[]
for v in [mp.mpf(0),mp.mpf('.5'),mp.mpf(1),mp.sqrt(8),mp.mpf(3)]:
    # Independent direct complex RK polynomial versus its real exact witness.
    a=abs(sum((mp.j*v)**k/mp.factorial(k) for k in range(5)))**2
    b=1+v**6*(v*v-8)/576
    errors.append(abs(a-b))
assert max(errors)<mp.mpf('1e-70')
rows=[]
for d in (12,16,32):
    # An anchored planning assumption, not a new two-body trajectory prediction.
    infall=306*(d/32)**1.5
    end=infall+100*2.002174604
    for dt in (.25,.375,.5):
        steps=math.ceil(end/(1.75*dt))
        rows.append(dict(separation_Mi=d,assumed_infall_Mi=infall,
            post_common_clock_Mi=100*2.002174604,planned_end_Mi=end,
            dt_multiplier=dt,coarse_steps=steps,
            limitation='306*(d/32)^(3/2) is an anchored planning assumption. '
                'A d32 CTT companion cannot initialize d12/d16; new initial data and qualification required. '
                'Stop remains actual first-qualified-common plus 100 initial M_ADM.'))
record=dict(claim='RK4 imaginary-axis interval and the 2D fourth-order scalar-wave CFL proxy',
    algebra='Q[y,q], real y and -1<=q=cos(theta)<=1',
    witness='exact RK polynomial magnitude residual; monotone fourth-order Laplacian symbol; CFL squared=3/4',
    status='PROVED for the frozen scalar wave only',cfl_proxy=float(cfl),
    numerical_check='five fixed y samples; independent complex RK evaluation; mpmath 80 digits',
    worst_residual=str(max(errors)),
    exclusions='Does not include native upwind advection, KO, Cartoon axis sources, CCZ4 full symbol, '
        'AMR transfer or nonlinear stability. Reported CFL margins are proxies, not a safe-dt certificate.',
    clocks=rows)
(HERE/'t18-budget.json').write_text(json.dumps(record,indent=2)+'\n')
print('PROVED frozen scalar-wave proxy; safe production dt still requires native evidence')
