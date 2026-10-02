#!/usr/bin/env python3
"""Small exact witnesses and independent numerical eigenvalue checks."""
import sys
sys.dont_write_bytecode=True
import json, math
from pathlib import Path
import numpy as np
import sympy as s
import mpmath as mp
HERE=Path(__file__).resolve().parent
alpha,chi,H=s.symbols('alpha chi H',positive=True);z=s.symbols('z')
blocks=[(s.Matrix([[0,-s.Rational(9,5)*alpha],[-chi*H,0]]),s.Rational(9,5)*alpha*chi*H),
        (s.Matrix([[0,s.Rational(3,4)],[s.Rational(4,3)*H,0]]),H)]
witnesses=[]
for matrix,speed2 in blocks:
    residual=s.cancel((z*s.eye(2)-matrix).det()-(z*z-speed2));assert residual==0
    witnesses.append(dict(matrix=str(matrix),speed_squared=str(speed2),residual=str(residual),status='PROVED'))
worst=0.
for av,cv,hv in [(1.,1.,1.),(.2,.003,1.1),(.0002,6e-7,1.002),(.99,.999,1.0001)]:
    for matrix,speed2 in blocks:
        num=np.array(matrix.subs({alpha:av,chi:cv,H:hv})).astype(float)
        expected=float(s.sqrt(speed2).subs({alpha:av,chi:cv,H:hv}))
        error=np.max(abs(np.sort(np.linalg.eigvals(num))-np.array([-expected,expected])))/expected
        worst=max(worst,float(error));assert error<16*np.finfo(float).eps
mp.mp.dps=80
high_precision_worst=mp.mpf(0)
for av,cv,hv in [('1','1','1'),('.2','.003','1.1'),('.0002','.0000006','1.002'),('.99','.999','1.0001')]:
    aa,cc,hh=map(mp.mpf,(av,cv,hv))
    for matrix,expected in [(mp.matrix([[0,-mp.mpf('1.8')*aa],[-cc*hh,0]]),mp.sqrt(mp.mpf('1.8')*aa*cc*hh)),
                            (mp.matrix([[0,mp.mpf('.75')],[mp.mpf(4)*hh/3,0]]),mp.sqrt(hh))]:
        eigenvalues=sorted(mp.eig(matrix,left=False,right=False))
        error=max(abs(eigenvalues[0]+expected),abs(eigenvalues[1]-expected))/expected
        high_precision_worst=max(high_precision_worst,error)
        assert error<mp.mpf('1e-75')
M=s.Rational('2.002174604');dt=s.Rational(7,16);V=s.Integer(2);boundary=s.Integer(2240);r=s.Integer(100)
rows=[]
for n in (100,200,300):
    end=s.Integer(306)+n*M;steps=int(s.ceiling(end/dt));clock=steps*dt
    incoming=(boundary-r)/V;roundtrip=(2*boundary-16-r)/V
    assert incoming>clock and roundtrip>clock
    rows.append(dict(post_common_Mf=n,Mref_Mi=float(M),estimated_end_Mi=float(end),coarse_steps=steps,
        rounded_end_Mi=float(clock),incoming_first_Mi=float(incoming),incoming_margin_Mi=float(incoming-clock),
        light_roundtrip_Mi=float(2*boundary-16-r),gauge_roundtrip_Mi=float(roundtrip)))
record=dict(claim='Eigenvalues of the declared scalar lapse and longitudinal shift principal blocks and exact travel-time budgets',
    algebra='Q[alpha,chi,H,z], alpha>0 chi>0 H>0, positive square roots',
    equivalence='For each real 2x2 block the determinant polynomial gives eigenvalues; z^2-speed_squared has positive roots under the printed assumptions.',
    witnesses=witnesses,numerical_samples=4,numerical_precision='Float64; independent numpy.linalg.eigvals',
    worst_relative_residual=worst,
    independent_high_precision_check='mpmath eig, 80 decimal digits, eight matrices',
    high_precision_relative_residual=str(high_precision_worst),
    limitation='This proves the two declared blocks, not a full nonlinear CCZ4 characteristic theorem or a universal V<=2 bound. Actual native h/lapse/shift must be measured and the future run must monitor V.',
    domain_budgets=rows)
(HERE/'t17-budget.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
