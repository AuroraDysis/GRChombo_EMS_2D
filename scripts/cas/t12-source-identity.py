#!/usr/bin/env python3
"""Exact rational witnesses for consult-7 identities (1)-(3), R>0."""
import json, random
from pathlib import Path
import sympy as sp
import mpmath as mp
R,r,rp,C,e,k,kp,a,ap,ak,akp,b,bp,bpp,c,m=sp.symbols('R r rp C e k kp a ap ak akp b bp bpp c m')
S=sp.Rational(4,3)*(bpp+2*bp/R-2*b/R**2)-4*k*ap-6*a*k*c-2*a*m
M=bp-b/R-3*ak*k;Mp=bpp-bp/R+b/R**2-3*(akp*k+ak*kp)
CM=2*kp+6*k/R-3*k*c-m
D=4*(kp*(ak-a)+k*(akp-ap))+12*k*(ak-a)/R
witnesses={}
def close(name,L,H):
 # CLOSE: unique rational normal form, excluded R=0 (and r=0 for map identity).
 num,den=sp.fraction(sp.cancel(sp.together(L-H)));assert num==0,(name,num)
 witnesses[name]={'numerator':str(num),'denominator':str(den)}
close('identity_3_off_constraint',S,D+sp.Rational(4,3)*(Mp+3*M/R)+2*a*CM)
constraints={bp:b/R+3*ak*k,bpp:3*(akp*k+ak*kp)+3*ak*k/R,m:2*kp+6*k/R-3*k*c}
close('identity_1_on_constraint',S.subs(constraints,simultaneous=True),D)
close('identity_2_matched_lapse',D.subs({a:ak,ap:akp}),sp.Integer(0))
br=-C/r**3+3*C*R*rp/r**4
close('map_shift_identity', (br+C/r**3).subs(rp,r*ak*e/R),3*ak*C*e/r**3)
close('metric_strain_with_Mbeta',sp.Rational(2,3)*(bp-b/R)-2*a*k,2*k*(ak-a)+sp.Rational(2,3)*M)
# Independent arbitrary-function falsifier of (3), not a shared residual evaluator.
mp.mp.dps=60;rng=random.Random(12012);worst=mp.mpf(0);points=[]
for i in range(24):
 x=mp.mpf(str(10**rng.uniform(-3,1)));points.append(str(x))
 bf=lambda z:z*z+mp.sin(z);kf=lambda z:mp.exp(-z)+z/7
 af=lambda z:1+z/3;akf=lambda z:z*z+mp.cos(z)/4
 cf=lambda z:2/z+z;mf=lambda z:mp.sin(2*z)
 mb=lambda z:mp.diff(bf,z)-bf(z)/z-3*akf(z)*kf(z)
 cm=2*mp.diff(kf,x)+6*kf(x)/x-3*kf(x)*cf(x)-mf(x)
 lhs=mp.mpf(4)/3*(mp.diff(bf,x,2)+2*mp.diff(bf,x)/x-2*bf(x)/x**2)-4*kf(x)*mp.diff(af,x)-6*af(x)*kf(x)*cf(x)-2*af(x)*mf(x)
 rhs=4*mp.diff(lambda z:kf(z)*(akf(z)-af(z)),x)+12*kf(x)*(akf(x)-af(x))/x+mp.mpf(4)/3*(mp.diff(mb,x)+3*mb(x)/x)+2*af(x)*cm
 worst=max(worst,abs(lhs-rhs)/max(1,abs(lhs),abs(rhs)))
assert worst<mp.mpf('1e-50')
result=dict(status='PROVED',scope='algebraic identity — production physics not certified',algebra='rational function field over QQ in R,r and independent jet symbols',domain='real R>0,r>0; e=exp(delta)>0; smooth radial fields for R>0; identity 1 assumes M_beta=M_beta_prime=C_M=0, identity 2 additionally alpha=alpha_K; identity 3 has no constraint assumption',normal_form='together/cancel numerator equality; no general simplifier or branch-crossing transformation',witnesses=witnesses,falsifier_seed=12012,precision_digits=60,falsifier_points=points,worst_normalized_falsifier=str(worst))
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
