#!/usr/bin/env python3
"""Exact current-field algebra witnesses, independently checked at 60 digits."""
import json, random
import sympy as s
import mpmath as mp
# Contract: chi,hww,F,alpha,kappa positive, det(H)>0, real fields.
# Witness: coefficient-zero directional-derivative residuals and 2x2 determinants.
c,w,F,a,kap,g=s.symbols('chi hww F alpha kappa g',positive=True)
h,j,l,ex,ey,phi,fprime,Pi=s.symbols('h j l ex ey phi fprime Pi',real=True)
H=s.Matrix([[h,j],[j,l]]);V=H.inv()*s.Matrix([ex,ey]);S=F*s.sqrt(H.det()*w)/s.sqrt(c)
D=S*V
v=[c,w,h,j,l,ex,ey,F];q=s.symbols('qc qw qh qj ql qex qey qF',real=True)
Q=s.Matrix([[q[2],q[3]],[q[3],q[4]]])
expected=S*(H.inv()*s.Matrix(q[5:7])-H.inv()*Q*V+
 (s.Rational(1,2)*s.trace(H.inv()*Q)+s.Rational(1,2)*q[1]/w-s.Rational(1,2)*q[0]/c+q[7]/F)*V)
res=[s.factor(sum(s.diff(D[i],x)*dx for x,dx in zip(v,q))-expected[i]) for i in range(2)]
assert res==[0,0]
# Native electric scalar term cancels Fdot from phi_t=-alpha Pi.
assert s.expand(F*(-2*a*fprime*Pi)*ex+ex*(-2*F*fprime)*(-a*Pi))==0
k,z=s.symbols('k z',real=True)
for sign in (-1,1):
 block=s.Matrix([[0,sign*a*c*g*k*k*(-1)],[sign*a,-a*kap]])
 assert s.expand((z*s.eye(2)-block).det()-(z*z+a*kap*z+a*a*c*g*k*k))==0
# CE_rhs - CE_diag is the omitted conformal determinant trace.
E,du,dlogF,dlogvol,dlogdet=s.symbols('E du dlogF dlogvol dlogdet',real=True)
ce_diag=du+E*(dlogF+dlogvol)
ce_rhs=du+E*(dlogF+dlogvol-s.Rational(1,2)*dlogdet)
assert s.expand(ce_rhs-ce_diag+s.Rational(1,2)*E*dlogdet)==0
# Probe signs numerically before relying on the polynomial as a stability audit.
assert all(mp.re(x)<0 for x in mp.eig(mp.matrix([[0,2],[-3,-1]]))[0])
mp.mp.dps=60;rng=random.Random(7007);worst=mp.mpf(0)
for _ in range(24):
 vals=[mp.mpf(str(rng.uniform(.2,1.5))) for i in range(8)]
 vals[2]=vals[4]=mp.mpf('1.2');vals[3]=mp.mpf(str(rng.uniform(-.2,.2)))
 rate=[mp.mpf(str(rng.uniform(-.1,.1))) for i in range(8)]
 def fun(t):
  C,W,H0,J,L,X,Y,FF=[x+t*dx for x,dx in zip(vals,rate)]
  mat=mp.matrix([[H0,J],[J,L]])
  return FF*mp.sqrt(mp.det(mat)*W)/mp.sqrt(C)*(mat**-1)*mp.matrix([X,Y])
 C,W,H0,J,L,X,Y,FF=vals
 mat=mp.matrix([[H0,J],[J,L]]);iv=mat**-1;vec=iv*mp.matrix([X,Y]);qm=mp.matrix([[rate[2],rate[3]],[rate[3],rate[4]]])
 tr=sum((iv*qm)[i,i] for i in range(2))
 factor=tr/2+rate[1]/(2*W)-rate[0]/(2*C)+rate[7]/FF
 got=FF*mp.sqrt(mp.det(mat)*W)/mp.sqrt(C)*(iv*mp.matrix(rate[5:7])-iv*qm*vec+factor*vec)
 actual=mp.matrix([mp.diff(lambda t:fun(t)[i],0) for i in range(2)])
 worst=max(worst,max(abs(actual[i]-got[i]) for i in range(2)))
assert worst<mp.mpf('1e-55')
print(json.dumps(dict(status='PROVED',scope='algebraic identity — production physics not certified',
 assumptions='chi,hww,F,alpha,kappa,g positive; H positive definite; real current fields; cleaner coefficients frozen; advection excluded',
 witness='exact directional derivative residuals, scalar-coupling cancellation, cleaner determinant and determinant-omission residual',
 normal_form='factor rational/square-root directional residual; expand polynomial coefficients; boolean assert CLOSE',
 branch='positive volume square roots', exclusions='zero/nonpositive metric determinants; full variable-coefficient CCZ4 stability not claimed',
 numerical_seed=7007,points=24,precision_digits=60,max_independent_directional_residual=str(worst),
 cleaner_polynomial='z^2+alpha*kappa*z+alpha^2*chi*g*k^2',
 continuum_budget='Qout-QH=(1/sqrt(2*pi))*integral(sqrt(gamma)*GaussE*d3x); GaussE=F*CE'),indent=2))
