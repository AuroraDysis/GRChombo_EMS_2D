#!/usr/bin/env python3
"""Exact Q-moments and Q[F0,F1,F2,F3,theta][H]/(H^4) stage witness."""
import json
from pathlib import Path
import sympy as s

def weights(q, first):
    return [s.prod((q-first-j)/s.Integer(k-j) for j in range(6) if j != k)
            for k in range(6)]
moments = []
for q,first,name in [(s.Rational(-1,4),-3,'prolong_left'),
                     (s.Rational(1,4),-2,'prolong_right'),
                     (s.Rational(1,2),-2,'restriction')]:
    w=weights(q,first)
    residuals=[s.expand(sum(w[k]*(first+k)**p for k in range(6))-q**p)
               for p in range(6)]
    assert residuals == [0]*6, (name,residuals)
    moments.append(dict(name=name,weights=[str(x) for x in w],residuals=[int(x) for x in residuals]))
# Generic vector ODE jet: B=F'F, C=F''(F,F), D=F'^2F are
# independent formal elementary differentials, so this does not assume that
# Jacobians commute or collapse vector relations to a scalar example.
H,theta,f,B,C,D=s.symbols('H theta f B C D')
def truncate(p):
    expanded=s.expand(p)
    return sum(expanded.coeff(H,j)*H**j for j in range(4))
k=[f,f+H*B/2+H**2*C/8,
   f+H*B/2+H**2*C/8+H**2*D/4,
   f+H*B+H**2*(C+D)/2]
a=[0,H*k[0],H*(-3*k[0]/2+k[1]+k[2]-k[3]/2),
   H*s.Rational(2,3)*(k[0]-k[1]-k[2]+k[3])]
diff=H*(k[2]-k[1]);r=s.Rational(1,2)
c=[(theta,theta**2,theta**3,0),
   (theta+r/2,theta*(r+theta),theta**2*(3*r/2+theta),0),
   (theta+r/2,r*r/2+theta*(r+theta),3*r**3/8+theta*(3*r*r/2+theta*(3*r/2+theta)),-r**3/4),
   (theta+r,r*r+theta*(2*r+theta),3*r**3/4+theta*(3*r*r+theta*(3*r+theta)),r**3/2)]
start=theta*H*f+theta**2*H**2*B/2+theta**3*H**3*(C+D)/6
fstart=f+theta*H*B+theta**2*H**2*(C+D)/2
Bstart=B+theta*H*(C+D)
fine=[start,start+r*H*fstart/2,
      start+r*H*fstart/2+r**2*H**2*Bstart/4+r**3*H**3*C/16,
      start+r*H*fstart+r**2*H**2*Bstart/2+r**3*H**3*(C+2*D)/8]
residuals=[]
for stage,(b1,b2,b3,d) in enumerate(c):
    residual=truncate(a[1]*b1+a[2]*b2+a[3]*b3+diff*d-fine[stage])
    coeffs=[s.expand(residual).coeff(H,j) for j in range(4)]
    assert all(x==0 for x in coeffs),(stage,coeffs)
    residuals.append([str(x) for x in coeffs])
# The original Chombo coefficients have exact nonzero counterexamples.
legacy_residual=[s.expand(truncate(diff*(-r*r/4+r**3/4))),
                 s.expand(truncate(diff*(r*r/2-r**3/2)))]
assert legacy_residual == [-D*H**3/128,D*H**3/64]
record=dict(status='PROVED',scope='algebraic identity — production physics not certified',
            assumptions='ratio=2; point values; six distinct consecutive real nodes; degree<=5 per coordinate; autonomous smooth vector ODE formal jet; theta real; equality modulo H^4',
            excluded_loci='none; stencil denominators are nonzero integer differences',
            simplifier_roles='no general simplifier; deterministic expand and coefficient extraction CLOSE; asserts exit nonzero',
            branches='none',moments=moments,stage_coefficients=residuals,legacy_stage_residuals=[str(x) for x in legacy_residual])
print(json.dumps(record,indent=2))
