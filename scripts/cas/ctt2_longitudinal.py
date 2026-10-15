"""Exact witness: g=P^2 I, P>0, arbitrary real W,dW,dlogP.
In Q[P,P^-1,W,dW,dlogP], g L_g(W) g = P^2 flat_strain(W).
Deterministic expand/cancel normal forms; excluded loci P=0, r=0, r=-1.
"""
import sympy as s
P=s.symbols('P', positive=True)
W=s.Matrix(s.symbols('w:3')); u=s.Matrix(s.symbols('u:3'))
D=s.Matrix(3,3,s.symbols('d:9'))
C=lambda i,j,k: int(i==j)*u[k]+int(i==k)*u[j]-int(j==k)*u[i]
DW=D+s.Matrix(3,3,lambda i,j:sum(C(i,j,k)*W[k] for k in range(3)))
L=(DW+DW.T-s.Rational(2,3)*s.trace(DW)*s.eye(3))/P**2
strain=D+D.T-s.Rational(2,3)*s.trace(D)*s.eye(3)
assert all(s.expand(v)==0 for v in P**4*L-P**2*strain)
# W derivative formula witness in independent (r,mu) coordinates, r>0.
r,mu,y=s.symbols('r mu y', real=True); rr=s.symbols('rr', positive=True)
a=s.Function('a')(rr,mu); b=s.Function('b')(rr,mu)
rj,muj,delta=s.symbols('rj muj delta')
assert s.cancel(s.diff(rr*a/(1+rr),rr)*rj+s.diff(rr*a/(1+rr),mu)*muj-
    (rr/(1+rr)*(s.diff(a,rr)*rj+s.diff(a,mu)*muj)+a*rj/(1+rr)**2))==0
assert s.cancel(s.diff(y*b/(1+rr),rr)*rj+s.diff(y*b/(1+rr),mu)*muj+s.diff(y*b/(1+rr),y)*delta-
    (y/(1+rr)*(s.diff(b,rr)*rj+s.diff(b,mu)*muj)+b*(delta/(1+rr)-y*rj/(1+rr)**2)))==0
print('PROVED: 9 longitudinal and 2 vector derivative residuals exactly zero; P>0, r>0; no branch crossing or general simplifier')
print('CAS scope: algebraic identity — production physics not certified')
