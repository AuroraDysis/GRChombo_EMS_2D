#!/usr/bin/env python3
"""Witness-first verification of the implicit-surface expansion geometry.

Claim: div(n/|n|) = (delta-ss):Hess(F)/|n| in Euclidean space,
F=x^2/c^2+(y^2+z^2)/a^2-1, a,c>0, n!=0, positive sqrt branch.
Witness: multiply the difference by |n|^3, a nonzero positive factor;
cancel its rational numerator. A generic curved metric is separately checked
against differentiating sqrt(det(g))*s^i directly at 50-digit precision.
This verifies geometry algebra, not the physical checkpoint or a horizon.
"""
import json,random
from pathlib import Path
import sympy as sp
import mpmath as mp
def main():
    x,y,z=sp.symbols('x y z',real=True)
    a,c=sp.symbols('a c',positive=True)
    q=sp.Matrix([x,y,z]);F=x*x/c**2+(y*y+z*z)/a**2-1
    n=sp.Matrix([sp.diff(F,v) for v in q]);N2=(n.T*n)[0];N=sp.sqrt(N2)
    direct=sum(sp.diff(n[i]/N,q[i]) for i in range(3))
    H=sp.hessian(F,q)
    projected=sum(((1 if i==j else 0)-n[i]*n[j]/N2)*H[i,j]/N for i in range(3) for j in range(3))
    # CLOSE: rational normal form after a branch-preserving positive factor.
    residual=sp.cancel(sp.together((direct-projected)*N**3))
    assert residual==0, 'NOT-CLOSED exact flat witness'
    mp.mp.dps=50;rng=random.Random(2501);errors=[]
    def metric(p):
        v=mp.matrix(p)
        return mp.exp(mp.mpf('.17')*p[0]+mp.mpf('.08')*sum(t*t for t in p))*(mp.eye(3)+mp.mpf('.2')*v*v.T)
    def normal(p,aa,cc):
        nn=mp.matrix([2*p[0]/cc**2,2*p[1]/aa**2,2*p[2]/aa**2]);gu=metric(p)**-1
        return gu*nn/mp.sqrt((nn.T*gu*nn)[0])
    def replace(p,k,v): return [v if j==k else p[j] for j in range(3)]
    for _ in range(16):
        p=[mp.mpf(str(rng.uniform(-.8,.8))) for j in range(3)]
        aa=mp.mpf(str(rng.uniform(.4,1.4)));cc=mp.mpf(str(rng.uniform(.4,1.4)))
        g=metric(p);gu=g**-1;nn=mp.matrix([2*p[0]/cc**2,2*p[1]/aa**2,2*p[2]/aa**2]);ss=normal(p,aa,cc)
        norm=mp.sqrt((nn.T*gu*nn)[0]);dg=[mp.diff(lambda v:metric(replace(p,k,v)),p[k]) for k in range(3)]
        lhs=sum(mp.diff(lambda v:mp.sqrt(mp.det(metric(replace(p,i,v))))*normal(replace(p,i,v),aa,cc)[i],p[i]) for i in range(3))/mp.sqrt(mp.det(g))
        rhs=mp.mpf(0)
        for i in range(3):
            for j in range(3):
                hess=(2/(cc**2 if i==0 else aa**2)) if i==j else mp.mpf(0)
                hess-=sum(gu[k,l]*(dg[i][l,j]+dg[j][l,i]-dg[l][i,j])*nn[k]/2 for k in range(3) for l in range(3))
                rhs+=(gu[i,j]-ss[i]*ss[j])*hess/norm
        errors.append(abs(lhs-rhs))
    assert max(errors)<mp.mpf('1e-45'), 'curved numeric falsifier'
    data=dict(exact_status='PROVED',witness=str(residual),domain='a,c>0; x,y,z real; grad F !=0; positive sqrt',
        curved_status='CORROBORATED',count=16,seed=2501,dps=50,region='x,y,z in [-.8,.8], a,c in [.4,1.4]',
        worst_residual=str(max(errors)),independent_definition='coordinate divergence of sqrt(det gamma)*s vs Hessian/Christoffel projector')
    path=Path(__file__).with_name('t25-geometry-verification.json');path.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data))
if __name__=='__main__':main()
