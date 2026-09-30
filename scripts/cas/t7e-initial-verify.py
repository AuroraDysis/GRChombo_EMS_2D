"""Radial conformal Ricci identity; r>0, c(r)>0, 0<theta<pi; exact rational witness."""
import sympy as s
import json,random,mpmath as m
from pathlib import Path
r,t,z=s.symbols('r t z',positive=True,real=True);c=s.Function('c')(r);x=(r,t,z)
g=s.diag(1/c,r*r/c,r*r*s.sin(t)**2/c);gi=s.diag(c,c/r**2,c/(r*r*s.sin(t)**2))
G={}
for k in range(3):
 for i in range(3):
  for j in range(3):G[k,i,j]=sum(gi[k,l]*(s.diff(g[l,j],x[i])+s.diff(g[l,i],x[j])-s.diff(g[i,j],x[l]))/2 for l in range(3))
R=0
for i in range(3):
 for k in range(3):
  R+=gi[i,i]*(s.diff(G[k,i,i],x[k])-s.diff(G[k,i,k],x[i])+sum(G[k,i,i]*G[l,k,l]-G[l,i,k]*G[k,i,l] for l in range(3)))
claim=2*(s.diff(c,r,2)+2*s.diff(c,r)/r)-s.Rational(5,2)*s.diff(c,r)**2/c
w=s.cancel(s.together(s.expand_trig(R-claim)))
assert w==0,w
weights=[-s.Rational(1,12),s.Rational(4,3),-s.Rational(5,2),s.Rational(4,3),-s.Rational(1,12)]
assert sum(map(abs,weights))==s.Rational(16,3)
assert 2*sum(map(abs,weights))==s.Rational(32,3)
assert all(sum(a*s.Integer(k)**n for a,k in zip(weights,range(-2,3)))==(2 if n==2 else 0) for n in range(6))
m.mp.dps=60;rng=random.Random(945825)
# Independent 3D Cartesian Christoffel contraction from radial metric jets.
errors=[]
for _ in range(12):
 v=[m.mpf(str(rng.uniform(.2,2))) for i in range(3)];rr=m.sqrt(sum(a*a for a in v));n=[a/rr for a in v]
 cc=1+rr**2+m.exp(-rr)/10;cp=2*rr-m.exp(-rr)/10;cpp=2+m.exp(-rr)/10
 f=1/cc;fp=-cp/cc**2;fpp=2*cp**2/cc**3-cpp/cc**2
 d=[fp*a for a in n];dd=[[fpp*n[i]*n[j]+fp/rr*(int(i==j)-n[i]*n[j]) for j in range(3)] for i in range(3)]
 ga=lambda k,i,j:(int(k==j)*d[i]+int(k==i)*d[j]-int(i==j)*d[k])/(2*f)
 dg=lambda k,i,j,a:(int(k==j)*dd[i][a]+int(k==i)*dd[j][a]-int(i==j)*dd[k][a])/(2*f)-ga(k,i,j)*d[a]/f
 direct=sum((dg(k,i,i,k)-dg(k,i,k,i)+sum(ga(k,i,i)*ga(l,k,l)-ga(l,i,k)*ga(k,i,l) for l in range(3)))/f for i in range(3) for k in range(3))
 radial=2*(cpp+2*cp/rr)-m.mpf('2.5')*cp**2/cc;errors.append(abs(direct-radial))
assert max(errors)<m.mpf('1e-50')
out=dict(status='PROVED',contract='r>0, c(r)>0, 0<theta<pi; literal equality over rational differential expressions',witness=str(w),excluded='r=0,c=0,sin(theta)=0',seed=945825,n=12,dps=60,worst_numeric_residual=str(max(errors)),sample_region='[0.2,2]^3 Cartesian, chi=1+r^2+exp(-r)/10',stencil_absolute_weight_sum='16/3',ricci_coordinate_sensitivity_factor='32/3',scope='algebraic identity — production physics not certified')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
