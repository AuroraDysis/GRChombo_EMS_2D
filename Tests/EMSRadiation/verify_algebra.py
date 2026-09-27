"""Exact rational stress-flux witnesses; no production code generation.

Domain: alpha>0, F>0, real E,B,beta,Pi,grad. Spatial orthonormal basis;
the identities are tensorial under spatial coordinate changes. Poles alpha=0
excluded. Witness: every component of direct -T^i_t minus the ADM formula,
cancelled over Q(alpha,F,E,B,beta,Pi,grad). NP witness uses outgoing B=s x E.
"""
import random
import sympy as s
import mpmath as mp

a,F=s.symbols('alpha F',positive=True)
beta=s.Matrix(s.symbols('b0:3',real=True))
E=s.Matrix(s.symbols('E0:3',real=True)); B=s.Matrix(s.symbols('B0:3',real=True))
p=s.Symbol('Pi',real=True); d=s.Matrix(s.symbols('d0:3',real=True))
g=s.eye(4); g[0,0]=-a*a+beta.dot(beta)
g[0,1:4]=beta.T; g[1:4,0]=beta
gi=s.zeros(4); gi[0,0]=-1/a**2
gi[0,1:4]=beta.T/a**2; gi[1:4,0]=beta/a**2
gi[1:4,1:4]=s.eye(3)-beta*beta.T/a**2
Fab=s.zeros(4)
for i in range(3):
    for j in range(3):
        Fab[i+1,j+1]=sum(s.LeviCivita(i,j,k)*B[k] for k in range(3))
for i in range(3):
    Fab[i+1,0]=a*E[i]+sum(Fab[i+1,j+1]*beta[j] for j in range(3))
    Fab[0,i+1]=-Fab[i+1,0]
T=2*F*(Fab*gi*Fab.T-g*sum((gi*Fab*gi)[i,j]*Fab[i,j]
                         for i in range(4) for j in range(4))/4)
rho=F*(E.dot(E)+B.dot(B)); S=2*F*E.cross(B)
stress=rho*s.eye(3)-2*F*(E*E.T+B*B.T)
J=a*S-beta*rho-stress*beta+beta*beta.dot(S)/a
direct=-(gi*T)[1:4,0]
assert all(s.cancel(x)==0 for x in direct-J), 'NOT-CLOSED EM flux'
dt=-a*p+beta.dot(d); grad4=s.Matrix([dt,*d])
Ts=2*grad4*grad4.T-g*(grad4.T*gi*grad4)[0]
assert all(s.cancel(x)==0 for x in -(gi*Ts)[1:4,0]+2*dt*(d-beta*p/a)), 'NOT-CLOSED scalar'
et,ep=s.symbols('Et Ep',real=True)
flat=Fab.subs(dict(zip([a,*beta,*E,*B],[1,0,0,0,0,et,ep,0,-ep,et])))
k=s.Matrix([1,-1,0,0])/s.sqrt(2); mb=s.Matrix([0,0,1,-s.I])/s.sqrt(2)
phi2=(mb.T*flat*k)[0]
assert s.expand(phi2*s.conjugate(phi2)-et**2-ep**2)==0, 'NOT-CLOSED Phi2'
u=s.Symbol('u',real=True); A,sigma=s.symbols('A sigma',positive=True)
assert s.integrate(2*A*A*u*u/sigma**4*s.exp(-u*u/sigma**2),(u,-s.oo,s.oo)) \
       == A*A*s.sqrt(s.pi)/sigma, 'NOT-CLOSED Gaussian pulse energy'

# Independent high-precision contractions at deterministic random points.
mp.mp.dps=30; rng=random.Random(1729); worst=mp.mpf(0)
for _ in range(24):
    av=mp.mpf(str(rng.uniform(.4,1.4))); fv=mp.mpf(str(rng.uniform(.2,3)))
    bv=mp.matrix([str(rng.uniform(-.3,.3)) for _ in range(3)])
    ev=mp.matrix([str(rng.uniform(-1,1)) for _ in range(3)])
    mv=mp.matrix([str(rng.uniform(-1,1)) for _ in range(3)])
    gv=mp.eye(4); gv[0,0]=-av*av+(bv.T*bv)[0]
    for i in range(3): gv[0,i+1]=gv[i+1,0]=bv[i]
    ff=mp.matrix(4)
    for i in range(3):
        for j in range(3):
            ff[i+1,j+1]=sum(int(s.LeviCivita(i,j,h))*mv[h] for h in range(3))
        ff[i+1,0]=av*ev[i]+sum(ff[i+1,j+1]*bv[j] for j in range(3))
        ff[0,i+1]=-ff[i+1,0]
    gu=gv**-1; fu=gu*ff*gu
    tt=2*fv*(ff*gu*ff.T-gv*sum(ff[i,j]*fu[i,j] for i in range(4) for j in range(4))/4)
    rr=fv*((ev.T*ev)[0]+(mv.T*mv)[0])
    sv=mp.matrix([2*fv*sum(int(s.LeviCivita(i,j,h))*ev[j]*mv[h]
                         for j in range(3) for h in range(3)) for i in range(3)])
    ss=rr*mp.eye(3)-2*fv*(ev*ev.T+mv*mv.T)
    jj=av*sv-bv*rr-ss*bv+bv*(bv.T*sv)[0]/av
    worst=max(worst,max(abs((gu*tt)[i+1,0]+jj[i]) for i in range(3)))
assert worst < mp.mpf('1e-14')
print('PROVED: scalar/EM rational flux components and outgoing Phi2 normalization')
print('Numerical transfer: 24 points, seed 1729, alpha [0.4,1.4], F [0.2,3], '
      'beta [-0.3,0.3], E/B [-1,1], 30 digits; residual < Float64 tolerance 1e-14')
print('CAS scope: algebraic identity; production extraction requires separate tests.')
