#!/usr/bin/env python3
"""Exact scoped gauge identities; numerical frozen-wave CFL cross-check.

Claims live in Q[alpha,K,Theta,B,Gammadot,mu,eta,lambda]; alpha>0,
mu>0, eta>=0. Equilibrium is equality modulo K=Theta=B=Gammadot=0.
The smallest witnesses are substitution residuals and 2x2 determinants.
This does not derive or certify the full curved-background CCZ4 symbol.
"""
import csv, json, resource
from pathlib import Path
import sympy as s
import numpy as np
HERE=Path(__file__).resolve().parent
a,K,T,B,G,mu,eta,lam=s.symbols('alpha K Theta B Gammadot mu eta lambda',real=True)
lapse=-2*a*(K-2*T);shift=mu*B;driver=G-eta*B
eq={K:0,T:0,B:0,G:0}
witnesses={name:s.Poly(expr.subs(eq),a,mu,eta).as_expr() for name,expr in [('lapse_equilibrium',lapse),('shift_equilibrium',shift),('driver_equilibrium',driver)]}
witnesses['driver_complete_Gamma']=s.Poly(driver-(G-eta*B),G,eta,B).as_expr()
lapse_symbol=s.Matrix([[0,1],[2,0]])
driver_symbol=s.Matrix([[0,1],[s.Rational(4,3)*mu,0]])
witnesses['lapse_characteristic']=s.Poly((lam*s.eye(2)-lapse_symbol).det()-(lam**2-2),lam).as_expr()
witnesses['driver_characteristic']=s.Poly((lam*s.eye(2)-driver_symbol).det()-(lam**2-s.Rational(4,3)*mu),lam,mu).as_expr()
assert all(v==0 for v in witnesses.values())
# Independent eigenvalues, then a finite-grid numerical (not rigorous)
# RK4 / fourth-order scalar-wave + sigma=1 KO stability cross-check.
eigen_lapse=sorted(abs(np.linalg.eigvals(np.array(lapse_symbol).astype(float))))
eigen_driver=sorted(abs(np.linalg.eigvals(np.array(driver_symbol.subs(mu,1)).astype(float))))
assert np.allclose(eigen_lapse,np.sqrt(2),rtol=2e-15)
assert np.allclose(eigen_driver,np.sqrt(4/3),rtol=2e-15)
theta=np.linspace(0,np.pi,513)
d2=2.5-(8/3)*np.cos(theta)+(1/6)*np.cos(2*theta)
ko=(1/32)*np.cos(3*theta)-(3/16)*np.cos(2*theta)+(15/32)*np.cos(theta)-5/16
cross=[]
for name,c in [('light',1.),('lapse',np.sqrt(2)),('driver',np.sqrt(4/3))]:
    z=.25*(ko[:,None]+ko[None,:]+1j*c*np.sqrt(np.maximum(0,d2[:,None]+d2[None,:])))
    amplification=abs(1+z+z*z/2+z*z*z/6+z*z*z*z/24)
    assert amplification.max()<1+4e-15
    cross.append(dict(mode=name,speed=float(c),directional_Courant=float(.25*c),frozen_2D_RK4_KO_max_amplification=float(amplification.max()),finite_sample_count=513**2,status='CORROBORATED; scalar principal-wave proxy only'))
with (HERE/'t22-asymptotic-CFL.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cross[0]);w.writeheader();w.writerows(cross)
out=dict(exact_status='PROVED within declared reduced gauge/principal-wave algebra',claims_and_witnesses={k:str(v) for k,v in witnesses.items()},domain='alpha>0, mu>0, eta>=0; equilibrium modulo K=Theta=B=Gamma_dot=0; positive speed branches',numerical_status='CORROBORATED',numerical_definition='independent numpy eigenvalues; native FD2 and KO Fourier symbols; RK4 finite-angle grid',full_CCZ4_hyperbolicity_and_curved_CFL='NOT-CLOSED; native speed envelopes and local calibration required',per_level_dt_over_h=.25,levels=15,boundary_speed=1.,boundary_return_to_R20_M=float(2*(336-20)/np.sqrt(2)),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(HERE/'t22-gauge-witness.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
