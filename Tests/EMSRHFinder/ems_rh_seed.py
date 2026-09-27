"""Known stationary horizon seed; algebra only, not a finder validation.

At t=0 the supplied Lorentz map is X=gamma*x, Y=y. A rest-frame
stationary horizon r_rest=R therefore satisfies gamma^2*x^2+y^2=R^2.
With x=f*cos(theta), y=f*sin(theta), choose the positive radial root.
Domain: R>0, gamma>=1, real theta; cos^2+sin^2=1 implies D>=1.
"""
import math
import sympy as s

R,g,c,z=s.symbols('R g c z',real=True)
D=g*g*c*c+z*z
f_squared=R*R/D
# CLOSE in Q(R,g,c,z); excluded D=0. No general simplifier is used.
assert s.cancel(g*g*f_squared*c*c+f_squared*z*z-R*R)==0
# D = 1+(g^2-1)c^2 on c^2+z^2=1, hence D>=1 for gamma>=1.
assert s.Poly(s.expand((D-1-(g*g-1)*c*c)-(c*c+z*z-1)),g,c,z).is_zero
# Positive R and D select f=R/sqrt(D), with no sign ambiguity.
worst=0.;count=0
for radius in (0.63593977642346233,0.5701934158668871):
    for eta in (0.,0.20273255,-0.20273255):
        for n in (48,96,192):
            for j in range(n):
                th=(j+0.5)*math.pi/n
                gamma=math.cosh(eta)
                f=radius/math.sqrt(gamma**2*math.cos(th)**2+math.sin(th)**2)
                # Independent Cartesian evaluation of the rest-frame radius.
                rest=math.hypot(gamma*f*math.cos(th),f*math.sin(th))
                worst=max(worst,abs(rest/radius-1));count+=1
assert worst<5e-14
print('Domain R>0, gamma>=1, theta real; D>=1; positive radial branch')
print('Witness: cancelled rational residual and zero positivity-identity polynomial')
print(f'Float64 transfer cross-check: {count} declared midpoint samples, no RNG, worst relative={worst}')
print('PROVED: boosted stationary surface seed solves the Lorentz-transformed sphere equation')
print('CAS scope: algebraic identity — production physics not certified')
