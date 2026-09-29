# T4b algebra evidence

Witness: over Q, eight consecutive-node Lagrange moments have zero residual for
degrees 0–7 at both quarter-cell positions. Tensor moments factor into these
one-dimensional witnesses. The log-root identity is
`f(x+delta)-f(x)+delta=0`, for `f=log(u)/nu-log(1-u)-x`.

Domain: distinct stencil nodes; 0<u<1, nu>0; x=log(R/ell), R>0, ell>0.
Excluded loci: coincident nodes, u=0 or 1, nu=0. No branch changes or general
simplifier: rational arithmetic and deterministic expansion close each witness.

Independent Float64 checks of the implemented offline P8 tensor interpolation
cover all 64 monomials of degree at most seven per coordinate. Maximum absolute
error is 2.168404344971009e-19. The log-coordinate numerical check is recorded in
the JSON transcript. Assertions fail with a nonzero exit status.

Reproduce: `PYTHONDONTWRITEBYTECODE=1 /Users/auroradysis/miniconda3/bin/python
scripts/cas/t4b-verify.py`. Transcript: `t4b-verify.json`.

Scope: algebraic identities and the offline probe. Production physics and
roundoff-floor saturation are not certified. No static solution is evaluated.
