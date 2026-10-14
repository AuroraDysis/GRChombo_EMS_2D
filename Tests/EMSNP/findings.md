# NPTET T1 exact-witness evidence card

Claim: at fixed perturbed fields the two defining X combinations have
O(a^2) change under an infinitesimal null rotation. Domain: characteristic
zero, normalized null tetrad, a real (a valid subset of complex rotations),
arbitrary complex NP scalar components. Counterexample components are
real, arbitrarily small perturbations of aligned Coulomb fields. Geometric
mass-one units; restoring M gives X weight M^-3. No poles or branch cuts.

Witness: coefficient of a^1 of the *unchanged* defining products after
binomial Weyl/Maxwell spinor transformation. A polynomial difference is
O(a^2) iff this coefficient vanishes (its constant term is zero by
construction). The four coefficients and a nonzero exact counterexample
are in `../../scripts/cas/np_null_rotation.wls` and `results/cas.log`.
No witness substitution or threshold is used for the exact refutation.

Simplifier roles: none. `Expand` provides deterministic polynomial CLOSE
in Q[a,Psi0..4,Phi0..2]; residual identities and the counterexample are
consumed by explicit tests exiting nonzero on failure. No general
simplifier, conditional identity, numerical rationalization, or branch
crossing occurs. All preconditions are discharged by this polynomial
domain and the stated normalized tetrad.

Independent falsifier/transfer check: direct rank-4/rank-2 tensor
contractions before/after tetrad rotation, rather than applying the
binomial rules, in `check_rotation.py`. Its deterministic electric/magnetic
Weyl fixture (no RNG) verifies both rotation families at a=.01,-.03;
worst Float64 residual 5.929e-17, tolerance 4096 eps=9.095e-13 for O(1)
fixture fields. The actual t=602 checkpoint jet provides the numerical
failed control with its arithmetic error budget, trace residual, and
scope limitations in the README and `results/rotation.json`.

Backend: Wolfram 15.0.1, via CLI-fallback (no Wolfram MCP evaluator exposed).
Fresh-kernel final adjudication completed in .75 s. Context7 resolved the
official `/websites/reference_wolfram_language` reference and confirmed
polynomial expansion/coefficient access; no external code was loaded.

Status: **REFUTED** for the universal fixed-perturbed-field O(a^2) claim.
The coefficient identities themselves are closed exact witnesses;
first-order Coulomb-background cancellation is retained. This status
does not assert a continuum error estimate for the supplied hierarchy.

CAS scope: algebraic identity — production physics not certified

## NPTET T2: corrected first-order invariance

The controller explicitly withdrew the T1 fixed-field O(a^2) claim.
This card covers the corrected joint perturbative statement, not a retry
of the refuted proposition. Artifact: `scripts/cas/np_first_order.wls`;
resident output in `results-t2/first-order-cas-mcp.txt`, cold-kernel final
output in `results-t2/first-order-cas.log`.

Contract/domain: polynomial ring Q[a,h,eps,Psi0..4,Phi0..2], h!=0;
arbitrary complex scalar components, real rotation direction. Count a
and every non-Coulomb component as order eps, with independent Coulomb
Psi2 and Phi1. No analytic branches, poles except h=0, or hidden kernel
assumptions. Witnesses: coefficients of eps^0 and eps^1 of the exact
rotated-minus-original products; both vanish in all four cases. The
four-point difference minus the first polynomial coefficient also
vanishes identically. A separate matrix contraction proves the residual
Weyl trace identity. No witness substitution or fitted threshold is used.

Normalizer roles: Expand and Cancel[Together[...]] are deterministic
polynomial/rational CLOSE operations. Every residual is consumed by an
explicit assertion exiting nonzero. h!=0 is the sole denominator
condition in the finite-difference witness. Numeric transfer: direct
tetrad/tensor contraction at 665 nonzero-weight nodes, both classes and
both X, for each Riemann form. All points, values, FD steps and scale-based
roundoff budgets are preserved in the two np_points.csv packets. No RNG.

Status: **PROVED** for the algebraic coefficient, joint first-order and
finite-difference identities. The runtime comparison does not close:
the matter tensor reconstructed with the required field-equation Ricci
has residual trace, and every sphere/class/combination group has points
outside the roundoff budget. The algebraic status does not certify that
runtime transfer or the discretized physics.

CAS scope: algebraic identity — production physics not certified

## NPTET T2: static metric oracle

Artifact: `scripts/cas/np_static_metric.wls`, with resident and fresh CLI
outputs in results-t2/static-cas-mcp.txt and static-cas.log. Contract:
static diagonal metric (-A,1/N,r^2,r^2 sin(theta)^2), r,A,N,F>0,
0<theta<pi, delta real; poles r=A=N=F=sin(theta)=0 excluded. The ring is
rational in the metric functions and their jets, with Exp(delta)
nonzero. The only square root is the positive constant sqrt(8 pi).

Witness: four directly differentiated metric Riemann components minus
the displayed orthonormal formulas, followed by the Ricci-subtracted
Psi2 formula and exact Schwarzschild/RN substitutions. Cancel[Together]
is a CLOSE on each residual; explicit boolean assertions exit nonzero.
No general simplifier, branch crossing or threshold tuning is used.
Independent numeric transfer: the Julia native spectral metric/potential
oracle and the production Cartesian 3+1 t=0 hierarchy, at all seven
static radii and two resolutions. RN is additionally checked against its
closed formula. Spatial and conditioning-based roundoff budgets, points,
precision (Float64), errors and orders are in validation.json; no RNG.

Status: **PROVED** for the static metric/NP identities. The 70 static
scalar/sphere checks corroborate the two finite-resolution profiles;
this is not a parametric continuum certificate. The outermost EMS
curvature has reached the Float64 derivative arithmetic floor.

CAS scope: algebraic identity — production physics not certified
