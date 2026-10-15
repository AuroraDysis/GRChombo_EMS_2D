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

## NPTET T3: authorized trace-free electric/magnetic construction

Artifact: `scripts/cas/np_tracefree.wls`, clean-kernel output in
results-t3/tracefree-cas.log. The contract is a real oriented orthonormal
frame of signature (-+++), arbitrary symmetric trace-free E and B, an
arbitrary real null-rotation parameter, and a general antisymmetric F.
Exact arithmetic is polynomial over Q extended by I and positive sqrt(2).
All denominators are nonzero constants; there are no parameter poles.

The witnesses are all 256 Riemann symmetry/Bianchi component residuals,
all 16 metric traces, and every coefficient of both NP rotation laws
computed from the rebuilt tensor and the rotated tetrad. A second generic
Gauss/Codazzi tensor retains independent Hamiltonian and momentum
residuals and field Ricci entries; direct contraction proves preservation
of Psi0 and Psi4 after the electric/magnetic STF construction. Setting
its Ricci entries to zero includes the consistent vacuum control.
Expand followed by literal zero assertions is the deterministic CLOSE;
failed assertions or Wolfram messages exit nonzero. No simplifier or
numerical rationalization is used.

Independent Float64 transfer consists of the C++ generic E/B symmetry
and Bianchi test (one general ten-component sample, tolerance 1e-14),
the vacuum/matter same-jet Weyl4 unit check (1e-13 absolute), and the
frozen checkpoint on 665 interior angular nodes per form. Four-point
finite differences rotate the actual tetrad and agree with the separate
CAS coefficients within the contraction roundoff budget. All point data,
budgets and both forms are retained. No random sampling is involved.
Psi0/Psi4 preservation is checked for both forms on every requested
sphere; the production-path comparison remains a separate failed gate.

Status: **PROVED** for the algebraic Weyl identities, both rotation laws
and radiative preservation under the stated Gauss/Codazzi construction.
Static and frozen results corroborate the finite-resolution transfer.
They do not close the requested pure-arithmetic comparison of distinct
spatial evaluation orders: the even production W modes fail that gate.
Dimensions remain [E_Weyl]=[B_Weyl]=[H]=[M]=L^-2,
[Phi]=L^-1, [X]=L^-3. STF and the reconstruction are linear in curvature
with dimensionless metric/frame factors; no new physical scale enters.

CAS scope: algebraic identity — production physics not certified

## NPTET T4: measured interpolation-order comparison

The T3 algebraic construction and its accepted static, rotation and default
controls are unchanged. T4 varies only the interpolation order of the two
extraction routes on the same authenticated t=602 hierarchy. Both forms,
all three requested radii and l=2,4,6 are measured at orders 2,4,6 in
Float64. Raw NP/radiation modes and the complete order changes are in
results-t4/study. No independent continuum solution is inferred from them.

Before the valid numerical study, contract.json fixed
E=abs(NP4-NP6)+abs(W4-W6), tolerance=2E. Each summand is a within-route
change; factor 2 was fixed as the safety margin for interpreting those
changes as errors of two approximations. It is an empirical estimate,
not a universal bound. All production W cell derivatives remain fourth
order. The two L4 routes match T3 and the ordinary uninstrumented L4
implementation within 4096 eps times each output column's scale.

The existing point-transfer law requires three ghosts. L6 second
derivatives can reach outside an owning FAB with that allocation. The
offline source therefore bounds the existing Lagrange stencil search to
the sampled FAB and checks every read; the unbounded negative control
exits before an invalid read. This shifts full-width stencils at box
edges without changing the hierarchy, transfer law, normalization or
production L4 interpolation. This verifies valid reads for the present
production-sphere study; earlier unbounded L6 artifacts remain unvalidated
and unused. It is not a continuum convergence certificate.

Measured gate: **16/18 PASS, 2/18 FAIL**. Both failures are l=6,r=75:
difference about 2.02193e-10 versus 1.37060e-10 allowed, relative to the
mode 3.92304e-4 versus 2.65928e-4. The needed multiplier of E is about
2.95, above the allowed 2. L2-to-L4 discrepancies decrease in all cases;
L6 discrepancies exceed L4 in all cases; cell-derivative order remains
fixed during interpolation variation. No larger L2 allowance
was substituted after seeing the failure. Readiness remains false.

Status: **CORROBORATED** for this finite checkpoint/order measurement;
no CAS proof or continuum error-bound claim is made. All sampled values,
controls and budgets are preserved; no RNG or Julia process was used.
