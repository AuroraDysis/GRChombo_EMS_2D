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
