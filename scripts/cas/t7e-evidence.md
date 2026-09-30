# T7-E radial initial-data curvature

**Claim contract.** For a smooth real chi(r)>0 and r>0, the scalar curvature of the three-dimensional metric gamma_ij=delta_ij/chi(r) is 2(chi_rr+2 chi_r/r) - 5 chi_r^2/(2 chi). The spherical-chart witness excludes sin(theta)=0; the resulting scalar is independent of theta and extends to the regular axis for r>0. Equality is literal in rational differential expressions. No statement is made at the puncture or about the evolution.

**Witness and method.** Construct the diagonal spherical metric, its inverse, all Christoffels and the Ricci contraction independently. Subtract the stated radial expression; expand trigonometric products and cancel the rational residual. The visible assert requires exactly zero. There is no general simplifier, numerical threshold substitution, or omitted denominator. Positivity of chi is a hypothesis and is checked on the evaluated checkpoint samples. No square-root branch is crossed in this exact witness.

**Independent falsifier.** A separate Cartesian implementation contracts Christoffels and their derivatives from radial metric jets for chi(r)=1+r^2+exp(-r)/10. Twelve positive-coordinate points, seed 945825, 60 decimal digits, give maximum disagreement 3.11150763893057e-60. The JSON records the contract and result. This is a numerical cross-check, not an extra proof of the production implementation.

**Status: PROVED.** CAS scope: algebraic identity — production physics not certified.

The runtime transfer is checked separately: the standalone T7EInitial program calls the actual EMSTRUMPET reader and the actual setter, verifies the same coupling, and evaluates analytic radial jets using the reader's derivative coefficients. The source profile SHA-256 equals exp-0019's manifest. The checkpoint's t=0 fields agree with the setter at ordinary Float64 accuracy. Only the t=0 branch may open this profile; the evolved checkpoint replay has no initial-data reader. Frozen coordinate controls retain the current t=0 fields other than chi and use the unchanged native constraint kernel. They are diagnostic operation comparisons, not evolution modifications or after-t=0 targets.

Reproduce: `OPENBLAS_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python scripts/cas/t7e-initial-verify.py`. The command and runtime conditions are registered in README T7-E. Exact closure certifies the geometric identity only; the CSVs carry the numerical evidence and its sampling limits.
