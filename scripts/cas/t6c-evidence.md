# T6 C conditional characteristic algebra

Claim contract: alpha, chi, g=conformal inverse-normal metric, and mu are positive;
beta_n is real. Coefficients are frozen in a homogeneous one-dimensional principal
gauge block, with Theta perturbation zero. Matter forcing, damping, variable
coefficients, cartoon terms and KO are excluded. This is not the full CCZ4 symbol.

Witness: exact coefficient-zero characteristic-polynomial residual in
QQ[alpha,chi,g,mu,beta_n,s], plus the transverse 2x2 determinant. For the source
principal block on (d_n alpha,K-2Theta,d_n beta_n,Gamma_n), the coordinate speeds
are eigenvalues of -P, since u_t=P u_n. The witness factors into
[(s+beta_n)^2-9 alpha chi g/5][(s+beta_n)^2-4 mu g/3].
The transverse factor is (s+beta_n)^2-mu g. Outward branches use positive square
roots and -beta_n+c. With mu=3/4 the longitudinal shift branch is -beta_n+sqrt(g).
The source correspondence is ExperimentalGauge.hpp:84-93 and the K, Theta,
Gamma principal terms of CCZ4Cartoon.impl.hpp. The correspondence and the physical
family assignment remain an interpretation of this restricted block.

Normal form: polynomial coefficient expansion and rational cancellation, with
every CLOSE consumed by an assert. No general-purpose simplifier is needed.
Excluded locus: nonpositive alpha/chi/g/mu; coincident characteristic speeds
cannot support an independent eigenvector interpretation. The polynomial
identity alone makes no diagonalizability or hyperbolicity claim there.

Independent numerical check: mpmath eigensolver versus direct positive-root
speed formulas, 24 pseudorandom points, seed 6019, 60 decimal digits; parameters
alpha,chi,g,mu in [0.2,1.5], beta_n in [-0.2,0.2]. Maximum residual
5.4295808299e-59, below 1e-55. The JSON retains the complete result.

Status: PROVED for these algebraic identities. Scope: **algebraic identity —
production physics not certified**. No characteristic projection or generated
production kernel is introduced. Matching a measured crest speed does not prove
a pure mode in the nonlinear coupled system.

The ray sampler reuses the T5 Lagrange witness and adds a runnable Float64 check
of all 36/64 tensor monomials, all 33 components, and even/odd axis parity at
six/eight nodes (tolerance 3e-13). A separate small assertion protects the local
index of the independently selected temporal peak's growth basin. These checks
do not certify physical pulse resolution or AMR evolution.

Reproduce:

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/auroradysis/miniconda3/bin/python scripts/cas/t6c-characteristics.py > scripts/cas/t6c-characteristics.json
```
