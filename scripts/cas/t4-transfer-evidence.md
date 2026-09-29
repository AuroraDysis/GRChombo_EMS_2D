# T4 finite transfer witness

CAS scope: algebraic identity — production physics not certified

Claim contract: 2:1 cell-centred point values, six consecutive distinct nodes, degree at most five in each coordinate; autonomous smooth vector ODE; real fine-step start fraction theta; RK stage equality modulo coarse timestep H^4. Arithmetic is exact Q for moments and Q[f,B,C,D,theta][H]/(H^4), with B=F'F, C=F''(F,F), D=F'^2F independent formal elementary differentials. There are no poles or branch choices.

Witness: all six monomial moments for each child offset and for midpoint restriction; all four H-coefficient residuals for all four stage values. Tensor exactness follows by the product of the two one-dimensional moment equalities. Degree-zero exactness also proves that accumulation anchored on a current source cell is algebraically identical. For the stage witness, the independent defining expressions are the fine-step RK4 Taylor stages, McCorquodale–Colella (2011), equations 43–49. Substituting the dense coarse jet gives an equivalent coefficient-vector equality; no threshold or condition is discarded.

Simplifier role: deterministic expansion and extraction of H^0 through H^3 are CLOSE operations, each consumed by an assert that exits nonzero on failure. No general simplifier, NORMALIZE transformation or branch-crossing operation is used.

Status: PROVED for the stated finite moments and truncated stage identities. The existing Chombo stage-2/stage-3 correction factors give the exact residuals -D H^3/128 and +D H^3/64 at refinement ratio 2. The local wrapper changes only those correction factors from r^2 to r^3.

Independent runtime falsifier: the compiled production operator is compared against direct point samples of all 36 tensor monomials at N=32,64,128, with even/odd parity at y=0, and against exp(2x) cos(3y) and exp(2x) sin(3y). The compiled stage wrapper is compared against separately constructed fine RK stages starting at exp(0.7 theta H), theta=0,1/2, H=1/8,1/16,1/32, in Float64. These deterministic checks have no RNG. Worst polynomial prolongation/restriction errors are 1.3877787807814457e-17 and 3.4694469519536142e-18; stage errors are 1.0732744670782779e-6, 6.6931720610341472e-8, 4.1786385462927456e-9. These checks exercise the transfer implementation; they do not certify the EMS constraints or production physics.

Command: `/Users/auroradysis/miniconda3/bin/python scripts/cas/t4-transfer-verify.py > scripts/cas/t4-transfer-verify.json`. Runtime command and all measured residuals are recorded in `Tests/EMSNative/t4-runs.csv` and `t4-operator-checks.csv`.
