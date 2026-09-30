# T5 sampling algebra

Claim contract: six/eight distinct consecutive integer stencil nodes, real target q, tensor polynomial degree at most five/seven per coordinate; ratio=3/2; psi>0, alpha>=0. The speed comparison uses chi=psi^-4 and positive square roots. Excluded locus: psi=0 and coincident stencil nodes.

Witness: coefficient-zero Lagrange moment residuals in QQ[q]; equality of squared speed expressions with nonnegative branches; exact zero of the scaled fourth-order difference residual. Equality of nonnegative squares is equivalent to equality of the speed expressions. No numerical acceptance threshold is fitted to the simulation data.

Normal form and simplifier roles: deterministic expand and cancel, each CLOSE consumed by an assert; no general simplifier and no unverified normalization.

Branch conditions: psi>0 ensures finite inverse powers; alpha>0 in the numerical samples. Alpha=0 is the continuous zero-speed boundary. These assumptions are checked against the saved positive chi/lapse values when speeds are evaluated.

Independent checks: two speed definitions at 60 decimal digits, 32 points, seed 1729, psi in [0.5,10], alpha in [0.0001,1]; maximum residual is in t5-sampling-verify.json. t5-analyze.py selfcheck independently exercises the Float64 tensor sampler on all 36/64 monomials with even/odd axis parity at both stencil orders, tolerance 3e-13.

Status: PROVED for the algebraic witnesses. Sampling checks certify neither gauge-mode identity nor AMR/evolution correctness.

CAS scope: algebraic identity — production physics not certified.

Reproduce: PYTHONDONTWRITEBYTECODE=1 /Users/auroradysis/miniconda3/bin/python scripts/cas/t5-sampling-verify.py > scripts/cas/t5-sampling-verify.json
