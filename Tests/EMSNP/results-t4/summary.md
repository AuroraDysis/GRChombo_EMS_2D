# NPTET T4 interpolation-order comparison

**Ready: False.**

The production order remains Lagrange<4>. The estimator is the sum of two within-route changes: E = |NP4-NP6| + |W4-W6|. The acceptance is |NP4-W4| <= 2 E. Safety factor 2 was fixed in contract.json before computing the study. All production W cell derivatives remain fourth order; only the interpolation order changes. This is a measured discretization estimate at one checkpoint, not a rigorous error bound.

The original three ghosts and point transfer are retained. The existing Lagrange stencil search is bounded to the sampled FAB, shifting full-width stencils near its edges. Every sample is bounds checked; the unbounded L6 negative control rejects an out-of-array stencil. All valid checkpoint data and hierarchy levels are unchanged. L4 reproduces T3, and checked/unchecked L4 agree within 4096 eps times each column scale.

| form | r | l | difference L2 | difference L4 | difference L6 | NP4-NP6 | W4-W6 | 2E | relative difference | relative budget | result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| vacuum | 50 | 2 | 1.207230e-07 | 4.094436e-10 | 5.842498e-10 | 6.572985e-10 | 4.824922e-10 | 2.279582e-09 | 1.423226e-06 | 7.923826e-06 | PASS |
| vacuum | 50 | 4 | 4.202409e-08 | 3.305362e-11 | 1.157639e-10 | 2.326968e-10 | 1.499865e-10 | 7.653666e-10 | 7.079980e-05 | 1.639391e-03 | PASS |
| vacuum | 50 | 6 | 1.687466e-09 | 3.669184e-11 | 4.419197e-11 | 2.205385e-10 | 1.396547e-10 | 7.203865e-10 | 3.388945e-04 | 6.653659e-03 | PASS |
| vacuum | 75 | 2 | 1.985833e-07 | 2.804404e-10 | 6.501002e-10 | 1.003404e-09 | 6.337439e-10 | 3.274295e-09 | 7.527839e-07 | 8.789163e-06 | PASS |
| vacuum | 75 | 4 | 8.474775e-09 | 1.915874e-10 | 4.746289e-10 | 5.408836e-10 | 2.578421e-10 | 1.597452e-09 | 4.367226e-05 | 3.641384e-04 | PASS |
| vacuum | 75 | 6 | 9.294215e-09 | 2.021934e-10 | 2.707232e-10 | 6.305133e-11 | 5.478430e-12 | 1.370595e-10 | 3.923037e-04 | 2.659283e-04 | FAIL |
| vacuum | 100 | 2 | 2.833972e-07 | 2.943231e-10 | 9.147719e-10 | 1.134018e-09 | 5.135695e-10 | 3.295176e-09 | 6.833395e-07 | 7.650516e-06 | PASS |
| vacuum | 100 | 4 | 1.145752e-07 | 4.825040e-11 | 2.064701e-10 | 2.128420e-10 | 5.462235e-11 | 5.349287e-10 | 1.597513e-05 | 1.771085e-04 | PASS |
| vacuum | 100 | 6 | 7.022713e-08 | 8.605493e-12 | 4.497275e-11 | 2.004593e-11 | 7.362418e-11 | 1.873402e-10 | 2.027263e-04 | 4.413320e-03 | PASS |
| matter | 50 | 2 | 1.207232e-07 | 4.094443e-10 | 5.842499e-10 | 6.572980e-10 | 4.824925e-10 | 2.279581e-09 | 1.423228e-06 | 7.923823e-06 | PASS |
| matter | 50 | 4 | 4.202409e-08 | 3.305355e-11 | 1.157639e-10 | 2.326971e-10 | 1.499867e-10 | 7.653676e-10 | 7.081185e-05 | 1.639676e-03 | PASS |
| matter | 50 | 6 | 1.688274e-09 | 3.669089e-11 | 4.419195e-11 | 2.205375e-10 | 1.396547e-10 | 7.203844e-10 | 3.388834e-04 | 6.653596e-03 | PASS |
| matter | 75 | 2 | 1.985833e-07 | 2.804405e-10 | 6.501002e-10 | 1.003404e-09 | 6.337440e-10 | 3.274295e-09 | 7.527841e-07 | 8.789163e-06 | PASS |
| matter | 75 | 4 | 8.474782e-09 | 1.915874e-10 | 4.746288e-10 | 5.408836e-10 | 2.578422e-10 | 1.597451e-09 | 4.367247e-05 | 3.641401e-04 | PASS |
| matter | 75 | 6 | 9.294137e-09 | 2.021935e-10 | 2.707233e-10 | 6.305149e-11 | 5.478299e-12 | 1.370596e-10 | 3.923040e-04 | 2.659285e-04 | FAIL |
| matter | 100 | 2 | 2.833972e-07 | 2.943231e-10 | 9.147720e-10 | 1.134018e-09 | 5.135695e-10 | 3.295176e-09 | 6.833395e-07 | 7.650516e-06 | PASS |
| matter | 100 | 4 | 1.145752e-07 | 4.825042e-11 | 2.064699e-10 | 2.128419e-10 | 5.462240e-11 | 5.349285e-10 | 1.597516e-05 | 1.771087e-04 | PASS |
| matter | 100 | 6 | 7.022713e-08 | 8.605531e-12 | 4.497271e-11 | 2.004595e-11 | 7.362419e-11 | 1.873403e-10 | 2.027264e-04 | 4.413304e-03 | PASS |

L2-to-L4 reduces every cross-route discrepancy, but the L6 discrepancy is larger than L4 in all 18 comparisons. Increasing only interpolation order leaves the fourth-order cell-Weyl derivative error in place. The l=6,r=75 pair in both forms needs a safety factor about 2.95 on the declared within-route estimate, exceeding the allowed maximum 2. No L2-based larger allowance is substituted after observing this failure.

## Extraction cost on this checkpoint

Uninstrumented production L4, 3 spheres, 97 theta nodes, serial Chombo, OMP_NUM_THREADS=8. Median and range of three calls after one warm-up; includes construction, interpolation, reconstruction, projection and CSV IO. Excludes restart loading, shared ghost refresh and shared cell-Weyl preparation. NP is the incremental cost when existing EMS radiation is already enabled.

| form | route | median seconds | min | max |
|---|---|---:|---:|---:|
| vacuum | production | 0.004103 | 0.003775 | 0.005316 |
| vacuum | np | 0.007996 | 0.006917 | 0.008209 |
| matter | production | 0.003290 | 0.003256 | 0.003298 |
| matter | np | 0.006323 | 0.006258 | 0.006412 |

All raw modes, per-call times, parameters, input hashes and logs accompany this report. T3 static, algebraic, rotation, normalization and default-output controls remain accepted; only the Weyl4 comparison tolerance is superseded.
