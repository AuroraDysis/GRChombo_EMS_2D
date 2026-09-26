# EMSTRUMPET C3 evolution results

Source profile: `reference-alpha20-qfile.trumpet`, SHA-256
`6f0820a576620f1f7230c701131af56a2312b130e31d5de2a752c32cefe6d24f`.
All runs used `OMP_NUM_THREADS=1`, CFL 0.25, the existing CCZ4/gauge/RK4 path,
and `min_chi = min_lapse = 1e-12`. Commands and masks are in the parent README;
the CSV files here contain unrounded norms and orders.

| Run | Grid and steps | Wall time | Result |
|:---|:---|---:|:---|
| Single | M/16, 16 | 0.74 s | finite, T=0.25 |
| Single | M/32, 32 | 5.31 s | finite, T=0.25 |
| Single | M/64, 64 | 39.88 s | finite, T=0.25 |
| Boosted single, η=+0.20273255 | M/32, 32 | 5.27 s | finite, T=0.25 |
| Binary, d=32, η=±0.20273255 | M/16, 16 | 3.29 s | finite, T=0.25 |
| Binary | M/32, 32 | 27.27 s | finite, T=0.25 |
| Binary | M/64, 64 | 200.83 s | finite, T=0.25 |
| AMR, levels 0–2 | base M/8, 4 coarse steps | 4.08 s | finite; levels 1 and 2 regridded |
| RHFinder, N=32 | 32 steps | 5.39 s | first output at t=1/128; T=0.25 |

The corrected t=0 axis strip is the **fixed** `0 < y ≤ 1/8` within
`0.25 ≤ R_rest ≤ 1.5`, giving 4, 8, 16 y rows at N=32, 64, 128.
[`fixed-axis-strip.csv`](fixed-axis-strip.csv) has all nine L2/L∞ norms and
successive orders. Its finer-pair nonzero L2/L∞ orders are:

| η | H | Mx | My | GaussE |
|---:|:---|:---|:---|:---|
| 0 | 3.99 / 3.94 | 3.97 / 3.85 | 3.98 / 3.97 | 3.96 / 4.19 |
| +0.20273255 | 3.97 / 3.85 | 3.91 / 3.63 | 4.01 / 3.87 | 4.36 / 4.30 |
| −0.20273255 | 3.97 / 3.85 | 3.91 / 3.63 | 4.01 / 3.87 | 4.36 / 4.30 |

## Uniform single evolution

Field ratios `Q = ||u_16 − I_32 u_32|| / ||I_32 u_32 − I_64 u_64||`
use the same 678 N=16 cell centres in `0.75 ≤ R ≤ 1.5`; `I` is tensor-product
degree-five Lagrange interpolation. These are L2 ratios, in field order
`chi, h11, K, A11, lapse, shift1, phi, Pi, Ex`:

| t | χ | hxx | K | Axx | lapse | shiftx | φ | Π | Ex |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1.01 | 1.03 | 1.42 | 1.00 | 1.00 | 1.01 | 1.01 | 1.01 | 0.99 |
| 1/16 | 16.21 | 15.85 | 18.93 | 15.68 | 15.41 | 16.24 | 16.21 | 15.95 | 14.40 |
| 1/8 | 16.77 | 15.94 | 17.79 | 17.42 | 15.76 | 15.60 | 16.42 | 15.65 | 14.86 |
| 1/4 | 16.98 | 15.86 | 14.96 | 16.27 | 15.53 | 15.59 | 16.19 | 15.56 | 15.19 |

[`field-convergence.csv`](field-convergence.csv) also gives both difference
norms, L∞ ratios and orders. At t=0 the initializer samples an analytic
profile: there is no evolution truncation error, and the two differences are
dominated by the same N=32 interpolation error. Its Q≈1 is not an evolution
order. For t>0 all L2 orders are 3.85–4.24; the one L∞ exception to 3.5–4.5
is K at t=1/16, order 4.57.

[`common-constraints.csv`](common-constraints.csv) evaluates constraint
diagnostics on those same physical points. Its L2 norms and finer-pair orders:

| t | Quantity | N16 | N32 | N64 | p32→64 |
|---:|:---|---:|---:|---:|---:|
| 0 | H | 1.560e-5 | 9.782e-7 | 6.118e-8 | 4.00 |
| 0 | Mx/My | 3.918e-5 | 2.430e-6 | 1.516e-7 | 4.00 |
| 0 | GaussE | 1.741e-7 | 1.048e-8 | 6.507e-10 | 4.01 |
| 1/16 | H | 2.203e-5 | 1.137e-6 | 5.520e-8 | 4.36 |
| 1/16 | Mx | 3.845e-5 | 2.388e-6 | 1.490e-7 | 4.00 |
| 1/16 | My | 3.943e-5 | 2.445e-6 | 1.519e-7 | 4.01 |
| 1/16 | GaussE | 1.608e-7 | 9.689e-9 | 6.034e-10 | 4.01 |
| 1/16 | Θ | 2.605e-7 | 1.455e-8 | 7.584e-10 | 4.26 |
| 1/8 | H | 2.637e-5 | 1.605e-6 | 1.060e-7 | 3.92 |
| 1/8 | Mx | 3.796e-5 | 2.361e-6 | 1.485e-7 | 3.99 |
| 1/8 | My | 4.536e-5 | 3.146e-6 | 2.148e-7 | 3.87 |
| 1/8 | GaussE | 1.504e-7 | 9.083e-9 | 5.683e-10 | 4.00 |
| 1/8 | Θ | 5.786e-7 | 3.342e-8 | 1.667e-9 | 4.33 |
| 1/4 | H | 3.513e-5 | 2.496e-6 | 1.713e-7 | 3.87 |
| 1/4 | Mx | 3.714e-5 | 2.312e-6 | 1.452e-7 | 3.99 |
| 1/4 | My | 5.848e-5 | 5.171e-6 | 3.331e-7 | 3.96 |
| 1/4 | GaussE | 1.415e-7 | 8.675e-9 | 5.479e-10 | 3.98 |
| 1/4 | Θ | 1.408e-6 | 9.281e-8 | 5.692e-9 | 4.03 |

GaussB is identically zero, and Θ is zero at t=0. The common-point L∞ order
exceeds 4.5 for H at t=1/16 (4.63) and Θ at t=1/8 (4.90). The
[`single-constraints.csv`](single-constraints.csv) retains norms on every
grid's own fixed-radius mask. Its L∞ maxima sometimes show lower apparent
orders: the maximizing cells change with resolution at the inner mask edge
and first axis row. Those raw rates are preserved in the CSV.

## Smokes

The boosted single at N=32 has H L2 `1.022e-6 → 2.672e-6`, Mx
`2.556e-6 → 2.342e-6`, My `2.464e-6 → 4.508e-6`, GaussE
`1.124e-8 → 1.085e-8`, and Θ `0 → 1.200e-7` from t=0 to T.
[`single-constraints.csv`](single-constraints.csv) includes the intermediate
times and L∞ norms.

The binary's union-of-shells H L2 is `7.331e-3 → 7.148e-3`, Mx
`1.149e-3 → 1.148e-3`, My `1.035e-3 → 1.018e-3`, and GaussE
`5.407e-4 → 5.242e-4`. This nonzero superposition residual is expected.
At t=0 and T, the largest parity-adjusted normalized mirror residual across
all 33 plotted fields is `2.22e-13` and `1.24e-12`; see
[`binary-symmetry.csv`](binary-symmetry.csv). Binary evolved-field L2
Richardson ratios are 14.75–16.90 at T, orders 3.88–4.08; see
[`binary-field-convergence.csv`](binary-field-convergence.csv). Its t=0 ratios
have the same interpolation limitation as the single-object data.

The AMR run logged repeated level-1 and level-2 regrids. On the exterior
shell at T, level 0/1/2 H L2 is `2.409e-3 / 4.945e-4 / 2.682e-6` and GaussE L2
is `3.381e-5 / 6.906e-6 / 9.013e-9`. These are separate level norms,
not a composite AMR norm; see [`amr-levels.csv`](amr-levels.csv).

All plotted components are finite. Across the 41 plotted level snapshots,
[`snapshot-status.csv`](snapshot-status.csv) counts zero interior χ and lapse
values below `1e-12`, including outside R=0.25. These are snapshot counts;
individual RK-substage clamp calls were not instrumented. The existing
`min_chi.dat` output is not a χ minimum: its code constructs
`AMRReductions<VariableType::diagnostic>` then passes `c_chi=0`, which selects
the diagnostic `mod_F` slot. The plots directly give positive χ minima.

The optional RHFinder first writes at t=1/128, not at t=0. Its areal radius
is `1.750533` then, versus the profile's `r_h=1.750590`; at T it reports
`1.745695`, with mean outgoing expansion `−4.93e-3` and mode `close`.
The latter is not a converged apparent-horizon measurement. See
[`rh-results.csv`](rh-results.csv).
