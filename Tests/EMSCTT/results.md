# Measured C++ handoff results

Decision: fixture parity and no-companion bit identity pass. The CTT constraints show fourth-order decay before the declared numerical floors; the density seed retains its nonzero continuum residual. No evolution or HPC claim is made.

Errors use `abs(C++ − Julia)/(1 + abs(Julia)) ≤ 5e-13`, the existing density-binary fixture tolerance. All 42 points pass in each table.

| Table | Group | Max scaled error | Worst point | Component |
|---|---|---:|---:|---|
| physical | B | 6.4988e-19 | 10 | B_z |
| physical | C | 6.4957e-19 | 10 | C_yy |
| physical | E | 2.2700e-15 | 3 | E_y |
| physical | Kij | 5.9062e-15 | 13 | K_xx |
| physical | gamma | 4.2597e-15 | 3 | gamma_yy |
| physical | logpsi | 1.2196e-19 | 30 | logpsi |
| physical | scalar | 2.0419e-17 | 15 | Pi |
| physical | shift | 9.1138e-17 | 31 | shift_x |
| ccz4 | B | 6.4988e-19 | 10 | Bz |
| ccz4 | E | 2.2700e-15 | 3 | Ey |
| ccz4 | auxiliary | 0.0000e+00 | all exact | all exact |
| ccz4 | curvature | 1.1141e-15 | 2 | A22 |
| ccz4 | gauge | 9.1138e-17 | 31 | shift1 |
| ccz4 | metric | 1.3972e-15 | 2 | h22 |
| ccz4 | scalar | 2.0419e-17 | 15 | Pi |

Each entry below is unweighted cell RMS (L2), with `log2(previous/current)` in parentheses after the first row. `M = hypot(Mx, My)`. GaussB is exactly zero at every level; its order is undefined. Full-precision L2 and L∞ values, both momentum components, orders, counts, staging errors and the reflection diagnostic are in [convergence.csv](convergence.csv).

The core spacings are 1/16, 1/32, 1/64, 1/128 M; the outer spacings are 8, 4, 2, 1 M. The 1/16 level was added to expose the asymptotic inter-hole H range before the finest-grid numerical floor. All masks and fields are fixed across the final exterior ladder.

## left

| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CTT | 0.0625 | 666 | 1.6194e-05 | 4.0384e-05 | 3.9942e-05 | 5.6800e-05 | 1.8074e-07 | 0.0000e+00 |
| density | 0.0625 | 666 | 3.6351e-03 | 1.1861e-03 | 1.0376e-03 | 1.5759e-03 | 1.8381e-07 | 0.0000e+00 |
| CTT | 0.03125 | 2668 | 1.0254e-06 (3.981) | 2.5036e-06 (4.012) | 2.4112e-06 (4.050) | 3.4759e-06 (4.030) | 1.0854e-08 (4.058) | 0.0000e+00 |
| density | 0.03125 | 2668 | 3.6237e-03 (0.005) | 1.1830e-03 (0.004) | 1.0378e-03 (-0.000) | 1.5737e-03 (0.002) | 1.1045e-08 (4.057) | 0.0000e+00 |
| CTT | 0.015625 | 10636 | 6.3669e-08 (4.009) | 1.5424e-07 (4.021) | 1.5093e-07 (3.998) | 2.1581e-07 (4.010) | 6.7774e-10 (4.001) | 0.0000e+00 |
| density | 0.015625 | 10636 | 3.6243e-03 (-0.000) | 1.1829e-03 (0.000) | 1.0381e-03 (-0.000) | 1.5738e-03 (-0.000) | 6.8968e-10 (4.001) | 0.0000e+00 |
| CTT | 0.0078125 | 42540 | 3.9910e-09 (3.996) | 9.6914e-09 (3.992) | 9.4695e-09 (3.994) | 1.3550e-08 (3.993) | 4.2470e-11 (3.996) | 0.0000e+00 |
| density | 0.0078125 | 42540 | 3.6256e-03 (-0.000) | 1.1837e-03 (-0.001) | 1.0394e-03 (-0.002) | 1.5753e-03 (-0.001) | 4.3218e-11 (3.996) | 0.0000e+00 |

## right

| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CTT | 0.0625 | 666 | 1.6194e-05 | 4.0384e-05 | 3.9942e-05 | 5.6800e-05 | 1.8074e-07 | 0.0000e+00 |
| density | 0.0625 | 666 | 3.6351e-03 | 1.1861e-03 | 1.0376e-03 | 1.5759e-03 | 1.8381e-07 | 0.0000e+00 |
| CTT | 0.03125 | 2668 | 1.0254e-06 (3.981) | 2.5036e-06 (4.012) | 2.4112e-06 (4.050) | 3.4759e-06 (4.030) | 1.0854e-08 (4.058) | 0.0000e+00 |
| density | 0.03125 | 2668 | 3.6237e-03 (0.005) | 1.1830e-03 (0.004) | 1.0378e-03 (-0.000) | 1.5737e-03 (0.002) | 1.1045e-08 (4.057) | 0.0000e+00 |
| CTT | 0.015625 | 10636 | 6.3669e-08 (4.009) | 1.5424e-07 (4.021) | 1.5093e-07 (3.998) | 2.1581e-07 (4.010) | 6.7774e-10 (4.001) | 0.0000e+00 |
| density | 0.015625 | 10636 | 3.6243e-03 (-0.000) | 1.1829e-03 (0.000) | 1.0381e-03 (-0.000) | 1.5738e-03 (-0.000) | 6.8968e-10 (4.001) | 0.0000e+00 |
| CTT | 0.0078125 | 42540 | 3.9910e-09 (3.996) | 9.6914e-09 (3.992) | 9.4695e-09 (3.994) | 1.3550e-08 (3.993) | 4.2470e-11 (3.996) | 0.0000e+00 |
| density | 0.0078125 | 42540 | 3.6256e-03 (-0.000) | 1.1837e-03 (-0.001) | 1.0394e-03 (-0.002) | 1.5753e-03 (-0.001) | 4.3218e-11 (3.996) | 0.0000e+00 |

## middle

| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CTT | 0.0625 | 10752 | 5.0347e-08 | 6.6049e-09 | 5.5372e-09 | 8.6188e-09 | 9.1377e-10 | 0.0000e+00 |
| density | 0.0625 | 10752 | 4.0040e-04 | 1.9429e-05 | 1.2106e-05 | 2.2892e-05 | 9.3169e-10 | 0.0000e+00 |
| CTT | 0.03125 | 43008 | 3.1535e-09 (3.997) | 4.1568e-10 (3.990) | 3.5010e-10 (3.983) | 5.4347e-10 (3.987) | 5.7187e-11 (3.998) | 0.0000e+00 |
| density | 0.03125 | 43008 | 4.0040e-04 (-0.000) | 1.9448e-05 (-0.001) | 1.2116e-05 (-0.001) | 2.2914e-05 (-0.001) | 5.8309e-11 (3.998) | 0.0000e+00 |
| CTT | 0.015625 | 172032 | 1.9775e-10 (3.995) | 2.6028e-11 (3.997) | 2.1947e-11 (3.996) | 3.4046e-11 (3.997) | 3.5754e-12 (4.000) | 0.0000e+00 |
| density | 0.015625 | 172032 | 4.0040e-04 (-0.000) | 1.9452e-05 (-0.000) | 1.2118e-05 (-0.000) | 2.2918e-05 (-0.000) | 3.6455e-12 (4.000) | 0.0000e+00 |
| CTT | 0.0078125 | 688128 | 5.6723e-11 (1.802) | 1.6278e-12 (3.999) | 1.3727e-12 (3.999) | 2.1293e-12 (3.999) | 2.2348e-13 (4.000) | 0.0000e+00 |
| density | 0.0078125 | 688128 | 4.0041e-04 (-0.000) | 1.9453e-05 (-0.000) | 1.2119e-05 (-0.000) | 2.2919e-05 (-0.000) | 2.2786e-13 (4.000) | 0.0000e+00 |

## outer

| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CTT | 8 | 120 | 5.4172e-09 | 5.0033e-09 | 2.6757e-09 | 5.6739e-09 | 1.6311e-09 | 0.0000e+00 |
| density | 8 | 120 | 2.4819e-07 | 5.1043e-09 | 2.8078e-09 | 5.8256e-09 | 1.6315e-09 | 0.0000e+00 |
| CTT | 4 | 496 | 4.6176e-10 (3.552) | 2.9060e-10 (4.106) | 1.7681e-10 (3.920) | 3.4016e-10 (4.060) | 9.6442e-11 (4.080) | 0.0000e+00 |
| density | 4 | 496 | 2.5727e-07 (-0.052) | 3.3974e-10 (3.909) | 1.3407e-09 (1.066) | 1.3830e-09 (2.075) | 9.6468e-11 (4.080) | 0.0000e+00 |
| CTT | 2 | 2002 | 3.1149e-11 (3.890) | 1.8295e-11 (3.990) | 1.0918e-11 (4.017) | 2.1305e-11 (3.997) | 6.0608e-12 (3.992) | 0.0000e+00 |
| density | 2 | 2002 | 2.5911e-07 (-0.010) | 1.4700e-10 (1.209) | 1.3748e-09 (-0.036) | 1.3827e-09 (0.000) | 6.0624e-12 (3.992) | 0.0000e+00 |
| CTT | 1 | 8038 | 2.0043e-12 (3.958) | 1.1531e-12 (3.988) | 6.8546e-13 (3.993) | 1.3414e-12 (3.989) | 3.8223e-13 (3.987) | 0.0000e+00 |
| density | 1 | 8038 | 2.5935e-07 (-0.001) | 1.4453e-10 (0.024) | 1.3762e-09 (-0.001) | 1.3838e-09 (-0.001) | 3.8233e-13 (3.987) | 0.0000e+00 |

## left_axis

| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CTT | 0.0625 | 46 | 3.5812e-05 | 5.7868e-05 | 1.9008e-05 | 6.0910e-05 | 1.9830e-07 | 0.0000e+00 |
| density | 0.0625 | 46 | 3.6104e-03 | 1.6995e-03 | 9.9808e-05 | 1.7024e-03 | 2.0255e-07 | 0.0000e+00 |
| CTT | 0.03125 | 190 | 2.4922e-06 (3.845) | 4.2860e-06 (3.755) | 1.6549e-06 (3.522) | 4.5943e-06 (3.729) | 1.5106e-08 (3.714) | 0.0000e+00 |
| density | 0.03125 | 190 | 3.6118e-03 (-0.001) | 1.7956e-03 (-0.079) | 1.2770e-04 (-0.356) | 1.8001e-03 (-0.080) | 1.5400e-08 (3.717) | 0.0000e+00 |
| CTT | 0.015625 | 754 | 1.5220e-07 (4.033) | 2.5933e-07 (4.047) | 9.6957e-08 (4.093) | 2.7686e-07 (4.053) | 9.4764e-10 (3.995) | 0.0000e+00 |
| density | 0.015625 | 754 | 3.6026e-03 (0.004) | 1.7862e-03 (0.008) | 1.2659e-04 (0.013) | 1.7906e-03 (0.008) | 9.6586e-10 (3.995) | 0.0000e+00 |
| CTT | 0.0078125 | 3020 | 9.5572e-09 (3.993) | 1.6331e-08 (3.989) | 6.1213e-09 (3.985) | 1.7440e-08 (3.989) | 5.9847e-11 (3.985) | 0.0000e+00 |
| density | 0.0078125 | 3020 | 3.6037e-03 (-0.000) | 1.7884e-03 (-0.002) | 1.2702e-04 (-0.005) | 1.7929e-03 (-0.002) | 6.0992e-11 (3.985) | 0.0000e+00 |

## right_axis

| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CTT | 0.0625 | 46 | 3.5812e-05 | 5.7868e-05 | 1.9008e-05 | 6.0910e-05 | 1.9830e-07 | 0.0000e+00 |
| density | 0.0625 | 46 | 3.6104e-03 | 1.6995e-03 | 9.9808e-05 | 1.7024e-03 | 2.0255e-07 | 0.0000e+00 |
| CTT | 0.03125 | 190 | 2.4922e-06 (3.845) | 4.2860e-06 (3.755) | 1.6549e-06 (3.522) | 4.5943e-06 (3.729) | 1.5106e-08 (3.714) | 0.0000e+00 |
| density | 0.03125 | 190 | 3.6118e-03 (-0.001) | 1.7956e-03 (-0.079) | 1.2770e-04 (-0.356) | 1.8001e-03 (-0.080) | 1.5400e-08 (3.717) | 0.0000e+00 |
| CTT | 0.015625 | 754 | 1.5220e-07 (4.033) | 2.5933e-07 (4.047) | 9.6957e-08 (4.093) | 2.7686e-07 (4.053) | 9.4764e-10 (3.995) | 0.0000e+00 |
| density | 0.015625 | 754 | 3.6026e-03 (0.004) | 1.7862e-03 (0.008) | 1.2659e-04 (0.013) | 1.7906e-03 (0.008) | 9.6586e-10 (3.995) | 0.0000e+00 |
| CTT | 0.0078125 | 3020 | 9.5572e-09 (3.993) | 1.6331e-08 (3.989) | 6.1213e-09 (3.985) | 1.7440e-08 (3.989) | 5.9847e-11 (3.985) | 0.0000e+00 |
| density | 0.0078125 | 3020 | 3.6037e-03 (-0.000) | 1.7884e-03 (-0.002) | 1.2702e-04 (-0.005) | 1.7929e-03 (-0.002) | 6.0992e-11 (3.985) | 0.0000e+00 |

## middle_axis

| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CTT | 0.0625 | 896 | 8.1991e-08 | 1.2347e-08 | 1.4813e-09 | 1.2436e-08 | 1.5143e-09 | 0.0000e+00 |
| density | 0.0625 | 896 | 4.3390e-04 | 2.6428e-05 | 1.4623e-06 | 2.6469e-05 | 1.5450e-09 | 0.0000e+00 |
| CTT | 0.03125 | 3584 | 5.1403e-09 (3.996) | 7.9065e-10 (3.965) | 1.0067e-10 (3.879) | 7.9704e-10 (3.964) | 9.5092e-11 (3.993) | 0.0000e+00 |
| density | 0.03125 | 3584 | 4.3388e-04 (0.000) | 2.6457e-05 (-0.002) | 1.5003e-06 (-0.037) | 2.6499e-05 (-0.002) | 9.7015e-11 (3.993) | 0.0000e+00 |
| CTT | 0.015625 | 14336 | 3.2187e-10 (3.997) | 4.9714e-11 (3.991) | 6.4205e-12 (3.971) | 5.0127e-11 (3.991) | 5.9503e-12 (3.998) | 0.0000e+00 |
| density | 0.015625 | 14336 | 4.3389e-04 (-0.000) | 2.6462e-05 (-0.000) | 1.5095e-06 (-0.009) | 2.6505e-05 (-0.000) | 6.0705e-12 (3.998) | 0.0000e+00 |
| CTT | 0.0078125 | 57344 | 6.0006e-11 (2.423) | 3.1123e-12 (3.998) | 4.0332e-13 (3.993) | 3.1383e-12 (3.998) | 3.7200e-13 (4.000) | 0.0000e+00 |
| density | 0.0078125 | 57344 | 4.3389e-04 (-0.000) | 2.6463e-05 (-0.000) | 1.5117e-06 (-0.002) | 2.6506e-05 (-0.000) | 3.7952e-13 (4.000) | 0.0000e+00 |

## Inter-hole Hamiltonian reflection diagnostic

The symmetric head-on configuration has an even Hamiltonian. The antisymmetric RMS is computed as `RMS((H(x,y) − H(−x,y))/2)` using the same grid. It measures symmetry-breaking evaluation noise; it is not subtracted from any reported norm.

| h/M | H RMS | Antisymmetric RMS |
|---:|---:|---:|
| 0.0625 | 5.0347e-08 | 3.6983e-13 |
| 0.03125 | 3.1535e-09 | 1.4875e-12 |
| 0.015625 | 1.9775e-10 | 5.9333e-12 |
| 0.0078125 | 5.6723e-11 | 2.3834e-11 |

Inter-hole H has orders 3.997 and 3.995 on h=1/16, 1/32, 1/64. At h=1/128 its RMS is 5.67e-11 and the antisymmetric RMS is 2.38e-11. The latter grows approximately fourfold per halving of h, consistent with Float64 field-evaluation noise amplified by second differences. The raw finest order 1.802 is retained; no residual is subtracted or threshold relaxed.

The collars are 0.75 ≤ R_rest ≤ 1.5 M (R_h=0.63593977642346233 M in the supplied profile). Inter-hole and outer masks and the fixed axis strips are defined in [README.md](README.md). The earlier 0.5-inner-radius collar control is archived separately. This ladder does not resolve the companion’s separate scaled D¹C interface plateau; that declared limitation is retained.

## Timings

| Mode | Refinement | Program seconds |
|---|---:|---:|
| ctt | 0.5 | 0.797 |
| ctt | 1 | 1.872 |
| ctt | 2 | 6.650 |
| ctt | 4 | 24.696 |
| density | 0.5 | 0.079 |
| density | 1 | 0.296 |
| density | 2 | 1.134 |
| density | 4 | 4.334 |

## Bit identity

Baseline and current fixture stdout compare byte-for-byte equal. The 13 fingerprints cover every computed numeric value in the reference, density-binary and B/E echo fixture tables. Frozen baseline was built from e1ee1b8 with only fingerprint printing added.

```text
Tests/EMSTrumpet/fixtures
jets.tsv all-output fingerprint=cd36dfcad48b7b46
radial.tsv all-output fingerprint=96d7d2dbbc8d287e
objects.tsv all-output fingerprint=e1107d452abc4e4b
single_ccz4.tsv all-output fingerprint=c30cf22126f14d61
binary_ccz4.tsv all-output fingerprint=3d5a69d4cda91686
Tests/EMSTrumpet/fixtures-echo/B
jets.tsv all-output fingerprint=82c55012bcb072e0
radial.tsv all-output fingerprint=cd3f68ffa3648881
objects.tsv all-output fingerprint=db6b48eb1c1fc5ad
single_ccz4.tsv all-output fingerprint=8cfc83d44aa754d9
Tests/EMSTrumpet/fixtures-echo/E
jets.tsv all-output fingerprint=71ee7188472435d3
radial.tsv all-output fingerprint=f9cd4bcbef959a70
objects.tsv all-output fingerprint=c5512734504c6e77
single_ccz4.tsv all-output fingerprint=59fb7d5a88e073fd
```
