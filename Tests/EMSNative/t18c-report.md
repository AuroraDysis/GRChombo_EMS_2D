READY-EXCEPT

Native t=0 audit and launch receipts pass the execution checks. ems_use_geometric_initial_lapse=true. Diagnostic LOW merger lapse rule: PASS. Convergence admission remains false.

LOW d16 uses h=1.75/4096, R_h/h=14.1381, box size24 and dt_multiplier=0.5. Every production cell and ghost enters the strict t0 bounds and bit comparison. The unused-array initialization fixture and translated-launch t0 checks are bit-identical. The maximum measured native-leg RSS is 2.8758 GB.

| Run | Quantity | Min | Max | d32 coarse min | d32 coarse max |
| --- | --- | --- | --- | --- | --- |
| full-sqrt | a_left | 0.000655491 | 0.999684 | 0.000190921 | 0.997955 |
| full-sqrt | a_right | 0.000655491 | 0.999684 | 0.000190921 | 0.997955 |
| full-sqrt | alpha_bin | 0.000655491 | 0.999367 | 0.000190921 | 0.995823 |
| full-sqrt | native_lapse | 0.00209741 | 0.99937 | 0.000839521 | 0.995836 |
| full-geometric | a_left | 0.000655491 | 0.999684 | 0.000190921 | 0.997955 |
| full-geometric | a_right | 0.000655491 | 0.999684 | 0.000190921 | 0.997955 |
| full-geometric | alpha_bin | 0.000655491 | 0.999367 | 0.000190921 | 0.995823 |
| full-geometric | native_lapse | 0.000655491 | 0.999367 | 0.000190921 | 0.995823 |

| Hole | R | alpha_bin/a_own d16 | d32 |
| --- | --- | --- | --- |
| left | 0.01 | 0.999781 | 0.999892 |
| left | 0.001 | 0.999998 | 0.999999 |
| left | 0.0001 | 1 | 1 |
| left | 1.0000000000000001e-05 | 1 | 1 |
| left | 9.9999999999999995e-07 | 1 | 1 |
| left | 9.9999999999999995e-08 | 1 | 1 |
| left | 1e-08 | 1 | 1 |
| left | 1.0000000000000001e-09 | 1 | 1 |
| left | 1e-10 | 1 | 1 |
| right | 0.01 | 0.999781 | 0.999892 |
| right | 0.001 | 0.999998 | 0.999999 |
| right | 0.0001 | 1 | 1 |
| right | 1.0000000000000001e-05 | 1 | 1 |
| right | 9.9999999999999995e-07 | 1 | 1 |
| right | 9.9999999999999995e-08 | 1 | 1 |
| right | 1e-08 | 1 | 1 |
| right | 1.0000000000000001e-09 | 1 | 1 |
| right | 1e-10 | 1 | 1 |

| Run | Hole | Source | Peak | RMS | d32 peak | d32 RMS |
| --- | --- | --- | --- | --- | --- | --- |
| full-sqrt | left | physical_Gamma_n | 357.189 | 203.391 | 379.531 | 293.187 |
| full-sqrt | left | KO_Gamma_n | 0.786041 | 0.420177 | 0.27212 | 0.209187 |
| full-sqrt | left | total_Gamma_n | 357.166 | 203.472 | 379.803 | 293.187 |
| full-sqrt | right | physical_Gamma_n | 357.189 | 203.391 | 379.531 | 293.187 |
| full-sqrt | right | KO_Gamma_n | 0.786041 | 0.420177 | 0.27212 | 0.209187 |
| full-sqrt | right | total_Gamma_n | 357.166 | 203.472 | 379.803 | 293.187 |
| full-geometric | left | physical_Gamma_n | 5.07636 | 4.11798 | 5.81954 | 4.57387 |
| full-geometric | left | KO_Gamma_n | 0.786041 | 0.420177 | 0.27212 | 0.209187 |
| full-geometric | left | total_Gamma_n | 5.8624 | 4.38819 | 6.09166 | 4.57912 |
| full-geometric | right | physical_Gamma_n | 5.07636 | 4.11798 | 5.81954 | 4.57387 |
| full-geometric | right | KO_Gamma_n | 0.786041 | 0.420177 | 0.27212 | 0.209187 |
| full-geometric | right | total_Gamma_n | 5.8624 | 4.38819 | 6.09166 | 4.57912 |

| Region | Constraint | d16 peak | d16 RMS | d32 peak | d32 RMS |
| --- | --- | --- | --- | --- | --- |
| far | Ham | 5.79894e-13 | 9.42371e-14 | 5.60165e-12 | 1.4841e-13 |
| far | Mom | 9.19347e-15 | 6.54199e-16 | 2.51469e-13 | 1.86996e-14 |
| far | GaussE | 1.34628e-14 | 1.75875e-15 | 6.96822e-13 | 4.73864e-14 |
| far | GaussB | 0 | 0 | 0 | 0 |
| far | C_Gamma | 0 | 0 | 0 | 0 |
| elsewhere-off-sheet | Ham | 0.000553099 | 1.31994e-10 | 1.19053e-05 | 2.2509e-11 |
| elsewhere-off-sheet | Mom | 0.00475509 | 2.12114e-09 | 5.22249e-05 | 3.96054e-11 |
| elsewhere-off-sheet | GaussE | 8.47324e-05 | 3.25805e-11 | 1.60815e-06 | 3.75366e-12 |
| elsewhere-off-sheet | GaussB | 0 | 0 | 0 | 0 |
| elsewhere-off-sheet | C_Gamma | 0 | 0 | 0 | 0 |
| interhole | Ham | 9.97801 | 4.32177e-06 | 1.00906e-09 | 9.82543e-11 |
| interhole | Mom | 49365.4 | 0.0205922 | 8.9205e-12 | 1.64887e-12 |
| interhole | GaussE | 363.705 | 0.000265001 | 1.73536e-10 | 1.34428e-11 |
| interhole | GaussB | 0 | 0 | 0 | 0 |
| interhole | C_Gamma | 0 | 0 | 0 | 0 |
| collar-left | Ham | 0.000950231 | 8.38598e-05 | 1.19053e-05 | 1.4573e-06 |
| collar-left | Mom | 0.0100997 | 0.00184096 | 6.65461e-05 | 1.65253e-05 |
| collar-left | GaussE | 0.000148074 | 2.1968e-05 | 1.60815e-06 | 3.7578e-07 |
| collar-left | GaussB | 0 | 0 | 0 | 0 |
| collar-left | C_Gamma | 0 | 0 | 0 | 0 |
| collar-right | Ham | 0.000950231 | 8.38598e-05 | 1.19053e-05 | 1.4573e-06 |
| collar-right | Mom | 0.0100997 | 0.00184096 | 6.65461e-05 | 1.65253e-05 |
| collar-right | GaussE | 0.000148074 | 2.1968e-05 | 1.60815e-06 | 3.7578e-07 |
| collar-right | GaussB | 0 | 0 | 0 | 0 |
| collar-right | C_Gamma | 0 | 0 | 0 | 0 |

Full sheet and AMR-face-excluded norms and panel crossing counts are retained in [interfaces](t18c-interfaces.csv). Sheet masks retain the same selector/stencil definition; panel surfaces move with the d16 companion, and the LOW reference spacing differs. These are not a fixed-separation convergence pair.

| Hole | Ray | Field | Norm | Ratio | five-spread interval | signal margins G / S | reduction margin | d32 endpoint ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| left | axis | Gamma | peak | 0.0886544 | 0.0865802 / 0.0907378 | 47.146 / 446.985 | 221.327 | 0.044044 |
| left | axis | Gamma | RMS | 0.0961972 | 0.0951469 / 0.0972506 | 106.002 / 666.379 | 375.31 | 0.0478794 |
| left | axis | metric_Gamma | peak | 0.401311 | 0.376539 / 0.427489 | 27.9115 / 36.2271 | 14.2607 | 0.188259 |
| left | axis | metric_Gamma | RMS | 0.32322 | 0.30394 / 0.343425 | 26.5842 / 42.6797 | 19.0167 | 0.0755707 |
| left | axis | shift | peak | 0.340177 | 0.337898 / 0.342466 | 215.959 / 479.929 | 180.337 | 0.351281 |
| left | axis | shift | RMS | 0.234703 | 0.233712 / 0.235695 | 311.279 / 989.011 | 433.571 | 0.293115 |
| left | axis | lapse | peak | 0.589841 | 0.53274 / 0.654962 | 26.6397 / 15.2389 | 4.67348 | 0.765521 |
| left | axis | lapse | RMS | 0.781286 | 0.746196 / 0.818742 | 72.6239 / 30.667 | 5.04341 | 0.794551 |
| left | diagonal | Gamma | peak | 0.0930309 | 0.0922778 / 0.0937862 | 151.74 / 659.362 | 425.865 | 0.0475481 |
| left | diagonal | Gamma | RMS | 0.0955331 | 0.0951035 / 0.0959636 | 295.886 / 891.455 | 626.088 | 0.0514909 |
| left | diagonal | metric_Gamma | peak | 0.374589 | 0.332841 / 0.421158 | 15.8852 / 18.3206 | 8.00124 | 0.0475441 |
| left | diagonal | metric_Gamma | RMS | 0.177966 | 0.159395 / 0.197753 | 13.1667 / 31.5353 | 18.1758 | 0.0519498 |
| left | diagonal | shift | peak | 0.209347 | 0.208593 / 0.210102 | 336.71 / 1575.94 | 629.357 | 0.222369 |
| left | diagonal | shift | RMS | 0.12769 | 0.127363 / 0.128018 | 439.754 / 3434.63 | 1500.05 | 0.174619 |
| left | diagonal | lapse | peak | 0.585755 | 0.559611 / 0.613398 | 55.5109 / 35.8908 | 10.7836 | 0.773988 |
| left | diagonal | lapse | RMS | 0.738076 | 0.72726 / 0.749099 | 191.345 / 104.526 | 19.5112 | 0.797868 |
| right | axis | Gamma | peak | 0.0886544 | 0.0865802 / 0.0907378 | 47.146 / 446.985 | 221.327 | 0.044044 |
| right | axis | Gamma | RMS | 0.0961972 | 0.0951469 / 0.0972506 | 106.002 / 666.379 | 375.31 | 0.0478794 |
| right | axis | metric_Gamma | peak | 0.401311 | 0.376539 / 0.427489 | 27.9115 / 36.2271 | 14.2607 | 0.188259 |
| right | axis | metric_Gamma | RMS | 0.32322 | 0.30394 / 0.343425 | 26.5842 / 42.6797 | 19.0167 | 0.0755707 |
| right | axis | shift | peak | 0.340177 | 0.337898 / 0.342466 | 215.959 / 479.929 | 180.337 | 0.351281 |
| right | axis | shift | RMS | 0.234703 | 0.233712 / 0.235695 | 311.279 / 989.011 | 433.571 | 0.293115 |
| right | axis | lapse | peak | 0.589841 | 0.53274 / 0.654962 | 26.6397 / 15.2389 | 4.67348 | 0.765521 |
| right | axis | lapse | RMS | 0.781286 | 0.746196 / 0.818742 | 72.6239 / 30.667 | 5.04341 | 0.794551 |
| right | diagonal | Gamma | peak | 0.0930309 | 0.0922778 / 0.0937862 | 151.74 / 659.362 | 425.865 | 0.0475481 |
| right | diagonal | Gamma | RMS | 0.0955331 | 0.0951035 / 0.0959636 | 295.886 / 891.455 | 626.088 | 0.0514909 |
| right | diagonal | metric_Gamma | peak | 0.374589 | 0.332841 / 0.421158 | 15.8852 / 18.3206 | 8.00124 | 0.0475441 |
| right | diagonal | metric_Gamma | RMS | 0.177966 | 0.159395 / 0.197753 | 13.1667 / 31.5353 | 18.1758 | 0.0519498 |
| right | diagonal | shift | peak | 0.209347 | 0.208593 / 0.210102 | 336.71 / 1575.94 | 629.357 | 0.222369 |
| right | diagonal | shift | RMS | 0.12769 | 0.127363 / 0.128018 | 439.754 / 3434.63 | 1500.05 | 0.174619 |
| right | diagonal | lapse | peak | 0.585755 | 0.559611 / 0.613398 | 55.5109 / 35.8908 | 10.7836 | 0.773988 |
| right | diagonal | lapse | RMS | 0.738076 | 0.72726 / 0.749099 | 191.345 / 104.526 | 19.5112 | 0.797868 |

All nine positive native clocks are retained in [launch ratios](t18c-launch-ratios.csv); early enhancements and unresolved entries remain in [history summary](t18c-history-summary.csv). Margins are T16b five-spread sampling sensitivity estimates, not rigorous error bounds. Native four-cell changes and KO peaks are in [puncture corroboration](t18c-native-puncture.csv).

The d16 endpoint is 0.001922607421875 M_i, 3.077% before d32 T16b 0.001983642578125. d32 T16b has R_h/h49.483/98.967 and dt=h/4; separation, grid phase, spatial resolution, dt and endpoint all differ. The table compares the actual endpoints without time alignment; it cannot attribute differences solely to d. No second-rung or same-grid temporal control was requested here.

The native constraints are identical between lapse variants. The documented d16 D2(logpsi) jump is 2.8922936531e-5 versus d32 about8.89e-6; both exceed the unchanged1e-9 gate. This initial-lapse choice cannot repair that defect. The requested LOW diagnostic lapse criterion passes at both holes; the companion may enter an explicitly diagnostic low-resolution merger with geometric lapse. The strict T16b two-rung/convergence criteria are not met.

Reproduce with t18c-work.py generate/census/prepare/build/smoke, then the measured detached pipeline plan, then t18c-analyze.py. Every native child has its own receipt; zero done.exit is execution evidence only. No SSH or commit.

The unchanged inter-hole box |x-midpoint|≤8 touches both punctures when d=16; at d=32 it excludes both. Its large d16 peaks therefore include the unresolved puncture neighborhood. The native extrema lie at r_near=0.000215792 M_i (Ham), 0.000215792 M_i (Mom), 0.000641595 M_i (GaussE). The mask is retained exactly and no near-hole samples are removed. These peaks are not a like-for-like exterior separation comparison. [Locations](t18c-constraint-hotspots.csv).

The selector census observes 24/44 declared geometric interfaces at native stencil spacing, versus 30/44 in T16b d32. Unobserved end sheets remain explicitly labelled without an inferred pass or order. [Coverage](t18c-interface-coverage.csv).

The cropped launch's entire refined t0 non-lapse state and base ±200-M_i neighborhood match 53,342,145 full-production values with zero bit mismatches. All 28 captured finest t0 fields also match. The recorded update counts give a deliberately generous coordinate-footprint estimate of 168.417 M_i, versus 440 M_i to the nearer cropped outer boundary; this is a scope estimate, not a full-production evolved-state bit comparison. [Domain check](t18c-domain-check.json).

The final analysis and its enclosing queue each have their own successful actual-child receipt; 640 cached launch amplitudes/spreads and every one of 288 five-spread ratio intervals reproduce exactly. Failed preliminary reporting receipts are retained. The authoritative completion marker is `/private/tmp/ems-t18c/completion-verified/done.exit = 0`. No native job remains pending. [Completion checks](t18c-completion-checks.csv).

Endpoint margins compare actual endpoints and include the minimum across peak/RMS; the d32 coarse and fine amplitudes, intervals and margins are also retained in [the endpoint table](t18c-launch-endpoint.csv).

| Hole | Ray | Field | d16 signal / 5 spread | d32 coarse signal / 5 spread | d16 reduction margin | d32 coarse reduction margin |
| --- | --- | --- | --- | --- | --- | --- |
| left | axis | Gamma | 47.146 | 8075.32 | 221.327 | 19941.8 |
| left | axis | metric_Gamma | 26.5842 | 55.8993 | 14.2607 | 68.1029 |
| left | axis | shift | 215.959 | 33730 | 180.337 | 19572.6 |
| left | axis | lapse | 15.2389 | 1529.76 | 4.67348 | 282.235 |
| left | diagonal | Gamma | 151.74 | 24202.9 | 425.865 | 22388.6 |
| left | diagonal | metric_Gamma | 13.1667 | 82.7623 | 8.00124 | 438.574 |
| left | diagonal | shift | 336.71 | 204404 | 629.357 | 192994 |
| left | diagonal | lapse | 35.8908 | 9311.95 | 10.7836 | 1654.98 |
| right | axis | Gamma | 47.146 | 8075.32 | 221.327 | 19941.8 |
| right | axis | metric_Gamma | 26.5842 | 55.8993 | 14.2607 | 68.1029 |
| right | axis | shift | 215.959 | 33730 | 180.337 | 19572.6 |
| right | axis | lapse | 15.2389 | 1529.76 | 4.67348 | 282.235 |
| right | diagonal | Gamma | 151.74 | 24202.9 | 425.865 | 22388.6 |
| right | diagonal | metric_Gamma | 13.1667 | 82.7623 | 8.00124 | 438.574 |
| right | diagonal | shift | 336.71 | 204404 | 629.357 | 192994 |
| right | diagonal | lapse | 35.8908 | 9311.95 | 10.7836 | 1654.98 |
