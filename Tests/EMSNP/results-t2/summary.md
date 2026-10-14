# NPTET T2 measured controls

**Not ready for a simulation:** (b)(i) fails and two (c) modes exceed their declared error estimates.

## Static metric/potential oracles

| case | coordinate r | areal R | exact Psi2 | exact native Phi1 |
|---|---:|---:|---:|---:|
| rn | 1 | 2.10062257 | -0.0827183136 | 0.0158216605 |
| rn | 2 | 3.05812987 | -0.0293625172 | 0.00746511094 |
| rn | 4 | 5.03091385 | -0.0070885228 | 0.00275838164 |
| ems | 0.0512 | 0.385750651 | -1.18753759 | 0.0349380132 |
| ems | 0.2622 | 0.96496845 | -0.157957718 | 0.0515062567 |
| ems | 0.6069 | 1.46646698 | -0.089369962 | 0.0386315184 |
| ems | 1.0429 | 1.96287774 | -0.0558914938 | 0.0250215527 |

## Fine-grid static residuals (70 scalar/sphere checks pass)

Psi2 and Phi1 cells show maximum absolute error / total tolerance. The other columns normalize maximum absolute residuals by their Coulomb field scale.

| case | r | dx | Psi2 error / tol | Phi1 error / tol | max radiative Psi / abs(Psi2) | max Phi0,2 / abs(Phi1) | max X / abs(Phi1 Psi2) |
|---|---:|---:|---:|---:|---:|---:|---:|
| rn | 1 | 0.03125 | 1.965e-07 / 1.458e-06 | 2.501e-08 / 2.002e-07 | 6.771e-06 | 1.281e-06 | 7.736e-06 |
| rn | 2 | 0.03125 | 9.663e-09 / 7.350e-08 | 9.456e-10 / 7.569e-09 | 5.085e-07 | 6.314e-08 | 8.276e-07 |
| rn | 4 | 0.03125 | 2.388e-10 / 1.875e-09 | 2.755e-11 / 2.204e-10 | 2.998e-08 | 6.698e-09 | 5.113e-08 |
| ems | 0.0512 | 0.00195312 | 1.488e-05 / 1.281e-04 | 1.943e-08 / 1.353e-07 | 1.428e-05 | 7.491e-07 | 2.718e-05 |
| ems | 0.2622 | 0.00195312 | 2.296e-09 / 2.593e-08 | 7.658e-11 / 6.096e-10 | 1.089e-08 | 1.544e-09 | 2.052e-08 |
| ems | 0.6069 | 0.00195312 | 3.759e-11 / 7.706e-09 | 3.700e-12 / 3.352e-11 | 4.134e-10 | 1.111e-10 | 5.231e-10 |
| ems | 1.0429 | 0.00195312 | 7.236e-11 / 7.496e-09 | 3.951e-13 / 3.169e-12 | 1.294e-09 | 1.631e-11 | 6.872e-10 |

## Matter checkpoint suppression: RMS of pointwise ratios

Both classes and both combinations are retained. k is the second null vector (called n in the task). Opposite-class denominators have no Coulomb leading term, so their O(1) ratios do not test the protected cancellation.

| r | areal R | Phi2/Phi1 | Psi3/Psi2 | l-fixed Xout | l-fixed Xin | k-fixed Xout | k-fixed Xin |
|---|---:|---:|---:|---:|---:|---:|---:|
| 50 | 51.9708451 | 3.293e-03 | 7.529e-03 | 4.160e-07 | 1.057e+00 | 9.995e-01 | 1.140e-06 |
| 75 | 76.9786032 | 4.127e-03 | 9.485e-03 | 3.079e-06 | 1.041e+00 | 9.996e-01 | 2.403e-06 |
| 100 | 101.986174 | 2.934e-03 | 6.558e-03 | 1.441e-06 | 1.022e+00 | 9.999e-01 | 1.506e-06 |
| 0.0512 | 0.519099906 | 2.164e-03 | 2.748e-03 | 3.271e-06 | 9.887e-01 | 1.009e+00 | 2.209e-06 |
| 0.2622 | 1.33301019 | 2.667e-03 | 9.111e-03 | 7.535e-06 | 9.898e-01 | 1.012e+00 | 1.377e-06 |
| 0.6069 | 2.07643525 | 2.007e-03 | 8.157e-03 | 3.180e-06 | 9.919e-01 | 1.965e+00 | 6.460e-07 |
| 1.0429 | 2.7160897 | 1.448e-03 | 5.428e-03 | 9.060e-07 | 9.920e-01 | 1.048e+00 | 7.020e-07 |

## Pointwise derivative roundoff gate (b)(i): FAIL

| r | maximum FD/CAS absolute error | maximum error / budget | failed comparisons / 380 | coordinate Weyl trace max |
|---|---:|---:|---:|---:|
| 50 | 2.623e-14 | 1.172e+04 | 380 | 3.247e-10 |
| 75 | 8.770e-15 | 1.998e+04 | 380 | 2.439e-10 |
| 100 | 1.736e-15 | 1.520e+04 | 380 | 8.462e-11 |
| 0.0512 | 3.489e-07 | 7.602e+00 | 380 | 2.164e-04 |
| 0.2622 | 1.261e-07 | 2.723e+02 | 364 | 2.499e-05 |
| 0.6069 | 5.105e-08 | 9.381e+02 | 380 | 6.325e-06 |
| 1.0429 | 2.080e-08 | 1.569e+03 | 380 | 2.177e-06 |

## Weyl4: l=2 matter modes pass; 28/30 total modes pass

| r | R Psi4 | old W | absolute error | tolerance | relative error |
|---|---:|---:|---:|---:|---:|
| 50 | -0.000287687025258 | -0.000287686615813 | 4.094e-10 | 4.641e-09 | 1.423e-06 |
| 75 | -0.000372537787728 | -0.000372537507288 | 2.804e-10 | 5.303e-09 | 7.528e-07 |
| 100 | -0.000430712881552 | -0.000430712587229 | 2.943e-10 | 7.292e-09 | 6.833e-07 |

Failed weak modes:

| form | r | l | absolute error | tolerance |
|---|---:|---:|---:|---:|
| vacuum | 100 | 3 | 3.343e-14 | 7.393e-15 |
| matter | 100 | 3 | 3.343e-14 | 7.392e-15 |

## Output checks

- Independent harmonic re-projection: 980 components pass; worst residual/budget 0.001805. Both pointwise X expressions pass.
- Default off and NP enabled: 42 old radiation rows unchanged within 1e-12 relative; maximum difference 0. No new file with NP off.
- Exact CAS witnesses and production syntax check pass. No production link/evolution or continuum certification is claimed.
- Complete budgets, observed convergence orders, incoming ratios, both Weyl forms and all modes are in validation.json.
