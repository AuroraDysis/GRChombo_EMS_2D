# NPTET T3 measured controls

**Ready for simulation under every requested gate: False**

## Removed components and same-point constraints

All entries are sphere RMS except the last three maximum residual columns. E trace uses the complete spatial metric. The production Hamiltonian uses independently evolved Gamma; H(metric) reconstructs Gamma and its derivative from the same metric jets. The latter agrees with the removed E trace after the measured K tr(A)/3 term. Magnetic trace is zero in this polar sector; its removed antisymmetric axial vector tracks minus the momentum constraint. The small tabulated axial residual includes derivatives of residual tr(A), which the production momentum formula assumes zero.

| r | trace E | H(production) | H(metric) | axial B | norm M | max E-Hmetric-KtrA/3 | max axialB+M | max projected C trace |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 1.057e-10 | 4.943e-10 | 1.057e-10 | 2.669e-10 | 2.669e-10 | 3.566e-17 | 5.021e-15 | 3.522e-20 |
| 75 | 1.807e-10 | 3.366e-10 | 1.807e-10 | 5.825e-11 | 5.825e-11 | 5.622e-17 | 1.250e-15 | 1.244e-20 |
| 100 | 2.703e-11 | 7.084e-11 | 2.703e-11 | 2.795e-11 | 2.795e-11 | 2.931e-17 | 3.356e-15 | 3.569e-21 |
| 0.0512 | 7.689e-06 | 1.391e-05 | 7.689e-06 | 2.112e-05 | 2.112e-05 | 1.228e-12 | 1.751e-10 | 2.285e-13 |
| 0.2622 | 6.305e-07 | 2.365e-06 | 6.305e-07 | 4.844e-06 | 4.844e-06 | 2.508e-13 | 2.334e-11 | 4.623e-15 |
| 0.6069 | 3.556e-07 | 9.981e-07 | 3.556e-07 | 1.804e-06 | 1.804e-06 | 5.982e-14 | 1.086e-12 | 1.028e-15 |
| 1.0429 | 3.252e-07 | 4.329e-07 | 3.252e-07 | 8.091e-07 | 8.091e-07 | 6.094e-14 | 2.687e-13 | 3.402e-16 |

## Both null-rotation classes: roundoff comparison

| r | max FD/CAS error | max error/budget | failed comparisons / 380 |
|---|---:|---:|---:|
| 50 | 1.903e-23 | 5.290e-06 | 0 |
| 75 | 2.864e-24 | 5.894e-06 | 0 |
| 100 | 7.238e-25 | 6.338e-06 | 0 |
| 0.0512 | 1.657e-16 | 3.495e-09 | 0 |
| 0.2622 | 2.682e-17 | 5.667e-08 | 0 |
| 0.6069 | 1.382e-17 | 2.523e-07 | 0 |
| 1.0429 | 8.964e-18 | 6.728e-07 | 0 |

## Protected suppression: RMS of pointwise ratios

X_out protects l fixed; X_in protects k fixed. Perturbation ratios are outgoing and incoming respectively. Opposite-class ratios remain in validation.json and have no Coulomb denominator term.

| r | areal R | Phi2/Phi1 | Psi3/Psi2 | l-fixed Xout | Phi0/Phi1 | Psi1/Psi2 | k-fixed Xin |
|---|---:|---:|---:|---:|---:|---:|---:|
| 50 | 51.9708451 | 3.293e-03 | 7.529e-03 | 7.243e-07 | 4.676e-04 | 7.208e-04 | 1.633e-06 |
| 75 | 76.9786032 | 4.127e-03 | 9.485e-03 | 6.797e-07 | 3.182e-04 | 5.176e-04 | 1.422e-06 |
| 100 | 101.986174 | 2.934e-03 | 6.558e-03 | 6.387e-07 | 1.530e-04 | 4.244e-04 | 4.987e-07 |
| 0.0512 | 0.519099906 | 2.164e-03 | 2.748e-03 | 2.620e-06 | 1.665e-03 | 2.738e-03 | 1.548e-06 |
| 0.2622 | 1.33301019 | 2.667e-03 | 9.111e-03 | 7.945e-06 | 9.049e-04 | 6.016e-03 | 1.719e-06 |
| 0.6069 | 2.07643525 | 2.007e-03 | 8.157e-03 | 3.716e-06 | 4.923e-04 | 3.775e-03 | 8.824e-07 |
| 1.0429 | 2.7160897 | 1.448e-03 | 5.428e-03 | 1.471e-06 | 3.989e-04 | 2.205e-03 | 6.549e-07 |

## Production Weyl4: requested pure sphere-scale tolerance

epsilon=2.2e-16, Ntheta=97, C=4096; tolerance=C epsilon Ntheta max_l abs(mode). C is 16 times the 256-term contraction count. No spatial-error allowance is included. Identical-jet Psi0/Psi4 preservation passes; the existing cell-then-interpolate production W differs from the NP interpolate-then-reconstruct route.

| form | r | l | NP mode | production W | sphere mode scale | error | tolerance | result |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| vacuum | 50 | 2 | -0.000287686955836 | -0.000287686546392 | 2.877e-04 | 4.094e-10 | 2.515e-14 | FAIL |
| vacuum | 50 | 3 | 1.7933057251e-11 | 1.79211568687e-11 | 2.877e-04 | 1.190e-14 | 2.515e-14 | PASS |
| vacuum | 50 | 4 | 4.66860295067e-07 | 4.66893348685e-07 | 2.877e-04 | 3.305e-11 | 2.515e-14 | FAIL |
| vacuum | 50 | 5 | 7.35352300946e-14 | 7.84648669772e-14 | 2.877e-04 | 4.930e-15 | 2.515e-14 | PASS |
| vacuum | 50 | 6 | -1.08269227071e-07 | -1.0823253523e-07 | 2.877e-04 | 3.669e-11 | 2.515e-14 | FAIL |
| vacuum | 75 | 2 | -0.000372537769251 | -0.00037253748881 | 3.725e-04 | 2.804e-10 | 3.256e-14 | FAIL |
| vacuum | 75 | 3 | -1.4539487033e-11 | -1.45336487136e-11 | 3.725e-04 | 5.838e-15 | 3.256e-14 | PASS |
| vacuum | 75 | 4 | 4.38693552019e-06 | 4.38674393281e-06 | 3.725e-04 | 1.916e-10 | 3.256e-14 | FAIL |
| vacuum | 75 | 5 | -8.96478137021e-14 | -8.90470567286e-14 | 3.725e-04 | 6.008e-16 | 3.256e-14 | PASS |
| vacuum | 75 | 6 | 5.15400273942e-07 | 5.15602467371e-07 | 3.725e-04 | 2.022e-10 | 3.256e-14 | FAIL |
| vacuum | 100 | 2 | -0.000430712878531 | -0.000430712584208 | 4.307e-04 | 2.943e-10 | 3.765e-14 | FAIL |
| vacuum | 100 | 3 | 7.0429547687e-12 | 7.00952437599e-12 | 4.307e-04 | 3.343e-14 | 3.765e-14 | PASS |
| vacuum | 100 | 4 | 3.02034391327e-06 | 3.02029566287e-06 | 4.307e-04 | 4.825e-11 | 3.765e-14 | FAIL |
| vacuum | 100 | 5 | -7.56396491209e-15 | -1.25571678957e-14 | 4.307e-04 | 4.993e-15 | 3.765e-14 | PASS |
| vacuum | 100 | 6 | 4.24488174033e-08 | 4.244021191e-08 | 4.307e-04 | 8.605e-12 | 3.765e-14 | FAIL |
| matter | 50 | 2 | -0.000287687025258 | -0.000287686615813 | 2.877e-04 | 4.094e-10 | 2.515e-14 | FAIL |
| matter | 50 | 3 | 1.79330551671e-11 | 1.7921154721e-11 | 2.877e-04 | 1.190e-14 | 2.515e-14 | PASS |
| matter | 50 | 4 | 4.66779909381e-07 | 4.6681296293e-07 | 2.877e-04 | 3.305e-11 | 2.515e-14 | FAIL |
| matter | 50 | 5 | 7.35367996082e-14 | 7.84664623595e-14 | 2.877e-04 | 4.930e-15 | 2.515e-14 | PASS |
| matter | 50 | 6 | -1.08269928702e-07 | -1.08233237816e-07 | 2.877e-04 | 3.669e-11 | 2.515e-14 | FAIL |
| matter | 75 | 2 | -0.000372537787728 | -0.000372537507288 | 3.725e-04 | 2.804e-10 | 3.256e-14 | FAIL |
| matter | 75 | 3 | -1.4539488328e-11 | -1.45336500547e-11 | 3.725e-04 | 5.838e-15 | 3.256e-14 | PASS |
| matter | 75 | 4 | 4.38691437106e-06 | 4.38672278367e-06 | 3.725e-04 | 1.916e-10 | 3.256e-14 | FAIL |
| matter | 75 | 5 | -8.96491561365e-14 | -8.90483469154e-14 | 3.725e-04 | 6.008e-16 | 3.256e-14 | PASS |
| matter | 75 | 6 | 5.15400047901e-07 | 5.15602241389e-07 | 3.725e-04 | 2.022e-10 | 3.256e-14 | FAIL |
| matter | 100 | 2 | -0.000430712881552 | -0.000430712587229 | 4.307e-04 | 2.943e-10 | 3.765e-14 | FAIL |
| matter | 100 | 3 | 7.04295521035e-12 | 7.00952477015e-12 | 4.307e-04 | 3.343e-14 | 3.765e-14 | PASS |
| matter | 100 | 4 | 3.02034044978e-06 | 3.02029219936e-06 | 4.307e-04 | 4.825e-11 | 3.765e-14 | FAIL |
| matter | 100 | 5 | -7.56346959302e-15 | -1.25567277402e-14 | 4.307e-04 | 4.993e-15 | 3.765e-14 | PASS |
| matter | 100 | 6 | 4.24489840551e-08 | 4.24403785241e-08 | 4.307e-04 | 8.606e-12 | 3.765e-14 | FAIL |

The l=3 values are symmetry-suppressed, but are not zero to machine roundoff: NP/scale and W/scale are retained in validation.json. The gate above compares the two measured modes, without replacing either by zero.

## Static backgrounds

All ten scalars on seven spheres are checked at both original T2 resolutions. Budgets are unchanged: 8/15 of the coarse/fine difference plus conditioning allowances (128 eps/dx² for curvature). Exact metric/potential oracle columns and all vanishing-field budgets are in validation.json.

| case | r | exact Psi2 | fine error / tolerance | exact Phi1 | fine error / tolerance |
|---|---:|---:|---:|---:|---:|
| rn | 1 | -0.0827183136 | 1.143e-07 / 8.380e-07 | 0.0158216605 | 2.501e-08 / 2.002e-07 |
| rn | 2 | -0.0293625172 | 6.508e-09 / 5.623e-08 | 0.00746511094 | 9.456e-10 / 7.569e-09 |
| rn | 4 | -0.0070885228 | 1.764e-10 / 1.288e-09 | 0.00275838164 | 2.755e-11 / 2.204e-10 |
| ems | 0.0512 | -1.18753759 | 1.261e-05 / 7.686e-05 | 0.0349380132 | 1.943e-08 / 1.353e-07 |
| ems | 0.2622 | -0.157957718 | 1.654e-09 / 2.394e-08 | 0.0515062567 | 7.658e-11 / 6.096e-10 |
| ems | 0.6069 | -0.089369962 | 2.969e-11 / 7.621e-09 | 0.0386315184 | 3.700e-12 / 3.352e-11 |
| ems | 1.0429 | -0.0558914938 | 2.512e-11 / 7.475e-09 | 0.0250215527 | 3.951e-13 / 3.169e-12 |

## Other controls

- Static: Counter({'PASS': 70})
- Default: {'status': 'PASS', 'rows': 42, 'relative_tolerance': 1e-12, 'absolute_tolerance': 1e-100, 'max_relative_difference': 0, 'new_files_with_switch_off': False}
- Independent projection: {'status': 'PASS', 'mode_components': 980, 'max_error_over_roundoff_budget': 0.0018045267922320635, 'pointwise_X_definitions': 'PASS', 'independent_harmonics': 'associated Legendre factorial polynomial'}
- Exact E/B, rotation and static CAS witnesses, complex Maxwell/flat/floor unit controls and prescribed production syntax check are retained in the logs.
- Full production link/evolution was not run. No HPC, SSH, transfer or push.
