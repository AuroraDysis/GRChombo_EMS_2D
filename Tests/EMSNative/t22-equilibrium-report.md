# T22 measured initialization: the exact-zero K premise is false

The E-high initialization audit **passes with explicitly measured K roundoff**, not exact-zero K. The consult's premise that the setter supplies exact zeros for K and Theta was false for K. Theta is exactly zero, and the new driver B is exactly zero on all valid cells and stored ghosts. K is nonzero in every valid cell of both gauges, not only at the puncture or in ghosts. It peaks at 4.7423391661e-15 M^-1 on level 14 at (x,y)=(-1.7801920575e-04,1.9582112630e-03) M, R=1.9662863953e-03 M. Its maximum scaled residual is 3.72038377 eps times the local absolute tensor-contraction scale. This is Float64 contraction/cancellation roundoff of the analytically traceless isolated, unboosted construction, not a physical nonzero mean curvature.

The initial setter builds K_ij=P^2 k(3 n_i n_j-delta_ij) at [EMSBH_trumpet_read.impl.hpp](../../Source/InitialConditions/EMSBH/EMSBH_trumpet_read.impl.hpp#L223). At zero rapidity the spatial metric is P^2 delta_ij, giving trace 3k(n.n-1)=0 for a unit direction. The exact rational witness and independent 50-digit numerical inversion are in [t22-trace-qualification.json](t22-trace-qualification.json): exact residual zero and worst scaled numerical residual 1.20490e-50. These are algebraic witnesses of the construction, not static-reference comparisons or positive-time reads. The setter initializes a trace accumulator at line 419, **computes** the metric contraction at lines 420–422, and stores it as K at line 426. It does not assign K as a literal zero. VarsTools::assign(vars,0) at line 424 supplies Theta's zero, which is never overwritten.

The original qualification assertion required every K/Theta cell to be exactly zero. Its exit-1 receipt is retained at `/private/tmp/ems-t22/evidence-pipeline/` and the original qualification log. The controller-authorized replacement explicitly requires Theta=0 and |K|<=64 eps S locally, where S=|h^11 A11|+2|h^12 A12|+|h^22 A22|+|Aww/hww| from the stored numerical state. K itself is never rounded, clipped, subtracted or overwritten. The 64-eps gate was stated before measuring the full normalized census; the measured maximum is 3.720384. Ghost copies are checked and labelled separately. This replacement concerns initialization qualification alone: the registered 10.5 M activity reading, masks, clocks, fivefold significance, temporal rule and budgets remain frozen.

## Every-level valid-cell K and Theta

| level | valid cells | K peak | K cell RMS | Theta peak / RMS | K / (eps scale) peak |
| --- | --- | --- | --- | --- | --- |
| 0 | 1492992 | 7.44007455e-19 | 7.29856170e-22 | 0 / 0 | 3.588787 |
| 1 | 165888 | 4.58991309e-18 | 1.22364578e-20 | 0 / 0 | 3.316082 |
| 2 | 165888 | 1.48961931e-17 | 5.21783120e-20 | 0 / 0 | 3.507084 |
| 3 | 165888 | 3.87186562e-17 | 1.25189258e-19 | 0 / 0 | 3.458140 |
| 4 | 165888 | 4.58416321e-16 | 1.15282423e-18 | 0 / 0 | 3.720384 |
| 5 | 165888 | 4.95453272e-16 | 1.59052959e-18 | 0 / 0 | 3.594048 |
| 6 | 165888 | 8.16610775e-16 | 3.29248726e-18 | 0 / 0 | 3.481005 |
| 7 | 165888 | 3.21814866e-15 | 1.01633111e-17 | 0 / 0 | 3.453051 |
| 8 | 165888 | 1.70165181e-15 | 1.52110726e-17 | 0 / 0 | 3.422600 |
| 9 | 165888 | 2.29871764e-15 | 3.22203217e-17 | 0 / 0 | 3.304040 |
| 10 | 165888 | 3.38582732e-15 | 6.25996156e-17 | 0 / 0 | 3.382850 |
| 11 | 165888 | 3.47584513e-15 | 1.29869847e-16 | 0 / 0 | 3.539205 |
| 12 | 165888 | 3.88968642e-15 | 2.51681471e-16 | 0 / 0 | 3.492669 |
| 13 | 165888 | 4.06515446e-15 | 4.69966085e-16 | 0 / 0 | 3.583897 |
| 14 | 165888 | 4.74233917e-15 | 7.72870859e-16 | 0 / 0 | 3.470597 |

The table uses unweighted native-cell RMS. [t22-KTheta-census.csv](t22-KTheta-census.csv) additionally reports cartoon coordinate-volume RMS (weights |y|), every level's peak location and radius, both gauge labels, the separate puncture-neighbourhood census (four rows by eight columns around the puncture), and stored ghost-copy counts. The two initial physical states agree on **113,549,904 Float64 values**, zero bit mismatches, excluding only driver B1/B2. Both hierarchies have the exact collected E-high box lists: levels 0–14, 546 boxes, 3,815,424 valid cells. The initial states and all sampled RHS pieces are finite.

## Every-level measured lapse rates

The new gauge's all-valid-cell pre-KO lapse rate is the native class's -2 alpha(K-2Theta), not an assumed zero. Its global peak is 1.5023657771e-16 M^-1; the formula check has zero residual. The old gauge's rate is native upwind beta.grad(alpha)-1.8 alpha(K-2Theta), with its historical coefficient unchanged. The new pre-KO shift rate is mu B=0 exactly. Units of the lapse rates are M^-1.

| level | MPG lapse peak | MPG lapse cell RMS | old lapse peak | old lapse cell RMS |
| --- | --- | --- | --- | --- |
| 0 | 3.88085539e-19 | 4.24208999e-22 | 9.27075232e-05 | 1.70272691e-07 |
| 1 | 1.64795954e-18 | 4.81501096e-21 | 1.56147230e-04 | 1.13799754e-06 |
| 2 | 3.81309001e-18 | 1.43676145e-20 | 2.85891011e-04 | 2.45047466e-06 |
| 3 | 7.17177402e-18 | 2.72101672e-20 | 6.18826205e-04 | 5.29170673e-06 |
| 4 | 5.91032446e-17 | 1.55209723e-19 | 1.46378101e-03 | 1.16310210e-05 |
| 5 | 3.99172532e-17 | 1.87022346e-19 | 2.97470911e-03 | 2.52851055e-05 |
| 6 | 7.16650734e-17 | 3.58398276e-19 | 3.93112912e-03 | 5.20904532e-05 |
| 7 | 6.47673008e-17 | 7.00944510e-19 | 4.68774381e-03 | 1.01625815e-04 |
| 8 | 8.00537126e-17 | 1.43763083e-18 | 4.71271033e-03 | 1.87819242e-04 |
| 9 | 1.11634311e-16 | 2.90068098e-18 | 4.74261795e-03 | 3.41702984e-04 |
| 10 | 1.08604776e-16 | 5.73088794e-18 | 4.74341182e-03 | 6.31221714e-04 |
| 11 | 1.22077290e-16 | 1.11006055e-17 | 4.74344112e-03 | 1.17202907e-03 |
| 12 | 1.22659893e-16 | 1.99010897e-17 | 4.74344511e-03 | 2.11303191e-03 |
| 13 | 1.50236578e-16 | 3.02760026e-17 | 4.74344748e-03 | 3.42849040e-03 |
| 14 | 1.39199660e-16 | 3.36817070e-17 | 4.74344750e-03 | 4.34323102e-03 |

[t22-initial-lapse-rates.csv](t22-initial-lapse-rates.csv) includes coordinate-volume RMS and cell counts for every level. The all-valid-cell rates use the checkpoint's stored stencil halos and actual native derivative/gauge classes; no new initialization or evolution is performed. At level 14 the new lapse cell RMS is 3.36817070e-17 versus old 4.34323102e-3. Ordinary KO(lapse), which is unchanged and can dominate this roundoff-level pre-KO gauge row, is reported separately below.

## Complete discrete Gamma and driver pieces

The new B pre-KO RHS equals the **full discrete Gamma RHS bit for bit** at every native sampled point, including the EMS matter momentum term. B=0 and advection is disabled. The geometric column is the native RHS with Newton coupling zero on a private current-state evaluation; the matter column is independently reconstructed from the native EMS stress tensor. Their sum agrees with full Gamma within the original arithmetic budget (measured budget fraction zero). Ordinary KO(B)=0 initially because B=0; KO(Gamma) remains an independent Gamma row and is never inserted into B.

| level | points | geometric Gamma peak / RMS | matter Gamma peak / RMS | full Gamma = MPG B peak / RMS | KO(Gamma) peak | KO(B) peak |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 32 | 7.43509199e-04 / 1.78963453e-04 | 9.64006948e-04 / 2.52755341e-04 | 1.70751615e-03 / 3.36834188e-04 | 5.51296891e-16 | 0.00000000e+00 |
| 1 | 32 | 3.02159552e-03 / 1.09274798e-03 | 6.74747880e-03 / 1.83231058e-03 | 9.76907432e-03 / 1.95191297e-03 | 2.22621666e-15 | 0.00000000e+00 |
| 2 | 32 | 3.84005932e-02 / 8.93396665e-03 | 3.77210089e-02 / 1.05421220e-02 | 7.61216021e-02 / 1.56676924e-02 | 8.35428665e-15 | 0.00000000e+00 |
| 3 | 32 | 3.52694042e-01 / 7.44966154e-02 | 1.75557294e-01 / 5.06112653e-02 | 5.28251337e-01 / 1.09858249e-01 | 2.61984312e-14 | 0.00000000e+00 |
| 4 | 32 | 1.97560737e+00 / 4.06959626e-01 | 6.18257999e-01 / 1.93101859e-01 | 2.59386537e+00 / 5.28797131e-01 | 5.65262123e-14 | 0.00000000e+00 |
| 5 | 32 | 6.05804686e+00 / 1.23090278e+00 | 1.34905950e+00 / 5.43531806e-01 | 7.40710636e+00 / 1.43381062e+00 | 3.87118908e-13 | 0.00000000e+00 |
| 6 | 32 | 9.88417291e+00 / 2.21601681e+00 | 2.83773071e+00 / 1.17318234e+00 | 1.14484862e+01 / 2.16065920e+00 | 2.41765791e-12 | 0.00000000e+00 |
| 7 | 35 | 1.03519708e+01 / 2.89093346e+00 | 3.25420447e+00 / 1.86293788e+00 | 1.13641044e+01 / 2.26431011e+00 | 1.03675436e-11 | 0.00000000e+00 |
| 8 | 40 | 8.35575657e+00 / 2.83290738e+00 | 3.34776387e+00 / 2.05796444e+00 | 8.80179033e+00 / 1.88899765e+00 | 4.97078297e-11 | 0.00000000e+00 |
| 9 | 44 | 5.93473662e+00 / 2.17748361e+00 | 3.38250003e+00 / 1.62167553e+00 | 6.09630567e+00 / 1.39971757e+00 | 1.39445027e-10 | 0.00000000e+00 |
| 10 | 50 | 4.56014724e+00 / 1.66438200e+00 | 3.38731188e+00 / 1.31536737e+00 | 4.50645314e+00 / 9.92146388e-01 | 7.14640102e-10 | 0.00000000e+00 |
| 11 | 52 | 3.58817860e+00 / 1.41063667e+00 | 3.38655493e+00 / 1.18905179e+00 | 3.57095554e+00 / 7.47679686e-01 | 2.33135046e-09 | 0.00000000e+00 |
| 12 | 55 | 3.38542059e+00 / 1.27754471e+00 | 3.38541216e+00 / 1.14246792e+00 | 2.82535487e+00 / 5.67774572e-01 | 8.17877031e-09 | 0.00000000e+00 |
| 13 | 53 | 3.38627002e+00 / 1.23195824e+00 | 3.38626949e+00 / 1.14449989e+00 | 2.23572834e+00 / 4.54664837e-01 | 3.51206019e-08 | 0.00000000e+00 |
| 14 | 48 | 3.38591205e+00 / 1.09457802e+00 | 3.38591201e+00 / 1.02736444e+00 | 1.76946008e+00 / 3.77251603e-01 | 1.26498314e-07 | 0.00000000e+00 |

These Gamma/B RMS values are explicitly **sample RMS**, not full-mask norms: 601 native points per gauge across all 15 levels, comprising 32 puncture-neighbourhood cells plus retained in-box axis/equator/diagonal anchors through R=0.02 M. A finite-difference Gamma residual is expected; it is neither asserted zero nor removed. On level 14 its peak/RMS are 1.76946008 / 0.377251603, while KO(Gamma) peaks at 1.26498314e-7. Geometric and matter peaks separately reach about 3.385912 and partly cancel. Full signed pointwise values and coordinates are in [t22-equilibrium-native.csv](t22-equilibrium-native.csv); all per-level sample counts, component-combined RMS, geometric/matter/KO and old-gauge rates are in [t22-equilibrium-summary.csv](t22-equilibrium-summary.csv).

| level | old shift peak | MPG shift peak | old offset B rate peak | KO(lapse) peak | KO(shift) peak |
| --- | --- | --- | --- | --- | --- |
| 0 | 1.49505203e-08 | 0 | 1.53323601e-05 | 3.00275977e-02 | 1.30471544e-04 |
| 1 | 3.59860414e-08 | 0 | 1.97741776e-05 | 3.55331932e-02 | 3.54885929e-04 |
| 2 | 2.59738719e-07 | 0 | 2.88392856e-05 | 4.38763892e-02 | 1.12378185e-03 |
| 3 | 1.56403444e-06 | 0 | 4.46782610e-05 | 6.50793725e-02 | 3.60263446e-03 |
| 4 | 4.81950086e-06 | 0 | 6.46273464e-05 | 1.14532323e-01 | 9.92681746e-03 |
| 5 | 1.36355473e-05 | 0 | 8.26054090e-05 | 2.04362552e-01 | 1.92978954e-02 |
| 6 | 4.73092386e-05 | 0 | 1.00637881e-04 | 2.93480249e-01 | 2.20871890e-02 |
| 7 | 8.30433237e-05 | 0 | 1.04418709e-04 | 2.98661543e-01 | 1.41328247e-02 |
| 8 | 1.19727345e-04 | 0 | 1.04717013e-04 | 2.36526929e-01 | 5.96970466e-03 |
| 9 | 1.25509596e-04 | 0 | 1.05796874e-04 | 1.75930023e-01 | 2.20049042e-03 |
| 10 | 1.20014182e-04 | 0 | 1.06140739e-04 | 1.33766732e-01 | 8.23233313e-04 |
| 11 | 1.18372776e-04 | 0 | 1.06267803e-04 | 1.03987986e-01 | 3.17273481e-04 |
| 12 | 1.17420258e-04 | 0 | 1.06223813e-04 | 8.17080476e-02 | 1.24197557e-04 |
| 13 | 1.16913397e-04 | 0 | 1.06250639e-04 | 6.44872023e-02 | 4.89392273e-05 |
| 14 | 1.17180781e-04 | 0 | 1.02995017e-04 | 5.09863586e-02 | 1.93354758e-05 |

The old offset driver and the new differential driver are different variables. The table reports their labelled source rows separately; it forms no cross-gauge B ratio. The initial Gamma, lapse KO and shift KO pieces coincide because their physical initial fields coincide.

## Stored ghosts and the native storage test

| level | stored ghost copies | K peak | K cell RMS | K / (eps scale) peak |
| --- | --- | --- | --- | --- |
| 0 | 148680 | 7.44007455e-19 | 2.29785292e-21 | 3.439124 |
| 1 | 28800 | 4.58991309e-18 | 5.06088097e-20 | 3.219725 |
| 2 | 28800 | 1.48961931e-17 | 2.16361461e-19 | 3.222619 |
| 3 | 28800 | 3.87186562e-17 | 5.17314031e-19 | 3.395017 |
| 4 | 28800 | 4.58416321e-16 | 4.78465596e-18 | 3.207745 |
| 5 | 28800 | 4.95453272e-16 | 6.53938498e-18 | 2.808733 |
| 6 | 28800 | 8.16610775e-16 | 1.29217458e-17 | 2.991330 |
| 7 | 28800 | 3.21814866e-15 | 3.95059623e-17 | 3.102722 |
| 8 | 28800 | 1.70165181e-15 | 4.28194430e-17 | 2.993220 |
| 9 | 28800 | 2.29871764e-15 | 7.05871481e-17 | 3.165336 |
| 10 | 28800 | 3.38582732e-15 | 9.92209391e-17 | 2.982914 |
| 11 | 28800 | 3.47584513e-15 | 1.49982889e-16 | 3.038738 |
| 12 | 28800 | 3.72395427e-15 | 2.54602328e-16 | 3.492669 |
| 13 | 28800 | 3.87052211e-15 | 4.68632084e-16 | 3.017319 |
| 14 | 28800 | 4.45645463e-15 | 7.63467095e-16 | 3.183948 |

Ghost copies are duplicated across boxes and are not physical-volume counts. Theta and new B are exactly zero throughout these copies; the largest ghost K is 4.45645463e-15 M^-1, scaled maximum 3.492669 eps. Their measured values pass the same numerical contraction bound, although the required valid-cell gate is kept distinct.

After replacing the K assertion, the old audit's internal total-RHS check exposed 92 one-component flags per gauge. [VarsTools.hpp](../../Source/utils/VarsTools.hpp#L41) maps both upper and lower symmetric tensor entries into the same evolved component; storage retains the last entry. The old test compared **both aliases** with that stored value. All flags are the overwritten upper A12 alias (component 7), not a mismatch of an evolved field. [t22-tensor-aliases.csv](t22-tensor-aliases.csv) retains every upper/lower discrepancy. The repaired Tests-only check compares exactly the 28 serialized, last-write values with the native stored output, still **bitwise with no tolerance**: 16,828 components per gauge, zero bit mismatches. [t22-native-storage-check.csv](t22-native-storage-check.csv) gives the counts. This does not change any production arithmetic or weaken the pinned default-path regression.

## Receipts and reproducibility

Both completed native initial runs are verified from their own actual child status, no gate event, finite native outputs, time-zero checkpoint and hierarchy census, not merely the launcher exit. Experimental: 45.709 s, peak RSS 1,642,774,528 bytes; moving_puncture: 43.159 s, peak RSS 1,369,243,648 bytes. Their logs say `T22_NATIVE_INITIAL_AUDIT_COMPLETE; no advances`. The sandboxed `time` sysctl warning has no scientific or child-exit meaning. [t22-equilibrium-resources.csv](t22-equilibrium-resources.csv) preserves these receipts. The offline trace/rate process stays below 0.6 GB, and the continuation qualification has its own receipt and done marker. All native runs use two OpenMP threads, below the four-thread limit. No completed evolution is repeated, no SSH/cluster operation or commit is made, and static data enter no evolution after t=0.

Regenerate with `t22-trace-audit.py compile-rates`, `t22-trace-audit.py analyze`, `t22-verify.py equilibrium`, then this script. Wrap the first three in `t13-run.py --measure` with the 6 GB per-process cap; run commands use absolute script paths because the measurement wrapper changes directory. Only stored numerical checkpoints are read. The design is updated only in its qualification note; `/private/tmp/ems-t22/design-frozen.md` and its registered SHA remain immutable.
