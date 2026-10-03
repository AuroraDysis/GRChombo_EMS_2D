# T19 — registered Stage D verdict of exp-0023

**Registered verdict: FAIL.** Resolved negative orders on fixed exterior masks at the registered positive clocks meet the contract's FAIL condition. This verdict does not classify positive orders above four, or narrow intervals excluding four, as failures of convergence. Those rows remain as measured. Four physically completed chains, small native constraint magnitudes and small horizon drifts do not establish the required field-and-constraint fourth-order convergence.

The production root is `/Users/auroradysis/SciData/EMS/.data/exp-0023/production` and the contract directory is `/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0023/submissions/exp-0023`. The script reads their tables and own numerical receipts only. No evolution, horizon solve, new interpolation, changed mask, clock, significance rule, tolerance or baseline is introduced. The supplied collection was SHA-256 verified by the controller (57,958 files); this analysis hashes every input it actually reads in `t19-manifest.txt`. It does not repeat or replace the controller's bulk collection verification.

## Registered reading and completeness

> PASS = on the smooth exterior masks at the synchronized diagnostic clocks t_k = 0.875*k M, k = 1…12, field and constraint orders consistent with 4 within the measured uncertainties, the sealed magnitudes met, the finest horizon drifts ≤ 1e-3 (area) / 2e-4 (charge) and decreasing under refinement; FAIL = negative or non-convergent orders on any mask outside measured uncertainty; anything else reported as measured, with the failing masks/times named.

The applied registration is `submit-contract.md` §8, including its controller-amended clocks and sealed magnitudes. The three spatial rung spacings are 7/8, 7/12 and 7/18 M, with ratio 3/2; max_level is 14. Positive acceptance clocks are exactly 0.875 k M for k=1…12. Numerical t=0 is reported separately and cannot trigger a positive-history verdict. All 13 fixed masks remain reported, including the outer boundary shell. For the exterior statement, the 12 radial masks have R≥0.005 M; the boundary shell is separately named, not silently used to establish an interior failure.

All four `run.done`/`chain.done` endpoint records agree at 10.5 M and steps 48/72/108/216. Evolution native exits, finite-state audits at all 13 common clocks, independent diagnostics and all 488 own N48/N96 qualification receipts were checked. Each accepted finder case has all three FOUND stages with squared expansion at most 1e-7/1e-10/1e-12 and negative inward expansion. The 13 E-high native-term replay receipts have unchanged current evolved bits, exactly matching saved constraint totals and term sums within their recorded roundoff allowance. Native pieces and signed-shell tables are produced only for the finest rung, as required by the sealed magnitude gate; no other-rung term tables are invented. The independent reduction receipt covers all 13 clocks, 42 quantities, 13 masks and I8/I10. Axis-parity receipts were checked independently. These are completion and diagnostic checks, not a convergence PASS.

The measured order is log(D_low,mid/D_mid,high)/log(1.5). Intervals are the collected I8–I10 plus Float64 bounds and are reproduced without widening or shrinking them. Both spatial differences must exceed five times their own spread-plus-floor. The half-step difference must be ≤0.2 of **each** spatial difference. Temporal uncertainty is reported separately in the source tables; a row exceeding that control cannot support a clean spatial-order statement, even if its displayed order looks fourth order. Undefined and sampling-limited entries remain unresolved. A negative central order whose interval includes zero is not established as negative outside uncertainty; its stored `non_convergent` label is preserved in the ledger rather than promoted to a strict FAIL witness.

## 1. Verdict and named exceptions

There are 344 qualified negative-order witnesses on radial exterior masks at positive clocks across the volume/geometry and ray tables, plus 7 negative direct constraint-to-zero pair orders outside their intervals and passing their significance/temporal checks. `t19-failing-orders.csv` names every such witness with mask, clock, field, geometry, differences, margins and order interval. `t19-order-exceptions.csv` also names every other non-qualifying self-difference entry; `t19-constraint-orders.csv` retains every constraint-to-zero exception. No exception is dropped by taking a median or restricting the report to category `all`.

For an explicit smooth-interior witness, Γ̃₁ on the horizon mask at t=1.75 M has p=-1.090154, interval [-1.090453, -1.089855], D_low,mid=8.070785e-10, D_mid,high=1.25569e-09, significance=13361.77, and temporal fractions 0.07027766/0.04517006. This resolves non-contraction at that clock under the registered uncertainty and temporal rules. It does not identify the physical or discretization mechanism from the tables alone.

| Fixed mask | Median p range over 28 fields, category all, k=1…12 | Qualified negative self-difference rows | Clock indices | Fields / constraints | Geometry / ray classes |
| --- | --- | --- | --- | --- | --- |
| puncture_excluded | 3.9885–4.7743 | 32 | 2;3;4;5;11;12 | A11;A12;A22;Aww;Ex;Ey;Gamma1;Gamma2;Pi;hww;phi | all;axis_junction;patch_edge;ray_axis;ray_diagonal;ray_equator;same_level_seam;smooth_interior |
| horizon | 2.8476–6.933 | 134 | 1;2;5;6;7;8;9;10;11;12 | A11;A12;A22;Aww;Bz;Ex;Ey;Gamma1;Gamma2;K;Pi;chi;h11;h12;h22;hww;phi | all;axis_junction;convex_corner;patch_edge;ray_axis;ray_diagonal;ray_equator;same_level_seam;smooth_interior |
| inside_inner_ring | 2.7033–6.6389 | 73 | 1;2;5;6;9;10;11;12 | A11;A12;A22;Aww;Bz;Ex;Ey;Gamma1;Gamma2;h11;h12;h22;hww;lapse;shift1;shift2 | all;axis_junction;convex_corner;patch_edge;ray_axis;ray_diagonal;ray_equator;same_level_seam;smooth_interior |
| cavity | 3.8181–5.2583 | 1 | 5 | A12 | convex_corner |
| between_rings | 3.9031–4.8867 | 9 | 4;5;10;11;12 | A12;Gamma1;Gamma2;chi;phi;shift2 | axis_junction;convex_corner;same_level_seam |
| outer_ring | 3.3575–5.863 | 22 | 4;5;6;7;8;10;11;12 | A11;A22;Aww;Bz;Ex;Gamma1;Gamma2;chi;phi;shift2 | axis_junction;convex_corner;ray_axis;ray_diagonal;same_level_seam |
| far | 3.9543–4.0225 | 9 | 5;6 | Gamma1;Gamma2;Ham;Mom1;Mom2 | patch_edge;ray_diagonal |
| cavity_core | 3.7974–5.7915 | 19 | 1;2;3;4;5;7;8;9;10;11;12 | Aww;Bz;Ex;Ey;Gamma1;Mom2;Pi;h22;hww;lapse;phi;shift1;shift2 | convex_corner;ray_axis;ray_diagonal;same_level_seam |
| outer_ring_core | 3.6471–5.0374 | 18 | 4;5;7;8;10;11;12 | A11;A22;Aww;Ex;Gamma1;Gamma2;Theta;Xi;chi;h11;phi;shift2 | axis_junction;ray_axis;ray_diagonal;same_level_seam |
| far_core | 3.9503–4.039 | 0 |  |  |  |
| exterior_wake | 3.96–5.1839 | 14 | 2;3;4 | Gamma1;Gamma2;Ham | all;axis_junction;patch_edge;ray_axis;ray_equator;same_level_seam;smooth_interior |
| receiving_side_3p5 | 3.9499–4.12 | 13 | 4;5;6 | Gamma1;Gamma2;Ham;Mom1 | all;patch_edge;ray_diagonal;same_level_seam;smooth_interior |
| outer_boundary_shell | 0.64566–2.1794 | 0 |  |  |  |

## 2. All evolved fields and geometry classes

The controller's full volume-table counts are independently reproduced: {'sampling_or_roundoff_limited': 9817, 'measured': 1774, 'consistent_with_4': 2101, 'measured_order_outside_4_interval': 24670, 'non_convergent': 404} over 38,766 rows, including t=0 and all 42 quantities. The independent ray table has {'sampling_or_roundoff_limited': 6619, 'consistent_with_4': 932, 'measured': 1160, 'measured_order_outside_4_interval': 10828, 'non_convergent': 117} over 19,656 rows. These counts are not restricted to the 28 evolved fields. The controller's reported median ranges reproduce when taking all 42 quantities together, including constraints (`t19-all-quantity-medians.csv`); they are not medians of the 28 evolved fields alone. The all-quantity early outer_ring/outer_ring_core medians are 3.349302/2.897805, and its boundary-shell range is 0.3518987–2.607843. `t19-field-summary.csv` supplies median/min/max defined order and every status count for exactly the 28 evolved fields at each mask × clock × existing geometry/ray class. It separately supplies the median after the spatial/temporal checks; neither median overrides a field's status. `t19-field-orders-outside.csv` lists every field with its actual order and interval excluding four, including temporally unqualified measurements, explicitly tagged. For a positive above-four example, K on puncture_excluded at t=0.875 M has p=4.054765, interval [4.00537, 4.104786]. It is as measured: neither a fourth-order pass nor a failure of convergence.

There are 91 absent common geometry-class records in `t19-absent-support.csv`; they remain absent rather than becoming zero-norm passes. The collected reducer explicitly skips the Cartesian outer boundary shell in its radial-ray loop (`production/reduce.py:170`); the 39 unavailable clock/ray combinations are recorded separately in that CSV. Patch-edge, convex-corner, axis-junction, same-level seam and smooth-interior classes follow the collected membership, and axis/equator/diagonal profiles retain their original ray weights. At k=1, outer_ring/outer_ring_core medians are 3.357503/3.647093. The broad median range near four elsewhere coexists with individually resolved negative and above-four orders.

## 3. Constraints and sealed magnitudes

`t19-constraint-orders.csv` contains every raw common-support low/mid/high constraint norm, the self-difference order and its interval, and both direct-to-zero orders and intervals. Direct-to-zero pair classifications at positive clocks, including the separately named boundary shell, are low–mid: {'consistent_with_4': 198, 'sampling_or_roundoff_limited': 595, 'measured_order_outside_4_interval': 1709, 'significance_or_temporal_unqualified': 80, 'negative_outside_uncertainty': 70}; mid–high: {'consistent_with_4': 260, 'sampling_or_roundoff_limited': 659, 'measured_order_outside_4_interval': 1380, 'significance_or_temporal_unqualified': 317, 'negative_outside_uncertainty': 36}. A direct-to-zero pair requires both norms to exceed five times their collected uncertainties; its registered temporal control is also retained. Raw native composite-grid constraint norms are reported in `t19-native-constraint-norms.csv`, and their native term-normalized RMS values are in `t19-magnitudes.csv`. The common-support constraint norms and native composite-grid norms use different supports and are not interchanged. The collected tables contain RMS orders, not peak-norm order intervals; maxima remain reported as raw values without inventing an uncertainty model or peak-order certificate.

Every applicable E-high normalized magnitude passes its sealed RMS threshold: 720 rows, comprising Hamiltonian, both momentum components, their joint norm and electric Gauss on all 12 radial masks and all 12 positive clocks. `t19-magnitudes.csv` records the actual winning individual term, its RMS, every other term's RMS, the raw residual and the recomputed ratio. `t19-native-term-scales.csv` contains the individual native pieces. The denominator is the largest individual term on that same native mask and clock, never a cancelled total or static reference.

| Constraint | Required rows | Largest normalized RMS | Sealed limit | Mask at largest ratio | t / M | Actual denominator | Denominator RMS | All met |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ham | 144 | 4.239974e-07 | 0.001 | receiving_side_3p5 | 4.375 | Ricci | 0.004321078 | True |
| Mom1 | 144 | 0.0001921884 | 0.001 | receiving_side_3p5 | 4.375 | cartoon_regularization | 5.18708e-06 | True |
| Mom2 | 144 | 0.0001129059 | 0.001 | far | 4.375 | cartoon_regularization | 1.960283e-06 | True |
| Mom | 144 | 0.0001434716 | 0.001 | receiving_side_3p5 | 4.375 | cartoon_regularization | 8.98435e-06 | True |
| GaussE | 144 | 5.20417e-09 | 0.001 | far_core | 0.875 | F_Eup_y_over_y | 0.0003771828 | True |

The largest absolute **signed** E-high shell ratio at positive clocks is 4.121209e-09 at t=10.5 M, N=96, below 2e-4. `t19-signed-GaussE-shell.csv` retains every produced E-high row and both angular resolutions, their own numerical Q_H(0), signed integrals, absolute integrals and two-cell boundary-band estimates. The signed physical-volume shell gate and the coordinate-volume normalized RMS gate both pass numerically. The measured angular/boundary-band scope is not a full shell quadrature certificate. Θ, Z, C_Γ, det(h)−1, tr(Ã), magnetic Gauss and Λ/Ξ have reported raw norms and orders, with no invented magnitude thresholds.

## 4. Temporal control

`t19-temporal-failures.csv` names all 4610 rows exceeding 0.2 in at least one pair fraction (1917 volume/geometry; 2693 ray), with both fractions, additive excess, factor above the limit, significance and original status. Of these, 599 have both spatial differences above five times uncertainty; these still cannot qualify a spatial-order claim. Zero-denominator fractions are separately counted as undefined in the JSON receipt and exceptions ledger.

The largest reported fraction is 155.9158 for det_h_minus_1, cavity_core, convex_corner, t=1.75 M, but its spatial significance is only 0.212224. Among rows exceeding the spatial significance rule, the largest fraction is 18.34119 for det_h_minus_1, inside_inner_ring, convex_corner, t=2.625 M. A large ratio on a tiny/floor-limited spatial difference is not evidence for a large absolute temporal error; both raw differences and floors remain in the ledger.

## 5. Every stored non_convergent row

`t19-non-convergent.csv` preserves all 521 stored labels with their exact differences, measured significance, I8/I10 spreads and Float64 floor. Per-pair floor flags explicitly test the five-times floor and five-times spread separately and together. No labelled non_convergent row has significance ≤5 under the table's own combined spread-plus-floor; near-roundoff absolute values alone do not justify relabelling a statistically resolved difference as roundoff-limited. Intervals straddling zero are identified and are not strict negative-order witnesses. Numerical t=0 and boundary-only cases remain reported outside the radial positive-time claim.

| Support | Stored non_convergent | Negative outside uncertainty, all clocks/masks | Interval includes zero | Significance ≤5 |
| --- | --- | --- | --- | --- |
| volume | 404 | 337 | 67 | 0 |
| ray | 117 | 109 | 8 | 0 |

## 6. Horizons and MQ

All 488 qualified numerical finder cases and stage residuals are retained in `t19-qualified-horizons.csv`; `t19-horizon-baselines.csv` selects the eight independent numerical-t=0 baselines. A positive-time numerical candidate is only a geometry seed when solving the t=0 fields. The native candidate `finder-values.csv` modes are not substituted for these qualified surfaces. `M_RN_legacy` is the emitted RN expression and is not promoted to an EMS equilibrium mass or a new horizon mass diagnostic.

| Leg | Nθ | Numerical A(0) | Numerical Q(0) | Final squared expansion | Own wall seconds |
| --- | --- | --- | --- | --- | --- |
| E-low | 48 | 0.3837463 | 1.048593 | 9.999772e-13 | 218.73932217573747 |
| E-low | 96 | 0.3836949 | 1.048453 | 9.976323e-13 | 1.3528791489079595 |
| E-mid | 48 | 0.3837463 | 1.048593 | 9.999132e-13 | 206.49662777315825 |
| E-mid | 96 | 0.3836949 | 1.048453 | 9.999939e-13 | 1.3044138234108686 |
| E-high | 48 | 0.3837463 | 1.048593 | 9.999788e-13 | 262.827638088027 |
| E-high | 96 | 0.3836949 | 1.048453 | 9.997517e-13 | 1.3507234170101583 |
| E-high-dt2 | 48 | 0.3837463 | 1.048593 | 9.999694e-13 | 278.23102558986284 |
| E-high-dt2 | 96 | 0.3836949 | 1.048453 | 9.997394e-13 | 1.2994970721192658 |

`t19-horizon-drifts.csv` reports the full spatial coarse histories (49/73/109 clocks per N) and only the registered 13 diagnostic clocks for the half-step control. No between-diagnostic-clock control maximum is claimed. The budgets and observed decrease under spatial refinement are met in the drift table.

| Leg | Quantity | Nθ | Largest signed relative drift | Budget | Required history complete | Within budget |
| --- | --- | --- | --- | --- | --- | --- |
| E-low | A | 48 | -6.344878e-07 | 0.001 | True | True |
| E-low | A | 96 | -6.344259e-07 | 0.001 | True | True |
| E-low | Q | 48 | -7.130233e-08 | 0.0002 | True | True |
| E-low | Q | 96 | -7.158267e-08 | 0.0002 | True | True |
| E-mid | A | 48 | -2.973145e-07 | 0.001 | True | True |
| E-mid | A | 96 | -2.97233e-07 | 0.001 | True | True |
| E-mid | Q | 48 | -7.515553e-09 | 0.0002 | True | True |
| E-mid | Q | 96 | -7.455081e-09 | 0.0002 | True | True |
| E-high | A | 48 | -2.510418e-07 | 0.001 | True | True |
| E-high | A | 96 | -2.510315e-07 | 0.001 | True | True |
| E-high | Q | 48 | -1.080348e-09 | 0.0002 | True | True |
| E-high | Q | 96 | -1.096404e-09 | 0.0002 | True | True |
| E-high-dt2 | A | 48 | -2.511057e-07 | 0.001 | True | True |
| E-high-dt2 | A | 96 | -2.511021e-07 | 0.001 | True | True |
| E-high-dt2 | Q | 48 | -1.073011e-09 | 0.0002 | True | True |
| E-high-dt2 | Q | 96 | -1.088157e-09 | 0.0002 | True | True |

All 96 positive-clock retention comparisons have sensitivities within the drift budgets and resolved **drift** differences under refinement, but 96 have absolute angular/stopping sensitivity larger than the minimum absolute inter-rung difference. `t19-horizon-sensitivities.csv` preserves the `interrung_sensitivity_unresolved` statuses. Stable correlated N48/N96 drifts do not resolve that absolute comparison, so the stronger horizon qualification clause is not established even though the budgets pass.

All 444 produced MQ rows, including off-ladder coarse clocks, are in `t19-mq.csv`. They are code-produced radius-20/50/100 mass/charge diagnostics, independent of the horizon baselines. No t=0 MQ row is fabricated and no new flux or loss threshold is imposed. The first and final reported values are:

| Leg | Radius / M | Produced rows | First t / M | Last t / M | First M | Final M | First Q | Final Q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| E-low | 20 | 48 | 0.21875 | 10.5 | 0.9949206 | 0.9949224 | 1.048406 | 1.048406 |
| E-low | 50 | 48 | 0.21875 | 10.5 | 0.9979967 | 0.9979967 | 1.048406 | 1.048406 |
| E-low | 100 | 48 | 0.21875 | 10.5 | 0.9990032 | 0.9990032 | 1.048406 | 1.048406 |
| E-mid | 20 | 72 | 0.1458333 | 10.5 | 0.9949206 | 0.9949225 | 1.048406 | 1.048406 |
| E-mid | 50 | 72 | 0.1458333 | 10.5 | 0.9979967 | 0.9979967 | 1.048406 | 1.048406 |
| E-mid | 100 | 72 | 0.1458333 | 10.5 | 0.9990032 | 0.9990032 | 1.048406 | 1.048406 |
| E-high | 20 | 108 | 0.09722222 | 10.5 | 0.9949206 | 0.9949225 | 1.048406 | 1.048406 |
| E-high | 50 | 108 | 0.09722222 | 10.5 | 0.9979967 | 0.9979967 | 1.048406 | 1.048406 |
| E-high | 100 | 108 | 0.09722222 | 10.5 | 0.9990032 | 0.9990032 | 1.048406 | 1.048406 |
| E-high-dt2 | 20 | 216 | 0.04861111 | 10.5 | 0.9949206 | 0.9949225 | 1.048406 | 1.048406 |
| E-high-dt2 | 50 | 216 | 0.04861111 | 10.5 | 0.9979967 | 0.9979967 | 1.048406 | 1.048406 |
| E-high-dt2 | 100 | 216 | 0.04861111 | 10.5 | 0.9990032 | 0.9990032 | 1.048406 | 1.048406 |

## 7. Figures and reproducibility

`t19-stageD-history.pdf` has 26 pages: ten fixed representative quantities on every one of the 13 sealed masks. Each row shows common-support coordinate-volume **RMS norms against time**, both raw spatial differences and the temporal difference, the finer spatial difference rescaled by (3/2)^p for p=2,4,6, and the observed order with its stored interval. These rescalings are displayed assumptions, not alternate verdicts. `t19-stageD-field-orders.pdf` has one page per mask summarizing all 28 fields' order range, medians and status counts. `t19-stageD-overview.pdf`/`.png` show the smooth-interior horizon Γ̃₁ witness; its negative order at 1.75 M is resolved while both spatial differences greatly exceed their diagnostic spread. Unqualified central orders are marked and t<0.875 M is shaded outside the claim. Raw maxima are available in the original and exported native norm tables but are not the plotted norm. `t19-figures.csv` indexes every page. No superseded run is plotted.

Regenerate with `python Tests/EMSNative/t19-stageD-verdict.py COLLECTED_PRODUCTION_ROOT CONTRACT_DIR [--output-dir DIR]`. The script checks duplicate keys, frozen mask/clock coverage, every stored self-difference order/interval/significance/status, direct constraint orders, native term denominators and numerical horizon stages before writing the verdict. Its own resource receipt is `t19-stageD-verdict.json`; SHA-256 input/output seals are in `t19-manifest.txt`. Only t19-* outputs and README §T19 are written in this worktree. EMS main remains read-only; no SSH, cluster work, new evolution or commit occurs.

## T20 ledger housekeeping

The five large ledgers default to `/Users/auroradysis/Workspace/EMS/.data/exp-0023/t19` in the generator; `--output-dir DIR` overrides that directory. The move is blocked by the session filesystem sandbox; all five worktree sources and their current references are retained. The per-file sizes, hashes, destinations and status are in [t20-housekeeping.csv](t20-housekeeping.csv). The original T19 numerical reading and its resource receipt are unchanged.
