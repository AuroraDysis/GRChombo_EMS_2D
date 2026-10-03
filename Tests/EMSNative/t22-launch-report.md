# T22 short E-mid launch: qualified endpoint reductions, mixed early history

All four previously unrun trajectories completed their own native clean stops and passed receipt-first analysis. This is a local early-window audit, **not** the registered 10.5 M decision or a convergence statement. No completed initialization/evolution was repeated. The old and new packages keep identical delivered physical initial data, hierarchy, KO, transfers and t2 static guard; only gauge package and driver initialization differ. B variables are not compared across gauges.

The hierarchy has levels 0–12, h0=7/12 M, finest h=7/(12*4096)=0.00014241536458333333 M and dt/h=0.25 (0.125 for each same-gauge control). The nominal native stop is t=56*dt=0.001993815104166667 M, within the registered 0.002 M early window. The nominal trajectories have 57 saved level-12 clocks including t=0; dt/2 has 113, compared at every second clock. No unsynchronized plot/checkpoint is written. The native recorder is the qualified serial T13 stopping audit, not production recording.

Signed Gamma_n = (Gamma1,Gamma2).n is sampled on the axis and diagonal over the frozen radius window 0.00075<=R<=0.0025 M. The disturbance is the current numerical value minus its own stored numerical t=0 value. No static solution is read, subtracted or used as a reference. Peak is max absolute disturbance; RMS is its radial-integral RMS over the same window. I8 is central, with |I8-I10| at both current and initial clocks and the unchanged 128-eps Float64 floor. Each gauge must exceed five times its sampling/floor bound plus its own nominal-minus-dt/2 difference. Ratio bounds propagate these separately; the joint reduction margin is their sum.

## Endpoint at the common native stop

| ray | norm | old disturbance | new disturbance | new / old [bound] | 5x joint sampling | joint dt/2 difference |
| --- | --- | --- | --- | --- | --- | --- |
| axis | peak | 9.40815192e-05 | 2.58098307e-05 | 0.274335 [0.268169, 0.280555] | 8.25132803e-07 | 5.80076435e-08 |
| axis | RMS | 7.36255675e-05 | 1.01739495e-05 | 0.138185 [0.136440, 0.139937] | 2.30098168e-07 | 3.09934037e-08 |
| diagonal | peak | 9.17723960e-05 | 3.45297429e-05 | 0.376254 [0.373071, 0.379455] | 3.55679749e-07 | 9.36412291e-08 |
| diagonal | RMS | 7.47678423e-05 | 1.31366304e-05 | 0.175699 [0.174317, 0.177085] | 1.51761659e-07 | 4.58570490e-08 |

Every endpoint is qualified; even the upper ratio bounds are below one half. These are reductions of the measured Gamma disturbance in this window, not accuracy certificates. At the endpoint the joint temporal difference is smaller than the fivefold joint sampling margin in all four reads. [t22-launch-endpoint.csv](t22-launch-endpoint.csv) records the individual gauge temporal and sampling bounds, significance margins and reduction bounds.

## Entire common-time history

| ray | norm | qualified / positive clocks | unqualified positive steps | ratio range | max joint dt/2 / old |
| --- | --- | --- | --- | --- | --- |
| axis | peak | 53 / 56 | 2 3 4 | 0.274335–1.846309 | 0.023400 |
| axis | RMS | 56 / 56 | none | 0.138185–1.391008 | 0.018515 |
| diagonal | peak | 56 / 56 | none | 0.376254–2.508489 | 0.031083 |
| diagonal | RMS | 56 / 56 | none | 0.175699–1.146165 | 0.024395 |

The new gauge is **not uniformly quieter in the early history**. Both peak reads and RMS reads have resolved intervals where new exceeds old; all such steps are listed in [t22-launch-history-summary.csv](t22-launch-history-summary.csv). Axis peak reaches ratio 1.846309 and diagonal peak 2.508489 before their later decrease. In the first four steps, axis RMS ratios are 1.0723, 1.1520, 1.2594, 1.3872; diagonal RMS ratios are 1.0356, 1.1430, 1.1462, 1.1367. This short record has finite, bounded transients and no observed first-step runaway, but cannot diagnose a slow instability or decide the long-time activity screen.

| ray | norm | step | time M | failed signal | old significance margin | new significance margin |
| --- | --- | --- | --- | --- | --- | --- |
| axis | peak | 2 | 7.1207682292e-05 | old new | -9.73517639e-07 | -5.42512785e-07 |
| axis | peak | 3 | 1.0681152344e-04 | old new | -2.49153374e-06 | -1.30478242e-06 |
| axis | peak | 4 | 1.4241536458e-04 | old new | -2.82542713e-06 | -1.43782935e-06 |

The three remaining positive-time unqualified reads are axis peak steps 2–4; their failed signal is reported explicitly. All four t=0 reads are initial zero-disturbance floor entries and are also retained in [t22-launch-unqualified.csv](t22-launch-unqualified.csv). None is relabelled as a pass. [t22-launch-ratios.csv](t22-launch-ratios.csv) has every one of the 228 nominal ray/norm/clock reads, with old/new sampling and temporal margins separately. [t22-launch-history.csv](t22-launch-history.csv) preserves both full nominal and half-step histories. Joint temporal/old disturbance is at most 3.11% over positive clocks; per-gauge maxima are in the summary, without suppressing early sampling-dominated entries.

## Floors, finite values and receipts

No chi/lapse floor activation and no nonfinite value was recorded at any inspected stage/level, including ghost copies. [t22-launch-floors.csv](t22-launch-floors.csv) lists every level and run. The smallest recorded chi is 5.49151383e-7 and lapse 1.90911861e-4, both above their 1e-12 floors. Every sampled state and I10 stencil is finite. Each run's own done marker, actual child returncode zero, empty gate reason, native stop clock, state count, matching temporal clocks, floor rows and `GRChombo finished.` log are checked. The sandboxed time-wrapper sysctl warning is not a failure.

| case | own child exit | wall seconds | peak RSS bytes |
| --- | --- | --- | --- |
| experimental | 0 | 50.800 | 1590312960 |
| experimental-half | 0 | 72.234 | 1596751872 |
| moving_puncture | 0 | 43.703 | 1611841536 |
| moving_puncture-half | 0 | 65.208 | 1544093696 |

The maximum native launch RSS is 1,611,841,536 bytes (1.612 GB), below the 6 GB per-process cap. Two OpenMP threads are used. The original evidence pipeline exit 1 and false exact-zero assertion remain retained; only the unfinished cases ran under `/private/tmp/ems-t22/continuation/plan.json`. The continuation's own receipt and marker confirm completion, but its exit alone is not used as a scientific pass.

[t22-launch-history.pdf](t22-launch-history.pdf) plots peak and RMS disturbance norms against physical time, with both same-gauge dt/2 histories. The registered activity reading remains frozen in [t22-gauge-control-design.md](t22-gauge-control-design.md); H is an engineering comparator. New-gauge 10.5 M and temporal-control legs, target compiler/MPI identity preflight and whole-node calibration remain for the submission worker/controller. No local 10.5 M evolution, cluster operation or commit is made.
