# Frozen-checkpoint Newton experiment — 2026-10-05

READY-EXCEPT: the implementation is opt-in, but the step-72 root-agreement gate fails.
Steps 40/54 both hit TIME_CAP; the default-off mode is not accepted as a flow replacement.
No evolution, cluster job, push, EMS edit, static-data read, or post-pipeline edit occurred.
Pristine flow is built from `24eae0c`; portable implementation commit: `d9c4fd4`.
Raw runs, binaries, checks and stdlib/h5py analysis: run directory `/private/tmp/rh-newton-run/`, archived as `/Users/auroradysis/Workspace/EMS/.data/exp-0027/perf/rh-newton-run.tgz` (SHA-256 `28d02ca7f17d53b98944a079284e96f6988a2131c87df12fed990913a1c449e2`).
`offline_solver = flow | newton` defaults to flow; both RHUnion/RHSurf flags default off.
Newton solves the existing Theta_plus with the existing even/odd polar ghosts and expansion_error.
df/d2f use -2..2; DivS/KSSmK and all interpolated values/Cartesian derivatives use fields at i only. Thus b=2.
Reflection maps -1→0, -2→1, N→N-1, N+1→N-2, preserving the nonperiodic band.
Five colors perturb j modulo 5 together, ghost-fill and re-interpolate all 23 fields once/color.
For point j, h = sqrt(eps)*sqrt(max(f_j,dx_j)*max(max(f_j,dx_j),|centre_x|,|centre_y|)).
Float64 forward truncation h/r balances eps*centre/h rounding: seed h≈4–16e-8 exceeds the ≈4.5e-13 coordinate ulp.
The relative prefactor, radial bound in cells and number of backtracks are driver parameters.
invert_banded5 is used only for b=2, with a solve-residual gate and no patched diagonal; another band requires banded LU.
Every fresh-field trial stays in the radial limits, obeys |delta_f_j|≤1 old/trial local cell, and passes Armijo on expansion_error.
One Newton iteration/update; re-centring is followed by fresh interpolation in Newton mode.
On failure the entire offline search rewinds to the original resampled seed and replays flow.
The original thresholds 1e-7/1e-10/1e-12, N48, quota400 and 1790 s row timer are retained.
Probe time consumes that timer: strict no-worse-than-flow behavior at a time cap is NOT guaranteed.
Colored-versus-48-column FD Jacobians, each with its own re-interpolation: all five seeds
have max absolute/relative difference 0 and outside-band maximum 0, including both poles.
Cold-kernel check_newton.wls PROVED the exact N48 reflected support/color witness over Z/Q; numeric J checks are separate.
Default-off first-update geometry matches pristine flow; a forced failed-J replay has identical
shape/A/Q/residual/cells to default flow. Valid restart passes; negative-radius restart aborts.
All five checkpoint SHA-256 values matched their readout receipts before execution and after all runs (input-checks.json).
Private params.txt/rh_surf_0.dat/rh_f0.dat copies use the FIRST listed attempt (warm for all five); only restart_file and opt-in controls change.
`cases.json`, receipts and `measurements.json` retain full precision and input identities.
F/F/F means all three stages FOUND; T/-/- means TIME_CAP at stage 0. Updates are stage-wise;
interpolation counts include startup and the 54 calls for the optional J check in Newton runs.
Wall includes setup; solve is cumulative stage time. Overlapping local runs: OMP_NUM_THREADS=1, serial Chombo, GCC16 -O3, CH_TIMER=1.

| Step | Solver | Stages | Stage updates (Newton iterations) | Interpolations | Wall/solve s | Final expansion squared | A | Q | Min cells |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| 4 | flow | F/F/F | 615/794/193 | 641922 | 201.550/199.517 | 9.999953e-13 | 0.224507802628 | 1.07190061062 | 14.247706 |
| 4 | newton | F/F/F | 1/1/1 (2) | 73 | 24.051/0.008 | 1.313452e-15 | 0.224507781806 | 1.07190056458 | 14.247704 |
| 34 | flow | F/F/F | 2629/1245/301 | 1673649 | 530.635/517.770 | 9.999499e-13 | 0.256295326518 | 0.815194591061 | 9.487590 |
| 34 | newton | F/F/F | 2/1/1 (3) | 81 | 13.098/0.011 | 2.398118e-16 | 0.25629534526 | 0.815194598407 | 9.487591 |
| 40 | flow | T/-/- | 13105 | 5255106 | 1818.200/1790.080 | 1.273786e+00 | 0.326623513796 | 1.02417439197 | 6.423127 |
| 40 | newton | T/-/- | 12824 (5) | 5142953 | 1812.487/1790.047 | 7.787868e-01 | 0.372621451771 | 0.921410790031 | 8.478292 |
| 54 | flow | T/-/- | 12615 | 5058616 | 1822.460/1790.050 | 6.544120e+03 | 0.00910794047013 | 0.000802148398121 | 0.238388 |
| 54 | newton | T/-/- | 13033 (8) | 5226780 | 1792.225/1790.014 | 2.482533e+00 | 1.3976467155 | 0.943162107133 | 13.461364 |
| 72 | flow | F/F/F | 2311/753/182 | 1300747 | 457.775/452.906 | 9.999477e-13 | 50.3012690595 | 1.01656591867 | 81.013092 |
| 72 | newton | F/F/F | 45/1/1 (46) | 426 | 23.436/0.136 | 2.554275e-15 | 15.3422214746 | 1.00512793107 | 69.520579 |

Both methods' best logged update residuals: step40 0.1334502427764; step54 0.1357984081572.
Squared residual 1e-12 permits RMS Theta 1e-6; the 4/34 A/Q differences are consistent with stopping and far below the
supplied N48/N96 angular uncertainties ΔA=3.0e-5 and ΔQ=1.43e-4. Step72 finds another root.
Delta_f uses the finer spacing at each point; centres differ <1.8e-7 cells. Capped rows have no two-root difference.

| Step | max abs delta_f, cells | Newton minus flow A | Newton minus flow Q | Agreement |
|---|---:|---:|---:|---|
| 4 | 2.41601e-6 | -2.08218e-8 | -4.60404e-8 | PASS |
| 34 | 1.84065e-6 | +1.87418e-8 | +7.34565e-9 | PASS |
| 72 | 581.580 | -34.959047585 | -0.0114379876015 | FAIL: different branch |

A/(16*pi) for M=1 and min(f/local_dx) are necessary checks, not horizon certificates; QUALIFIED alone is insufficient.

| Step | Flow ratio; cells≥10 | Newton ratio; cells≥10 |
|---|---|---|
| 4 | 0.00446644082; yes | 0.00446644040; yes |
| 34 | 0.00509883351; NO (9.48759) | 0.00509883389; NO (9.48759) |
| 72 | 1.00071195183; yes, exceeds M=1 bound | 0.30522379821; yes, inner root |

The 57 readouts contain 56 active RH rows (row0 is sealed) and 74 recorded finder attempts.

| Cost model | 17 chunks, readouts / whole post | 57 row processes on 112 cores, readouts / whole post |
|---|---:|---:|
| Mixed nearest-sample, 74 attempts | 251.017 / 311.102 min | 64.002 / 124.087 min |
| Optimistic first-attempt success, 56 attempts | 13.657 / 73.742 min | 3.987 / 64.072 min |

Mixed assigns each row the full Newton wall of its nearest sampled step (ties earlier), including
failed warm costs for untested reseeds; all recorded retries/chunk assignments remain. Critical chunk: 011.
Optimistic assumes every active row finds the accepted branch first time, at the slowest successful cost
24.050586 s. This is unestablished: 40/54 fail and 72 selects another branch. Neither model is a forecast.
Both keep receipt overhead=row wall−all finder walls (6172.917039 s total), per-chunk overhead and
1.245304 s parent closure; whole post adds the unchanged serial N96 tail, 3605.099940 s.
Assume sampled Mac costs transfer unchanged to Leonardo, one thread/rank per process, serial retries
within each row, and no extra CPU/I/O/memory contention at 57 processes; that scaling is unmeasured.
Reproduction (fresh output directories; retained binaries reproduce the measured source states):

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
make -C Tests/EMSRHFinder all DIM=2 -j 8
R=/private/tmp/rh-newton-run
cp Tests/EMSRHFinder/EMSRHCheckpoint2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex "$R/bin/final.ex"
python3 Tests/EMSRHFinder/check_newton.py flow 4 34 --executable "$R/bin/main-flow.ex" --label main-flow
python3 Tests/EMSRHFinder/check_newton.py flow 40 54 --executable "$R/bin/main-flow.ex" --label measured-flow
python3 Tests/EMSRHFinder/check_newton.py flow 72 --executable "$R/bin/main-flow.ex" --label independent-flow
python3 Tests/EMSRHFinder/check_newton.py newton 4 34 40 54 --executable "$R/bin/final.ex" --label final-newton
python3 Tests/EMSRHFinder/check_newton.py newton 72 --executable "$R/bin/final.ex" --label measured-newton
uv run --with h5py --with numpy python "$R/analyze.py"
/opt/homebrew/bin/gtimeout -k 5s 60 /Applications/Wolfram.app/Contents/MacOS/wolfram -script Tests/EMSRHFinder/check_newton.wls
```
Pristine baseline used the same make command in `$R/commit-clone` at 24eae0c; retained as main-flow.ex. Use --output with copied cases.json
and readout-overheads.json, then analyze.py --root, for a repeat without overwriting this evidence.
The core step can use current evolved fields, but offline seed rewind assumes a frozen hierarchy;
live integration would require explicit flag wiring and validation. No live integration is claimed.
Worktree index-lock writes were refused; final portable commits/bundle/patch use `$R/commit-clone`.
