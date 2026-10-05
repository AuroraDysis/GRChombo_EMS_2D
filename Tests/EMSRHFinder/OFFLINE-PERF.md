# Frozen checkpoint RH cost — 2026-10-05

Leg X step 4; baseline `13b5010bc84ba8d4e3d115a08d138eb8e97674f6`; relevant sources match deployed `758f0d2`.
Checkpoint: `/Users/auroradysis/Workspace/EMS/.data/exp-0027/perf/EMS_000004.2d.hdf5`.
Verified SHA-256: `9bf953c0109cfce3bdd5678e547964b466402a5035afaee49cb06ebb888c1f33`.
macOS arm64, GCC 16.2, -O3, serial Chombo, OMP_NUM_THREADS=1, CH_TIMER=1.
All three executables were clean-built in this worktree and run sequentially.

Correction to the initial cost model: resample() overrides RH_time_step_freq=1
with quota 400. Each update permits 402 interpolations of 23 fields at 48 points.
Quota/chase settings, the three thresholds, 1790 s cap and threads are unchanged.
RHUnion's skip flag defaults false; only EMSRHCheckpoint enables it after the
one-time refresh. All live callers retain refresh. The extra driver interpolation
was removed only after its separate identity gate passed.

| Measurement | Unchanged baseline | Flag only | Flag + driver deletion |
|---|---:|---:|---:|
| Flow seconds/update | 0.40850945 | 0.12198192 | 0.11331486 |
| Refresh/ghost fill seconds/update | 0.29349224 | 0 | 0 |
| Interpolation seconds/update | 0.10148347 | 0.10760789 | 0.09955407 |
| Other seconds/update | 0.01353374 | 0.01437404 | 0.01376079 |
| Full process wall, s | 667.972469 | 208.527695 | 194.623434 |
| Stage updates | 615/794/193 | 615/794/193 | 615/794/193 |
| Total updates | 1602 | 1602 | 1602 |
| Refresh timer s/calls | 470.46806/1603 | 0.93088/1 | 0.80501/1 |
| Interp timer s/calls | 162.57677/643524 | 172.38810/643524 | 159.48587/641922 |

Split estimates subtract the single setup call at mean timer cost; unchanged
source does not separately time that call. Other = flow wall/update − refresh −
interp. calculateAnswers is included in interp. Exit time.table used once.
Final process speedup: 3.432x; baseline refresh was 71.84% of flow.
Deletion saves 1602 interpolation calls (0.25%); timings are single observations,
so the flag-only/final timing difference is not attributed solely to deletion.

Leonardo: 855.468662 s process, 853.348381 s flow, 1602 updates (615/794/193).
All local runs give A=0.22450780262796724, Q=1.07190061061694.
Mac − Leonardo: ΔA=+1.4156176231239215e-12, ΔQ=+4.239941731043473e-12.
No local run approached 30 minutes; no capped substitute was needed.

| Identity file | Baseline vs flag | Flag vs deletion; baseline vs final |
|---|---|---|
| progress.csv | byte-identical | byte-identical |
| offline_rh_f0.dat | byte-identical | byte-identical |
| offline_rh_surf_0.dat | byte-identical | byte-identical |
| shape-0--1.dat | byte-identical | byte-identical |
| shape-0-0.dat | byte-identical | byte-identical |
| shape-0-1.dat | byte-identical | byte-identical |
| shape-0-2.dat | byte-identical | byte-identical |
| surfaces.csv excluding only seconds | byte-identical | byte-identical |

Inventory reverified against every progress.csv: 230029 updates / 74 attempts,
including 18 reseeds (79800 updates); all are N=48. Recursive readouts/*/horizon/**
receipt coverage includes all locally present attempts; no attempt was omitted.
Fixed-work finder = 0.11331486142 × 230029 + 13.09302583 × 74 = 27034.588 s (7.51 h work).
Recorded non-finder/row overhead adds 6172.917 s: 33207.505 s (9.22 h serialized readouts).
Preserving the recorded 17-worker chunk assignments, critical chunk 011 has
31423 updates / 8 attempts: projected readout wall = 4388.661 s (73.14 min).
Holding the outside-inventory N=96 final-angular tail at its recorded 3605.100 s
gives 7993.761 s ≈ 2 h 13 min, plus unmodelled launch/closure overhead.
Assume Mac costs transfer to the recorded concurrent Leonardo workers/overheads and counts stay fixed.
Capped rows may now perform more updates/retries, so this is not a completion forecast.

All raw evidence/stdlib helpers: run directory `/private/tmp/rh-offline-perf-run/`, archived at `/Users/auroradysis/Workspace/EMS/.data/exp-0027/perf/rh-offline-perf-run/`. Inputs copied from X/000004/horizon/warm;
the ONLY parameter edit was restart_file. Full ordered source/build/run reproduction: `commands.sh` there.
Per-variant commands (fresh directories; baseline, then flag-only skip, then tested deletion candidate):
```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
make -C Tests/EMSRHFinder clean DIM=2 && make -C Tests/EMSRHFinder all DIM=2 -j 8
R=/private/tmp/rh-offline-perf-run; code=baseline # repeat with skip/candidate source states
cp Tests/EMSRHFinder/EMSRHCheckpoint2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex "$R/bin/$code.ex"
nohup python3 "$R/run.py" "$code" > "$R/$code-launch.log" 2>&1 < /dev/null & wait $!
python3 "$R/measure.py" "$code" # run.py supplies OMP_NUM_THREADS=1 CH_TIMER=1
python3 "$R/compare.py" baseline skip
python3 "$R/compare.py" baseline candidate
python3 "$R/compare.py" skip candidate
python3 "$R/projection.py" candidate
```
Worktree git refused the ref lock; final commit uses /private/tmp/rh-offline-perf-r2.git.
Final bundle/patch: archived as `/Users/auroradysis/Workspace/EMS/.data/exp-0027/perf/rh-offline-perf-r2.{bundle,patch}`.
