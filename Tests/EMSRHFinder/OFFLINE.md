# Strict horizons on frozen checkpoints

`offline.py` loads a production checkpoint with its production parameter file,
restores the saved centre and full angular shape through the existing RHUnion
restart reader, and re-finds each active surface without advancing the hierarchy.
The frozen level throws if initial-data construction, evolution or regridding is
requested. RHUnion and RHSurf are unchanged. No mass code is added to the evolution.

The command must run where the checkpoint is already accessible. It does not
connect to another host or submit a job. From a built fork, one invocation handles
several times:

```sh
uv run --with h5py --with numpy python Tests/EMSRHFinder/offline.py \
  /data/run/params.txt /data/run /data/offline-new \
  /data/run/chk/EMS_000160.2d.hdf5 /data/run/chk/EMS_000320.2d.hdf5 \
  --ems-root /data/EMS --points 96 192 --threads 8 --cap 900
```

The output directory must be new. Add `--dry-run` to check times and seed rows and
print the resource plan without creating output. `--executable PATH` selects a
specific build. The default cap is 240 seconds per checkpoint/resolution; the
maximum permitted cap is 1800 seconds. A 5–30 minute case needs a short pilot and
an explicit cap; a measured projection above 30 minutes is blocked. A cap is a
resource limit, not evidence of a residual floor. Batch wall time scales with the
number of checkpoint/resolution pairs. The driver exits nonzero if any active
surface misses the final target, while retaining completed measurements.

The build uses the same dimension, variable layout, Chombo and HDF5 toolchain as
production, with OpenMP enabled:

```sh
export CHOMBO_HOME=/path/to/Chombo/lib
make -C Tests/EMSRHFinder all DIM=2 -j 8
```

Runtime resources are one process/rank and eight CPU threads, Python >=3.11 with
NumPy and h5py, Julia with the existing EMS `--project=test` environment already
instantiated, and the unchanged certified table
`EMS/artifacts/horizon-family/production/alpha20-positive.hfamily`. The pinned
payload is `dd4013ca19a746b5bab847fcaaf7e9ece4d5ddf7f4ba5b29d84cb86386b23a07`.
That same family contains the exact q=0.7 member; no separate table is fitted or
regenerated. The input checkpoint, original parameter file, and both
`rh_surf_ID.dat` and `rh_fID.dat` histories are required. Initial-data files and
plot files are not read by the offline executable. MPI scaling is untested; the
Python driver is launched once, not once per MPI rank. Budget memory as 2 GiB
plus six times the checkpoint field payload for a first allocation; this is a
conservative planning estimate, not a measured production-memory bound.

## Measurement and acceptance

The driver selects the last saved pair at the checkpoint time, matching the
native writer's nine significant digits. It refuses a missing or inconsistent
active seed instead of silently taking an older tracking row. Future rows are
not consumed. A previously dead surface remains dead; a not-yet-started surface
remains dormant. Only private copies reach the restart reader's pruning path.
Input hashes are checked again after processing. At t=0 a positive adjacent
Float64 restart time enables the same reader to restore the zero-time row.

The original-resolution seed is first replayed against the checkpoint. Each
requested angular count starts from a linear resampling of that saved shape
with even polar ghosts. The unchanged finder then runs at chase 0.125 and quota
400 through thresholds 1e-7, 1e-10 and 1e-12. Acceptance uses freshly interpolated
fields at the final centre and surface, including after the finder's re-centring.
A floor of the worst active-surface residual is reported only when its window minimum and maximum change
by less than 1% between consecutive 64-update windows. This is an empirical
stagnation test, not a proof that further iteration cannot improve the result.
It stops the current group: a stalled search can prevent the remaining searches
from proceeding to a tighter threshold. Individual fresh residuals and achieved
statuses remain explicit; this case would need separate per-surface stopping.
Time/update caps and floors have separate statuses. `progress.csv` retains the
fresh residual, A, Q and centre at every update so that this decision is auditable.

Area and charge are RHSurf's own `Area()` and `Q_charge()`. The area weights and
interpolated scalar already held by RHSurf give the area-weighted scalar mean
and RMS fluctuation about that mean. The latter is not the RMS field amplitude.
The output preserves all native tracking A/Q rows used as seeds, their fresh
checkpoint replay, and all threshold-stage measurements. Native tracking rows
are labelled with their original solver status. Their scalar columns refer to
the replay of their saved shape on the checkpoint.

The existing Julia post-processor reads the unchanged table, binds the coupling,
branch and scalar boundary condition, and calls `horizon_mass` for M_eq and its
first-law uncertainty. The requested angular error estimate is the larger of
the analytic midpoint-sphere bias and the measured adjacent-resolution estimate
assuming second order (one third of the coarse/fine difference for the fine
value, four thirds for the coarse value). The stopping estimate is the absolute
change from the preceding threshold, separately for A and Q. They are added
before the existing first-law propagation; the table allowance remains 1e-8
relative. For the final target this uses the measured 1e-10 to 1e-12 change.
The lower-threshold rows and seed replays are diagnostic, not final precision
certificates. Duplicate geometries use the existing mutual-enclosure test;
separated and nested surfaces remain distinct, and duplicate spreads stay explicit.

This estimates angular and stopping errors on a particular saved grid. It does
not include evolution/discretization error, interpolation error of an unknown
binary grid, branch ambiguity or departure from equilibrium. The static control
below checks the combined result against an exact member. The binary study is
angular convergence on one frozen hierarchy, not a continuum evolution or
constraint-convergence certificate. Every mass retains
`BRANCH_UNVERIFIED_STATIC_PROJECTION`; `M_settled` is unset.

Each checkpoint has `masses.csv`, the pre-Julia `measurements.csv`, `inputs.json`
with SHA-256 identities, and `status.json` with process duration/exit/peak child
RSS. Each angular directory contains the effective parameter file, copied seeds,
`surfaces.csv`, `progress.csv`, native-format `offline_rh_surf_*`/`offline_rh_f*`
histories, and `shape-ID-STAGE.dat` with 17-digit time, centre and radii. Stages
-1/0/1/2 mean seed replay / 1e-7 / 1e-10 / 1e-12. Skipped surfaces are listed in
`skipped.csv`; their geometry is never revived. The mass CSV additionally uses
stage -2 for original in-run A/Q. No plot file is needed.

## Local reproduction

Raw data are in `/private/tmp/ems-offline-f`. The fork is based on
`1e8684cd6b4738f6c5c3220de59979e46c13113a`; no commit, SSH operation, HPC run,
author-finder change or EMS.jl change is part of this work.

```sh
uv run --with h5py --with numpy python Tests/EMSRHFinder/check_offline.py \
  prepare /private/tmp/offline-repeat --ems-root /Users/auroradysis/Workspace/EMS
uv run --with h5py --with numpy python Tests/EMSRHFinder/check_offline.py \
  check /private/tmp/offline-repeat
# Then run offline.py on binary/chk/EMS_000003.2d.hdf5 and
# static32/chk/EMS_000000.2d.hdf5, using --points 48 96 192.
# For the static case also pass:
# --expected-mass 1 --expected-charge .7 --expected-area 38.530469490246766
# Finally run check_offline.py results CHECKPOINT_OUTPUT_DIRECTORY.
```

The pilot uses L=96, 192x96 coarse cells, max_level=3, max box/block factor 16,
separation 32, the q=0.7 trumpet and degree-28 CTT companion, exact rapidity
0.2027325540540822 and bh_charge=0.7. Its finest spacing is M/16, as in the reduced
performance pilot. It takes three coarse steps to t=0.375, writes checkpoints at
every coarse step, and tracks at chase 1, threshold 1e-7. The static control is
an isolated, unboosted q=0.7 member at t=0, uniform dx=M/32, initially found from
radius .64. All new checkpoint finds use eight OpenMP threads.

## Local results

Both binary surfaces passed all three thresholds at all three angular counts.
The table gives either individual hole at t=0.375; the two results agree to
1.3e-13 in area and 1.3e-15 in mass. The in-run mass below is the Julia projection
of the saved tracking A/Q, not the author's legacy RN mass column.

| Binary measurement | N_theta | A | Q | M_eq | estimated delta_M |
|---|---:|---:|---:|---:|---:|
| In-run tracking | 96 | 38.4323457 | 0.700036401 | 0.9989752546 | diagnostic only |
| Offline, 1e-12 | 48 | 38.4422634469 | 0.700126736877 | 0.9991041473 | 1.10231e-4 |
| Offline, 1e-12 | 96 | 38.4373016478 | 0.700036250820 | 0.9990274965 | 3.07539e-5 |
| Offline, 1e-12 | 192 | 38.4360574115 | 0.700013647385 | 0.9990082982 | 1.08886e-5 |

The area differences have observed angular order **1.99560**, against the
predeclared midpoint expectation of two. At N=192 the native scalar mean is
0.0550196736207, with RMS fluctuation 3.11160604e-5. The fresh mean-square outgoing
expansion is 9.99975e-13 or less for both holes, and the inward expansion is
negative. Restoring the original tracking geometry gives A=38.4323456542253 and
Q=0.700036400595933, consistent with its nine-digit saved row.

The isolated static reference has M=1, Q=0.7, A=38.530469490246766 and
phi_H=0.052274528945421062. Its initial tracking mass was 1.0001747498. The strict
checkpoint results are:

| Static t=0 | A | Q | M_eq | absolute mass error | estimated delta_M |
|---:|---:|---:|---:|---:|---:|
| N=48 | 38.5374382972 | 0.700125283389 | 1.0001070803 | 1.07080e-4 | 1.10342e-4 |
| N=96 | 38.5322754120 | 0.700031470130 | 1.0000274822 | 2.74822e-5 | 3.07216e-5 |
| N=192 | 38.5309892406 | 0.700008081064 | 1.0000076473 | 7.64726e-6 | 1.08401e-5 |

The static area order is **2.00509**. All three mass errors are covered by the
stated combined allowance on this dx=M/32 snapshot. At N=192,
phi_mean=0.0522744836964 and phi_rms=2.47807e-9. The residual is 9.99976e-13.
Quadrature, finite-grid interpolation and the remaining stopping displacement
are still present; the bare table allowance is not the complete error budget.
The measured threshold difference is retained rather than correcting the mass.

| Per frozen checkpoint process | N=48 | N=96 | N=192 |
|---|---:|---:|---:|
| Two-hole reduced binary, seconds | 25.58 | 116.67 | 784.80 |
| Single static member, seconds | 4.89 | 27.16 | 206.13 |

The binary N=192 run was preceded by the shorter resolutions and bounded at
900 seconds. The static cases were bounded at 600 seconds. Default N=96/192
processing therefore costs about **15.0 minutes per reduced binary checkpoint**
and **3.9 minutes per static checkpoint**, plus the short Julia/hash I/O pass.
These are single local observations, not MPI or production-domain benchmarks.
Peak child RSS before Julia was 183 MB for the binary finder and 21 MB for the
static finder. The driver's RSS field is the cumulative maximum over completed
children, so later checkpoint records can include an earlier Julia process.
The binary checkpoint field payload is 54.4 MB; the planning memory estimate is
2.30 GiB, rounded up to a 3 GiB allocation for this reduced case. Full exp-0005
checkpoints need their own dry run and capped timing pilot; no production-size
cost or MPI scaling is established here.

The 17-digit shapes at 48/96/192 points were replayed through the checkpoint
reader. A, Q, scalar statistics and inward expansion agree within 1e-12 absolute;
the outgoing mean-square residual agrees within 1e-20 absolute. Input hashes
confirm that checkpoints, original histories and the static-family table were
unchanged. Missing-time, mismatched-time, negative-radius and NaN seeds are
rejected. Dead/dormant status tests pass. All 24 archived native Julia records
retain exactly the same A, Q, M_eq, error components, scalar diagnostics and table
hash with the extended post-processor.

Compact data and full-precision accepted shapes are in [results/offline](results/offline).
Raw checkpoints, iteration histories, logs and builds remain in
`/private/tmp/ems-offline-f`; generated executables/objects were moved into its
`build-products/` subtree to leave the worktree clean of build products. Rebuild
before using the default executable lookup, or supply `--executable` with that
retained binary path. The second-order claims above concern only the
angular observable on each fixed snapshot. Grid/evolution error for the binary,
late common-horizon behavior, branch identification and remnant settling remain
unresolved for the eventual energy balance.

Files changed: `GNUmakefile`, `harness/EMSRHSurfaceTest.cpp` (optional t=0 checkpoint),
`harness/EMSRHCheckpoint.cpp`, `offline.py`, `check_offline.py`, `postprocess.jl`,
this report, the README pointer and compact results. The radiation test driver's
two stale pilot assertions were aligned with the already checked-in quota 2000
and threshold 1e-10; the production parameter file was not changed.

The batch command was also exercised on checkpoints t=0.25 and 0.375 at N=48,
with the original future rows still present when selecting t=0.25. Both targets
were reached, and the repeated t=0.375 final observables are bit-identical to the
independent N=48 run. These overlapping regression runs took 98.5 and 40.0 s and
are not the isolated cost measurements above.

The existing EMS regression suite passes: reference/echo trumpet and CTT
fixtures/rejections/fingerprints, 26 t=0 grid jobs, four native RH pilots,
restart valid/invalid controls, radiation tensor/analytic and pulse-grid tests,
stationary and production-gauge evolutions, AMR, threshold and parameter checks,
and default-off/explicit-off byte identity against the archived fork. The
existing radiation collector passed with its output directory redirected to the
scratch tree; its assertions were unchanged. Its archived algebra witness was
reused because no radiation algebra changed; no new symbolic verification is
claimed. Compact evidence is in `results/offline/regression-gates.json` and
`regression-summary.json`. There are no remaining running jobs.
