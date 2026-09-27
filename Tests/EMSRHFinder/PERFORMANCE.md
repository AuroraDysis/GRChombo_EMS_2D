# Finder/interpolation performance port — 2026-09-27

Implemented without committing or using SSH. Reference: fork `e478b5fc08cd038c45f8bd1499b22a884d851a1a`.
All measurements here used the same local macOS arm64 machine, GCC 16.2, `-O3`,
and **serial Chombo: one MPI rank**. MPI correctness and scaling remain a cluster gate.
Raw builds, parameter snapshots, output files and logs are in
`/private/tmp/ems-finder-perf/`; compact results are in [perf-results.json](perf-results.json).

## Changes and decisions

| Change | Status | Reason / scope |
|---|---|---|
| Trinary `parallel-answers.patch` | Ported verbatim, then extended with the layout hoist | The original interpolator files compared byte-for-byte. Point-local scratch and disjoint answer writes; arithmetic/order within each point unchanged. Verbosity >=2 stays serial. |
| `layoutIterator()` outside the parallel loop | Ported | One evolution/diagnostic iterator per level per call removes refcount traffic inside the point loop. The identity gates below pass. |
| One Lagrange object per thread | Rejected for this port | The object holds a reference to a particular level's `InterpSource`; answers can cross levels. Its stencil cache also grows and is searched linearly. A safe, bounded reuse policy would be a separate change, not a simple hoist. |
| Trinary `root-query.patch` | Adapted | Root queries; all ranks answer for owned boxes; one result broadcast. One descriptor list defines all 23 value/derivative fields, including Ex, Ey and phi, and both query and scatter use its derived size. No hard-coded component count. |
| Empty MPI buffers | Fixed | Query/answer exchanges use `.data()` instead of forming `&vector[0]` on zero-point participants. MPI calls remain outside OpenMP. |
| Trinary restart centre read | Adapted and extended | Restore the centre **and matching full angular shape**, plus dead status, on root and broadcast before pruning. Trinary's centre-only read still reinitializes radii. |
| Pilot | Changed | Level-3 call cadence retained; chase quota 400; max box 16, block factor 16. See measured layout and limitations below. |
| Manual AMR cap, diagnostic-cadence patches, box-search.patch | Not ported | Explicit archive-only scope; box-search had an unaccepted last-digit gate and cadence had a Weyl correctness failure. No changes to equations, tagging thresholds, precision, or diagnostic scheduling. |

Production files changed: `Source/AMRInterpolator/AMRInterpolator.impl.hpp`,
`Source/RHFinder/RHUnion.hpp`, `Examples/EMS/params-radiation-pilot.txt`.
Test/report files: this report, `README.md`, `perf-results.json`, `compare_perf.py`, `run_perf.py`,
`harness/EMSRHRestartTest.cpp`, the RH test `GNUmakefile`, and the radiation
`README.md` / `run_local.py` (their old pilot assertions/documentation were stale).

## Identity and tests

Both reference executables were clean builds of the archived reference revision.
Candidate application and test objects were also rebuilt from scratch, including
PunctureTracker; no weak template definition was inherited from an old object.
After the final cleanup, another clean rebuild reproduced the final-box-16
baseline exactly (26 files / 27,697,770 numeric values), passed the frozen
full-precision shape replay and both restart unit controls.

* Frozen harness: N_theta=96, dx=1/32, initial sphere 0.64, chase speed .125,
  quota 400 and threshold 1e-12. All four baseline/candidate × OMP 1/8 runs found
  the surface after **787 calls**. Each comparison checks 361,700 numeric values,
  including the t=0 plot. The complete RH trajectory and shape files are byte-identical.
* Reduced binary pilot: N1/N2=192/96, L=96, centre (48,0), separation 32,
  certified trumpet and degree-28 CTT input, max_level=3, RH_level=3,
  extraction level 1 at radii 20/24/28. Three coarse steps through t=.375,
  24 finder calls, regrids and checkpoints 0/1/2/3. The timing matrix used max box 32.
  Baseline/candidate and OMP 1/8 comparisons have **zero numeric differences**.
* Final max-box-16 repeat, with raw Weyl output enabled: baseline OMP 1,
  candidate OMP 1 and candidate OMP 8 agree over **26 files / 27,697,770 numeric values**.
  This includes all checkpoint datasets and attributes, evolution ghost data,
  RH histories/shapes, native radiation CSV, legacy GW/MQ modes, raw Weyl and
  the three legacy scalar diagnostics. HDF5 creation timestamps are excluded
  by comparing its contents rather than its container bytes.
* Restart from step 2 at t=.25, with deliberately retained future rows through
  t=.375: both active surfaces restore 96 saved radii and their saved centres;
  the t=70 surface stays dormant. The old prefix is retained exactly, future
  RH rows are pruned, and continuation reaches .375 without duplicate RH times.
  OMP 1/8 restart outputs have zero numeric differences.
* The new restart unit test checks a displaced, nonspherical surface, ghost
  reconstruction, the dormant surface, saved dead status and future pruning.
  A negative-radius input is rejected before its history is truncated.
* Existing fork regressions pass: reference and echo B/E trumpet fixtures,
  CTT fixtures/rejection controls, all 26 grid regressions, all four native RH
  pilot cases, radiation analytic/tensor tests, pulse grids, stationary and
  production-gauge evolutions, threshold tests, AMR/regrid smoke, parameter
  rejection tests, and default-off identity against the earlier fork revision.
  Fingerprints remain `e1107d452abc4e4b`, `e68444ed5f33ba7b`,
  `024d0092dba552bc`, `b3b32ee0d0a729dc`.
  The radiation collector returned PASS. Its unchanged symbolic witness was
  reused from the archive; no new symbolic derivation is claimed.

`compare_perf.py` checks complete file/row coverage, finite values and absolute
error <=1e-12, retaining repeated-time rows. It does not silently intersect
output times. `--mpi` only permits differing checkpoint `Processors` datasets.

These are **same-parameter code identity gates**. Reducing the quota changes
intermediate finder trajectories/output opportunities and is not claimed to
reproduce the old 4000-quota history at every physical time.

## Local timings

Seconds; single observations, sequential benchmark runs, no performance extrapolation
to MPI. `RH` is the sum of the existing update messages, rounded to milliseconds
per call. Interpolation timers are Chombo wall times, summed across call sites
from the last timer report (not summed again across repeated reports).

| Case | OMP | Code | Whole process | RH updates | interp | calculateAnswers |
|---|---:|---|---:|---:|---:|---:|
| Frozen | 1 | baseline | 117.215 | 116.674 | 100.838 | 98.672 |
| Frozen | 1 | candidate | 113.375 | 112.921 | 97.066 | 94.866 |
| Frozen | 8 | baseline | 116.389 | 116.247 | 99.909 | 97.703 |
| Frozen | 8 | candidate | 66.660 | 66.468 | 50.697 | 48.110 |
| Three-step binary, box 32 | 1 | baseline | 27.458 | 12.311 | 10.535 | 9.879 |
| Three-step binary, box 32 | 1 | candidate | 23.971 | 10.389 | 8.802 | 8.256 |
| Three-step binary, box 32 | 8 | baseline | 52.376 | 19.460 | 10.396 | 9.727 |
| Three-step binary, box 32 | 8 | candidate | 35.692 | 9.838 | 2.985 | 2.404 |

At eight threads the frozen RH workload is 1.75x faster; the short binary's RH
updates are 1.98x faster, interpolation 3.48x faster, and process time 1.47x faster.
Eight threads still lose to one thread for this small whole binary run. The
local result establishes a benefit at fixed thread count, not ideal scaling.
Root-query's removal of duplicate **MPI** queries cannot be timed here.
Later validation jobs overlapped; their durations are not controlled benchmarks.

For comparison, rereading exp-0004 pout logs gives summed RH fractions
72.26% at n32 and 62.46% at n64, consistent with the supplied bottleneck evidence.

## Restart finding and limits

Before this change, setup always filled the centre/radii from parameters,
then pruned rows beyond the checkpoint time and opened files in append mode.
Neither a checkpoint nor the retained RH output was used to restore geometry.
The first post-step callback lazily constructs the union; it receives Chombo's
restart time even when it first runs on a finer level. Start times are absolute
simulation times, not offsets from restart.

The new reader takes the last centre and shape rows at or before the checkpoint,
requires matching times and unchanged point counts, validates finite positive
radii, restores ghosts and broadcasts the same state to every rank. It then
uses the existing pruning/appending path. An active surface with missing or
inconsistent history fails explicitly; copy `rh_surf_*.dat` **and** `rh_f*.dat`
with each checkpoint. A not-yet-started surface may have no rows and retains its
parameter seed. Dead surfaces are not revived. Active solver regimes are
reclassified using fresh fields on the next update.

The legacy RH writer uses nine significant digits. Restart now continues from
that saved geometry, but **not** from the original full-precision in-memory
shape. This change leaves output precision intact so the performance gates are
not confused with a serialization change. Exact uninterrupted/restarted RH
identity is therefore not promised; a dedicated full-precision checkpoint record
would be needed for that stronger contract. The evolved grid does not depend on RH.

## Quota and box choice

`RH_level = 3 3 3` controls call cadence: .0625 M for the full pilot.
`RH_time_step_freq = 400 400 400` is the **per-call chase quota**, ten times below
4000. The frozen test demonstrates convergence over successive calls at 1e-12.
Its deliberately displaced seed took 787 calls: at production cadence that would
be about 49 M. This is not a moving-binary tracking certificate. The short binary
runs do not establish strict FOUND status for the individual surfaces; tracking
lag, individual-horizon loss and common-horizon activation at t=70 must be watched
in the cluster pilot. One non-converging active surface still keeps a call at its
full quota; the solver's stopping/chasing rules were intentionally preserved.

Full production-domain initialization, with the actual CTT data, all seven levels
and the radiation refinement floors, produced:

| Level | Valid cells | max box 64 | max box 32 | selected max box 16 |
|---:|---:|---:|---:|---:|
| 0 | 294912 | 72 | 288 | 1152 |
| 1 | 41472 | 18 | 50 | 162 |
| 2 | 147968 | 50 | 162 | 578 |
| 3 | 12288 | 4 | 16 | 48 |
| 4 | 35840 | 16 | 48 | 140 |
| 5 | 65536 | 16 | 64 | 256 |
| 6 | 82944 | 36 | 100 | 324 |

All use block factor 16. The covered cell sets and all valid initial field values
agree exactly after matching coordinates across decompositions. Smaller boxes
here split the same mesh; no refinement cap was applied.

Use **max_box_size=16, block_factor=16, about 256 MPI ranks** as the starting
configuration. Levels 5/6 dominate subcycled updates and have 256/324 boxes.
At 512 ranks level 6 can occupy at most 324 ranks (63%), and level 5 at most half;
coarser sparse levels also leave ranks idle. Thus **512-rank efficiency is not
established and should not be assumed**. Measure 128/256/512 ranks inside the
smoke before choosing production topology. Smaller-than-16 boxes would require
changing block factor and revalidating mesh coverage, with substantial ghost and
communication overhead; that change is deferred. Even size-16 boxes have 22x22
stored cells for 16x16 valid cells with three ghosts, so decomposition has a cost.

## Read-only radiation audit

The generic Trinary hazard remains in shared infrastructure:
`GRAMRLevel::regrid` reallocates diagnostics without recomputing them;
`postRegrid` only propagates restart time; `refresh()` fills ghosts rather than
computing Weyl. A fine-level callback happens before its coarser post-step callback,
and interpolation searches down to level 0 without a minimum-level/time assertion.
Consequently a query or diagnostic ghost stencil that reaches an unprepared
coarser level can read stale/uninitialized Weyl data. No cadence patch was applied.

There are relevant protections in this fork that Trinary lacked: native radiation
activation enforces refinement through every requested extraction sphere, plus the
pilot's radius-132 level-2 wave zone. Main also computes Weyl on **all levels**
before the initial extraction, including on restart. Pilot radii 50/75/100 are
inside that buffered region. This mitigates the specific missing-coarse-coverage
failure; it is not a general timestamp guarantee for arbitrary radii or AMR interfaces.

A separate max-box-16 reduced pilot with three coarse steps, actual regridding,
no initial plot and raw Weyl output had 6 files / 1170 raw samples per build.
Every sample was finite and raw Im was exactly zero. Baseline OMP 1 and candidate
OMP 8 matched exactly, including the real data, modes and checkpoints. This is
bounded positive evidence for the tested layout, not proof that all future
coarse/fine ghost sources have the right time. Keep the raw-Im and stencil/time
coverage checks in the cluster smoke; a successful RH/checkpoint gate alone
cannot certify Weyl correctness.

## Reproduction and cluster gate

Local build settings are in `/Users/auroradysis/Workspace/EMS-deps/BUILD.md`.
Use separate clean source trees for baseline and candidate; clean every application
object before building. Run long commands in the background and poll with backoff.
The drivers cap each process at 240 s. The complete regression batch had an
1800-second outer cap; its individual evolutions were under five minutes.
No >30-minute run or SSH/HPC operation was launched.

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
# In each source tree, independently:
make -C Examples/EMS clean DIM=2
make -C Examples/EMS all DIM=2 -j 8
make -C Tests/EMSRHFinder clean DIM=2
make -C Tests/EMSRHFinder all DIM=2 -j 8

python3 Tests/EMSRHFinder/run_perf.py /private/tmp/ems-perf-repeat BASELINE_TREE identity
python3 Tests/EMSRHFinder/run_perf.py /private/tmp/ems-perf-repeat BASELINE_TREE restart
python3 Tests/EMSRHFinder/run_perf.py /private/tmp/ems-perf-repeat BASELINE_TREE boxes
python3 Tests/EMSRHFinder/run_perf.py /private/tmp/ems-perf-repeat BASELINE_TREE final
uv run --with numpy --with h5py python Tests/EMSRHFinder/compare_perf.py REFERENCE CANDIDATE
# Run EMSRHRestartTest2d.*.ex in an empty directory; "invalid" must fail.
```

For the **cluster smoke, not executed here**, build both e478b5f and this candidate
cleanly against the same MPI Chombo, compiler and flags. Inside an existing
allocation, use the generated reduced parameter snapshot
`runs/final-candidate-8/params.txt` from the local evidence (or regenerate with
`run_perf.py final`), replacing its two local input paths with the staged certified
trumpet and CTT paths. It has max_box_size=16/block_factor=16, N=192x96, L=96,
max_level=3, RH_level=3, three coarse steps, checkpoint_interval=1,
plot_interval=-1 and raw extraction on. All compared runs must read identical
physics/input parameters and hashes. Preserve the site's already validated MPI
transport configuration.

The following is a bounded **first pass**: quota 4 on both builds to exercise all
23 fields, zero-point participants and multiple regrids without committing to the
unoptimized 256/512-rank cost. `SMOKE`, `BASE_EXE`, `CAND_EXE`, `PROFILE`, and `CTT`
are absolute paths. Use `RANKS=256 THREADS=2` (512 allocated cores), or
`RANKS=512 THREADS=1`. The map/bind syntax follows the official
[Open MPI mpirun manual](https://github.com/open-mpi/ompi/blob/main/docs/man-openmpi/man1/mpirun.1.rst).

```bash
export OMP_NUM_THREADS="$THREADS" CH_TIMER=1
for code in baseline candidate; do
    exe="$BASE_EXE"
    [ "$code" = candidate ] && exe="$CAND_EXE"
    for ranks in 1 "$RANKS"; do
        directory="$SMOKE/$code-$ranks"
        mkdir -p "$directory"/chk "$directory"/plt "$directory"/data/extraction
        (
            cd "$directory"
            timeout -k 10s 240s mpirun -np "$ranks" \
                --map-by "slot:PE=$THREADS" --bind-to core --report-bindings \
                -x OMP_NUM_THREADS -x CH_TIMER "$exe" "$SMOKE/params.txt" \
                "ems_data_path=$PROFILE" "ems_ctt_data_path=$CTT" \
                'RH_time_step_freq=4 4 4' > run.log 2>&1
        ) || exit 1
    done
done
uv run --with numpy --with h5py python Tests/EMSRHFinder/compare_perf.py \
    "$SMOKE/baseline-1" "$SMOKE/candidate-1"
uv run --with numpy --with h5py python Tests/EMSRHFinder/compare_perf.py \
    "$SMOKE/baseline-$RANKS" "$SMOKE/candidate-$RANKS"
uv run --with numpy --with h5py python Tests/EMSRHFinder/compare_perf.py \
    "$SMOKE/candidate-1" "$SMOKE/candidate-$RANKS" --mpi
```

Also require: no hang when ranks have zero answers or queries; finite all outputs;
raw `Weyl4_Im == 0`; unchanged RH iteration/row counts for equal parameters;
checkpoint fields/boxes/offsets/attributes equal to <=1e-12 (only Processors may
differ across rank counts). Repeat OMP 1 versus the production thread count on
the same MPI layout. Run the new restart unit with MPI; copy checkpoint 2 and
both RH histories, including future rows, then restart through step 3 at 1 and N
ranks and compare the complete resulting histories/fields with the same gate.

Then repeat with the actual quota **400 400 400** in fresh directories. Use the
same `RH_start_times=0 0 0` override on both builds for an additional bounded
common-surface/non-convergence stress test. For the quota-400 repeat, use the
short measurement to project cost: under five minutes run; 5–30 minutes use a
bounded pilot and explicit cap; a projection beyond 30 minutes is BLOCKED.
Record RH and interpolation wall timers, actual box/rank ownership, transport,
OMP binding and full input/executable hashes. Do not infer production scaling
or moving-surface convergence from the one-rank local results.
