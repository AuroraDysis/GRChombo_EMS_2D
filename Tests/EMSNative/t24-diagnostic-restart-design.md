# T24: diagnostic restart of exp-0024 (submission-worker design only)

This packet is **not submitted** and makes no evolution fix. It addresses the
first numerical loss of finiteness, rather than trying to infer it from a final
NaN dump. The native source is fork `d50743867e48bcb466a60ca16869b07fe1546c27`.
Use the production executable/build and layout: four exclusive Leonardo nodes,
224 MPI ranks × two OpenMP threads, Float64, no fast-math, FMA contraction off.
Keep the supplied runtime parameter file, box size 24, levels 0–12,
`dt_multiplier=0.5`, point transfers, KO, ExperimentalGauge, cleaning, tracking,
tagging and all physical diagnostics. No static or CTT file is present on restart.
The sealed input is `chk/EMS_000152.2d.hdf5`, t=133, 1,500,693,480 bytes,
whole-file SHA-256 `f19f9e3a7ee6d0c9906d1aeb5042488ab0d3570e84aac150061f2186bb7bc6c7`.
Copy the numerical puncture history into the isolated diagnostic output directory;
use the production restart-history validation and remove only duplicate times.
Neither current tracked centres nor driver variables are reinitialized.

## Existing output and exact native overrides

The native checkpoint writer can produce the needed synchronized step-154 state.
The standard plot writer cannot capture an asynchronous mid-subcycle regrid.
The existing Tests-only `t24-observer.hpp` captures *post-regrid box coordinates*
only; it does not provide pre-regrid fields or the first failing RK stage.
The T13 recorder is serial, stops the evolution, and is not the production recorder.

Apply these overrides to an isolated diagnostic directory (replace each placeholder
with the submit worker's actual absolute path; leave every other key unchanged):

```text
restart_file = <sealed-input>/EMS_000152.2d.hdf5
max_steps = 155
stop_time = 1e100
checkpoint_interval = 2
plot_interval = -1
plot_period = -1
hdf5_path = <diagnostic-output>/
data_path = <diagnostic-output>/
ems_data_path = <absent-input>/ABSENT.trumpet
ems_ctt_data_path = <absent-input>/ABSENT.ctt
t2_guard_initial_data_after_t0 = true
nan_check = 1
```

The native guard parameter name is retained from the production file; do not add
an unrecognized alias. The production driver must stop expecting only a step-156
endpoint: seal the completed step-154 checkpoint independently once HDF5 has
closed it, even if step 155 later aborts. Store its per-level clocks, layouts,
partition-independent valid-field digest and whole-file digest. A failed segment
must not emit a successful `run.done` or `chain.done`. Its nonzero exit, snapshot
completion receipts and abort witness are separate artifacts. The production
driver currently restarts in four-step segments; this diagnostic is a distinct
three-step attempt. External finder/post-segment reductions must not conceal the
native failure or prevent collecting the already closed step-154 checkpoint.

## Required read-only hook (to implement and qualify before execution)

Use an opt-in, default-off **diagnostic variant**, not a change to the in-flight
production source. Proposed parameter contract:

```text
ems_regrid_audit = false                    # default
ems_regrid_audit_levels = 8 9 10 11 12
ems_regrid_audit_region_lo = 2239.60 0.0
ems_regrid_audit_region_hi = 2240.40 0.15
ems_regrid_audit_time_window = 135.15 135.19
ems_regrid_audit_fields = chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi
ems_regrid_audit_saved_ghosts = true
ems_regrid_audit_directory = <diagnostic-output>/regrid-audit/
ems_nonfinite_audit = false                 # default
ems_nonfinite_audit_time_window = 133.0 135.625
ems_nonfinite_audit_max_abs = 1e20           # the current NanCheck trigger
ems_safe_nan_abort = false                  # proposed failure-control option
```

For this diagnostic set the three switches true. The field/radius/level/time
settings are parameters, not baked-in physics constants. The audit saves raw
numerical buffers as they stand. It must never fill ghosts, project the fields,
call tracking, retag, alter time steps or read either static input to obtain an
output. Gauge-driver columns are labelled `ExperimentalGauge:B1/B2`.

1. Override the already available `AMRLevel::preRegrid(base,new_grids)` callback
   in the diagnostic variant: it is invoked in Chombo `AMR.cpp:1477–1480` before
   redefining any level. Save selected current buffers, old layouts, proposed
   layouts and each level's **own** full-precision `m_time`, `m_dt`, h, coarse-step
   counter and regrid sequence ID. This is after tagging (and its ghost fills).
   If an input first becomes bad during tagging, the separate phase scan below
   records that fact. Do not describe these asynchronous states as synchronized.
2. Immediately after each `GRAMRLevel::regrid` returns, retain its local raw
   fields and the new layout. Also snapshot at the existing `postRegrid(base)`
   callback (`AMR.cpp:1495–1497`), after the whole hierarchy is reconstructed.
   Distinguish `post_level_regrid` from `post_hierarchy_regrid`. The latter is
   the before/after comparison counterpart of `preRegrid`; the former locates
   interpolation/old-overlap-copy faults. Include both the instantaneous
   coarse-source clock and the receiving fine clock.
3. Rank 0 writes the complete global box lists of levels 0–12, coordinates,
   domain, owner rank, per-level clocks and refinement ratios, once per phase.
   These are global metadata already available without a new collective. Every
   rank writes only its owned FABs intersecting the fixed region. Save all 28
   evolved components as Float64 bit patterns plus box coordinates and the saved
   ghost width. Keep valid and ghost cells explicitly distinguishable. The
   overlapping ghost copies are separate records, not extra physical cells.
   Use per-rank files, with no MPI/HDF5 collective inside a failure path.
4. Compare the *unions* of boxes, not their partition IDs: old minus new and
   new minus old, for each level and for child coverage of parent valid cells.
   Check nesting with the registered eight-parent-cell buffer. Mark the L10
   bad cell and its mirrored cell as retained/new/exposed/removed; save the
   complete 6×6 donor stencil (and its clock) when a donor is involved. A new
   box or a different owner is not by itself a new physical cell.

The nonfinite observer scans all locally valid evolved fields on all levels,
including covered coarse cells, throughout the scan window; ghost results are
reported separately. Scan after restart read, before/after tagging, pre/post
regrid, before/after `postTimeStep` point restriction, at each RK input before
ghost filling, after exchange/time interpolation, before/after the existing
trace/floor projection, at RHS output before/after KO, after each ODE update,
and before/after final-advance projection. The insertion sites are
`GRAMRLevel.cpp:1024–1118` (RHS/update), `:247–257` (restriction), `:380–474`
(regrid), and `EMSBH2DLevel.cpp:130–154,326–345` (projection and existing check).
RHS-buffer witnesses are explicitly labelled as RHS, not evolved-state values.
Track the actual RK evaluation index and stage argument; do not infer a stage
from a rounded progress line. Keep phase extrema, floor counts and negative
chi/lapse counts in bounded memory, including the immediately preceding phase
summary. Emit them at selected regrids, coarse endpoints and the first witness,
not at every fine RK call on every rank. Record negative-metric or large finite
values as observations; do not add an evolution projection. The stage observer
aborts on nonfiniteness; a finite threshold crossing is recorded without ending
the stage audit prematurely. At the native end-of-advance check retain the
existing NaN/Inf/1e20 threshold trigger and label that termination accurately.

A witness includes source/build hashes, input seal, rank/thread, level, phase,
regrid ID, RK index, exact stage/level clocks (decimal and hex), dt/h, valid/ghost
class, box bounds, integer and physical cell coordinates, component, IEEE bits,
`nan`/`inf`/`finite_threshold` class, all 28 cell values and its stored stencil.
Write and close the per-rank witness before nonzero termination. The first
*detected* local phase is established; without synchronizing the entire
subcycling schedule do not claim a unique globally earliest physical time.

## Failure-control proposal and qualification

The current checker uses conditional `omp single` at `NanCheck.hpp:51`, then
`MayDay::Error` inside an OpenMP work region. Propose collecting a stop flag and
a bounded witness without an implicit barrier, joining the worker region, then
writing the witness and directly calling `MPI_Abort(MPI_COMM_WORLD, nonzero)`
from the calling/main thread. In a serial build exit nonzero. Do not use an
Allreduce or Barrier to agree on the failure; other ranks may be blocked in
existing exchanges. Do not rely on `MayDay::Error` to reach MPI_Abort: its default
branch calls `abort()` first (`MayDay.cpp:67–78`), and timer/debugger SIGABRT
handlers may intervene. The proposed switch leaves legacy failure control
available by default; it does not alter any finite-state arithmetic.

Before cluster execution, compare omitted/off/on recorder builds on a small
fixture through initialization, an advance and an actual regrid. Require zero
evolved-bit mismatches; compare defined diagnostics as established in T23.
Inject one NaN under a Tests-only build on a non-master worker and require a
closed witness and nonzero termination of **every** MPI rank. Also test Inf
and the unchanged finite threshold. This MPI test must use an exclusive whole
node under the standing orders. The completed local two-thread collector test
is not the MPI qualification. Use the step-152 production-layout diagnostic
as the subsequent four-node reproduction, not as a long physics run.

## Size, duration and receipts

[t24-recording-size.csv](t24-recording-size.csv) counts the chosen region directly
in step 152. Full intersecting FABs, including saved ghosts, total **86,798,208
bytes per field snapshot** across levels 8–12 (all 28 fields). Selected valid
cells alone total **45,257,856 bytes**. Regrids are separated by dt6=0.013671875;
the window [135.15,135.19] contains three lattice clocks, 135.16015625,
135.173828125 and 135.1875 (the last may not be reached). Three field phases
at all three clocks cost about **781 MB** at the saved layout; budget **1.2 GB**
for movement/halos plus metadata/witnesses, and **1.6 GB** for step 154. Allocate
**4 GB** total incremental output. Record actual sizes; cap the diagnostic
at 8 GB, retaining partial/completion flags if the cap fires. Global coordinate
metadata are sub-MB per phase with 2,704 nonzero-level boxes in step 152.
Do not write whole ROI fields or per-rank extrema rows at every RK evaluation.
Scanning all phases is necessary, but the bounded phase summaries are flushed
only at the stated regrids/endpoints/abort. Otherwise millions of fine-stage
rows across 224 ranks would invalidate the storage budget. A first-witness
stencil and its immediately preceding input/RHS summary are retained.

The controller measures roughly two minutes to the NaN with this exact layout;
record that as the baseline, not a measured cost for the added observer. Give
the diagnostic a 15-minute native phase cap. Preserve all per-rank progress,
abort, scheduler/native exit and resource receipts on termination. A resource,
storage or phase-cap stop is distinct from a numerical abort. The failed local
full replay is not repeated: its footprint exceeded 7.5 GB before the failure
time even though peak native wait4 RSS was 5.119 GB.

## Reading that distinguishes the mechanisms

If negative lapse/metric/chi first appears at point restriction or a newly
interpolated fine cell, and its exact donor arithmetic reproduces it, identify
that operation and stencil; then see whether the first subsequent RHS/floor
event amplifies it. If old valid overlaps change across a same-clock regrid,
check the copy/communication path and nesting before attributing it to gauge.
If all regrid states are unchanged on overlaps and finite, but a particular
RK/RHS phase fails in an old valid cell, the immediate cause is in that
advance, with regrid at most an antecedent. If nonfiniteness predates regrid,
the last regrid is not the initiating event. A normal large-volume constraint
norm does not exclude an inner local failure. No transfer, floor, gauge,
refinement or puncture-coincidence rule is changed by this diagnostic.
