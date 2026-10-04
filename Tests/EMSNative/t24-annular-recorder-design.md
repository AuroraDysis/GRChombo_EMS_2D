# T24 annular RK-stage recorder: build cancelled; historical design only

The controller cancelled this build after the exp-0026 G/D interventions.
No implementation was started. This retained design and size estimate are
not an active submission plan or a qualified recorder.

This replaces the earlier fixed box around the printed L10 cell. It is an
historical implementation contract. The cancelled build produced no installed
hook, tested MPI recorder or evolution remedy. The abort hook is a separate
implemented option.

The regions follow the two **current numerical tracked punctures**, retaining
their full precision coordinates. A cell is selected when its distance in the
cartoon half-plane from either centre is in the configured closed interval.
Use the union, write each owning valid cell once, and give it a hole-membership
bit mask. Include covered coarse valid cells: the original failure is first
reported on a covered L10 advance. Ghost records, if requested, must be
separately labelled by owning box and coverage; they cannot enter a valid-cell
checkerboard estimate. No puncture separation denominator, static seed, static
field, lapse clipping or modification of a right-hand side is introduced.

The proposed parameters are `ems_stage_record=false`, `ems_stage_record_levels`,
`ems_stage_record_time_window`, `ems_stage_record_step_stride=1`,
`ems_stage_record_radii=0.04 0.09`, `ems_stage_record_fields`,
`ems_stage_record_prefix`, and `ems_stage_record_max_bytes`. Times, levels, radii,
cadence, fields and quota are parameters; there are no embedded cell locations
or merger times. An enabled exhausted output quota writes a closed
`RECORDING_INCOMPLETE` receipt and terminates the diagnostic job with a nonzero
exit rather than silently producing a purported complete record.

Capture the stage input before RHS projections, then perform read-only
nonfinite scans on the input, completed RHS and each ODE update. A failure
witness names phase, level, native cell and box, rank, field, stage ordinal,
fine-step ordinal, RK evaluation time, step start, dt and ODE update weight.
These distinguish the two RK evaluations at half time and do not label an
intermediate state as an accepted physical-time solution. Save a whole selected
annulus at stage input; additionally save the corresponding small witness
neighbourhood at the first bad RHS/update. The "first" is the first sampled
phase on each rank, not a globally ordered arithmetic operation across MPI.
The first nonfinite scan covers all evolved components on the selected levels,
including cells outside the output annulus; it records an outside-region flag.
Derived S failures are labelled separately from nonfinite evolved fields.

Store beta1/2, driver B1/2, Gamma1/2, chi, lapse, h11/h12/h22/hww, K, Theta and
`S=c_Gamma*(H11+4 H22/3)` (equal to `H22+0.75 H11` for this run, with the
coefficient read from its gauge parameters), as binary Float64 values. The
inverse metric is formed from the **current numerical metric**. Also record
`S_effective=c_Gamma*(trace(H)+lambda_max(H)/3)`, or compute it offline from the
stored h fields, to retain the h12 coupling. B has a gauge-class label and is
never compared to a different gauge's driver. Store global integer coordinates,
dx, centres and the complete component list, so cell parity and staggering
cannot be lost to decimal rendering or interpolated profiles.

Maintain an opt-in per-level regrid generation counter and the generation at
the last completed fine step. Every stage header states whether a regrid has
intervened since that step, and whether the box list actually changed. Emit
pre/post-regrid coordinate box lists and current tracked centres inside the
same time window, even when the box count is unchanged. The proposed hook
locations are `GRAMRLevel::regrid`, `advance`, `evalRHS`, and `updateODE`, with
an EMS diagnostic callback for the tracked-centre region. Existing
`m_rk_stage` alone is insufficient: it is advanced on the point-transfer branch
and cannot serve as a gauge/transfer-independent production stage counter.
Do not reuse T13's serial early-stop recorder.

After the BoxLoops workers join, each MPI rank writes only its own chunks and
receipts, with an atomic final rename and no recording collectives. A node
reducer verifies every rank's sequence, field count and checksum against the
saved box lists. It must recognize ranks with zero selected cells and an
aborted partial final step. An incomplete capture cannot exclude the proposed
signature. The disabled branch allocates no recording buffer, scans no fields,
updates no recording counter and calls the original execution path.

Before use, qualify omitted versus false versus true on at least two MPI ranks
and two OpenMP threads, with a real regrid and a restart, comparing all evolved
checkpoint/plot values and all defined valid diagnostic values bit for bit.
Exclude only the already demonstrated uninitialized plot-diagnostic ghost
slots, with their counts reported. A second tests-only input contains a
checkerboard mode with a known per-step factor and verifies cell parity,
stage identity, growth factor and hole labels in the emitted record. An injected
worker/rank nonfinite value verifies the joined abort and closed record. None
of these MPI qualifications has been claimed for this designed hook.

The offline reader should report the signed adjacent-cell alternation of beta2
and Gamma2 on each unchanging support, the native Fourier projection
`sum((-1)^(i+j)*u)` together with directional checkerboards, RMS and peak
amplitudes, and consecutive-step amplitude ratios. Compare these with the
frozen-operator growing eigenmode and its predicted factor; a fit over tens of
fine steps must name its interval and uncertainty. Report sampling support,
regrid generations, lapse/chi floor sets, h positivity and S, and retain the
unprojected raw data. A changing ROI or box support cannot by itself be called
exponential growth. No static solution is subtracted.

## Storage and bounded diagnostic window

[t24-annulus-size.py](t24-annulus-size.py) and
[t24-annulus-size.csv](t24-annulus-size.csv) count the exact sealed step-152 box
support, using its numerical tracked centres. A packed cell record uses 120
bytes for 15 doubles and 12 bytes for two int32 coordinates and a uint32 hole
mask. Per-stage headers and box lists add a small amount.

| Level | Native valid annular cells | Payload per stage | Four-stage payload per coarse step |
| --- | ---: | ---: | ---: |
| 11 | 27,972 | 3,692,304 B | 30,247,354,368 B |
| 12 | 46,442 | 6,130,344 B | 100,439,556,096 B |

The sum is **130.687 GB**, uncompressed, at stride 1; it exceeds the inherited
40 GB diagnostic quota. Shapes/counts may change under tracking/regridding,
so these are measured estimates on the retained layout, not hard bounds for
future output. A payload with full RHS and update arrays at every stage would
multiply this estimate and is not part of the proposed default capture.

For a diagnostic repeat with the original dt, use levels 10–12 (to retain the
original covered-L10 witness) and a declared window such as
`135.146484375 <= t <= 135.1875` M_i, with stride 1. On levels 11–12 this
0.041015625 interval contains 96/192 fine advances and costs **6.126 GB** for
the stage-input payload, plus endpoint-inclusive RK records, L10 data, layouts
and first-failure neighbourhoods. It retains tens of consecutive fine steps
before the observed event. Keep a global 40 GB quota and verify actual sizes;
do not record a whole coarse interval merely because its coarse step is 155.
For half dt the same physical window doubles stage counts and payload, about
12.252 GB on levels 11–12. The submit worker must calculate the selected L10
and header allowance at runtime and register it before starting.

Existing checkpoint output can save coarse step 154. At dt_multiplier=0.25,
the same physical t=134.75 is step **156**, not 154; the old failure time falls
within new coarse step 157. Saved clocks and comparisons must be in physical
time. Nothing in this design changes the production transfer, gauge, KO,
initialization or evolution equations.
