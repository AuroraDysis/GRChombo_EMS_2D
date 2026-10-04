#!/usr/bin/env python3
"""Seal completed T24 read-only evidence; no job launch, no physics changes."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path('/private/tmp/ems-t24')


def rows(name):
    return list(csv.DictReader((HERE/name).open()))


def table(headers, data):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join(
        '| '+' | '.join(map(str,row))+' |\n' for row in data)


def main():
    audit = json.loads((HERE/'t24-offline-receipt.json').read_text())
    native = json.loads((ROOT/'restart152-rss/native.receipt.json').read_text())
    probe = json.loads((HERE/'t24-nan-barrier-qualification.json').read_text())
    assert audit['returncode']==0 and audit['advances']==0
    assert len(audit['witnesses'])==2 and len(set(w['dump_sha256'] for w in audit['witnesses']))==1
    assert native['gate_reason']=='AGGREGATE_FOOTPRINT_GATE_7.5GB' and not native['complete_to_step156']
    extrema = rows('t24-checkpoint-extrema.csv')
    nearest = rows('t24-checkpoint-nearest.csv')
    def e(level,region,name):
        return next(r for r in extrema if int(r['level'])==level and r['region']==region and r['field']==name)
    def fmt(x):
        return f'{float(x):.6g}'
    def range_(r):
        return fmt(r['minimum'])+' … '+fmt(r['maximum'])
    pts = [r for r in nearest if r['region']=='failing_cell']
    nearest_table = table(['L','χ','lapse','K','h11 / h12 / h22 / hww','A11 / A12 / A22 / Aww','β1 / β2','Θ'],
        [[r['level'],fmt(r['chi']),fmt(r['lapse']),fmt(r['K']),
          ' / '.join(fmt(r[n]) for n in ('h11','h12','h22','hww')),
          ' / '.join(fmt(r[n]) for n in ('A11','A12','A22','Aww')),
          ' / '.join(fmt(r[n]) for n in ('shift1','shift2')),fmt(r['Theta'])] for r in pts])
    hole_table = table(['L','χ range, P1','χ range, P2','lapse minimum P1 / P2','negative lapse P1 / P2','lapse floor P1 / P2','sampled cells per hole'],
        [[l,range_(e(l,'puncture1','chi')),range_(e(l,'puncture2','chi')),
          ' / '.join(fmt(e(l,p,'lapse')['minimum']) for p in ('puncture1','puncture2')),
          ' / '.join(e(l,p,'lapse')['nonpositive'] for p in ('puncture1','puncture2')),
          ' / '.join(e(l,p,'lapse')['floor_cells'] for p in ('puncture1','puncture2')),
          e(l,'puncture1','lapse')['valid_cells']] for l in range(8,13)])
    floor_table = table(['L','χ minimum, all valid','lapse minimum, all valid','negative lapse, all valid','lapse floor, active valid','lapse floor, all valid'],
        [[l,fmt(e(l,'all_valid','chi')['minimum']),fmt(e(l,'all_valid','lapse')['minimum']),
          e(l,'all_valid','lapse')['nonpositive'],e(l,'active_valid','lapse')['floor_cells'],
          e(l,'all_valid','lapse')['floor_cells']] for l in range(8,13)])
    probe_table = table(['Case','Native return','Runner return','RSS MB','Meaning'],
        [[r['case'],r['native_returncode'],r['runner_returncode'],f"{r['peak_RSS_bytes']/1e6:.2f}",
          ('SIGABRT locally; no hang reproduced' if r['case'].startswith('legacy') else
           'closed witness, unchanged state, nonzero' if r['case']=='safe' else 'finite state unchanged, zero')]
         for r in probe['cases']])
    report = f'''# T24: merger blow-up and missing termination — completed local evidence

**READY-EXCEPT: a numerical blow-up is established, but its initiating operation
and the exact cluster abort wait remain unestablished.** No evolution fix is
applied. The full local replay is ended, not pending, and will not be retried.
The submission-worker diagnostic design is [t24-diagnostic-restart-design.md](t24-diagnostic-restart-design.md).

All source line references below refer to fork
`d50743867e48bcb466a60ca16869b07fe1546c27`, archived in
`/private/tmp/ems-t24/source`, and its strict private Chombo dependency. The
dirty worktree's later gauge changes were not used. Both static paths were
absent in the restart, with `t2_guard_initial_data_after_t0=true`.
The standing KO/gauge/equations/transfers and initial-data rules are intact.

The sealed checkpoint is step 152, t=133 M_i, **1,500,693,480 bytes**. Whole SHA-256
`{audit['checkpoint_sha256']}` matches the previous audit; the 28-field
partition-independent valid-state seal is
`7146c763475ad64b73a62f19ab2e17f5e091a3da6891698ccda8a18c03b389cc`.
Both supplied manifests now verify: **1,855 ev1** and **692 ev2** files.
[t24-checkpoint-read.py](t24-checkpoint-read.py) regenerates the checkpoint/log
tables directly from the collected root, using no static solution. Its own
receipt is [t24-offline-receipt.json](t24-offline-receipt.json): return 0,
{audit['elapsed_s']:.2f} s, **{audit['peak_RSS_bytes']:,} bytes peak RSS**.

## A1. The saved state and the actual nonfinite witness

The L10 coordinate (1310826,27) is at relative (x,y)=(0.1820068359375,
0.0469970703125) M_i, r=0.187976629 M_i, **0.0643615 M_i from puncture 2**
at t=133. The tracked positions then are x=±0.2259805 relative to the midpoint;
the tiny tracked y offsets are retained in the CSVs. The point and its mirror,
and both punctures, are covered by **every level 8–12** at step 152. In
particular the failing L10 cell is covered coarse data, inside the actual
level-11 and level-12 refinement footprint; it is not a saved fine-face hole.
[t24-cell-coverage.csv](t24-cell-coverage.csv) gives the native cells/boxes.

The following are the nearest *native valid cell* on each level to the bad
point. Different levels have different cell-centre offsets; these values are
not an interpolation to a common point. The L10 row is the exact requested cell.

{nearest_table}

All corresponding fields are finite. The mirrored rows show the expected
opposite signs of off-diagonal/longitudinal odd components, with close but not
bit-identical even fields. The hole-centred native cells have lapse 1e-12 on
every L8–12, |K|≈0.81–1.15, |Θ|≈0.51–0.68, h11≈3.87–4.36 and
|A11|≈1.85–3.81; both holes have the same qualitative state. Full 28-field
values are in [t24-checkpoint-nearest.csv](t24-checkpoint-nearest.csv), and
366 nearby stored valid-cell samples in
[t24-checkpoint-neighbourhood.csv](t24-checkpoint-neighbourhood.csv). These
samples are clipped at their owning boxes/axis; the region statistics below
include all valid boxes intersecting each specified window.

For each hole the region is the stored portion of a radius-0.075 M_i circle,
y≥0. L12 does not cover that entire circle: counts quantify the retained support.
Counts below include covered coarse valid cells, not ghosts.

{hole_table}

For the bad-point window |x−x_bad|,|y−y_bad|≤6h10=0.01025390625 M_i,
χ ranges from 0.00154 to 0.00292 across L8–12; all are positive and above
the floor. L12 lapse ranges from 1e-12 to 8.65e-5, with **116 floor cells
of 2,304**. L10 has **13 negative lapse cells of 169** in this window, and
L11 has **14 of 576**, in covered coarse data. Metric components there remain
finite (h11≈3.83–4.29, h22≈0.489–0.517, hww≈0.504–0.538); |A11|≤2.853,
|K|≤0.206, |Θ|≤0.141, |β1|≤0.01893, |β2|≤0.00227. The full
minimum/maximum/**RMS**, location and floor tables for both holes and the
bad/mirrored windows are [t24-checkpoint-extrema.csv](t24-checkpoint-extrema.csv).
The fixed windows are coordinate-volume cell summaries, not horizon masks.

Both attempts' rank-96 logs contain **the identical printed dump**, SHA-256
`{audit['witnesses'][0]['dump_sha256']}`. This is a detected numerical
failure, not merely a silent regrid: χ=lapse=1e-12, h11=−1.26742e12,
h22=−3.53872e11, hww=−7.98873e10, K and A are NaN, Θ and Γ1 are Inf,
β2=1.85437e11, φ=−1.0221e7 and Pi is Inf. Printed precision is six
significant digits; identical text is not a proof of identical complete
pre-failure states. The check is `EMSBH2DLevel.cpp:151–154`, **valid cells**,
after the final RK trace/floor projection. It does not identify the first
nonfinite RK stage or whether a field already failed before regrid.

## A2. `min_chi.dat`, actual positivity, and the floor

**The filename is misleading: it records min(mod_F), not min(χ).** At
`EMSBH2DLevel.cpp:608` the reducer is
`AMRReductions<VariableType::diagnostic>`; line **641** passes `c_chi=0`
to that reducer. `DiagnosticVariables.hpp:12` assigns diagnostic index 0
to `c_mod_F`; `AMRReductions.impl.hpp:25–35,70–80` selects the diagnostic
buffers and that component. The Lorentz scalar at
`Source/Cartoon/EMSCartoonLorentzScalars.hpp:86–112,362` is the code's
`2(B²−E²)` contraction. An electric-dominated state makes it negative.
The −0.019495763002 entry at t=0.875, and subsequent negative values,
therefore do **not** date a χ positivity failure. Hamiltonian/momentum norms
use the correct diagnostic indices; this component mix-up is specific to
the purported χ minimum. A separate evolution reducer would be needed to
repair that output; no repair is applied here.

Every saved valid χ on **all 13 levels** is positive and >1e-12; the global
minimum is **3.976130075854864e-6**, L12 cell (5242349,0) near the left
puncture. On L8–12 saved ghost copies χ is also positive and above the floor.
There are **no negative or floor χ cells** in this checkpoint. There are
**1,002 negative lapse valid cells**, all covered by a finer level; all active
valid cells have lapse≥1e-12. L12 has 29,762/73,984 valid cells at that floor.

{floor_table}

Ghosts are recorded separately, with duplicates: for example L12 has 60
negative saved lapse ghosts, minimum **−6.4982174632904995e-6** at
(5243281,136), beyond its y=135 valid face. These are not 60 physical
negative-lapse cells. L11's ghost minimum is −2.11526e-6. No negative
valid χ claim can be made from these ghost/lapse observations.

`PositiveChiAndAlpha.hpp:28–37` **writes** max(χ,min_chi) and
max(lapse,min_lapse) back into the state. It runs immediately before RHS
construction (`EMSBH2DLevel.cpp:326–345`) and after full RK advance
(`:130–147`). It is not a clamp on denominators alone. `postTimeStep`
then performs point restriction (`GRAMRLevel.cpp:247–257`) **without
another floor projection**; the coarse state can therefore be below the
lapse floor when checkpointed. Newly interpolated fine/ghost data likewise
have no positivity guarantee before the next prescribed projection.

[t24-point-restriction.csv](t24-point-restriction.csv) replays the native
36-point, anchored, signed-weight arithmetic from saved numerical children.
All **14 requested L10 fields**, and the L8–11 lapse-minimum witnesses,
match their stored coarse values **bit for bit**. Most decisively, L11
lapse at (2621639,41)=**−1.3426383648371819e-6** comes from **36 positive
L12 values**, range [2.818563421019781e-6,2.1094440683440722e-4]. Thus
point restriction producing negative covered coarse lapse is established;
its causal connection to the later blow-up is not.

If χ becomes negative in a later RK input/new cell, the next prescribed
floor replaces it by 1e-12 and χ-inverse factors can be **1e12** (and
χ^(-3/2) **1e18**); this prevents division by zero but neither bounds
derivatives nor ensures a regular metric. Such factors appear in
`CCZ4Cartoon.impl.hpp:55,393,431–433,532`. At the saved bad cell 1/χ is
about 453; at the final dump it is 1e12. This huge change may be a cause,
an amplifier or an effect; the two endpoint observations do not choose.

To date onset from retained checkpoints, inspect the numerical t=0 state
and steps **4,8,…,148,152** with active/covered/ghost classes (step 4 is
t=3.5, so it cannot settle t=0.875). Step **154** is needed to bracket the
final failure; it is not supplied because production checkpoints every four
steps. Steps 148 and 152 alone give the nearest available pre-failure interval.
A negative χ that is projected away between checkpoints requires the proposed
stage observer; coarse checkpoint cadence cannot establish its first appearance.

## A3. Regrid, overlapping footprints and numerical transfers

Both attempts show L9 **144→128** at printed t=**134.176**, and L10
**256→240** at **134.436**. L8 remains 128 boxes and L11/L12 remain 256.
See [t24-box-count-transitions.csv](t24-box-count-transitions.csv) and
[t24-level10-count-history.csv](t24-level10-count-history.csv). The L10
drop occurred before the final regrid at **135.174**, not first at that
event. Complete LB records show all 224 ranks have 1–2 boxes in the
240-box layout. A box-count reduction does not establish a hole in the
physical union or an idle/absent MPI rank.

The two-centre criterion at `EMSBH2DLevel.cpp:430–465` is a finite
cellwise maximum, with forcing `r<1.2*r_extraction` at
`EMSExtractionTaggingCriterion.hpp:93–106`. Creating L10, L11, L12 uses
forcing radii **0.21875**, **0.109375**, **0.0546875 M_i**, respectively.
The L10 forced circles first overlap once d<0.4375; the retained tracks
cross that threshold between t=133 (d=0.451961) and t=133.875
(d=0.4149546). The L11/L12 nominal circles remain disjoint at the
failure's roughly 0.365 separation. Their actual rectangular footprints
also include tag buffer 3, grid/nesting buffer 8 and box quantization.

The last logged regrid changes descendants **7–12** of base level 6;
levels 0–6 are not redefined by that callback. The base-6 interval is 1;
zero descendant intervals do not freeze those descendants when their
ancestor regrids. Tracking runs on L6 after its children return
(`EMSBH2DLevel.cpp:486–488`), interpolating shift from minimum L11
(`Main_EMSBH2DBH.cpp:44–45`). There is no separation-dependent infinite
tag loop, division by separation or automatic transition to one label.

For the L11 parent cell containing the bad physical point, its centre is
(0.18243408203125,0.04742431640625). Its distance to the observed right
puncture at t=134.75 is **0.0481357**, less than the L12 forcing radius.
Linear extrapolation of the **last two retained coarse track rows**, explicitly
an inference and not a recovered tracking sample, gives x2≈0.18254262 at
the likely lattice regrid clock 135.173828125, distance **0.04742444**.
The same cell remains in the nominal L12 forced core under that extrapolation.
It was already covered at step 152. Thus the tracks/tag rule do not support
the specific claim that this cell loses L11/L12 coverage as the holes close.

**Actual gained/lost coverage at 135.174 cannot be reconstructed from the
collected evidence.** The retained tracker writes only coarse-clock rows,
not the intervening L6 centres/tag sets; the current fields at that regrid
are also missing. LB omits coordinates (`LoadBalance.cpp:443`), and the
supplied `checkpoint-layout.txt` is Lustre striping metadata. Stable box
counts do not imply stable unions or partitions. All **36** saved
step-152 nesting tests (0/4/8 parent-cell support) pass; they certify only
that checkpoint. No post-135.174 nesting or coverage claim is substituted
for the missing layout. The restart design writes the required unions and
actual tag-time centres.

`GRAMRLevel.cpp:387–450` defines the new grids, spatially interpolates
**all** fine data from the coarser **current** `m_state_new`, then copies
old same-level valid overlaps back. New fine cells receive tensor-product
six-point-per-direction Lagrange interpolation (degree five), not injection
or a positivity-preserving average (`PointAMRTransfer.hpp:152–173,192–224`).
The regrid fill uses **no temporal interpolation**. Ordinary RHS
coarse–fine ghost filling instead uses the point RK4 stage interpolator
(`GRAMRLevel.cpp:1024–1091`, `PointAMRTransfer.hpp:113–118`). During
synchronized base-level regrid the descendants normally share that base
clock; no clock mismatch at this event is established.

Coarse restriction also uses **six by six fine point values**, with signed
weights and an anchor (`PointAMRTransfer.hpp:123–144`), not injection or
volume averaging. A coarse cell exposed when its child footprint withdraws
retains its last restricted value; there is no special exposure floor or
new injection. If that coarse level itself gains a cell, its own coarse
interpolation supplies it. Old valid overlaps are copied without arithmetic,
while their ghost support can change. This is why the observer must record
the donor stencil, retained overlap and covered-parent state independently.

## A4. Detected NaN versus job termination

The dump is inside `NanCheck.hpp:51–68`'s conditional `omp single`;
that construct has an implicit barrier. It is called by the y-row parallel
loop at `BoxLoops.impl.hpp:82–94`, on valid cells after full RK projection.
The stop condition is NaN or |value|>1e20 (including Inf), not NaN alone.
After the single block, line 68 calls `MayDay::Error` **inside the worker
region**. In external Chombo `MayDay.cpp:67–78`, default `abort()`
precedes `MPI_Abort`. `SetupFunctions.hpp:52` uses `MPI_Init`, without
negotiating worker-thread MPI calls. The conditional work-sharing barrier
is invalid parallel control; it is a concrete failure-path defect even
if a given runtime happens to release it. A worker must not wait there
for finite-cell workers to execute that same conditional construct.

The rank-96 dump completes in both cluster attempts, but no MayDay error
line is present in the collected stdout/stderr. That is consistent with
not passing the single barrier, but does **not** prove that wait site:
MayDay's `fflush(NULL)` or runtime/signal machinery could also intervene.
Other ranks can be waiting for that rank in an evolved-field exchange,
coarse–fine copy or later tracker collective; their last regrid messages
are emitted at regrid return, not their current stacks. No specific
collective or stack is established from progress messages alone.

The controlled Tests-only probe uses the **actual pinned checker** and
two OpenMP threads. Its own native receipts, rather than the outer runner
code, give:

{probe_table}

The earlier timeout hypothesis is **not reproduced** locally: both isolated
legacy cases abort normally under arm64 GCC16/libgomp. The safe collector
returns after the parallel join, closes its record, and preserves all
**112 Float64 cell values** bit for bit, including the injected NaN bits;
the finite case also preserves 112. Probe/build recipes and evidence are
[t24-nan-barrier-probe.cpp](t24-nan-barrier-probe.cpp),
[t24-nan-barrier-qualification.json](t24-nan-barrier-qualification.json)
and `/private/tmp/ems-t24/nan-barrier/`. This is a local non-MPI fixture,
not a repair applied to production or a reproduction of the cluster hang.

The minimal correct **proposal** is a read-only stop-flag/witness collector
with no conditional barrier, followed by a worker join, a closed per-rank
abort record and a direct nonzero `MPI_Abort` on the caller/main thread
(serial exit in non-MPI builds). No failure-path collective is introduced.
Use a default-off `ems_safe_nan_abort` if preserving existing failure
behaviour is required. Preserve the existing detection threshold. The
record must identify rank, level, phase/RK stage, exact clock, coordinates,
field class and IEEE bits. The MPI termination qualification and exact
cluster wait diagnosis remain for the submit worker. Neither disabling
NaN checks nor changing an evolution floor is an acceptable remedy for
the missing termination.

## B. The diagnostic restart packet

[t24-diagnostic-restart-design.md](t24-diagnostic-restart-design.md) fixes
the production-layout replay from 152 to 155, checkpoint interval **2**
(to preserve 154), absent static paths, unchanged physical parameters,
the default-off recording contract, native phase insertion sites and
qualification. The region is midpoint-relative **x∈[−0.4,0.4],
y∈[0,0.15] M_i**, levels **8–12**, recording all 28 evolved fields and
separate valid/ghost classes. Regrid field window **[135.15,135.19]**
covers the final event and its neighbours; phase/nonfinite scans run from
the restart to **135.625**, on all levels and RK phases.

Existing native checkpoints suffice for step 154. Existing plots and the
post-layout observer do not provide the required asynchronous pre/post
fields: a read-only default-off diagnostic variant must be implemented
and qualified before the submit worker uses this design. The design names
every callback and every required datum, including actual centres,
old/proposed/new boxes, restriction/prolongation donors and the first
local phase that becomes nonfinite. No implementation/qualification of
those new production MPI field hooks is claimed in this turn.

Saved-layout measurement is **86.80 MB per full intersecting-FAB snapshot**
including ghosts; three field phases at three selected regrids are about
781 MB. Budget **1.2 GB** for snapshots/metadata and **1.6 GB** for the
step-154 checkpoint, **4 GB** incremental allocation, an **8 GB** disk gate
and **15-minute** native phase cap. The controller's approximately
two-minute reproduction is a baseline, not a measured observer cost.
Save per-rank witnesses/completion receipts even after nonzero abort.

## Candidate mechanisms, ranked without an evolution remedy

| Rank | Mechanism | Evidence for | Evidence against / still missing |
| --- | --- | --- | --- |
| 1 | Point transfer and subsequent floor/RHS amplification of an inner state | Signed restriction demonstrably creates negative coarse lapse from positive fine lapse; negative fine ghosts and widespread lapse-floor activity exist; failed χ/lapse are floored and metric/shift have exploded | χ is positive at 152; the specific bad cell and requested fields are then regular and correctly restricted. No first failing operation/stage or donor at 135.174 is retained. Floor may be an effect. |
| 2 | Inner gauge/dynamical instability or insufficient resolution near the moving holes, amplified in a covered coarse advance | Very small lapse, inner Θ magnitude≈0.69 and sizeable K are already present; L10 advances even its fine-covered cells and detects the failure there | No short-time spatial/RK growth history distinguishes gauge, truncation, KO or transfer effects. No local replay reaches the event. |
| 3 | Regrid-created cell/support or overlap-copy/communication defect | The final dump follows a regrid; moving regions alter support/partition and new fine interpolation is not positivity preserving | Bad cell is already covered at 152 and predicted to remain forced; L10 count drops at 134.436, not 135.174; saved nesting passes. Actual later unions, donors and overlap identity missing. |
| 4 | Tracker failure or a separation-dependent tagging nontermination | Moving centres determine all fine footprints; tracker has collective exchanges | Both centres are finite and distinct; two-centre max is finite even at coincidence. An explicit numerical NaN is recorded, so a pure no-progress tag loop does not explain the primary failure. |
| 5 | Pure load imbalance / missing rank / disappearance of the finest level | Repartition and native collective waits could delay progress | Repeated same-cell NaN, all 224 ranks in complete 240-box LB records, and L11/L12 still 256 strongly weigh against a purely performance/layout-count explanation. |

The **abort-control defect** is separate from these numerical mechanisms.
No expensive global refinement, gauge or transfer change is admitted from
this evidence. The proposed single diagnostic restart distinguishes whether
the initiating event is restriction, regrid interpolation, a changed ghost,
projection or an old valid-cell RHS/update. It also supplies a clean failure
record instead of interpreting an external phase-cap stop as the physics.

## Completed replay and resource receipts

The local native receipt is `/private/tmp/ems-t24/restart152-rss/native.receipt.json`:
return **−15**, `{native['gate_reason']}`, **{native['elapsed_s']:.2f} s**.
It ends at the first coarse endpoint t=133.875 and does not reach the
failure. Sampled RSS peaks at **5,050,679,296 bytes**, native wait4 RSS
at **5,119,148,032 bytes**, footprint **7,948,230,600 bytes**. Its own
done marker is 125. No backtrace of the blow-up was obtained. No full replay
is retried. This is a resource stop, not numerical success or failure.

The prior small default-off/on post-layout observer regression remains
qualified: **2,588,800 defined Float64 comparisons**, zero bit mismatches,
through actual regrid. Its scope does not qualify the proposed field/abort
hooks. The old sparse probe covers only saved-box grow/subtract/iteration,
not new tag unions or the later layout. All receipts, including failed
preparatory jobs, remain preserved; [t24-resources.csv](t24-resources.csv)
and [t24-manifest.txt](t24-manifest.txt) seal completed local evidence.
No commit, SSH, cluster execution, positive-time static read or evolution
fix was made.
'''
    (HERE/'t24-merger-stall.md').write_text(report)
    resources = []
    for p in sorted(ROOT.rglob('*.resources.json')):
        q = json.loads(p.read_text())
        resources.append(dict(receipt=str(p),process=q['process'],returncode=q['returncode'],
            peak_RSS_bytes=q['peak_rss_bytes'],wall_seconds=q['wall_seconds'],gate_reason=q['gate_reason']))
    resources += [dict(receipt=str(ROOT/'restart152-rss/native.receipt.json'),process='native-replay-wait4',
        returncode=native['native_returncode'],peak_RSS_bytes=native['native_wait4_peak_RSS_bytes'],
        wall_seconds=native['elapsed_s'],gate_reason=native['gate_reason']),
        dict(receipt=str(HERE/'t24-offline-receipt.json'),process='checkpoint-reader',returncode=0,
        peak_RSS_bytes=audit['peak_RSS_bytes'],wall_seconds=audit['elapsed_s'],gate_reason='')]
    build = json.loads((ROOT/'nan-probe-build.json').read_text())
    resources.append(dict(receipt=str(ROOT/'nan-probe-build.json'),process='nan-probe-compile',
        returncode=build['returncode'],peak_RSS_bytes=build['peak_RSS_bytes'],wall_seconds='',gate_reason=''))
    with (HERE/'t24-resources.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(resources[0]));w.writeheader();w.writerows(resources)
    status = dict(status='READY-EXCEPT',numerical_blowup='ESTABLISHED_BY_TWO_PRINTED_WITNESSES',
        initiating_operation='UNESTABLISHED',cluster_abort_wait='UNESTABLISHED',
        mislabelled_min_chi='DIAGNOSTIC_MOD_F',negative_valid_chi_at152=False,
        negative_valid_lapse='1002, all covered',safe_abort_proposal='NOT_APPLIED_TO_PRODUCTION',
        new_MPI_field_hooks='DESIGN_ONLY_NOT_YET_IMPLEMENTED_OR_QUALIFIED',
        full_local_replay='RESOURCE_STOP_NO_RETRY',native_receipt=native,
        peak_native_RSS_bytes=native['native_wait4_peak_RSS_bytes'],offline_peak_RSS_bytes=audit['peak_RSS_bytes'],
        input_manifests=audit['input_manifests'],local_abort_probe=probe,
        running_jobs=[],evolution_fixes_applied=False)
    (HERE/'t24-status.json').write_text(json.dumps(status,indent=2)+'\n')
    readme=HERE/'README.md'
    text=readme.read_text()
    begin=text.index('## T24 — binary merger stall, d507438')
    end=text.find('\n## ',begin+1)
    if end<0:end=len(text)
    replacement='''## T24 — binary merger stall, d507438

**READY-EXCEPT: numerical blow-up established; first failing operation and cluster
abort wait unestablished.** [Report](t24-merger-stall.md),
[offline reader](t24-checkpoint-read.py), [native/region values](t24-checkpoint-extrema.csv),
[coverage](t24-cell-coverage.csv), [point-restriction witnesses](t24-point-restriction.csv),
[abort probe](t24-nan-barrier-qualification.json),
[diagnostic restart design](t24-diagnostic-restart-design.md), [status](t24-status.json),
[resources](t24-resources.csv), [manifest](t24-manifest.txt).
Both manifests (1,855 + 692 files) and the 1.5 GB step-152 checkpoint verify.
The negative `min_chi.dat` history is **min(mod_F)** due to a diagnostic/evolution
component-index mix-up. Saved valid chi is positive everywhere (minimum
3.97613e-6); 1,002 negative valid lapse cells are all covered coarse cells.
Point restriction reproduces an L11 negative lapse from 36 positive L12 values
bit for bit. The bad L10 cell is already covered by L11/L12 at step 152.
L9 144→128 and L10 256→240 occur at logged times 134.176 and 134.436,
before the final regrid. Both attempts dump the same nonfinite state on rank 96.
The conditional OpenMP single/MayDay failure path is unsafe, but a local actual
checker probe aborts; the exact cluster wait is not proved. The proposed joined
witness/direct MPI abort and pre/post-regrid/stage hooks are **design only**.
The full local replay ended at t=133.875 under the 7.5 GB footprint gate,
native −15, peak wait4 RSS 5.119 GB; it is not retried. Offline analysis peak
RSS 0.141 GB. No job remains running, no evolution fix is applied, and no
commit, SSH or cluster operation is performed.
'''
    readme.write_text(text[:begin]+replacement+text[end:])
    paths=[p for p in HERE.glob('t24-*') if p.is_file() and p.name!='t24-manifest.txt']
    paths += [readme,HERE/'t13-run.py',HERE/'t23-compare.py',ROOT/'build/build-spec.json',
              ROOT/'restart152-rss/native.receipt.json',ROOT/'restart152-rss/done.exit',
              ROOT/'nan-probe-build.json']
    paths += sorted(ROOT.glob('nan-barrier/*-measured/*.resources.json'))
    paths += sorted(ROOT.glob('nan-barrier/*.record.json'))
    paths += [ROOT/'source'/p for p in ('Examples/EMS/EMSBH2DLevel.cpp',
              'Examples/EMS/DiagnosticVariables.hpp','Examples/EMS/UserVariables.hpp',
              'Source/GRChomboCore/GRAMRLevel.cpp','Source/GRChomboCore/PointAMRTransfer.hpp',
              'Source/BoxUtils/AMRReductions.impl.hpp','Source/BoxUtils/NanCheck.hpp',
              'Source/BoxUtils/BoxLoops.impl.hpp','Source/CCZ4/PositiveChiAndAlpha.hpp',
              'Source/Cartoon/EMSCartoonLorentzScalars.hpp')]
    paths += [ROOT/'Chombo/lib/src/BaseTools/MayDay.cpp',ROOT/'Chombo/lib/src/AMRTimeDependent/AMR.cpp']
    entries=[]
    for p in sorted(set(paths)):
        with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        entries.append(f'{digest}  {p.stat().st_size}  {p}\n')
    (HERE/'t24-manifest.txt').write_text('# T24 completed read-only evidence; no commit; no evolution fix\n'
        '# Numerical blow-up established; initiating operation/cluster wait UNESTABLISHED\n'
        '# New production MPI diagnostic/abort hooks are DESIGN ONLY\n'
        '# Native peak RSS 5119148032; offline peak RSS '+str(audit['peak_RSS_bytes'])+' bytes\n'
        '# SHA256 bytes absolute_path\n'+''.join(entries))
    print('T24_PACKET_SEALED',len(entries),'files; no pending job')


if __name__=='__main__':
    main()
