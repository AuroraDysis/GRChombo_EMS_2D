# T14 Stage C — registered before evolution

Baseline: fork main 949a744, unchanged T13 production numerical objects. Consult 7
Stage C and the prepared t13-stageC-reading.md govern the decision. No production
C++ or numerical-method change; sigma=1, point transfers, initial-data guard on.
The new roots and markers do not overwrite any T13 result.

Pre-evolution initialization amendment: the real max14 census with the original
factory crashed in TreeIntVectSet::clearTree/grow during level-13 tagging.
Chombo's traversal stacks have 24 entries; these global indices require a deeper
tree. Peak RSS was 0.730 GB, not an OOM. The Tests-only T14Launch.cpp includes the
production main verbatim and links the existing production numerical objects.
One option, t14_dense_initial_tags (default false), selects a Tests factory that
calls the ORIGINAL tagging routine below level 13 and uses the library's dense
set representation at level 13. Native criteria, selected cells, dilation and
clipping are preserved; production source and Chombo remain untouched. The test
compiler flag -fno-access-control permits this initialization-only override
without editing production access specifiers. Real old-box and t=0 retained
state identity, plus default/off and enabled two-step E/reference controls in
legacy/point modes, must pass before this Tests launcher is used for evolution.

First require a real initialization-only census with the production AMR factory
and level implementation, proving exact old boxes on levels 0–12, the two new
faces, eight-cell block alignment, positive floor margins and complete W stencil
support. The Tests-only census executable links the existing production objects.
It never advances. Initial positivity/floor or old-face failure stops this card.

Run alpha_K at max13 and max14 with dt0=7/48 M, and max14 with dt0=7/96 M.
Each stops at its last native step inside 0.002 M. Expected endpoints: nominal
112/224 steps at 0.001993815104166667 M; the half-step control has 449 steps and
ends at 0.001998265584309896 M. Use its first 448 steps for the common endpoint.
Common history times are the T13 max12 nominal times (every 2 max13, 4 max14,
8 max14-control steps). Verify actual timestamps to floating-point tolerance;
no temporal interpolation is expected. Report any measured interpolation error
if this expectation fails. No pulse alignment or background fitting.

On W=[0.00075,0.0025] M, axis and diagonal, use the SAME fixed 258 physical sample
points and anchor-subtracted tensor P6/P8 current-state sampling as T13. Compute
native metric Gamma with T13Replay's fourth-order cartoon derivatives, including
ww/hww. Compare raw signed longitudinal Gamma, metric Gamma, shift and lapse.
Report both peak and unweighted radial RMS norms of each resolution difference,
rho_D at the common endpoint and at every common nonzero time, plus the history
supremum of each difference norm. Initial Gamma zero/floor differences are
reported as uncertainty dominated, never divided into a spurious order.

For each difference, use the sum of the two profiles' P6/P8 difference norms as
the conservative interpolation bound. Require both spatial differences >5 times
their respective bounds. The finest dt/2 profile difference is the measured
temporal error; call it subdominant only if it is <=0.2 of BOTH spatial difference
norms (a conservative quantified reading of the expert's qualitative criterion).
Carry the T13 baseline's measured temporal error separately. Report all failed
qualifications instead of suppressing them. Do not infer fourth order from
rho_D<=0.8, whose equivalent dyadic order is only about 0.322.

Separately define Gamma disturbance as each rung's change from its own numerical
t=0 field on the same native cells. Report its P6/P8 peak/RMS, native-ray
corroboration (axis first row y=h/2, diagonal native centers), C_Gamma, Ham/Mom,
probe history at 0.0015 M and floors/nonfinites. A dominant disturbance lobe uses
the largest absolute native value, contiguous same-sign cells and the width
above half that value, with linearly interpolated half-height crossings. Flag
lobes clipped by W and distinguish native cell width from physical ray spacing.
There is NO nonzero-amplitude-stability requirement.

Replay every recorded finest-level stage for the first four steps using the
unchanged native T13Replay harness. Require actual total RHS and projected inputs
bit-identical, RK additions within the recorded raw-operation roundoff budget,
and the Gamma term sum within the native raw-stencil budget. Stream stages,
never hold the full max14 recorder in memory. Report geometric/advection/direct
KO, trace-projection and floor updates separately, at puncture cells and W.
If feedback attribution requires it, reuse T13's offline first-predictor
dependency test with and without A-only/all KO; this is a diagnostic replay,
not an altered evolution or a claimed additive outgoing-amplitude fraction.

Physical characteristic bounds use the frozen longitudinal shift family
|beta|+sqrt(lambda_max(inverse conformal metric)), compared with the current
light/lapse speeds. Use a conservative 1.1 upper bound and a 7h combined sampling
plus RHS coordinate support buffer, verify the bound on evolved recorder states.
These are characteristic separation estimates; FD/KO numerical tails do not have
strictly compact support. W reads no coarse–fine ghost directly.

Decision: significant rho_D<=0.8 with a resolved stabilizing pulse or disturbance
decreasing toward zero advances to a DRAFT Stage D; 0.8<rho_D<1, floors or unresolved
uncertainty gives inconclusive, with no automatic extra rung. Significant rho_D>=1
means lapse alone is insufficient at these resolutions; identify KO/projection
regeneration or gauge readjustment from the stages. Report mixed field/ray/norm
results explicitly; a mixed screen does not admit production. This 0.002 M window
cannot exclude delayed emission or establish exterior/100 M behavior.

Only after a pass draft fresh alpha_K E three-grid exterior chains through
2.625/4.375/10.5 M, exp-0020 common transport faces with the justified puncture
levels, full constraint battery and a predeclared puncture-excluded domain.
Design order is four on a smooth domain; assume no regularity exemption. Keep
numerical-t=0 horizon budgets |Delta A/A|<=1e-3, |Delta Q/Q|<=2e-4 and independent
sphere charges. Estimate cluster cost from the recorded exp-0020 rates. Do not
launch Stage D under this task.

Serial queue wall estimate: 536 s, upper 803 s. Peak per-process estimates
2.15/2.33/2.33 GB, measured 3 GB process AND active-tree gates; OMP=2, BLAS=1,
two single-thread compression streams, <=4 compute threads. New output root
/private/tmp/ems-t14, hard 5.7 GB gate. Every numerical process invokes
/usr/bin/time -l, with T13's explicit child/wait4 RSS fallback and libproc live
RSS/footprint monitoring where sandbox kern.clockrate is denied. Record every
measured peak and stop on the cap. Any long queue is detached, then end the turn.
