#!/usr/bin/env python3
"""Seal the T24 exp0026/CFL/abort/dt follow-up from its own completed receipts."""
import csv
import difflib
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
ROOT=Path('/private/tmp/ems-t24')


def read(name):return json.loads((HERE/name).read_text())
def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,row))+' |\n' for row in rows)


def revised_reading(report,g,nn):
    rule=read('t24-stability-rule.json');cluster=read('t24-cluster-check.json');hs=read('t24-horizon-summary.json')
    assert cluster['status']=='VERIFIED' and cluster['manifest_files']==147
    assert rule['peak_RSS_bytes']<6e9 and len(rule['rows'])==48
    assert all(q['native_exit']==1 and q['peak_RSS_bytes']<6e9 for q in hs['receipts'])
    begin=report.index('**READY-EXCEPT.**');end=report.index('## 1.',begin)
    report=report[:begin]+'''**READY-EXCEPT.** The collected exp-0026 interventions and the code-stencil
audit support the shift/Gamma discrete-stability candidate under the registered
reading. They do not uniquely prove its mechanism or admit a production restart.
The expanded RK4 rule, all local bounded numerical horizon probes and the
serial-qualified opt-in safe abort are complete. No strict individual or common
horizon is demonstrated. Original common-finder iteration logs are absent from
the collected packet; MPI abort execution is blocked at local MPI_Init.
No full evolution replay is repeated. No job remains pending. No gauge, KO,
transfer, initial-data or evolution-equation remedy is applied; no commit,
SSH or cluster operation is performed.

'''+report[end:]
    begin=report.index('**Established:**');end=report.index('## 2.',begin)
    limits=table(['dt/dx','sigma','weighted S limit','metric U limit, c=.75','U limit, c=.375','U limit, c=1','sampled S, measured cell/c=.75'],[
        [r['dt_dx'],r['sigma'],f"{r['weighted_S_Nyquist_limit']:.8f}",f"{r['weighted_S_Nyquist_limit']:.8f}",
         f"{2*r['weighted_S_Nyquist_limit']:.8f}",f"{.75*r['weighted_S_Nyquist_limit']:.8f}",
         f"{r['measured_shape_shift_sampled_S_lo']:.8f}"] for r in rule['rows'] if r['cGamma']==.75])
    measured=table(['Leg/step','t','max weighted S','cells >4.9','native evidence'],[
        [r['leg']+'/'+str(r['step']),r['time_M'],f"{float(r['max_S']):.9f}",r['cells_S_gt_4p9'],r['native_rank_exits']] for r in cluster['summary']])
    report=report[:begin]+f'''### RK4 rule, expanded before use in a run design

For arbitrary positive c_Gamma define **U=H22+0.75 H11** and the weighted
rho-mode coefficient **S=c_Gamma*(H11+4 H22/3)=(c_Gamma/0.75) U**.
Earlier columns named S in the single-cell script mean U; the CSV now also
records explicit `metric_U` and `weighted_S` columns. With diagonal H and
H22 >= H11, zero shift, alpha=0 and negligible eta*dx, the two-direction
Nyquist branch has z=−2 sigma nu ± i nu sqrt(16 S/3), nu=dt/dx.
Let y_RK(a) be the first positive edge of the RK4 stable segment on Re(z)=a.
Then **S_lim=3 y_RK(−2 sigma nu)^2/(16 nu^2)**. The weighted limit is
independent of c_Gamma; the permitted metric U is 0.75 S_lim/c_Gamma.
For general off-diagonal H use S_eff=c_Gamma*(trace(H)+lambda_max(H)/3)
at Nyquist, and scan all wave numbers with the actual shift and metric.
The restriction on H's ordering is part of the rho-mode claim, not an implicit
assumption about every cell in a simulation.

{limits}

[t24-stability-rule.py](t24-stability-rule.py) supplies all **48** requested
(nu,sigma,c_Gamma) combinations in [t24-stability-limits.csv](t24-stability-limits.csv).
For each real RK argument it closes the exact conjugate-polynomial witness
and independently checks the positive root with a 70-digit four-stage recurrence.
The last column is a 257×257 **decaying principal-branch** scan at the measured
L12 metric shape, shift, Gamma and dx, scaling the inverse metric to each S.
It is a reproducible sampled threshold, not a universal table depending only
on three parameters. `t24-stability-rule.json` records the narrow numerical
brackets, modes, witness statuses and resource receipt. Shift, h12, radial
source terms and background gradients prevent a complete universal criterion
in just (nu,sigma,c_Gamma). A safety factor is not invented after these data.
Sigma variants are offline calculations only; evolution retains sigma=1.

At nu=0.25 the measured-cell h/chi advection and B-decay branches have max
RK amplification <=1 across the full sampled domain. The full 11-row cartoon
scan at weighted S=25 has no decaying-mode RK violation; at S=26 its beta2/
Gamma2 wave is unstable (max |R|=1.10025). No other **decaying** alpha-independent
mode crosses the RK boundary first in this audit. The positive-real,
low-frequency homogeneous cylindrical modes reported above remain present
even at smaller dt: they are source/transport growth in this frozen model,
not a competing numerical CFL threshold. Their interpretation in the actual
inhomogeneous radial problem is **not closed**. K/Theta/A and lapse perturbations
do not add a fast feedback loop at fixed alpha=0: their advected/zero-order
blocks are triangular with respect to this principal subsystem. In the actually
compiled CCZ4 file, `COVARIANTZ4` is defined and the damping uses
kappa1*alpha/(0.005+alpha), which also vanishes as alpha approaches zero
(`CCZ4Cartoon.impl.hpp:13,346–351`). No finite-alpha, matter, floor-projection,
AMR-interface or boundary stability guarantee follows from this limit.

### Growth and the exp-0026 interventions

At nu=.5, sigma=1, **S=4.96 / 5.03** gives scalar Nyquist gains
**1.015262958 / 1.035814535 per fine step**, e-folds **66.01684 / 28.41869
steps**. On L10/L11/L12 these are e-fold times
0.056411/0.028205/0.014103 and 0.024284/0.012142/0.006071 M_i.
The measured-shift full cartoon variants instead give Nyquist gains
**1.052098413 / 1.073316283** (whole-grid maxima 1.052878456/1.074156157).
They are a sensitivity scan with the step-152 background held fixed, not
measured growth of a trajectory. Both sets are exposed rather than using
the scalar estimate as the full code's amplification.

The collected S153/S154 slope gives the scalar crossing **t={rule['linear_ramp']['crossing_M']:.9f}**
and extrapolated S={rule['linear_ramp']['end_S']:.9f} at 135.173828125.
The resulting quasi-frozen Nyquist ramp predicts log10 amplitude gains
**6.903/13.806/27.611** on L10/L11/L12 over the **0.765329 M_i** interval.
This accounts for how a tiny grid-scale seed can become very large after
hundreds to thousands of fine steps; it is consistent with a delayed NaN.
It is not a fit to measured checkerboard amplitudes: the cell of the spatial
maximum can migrate, the coefficients change, and the scalar model omits
shift and AMR effects. See [t24-stability-growth.csv](t24-stability-growth.csv).

{measured}

[t24-cluster-evidence.py](t24-cluster-evidence.py) checked **147** locally
collected files against their SHA256 manifests. It then checked each G/D/D0
**native** receipt: 224/224 rank exits zero, the independent finite-state
reduction has no unreadable checkpoints or valid non-finite values, and the
dt readback and absent-static-input proofs hold. The own reduction receipt
also reports all active metric inversions valid. The summary does not mistake
the scheduler exit for a numerical pass. Raw intervention checkpoints remain
on the cluster; no bulk file was copied into this worktree.

Both one-parameter interventions pass coarse 155 and the original event.
G lowers c_Gamma to .375, giving weighted S=2.658946 at t=138.25;
D keeps c_Gamma=.75 and halves dt, giving S=5.244598 at t=136.5 against the
scalar limit 25.448682. Restart already supplies the requested per-level dt.
The checkpoint-state gauge kick in G leaves B untouched, so G and D are
different gauge trajectories and do not certify each other's physical accuracy.
**Supported:** the collapsed-lapse beta/Gamma loop approaches/crosses its
discrete RK boundary and two targeted interventions prevent the original NaN.
**Not established:** a unique causal mechanism, the exact first unstable
stage, measured checkerboard growth, convergence or production admission.
The appropriate run-design rule evaluates the complete discrete symbol with
the current numerical metric/shift and leaves an explicit margin; S<4.906
alone is too permissive at the measured shift. The existing half-dt control
has a large measured margin through its stated bound, not through merger.

'''+report[end:]
    begin=report.index('There is therefore **no qualified apparent horizon');end=report.index('## 3.',begin)
    commonrows=list(csv.DictReader((HERE/'t24-horizon-all-probes.csv').open()))
    common=table(['Probe','Seed','Squared expansion','A (trial)','Q (trial)','rho_max (trial)','stop'],[
        [r['probe'],r['seed_radius'],f"{float(r['expansion_squared']):.8g}",f"{float(r['A']):.8g}",
         f"{float(r['Q']):.8g}",f"{float(r['rho_max']):.8g}",r['status']] for r in commonrows if r['hole']=='common'])
    receipts=table(['Probe','Native exit','Wall s','RSS GB','Best final squared residual'],[
        [r['probe'],r['native_exit'],f"{r['wall_s']:.3f}",f"{r['peak_RSS_bytes']/1e9:.3f}",f"{r['best_final_squared_residual']:.8g}"] for r in hs['receipts']])
    report=report[:begin]+f'''### Expanded common search and the production flow

The same unchanged author's finder searched fresh spheres at numerical midpoint
x=2240 with coordinate radii **.3,.5,.75,1,1.5,2,3,4 M_i**. N48/chase1/quota1
was bounded by 128 updates; N96/chase1/quota1 by 90 s and 2048 updates.
An additional N48/chase.125/quota1 control tried .75 and 3 for 90 s.
All queries used the saved t=133 state and native interpolation, never a static
seed, target or subtraction. All retain the same first-stage squared threshold
1e-7. The N96 run ended after 587 updates; the slow N48 run after 574.

{common}

{receipts}

Every native probe exits **1**, with no resource gate: these are bounded
finder failures, not numerical evolution failures. None qualifies even stage
one, so later 1e-10/1e-12 qualification and N48/N96 A/Q agreement cannot be
claimed. **No qualified individual or common horizon, area or charge is
established at t=133.** Trial areas/charges/extents in the table describe
unconverged surfaces. Thus whether rho=.047–.08 above either puncture is
inside an apparent horizon is **undetermined**. Failure of these finite
searches does not prove the absence of a trapped surface or of a non-star-shaped
surface. All trials are reproduced by [t24-horizon-search.py](t24-horizon-search.py),
[t24-horizon-all-probes.csv](t24-horizon-all-probes.csv) and their own receipts
in [t24-horizon-summary.json](t24-horizon-summary.json).

The collected run summary, [t24-common-horizon-history.csv](t24-common-horizon-history.csv),
has **18** common rows from steps 96–128, all TIME_CAP/UNRESOLVED: areas
1016.6–1279.9, squared residuals .0100433–563.326, and about 234.8 s/call.
The cap is the wrapper's `offline_seconds`; the author's stopping test is
area-average Theta_plus^2 (`RHSurf.hpp:228`) against the specified threshold,
not stabilization of area. Its update step is
delta_f=−.25*chase*dtheta^2*f^2*(Theta_plus+.01*Theta_plus/sqrt(Theta_plus^2+1e-13)*sqrt(abs(Theta_plus)))
(`RHSurf.hpp:823–839`). FAR uses one fresh plus **five stale** steps before
another interpolation (`RHUnion.hpp:297–313`). Per-node radii are clamped
to [.0001,10]; a mean radius outside those bounds resets to 5
(`RHUnion.hpp:220–229`). Newton polishing is commented out, so changing its
parameter cannot rescue this build. The residual/time tests, step and bounds
are the unchanged author's machinery, not a new horizon solver.

The actual driver initially seeds common surfaces at d/2+.5 and d+1 and
midpoint. At later calls it copies the preceding `common.dat` into **both**
seeds (`ev1/.../runtime/production.py:231–236`). It saves the lowest-error
finite trial as `common.dat` even when UNRESOLVED (`:259–265`); thus subsequent
paired searches are not independent fresh-radius trials. Already at t=133,
fresh radius .3 has numerical area **370.837**, and fresh radii 1–4 have
areas about 777–828. A large area is neither an admission nor a measure of
proximity to a horizon on this numerical geometry.

The controlled flow-step comparison is informative: for seed 3 at accumulated
chase 72, N48/chase1 has squared residual **2.37627**; at accumulated chase
71.75, N48/chase.125 gives **.0508586** (about 46.7 times smaller).
Thus large angular-flow excursions are sensitive to the step, beyond merely
allowing less chase time. Yet the smaller step remains more than five orders
above the strict threshold and does not demonstrate a horizon. The numerical
pilot with larger seeds develops very distorted and near-bound trial radii.
This supports investigating bounded, independently reseeded flow controls;
it does not prove the area-1000 history was caused by a reset or a particular
seed. The collected **ev1 packet lacks `diagnostics/` and `shapes/`**, including
the original finder run/progress/rh_f logs. Its available rank logs are the
evolution logs (`RH_activate=false`). Therefore the exact original reset,
stale-step excursion or seed basin cannot be read from the supplied finder
log: that requested evidence is missing. The report separates these measured
flow and policy facts from an unproved diagnosis of each historical call.

'''+report[end:]
    report=report.replace('cell, level, time, dx, dt, phase, threshold, exit code and all evolved IEEE\nhexadecimal values,','cell, level, time, dx, dt, phase, threshold, exit code, finite values in\nexact hexadecimal notation and non-finite labels (not NaN payload bits),')
    begin=report.index('## 4.');report=report[:begin]+'''## Cancelled builds and retained evidence

The controller cancelled the restart-dt feature and the annular stage-recorder
build after the interventions. **Neither was started.** Restart already
re-derives every level's dt from dt_multiplier; the completed local dyadic
readback regression is retained in [t24-restart-dt.csv](t24-restart-dt.csv)
and agrees with the cluster's 13-level readback. The annular recorder remains
a historical **design only**, explicitly cancelled as a build in
[t24-annular-recorder-design.md](t24-annular-recorder-design.md); no hook code
exists. Its old size estimate remains evidence, not a submission plan.

## Mechanism ranking after both interventions

1. **Collapsed-lapse shift/Gamma discrete CFL crossing:** supported by the
   code's alpha-independent loop, beta2/Gamma2 unstable mode, active-S history
   and two successful targeted interventions. Exact stage/checkerboard growth
   and uniqueness remain unproved.
2. **Inherited point-transfer/resolution damage:** remains possible as a
   contributor to the large inner metric and S; a switch only at step 152 does
   not exonerate the preceding 152 steps. Negative covered lapse is excluded
   as the proximate requirement by legacy's failure with positive valid lapse.
3. **Regrid/interpolation amplification:** remains a possible seed or contributor;
   frequent descendant regrids coincide with every fine-time neighbourhood.
   The named final regrid does not change count, and a lower-CFL run with the
   same transfer/regrid policy passes it. A unique footprint-change trigger
   is unsupported.
4. **Tracker/tag nontermination:** disfavoured as the initiating cause by the
   reproduced numerical non-finite witnesses and successful interventions.
   The unsafe abort is a separate failure to terminate the MPI job cleanly.

The gauge control does not establish accuracy, convergence, horizon retention
or a unique cause. No evolution fix beyond the optional termination mechanism
is applied in this turn.
'''
    return report


def main():
    g=read('t24-gauge-stability.json');h=read('t24-horizon-search.json');c=read('t24-controls-qualification.json');mpi=read('t24-mpi-abort-qualification.json')
    assert g['scan_grid']==257 and g['peak_RSS_bytes']<6e9
    assert h['resources']['returncode']==h['native']['returncode']==1 and not h['resources']['gate_reason'] and len(h['rows'])==28 and h['advances']==0
    assert all(r['bit_mismatches']==0 for r in c['comparisons'])
    for r in c['resources']:
        q=json.loads(Path(r['receipt']).read_text());assert q['returncode']==q['child_measurement']['returncode']==r['native_returncode'] and not q['gate_reason'] and q['peak_rss_bytes']<6e9
    assert (ROOT/'controls/done.exit').read_text().strip()=='0'
    for r in mpi['cases']:
        q=json.loads(Path(r['receipt']).read_text())
        assert q['returncode']==q['child_measurement']['returncode']==r['native_returncode'] and not q['gate_reason'] and q['peak_rss_bytes']<6e9
    size=list(csv.DictReader((HERE/'t24-annulus-size.csv').open()))
    payload=sum(int(r['payload_bytes_per_coarse_step']) for r in size);window=sum(int(r['payload_bytes_per_window']) for r in size)
    nn=sum(r['defined_values'] for r in c['comparisons']);rr=[r for r in g['results'] if r['operator']=='full_cartoon' and r['case']!='no_shift_scalar']
    waves=table(['Case','metric U','weighted S','max |R|, full','max |R|, Re(lambda) <= 0','limiting decaying k/pi','Nyquist wave |R|'],[
        [r['case'],f"{r['metric_U']:.8f}",f"{r['weighted_S']:.8f}",f"{r['max_amplification']:.9f}",f"{r['max_decaying_amplification']:.9f}",
         f"({r['kx_over_pi']:.7g}, {r['ky_over_pi']:.7g})",f"{r['nyquist_wave_amplification']:.9f}"] for r in rr])
    hr=[r for r in h['rows'] if r['stage']=='0']
    horizons=table(['Hole','Seed radius','Squared expansion after 64 updates','Trial rho_max','Status'],[
        [r['hole'],r['seed_radius'],f"{float(r['expansion_squared']):.8g}",f"{r['rho_max']:.8g}",r['status']] for r in hr])
    s=g['sampled_S_transition_interval'];p=g['parameters']
    report=f'''# T24 follow-up: gauge stability, numerical horizons and safe abort

**READY-EXCEPT.** The frozen operator and native serial controls are complete;
the numerical blow-up's initiating mechanism remains unproved. The safe abort
is implemented, default off, with native serial regression. MPI qualification
is blocked before the checker runs by the local MPI shared-memory bootstrap.
The annular recorder is explicitly **designed, not built**, using the card's
permitted fallback. No evolution remedy, commit, SSH or cluster operation is
applied. No long replay is repeated and no job remains pending.

## Corrections and the completed transfer A/B

Both exp-0026 legs complete steps 153/154 and encounter a nonfinite state
inside coarse step 155. Point reports the original covered-L10 witness;
legacy reports its mirror on L12 near the other puncture. Legacy has no negative
valid lapse from step 153 onward and still fails. Therefore the switch at 152
does not rescue the inherited state, and negative covered lapse is not the
proximate explanation. This does not exonerate the previous 152 point-transfer
steps. The controller's collected report is
`stage-exp-0026/state/threads/ems-spectral-solver/runs/exp-0026/submissions/exp-0026/transfer-diagnostic-report.md`;
its two native stops, complete rank-log census and corrected diagnostic receipt
are distinguished from the scheduler/wrapper result.

Levels 7–12 are regridded **every L6 step, 0.013671875 M_i**. The final printed
regrid at 135.174 (lattice time 135.173828125) changed no box count. L9's
144→128 change is earlier, at 134.176; L10's 256→240 is at 134.436. A final
regrid is not a unique trigger, and constant counts alone do not fix coordinates.

At the failure the punctures are approximately x=±0.1826, given the measured
0.0192/M_i motion; the printed x=+0.1820 and −0.1826 cells at rho=0.047–0.050
are directly above them. The old saved-cell distance 0.06436 refers to the
**t=133** puncture position ±0.2259805; it is not their distance at failure.
The checkpoint values below remain numerical observations at their actual clock.

## 1. Code-level frozen-coefficient audit

All equation references here are to archived **d507438**, not the later gauge
control. `ExperimentalGauge.hpp:77–90` supplies beta_t=advec(beta)+c_Gamma Gamma
−eta beta−B and B_t=−0.1 B. Its lapse coefficient 1.8 and B decay 0.1 are
existing literals; this tranche changes neither. `CCZ4Cartoon.impl.hpp:369–408`
supplies the pure inverse-metric shift Laplacian and **1/3** longitudinal second
derivative, the cartoon reg03/reg04 contributions and the evolved Gamma row.
With the production CCZ4 formulation and kappa3=1, contracted Christoffel plus
2 Z/chi is **evolved Gamma** exactly. No contracted-Christoffel replacement is
silently made. Matter, curvature, A and damping contributions proportional to
lapse vanish in the specified alpha→0 limit. This does not remove the finite
alpha/chi terms from the actual evolution, especially once chi hits a floor.

Write H=inverse([[h11,h12],[h12,h22]]), Hww=1/hww, rho the cartoon radius,
d=partial_i beta^i+beta_rho/rho, and dw=partial_i beta^i−2 beta_rho/rho.
    The alpha=0 subsystem implemented by the audit is

```
beta_t = adv(beta) + c_Gamma Gamma - eta beta - B + KO(beta)
B_t = -0.1 B + KO(B)
Gamma_i,t = adv(Gamma_i) + (2/3) d Gamma_i - Gamma_j partial_j beta_i
            + H_jk partial_j partial_k beta_i
            + (1/3) H_ij partial_j partial_k beta_k
            + Hww[(partial_rho beta_i)/rho - delta_i,rho beta_rho/rho^2]
            + (1/3) H_ij[(partial_j beta_rho)/rho - delta_j,rho beta_rho/rho^2]
            + KO(Gamma_i)
h_ij,t = adv(h_ij) -(2/3) h_ij d
          + h_ki partial_j beta_k + h_kj partial_i beta_k + KO(h_ij)
hww_t = adv(hww) -(2/3) hww dw + KO(hww)
chi_t = adv(chi) -(2/3) chi d + KO(chi)
```

The script differentiates these expressions analytically around a homogeneous
Cartesian background, retains cylindrical source terms and delta(H)=−H delta(h) H,
and assembles an **11×11** complex matrix. Gamma and B perturbations are scaled
by dx, chi by its numerical background. The damping rates are physical rates:
eta enters as eta*dx, not eta=1 in grid units. Actual numerical Gamma and shift
are supplied. K/Theta/A have advective alpha-independent blocks; at alpha=0
they do not close an extra fast feedback loop into this subsystem. A lapse
perturbation can feed those rows, but its coupling back through alpha K vanishes
at alpha=0; this audit freezes the collapsed lapse rather than claiming a
complete finite-lapse CCZ4 spectrum.

`FourthOrderDerivatives.hpp:114–190,258–365` gives the pure second stencil
(-1,16,-30,16,-1)/(12 dx^2); mixed derivatives are the tensor product of
centered fourth-order first derivatives. The shift-dependent one-sided
advection stencil is (-3,-10,18,-6,1)/(12 dx) at offsets (-1,0,1,2,3) for
positive shift. Its Nyquist symbol is **−8/(3 dx)**, not a centered imaginary
symbol. KO is (1,-6,15,-20,15,-6,1)*sigma/(64 dx), **on every evolved row,
including B** (`CCZ4Cartoon.impl.hpp:162`). Each direction contributes
−sigma/dx at Nyquist. The compiled Chombo `LevelRK4.H:107–162` has the classical
RK4 weights and stage times; its stability polynomial is
R(z)=1+z+z^2/2+z^3/6+z^4/24.

The exact CAS witnesses in [t24-gauge-stability.json](t24-gauge-stability.json)
verify the stencil symbols, two-component beta/Gamma characteristic polynomial,
and |R(−1+i y)|^2−1. A separately evaluated 70-digit **four-stage recurrence**
checks the boundary and both sides. The reduced zero-shift, diagonal-metric,
two-direction Nyquist criterion is indeed
**S <= 4.9062387672032055**, where S=H22+0.75 H11. Without KO it is 6;
the single-rho-direction KO limit on H22 is **6.36217048021**.
These are reduced-model limits, not a complete AMR or variable-background proof.

For nonzero h12 the fastest Nyquist coefficient is
c_Gamma*(trace(H)+lambda_max(H)/3); H22+0.75 H11 is only its diagonal rho-mode
special case. Shift upwinding and eta move the RK4 argument left. At the saved
cell the Nyquist wave has Re(z)≈−1.03274, changing the bound. The resolved
129×129 wave-number scan brackets the decaying-branch transition, when all four
h components are scaled together, at **{s[0]:.9f} < S_limit < {s[1]:.9f}**.
This is a sampled numerical bracket, not a rigorous bound between wave numbers.
Thus **4.906 is the correct scalar reduction but is too permissive as this
cell's complete discrete criterion**. The saved S is about
{100*4.662576980611591/s[0]:.3f}% of that sampled transition.

The independent blockwise L12 active-valid maximum is S=4.662576980611591 at
(5243409,135), rho=0.057891845703125, dx=0.00042724609375. Its h11/h12/h22/hww
are 7.75956300347 / 0.0422448889773 / 0.219248673040 / 0.588659453868;
beta=({p['beta'][0]:.12g},{p['beta'][1]:.12g}), chi=0.000455879319582,
lapse=1e-12. All 28 values are in [t24-gauge-cell.json](t24-gauge-cell.json).
The prescribed variants use the measured cell, scaling all h entries by
S_original/S_target; this choice preserves metric shape but not unit determinant.
They are operator sensitivity tests, not physically initialized states.

{waves}

The unstable high-frequency branches at S=4.9/5.5/6.4 are dominated by **beta2
and Gamma2**, with ky at Nyquist and kx near Nyquist. The exact modes and complex
RK arguments are in [t24-gauge-stability.csv](t24-gauge-stability.csv).
The neutral near-zero-wave branches mean a stable principal scan's overall
maximum is 1; the Nyquist column retains the useful margin. Half dt, half
c_Gamma, sigma=0.5 and sigma=0 have no decaying-branch RK4 violation at the
saved metric. Sigma=0 is **an offline operator test only**; KO remains 1 in
every production/evolution test and no KO removal is proposed.

The **full** frozen cartoon matrix has positive-real low-frequency eigenvalues:
its overall maximum is 1.001982115 at k/pi≈(−0.0078125,−0.0546875) for the saved
cell. These are already growing eigenvalues of that homogeneous cylindrical
semi-discrete operator, not an RK4 amplification error of a decaying mode.
Freezing Cartesian background derivatives at zero loses the actual numerical
gradients and radial wave envelopes. Their interpretation as radial amplitude
transport, physical/local source growth or a real instability remains open.
S alone cannot bound them, and they are not silently removed from the table.

**Established:** the code has an alpha-independent fast shift/Gamma loop, its
actual discrete stability margin at step 152 is small, and the proposed larger-S
checkerboard is numerically unstable. **Not established:** a violation at the
first bad RK stage, growing checkerboard in the trajectory, or gauge CFL as the
cause. The step-153/154 and stage history is still needed. No frozen-grid scan
proves stability of subcycled AMR, exchanges, regrids or finite-alpha matter.

## 2. Numerical apparent-horizon searches

[t24-horizon-search.py](t24-horizon-search.py) loaded only the numerical step-152
checkpoint into the frozen T17 harness, with unchanged RHUnion/RHSurf and native
Lagrange<4> interpolation. Static and CTT paths were deliberately absent; no
advance or initialization was executed. Seeds 0.002,0.004,0.008,0.016,0.032,0.064,
0.096 M_i were chosen as a numerical radius sweep, independently at both
numerically tracked centres, not from a static trumpet. N48, chase 1, quota 1,
64 updates and the inherited stage-0 squared-expansion threshold 1e-7 were used.

{horizons}

None qualifies even the first stage; the minimum squared residual is
**0.2920168877**, and every search reports UPDATE_CAP. The trial extents in the
table are **unconverged surfaces**, never apparent-horizon extents. N96 and
later strict stages were consequently not run. The native nonzero exit **1**
is the bounded finder failure, not a launcher success, memory failure or NaN.
Cost including checkpoint load: **{h['resources']['wall_seconds']:.3f} s**,
peak RSS **{h['resources']['peak_rss_bytes']:,} B**; actual chase iteration
time is about 18.675 s across all 14 surfaces. See
[t24-horizon-search.csv](t24-horizon-search.csv) and its own JSON receipt.

There is therefore **no qualified apparent horizon from which to classify
rho=0.047–0.058 as inside or outside**. A failed bounded search does not prove
the absence of an individual or common horizon. The at-failure moving centres
also cannot be substituted for the step-152 centres during this frozen search.

## 3. Implemented opt-in joined NaN abort

`ems_safe_nan_abort=false` preserves the old NanCheck path.
`ems_safe_nan_abort=true` uses [SafeNanAbort.hpp](../../Source/BoxUtils/SafeNanAbort.hpp)
in `EMSBH2DLevel::specificAdvance`, after the existing projections and on valid
cells only. Parameters are `ems_nan_max_abs` (default 1e20, the old threshold),
`ems_nan_abort_exit_code` (default 86) and `ems_nan_abort_prefix`.
Detection collects a first local witness without a conditional OpenMP barrier.
After all workers join it closes/renames a per-rank JSON record with field,
cell, level, time, dx, dt, phase, threshold, exit code and all evolved IEEE
hexadecimal values, then directly calls MPI_Abort on Chombo's communicator.
The non-MPI fallback uses the same nonzero code. There is no abort collective
that could wait for a rank stranded in an exchange. If the record cannot be
written it prints an explicit write failure and still terminates nonzero.
The hook diagnoses the post-complete-RK check; it does not claim to identify
the first generating RK stage.

The isolated patch [t24-safe-abort.patch](t24-safe-abort.patch) applies only this
feature to d507438. It contains no later gauge or recording changes. Every
native unit was rebuilt against the new parameter layout: GRAMRLevel stores
SimulationParameters by value. The initial partial-object build was rejected
and corrected, not accepted as a regression.

Omitted/off/on native checks at t=0, after two regridding coarse steps and after
restart have **{nn:,} defined Float64 comparisons, zero bit mismatches**.
All evolved checkpoint/plot values, including saved evolution ghosts, and
valid diagnostic cells are compared. Only the previously demonstrated undefined
plot-diagnostic ghost slots are excluded. A NaN on the non-master OpenMP worker
preserves all 112 input Float64 bit patterns, closes the registered record,
and exits **86**; a finite input is unchanged and exits 0. No RHS, floor,
transfer, gauge, tagging or tracker operation is changed. Evidence:
[t24-controls-qualification.json](t24-controls-qualification.json),
[t24-control-bits.csv](t24-control-bits.csv),
[t24-control-resources.csv](t24-control-resources.csv).

The CH_MPI probe and the dependency's actual MPI SPMD unit compile with the
strict flags. **Every local one/two-rank attempt fails in MPI_Init**, exit 15,
at POSIX shared-memory bootstrap before the checker; the direct-abort path
is **not MPI-qualified**. [t24-mpi-abort-qualification.json](t24-mpi-abort-qualification.json)
states this failure rather than treating its nonzero code as an abort pass.
Changing the runtime SHM selector did not cure it. The controller's submission
worker still needs the two-rank/two-thread finite/default identity test and a
worker/rank NaN injection while the healthy rank waits in MPI, on an exclusive
node. No such cluster test was run here. The old cluster wait remains unproved.

The MPI choice is supported by the [MPICH API source](https://github.com/pmodels/mpich/blob/main/src/binding/mpi_standard_api.txt)
and its [abort notes](https://github.com/pmodels/mpich/blob/main/doc/mansrc/funcnotes.txt):
terminate the communicator's processes and perform the request outside worker
threads. That API argument is separate from the missing execution qualification.

## 4. Restart dt already honours the parameter

At `GRAMRLevel.cpp:747–749`, readCheckpointLevel reads dx and calls
computeInitialDt; `GRAMRLevel.cpp:543–549` computes dt=dt_multiplier*dx.
It does **not** use the serialized dt value. Chombo `AMR.cpp:674–686` then
populates its current/new level-dt arrays and calls assignDt. No new option is
needed or added. The existing parameter already implements the requested action.

[t24-restart-dt.csv](t24-restart-dt.csv) measures two native levels after restart
from t=0.03125 with initial coarse/fine dt=0.03125/0.015625. Omitted/off/on safe
abort with the same dt_multiplier=0.25 is bit-identical to the pinned restart.
With dt_multiplier=0.125 the stored dt is **0.015625/0.0078125**, exactly half;
coarse clocks are **0.046875 then 0.0625**, aligned with the old endpoint 0.0625
after two half steps. Static paths remain absent throughout these restarts.
The same source formula applies to every checkpoint level: on the binary
13-level grid half dt_multiplier=0.25 gives dt0=**0.4375** and dt12=
**0.0001068115234375** M_i.

A time-step diagnostic, if the controller admits it, must compare common
physical clocks: original steps 153/154 are half-dt steps 154/156; the old
event lies in half-dt step 157. It is not a repetition with identical step
numbers. No production parameter is changed silently in this turn.

## 5. Annular recorder contract and output budget

[t24-annular-recorder-design.md](t24-annular-recorder-design.md) fixes the
parameter interface, tracked annuli, stage/phase witness, per-level regrid
generation, gauge-labelled B rows, disabled-path and MPI qualification gates,
and checkerboard reading. It is **designed, not built**. No claim of MPI
qualification is made for an absent implementation.

On the sealed step-152 layout the .04–.09 annuli contain **27,972/46,442**
valid cells on L11/L12. The 15-double plus native coordinates/hole-mask format
costs **{payload/1e9:.6f} GB per coarse step**, uncompressed, at four stages and
stride 1. A declared 135.146484375–135.1875 window retains 96/192 consecutive
fine advances and costs **{window/1e9:.6f} GB** on those two levels, plus headers,
L10 witnesses and endpoint-inclusive records. Half dt doubles this payload.
Native layout counts, rather than a hypothetical complete radial shell, are
used; moving layouts can change the estimate. The 40 GB limit is retained.

## Mechanism ranking after exp-0026 and this audit

1. **Inner gauge/discretization growth, especially a shift/Gamma CFL crossing:**
   strongest new candidate; the alpha-independent beta2/Gamma2 mode and small
   sampled margin fit the failure's dominant beta2. Against: no measured stage
   crossing or checkerboard growth; the saved cell remains below the limit,
   and freezing does not resolve actual background/AMR effects.
2. **Inherited point-transfer or under-resolved inner-state damage:** the
   signed restriction and covered negative lapse are measured, and both legs
   inherit 152 point-transfer steps. Against: switching to legacy eliminates
   negative valid lapse without rescuing the evolution; it is not a proximate
   negative-lapse explanation. A clean earlier transfer-history control would
   be needed to distinguish inherited damage.
3. **Regrid/interpolation/interface amplification:** remains possible from
   frequent descendant regrids and unresolved instantaneous ghost/donor states.
   Against: the named final regrid changes no box count, and failure cells
   remain in nominal forced cores. Count-only evidence cannot prove or exclude
   a changed union, interpolation error or exchange defect.
4. **Tracker/tag nontermination or pure layout imbalance:** strongly disfavoured
   as the numerical initiating cause by distinct finite centres, finite max
   tagging even at coincidence, and explicit reproduced nonfinite states.
   The abort-control issue is a separate failure to terminate cleanly.

The next discriminating observations are step/stage S, both checkerboard rows,
floor sets and pre/post-regrid support, plus the existing half-dt restart
parameter. None requires removing KO, changing equations or using static data.
'''
    report=revised_reading(report,g,nn)
    (HERE/'t24-gauge-control-followup.md').write_text(report)
    # Keep the earlier numerical tables/report, but correct its current reading
    # and replace the outdated ranking rather than silently retaining rank 1.
    old=(HERE/'t24-merger-stall.md').read_text()
    start=old.index('**READY-EXCEPT');end=old.index('All source line references below',start)
    old=old[:start]+'''**READY-EXCEPT.** The shift/Gamma discrete-CFL candidate is supported by
exp-0026's history and two interventions, with the numerical and code-level
limits in [the current report](t24-gauge-control-followup.md). A unique initiating
operation and the exact original abort wait are unproved; strict horizon
membership and MPI safe-abort execution remain unqualified. The full local
replay is ended, not pending, and will not be retried. The original
[diagnostic restart design](t24-diagnostic-restart-design.md) below is retained
as historical evidence; the intervention has superseded its recorder build.
No evolution remedy is applied.

'''+old[end:]
    marker='## Exp-0026 and frozen-gauge follow-up — current reading\n'
    if marker in old:
        begin=old.index(marker);end=old.index('## A1.',begin);old=old[:begin]+old[end:]
    insert=f'''{marker}

The current reading is [t24-gauge-control-followup.md](t24-gauge-control-followup.md).
Both point/legacy restarts fail inside step 155; legacy's negative valid lapse
is gone from step 153, so it is not the proximate mechanism. Descendant levels
7–12 regrid every L6 step (0.013671875); the final 135.173828125 regrid changes
no box count. At failure the printed cells lie directly above x_p≈±0.1826,
not at the distance inferred using the t=133 positions. The active-S history
crosses the scalar limit 4.90623877 near t=134.4085, and both targeted G/D
interventions pass the original event. The shift/Gamma CFL candidate is
**supported**, not uniquely proved; measured-shift discrete limits are lower
than the reduced scalar rule. The expanded 48-row run-design table is in
t24-stability-limits.csv. The unchanged finder finds no qualified individual
or common surface in the numerical t=133 sweeps (N48/N96 midpoint radii .3–4);
horizon membership and qualified A/Q remain undetermined. The original
finder iteration logs were not included in ev1, limiting the historical
area-1000 diagnosis; flow-step sensitivity is locally measured.
The safe NaN abort is implemented/default off and serial-bit-qualified;
MPI qualification is blocked at MPI_Init. Restart already honours changed
dt_multiplier; no restart feature is added. The annular recorder build is
**cancelled**, its design only retained. These statements supersede earlier
design-only abort wording below.

'''
    idx=old.index('## A1.');old=old[:idx]+insert+old[idx:]
    title='## Candidate mechanisms, ranked without an evolution remedy'
    if title not in old:title='## Candidate mechanisms — ranking updated after exp-0026'
    begin=old.index(title)
    end=old.index('## ',begin+3) if '## ' in old[begin+3:] else len(old)
    old=old[:begin]+'''## Candidate mechanisms — ranking updated after exp-0026

The current ranked evidence is in [t24-gauge-control-followup.md](t24-gauge-control-followup.md).
Shift/Gamma discretization and the supported CFL crossing rank first; inherited
transfer/resolution damage remains, while negative covered lapse as the
proximate cause is excluded by the failed legacy control. Frequent regrid/
interpolation effects remain open. A unique final regrid trigger and pure
tagging nontermination are not supported. A unique cause is not claimed proved.

'''+old[end:]
    (HERE/'t24-merger-stall.md').write_text(old)
    # Isolate the implemented feature from the unrelated dirty T21/T22 work.
    patch=[]
    for rel in ('Source/BoxUtils/SafeNanAbort.hpp','Examples/EMS/SimulationParameters.hpp','Examples/EMS/EMSBH2DLevel.cpp'):
        before=ROOT/'source'/rel;after=ROOT/'controls/source'/rel
        b=before.read_text().splitlines(True) if before.exists() else []
        patch.extend(difflib.unified_diff(b,after.read_text().splitlines(True),fromfile='a/'+rel if before.exists() else '/dev/null',tofile='b/'+rel))
    (HERE/'t24-safe-abort.patch').write_text(''.join(patch))
    hs=read('t24-horizon-summary.json')
    peak=max(r['peak_RSS_bytes'] for r in hs['receipts'])
    status=dict(status='READY-EXCEPT',causal_verdict='GAUGE_CFL_CANDIDATE_SUPPORTED_NOT_UNIQUE_PROOF',safe_abort='implemented, default off, serial qualified, MPI blocked',
        annular_recorder='build cancelled; historical design only',restart_dt='feature cancelled; existing parameter already honoured',running_jobs=[],
        horizon='no strict surface found in bounded individual and common searches; membership undetermined',
        original_finder_logs='absent from collected ev1 diagnostics/shapes',
        default_comparisons=nn,bit_mismatches=0,peak_this_followup_RSS_bytes=peak,
        new_evolution_remedy_applied=False,initial_data_reads_after_t0=0)
    (HERE/'t24-status.json').write_text(json.dumps(status,indent=2)+'\n')
    readme=HERE/'README.md';text=readme.read_text();begin=text.index('## T24 — binary merger stall, d507438');end=text.find('\n## ',begin+1)
    if end<0:end=len(text)
    replacement=f'''## T24 — binary merger stall, d507438

**READY-EXCEPT.** [Current follow-up](t24-gauge-control-followup.md),
[saved-state report](t24-merger-stall.md), [stability script](t24-gauge-stability.py),
[amplification table](t24-gauge-stability.csv), [numerical horizon searches](t24-horizon-search.csv),
[native controls](t24-controls-qualification.json), [restart dt](t24-restart-dt.csv),
[isolated safe-abort patch](t24-safe-abort.patch), [MPI qualification](t24-mpi-abort-qualification.json),
[annular recorder design](t24-annular-recorder-design.md), [storage](t24-annulus-size.csv),
[status](t24-status.json), [manifest](t24-manifest.txt).

Both exp-0026 transfers complete steps 153/154 and fail inside 155; removing
negative valid lapse via the legacy switch does not rescue the inherited state.
L7–12 regrid every L6 step and the final regrid changes no count. Failure cells
lie above the contemporaneous moving punctures. The scalar Nyquist bound is
4.90623877; with the measured shift/metric the sampled transition is near
{s[0]:.8f}, against saved S=4.66257698. Larger-S beta2/Gamma2 checkerboards grow.
The full frozen cartoon operator also has positive-real source modes, so this
is not a complete AMR/continuum stability or causal proof. No numerical horizon
qualifies in the bounded 14-seed native search (minimum squared residual .2920).
The implemented safe abort is default off and gives {nn:,} bit-identical
defined Float64 comparisons, including restart; its MPI qualification is
blocked during local MPI_Init. Restart already recomputes dt from the parameter;
native half-step clocks align exactly. The annular hook is **designed, not built**;
full L11–12 stage payload is {payload/1e9:.3f} GB/coarse step, or {window/1e9:.3f} GB
over the declared short window on the retained layout. Follow-up peak RSS
{h['resources']['peak_rss_bytes']/1e9:.3f} GB; earlier aborted full replay peak
5.119 GB is retained. No full replay is repeated, no job remains pending,
no evolution remedy is applied and nothing is committed or run remotely.
'''
    replacement=f'''## T24 — binary merger stall, d507438

**READY-EXCEPT.** [Current report](t24-gauge-control-followup.md),
[saved-state report](t24-merger-stall.md), [frozen operator](t24-gauge-stability.py),
[RK4 rule](t24-stability-rule.py), [48-row limits](t24-stability-limits.csv),
[growth](t24-stability-growth.csv), [verified cluster evidence](t24-cluster-check.json),
[all horizon probes](t24-horizon-all-probes.csv), [finder receipts](t24-horizon-summary.json),
[native controls](t24-controls-qualification.json), [safe-abort patch](t24-safe-abort.patch),
[MPI qualification](t24-mpi-abort-qualification.json), [status](t24-status.json),
[manifest](t24-manifest.txt).

The shift/Gamma CFL candidate is **supported** by active S=4.67762/4.81877/4.96223
at steps 152–154 and successful one-parameter G/D interventions; it is not a
unique-mechanism proof or production admission. The scalar Nyquist weighted-S
limit is 4.90623877 at dt/dx=.5,sigma=1, and 25.44868192 at .25; measured shift
and metric lower the local discrete limit. At S=4.96/5.03 scalar gains are
1.015263/1.035815 per fine step, consistent with a delayed NaN. The complete
table covers all requested dt/sigma/c_Gamma values and reports its frozen-model
scope; no other decaying alpha-independent branch crosses first at .25 in
the sampled scan. Homogeneous cartoon source growth is reported separately.

No strict individual or common horizon is found on the numerical t=133 state:
individual 14-seed N48, common eight-seed N48/N96 .3–4 M_i, plus a smaller-chase
control. Best final squared common residual .0508586 remains above 1e-7;
trial A/Q are not qualified horizon measurements. The collected ev1 packet
lacks original common-finder iteration logs; its summary and warm-seed policy
are available, and the local flow-step sensitivity is measured. Membership
of rho=.047–.08 remains undetermined. Follow-up peak RSS {peak/1e9:.3f} GB.

The safe NaN abort is implemented, default off, with {nn:,} defined Float64
comparisons and zero mismatches, including regrid/restart. A worker NaN closes
the record and exits 86. MPI execution is blocked during local MPI_Init and
must be qualified separately before cluster use. The isolated d507438 patch
contains no later gauge or recording changes. Restart already honours changed
dt_multiplier: its proposed feature and the annular-recorder build are
**cancelled**; the old recorder design is retained as design only. L7–12 regrid
every L6 step; the named final regrid changes no count. Negative covered lapse
is not the proximate requirement, since legacy removes it but still fails.
No full replay is repeated, no job remains pending, no evolution remedy is
applied, nothing is committed and no remote operation is performed.
'''
    readme.write_text(text[:begin]+replacement+text[end:])
    # Preserve earlier manifest evidence and add all new packet/source/receipts.
    paths=set()
    for line in (HERE/'t24-manifest.txt').read_text().splitlines():
        if not line or line.startswith('#'):continue
        paths.add(Path(line.split(None,2)[2]))
    paths.update(p for p in HERE.glob('t24-*') if p.is_file() and p.name!='t24-manifest.txt')
    paths.update([readme,REPO/'Source/BoxUtils/SafeNanAbort.hpp',REPO/'Examples/EMS/EMSBH2DLevel.cpp',REPO/'Examples/EMS/SimulationParameters.hpp',ROOT/'controls/build-spec.json'])
    paths.update(ROOT.glob('controls/**/*.resources.json'));paths.update(ROOT.glob('controls/**/abort*.json'))
    paths.update(ROOT.glob('mpi-abort/**/*.resources.json'));paths.update(ROOT.glob('horizon/*.resources.json'))
    paths.update(ROOT.glob('horizon/params.txt'))
    for rel in ('horizon-common-pilot','horizon-common-N96','horizon-common-slow'):
        paths.update(ROOT.glob(rel+'/*.resources.json'));paths.update(ROOT.glob(rel+'/params.txt'))
        paths.update(ROOT.glob(rel+'/run.log'));paths.update(ROOT.glob(rel+'/progress.csv'))
        paths.update(ROOT.glob(rel+'/shape*.dat'));paths.update(ROOT.glob(rel+'/done.exit'))
    paths.update(Path(p) for p in read('t24-cluster-check.json')['verified_files'])
    paths.update(Path(p) for p in hs['own_receipts'])
    entries=[]
    for p in sorted(paths):
        assert p.is_file(),p
        with p.open('rb') as f:d=hashlib.file_digest(f,'sha256').hexdigest()
        entries.append(f'{d}  {p.stat().st_size}  {p}\n')
    (HERE/'t24-manifest.txt').write_text('# T24 current follow-up: no commit, no evolution remedy\n'
        '# Safe abort default off: serial qualified; MPI qualification blocked at Init\n'
        '# Annular recorder build cancelled; gauge CFL supported, unique cause not proved\n'
        '# No strict horizon found; membership undetermined; original finder logs missing\n'
        f'# Follow-up peak RSS {peak} B; earlier native replay peak 5119148032 B\n'
        '# SHA256 bytes absolute_path\n'+''.join(entries))
    print('T24_FOLLOWUP_SEALED',len(entries),flush=True)


if __name__=='__main__':main()
