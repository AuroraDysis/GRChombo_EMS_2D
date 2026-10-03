#!/usr/bin/env python3
"""Write the measured t=0 report from the qualified numerical ledgers only."""
import csv,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def read(name):return list(csv.DictReader((HERE/name).open()))
def number(x):return f'{float(x):.8e}'
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(x) for x in row)+' |' for row in rows])
q=json.loads((HERE/'t22-equilibrium-qualification.json').read_text())
assert q['status']=='E_HIGH_NATIVE_EQUILIBRIUM_AUDIT_PASS_WITH_MEASURED_K_ROUNDOFF'
census=read('t22-KTheta-census.csv');rates=read('t22-initial-lapse-rates.csv');summary=read('t22-equilibrium-summary.csv')
valid=[a for a in census if a['gauge']=='moving_puncture' and a['scope']=='valid']
def rate(g,l):return next(a for a in rates if a['gauge']==g and int(a['level'])==l)
kt=table(['level','valid cells','K peak','K cell RMS','Theta peak / RMS','K / (eps scale) peak'],[[a['level'],a['cells'],number(a['K_peak']),number(a['K_cell_RMS']),'0 / 0',f"{float(a['K_over_eps_contraction_scale_peak']):.6f}"] for a in valid])
lt=table(['level','MPG lapse peak','MPG lapse cell RMS','old lapse peak','old lapse cell RMS'],[[l,number(rate('moving_puncture',l)['preKO_lapse_peak']),number(rate('moving_puncture',l)['preKO_lapse_cell_RMS']),number(rate('experimental',l)['preKO_lapse_peak']),number(rate('experimental',l)['preKO_lapse_cell_RMS'])] for l in range(15)])
mt=table(['level','points','geometric Gamma peak / RMS','matter Gamma peak / RMS','full Gamma = MPG B peak / RMS','KO(Gamma) peak','KO(B) peak'],[[a['level'],a['native_samples'],number(a['Gamma_geometric_peak'])+' / '+number(a['Gamma_geometric_RMS_sampled']),number(a['Gamma_matter_peak'])+' / '+number(a['Gamma_matter_RMS_sampled']),number(a['Gamma_full_preKO_peak'])+' / '+number(a['Gamma_full_preKO_RMS_sampled']),number(a['KO_Gamma_peak']),number(a['KO_B_peak'])] for a in summary if a['gauge']=='moving_puncture'])
st=table(['level','old shift peak','MPG shift peak','old offset B rate peak','KO(lapse) peak','KO(shift) peak'],[[l,number(next(a for a in summary if a['gauge']=='experimental' and int(a['level'])==l)['preKO_shift_peak']),'0',number(next(a for a in summary if a['gauge']=='experimental' and int(a['level'])==l)['preKO_B_peak']),number(next(a for a in summary if a['gauge']=='moving_puncture' and int(a['level'])==l)['KO_lapse_peak']),number(next(a for a in summary if a['gauge']=='moving_puncture' and int(a['level'])==l)['KO_shift_peak'])] for l in range(15)])
ghosts=[a for a in census if a['gauge']=='moving_puncture' and a['scope']=='stored_ghost_copies']
gt=table(['level','stored ghost copies','K peak','K cell RMS','K / (eps scale) peak'],[[a['level'],a['cells'],number(a['K_peak']),number(a['K_cell_RMS']),f"{float(a['K_over_eps_contraction_scale_peak']):.6f}"] for a in ghosts])
peak=max(valid,key=lambda a:float(a['K_peak']))
text=f'''# T22 measured initialization: the exact-zero K premise is false

The E-high initialization audit **passes with explicitly measured K roundoff**, not exact-zero K. The consult's premise that the setter supplies exact zeros for K and Theta was false for K. Theta is exactly zero, and the new driver B is exactly zero on all valid cells and stored ghosts. K is nonzero in every valid cell of both gauges, not only at the puncture or in ghosts. It peaks at {q['K_measured_peak']:.10e} M^-1 on level {peak['level']} at (x,y)=({float(peak['K_max_x_M']):.10e},{float(peak['K_max_y_M']):.10e}) M, R={float(peak['K_max_radius_M']):.10e} M. Its maximum scaled residual is {q['K_scaled_eps_peak']:.8f} eps times the local absolute tensor-contraction scale. This is Float64 contraction/cancellation roundoff of the analytically traceless isolated, unboosted construction, not a physical nonzero mean curvature.

The initial setter builds K_ij=P^2 k(3 n_i n_j-delta_ij) at [EMSBH_trumpet_read.impl.hpp](../../Source/InitialConditions/EMSBH/EMSBH_trumpet_read.impl.hpp#L223). At zero rapidity the spatial metric is P^2 delta_ij, giving trace 3k(n.n-1)=0 for a unit direction. The exact rational witness and independent 50-digit numerical inversion are in [t22-trace-qualification.json](t22-trace-qualification.json): exact residual zero and worst scaled numerical residual 1.20490e-50. These are algebraic witnesses of the construction, not static-reference comparisons or positive-time reads. The setter initializes a trace accumulator at line 419, **computes** the metric contraction at lines 420–422, and stores it as K at line 426. It does not assign K as a literal zero. VarsTools::assign(vars,0) at line 424 supplies Theta's zero, which is never overwritten.

The original qualification assertion required every K/Theta cell to be exactly zero. Its exit-1 receipt is retained at `/private/tmp/ems-t22/evidence-pipeline/` and the original qualification log. The controller-authorized replacement explicitly requires Theta=0 and |K|<=64 eps S locally, where S=|h^11 A11|+2|h^12 A12|+|h^22 A22|+|Aww/hww| from the stored numerical state. K itself is never rounded, clipped, subtracted or overwritten. The 64-eps gate was stated before measuring the full normalized census; the measured maximum is 3.720384. Ghost copies are checked and labelled separately. This replacement concerns initialization qualification alone: the registered 10.5 M activity reading, masks, clocks, fivefold significance, temporal rule and budgets remain frozen.

## Every-level valid-cell K and Theta

{kt}

The table uses unweighted native-cell RMS. [t22-KTheta-census.csv](t22-KTheta-census.csv) additionally reports cartoon coordinate-volume RMS (weights |y|), every level's peak location and radius, both gauge labels, the separate puncture-neighbourhood census (four rows by eight columns around the puncture), and stored ghost-copy counts. The two initial physical states agree on **{q['physical_initial_identity']['Float64_values']:,} Float64 values**, zero bit mismatches, excluding only driver B1/B2. Both hierarchies have the exact collected E-high box lists: levels 0–14, 546 boxes, 3,815,424 valid cells. The initial states and all sampled RHS pieces are finite.

## Every-level measured lapse rates

The new gauge's all-valid-cell pre-KO lapse rate is the native class's -2 alpha(K-2Theta), not an assumed zero. Its global peak is {q['preKO_MPG_lapse_measured_peak']:.10e} M^-1; the formula check has zero residual. The old gauge's rate is native upwind beta.grad(alpha)-1.8 alpha(K-2Theta), with its historical coefficient unchanged. The new pre-KO shift rate is mu B=0 exactly. Units of the lapse rates are M^-1.

{lt}

[t22-initial-lapse-rates.csv](t22-initial-lapse-rates.csv) includes coordinate-volume RMS and cell counts for every level. The all-valid-cell rates use the checkpoint's stored stencil halos and actual native derivative/gauge classes; no new initialization or evolution is performed. At level 14 the new lapse cell RMS is 3.36817070e-17 versus old 4.34323102e-3. Ordinary KO(lapse), which is unchanged and can dominate this roundoff-level pre-KO gauge row, is reported separately below.

## Complete discrete Gamma and driver pieces

The new B pre-KO RHS equals the **full discrete Gamma RHS bit for bit** at every native sampled point, including the EMS matter momentum term. B=0 and advection is disabled. The geometric column is the native RHS with Newton coupling zero on a private current-state evaluation; the matter column is independently reconstructed from the native EMS stress tensor. Their sum agrees with full Gamma within the original arithmetic budget (measured budget fraction zero). Ordinary KO(B)=0 initially because B=0; KO(Gamma) remains an independent Gamma row and is never inserted into B.

{mt}

These Gamma/B RMS values are explicitly **sample RMS**, not full-mask norms: 601 native points per gauge across all 15 levels, comprising 32 puncture-neighbourhood cells plus retained in-box axis/equator/diagonal anchors through R=0.02 M. A finite-difference Gamma residual is expected; it is neither asserted zero nor removed. On level 14 its peak/RMS are 1.76946008 / 0.377251603, while KO(Gamma) peaks at 1.26498314e-7. Geometric and matter peaks separately reach about 3.385912 and partly cancel. Full signed pointwise values and coordinates are in [t22-equilibrium-native.csv](t22-equilibrium-native.csv); all per-level sample counts, component-combined RMS, geometric/matter/KO and old-gauge rates are in [t22-equilibrium-summary.csv](t22-equilibrium-summary.csv).

{st}

The old offset driver and the new differential driver are different variables. The table reports their labelled source rows separately; it forms no cross-gauge B ratio. The initial Gamma, lapse KO and shift KO pieces coincide because their physical initial fields coincide.

## Stored ghosts and the native storage test

{gt}

Ghost copies are duplicated across boxes and are not physical-volume counts. Theta and new B are exactly zero throughout these copies; the largest ghost K is 4.45645463e-15 M^-1, scaled maximum 3.492669 eps. Their measured values pass the same numerical contraction bound, although the required valid-cell gate is kept distinct.

After replacing the K assertion, the old audit's internal total-RHS check exposed 92 one-component flags per gauge. [VarsTools.hpp](../../Source/utils/VarsTools.hpp#L41) maps both upper and lower symmetric tensor entries into the same evolved component; storage retains the last entry. The old test compared **both aliases** with that stored value. All flags are the overwritten upper A12 alias (component 7), not a mismatch of an evolved field. [t22-tensor-aliases.csv](t22-tensor-aliases.csv) retains every upper/lower discrepancy. The repaired Tests-only check compares exactly the 28 serialized, last-write values with the native stored output, still **bitwise with no tolerance**: 16,828 components per gauge, zero bit mismatches. [t22-native-storage-check.csv](t22-native-storage-check.csv) gives the counts. This does not change any production arithmetic or weaken the pinned default-path regression.

## Receipts and reproducibility

Both completed native initial runs are verified from their own actual child status, no gate event, finite native outputs, time-zero checkpoint and hierarchy census, not merely the launcher exit. Experimental: 45.709 s, peak RSS 1,642,774,528 bytes; moving_puncture: 43.159 s, peak RSS 1,369,243,648 bytes. Their logs say `T22_NATIVE_INITIAL_AUDIT_COMPLETE; no advances`. The sandboxed `time` sysctl warning has no scientific or child-exit meaning. [t22-equilibrium-resources.csv](t22-equilibrium-resources.csv) preserves these receipts. The offline trace/rate process stays below 0.6 GB, and the continuation qualification has its own receipt and done marker. All native runs use two OpenMP threads, below the four-thread limit. No completed evolution is repeated, no SSH/cluster operation or commit is made, and static data enter no evolution after t=0.

Regenerate with `t22-trace-audit.py compile-rates`, `t22-trace-audit.py analyze`, `t22-verify.py equilibrium`, then this script. Wrap the first three in `t13-run.py --measure` with the 6 GB per-process cap; run commands use absolute script paths because the measurement wrapper changes directory. Only stored numerical checkpoints are read. The design is updated only in its qualification note; `/private/tmp/ems-t22/design-frozen.md` and its registered SHA remain immutable.
'''
(HERE/'t22-equilibrium-report.md').write_text(text)
print('T22_MEASURED_EQUILIBRIUM_REPORT_WRITTEN')
