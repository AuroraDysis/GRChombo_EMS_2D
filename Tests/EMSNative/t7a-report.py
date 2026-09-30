#!/usr/bin/env python3
"""Append completed T7 A evidence without changing its registration."""
import csv
from pathlib import Path
HERE=Path(__file__).resolve().parent;TMP=Path('/private/tmp/ems-t7a')
def rows(name):return list(csv.DictReader((HERE/name).open()))
def num(x):return f'{float(x):.6g}'
def table(headers,rs):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rs])
pk=[r for r in rows('t7a-common-packet-ratios.csv') if r['ray'] in ('axis_plus','diagonal') and int(r['level'])==int(9-float(r['face_M']))]
packet=table(['face / ray','field','T6 peak','clock peak','space peak','clock/T6','space/clock','p_h'],[[f"{r['face_M']} / {r['ray']}",r['field'],*[num(r[k]) for k in ('T6_peak','clock_peak','space_peak','clock_over_T6','space_over_clock','spatial_apparent_order')]] for r in pk])
wake=[r for r in rows('t7a-ratios.csv') if r['field']=='Ham' and r['band'] in ('wake_inner','wake_outer')]
waketable=table(['face','t/M','collar','T6 RMS','clock RMS','space RMS','clock/T6','space/clock','p_h'],[[r['face_M'],r['time_M'],r['band'],*[num(r[k]) for k in ('T6_face3','clock','space','clock_over_T6','space_over_clock','spatial_apparent_order')]] for r in wake])
width=[r for r in rows('t7a-width-summary.csv') if (r['level']=='6' and r['ray'] in ('axis_plus','diagonal')) or (r['level']=='5' and ((r['ray']=='axis_plus' and abs(float(r['time_M'])-3.75)<1e-8) or (r['ray']=='diagonal' and abs(float(r['time_M'])-5.5)<1e-8)))]
widthtable=table(['run','face / ray','t/M','half-prominence width/M','fine / receiving cells','end-limited prominence?','closed half-height width/M','receiving cells (half-height)'],[[r['case'],f"{r['face_M']} / {r['ray']}",num(r['time_M']),num(r['width_M']),num(r['fine_cells'])+' / '+num(r['receiving_cells']),r['domain_end_baseline'],num(r['halfheight_width_M']),num(r['halfheight_receiving_cells'])] for r in width])
contr=[r for r in rows('t7a-profile-contrasts.csv') if r['method']=='tensor' and r['ray']=='axis_plus' and r['time_M']=='3.0' and r['field'] in ('lapse','chi','K','Theta','Gamma1','shift1','B1','phi','Pi','Ex','Xi')]
contrast=table(['field','clock range','clock−T6 RMS','space−clock RMS','space contrast / range','P6/P8 max difference'],[[r['field'],*[num(r[k]) for k in ('clock_range','clock_vs_T6_RMS','space_vs_clock_RMS','space_vs_clock_normalized_RMS','P6_P8_max')]] for r in contr])
sources=[r for r in rows('t7a-source-summary.csv') if r['ray'] in ('axis_plus','diagonal') and ((r['level']=='6') or (r['level']=='5' and r['face_M']=='4.0'))]
sourcetable=table(['run / face / ray','first 1% source t/M','peak source t/M','pre-KO Theta RHS','KO Theta RHS','total Theta RHS','0.5 alpha Ham','RK physical ΔTheta','RK KO ΔTheta','RK total ΔTheta'],[[f"{r['case']} / {r['face_M']} / {r['ray']}",*[num(r[k]) for k in ('first_1pct_source_time_M','peak_source_time_M','Theta_pre_KO_at_peak','Theta_KO_at_peak','Theta_total_at_peak','half_alpha_Ham_at_peak','net_physical_increment_at_peak','net_KO_increment_at_peak','net_total_increment_at_peak')]] for r in sources])
phase=table(['run / ray','curvature peak r (P6)/M','peak r (P8)/M','P6/P8 phase difference/M'],[[r['case']+' / '+r['ray'],*[num(r[k]) for k in ('radius_P6_M','radius_P8_M','P6_P8_phase_difference_M')]] for r in rows('t7a-phase-summary.csv') if r['ray'] in ('axis_plus','diagonal')])
restrict=table(['run / face / ray','peak jump t/M','covered-coarse ΔTheta','fine Theta support max before copy'],[[f"{r['case']} / {r['face_M']} / {r['ray']}",*[num(r[k]) for k in ('peak_jump_time_M','Theta_jump_at_peak','fine_Theta_support_absmax_at_peak')]] for r in rows('t7a-restriction-summary.csv') if r['ray'] in ('axis_plus','diagonal')])
sens=[r for r in rows('t7a-ratios.csv') if r['band']=='packet' and ((r['face_M']=='3.0' and r['time_M']=='4.0') or (r['face_M']=='4.0' and r['time_M']=='6.0')) and r['field'] in ('Ham','Theta')]
sens_table=table(['face','t/M','field','T6 collar RMS','clock collar RMS','space collar RMS','clock/T6','space/clock'],[[r['face_M'],r['time_M'],r['field'],*[num(r[k]) for k in ('T6_face3','clock','space','clock_over_T6','space_over_clock')]] for r in sens])
audit=rows('t7a-run-audit.csv');audit_table=table(['run','exit','sigma','stage frames','wall / min','evolution RSS / GB','captured min chi / lapse','captured floor cells','GaussB nonzero'],[[r['case'],r['exit'],r['sigma'],r['stage_frames'],num(float(r['wall_seconds'])/60),num(float(r['evolution_peak_rss_bytes'])/1e9),num(r['chi_min_captured'])+' / '+num(r['lapse_min_captured']),r['floor_cells_captured'],r['GaussB_nonzero']] for r in audit])
seams=rows('t7a-common-seams.csv');seamrows=[]
for ray in ('axis_plus','diagonal'):
 for seam in (1.5,2.25):
  for case in ('T6_face3','clock','space'):
   rr=[r for r in seams if r['case']==case and r['ray']==ray and float(r['seam_coordinate_M'])==seam and r['field']=='Ham' and float(r['time_M'])<=3.0+1e-9]
   seamrows.append([case,ray,seam,num(max(float(r['peak_abs']) for r in rr))])
seamtable=table(['run','ray','common plane/M','max Ham through t=3 M'],seamrows)
text=f'''
### A — completed run analysis (controller resumption)

**READY-EXCEPT — mixed decision.** Receiving resolution is the dominant improvement for the first packet peaks and the interface-excluded bulk wake. Halving the clock changes those peaks by about 1–3%, whereas halving h reduces them by factors about 2–4. Bulk collar Ham orders are about 1.6–2.4 after passage, and a fixed diagonal wake gives about fourth order. The strict resolution-only acceptance is **not** certified: some later interface collars have material clock sensitivity, the incident curvature width changes appreciably with h, and the recorder establishes the direct evolved-Theta source without uniquely separating the origin of the Hamiltonian input error into ghost space/time closure versus earlier metric evolution. Temporal coupling remains a priority for the clock-sensitive collars. The finite-width premise is neither proved nor killed by this single spatial contrast.

The registration above is preserved verbatim, including its historical pending statuses. Both controller-supplied markers now read zero. No new evolution, static-profile evaluation, gauge/equation change, production-source modification or commit occurs in this analysis. Every evolved comparison uses current saved fields or another numerical run. The standalone [T7AReplay.cpp](T7AReplay.cpp) invokes the unchanged native kernel on recorded operands; it is never called by evolution.

#### Audit and actual hierarchy

{audit_table}

Both runs have point transfer, sigma=1, fixed regrid intervals, the t>0 initial-data guard and global runtime NaN checks enabled. All 7,399,296 stage frames and all 438 per-level snapshot times decode to finite Float64 values. All 28 evolved variables and five constraints are present. No floor activation occurs in the captured face regions; global runtime NaN checks completed. These face recordings cannot count floor hits in unsaved puncture cells. Parameter-default messages are benign defaults or inactive features and are retained in [t7a-default-parameters.csv](t7a-default-parameters.csv), with logs/parameters hashed in [t7a-run-audit.csv](t7a-run-audit.csv).

Actual unions, unchanged from the registration, are [-a,a]×[0,a], a=(256,64,32,16,8,4,3) M. The 3 M face is level 6→5; the 4 M face is level 5→4. Clock spacings there are (1/48,1/24) and (1/24,1/12) M; space halves each at the same corresponding stage times. [t7-layout.csv](t7-layout.csv) / [t7-boxes.csv](t7-boxes.csv) retain actual boxes. The different same-level seams are compared on common physical bands below, rather than by per-box cell classes.

#### Independent events and sampling qualification

[t7a-events.csv](t7a-events.csv) searches each run, face, ray, level and Ham/Theta field over the whole interval for log-prominence ≥0.5 decades. No expected propagation time enters selection. The growth time maximizes log growth within its basin with both endpoints ≥1% of that peak; this excludes the trivial zero-Theta startup jump. Rising endpoint candidates meeting the same prominence threshold are separately marked **censored**. Thus the 4 M corner near t=6 is a creation event, not a completed peak certificate. Replay windows span independently selected Ham growth through the selected peak, including the censored endpoint.

[t7a-analyze.py](t7a-analyze.py) reuses the lossless T7 decoder and T6 readers. Six/eight-node radial point interpolation is exact to degrees 5/7 on smooth test polynomials. On the diagonal, saved native points lie on the identical physical ray x=y, making this a qualified one-dimensional comparison. Supported face collars additionally use six/eight-node tensor sampling at fixed physical (x,y), with current parity-filled axis values. Polynomial/parity checks pass. No static target is subtracted.

Outside complete collars, axis/equator tubes contain too few transverse rows for a full six/eight-node tensor stencil. Their radial native-first-row profiles have an h/2 transverse offset, explicitly labelled in [t7a-profile-contrasts.csv](t7a-profile-contrasts.csv); radial P6/P8 differences do not bound that transverse error. At the receiving side of the clock 4 M face the saved collar is also too narrow for a full tensor-eight stencil. No tensor qualification is claimed there. Native operand replay itself has complete support and is unaffected by these profile-sampling limits. Counts/missing support and interpolation sensitivity are in [t7a-sampling-sensitivity.csv](t7a-sampling-sensitivity.csv). One-sided stencils stay on the selected level; covered coarse values are excluded.

The following fixed-axis tensor comparison is at t=3 M, 2.825<r<2.98 M, with 60 common samples. Range measures variation of the current clock profile, not variation from an equilibrium target.

{contrast}

The leading lapse, chi, K, shift, driver-B and EMS profiles are preserved closely in field value. Gamma changes by RMS 0.087% of this profile range, well above the interpolation uncertainty. Theta and Xi are small constraint/cleaning signals and change by large relative factors; they are not used to normalize the incident gauge amplitude. [t7a-profile-metrics.csv](t7a-profile-metrics.csv) retains amplitude ranges and P6/P8 sensitivity for every saved variable on all three rays, while the private ray caches preserve full profiles at every snapshot time.

![Incident fields on a fixed diagonal ray](figures/t7a-incident-profiles.png)

The displayed incident values nearly coincide, while the curvature flank sharpens under refinement and the small incident Theta error drops. These are distinct measurements; field-value agreement alone does not establish resolved derivatives.

The curvature-peak phase at t=2.5 M is sampled on the same radial mesh (step 1/384 M). Clock and T6 agree in peak position. Space shifts the diagonal peak by −0.007813 M with P6 and −0.018229 M with P8; the coarse peak itself has a 0.013021 M interpolation sensitivity. This prevents a stronger claim of phase equality. The axis rows retain the transverse offset limitation stated above. [t7a-phase-summary.csv](t7a-phase-summary.csv) also records both curvature amplitudes and sampled widths; policy widths below use native fourth-order derivatives rather than differentiation of an interpolated profile.

{phase}

#### Identity, physical widths, and receiving cells

Before the first face, native Gamma-curvature tracking gives clock axial speed 0.911765±0.003389 and space 0.911458±0.002121 M/M. Their local frozen longitudinal shift speeds average 0.906587/0.907118, versus lapse 0.551341/0.558179, transverse shift 0.774007/0.774567 and light 0.291783/0.297056 on the same samples. Diagonal measured speeds are 0.908147/0.908580. The incoming feature is carried principally by Gamma and shift curvature, with lapse/K response and smooth driver-B forcing; this supports the longitudinal shift/gauge-family interpretation. It is not a full characteristic projection or a claim of a pure eigenmode. Characteristic formulas retain T6's conditional [CAS witness](../../scripts/cas/t6c-evidence.md): algebraic identity — production physics not certified.

Between faces, the first contiguous outward ridge gives axial speeds 0.910714±0.023053 / 0.954167±0.011393 and local longitudinal speeds about 0.95256/0.95246. The later strongest peak can instead be a stationary wake; it is excluded from that propagation fit. [t7a-speed-fits.csv](t7a-speed-fits.csv) records fit intervals and uncertainties.

Widths below use native fourth-order radial second derivatives, with linear crossing locations. Half-prominence is the T6 definition. When its baseline is limited by the saved level end, a separate **closed half-height crossing relative to zero curvature** is retained as a sensitivity diagnostic; it does not silently replace the T6 definition. Counts divide radial width by isotropic local h. Along the diagonal, counts per native diagonal step/normal face direction are smaller by sqrt(2).

{widthtable}

At t=2.5 M, T6's axial half-prominence width was 0.11726318 M; clock is 0.11726311 M. Space is 0.08575413 M, ratio 0.73130, with receiving count increasing 2.814→4.116. Diagonal widths are 0.10689131→0.08042531 M, ratio 0.75240, counts 2.565→3.860. Thus width does **not** simply halve with h and keep a constant cell count, but it also has not settled to a resolution-independent value. The curvature peak steepens by about 34–36%. Further refinement is required to establish a finite-width asymptote.

At the 4 M axial face, the closed half-height width at t=3.75 M gives only 1.384/1.932 receiving cells, and the diagonal at t=5.5 gives 1.079/1.532. The associated half-prominence estimates are end-limited. Even this refined run has a poorly sampled receiving side at the next face. The data do not provide any qualified resolution at the 8 M face, which the principal front has not crossed by 6 M.

#### Packet amplitudes and wake ratios

The next table uses identical physical windows and the matching 1/12 M pre-parent-restriction native sampling cadence: half-width 1/24 M around the 3 M face, 1/12 M around 4 M, including the corresponding axial or corner tube. Peaks are selected independently over the full available interval. These fixed windows can clip a packet that shifts with h; broader-window and tensor sampling sensitivities remain in the other tables. The 4 M corner values at 6 M are censored endpoint maxima. p_h=log2(clock/space) is an **apparent two-run constraint amplitude order**, not a three-grid field self-convergence proof.

{packet}

For the bulk wake, define e=max(|x|,y). Inner and outer collars are a−0.1875<e<a−1/12 and a+1/12<e<a+0.1875 M respectively. The exclusion is a fixed physical distance from the interface on every grid. Every uncovered cell in the saved full collar is used, with the physical coordinate-volume weight 2πy h_level². The box decomposition therefore changes quadrature resolution, not the selected physical region. These are local wake norms, not the full exp-0019 far-mask norms.

{waketable}

The 4 M rows at t=4 precede the main crossing and their high apparent orders do not qualify a post-passage wake. At t=6, receiving-side outer collar orders are 1.680 at 3 M and 2.429 at 4 M, while clock changes are about 1.4–1.5%. Inner-side errors are reduced with orders 2.029/1.637. This is evidence that the returned/upstream error is reduced too, but that norm does not isolate a reflected characteristic mode. A pure reflection coefficient is not claimed.

A separate exact-diagonal ray wake, 3.125<r<3.5 M at t=4 and 3.125<r<5 M at t=6, has clock/T6 Ham ratios 1.000242 / 0.991497 and space/clock 0.070742 / 0.056702, apparent orders 3.8213 / 4.1405. [t7a-ray-wake-ratios.csv](t7a-ray-wake-ratios.csv) also retains P8 and other fields. Angular and interface sampling therefore matter: neither a universal low-order wake nor universal fourth-order recovery follows from one norm.

![Packet and local wake histories](figures/t7a-packet-wake.png)

Clock and T6 nearly overlap at the first peaks and in the bulk collars. Space suppresses the packets and wake but does not remove continued error production near the faces.

There is a material clock effect in later **interface-including** collars, even though the dominant peaks and bulk norms hardly change. All three columns below have the same broad physical collar |e−a|<1/12; T6 has 1/4 M full plots, but t=4 and 6 are exact matches, while the new snapshots are level post-step recordings. Diagnostic ghost-refresh ordering can affect Ham, so Theta is reported alongside it.

{sens_table}

These rows prevent a blanket declaration of the registered resolution-only win. Absolute clock/spatial RMS-change ratios are 0.988/0.737 at the 3 M face at t=4 (Ham/Theta), and 1.180/0.800 at the 4 M face at t=6. They prioritize temporal coupling/order for the remaining interface-local effect, without demonstrating a particular dense-output bug. Their absolute signal is smaller than the main packet, but the clock effect is comparable to the spatial effect there. No threshold was changed to call this a clean pass.

#### Common physical seam regions

The shared 1.5 M planes and space-only 2.25 M planes are compared in **the same** physical ray bands of half-width 1/12 M in every run. No per-box sample is compared with a differently placed box. The following prefix is before the independently selected first-face large source; all-time histories, including later returned errors, are retained separately.

{seamtable}

The changed internal seams do not retain a first burst comparable to the 1e-4 face packet. Later errors can travel back through them; their later presence is not proof of creation at a same-level seam. [t7a-common-seams.csv](t7a-common-seams.csv) records common-region histories.

#### Actual operation replay and causal limit

All captured actual pre-KO/KO RHS values in the independent event windows replay **bit-identically** with the standalone native production kernel. Spatial point prolongation from the saved actual dense parent operands agrees within 1.77e-15 in max(1,|value|) units, including parity-filled supports; stage times agree within 2.97e-16 M. Start fractions are 0 or 0.5; stage-time fractions are 0,0.25,0.5,0.75,1, with the two half-time stages retained separately. This verifies capture and arithmetic, not truncation accuracy on a poorly sampled incident feature.

The first **direct evolved-Theta creation** is the physical RHS reading the ghost-filled stencil, followed by the RK update. The already-filled stage input has a large Hamiltonian residual before trace projection. Near the positive source lobe, the physical Theta RHS closely equals 0.5 alpha Ham; KO opposes it by about 13–16%. The table reports physical/KO/total at the same source-peak cell and time, plus the largest net actual RK4-step increment. The two peak selections need not coincide; their times are retained in [t7a-source-summary.csv](t7a-source-summary.csv). Rows for the 4 M corner are censored by the 6 M stop.

{sourcetable}

The actual four-stage Theta increment equals dt/6 times the native (1,2,2,1)-weighted physical and KO rates, with maximum error 1.50e-21. This is an evolved increment, not an endpoint diagnostic-ghost change. Stage/update/end projections change Theta by exactly zero; stage trace removal changes Ham by at most 1.74e-18. Captured chi/lapse floors never activate. Restriction changes covered coarse cells only, and replays even the outermost ghost-dependent row/corner within 1.55e-15; the fine Theta packet already exists in current fine support before copying it. It cannot be the contemporaneous first creator of that fine packet, although earlier restriction can feed later coarse ghost history. Actual covered-cell jumps are below; these are changes per copy operation, not uncovered-cell evolution increments. [t7a-restriction-summary.csv](t7a-restriction-summary.csv) retains all rays, before/after values and replay errors.

{restrict}

![Physical RHS, KO, and actual RK4 increments](figures/t7a-operation-replay.png)

**Measured:** pre-KO physical RHS supplies the dominant direct Theta emission; direct Theta KO is dissipative at its positive source peak, and projections/floors do not create Theta. **Still unresolved:** the original cause of the Ham error already present in that stage input. The recorder starts after ghost fill and stores dense parent operands, not a pre-fill ghost state or the raw parent RK stage vectors. It also stores RHS on the target core, not RHS on every neighboring cell needed for a Hamiltonian derivative of the update. Therefore an exact before/after ghost-fill Ham budget, a unique spatial-versus-temporal ghost-error split, and KO-on-metric versus physical-metric contributions to creating Ham cannot be recovered. The data rule out those overclaims; agreement with the interpolation formula is not a proof that its boundary closure is accurate on this wave.

#### Implication for exp-0019 and the 100 M E layout

The measured reference mechanism explains why adding puncture-only levels or moving one face is insufficient: a sharp outgoing gauge feature loses receiving cells at each outward factor-two transition and excites constraints there. It supports receiving resolution as a remedy candidate for the reference far-mask order 1–1.7, while the remaining stage-local sensitivity needs its own qualification. It does not establish a universal roundoff floor or a proven asymptotic fix.

E's exp-0019 far Ham/Mom are real current-field residuals, but low/mid/high have different physical faces. In particular level-4 faces are 5.6875 / 5.541667 / 5.444444 M, with different corners inside the far shell; its order near zero is confounded by geometry/phase and cannot be assigned solely to this reference mechanism. Exp-0020's common faces remove that confound; no running exp-0020 data are accessed here. E-specific incident widths and charge budgets remain necessary.

For a prospective 100 M chain, use **E's measured evolved** width at each significant face, h_receiving≤w_E/N_qualified, and common physical face unions across rungs. This tranche qualifies no universal N: the tested first-face half-prominence counts are only about 2.6–2.8→3.9–4.1, and the next receiver has fewer cells still. Ten/twelve cells is a test target, not an admission certificate. Audit both faces and diagonal corners at every rung the disturbance crosses; for exp-0019 this includes the inner rungs and exterior faces near 2.7–2.84, 5.44–5.69, 10.89–11.38, 21.78–22.75 and 43.56–45.5 M. Their receiving levels are 4,3,2,1,0 respectively. The common exp-0020 layout must use its own actual unions. Retain the near-hole hierarchy but qualify a transport region and its ancestors, rather than using the reference width as an E constant. No 100 M production admission follows from the present pair.

**Cheapest next disambiguating test, proposed only:** one more clock leg at h0=4/3 M, dt0=1/48 M (CFL=1/64), the same fixed faces and sigma=1, to 6 M. Approximate local wall is 34 min from the measured clock wall; RAM remains roughly the current clock bound. Save pre-fill and post-fill ghosts, raw parent RK stage vectors, and RHS on the full Hamiltonian dependency halo, while preserving the opt-in bit controls. Freeze event selection, physical masks and profile sampling before launch. Pass the clock-negligible condition only if every relevant packet/wake/interface-collar norm changes by ≤5% versus the present clock and by <20% of the clock→space contrast, with incident profiles and phase within their sampling uncertainty. Kill a **resolution-only** interpretation if a remaining clock change is ≥50% of the clock→space contrast in any registered region; 5% alone is not a production accuracy budget. Intermediate results remain unqualified. This cheap test addresses the material temporal exception; it cannot by itself prove a finite-width continuum limit or admit E production.

If that clock condition passes, the next necessary spatial qualification is the fixed-face h0=1/3 M, dt0=1/48 M rung. Require stable incident width (ratio ≥0.85), rising receiving count, and packet plus interface-excluded wake apparent order ≥3 on the new pair, with P6/P8 uncertainty <10% of each measured contrast. Kill the finite-width layout premise if width again approaches halving with h (ratio ≤0.6) and cell count stays nearly constant; kill resolution alone as sufficient if a demonstrably resolved incident feature still gives order <2. The extrapolated full 6 M local wall is about 7.5 h, so the clock test is cheaper. No test is launched. E's common-face admission analysis and qualified horizon budgets must pass independently before applying either result to a 100 M chain.

#### Retention, load, controls, and reproduction

The initial registration and all original completed outputs are retained. All 13 checked production/recorder/harness/author-finder source hashes match the preceding controls, recorded in [t7a-source-controls.csv](t7a-source-controls.csv). Existing default-off and recording-on bit controls remain the applicable evolution controls; this analysis adds no evolution path. The analysis uses one thread per stream, at most two streams/native replays concurrently (≤4 threads total), with measured per-process peak RSS ≤1.37 GB; even a conservative sum remains below 8 GB. Each largest stream analysis completed in under 90 s. No simulation or long detached job was started.

Small consolidated CSVs are in this directory. Detailed per-level tables, composite six/eight-node ray caches and losslessly compressed selected operation frames remain under `/private/tmp/ems-t7a/`; their hashes are recorded in [t7a-archived-tables.csv](t7a-archived-tables.csv) and the manifest. Composite caches were losslessly compressed with every array verified identical, including NaNs; redundant derived per-level ray caches were removed after hashing. [t7a-cleanup.csv](t7a-cleanup.csv) records the original and retained hashes. Native-stencil/replay scratch was also removed after hashing. All original streams, every cited composite profile and selected operation frame, and all detailed CSVs remain available. The analysis directory is now about 1.0 GB, alongside about 4.7 GB of original completed run output.

One bounded operation per command, with installed local Python and single-thread BLAS:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7a-analyze.py snapshots clock 6
# Repeat snapshots per case/level, then baseline, combined, profiles, common_packet_peaks.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7a-analyze.py operations clock 6
# Repeat operations per case/level; rk_budget per case/level; restriction for coarse levels 4 and 5.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7a-analyze.py summaries
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/private/tmp/ems-t7a/mpl /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7a-analyze.py figures
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7a-check.py
```

Native widths and ray_wakes are separate bounded operations. The native replay build uses the already installed local serial Chombo and T7E's one-compiler build helper; no library or production source is modified. [COMMIT-MANIFEST-T7-A.txt](COMMIT-MANIFEST-T7-A.txt) records this analysis, and the shared T7 manifest is refreshed. No commit.
'''
p=HERE/'README.md';old=p.read_text();heading='### A — completed run analysis (controller resumption)';boundary='## T7-E —'
if heading in old:
 start=old.index(heading);end=old.index(boundary,start);old=old[:start]+old[end:]
pos=old.index(boundary);p.write_text(old[:pos]+text.strip()+'\n\n'+old[pos:])
print('README T7 A APPENDED')
