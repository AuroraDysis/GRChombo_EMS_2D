from pathlib import Path
import csv,math,json
H=Path('Tests/EMSNative')
def rows(n):return list(csv.DictReader((H/('t7e-'+n+'.csv')).open()))
def f(x):return f'{float(x):.6e}'
def table(headers,data):
 headers=[x.replace('|',r'\|') for x in headers]
 return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in data])+'\n'
q=rows('qualified-horizons');drift=rows('horizon-budgets');orders=rows('constraint-orders');local=rows('checkpoint-localization');control=rows('coordinate-controls');layout=rows('layout');coord=rows('coordinate-errors');rich=rows('horizon-richardson');perf=rows('analysis-performance');runs=rows('finder-runs')
out=[]
out.append('''## T7-E — exp-0019 E, three-grid qualification and constraint diagnosis

**READY-EXCEPT:** every supplied E checkpoint has a qualified 48/96-point numerical horizon. High passes both budgets on those checkpoints; the early interval between t=0 and the first late checkpoint is not independently qualified. The t=0 Hamiltonian loss of order is dominated on the finer grids by absolute-coordinate rounding in the initial setter, amplified by second derivatives. The evolved far-mask signal is a time-dependent discretization packet/wake above the instantaneous Float64 floor; its exact source operation is not established by the saved E states alone.

All EMS and cluster inputs are read-only. This tranche adds standalone diagnostics and Python analysis; it changes no production evolution or finder C++ source, equations, gauge, floors, transfer operator, parameters, or checkpoint. No commit is made. Static-file evaluation is guarded by checkpoint time exactly zero. Evolved diagnostics use the saved current fields and other resolutions only. T7 A's space process, directory, parameters, and marker were not accessed. No new long simulation was started.

### Qualified horizons and numerical baselines

The four checkpoints supplied for each grid are used. Low/mid reuse the eight qualified T6 checkpoint rows; high adds eight independent frozen solves, one at each N_theta=48 and 96 for each checkpoint. T6's frozen-state setup, native Lagrange-4 field sampling, parity treatment, chase multiplier 2, 3000-update ceiling, and mean-square expansion thresholds 1e-7, 1e-10, 1e-12 are unchanged. The final condition is FOUND at stage 2 and sqrt(err)<=1e-6, independently at both angular resolutions. A finder label such as close or found from the in-run history is not substituted for this qualification. `err` is the surface-area-weighted mean square of outgoing expansion, whereas `<Theta+>` is its area-weighted mean. Numbers below are N_theta=96; the CSV also retains N_theta=48, angular and stopping sensitivities, and the native expansion RMS.

The current executable includes disabled T7 diagnostic branches; its executable hash differs from the earlier T6 builds. Both author's finder headers have identical hashes. A bounded control replays E-low t=0 at N_theta=48: all four surface rows and all 15 non-wall-clock columns exactly match the cached T6 result, including intermediate stages, area, charge, residual, updates and status. See [t7e-controls.csv](t7e-controls.csv). Thus the qualification is the same numerical frozen harness; executable byte identity is not claimed. [t7e-finder-runs.csv](t7e-finder-runs.csv) and [t7e-finder-input-audit.csv](t7e-finder-input-audit.csv) retain the eight high solves and hashes.

There is no in-run t=0 finder row. Each t=0 qualification uses the earliest saved numerical shape as a seed but solves on t=0 checkpoint fields, supplying the numerical baseline. High's 9.43055556 M checkpoint has no same-time saved shape; it uses the preceding 9.33333333 M shape. The independent root solve still satisfies the same criterion. The native values beside this row are linearly interpolated between the native 9.33333333 and 9.52777778 M rows. All other late native values are exact printed-time rows.
''')
data=[]
for r in q:data.append([r['scale'][2:],f"{float(r['time_M']):.8f}",f"{float(r['A96']):.12f}",f"{float(r['Q96']):.12f}",f(r['theta_plus_rms96']), '—' if r['inrun_row']=='NO_T0_ROW' else f"{float(r['inrun_A']):.9f}",'—' if r['inrun_row']=='NO_T0_ROW' else f"{float(r['inrun_Q']):.9f}"])
out.append(table(['grid','t/M','qualified A/M²','qualified Q/M','expansion RMS × M','in-run A/M²','in-run Q/M'],data))
out.append('Full values and sensitivities: [t7e-qualified-horizons.csv](t7e-qualified-horizons.csv). Relative curves use each grid\'s own qualified numerical t=0 value, with no static-profile target.\n\n')
data=[]
for r in drift:data.append([r['scale'],f"{float(r['terminal_time_M']):.8f}",f(r['terminal_relative_A']),f(r['terminal_relative_Q']),r['checkpoint_area_pass'],r['checkpoint_charge_pass']])
out.append(table(['grid','last t/M','ΔA/A(0)','ΔQ/Q(0)','|ΔA/A|≤1e-3 on checkpoints','|ΔQ/Q|≤2e-4 on checkpoints'],data))
out.append('''Absolute area and charge drifts both decrease with resolution. High's maximum qualified-sample drifts are 4.634627e-7 and 4.676447e-5. The 48/96 difference in **relative drift** is retained separately from absolute angular quadrature bias; root and angular sensitivities do not threaten these checkpoint budget decisions. High's area drift is smaller than the absolute stopping sensitivity, so its tiny sign should not be assigned physical significance. Low/mid charge drifts survive independent qualification and exceed the charge budget. The high charge drift changes sign.

The dense native high history has maximum |ΔA/A(0)|=4.478462e-3 and |ΔQ/Q(0)|=4.675473e-5. Its area excursion exceeds the area budget but occurs on unqualified surfaces. These rows neither prove a true horizon budget violation nor certify the entire interval. The available checkpoints leave 0<t<4.958333 M unqualified for high; connecting qualified points in the figure is not a proof between them.

![Qualified E horizon drifts; pale curves are unqualified native finder values](figures/t7e-horizons.png)

For unequal checkpoint times, [t7e-horizon-richardson.csv](t7e-horizon-richardson.csv) gives conditional linear interpolation of qualified A and Q to common 0, 4.5, 9 and 10 M. An alternative adds the nearby native-history trend to the qualified checkpoint value and records its discrepancy. Native shape lag makes area time correction particularly unreliable. Richardson uses p=log(|low-mid|/|mid-high|)/log(1.5) and extrapolates only for monotone triples with p>0. These are observable-difference fits, not constraint-RMS orders or a registered asymptotic proof.
''')
data=[]
for r in rich:data.append([r['time_M'],r['observable'],f"{float(r['p']):.6f}",'non-monotone' if r['monotone']=='False' else (f(r['extrapolated']) if math.isfinite(float(r['extrapolated'])) else 'p≤0'),f"{float(r['relative_p']):.6f}" if math.isfinite(float(r['relative_p'])) else '—',f(r['relative_extrapolated']) if math.isfinite(float(r['relative_extrapolated'])) else '—',f(r['time_matching_sensitivity_max'])])
out.append(table(['common t/M','observable','apparent p','extrapolated absolute value','relative-drift p','extrapolated relative drift','max absolute time-matching sensitivity'],data))
out.append('At 10 M, the conditional charge fit gives p=4.903729 and ΔQ/Q(0)=2.015876e-4 in the extrapolation, marginally above the charge budget; high\'s own measured drift is below it. This extrapolation should not be treated as a certified continuum limit. Late area triples are non-monotone and have no Richardson extrapolation.\n\n### t=0 Hamiltonian floor: operation localization\n\n')
out.append('''[t7e-coordinate-errors.csv](t7e-coordinate-errors.csv) enumerates **all** uncovered cells of the four masks and compares the ordinary global-coordinate expression `(i+0.5)*h-224` with correctly rounded fused evaluation. Low's dyadic spacings give exactly zero difference. Mid/high give coordinate RMS around 9.1e-15/8.2e-15 M and maxima 1.421085e-14/1.304512e-14 M. The error is in initial point coordinates; it is not a stored continuum Hamiltonian residual, an evolution damping target, a coordinate stretch, or a change to the evolution mesh.

The analytic diagnostic loads the checksum-identical exp-0019 E.trumpet through the actual C++ EMSTRUMPET reader, using its inversion, Clenshaw evaluation, and derivative-coefficient construction. It differentiates the reconstructed radial metric analytically and evaluates R-6k²-16πrho with the code's EMS normalization. It does not finite-difference a radial table, use the ODE to force the residual to zero, or substitute stored static targets for evolved data. The curvature identity is verified independently in [CAS evidence](../../scripts/cas/t7e-evidence.md); algebraic identity — production physics not certified. The profile SHA-256 is `2a8de074ae17c0b11d323d4b0933a6bdb7430a5305473c8cc7d4ce37f39fa793`, equal to the cluster input manifest.

The first three RMS columns below are the exact published cylindrical reductions. Continuum and Float64 columns are a deterministic checkpoint stencil census: stride 4 in each direction in bulk, every axis/interface/near-finer strip cell retained, with multiplicity weighting. They are sampled estimates, not replacements for the raw reductions. The entire far mask is evaluated without subsampling. Matching current same-level/parity ghosts and six-point coarse interpolation are reconstructed; covered coarse cells are removed. T4b's native arithmetic sensitivity audit is reused alongside a separate call to the actual constraint kernels. The continuum residual lies at its own Float64 term-cancellation scale, around 1e-15, many orders below the finite-difference residual.
''')
data=[]
for mask in ('inside_inner_ring','cavity','between_rings','cavity_core'):
 r=next(x for x in orders if float(x['time_M'])==0 and x['mask']==mask and x['constraint']=='Ham');v=[next(x for x in local if x['scale']==s and float(x['time_M'])==0 and x['mask']==mask and x['subset']=='all') for s in ('low','mid','high')]
 data.append([mask,' / '.join(f(r[s+'_rms']) for s in ('low','mid','high')),f"{float(r['order_low_mid']):.3f}/{float(r['order_mid_high']):.3f}",' / '.join(f(x['continuum_Ham_rms']) for x in v),' / '.join(f(x['eps_chi_h2']) for x in v),' / '.join(f(x['coordinate_Ham_sensitivity']) for x in v)])
out.append(table(['mask','raw Ham low / mid / high','p low→mid / mid→high','analytic Ham RMS low / mid / high','ε|chi|/h² low / mid / high','coordinate-sensitive Ham scale low / mid / high'],data))
out.append('''The last column is the conservative first-order scale (32/3) ε·224·|chi_r x/r|/h², using actual analytic jets and local spacing; 32/3 is twice the absolute fourth-order second-derivative stencil weight sum. It is a sensitivity scale, not an interval bound. The ε|chi|/h² column alone misses the loss of absolute precision from global centring. [t7e-checkpoint-localization.csv](t7e-checkpoint-localization.csv) retains spacings, level membership, counts, native-operation budget, setter replay differences and metric/Gamma terms. Setter-versus-checkpoint chi differences on ordinary coordinates are at Float64 precision. Metric Ricci at t=0 is far too small to account for the 1e-7 residual.

The decisive paired control evaluates **only chi** with fused initial coordinates, keeping all other t=0 checkpoint fields and the native constraint stencil unchanged. No current evolution field is modified. On the same bulk cells within each grid:
''')
data=[]
for mask in ('inside_inner_ring','cavity','between_rings','cavity_core'):
 rr=[next(x for x in control if x['scale']==s and x['mask']==mask) for s in ('low','mid','high')]
 data.append([mask,' / '.join(x['bulk_cells'] for x in rr),' / '.join(f(x['native_Ham']) for x in rr),' / '.join(f(x['chi_fused_Ham']) for x in rr)])
out.append(table(['mask','paired bulk cells low / mid / high','original Ham RMS','fused-coordinate chi Ham RMS'],data))
out.append('''High's reductions are factors 10.91, 15.51, 14.67 and 15.59 respectively. Low scarcely changes and remains dominated by genuine fourth-order truncation. This directly localizes the dominant finer-grid floor to coordinate evaluation feeding the initial chi values and then second-derivative stencils. Residual field-value/reader rounding and finite-resolution truncation remain; the sparse paired control is not a new all-mask convergence gate. The effect occurs in bulk, so it is not tied to a single patch or seam. The largest remaining analytic-data residual is negligible on these masks. The raw CSV additionally flags t=0 outer_ring Ham at 3.75/2.65; the blanket statement that every other mask passes is therefore not used.

### Evolved norms, matched times, and far-mask attribution

The published reduction is **legacy cylindrical RMS**: sqrt(sum y C²/sum y) on uncovered cells, with no h_level² factor. It is not T6's physical AMR-volume norm. The aggregate CSV has no level breakdown, so its full time history cannot be converted to a volume norm. The checkpoint replay supplies both legacy and correct volume norms separately; its far-mask evaluation includes every cell. Claims about the full time history below are explicitly conditional on the published legacy norm.

[t7e-constraint-rms.csv](t7e-constraint-rms.csv) retains all native-time RMS and maxima. [t7e-constraint-orders.csv](t7e-constraint-orders.csv) gives **all ten E masks**, all three constraints, and 24 requested times: every low/mid output through 9.625 M plus 10 M. Low/mid coincide to roundoff at their 0.4375 M cadence; high has a 0.486111 M cadence plus its special checkpoint output. All three coincide exactly at 0, 4.375 and 8.75 M. Other rows linearly interpolate RMS within each grid; bounds, fractions and endpoint counts are recorded. No extrapolation past a run endpoint is made. Alternative nearest-time and squared-RMS interpolation orders are also retained: transient rows sensitive to time pairing are not order qualifications. Orders are log(RMS_coarse/RMS_fine)/log(1.5), not self-convergence of evolved fields. High near-horizon apparent orders up to 14 indicate strong pre-asymptotic differences and should not be called established fourteenth-order convergence.

In the tables, `!` means at least one pair is below 3, and `N` means at least one pair is negative. The CSV retains flags at every output time. Times 0, 4.375 and 8.75 below are exact; 10 M is interpolated.
''')
for t in (0.,4.375,8.75,10.):
 out.append(f'**Orders at t={t:g} M**\n\n');data=[]
 for mask in ('horizon','inside_inner_ring','cavity','cavity_core','between_rings','outer_ring','outer_ring_core','far','far_core','outer_boundary_shell'):
  vals=[]
  for field in ('Ham','Mom','GaussE'):
   r=next(x for x in orders if float(x['time_M'])==t and x['mask']==mask and x['constraint']==field);vals.append(f"{float(r['order_low_mid']):.3f}/{float(r['order_mid_high']):.3f}"+(' N' if r['flag']=='NEGATIVE' else (' !' if r['flag']=='BELOW_3' else '')))
  data.append([mask,*vals])
 out.append(table(['mask','Ham low→mid / mid→high','Mom low→mid / mid→high','GaussE low→mid / mid→high'],data)+'\n')
out.append('Exact common-time RMS triples are retained below for the questioned far rows.\n\n');data=[]
for t in (4.375,8.75):
 for field in ('Ham','Mom','GaussE'):
  r=next(x for x in orders if float(x['time_M'])==t and x['mask']=='far' and x['constraint']==field);data.append([t,field,*[f(r[s+'_rms']) for s in ('low','mid','high')]])
out.append(table(['t/M','constraint','low RMS','mid RMS','high RMS'],data))
out.append('''![All E mask constraint histories at their native times](figures/t7e-constraint-rms.png)

**Measured:** far Ham stays small and convergent through roughly 3.4 M, then changes sharply at 3.888889/3.9375 M on high/mid. High grows from 3.808841e-11 at 1.944444 M to 4.146612e-9 at 5.833333 M and falls to 3.053811e-10 at 10.013889 M. Mid's large pulse reaches 2.522357e-9 at 5.6875 M. A constant continuum-data or instantaneous derivative-roundoff floor does not describe that history. Far GaussE continues to converge near fourth order.

Far cells lie on uncovered levels 3 and 4. The actual level-4 square faces are ±5.541666667 M and y=5.541666667 M on mid, ±5.444444444 M and y=5.444444444 M on high. Their corners are at 7.837100158 and 7.699607 M. The inward parent/child interface (level-5 face) is at 2.770833333 M on mid and 2.722222222 M on high, with corners near 3.918550079/3.849803 M immediately inside the far mask. [t7e-layout.csv](t7e-layout.csv) records every actual union and same-level seam. These are square faces, not spherical radii; a radial shell can contain both levels. The geometry changes by block alignment between grids, so sharp transient norms also compare different physical face positions.

The full-cell checkpoint replay gives the following far subdivisions. An interface is the outer two-cell layer of the fine union; near-finer is the two coarse-cell layer outside the covered region; bulk excludes both and the two-cell cartoon-axis strip. No covered coarse cell enters the norm. Native-versus-sensitivity replay differences are much smaller than the observed far residual.
''')
data=[]
for scale in ('mid','high'):
 for t in sorted(set(float(x['time_M']) for x in local if x['scale']==scale and float(x['time_M'])>0)):
  v=[next(x for x in local if x['scale']==scale and float(x['time_M'])==t and x['mask']=='far' and x['subset']==ss) for ss in ('all','bulk','interface','near_finer')]
  data.append([scale,f'{t:.8f}',*[f(x['Ham_rms']) for x in v],f(v[0]['Ham_roundoff']),f(v[0]['volume_Ham_rms'])])
out.append(table(['grid','t/M','all Ham','bulk Ham','fine-face Ham','receiving-coarse Ham','instant Float64 Ham scale','volume Ham'],data))
out.append('''At high 4.958333 M, bulk Ham is about 14 times the fine-face value. At 9.430556 M, the level-4 radial band around 6.44–6.56 M has Ham about 6.6e-10 and Mom about 5.7–6.2e-10, while a separate inner band near 4.06 M has Ham about 6.5e-10. Near 10 M that inner band persists at about 8.1e-10. These features are present in bulk cells, not solely in first ghost-reading cells. This is a real error in the current evolved state at these spacings, rather than just roundoff in printing the diagnostic. The instantaneous conservative Ham sensitivity is about 1.8e-11 on mid and 3.9e-11 on high; direct ε|chi|/h² is smaller still. Accumulated evolution roundoff is not bounded by this instantaneous estimate, so it is not excluded solely by the ratio.

A separate frozen nine-point-support control recomputes Gamma from the **current** metric before evaluating the same Ham kernel. It does not uniformly remove the residual and can increase it: at high 4.958333 M on sampled level-4 cells Ham changes from 5.114030e-10 to 3.290200e-9. At high 9.430556 M it changes from 3.401089e-10 to 2.855359e-10. This rules out assigning all of the observed signal to a removable stored-Gamma/instantaneous-evaluation floor. The control is recorded in [t7e-gamma-controls.csv](t7e-gamma-controls.csv); it is an algebraic frozen-state comparison and is never an evolution projection.

**Inferred:** the delayed onset, pulse-shaped growth/decay, and bulk bands are consistent with an outgoing constraint error and wake of the kind measured in T5/T6. The E data do not independently establish which E interface or which transfer/RHS operation emitted it. In particular, an error now in level-4 bulk may have been emitted at the inner level-5 face; its location at a later time does not identify its source. The four saved E times cannot isolate stage-time filling, restriction or KO. A claim that the E floor is proven to be the T6 face source is therefore withheld. The outer domain boundary is about 224 M from the puncture, far outside the 4–8 M mask; an outer-boundary packet cannot traverse that separation by 10 M at the light/gauge speeds in this run. The separately nonconvergent outer_boundary_shell is reported and flagged, not confused with the far wake.

![Far constraint profiles split by current AMR level](figures/t7e-far-profiles.png)

### Artifacts, checks, and load

[t7e-analysis-performance.csv](t7e-analysis-performance.csv) records the checkpoint replay. Peak parent RSS was 1,207,238,656 bytes, with one analysis thread; the largest high finder peak was 974,831,616 bytes per process. At most two two-thread finds ran concurrently. Builds used one compiler process; Python/BLAS and replay used one thread. The analysis stayed below four threads and 6 GB RAM. Each frozen find was capped at 115 solver seconds/120 process seconds; the slowest completed in 108.220 s. Large stencil scratch files were deleted immediately after use; final CSVs and figures are small. No simulation was evolved or monitored, and T7 A remains under the controller's existing run ledger.

Reproduction (installed local environment; every command is bounded analysis):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7e-localize.py build
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7e-localize.py layout
# For each supplied checkpoint: current-state replay; static evaluation occurs only if its time is exactly zero.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7e-localize.py high EMS_000103.2d.hdf5
# Each high checkpoint has separate 48/96 directories. Use a fresh OUT root for a deliberate repeat:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7e-horizons.py EMS_000103.2d.hdf5 96
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/private/tmp/ems-t7e/mpl /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7e-analyze.py
OPENBLAS_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python scripts/cas/t7e-initial-verify.py
OPENBLAS_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7e-report.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t7e-check.py
```

The qualification runner refuses an existing output directory; preserve cached results or set its OUT constant to a fresh output root before a deliberate repeat. Per-checkpoint scratch tables live in /private/tmp/ems-t7e/tables; consolidated deliverables are in this directory. [COMMIT-MANIFEST-T7-E.txt](COMMIT-MANIFEST-T7-E.txt) freezes this tranche, and [COMMIT-MANIFEST-T7.txt](COMMIT-MANIFEST-T7.txt) refreshes shared artifact hashes while retaining the old detached-run records without opening the running space directory.
''')
text=(H/'README.md').read_text();start=text.find('\n## T7-E')
section='\n'+ '\n'.join(out)
if start<0:text+=section
else:
 end=text.find('\n## ',start+4);text=text[:start]+section+(text[end:] if end>=0 else '')
(H/'README.md').write_text(text)
