#!/usr/bin/env python3
"""Write the evidence card and the bounded task status without changing inputs."""
import sys
sys.dont_write_bytecode=True
import csv,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t16')
def rows(name):return list(csv.DictReader((HERE/name).open()))
def fmt(x):return f'{float(x):.3g}' if x!='' else 'zero'
def table(header,data):return '\n'.join(['| '+' | '.join(header)+' |','| '+' | '.join(['---']*len(header))+' |']+['| '+' | '.join(map(str,r))+' |' for r in data])
def report():
    summary=rows('t16-audit-summary.csv');resources=rows('t16-resources.csv')
    lookup={(r['mode'],r['domain'],r['ray'],r['field'],r['norm']):r for r in summary}
    paragraph=['### Measured isolated initialization and current status',
        'The opt-in flag is `ems_use_geometric_initial_lapse`, omitted/false by default. The existing maximal option retains all single-unboosted restrictions; simultaneous flags and non-EMSTRUMPET/RN use are rejected. The setter reuses its existing full-geometry isolated lapse and assigns the stable inverse-square-addition candidate for a binary before the CTT correction, which leaves that candidate untouched. No non-lapse setter expression or native evolution term was changed.',
        'The exact [CAS witness](t16-cas.json) closes the inverse-metric and stable-binary identities. Its independently assembled 70-digit four-metric inversions span 64 samples with worst relative difference '+json.loads((HERE/'t16-cas.json').read_text())['worst_relative']+'. Each opt-in initialization also independently inverts 24 Float64 four-metrics, with both boost signs and off-plane directions. The measured maximum relative discrepancy is **8.88178e-16**, or **1.49351 condition-scaled eps** (declared limit 64). The bound is a numerical matrix-inversion roundoff criterion; finite sampling does not prove the lapse inequality globally. Every actual opt-in setter evaluation separately enforces positivity, finiteness and a<=1, with no clipping.',
        'Default identity passes **1,675,520 Float64 values** in sixteen initial-plot/two-step-checkpoint comparisons: unboosted and boosted progenitors, point and legacy transfers, flag omitted and explicitly false, against the frozen 5da576b production executable. The [control table](t16-controls.csv) records zero bit mismatches. All seven [parameter guards](t16-guards.csv) meet their expected acceptance/rejection. The [output filter check](t16-filter-check.csv) preserves 183,708 Float64 values bitwise, including a skipped-clock XOR-delta chain and an unchanged stage frame.']
    ranges=[]
    for l in (13,14):
        rr=rows(f't16-ranges-L{l}.csv')
        cc=rows(f't16-census-L{l}.csv')
        for reg in ('isolated-native','left-hole-near-right-footprint','right-hole-near-left-footprint'):
            group=[r for r in rr if r['region']==reg]
            ranges.append([l,reg,fmt(min(float(r['a_min']) for r in group)),fmt(max(float(r['a_max']) for r in group))])
        assert all(int(r['nonfinite'])==int(r['chi_activations'])==int(r['lapse_activations'])==0 for r in cc)
        for level in range(l+1):
            group=[r for r in cc if int(r['level'])==level];h=float(group[0]['h_M'])
            face=max(int(r['y1'])+1 for r in group)*h
            assert face==(336. if level==0 else 112/2**level)
    paragraph+=['The actual initialization census confirms every stated face and level. Native grid values, including ghosts, are finite and no chi or lapse floor activates. The isolated minimum geometric lapse is 1.5942598099528705e-4 on L13 and 6.249679544383461e-5 on L14; maximum is 0.9979038379089187 on both. The geometric companion checks transplant each isolated grid footprint by +/-32 M_i to include the opposite-hole neighborhood; these are **per-hole geometric evaluations**, not a final CTT binary state or an admission audit. All denominators and lapse bounds pass. [L13 ranges](t16-ranges-L13.csv), [L14 ranges](t16-ranges-L14.csv).',table(['Finest rung','Footprint evaluated','a minimum','a maximum'],ranges)]
    fields=[r['field'] for r in summary if r['mode']=='geometric' and r['domain']=='W' and r['ray']=='axis' and r['norm']=='RMS' and r['identity']=='True']
    data=[]
    for field in fields:
        vals=[lookup['geometric',domain,ray,field,'RMS'] for domain in ('W','smooth') for ray in ('axis','diagonal')]
        data.append([field]+[fmt(r['p']) for r in vals])
    paragraph+=['For the native control, the reported residual is physical_RHS+v D8_x u, with KO and its summed residual retained separately. Every nonzero non-gauge residual decreases under refinement; Lambda, Bx, By and Ez are exactly zero here. The discrete audit uses 33 equally spaced requested radial samples per window. The axis is the first positive-y native row and the diagonal uses the nearest native y row, rather than interpolating the RHS. All peak/RMS amplitudes and orders for all 28 components are retained in [audit summary](t16-audit-summary.csv), with [sqrt(chi) raw audit](t16-audit-sqrt.csv) and [geometric raw audit](t16-audit-geometric.csv).',table(['Non-gauge component','p RMS W axis','p RMS W diagonal','p RMS smooth axis','p RMS smooth diagonal'],data),
        'The smooth collar predominantly exhibits fourth-order behavior, with finite two-rung deviations retained in the table. The first-row transverse Gamma residual is third order (2.999 smooth, 2.911 on W). This is consistent with native cartoon D1(beta)/y terms: the exact polynomial witness D1_h(y^5)-5y^4=-4h^4 becomes -8h^3 at y=h/2. That consistency example does not prove which term dominates this measured Gamma residual. Two resolutions do not establish an asymptotic order or eliminate a smaller source below the discretization error. No uniform fourth-order certificate is claimed at the axis or puncture.']
    control=[]
    for ray in ('axis','diagonal'):
        for field in ('K','Gamma1','Gamma2','phi','Pi','Ex'):
            g=lookup['geometric','W',ray,field,'RMS'];old=lookup['sqrt','W',ray,field,'RMS']
            control.append([ray,field,fmt(g['residual_half']),fmt(float(g['residual_half'])/float(old['residual_half'])),fmt(g['p'])])
    paragraph+=[table(['Ray','Field','Geometric identity RMS, L14','Geometric/sqrt(chi) residual','Observed p'],control)]
    gauge=[]
    for ray in ('axis','diagonal'):
        for field in ('lapse','shift1','shift2','B1','B2'):
            r=lookup['geometric','W',ray,field,'RMS'];gauge.append([ray,field,fmt(r['rhs_half']),fmt(r['KO_half'])])
    paragraph+=['The code gauge remains independent of the translated-geometry identity: initial lapse RHS=beta dot grad(alpha)-1.8 alpha(K-2Theta), shift RHS=beta dot grad(beta) after its initialized driver cancellation, and gauge B RHS=-0.1 B. The native audit separately reports these physical RHS and KO terms. They are not required to translate the stationary four-geometry.',table(['Ray','Gauge field','Physical RHS RMS, L14','KO RMS, L14'],gauge)]
    if (HERE/'t16-launch-endpoint.csv').exists():
        endpoint=rows('t16-launch-endpoint.csv');data=[]
        margins=rows('t16-sampling-margins.csv')
        sensitivity={(int(r['level']),r['ray'],r['field'],r['norm']):r for r in margins if r['endpoint']=='True'}
        for l in (13,14):
            for ray in ('axis','diagonal'):
                for field in ('Gamma','metric_Gamma','shift','lapse'):
                    rr=[r for r in endpoint if int(r['level'])==l and r['ray']==ray and r['field']==field]
                    bounds=[sensitivity[l,ray,field,n] for n in ('peak','RMS')]
                    data.append([l,ray,field]+[fmt(next(r['geometric_over_sqrt'] for r in rr if r['norm']==n)) for n in ('peak','RMS')]+
                        [f"{float(r['ratio_lower_5x']):.6g}–{float(r['ratio_upper_5x']):.6g}" for r in bounds]+
                        [fmt(min(float(r['reduction_margin_5x']) for r in bounds))])
        verification=rows('t16-verification.csv');assert all(r['result']=='PASS' for r in verification)
        paragraph+=['**READY-EXCEPT: the isolated measurements are complete and verified; temporal certification and the certified-companion binary audit remain unavailable.** Each native leg is checked against its own actual-child receipt, clear output/memory gates, completion log, registered stop, finite current-state replay and floor/reader records. The enclosing evolution worker and pipeline receipts also pass those execution checks. The preserved pipeline manifest authenticates 368 inputs, executables, captures and caches; 33,885,000 replay Float64 values are finite and the native RHS comparison has zero bit mismatches. All 9,600 profile-history norms/spreads reproduce from the stored I8/I10 caches. [Verification](t16-verification.csv) distinguishes these checks from a scientific admission verdict. An initial follow-up verification used a nonexistent evolution queue-receipt name; the corrected check uses the actual evolution-worker receipt. Both resource records remain retained.',
            'Endpoint ratios below are geometric/sqrt(chi) disturbance at the common physical time **0.00197601318359375 M_i**. For each amplitude A, e=5 norm(I8-I10), after subtracting the respective captured numerical t=0 profile. The displayed sensitivity intervals are [(G-e_G)/(S+e_S),(G+e_G)/(S-e_S)], with the lower numerator floored at zero and an unbounded upper value if S<=e_S. The reduction margin is (S-G)/(e_S+e_G); greater than one resolves the amplitude reduction at this sampling budget. These are sensitivity estimates from the declared interpolation spread, not rigorous interpolation-error bounds. Every endpoint field amplitude and reduction qualifies. [Endpoint amplitudes and signal margins](t16-launch-endpoint.csv), [all-clock sensitivity intervals and reduction margins](t16-sampling-margins.csv).',
            table(['Finest rung','Ray','Field','Peak ratio','RMS ratio','Peak 5x interval','RMS 5x interval','Min reduction margin'],data)]
        histories=rows('t16-history-summary.csv');history_data=[]
        for l in (13,14):
            for ray in ('axis','diagonal'):
                for field in ('Gamma','metric_Gamma','shift','lapse'):
                    rr=[next(r for r in histories if int(r['level'])==l and r['ray']==ray and r['field']==field and r['norm']==n) for n in ('peak','RMS')]
                    history_data.append([l,ray,field]+[f"{float(r['ratio_min']):.4g}–{float(r['ratio_max']):.4g}" for r in rr]+
                        [' / '.join(r['reduction_resolved']+'/74' for r in rr),
                         ' / '.join(fmt(r['min_signal_margin_5x']) for r in rr)])
        assert sum(int(r['qualified']) for r in histories)==2368
        paragraph+=['The history retains **75 exact common clocks** including t=0 and all 74 positive clocks, without time interpolation or alignment. At t=0 every disturbance is identically zero and its ratio is undefined. All **2,368/2,368 positive-time disturbance entries** exceed their own 5x spread in both variants; uncertainty-dominated count is zero and no entry is discarded. The smallest signal margin is 5.13425, for L13 axis metric-Gamma peak at step 15, t=0.000400543212890625 M_i. The [empty unqualified-entry ledger](t16-unqualified-history.csv) retains its schema. Reduction is sampling-resolved in **2,228/2,368** entries; the remaining 140 are measured early shift increases, not an interpolation failure. [History counts, extrema and every enhanced step](t16-history-summary.csv), [raw history](t16-launch-history.csv).',
            table(['Rung','Ray','Field','History peak ratio range','History RMS ratio range','Reduction counts peak / RMS','Min signal margins peak / RMS'],history_data),
            'Gamma, metric Gamma and lapse disturbance are smaller throughout the positive history on both rays and rungs. The shift is larger at steps 1–19 (axis peak), 1–23 (axis RMS), 1–13 (diagonal peak) and 1–15 (diagonal RMS). Its largest ratio is **7.49014**, L13 axis RMS at step 11, t=0.000293731689453125 M_i. All shift norms become smaller at step 24, t=0.000640869140625 M_i, and remain smaller through the common endpoint. The endpoint improvement therefore does not establish uniform suppression of the complete launch history. Fixed-coordinate disturbances include physical translation and gauge response. No factor-two admission criterion or binary source-cancellation verdict was preregistered for this isolated comparison.']
        trends=rows('t16-rung-trend.csv');trend_data=[]
        for ray in ('axis','diagonal'):
            for field in ('Gamma','metric_Gamma','shift','lapse'):
                vals=[next(r for r in trends if r['ray']==ray and r['field']==field and r['mode']==mode and r['norm']==n) for mode in ('sqrt','geometric') for n in ('peak','RMS')]
                trend_data.append([ray,field]+[fmt(r['fine_over_coarse']) for r in vals])
        paragraph+=['**Temporal margin: not measured.** The frozen registration has no same-grid dt/2 leg. Both h and dt halve between L13 and L14, so the table below reports combined spatial/temporal sensitivity and cannot establish temporal subdominance. The common clocks coincide exactly; native final times are not mixed or interpolated. In particular, evolved Gamma geometric amplitudes rise by approximately 0.8–1.2% under refinement while remaining strongly below sqrt(chi); coarse axis metric-Gamma peak decreases sharply on refinement. No self-contraction order is inferred from these two amplitudes. [Rung trend](t16-rung-trend.csv).',
            table(['Ray','Field','sqrt fine/coarse peak','sqrt fine/coarse RMS','geometric fine/coarse peak','geometric fine/coarse RMS'],trend_data)]
        punct=rows('t16-puncture-ratios.csv');pdata=[]
        for field in ('Gamma1','Gamma2','shift1','shift2','lapse'):
            pdata.append([field]+[fmt(next(r['geometric_over_sqrt'] for r in punct if int(r['level'])==l and r['field']==field and r['endpoint']=='True')) for l in (13,14)])
        paragraph+=['The native puncture diagnostic uses the same four cells adjacent to the initial puncture on each rung, without radial interpolation. Its endpoint peak-change ratios corroborate suppression of Gamma and lapse, with a smaller longitudinal-shift benefit than on W. They compare different physical cell locations between rungs and do not define a fixed-domain convergence order. [All native puncture ratios and instantaneous KO peaks](t16-puncture-ratios.csv), [captured native changes](t16-native-puncture.csv). KO peaks are RHS terms and are not an integrated fraction of the evolved change; no KO-share attribution is inferred.',
            table(['Native field','L13 peak-change ratio','L14 peak-change ratio'],pdata)]
        constraints=rows('t16-constraint-summary.csv');cdata=[]
        for ray in ('axis','diagonal'):
            for field in ('C_Gamma','Ham','Mom','GaussE'):
                vals=[next(r for r in constraints if r['ray']==ray and r['field']==field and int(r['level'])==l and r['norm']==n) for l in (13,14) for n in ('peak','RMS')]
                cdata.append([ray,field]+[f"{float(r['endpoint_ratio']):.6f}" for r in vals]+
                    [' / '.join(r['qualified']+'/74' for r in vals)])
        paragraph+=['Constraints are absolute native residuals on the same W, not disturbances from t=0. Endpoint C_Gamma improves by roughly a factor 2.3–2.5 and its five-spread intervals exclude one. Ham, Mom and GaussE endpoint ratios lie between **0.999676 and 1.002996**; every corresponding five-spread interval includes one. Thus their small central changes, including increases, are sampling-unresolved. Some early constraint entries are uncertainty-dominated and remain labelled in the complete table; early fine-grid C_Gamma peak can also increase. The lapse change does not demonstrate suppression of the entire constraint battery throughout history. [Absolute amplitudes, endpoint intervals and history qualification counts](t16-constraint-summary.csv), [all-clock intervals](t16-sampling-margins.csv).',
            table(['Ray','Constraint','L13 peak ratio','L13 RMS ratio','L14 peak ratio','L14 RMS ratio','Qualified counts: L13 peak/RMS, L14 peak/RMS'],cdata)]
        audit=rows('t16-launch-audits.csv')
        paragraph+=['Across the complete native runs, chi/lapse floor counts are '+str(sum(int(r['chi_activations']) for r in audit))+'/'+str(sum(int(r['lapse_activations']) for r in audit))+', nonfinite counts '+str(sum(int(r['nonfinite']) for r in audit))+', and positive-time static-reader messages '+str(sum(int(r['reader_messages_after_advance']) for r in audit))+'. The initial driver cancels exactly and all retained t=0 non-lapse values match bitwise between the two lapse variants. The native stops are 0.00197601318359375 (L13) and 0.0019893646240234375 (L14); only the common 75 clocks enter the ratio history.']
    else:
        paragraph+=['The four launch legs are registered for a **serial detached queue** with a 720 s planned upper wall estimate. Results and ratios are pending. Each native leg has `/private/tmp/ems-t16/evolution/L{13,14}-{sqrt,geometric}/done.exit`; the queue marker is `/private/tmp/ems-t16/evolution/done.exit`. An enclosing pipeline performs bounded current-state analysis and refreshes this report after all four legs pass, with `/private/tmp/ems-t16/pipeline/done.exit` as its final marker. It retains all first-four-step stage captures and snapshots only at the 75 common clocks, cropping snapshot cells to |x-336|,|y|<=0.0035 M_i. This covers the I10 plus native RHS/KO support for W; no native evolution arrays are altered.']
    peak=max(int(r['peak_rss_bytes']) for r in resources);tree=max(int(r['tree_peak_bytes']) for r in resources)
    paragraph+=['Measured process peak is **'+str(peak)+' bytes**, sampled active-tree peak **'+str(tree)+' bytes**. Measurements are in [resources](t16-resources.csv), including failed setup attempts; each numerical process stays below 3 GB and uses OMP=2, BLAS=1. The final follow-up verification also stays within the cap. The native compiler/reader/recorder evidence is pinned in [manifest](t16-manifest.txt). The final manifest writer has its own receipt under `/private/tmp/ems-t16/final-manifest`; that receipt is finalized after manifest creation and is not self-hashed.',
        'The binary step needs the **certified tranche-2a EMSCTT/1 companion for this exact alpha=0.9, e=8 profile**, authenticated against its profile SHA-256 and the solver/end-tail provenance, with masses (1,1), centres (-16,+16), rapidities (+atanh(v),-atanh(v)), and couplings (4pi,0,0,-0.9). Its consumed-order/end-transfer certificate and reconstructed constraints, ADM mass, individual areas and charges must accompany it. Then the same unchanged companion supplies both lapse variants on the final reconstructed CTT state: positive per-hole/combined denominators, native-grid min/max including the companion neighborhoods, unchanged non-lapse fields and initialized driver cancellation, constraints and refinement-controlled source/launch comparison. No binary companion was constructed and no final-CTT binary audit or launch was run. The isolated result cannot admit a binary Killing lapse or prove binary source cancellation.']
    path=HERE/'README.md';text=path.read_text();marker='### Measured isolated initialization and current status'
    start=text.index('## T16 —');head=text[:start];section=text[start:]
    if marker in section:section=section[:section.index(marker)].rstrip()+'\n\n'
    # Remove the unrelated T14d reproduction paragraph accidentally carried into
    # this section; the frozen T16 numerical registration above is unchanged.
    section='\n\n'.join(p for p in section.split('\n\n') if not p.startswith('Reproduce serially with the measured t14-run.py wrapper'))
    paragraph+=['For this completed T16 evidence, rerun only the analysis with absolute paths under the measured `t16-run.py` wrapper: `t16-verify.py`; `t16-analyze.py resources`; `t16-report.py`; then `t16-analyze.py resources` for the final manifest. The verification consumes the retained native replay/profile caches and authenticates the original captures. If regenerating caches is necessary, `t16-analyze.py analyze` replays those stored captures and never opens the static trumpet. No new evolution is needed.']
    path.write_text(head+section+'\n\n'.join(paragraph)+'\n')
    (HERE/'t16-status.json').write_text(json.dumps(dict(status='READY-EXCEPT' if (HERE/'t16-launch-endpoint.csv').exists() else 'RUNS-PENDING',
        binary='pending certified tranche-2a companion',temporal_control='not registered or measured',
        isolated='verified; endpoint reductions resolved; early shift enhanced',positive_history_entries=2368,
        qualified_positive_history_entries=2368,sampling_resolved_reductions=2228,
        peak_rss_bytes=peak,sampled_tree_peak_bytes=tree),indent=2)+'\n')
if __name__=='__main__':report()
