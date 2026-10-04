#!/usr/bin/env python3
"""Verify completed post-flow maps and write the final T25 continuation.

No native probes, field reads, searches or advances. Run after t25-early-audit.py
and t25-early-postmap.py. The original seed/finder table remains separate from
the refined mean-root trials; a capped contour never becomes a measured MOTS.
"""
import importlib.util,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('early',HERE/'t25-early.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)

def main():
    audit=json.loads((HERE/'t25-early-audit.json').read_text())
    assert audit['all_native_receipts_verified'] and audit['physical_probes']==10
    assert audit['all_three_stage_passes']==0
    assert int((e.ROOT/'postmap-job/done.exit').read_text())==0
    assert 'T25_POSTMAP_COMPLETE; no searches, no advances' in (e.ROOT/'postmap-job/run.log').read_text()
    post=json.loads((HERE/'t25-early-postmap.json').read_text())
    assert post['status']=='COMPLETE' and post['new_finder_probes']==0 and len(post['rows'])==10
    bykey={(float(x['time']),int(x['hole'])):x for x in post['rows']}
    original={(float(x['time']),int(x['hole'])):x for x in audit['rows']}
    rows=[];brackets=[];spreads=[];peaks=[];native_maps=0
    for p in sorted(e.ROOT.glob('step*/post-flow-*/N*/early-map.resources.json')):
        receipt=e.verified(p.parent,'early-map',0,'T25_EXPANSION_MAP_COMPLETE; no advances')
        peaks.append(receipt['peak_rss_bytes']);native_maps+=1
    assert native_maps==20
    for s in e.samples():
        d=e.ROOT/s['tag']
        for hole in range(2):
            old=original[s['time'],hole];finder=old['finder'];x=bykey.get((s['time'],hole))
            out=dict(time=s['time'],hole=hole,configuration=s['configuration'],checkpoint=str(s['cp']),
                bracket_status=x['status'] if x else 'NO_RESOLVED_BRACKET',root_kind='MEAN_LINEAR_TRIAL' if x else 'NONE',
                pointwise_width_fraction=x['fractional_width'] if x and x['status']=='POINTWISE_NEGATIVE_TO_POSITIVE' else '',
                mean_pair_width_fraction='',r_cyl='',r_axial='',cells='',r_cyl_finest_cells='',r_axial_finest_cells='',
                local_dx_min='',local_dx_max='',A_trial='',Q_trial='',lapse_mean='',chi_mean='',K_mean='',K_minus_2Theta_mean='',
                theta_min='',theta_max='',theta_mean='',theta_rms='',centre='',
                finder_status=finder['status'] if finder else 'NOT_RUN_NO_BRACKET',
                finder_squared=finder['squared'] if finder else '',finder_A_trial=finder['A'] if finder else '',
                finder_Q_trial=finder['Q'] if finder else '',finder_seconds=finder['seconds'] if finder else '',
                finder_seed='PRELIMINARY_MAP' if finder else 'NONE',all_three_stages_qualified=0)
            if x:
                br=x['bracket'];pair=x['mean_pair'];r=x['root'];r96=x['root_N96']
                assert r['status']==r96['status']=='RESOLVED' and float(r['min_extent_over_dx'])>=3
                assert x['finder_seed_was_preliminary'] and x['new_finder_probes']==0
                mapped96={z['id']:z for z in e.read(d/'post-flow-centres/N96/map.csv')}
                inner,outer=mapped96[br['inner']],mapped96[br['outer']]
                if br['status']=='POINTWISE_NEGATIVE_TO_POSITIVE':
                    assert float(inner['theta_max'])<0<float(outer['theta_min'])
                    assert br['theta_inner_max']<0<br['theta_outer_min']
                else:
                    assert br['status']=='AVERAGE_ONLY_WEAKER'
                    assert float(inner['theta_mean'])<0<float(outer['theta_mean'])
                assert pair['theta_inner_mean']<0<pair['theta_outer_mean']
                assert float(pair['a_inner'])<=float(r['a'])<=float(pair['a_outer'])
                assert float(pair['c_inner'])<=float(r['c'])<=float(pair['c_outer'])
                for k in ('A','Q','lapse_mean','chi_mean','K_mean','K_minus_2Theta_mean','theta_mean','theta_rms'):
                    assert math.isfinite(float(r[k])) and math.isfinite(float(r96[k]))
                    spreads.append(dict(time=s['time'],hole=hole,quantity=k,N96=r96[k],N192=r[k],
                        absolute_spread=abs(float(r[k])-float(r96[k]))))
                out.update(r_cyl=r['a'],r_axial=r['c'],cells=r['min_extent_over_dx'],
                    r_cyl_finest_cells=float(r['a'])/e.H,r_axial_finest_cells=float(r['c'])/e.H,
                    local_dx_min=r['dx_min'],local_dx_max=r['dx_max'],A_trial=r['A'],Q_trial=r['Q'],
                    **{k:r[k] for k in ('lapse_mean','chi_mean','K_mean','K_minus_2Theta_mean','theta_min','theta_max','theta_mean','theta_rms','centre')},
                    mean_pair_width_fraction=float(pair['a_outer'])/float(pair['a_inner'])-1)
                brackets.append(dict(hole=hole,**br,mean_pair_inner=pair['inner'],mean_pair_outer=pair['outer'],
                    mean_pair_a_inner=pair['a_inner'],mean_pair_a_outer=pair['a_outer'],
                    mean_pair_c_inner=pair['c_inner'],mean_pair_c_outer=pair['c_outer'],
                    endpoint_N96_signs_agree=1,theta_inner_max_N96=inner['theta_max'],theta_outer_min_N96=outer['theta_min'],
                    theta_inner_mean_N96=inner['theta_mean'],theta_outer_mean_N96=outer['theta_mean']))
            else:assert old['root'] is None and finder is None
            rows.append(out)
    e.write(HERE/'t25-early-refined-time-table.csv',rows)
    e.write(HERE/'t25-early-final-brackets.csv',brackets)
    e.write(HERE/'t25-early-refined-angular-spreads.csv',spreads)
    midpoint=e.read(HERE/'t25-early-midpoint-t28.csv')
    assert len(midpoint)==20 and not e.strong(midpoint) and not e.b.brackets(midpoint)
    data=dict(status='PASS',postmap_native_receipts_verified=native_maps,endpoint_N96_sign_agreements=len(brackets),
        no_additional_finder_probes=True,no_evolution=True,physical_searches=10,strict_finder_passes=0,
        observed_postmap_peak_rss_bytes=max(peaks),rows=rows,
        midpoint=dict(surfaces=20,pointwise_pairs=0,average_pairs=0,
            minimum_theta=min(float(r['theta_min']) for r in midpoint),
            minimum_mean_theta=min(float(r['theta_mean']) for r in midpoint)))
    e.atomic_json(HERE/'t25-early-final-audit.json',data)
    report(rows,post,spreads,audit,data)
    print('T25_EARLY_FINAL_PACKET_VERIFIED',len(rows),native_maps,len(brackets),'strict passes 0')

def fmt(v):return f'{float(v):.7g}' if v!='' else '—'

def report(rows,post,spreads,audit,final):
    text=['## Time-resolved exp-0024 map, t = 0 to 28 M_i','',
        '**READY-EXCEPT: the bounded analysis is complete; no binary horizon search qualifies.**','',
        "The matching time-zero state uses exp-0024's d=16, rapidity 0.05778205303580913, geometric initial lapse and compiled ExperimentalGauge. "
        'Its companion SHA-256 is `107370ae00d8b69dc3122dc023e34bbff23b90afc8401bb233743c1182e6083c`; both initial inputs match the sealed submission. '
        'The d=32 T17 control above is a different initialization and is not substituted into this history. '
        'All eight collected checkpoints match their cluster SHA-256. Each uses the contemporaneous centres from `punctures.dat`. '
        'Static and CTT paths are absent from every positive-time map and search.','',
        'The native finest spacing is h=0.00042724609375 M_i. Spheres span 3h to 0.05 M_i with 4% initial radial spacing. '
        'Sign boundaries are bisected up to four times to at most 0.5% local spacing. '
        'If no spherical pointwise barrier exists, or it remains wider than 3%, the fixed fallback uses offsets −4, −2, 0, 2, 4 h '
        'and c/a=0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 2, with 8% initial radial spacing and the same refinement. '
        'Extents under three local cells are ineligible. N96/N192 and actual AMR ownership are retained. '
        'The selected endpoint signs agree at both angular resolutions. “Pointwise” here means every sampled angular point; it is not an interval proof between those points.','',
        'Before each search, the narrowest available pointwise barrier was selected, with endpoint RMS as the tie-breaker; '
        'a negative-to-positive mean pair inside it supplies the linearly interpolated seed. If no pointwise pair exists, '
        'one weaker mean pair may seed through `map_find_allow_average=true` (Tests-only, default false). '
        'There are eight pointwise-seeded and two mean-seeded searches. No bracket at t≥17.5 means no search there. '
        'All ten use the unchanged author finder: N48, chase=quota=1, 235 s cap, update cap 100000, inert floor window 100001, '
        'stages 1e-7/1e-10/1e-12. Every search stops at TIME_CAP in stage 1. The native squared residual is area-weighted mean θ₊², not its square root.','',
        'After those capped probes, additional maps use only their numerical final centres and a spheroid aspect ratio estimated from their numerical contours. '
        'These add no finder searches. They improve the t=0 pointwise width to 1.0025%, but cannot meet the requested few-percent pointwise width '
        'at positive time: best post-flow widths are about 38.1%, 30.5% and 22.1–22.7% at t=3.5, 7 and 10.5. '
        'The adjacent mean-root interval is narrower (about 0.5%), which does not make its angular expansion vanish. '
        'At t=14 only a weaker mean crossing exists. The original seed and receipts remain separate from these post-flow maps.','',
        'The table below reports the refined **mean-root trial** semiaxes a/c and physical surface integrals. '
        'Even inside a pointwise barrier, a mean root is not a solved marginally outer trapped surface (MOTS). '
        'Its finder residual belongs to the single earlier bracket-seeded probe, not to a new search at the refined centre. '
        'Hole 0 is the lower-x hole; hole 1 the upper-x hole. All entries are numerical-rung/checkpoint quantities.','',
        '| t / M_i | hole | bracket | a / c (M_i) | min extent / h | A trial | Q trial | mean lapse | finder θ² |',
        '|---:|---:|:---|:---|---:|---:|---:|---:|---:|']
    labels={'POINTWISE_NEGATIVE_TO_POSITIVE':'pointwise','AVERAGE_ONLY_WEAKER':'mean only','NO_RESOLVED_BRACKET':'none sampled'}
    for r in rows:
        text.append('| '+' | '.join([fmt(r['time']),str(r['hole']),labels[r['bracket_status']],fmt(r['r_cyl'])+' / '+fmt(r['r_axial']),
            fmt(r['cells']),fmt(r['A_trial']),fmt(r['Q_trial']),fmt(r['lapse_mean']),fmt(r['finder_squared'])])+' |')
    text+=['','The full 18-row table, including χ, K, K−2Θ, signed expansion range/RMS, local spacing and native finder A/Q, is '
        '[t25-early-refined-time-table.csv](t25-early-refined-time-table.csv). '
        '[t25-early-time-table.csv](t25-early-time-table.csv) retains the preliminary roots and actual finder seed centres. '
        '[t25-early-final-brackets.csv](t25-early-final-brackets.csv) gives both semiaxes of each refined barrier and the adjacent mean pair, '
        'signed endpoint bounds and the N96 confirmation. Per-checkpoint `t25-early-step*.csv` packets retain the spherical/selected-family map; '
        '`*-post-centres.csv`, `*-post-brackets.csv`, `*-post-barriers.csv` and `*-post-roots.csv` retain the additional numerical-centre maps. '
        'All remaining families and angular rows stay under `/private/tmp/ems-t25/early/step*/`, indexed by the manifest.','',
        '| t / M_i | χ mean (hole 0) | K mean | K−2Θ mean | θ₊ RMS on trial |',
        '|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['hole']==0 and r['r_cyl']!='':text.append('| '+' | '.join(fmt(r[k]) for k in ('time','chi_mean','K_mean','K_minus_2Theta_mean','theta_rms'))+' |')
    text+=['','Angular N96/N192 spreads are reported as sampling differences, not full spatial error bounds '
        '([t25-early-refined-angular-spreads.csv](t25-early-refined-angular-spreads.csv)). '
        f"The largest A/Q spreads on these ten trials are {max(r['absolute_spread'] for r in spreads if r['quantity']=='A'):.3g} / "
        f"{max(r['absolute_spread'] for r in spreads if r['quantity']=='Q'):.3g}. "
        'The signed θ₊ RMS grows from about 0.008 at t=0 to 0.137/0.276/0.46/0.57; its tiny mean is cancellation, not a small pointwise residual. '
        'The map uses the actual metric inverse, while the author finder uses unit-determinant conformal cofactors. '
        'Their maximum θ₊ discrepancy on the refined hole-0 trial grows from 8.7e-10 to 0.00449, 0.00671, 0.0604 and 0.204. '
        'This measured geometry-prescription difference is an additional limitation on late finder/map comparisons, not a silently corrected determinant.','',
        '### What the time series establishes','',
        'The initial narrow, grid-supported barrier and repeated capped flow support the **early finder/seeding/stopping part of case (i)**. '
        'Pointwise barriers remain sampled through t=10.5, but their angular zero bands are broad and no solved near-zero surface is obtained. '
        'At t=14 the evidence weakens to an average-only crossing; at t=17.5, 21, 24.5 and 28 none of the eligible sampled individual families has '
        'a negative-to-positive pointwise or mean pair. This is an absence within the tested families, not proof that no horizon exists. '
        'Thus a single one of the consult’s three cases cannot be asserted for the entire window: case (iii) “even early” is contradicted by the early barriers, '
        'and case (ii) is not established by qualified horizons.','',
        'The mean-root trial first falls below ten cells at the sampled t=14 (9.687/9.689 cells). It never supplies a resolved below-three-cell root. '
        '**The times at which an actual horizon falls below ten or three cells are unknown.** '
        'There is no accepted horizon whose area stays about 0.2245 while its coordinate radius squeezes; absence of later brackets cannot supply that claim.','',
        'The initial refined trial gives A=0.22208851 and Q=1.06932399. Its area is about 1.07% below the isolated progenitor value 0.2245, '
        'but agrees with the independent d=16 CTT initial reference A_i=0.2220848223, Q_i=1.06931005935 within 1.7e-5/1.3e-5 relative. '
        'Later trial A/Q do not remain at either initial pair: A is about 0.21384/0.19410 at t=7/10.5 and Q about 0.90922/0.71384. '
        'Those changes greatly exceed angular quadrature spreads, but the trial is not a MOTS and its pointwise residual increases. '
        '**Neither conservation nor a physical loss of horizon area/charge is certified by this series.** '
        'Nearest-cell coordinate-sphere areas and warm-tracking output are not substituted for these missing horizon measurements.','',
        '### Midpoint at t=28 and receipts','',
        'The 20 enclosing midpoint surfaces are spheres and prolates with c/a=1,2,4,8 and axial semiaxes 1.05,1.25,1.5,2,4 times '
        'the measured half separation (c=7.56963–28.83668 M_i). All enclose both punctures and exceed 16.47 cells per smallest extent. '
        'Every mean expansion is positive (minimum 0.06068095), but some angular sectors are negative (global minimum −0.06285768). '
        'No enclosing pointwise or mean barrier pair exists in this sampled set; no midpoint finder is run '
        '([t25-early-midpoint-t28.csv](t25-early-midpoint-t28.csv)). This is the requested sanity test, not a global horizon-exclusion proof.','',
        'Each native map/finder receipt was verified independently: its done marker, child return code, no gate, measured RSS and completion text. '
        'All ten finder return codes are 1 with finite TIME_CAP stage data; the workflow marker 0 only records completion. '
        f"The continuation peak RSS is {audit['peak_analysis_rss_bytes']/1e9:.6f} GB; initialization peaks at 2.867 GB. "
        'Maps have a 3 GB process cap, initializers/finders a 6 GB cap, with one process/thread and one checkpoint at a time. '
        'The regression preserves all 55,944 existing map/angle values bit for bit and verifies the weaker-seed gate off/on. '
        'See [t25-early-audit.json](t25-early-audit.json), [t25-early-final-audit.json](t25-early-final-audit.json), '
        '[t25-early-finder-stages.csv](t25-early-finder-stages.csv), [t25-early-receipts.csv](t25-early-receipts.csv) and '
        '[t25-early-regression.csv](t25-early-regression.csv).','',
        'To regenerate the final ledgers/report from completed receipts without probes, run `python3 t25-early-audit.py` then '
        '`python3 t25-early-finalize.py` then `python3 t25-manifest.py`. '
        'The serial pipeline and post-flow map markers are both 0. No job remains pending. '
        'No positive-time static input, evolution, SSH, source change, NaN-abort build or commit is part of this continuation.','']
    p=HERE/'t25-expansion-map.md';old=p.read_text().split(text[0])[0].rstrip()
    old=old.replace('time-resolved continuation below is a separate, pending bounded analysis.',
        'time-resolved continuation below is a completed bounded analysis, with no qualified binary horizon.')
    p.write_text(old+'\n\n'+'\n'.join(text)+'\n')
    heading='### T25 early-checkpoint continuation';p=HERE/'README.md';old=p.read_text().split(heading)[0].rstrip()
    p.write_text(old+'\n\n'+heading+'\n\n'+
        '**READY-EXCEPT; completed offline.** Eight SHA-256-verified early checkpoints and a matching d=16 t=0 regeneration are processed serially. '
        'The initial numerical-centre/spheroid pointwise bracket is 1.0025% wide; at t=3.5/7/10.5 the refined widths remain about 38/31/22%. '
        'Only a weaker mean crossing remains at t=14; no pair is found in eligible individual families at t=17.5–28. '
        'All ten unchanged N48, chase=quota=1, 235 s searches hit TIME_CAP in stage 1. '
        'The trial extent first falls below ten cells at t=14; an actual horizon’s support-loss times and A/Q preservation are unestablished. '
        'This supports an early finder problem, with later map ambiguity; it does not establish coordinate squeezing. '
        'The 20 enclosing midpoint t=28 surfaces give no barrier pair and have no finder.\n\n'
        'See [t25-early-refined-time-table.csv](t25-early-refined-time-table.csv), '
        '[t25-early-final-brackets.csv](t25-early-final-brackets.csv), per-checkpoint t25-early-step*.csv, '
        '[t25-early-finder-stages.csv](t25-early-finder-stages.csv), '
        '[t25-early-final-audit.json](t25-early-final-audit.json) and '
        '[t25-expansion-map.md](t25-expansion-map.md). The 55,944-value disabled/unchanged-output regression passes; '
        'peak RSS 3.322085 GB, within 6 GB (maps capped at 3 GB), one thread, one checkpoint at a time. '
        'Pipeline and postmap markers are 0, no job remains pending. No Source/Examples edit, NaN-abort build, evolution, SSH or commit was made.\n')

if __name__=='__main__':main()
