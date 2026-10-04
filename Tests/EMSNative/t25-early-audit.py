#!/usr/bin/env python3
"""Final receipt/geometry audit, without new probes or checkpoint field reads."""
import csv,importlib.util,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('early',HERE/'t25-early.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def main():
    assert (e.ROOT/'finders-complete.json').exists()
    assert int((e.ROOT/'pipeline/done.exit').read_text())==0
    e.summarize();summaries=[];angles=[];stage_rows=[];counts=[];peaks=[]
    for s in e.samples():
        d=e.ROOT/s['tag'];selection=json.loads((d/'selection.json').read_text())
        full=e.read(d/'map-N192.csv')
        assert len({r['id'] for r in full})==len(full)
        statuses={k:sum(r['status']==k for r in full) for k in ('RESOLVED','UNRESOLVED','INVALID_GEOMETRY')}
        for r in full:
            if r['status']=='RESOLVED':
                assert float(r['min_extent_over_dx'])>=3
                assert all(math.isfinite(float(r[k])) for k in ('theta_min','theta_max','theta_mean','theta_rms','A','Q','K_mean','K_minus_2Theta_mean'))
        counts.append(dict(time=s['time'],**statuses,pointwise_families=len(e.strong(full)),
            mean_only_pairs=sum(x['status']=='AVERAGE_ONLY_WEAKER' for x in e.b.brackets(full))))
        for item in selection['selections']:
            root=item['root'];br=item['bracket'];n48=[];result=None
            if br:
                sd=Path(item['seed_dir']);fd=sd/'finder'
                rc=int((fd/'done.exit').read_text());receipt=e.verified(fd,'early-finder',rc,'T25_FINDER_COMPLETE',6000000000)
                peaks.append(receipt['peak_rss_bytes']);n48=e.read(fd/'finder.csv');last=n48[-1]
                assert rc in (0,1) and int(last['N'])==48
                assert all(math.isfinite(float(r[k])) for r in n48 for k in ('expansion_squared','A','Q','r_min','r_max'))
                if rc==0:assert len(n48)==3 and all(r['status']=='FOUND' and float(r['expansion_squared'])<=float(r['threshold']) for r in n48)
                else:assert last['status']=='TIME_CAP' and float(last['seconds'])>=235
                for stage in n48:stage_rows.append(dict(hole=item['hole'],configuration=s['configuration'],**stage))
                root96=item['angular_root_N96']
                for k in ('A','Q','lapse_mean','chi_mean','K_mean','K_minus_2Theta_mean','theta_mean','theta_rms'):
                    angles.append(dict(time=s['time'],hole=item['hole'],quantity=k,N96=root96[k],N192=root[k],
                        absolute_spread=abs(float(root96[k])-float(root[k]))))
                result=dict(status=last['status'],squared=float(last['expansion_squared']),A=float(last['A']),Q=float(last['Q']),
                    r_min=float(last['r_min']),r_max=float(last['r_max']),finest_cells_lower_bound=float(last['r_min'])/e.H,
                    seconds=float(last['seconds']),peak_rss_bytes=receipt['peak_rss_bytes'],stages_qualified=sum(r['status']=='FOUND' for r in n48),
                    under_3_finest_cells=float(last['r_min'])/e.H<3)
            summaries.append(dict(time=s['time'],hole=item['hole'],bracket_status=br['status'] if br else 'NO_RESOLVED_BRACKET',
                bracket_fractional_width=item.get('barrier_fractional_width'),root=root,finder=result,
                interpretations='Trial data only unless finder qualifies; no absent-bracket inference below 3 cells'))
        packet=e.read(HERE/f't25-early-{s["tag"]}.csv')
        existing={r['id'] for r in packet}
        packet += [x['root'] for x in selection['selections'] if x['root'] and x['root']['id'] not in existing]
        e.write(HERE/f't25-early-{s["tag"]}.csv',packet)
    for p in e.ROOT.rglob('early-map.resources.json'):
        r=json.loads(p.read_text());assert r['returncode']==0 and not r['gate_reason'] and r['peak_rss_bytes']<3e9
        e.verified(p.parent,'early-map',0,'T25_EXPANSION_MAP_COMPLETE; no advances');peaks.append(r['peak_rss_bytes'])
    e.write(HERE/'t25-early-map-counts.csv',counts)
    e.write(HERE/'t25-early-angular-spreads.csv',angles)
    e.write(HERE/'t25-early-finder-stages.csv',stage_rows)
    midpoint=e.read(HERE/'t25-early-midpoint-t28.csv')
    assert len(midpoint)==20 and all(r['encloses_0']==r['encloses_1']=='1' for r in midpoint)
    data=dict(status='COMPLETE',all_native_receipts_verified=True,no_evolution=True,inputs=json.loads((e.ROOT/'inputs.json').read_text()),
        rows=summaries,counts=counts,physical_probes=sum(x['finder'] is not None for x in summaries),
        all_three_stage_passes=sum(x['finder'] and x['finder']['stages_qualified']==3 or False for x in summaries),
        peak_analysis_rss_bytes=max(peaks),midpoint_pointwise_pairs=len(e.strong(midpoint)),
        midpoint_average_pairs=len(e.b.brackets(midpoint)))
    e.atomic_json(HERE/'t25-early-audit.json',data)
    print('T25_EARLY_ALL_RECEIPTS_VERIFIED',data['physical_probes'],data['all_three_stage_passes'],max(peaks))
if __name__=='__main__':main()
