#!/usr/bin/env python3
"""Qualify each frozen numerical finder probe from its own native receipt."""
import csv
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t24')
E=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall/ev1/runs/exp-0024/merger')


def main():
    allrows=[];receipt_rows=[];inputs=[]
    for label,name,directory in [('individual','t24-horizon-search','horizon'),
        ('common N48 chase1','t24-horizon-common-pilot','horizon-common-pilot'),
        ('common N96 chase1','t24-horizon-common-N96','horizon-common-N96'),
        ('common N48 chase.125','t24-horizon-common-slow','horizon-common-slow')]:
        q=json.loads((HERE/(name+'.json')).read_text());r=q['resources']
        assert (ROOT/directory/'done.exit').read_text().strip()==str(r['returncode'])=='1'
        own=json.loads((ROOT/directory/'horizon.resources.json').read_text())
        assert own==r and r['child_measurement']['returncode']==1 and r['peak_rss_bytes']<6e9 and not r['gate_reason']
        assert q['advances']==q['static_data_reads']==0
        inputs.extend([str(ROOT/directory/'horizon.resources.json'),str(ROOT/directory/'params.txt'),str(ROOT/directory/'run.log')])
        rows=[row for row in q['rows'] if row['stage']=='0']
        for row in rows:
            row=dict(row,probe=label,qualification='UNQUALIFIED; extents/A/Q are trial-surface values')
            assert float(row['expansion_squared'])>1e-7 and row['status'] in ('UPDATE_CAP','TIME_CAP','FLOOR')
            allrows.append(row)
        receipt_rows.append(dict(probe=label,native_exit=1,gate='',seeds=len(rows),wall_s=r['wall_seconds'],
            peak_RSS_bytes=r['peak_rss_bytes'],best_final_squared_residual=min(float(row['expansion_squared']) for row in rows)))
    with (HERE/'t24-horizon-all-probes.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=allrows[0]);w.writeheader();w.writerows(allrows)
    with (HERE/'t24-horizon-probe-receipts.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=receipt_rows[0]);w.writeheader();w.writerows(receipt_rows)
    history=[r for r in csv.DictReader((E/'horizon-history.csv').open()) if r['kind']=='common']
    assert len(history)==18 and all(r['status']=='TIME_CAP' and r['label']=='UNRESOLVED' for r in history)
    with (HERE/'t24-common-horizon-history.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=history[0]);w.writeheader();w.writerows(history)
    # Matched accumulated chase factor (72 vs 574*.125=71.75).
    def at(d,idx,step):
        return next(q for q in csv.DictReader((ROOT/d/'progress.csv').open()) if int(q['search_index'])==idx and int(q['update'])==step)
    fast=at('horizon-common-pilot',6,72);slow=at('horizon-common-slow',1,574)
    # The collected ev1 packet has no original finder diagnostics/shapes/log.
    files=[str(p.relative_to(E)) for p in E.rglob('*') if p.is_file()]
    missing=[p for p in ('diagnostics','shapes') if not (E/p).exists()]
    assert missing==['diagnostics','shapes']
    data=dict(status='NO_STRICT_INDIVIDUAL_OR_COMMON_SURFACE_IN_BOUNDED_PROBES',receipts=receipt_rows,
        qualified_area=None,qualified_charge=None,horizon_membership='UNDETERMINED at rho=.047-.08',
        missing_original_ev1_finder_directories=missing,original_common_rows=len(history),
        matched_chase_comparison=dict(seed_radius=3,N=48,fast_update=72,fast_chase=1,fast=fast,slow_update=574,slow_chase=.125,slow=slow,
            conclusion='flow-step sensitivity measured; neither is a horizon; accumulated chase differs by 0.35%'),
        own_receipts=inputs,author_finder_unchanged=True,checkpoint_time_M=133,static_inputs_read=False,
        native_binary_sha256=hashlib.file_digest(Path('/private/tmp/ems-t17/rh.ex').open('rb'),'sha256').hexdigest(),
        limitation='bounded failures do not prove absence; no qualified A/Q or inside/outside classification')
    (HERE/'t24-horizon-summary.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(receipt_rows,indent=2))
if __name__=='__main__':main()
