#!/usr/bin/env python3
"""Read-only verification of collected exp0026 receipts/tables/manifests."""
import csv
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0026/submission')


def main():
    verified=[];inputs=[]
    for name in ('S-history-59294736','G-D-59294962'):
        root=BASE/name
        for line in (root/'manifest.sha256').read_text().splitlines():
            digest,rel=line.split(None,1);path=root/rel.strip().lstrip('*')
            assert path.is_file(),path
            actual=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
            assert actual==digest,(path,actual,digest)
            verified.append(str(path))
    history=BASE/'S-history-59294736'
    reduction=json.loads((history/'reduction.receipt.json').read_text())
    assert reduction['status']=='PASS' and reduction['rows']==112 and reduction['invalid_metric_active']==0
    assert not reduction['static_inputs_read'] and not reduction['evolution_executed']
    h=list(csv.DictReader((history/'S-history.csv').open()));assert len(h)==112
    summary=[]
    for label,step in [('exp-0024',152),('P',153),('P',154)]:
        r=[q for q in h if q['label']==label and int(q['step'])==step]
        assert len(r)==7 and all(int(q['invalid_metric_active'])==0 for q in r)
        maximum=max(r,key=lambda q:float(q['max_S']))
        summary.append(dict(leg=label,step=step,time_M=maximum['time_M'],max_S=maximum['max_S'],
            max_level=maximum['level'],cells_S_gt_4p9=sum(int(q['cells_S_gt_4p9']) for q in r),
            nonfinite_valid='not in this table',native_rank_exits='reduction only',source='d507438'))
    for name in ('G','D','D-dt-readback'):
        root=BASE/'G-D-59294962'/name
        n=json.loads((root/'native.receipt.json').read_text());r=json.loads((root/'intervention.receipt.json').read_text())
        assert n['reason']=='NATIVE_EXIT' and n['launcher_returncode']==0 and not n['wrapper_error']
        assert len(n['native_rank_exits'])==224 and all(v=='0' for v in n['native_rank_exits'].values())
        assert r['receipt_complete'] and r['native_dt_verified']
        inputs.extend([str(root/'native.receipt.json'),str(root/'intervention.receipt.json')])
        if name=='D-dt-readback':
            assert r['zero_advance_verified'] and len(r['level_dt_readback'])==13
            continue
        assert r['source']=='d507438'
        assert not r['unreadable_checkpoints'] and r['static_inputs_absent']
        assert all(q['valid_nonfinite']==0 for q in r['checkpoints'])
        assert r['completed_coarse155'] and r['passed_135p174'] and r['last_finite_step']==158
        t=list(csv.DictReader((root/'S-per-checkpoint.csv').open()))
        endpoint=[q for q in t if int(q['step'])==158]
        maximum=max(endpoint,key=lambda q:float(q['max_S']))
        summary.append(dict(leg=name,step=158,time_M=r['last_finite_time_M'],max_S=maximum['max_S'],
            max_level=maximum['level'],cells_S_gt_4p9=sum(int(q['cells_S_gt_4p9']) for q in endpoint),
            nonfinite_valid=0,native_rank_exits='224/224 zero',source=r['source']))
    with (HERE/'t24-cluster-summary.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
    (HERE/'t24-cluster-check.json').write_text(json.dumps(dict(status='VERIFIED',manifest_files=len(verified),
        verified_files=verified,native_receipts=inputs,summary=summary,read_only=True,
        limitation='raw G/D checkpoints remain on cluster; collected native reductions and their receipts checked'),indent=2)+'\n')
    print('CLUSTER_COLLECTED_EVIDENCE_VERIFIED',len(verified),summary)
if __name__=='__main__':main()
