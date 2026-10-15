"""Measured within-route interpolation-order error budget, fixed before the run."""
import collections
import csv
import json
import statistics
from pathlib import Path
import report_t2 as base

HERE=Path(__file__).resolve().parent
OUT=HERE/'results-t4'


def main():
    contract=json.loads((OUT/'contract.json').read_text())
    assert contract['safety_factor']==2 and contract['target_order']==4
    report=[]
    baseline=[]
    for form in ('vacuum','matter'):
        data={}
        for order in contract['comparison_orders']:
            folder=OUT/f'study/order{order}'/form
            for route,file,name in [('np','ems_np.csv','Psi4'),('W','ems_radiation.csv','W')]:
                rows=base.rows(folder/file)
                assert len(rows)==21 and all(d['t']==602 and d['m']==0 for d in rows)
                data[route,order]={(d['r'],int(d['l'])):base.value(d,name) for d in rows}
                if order==4:
                    previous=base.rows(HERE/'results-t3/checkpoint'/form/file)
                    plain=base.rows(OUT/'study/production-cost'/form/file)
                    for label,reference in [('T3',previous),('unchecked',plain)]:
                        assert len(reference)==len(rows)
                        maximum=0
                        for col in rows[0]:
                            scale=max(abs(d[col]) for d in reference)
                            tolerance=4096*base.EPS*scale+1e-100
                            for a,b in zip(rows,reference):
                                maximum=max(maximum,abs(a[col]-b[col])/tolerance)
                        baseline.append({'form':form,'route':route,'reference':label,
                            'max_error_over_budget':maximum,'status':'PASS' if maximum<=1 else 'FAIL'})
        for r in contract['radii']:
            for l in contract['modes']:
                n={o:data['np',o][r,l] for o in (2,4,6)}
                w={o:data['W',o][r,l] for o in (2,4,6)}
                error={o:abs(n[o]-w[o]) for o in (2,4,6)}
                np46=abs(n[4]-n[6]);w46=abs(w[4]-w[6])
                estimate=np46+w46
                tolerance=contract['safety_factor']*estimate
                amplitude=abs(n[4]); assert amplitude>0
                report.append({'form':form,'r':r,'l':l,
                    'NP':{o:[v.real,v.imag] for o,v in n.items()},
                    'production':{o:[v.real,v.imag] for o,v in w.items()},
                    'route_difference':error,'NP_4_to_6':np46,'W_4_to_6':w46,
                    'NP_2_to_4':abs(n[2]-n[4]),'W_2_to_4':abs(w[2]-w[4]),
                    'estimate':estimate,'tolerance':tolerance,'mode_amplitude':amplitude,
                    'required_safety_factor':error[4]/estimate if estimate else None,
                    'relative_difference':error[4]/amplitude,
                    'relative_tolerance':tolerance/amplitude,
                    'error_over_tolerance':error[4]/tolerance if tolerance else None,
                    'status':'PASS' if error[4]<=tolerance else 'FAIL'})
    with (OUT/'study/cost.csv').open() as stream:
        costs=list(csv.DictReader(stream))
    timing=[]
    for form,matter in [('vacuum','0'),('matter','1')]:
        for route in ('production','np'):
            values=[float(d['seconds']) for d in costs if d['matter']==matter and
                d['checked']=='0' and d['route']==route and int(d['sample'])>=0]
            assert len(values)==3 and all(v>0 for v in values)
            timing.append({'form':form,'route':route,'samples':len(values),
                'median_seconds':statistics.median(values),'min_seconds':min(values),'max_seconds':max(values)})
    accepted=json.loads((HERE/'results-t3/validation.json').read_text())
    statuses=[v['status'] for d in accepted['static'] for v in d['scalars'].values()]
    statuses += [v['status'] for d in accepted['rotation'] for v in d['rotations'].values()]
    statuses += [accepted[k]['status'] for k in ('default','projections')]
    statuses += [accepted['radiative_preservation_status']]
    assert all(s=='PASS' for s in statuses)
    ready=all(d['status']=='PASS' for d in report+baseline)
    result={'contract':contract,'comparison':report,'baseline_reproduction':baseline,
        'timing':timing,'accepted_T3_controls':'PASS; old Weyl4 tolerance superseded by T4',
        'L6_difference_larger_than_L4':sum(d['route_difference'][6]>d['route_difference'][4] for d in report),
        'ready_for_simulation':ready}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# NPTET T4 interpolation-order comparison','',f'**Ready: {ready}.**','',
        'The production order remains Lagrange<4>. The estimator is the sum of two within-route changes: E = |NP4-NP6| + |W4-W6|. The acceptance is |NP4-W4| <= 2 E. Safety factor 2 was fixed in contract.json before computing the study. All production W cell derivatives remain fourth order; only the interpolation order changes. This is a measured discretization estimate at one checkpoint, not a rigorous error bound.','',
        'The original three ghosts and point transfer are retained. The existing Lagrange stencil search is bounded to the sampled FAB, shifting full-width stencils near its edges. Every sample is bounds checked; the unbounded L6 negative control rejects an out-of-array stencil. All valid checkpoint data and hierarchy levels are unchanged. L4 reproduces T3, and checked/unchecked L4 agree within 4096 eps times each column scale.','',
        '| form | r | l | difference L2 | difference L4 | difference L6 | NP4-NP6 | W4-W6 | 2E | relative difference | relative budget | result |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for d in report:
        vals=[d['route_difference'][o] for o in (2,4,6)]+[d['NP_4_to_6'],d['W_4_to_6'],d['tolerance'],d['relative_difference'],d['relative_tolerance']]
        lines.append(f"| {d['form']} | {d['r']} | {d['l']} | "+' | '.join(f'{v:.6e}' for v in vals)+f" | {d['status']} |")
    lines += ['',
        'L2-to-L4 reduces every cross-route discrepancy, but the L6 discrepancy is larger than L4 in all 18 comparisons. Increasing only interpolation order leaves the fourth-order cell-Weyl derivative error in place. The l=6,r=75 pair in both forms needs a safety factor about 2.95 on the declared within-route estimate, exceeding the allowed maximum 2. No L2-based larger allowance is substituted after observing this failure.']
    lines += ['','## Extraction cost on this checkpoint','',
        'Uninstrumented production L4, 3 spheres, 97 theta nodes, serial Chombo, OMP_NUM_THREADS=8. Median and range of three calls after one warm-up; includes construction, interpolation, reconstruction, projection and CSV IO. Excludes restart loading, shared ghost refresh and shared cell-Weyl preparation. NP is the incremental cost when existing EMS radiation is already enabled.','',
        '| form | route | median seconds | min | max |','|---|---|---:|---:|---:|']
    for d in timing:lines.append(f"| {d['form']} | {d['route']} | {d['median_seconds']:.6f} | {d['min_seconds']:.6f} | {d['max_seconds']:.6f} |")
    lines += ['','All raw modes, per-call times, parameters, input hashes and logs accompany this report. T3 static, algebraic, rotation, normalization and default-output controls remain accepted; only the Weyl4 comparison tolerance is superseded.']
    (OUT/'summary.md').write_text('\n'.join(lines)+'\n')
    print('Comparison:',collections.Counter(d['status'] for d in report))
    print('Baseline:',collections.Counter(d['status'] for d in baseline))
    print('Ready:',ready)
    return 0 if ready else 1


if __name__=='__main__':raise SystemExit(main())
