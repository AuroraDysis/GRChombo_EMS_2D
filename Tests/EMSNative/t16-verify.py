#!/usr/bin/env python3
"""Verify completed receipts/caches and summarize the frozen T16 comparison."""
import sys
sys.dont_write_bytecode = True
import csv, hashlib, json, math, shutil
from pathlib import Path
import numpy as np
import importlib.util
HERE = Path(__file__).resolve().parent
ROOT = Path('/private/tmp/ems-t16')
spec = importlib.util.spec_from_file_location('t16analysis', HERE/'t16-analyze.py')
a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
FIELDS = ('Gamma', 'metric_Gamma', 'shift', 'lapse')
def rows(name):
    return list(csv.DictReader((HERE/name).open()))
def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()
def receipt(path):
    r = json.loads(path.read_text())
    assert r['returncode'] == 0 and not r['gate_reason'], path
    assert r['child_measurement']['returncode'] == 0, path
    assert r['peak_rss_bytes'] < 3_000_000_000 and r['tree_peak_bytes'] < 3_000_000_000, path
    assert all(max(v) < 3_000_000_000 for v in r['per_pid_peaks'].values()), path
    return r
def main():
    checks = []
    def checked(check, count, evidence):
        checks.append(dict(check=check, count=count, result='PASS', evidence=evidence))
    # Verify the pipeline's pinned numerical inputs before refreshing its manifest.
    pinned = ROOT/'verification/pipeline-manifest.txt'
    pinned.parent.mkdir(exist_ok=True)
    if not pinned.exists(): shutil.copyfile(HERE/'t16-manifest.txt', pinned)
    count = 0
    for line in pinned.read_text().splitlines()[2:]:
        sha, size, name = line.split(' ', 2); path = Path(name)
        if path.name in ('t16-manifest.txt', 't16-analyze.py', 't16-report.py', 'README.md', 't16-status.json', 't16-resources.csv'):
            continue  # report/manifest refresh only; numerical artifacts must remain pinned
        assert path.is_file() and path.stat().st_size == int(size) and digest(path) == sha, path
        count += 1
    checked('pipeline pinned inputs, executables, captures and caches', count, str(pinned))
    for directory, label in [('evolution', 'evolution-worker'), ('pipeline', 'queue'), ('pipeline/job', 'isolated-pipeline')]:
        d = ROOT/directory
        assert (d/'done.exit').read_text().strip() == '0'
        r = receipt(d/(label+'.resources.json'))
        checked('completed '+directory, 1, f"child rc=0, no gate; RSS={r['peak_rss_bytes']}; tree={r['tree_peak_bytes']}")
    for directory, label in [('evolution','evolution-worker'), ('analysis','analysis'), ('report','report'), ('manifest','manifest')]:
        receipt(ROOT/directory/(label+'.resources.json'))
    checked('pipeline component receipts', 4, 'successful child receipts; timer rc=1 is sandbox kern.clockrate denial')
    expected = np.arange(75) * (.875/8192/4)
    registration = {r['run']:r for r in rows('t16-registration.csv')}
    audit = {r['run']:r for r in rows('t16-launch-audits.csv')}
    hist = {(r['run'],r['ray'],r['field'],r['norm'],int(r['step'])):r for r in rows('t16-launch-history.csv')}
    replay_values = 0
    for name, reg in registration.items():
        d = ROOT/'evolution'/name; r = receipt(d/(name+'.resources.json'))
        assert (d/'done.exit').read_text().strip() == '0'
        assert json.loads(r['command'])[-1] == reg['parameters']
        log = (d/'run.log').read_text()
        assert 'GRChombo finished.' in log and 'T13 clean native stop' in log
        assert 'MayDay' not in log and 'nan' not in log.lower() and 'Read EMSTRUMPET' not in log[log.index('GRAMRLevel::advance level 0 at time 0'):]
        stop = next(csv.DictReader((d/'t13-stop.csv').open()))
        assert float(stop['actual_time_M']) == float(reg['native_t_stop_M'])
        meta = np.load(ROOT/'analysis'/(name+'-meta.npz')); times = meta['times']; cells = meta['cells']
        assert np.array_equal(times, expected)
        q = np.memmap(ROOT/'analysis'/(name+'-q.bin'), mode='r', dtype='f8', shape=(75,len(cells),180))
        assert np.isfinite(q).all() and np.count_nonzero(q[:,:,144]) == 0
        replay_values += q.size
        floor = [v for f in d.glob('t13-floors-*.csv') for v in csv.DictReader(f.open())]
        assert floor and all(int(v[k]) == 0 for v in floor for k in ('chi_activations','lapse_activations','nonfinite'))
        assert int(audit[name]['snapshots']) == 75 and float(audit[name]['last_retained_time_M']) == expected[-1]
        for ray, vec in a.a.NV.items():
            p = np.load(ROOT/'analysis'/(name+'-'+ray+'-profiles.npz'))
            assert np.array_equal(p['times'], expected) and np.array_equal(p['r'], a.a.R)
            f8 = [a.a.fields(frame, vec) for frame in p['p8']]
            f10 = [a.a.fields(frame, vec) for frame in p['p10']]
            for step in range(75):
                for field in (*FIELDS, 'C_Gamma', 'Ham', 'Mom', 'GaussE'):
                    u = f8[step][field] - (f8[0][field] if field in FIELDS else 0)
                    v = f10[step][field] - (f10[0][field] if field in FIELDS else 0)
                    for norm, fn in [('peak',lambda x:float(np.max(abs(x)))), ('RMS',a.a.rms)]:
                        h = hist[name,ray,field,norm,step]
                        assert fn(u) == float(h['amplitude']) and fn(u-v) == float(h['interpolation_spread'])
        checked(name+' native completion', 75, f"stop={stop['actual_time_M']}; finite replay; Base RHS bit mismatches=0; floors=0; reader after advance=0; gate clear")
    checked('all replay Float64 values and profile-derived CSV norms', replay_values, 'four independent native legs; 9600 history rows')
    margins=[]; unqualified=[]; summary=[]; trend=[]
    for r in rows('t16-launch-ratios.csv'):
        level=int(r['level']); step=int(r['step']); key=(r['ray'],r['field'],r['norm'],step)
        g=hist[f'L{level}-geometric',*key]; s=hist[f'L{level}-sqrt',*key]
        G=float(g['amplitude']); S=float(s['amplitude']); eg=5*float(g['interpolation_spread']); es=5*float(s['interpolation_spread'])
        assert G/S == float(r['geometric_over_sqrt'])
        low=max(0,G-eg)/(S+es); high=(G+eg)/(S-es) if S>es else math.inf
        margin=(S-G)/(eg+es) if eg+es else math.inf
        qualified=G>eg and S>es
        assert qualified == (r['qualified']=='True')
        m=dict(level=level,ray=r['ray'],field=r['field'],norm=r['norm'],step=step,time_M=r['time_M'],
               ratio=G/S,ratio_lower_5x=low,ratio_upper_5x=high,reduction_margin_5x=margin,
               geometric_margin_5x=r['geometric_margin_5x'],sqrt_margin_5x=r['sqrt_margin_5x'],qualified=qualified,
               reduction_resolved=margin>1,endpoint=step==74)
        margins.append(m)
        if r['field'] in FIELDS and not qualified: unqualified.append(m)
    a.save('t16-sampling-margins.csv',margins)
    for level in (13,14):
        for ray in ('axis','diagonal'):
            for field in FIELDS:
                for norm in ('peak','RMS'):
                    rs=[r for r in margins if (r['level'],r['ray'],r['field'],r['norm'])==(level,ray,field,norm)]
                    assert len(rs)==74
                    summary.append(dict(level=level,ray=ray,field=field,norm=norm,positive_clocks=74,
                        qualified=sum(r['qualified'] for r in rs),uncertainty_dominated=sum(not r['qualified'] for r in rs),
                        reduction_resolved=sum(r['reduction_resolved'] for r in rs),central_improvement=sum(r['ratio']<1 for r in rs),
                        ratio_min=min(r['ratio'] for r in rs),ratio_max=max(r['ratio'] for r in rs),
                        enhanced_steps=' '.join(str(r['step']) for r in rs if r['ratio']>1),
                        geometric_history_sup=max(float(hist[f'L{level}-geometric',ray,field,norm,i]['amplitude']) for i in range(75)),
                        sqrt_history_sup=max(float(hist[f'L{level}-sqrt',ray,field,norm,i]['amplitude']) for i in range(75)),
                        min_signal_margin_5x=min(min(float(r['geometric_margin_5x']),float(r['sqrt_margin_5x'])) for r in rs),
                        min_reduction_margin_5x=min(r['reduction_margin_5x'] for r in rs),
                        endpoint_lower_5x=rs[-1]['ratio_lower_5x'],endpoint_upper_5x=rs[-1]['ratio_upper_5x']))
    a.save('t16-history-summary.csv',summary)
    constraints=[]
    for level in (13,14):
        for ray in ('axis','diagonal'):
            for field in ('C_Gamma','Ham','Mom','GaussE'):
                for norm in ('peak','RMS'):
                    rs=[r for r in margins if (r['level'],r['ray'],r['field'],r['norm'])==(level,ray,field,norm)]
                    e=rs[-1]
                    constraints.append(dict(level=level,ray=ray,field=field,norm=norm,
                        geometric_endpoint=hist[f'L{level}-geometric',ray,field,norm,74]['amplitude'],
                        sqrt_endpoint=hist[f'L{level}-sqrt',ray,field,norm,74]['amplitude'],
                        endpoint_ratio=e['ratio'],endpoint_lower_5x=e['ratio_lower_5x'],endpoint_upper_5x=e['ratio_upper_5x'],
                        endpoint_reduction_margin_5x=e['reduction_margin_5x'],qualified=sum(r['qualified'] for r in rs),
                        history_ratio_min=min(r['ratio'] for r in rs),history_ratio_max=max(r['ratio'] for r in rs)))
    a.save('t16-constraint-summary.csv',constraints)
    with (HERE/'t16-unqualified-history.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=margins[0]);w.writeheader();w.writerows(unqualified)
    for ray in ('axis','diagonal'):
        for mode in ('sqrt','geometric'):
            for field in FIELDS:
                for norm in ('peak','RMS'):
                    coarse=hist[f'L13-{mode}',ray,field,norm,74]; fine=hist[f'L14-{mode}',ray,field,norm,74]
                    trend.append(dict(ray=ray,mode=mode,field=field,norm=norm,coarse_amplitude=coarse['amplitude'],fine_amplitude=fine['amplitude'],
                        fine_over_coarse=float(fine['amplitude'])/float(coarse['amplitude']),interpretation='combined spatial and temporal sensitivity; not a temporal error bound'))
    a.save('t16-rung-trend.csv',trend)
    native={(r['run'],int(r['step']),r['field']):r for r in rows('t16-native-puncture.csv')}
    punct=[]
    for level in (13,14):
        for step in range(1,75):
            for field in ('Gamma1','Gamma2','shift1','shift2','lapse'):
                s=native[f'L{level}-sqrt',step,field];g=native[f'L{level}-geometric',step,field]
                S=float(s['puncture_peak_change']);G=float(g['puncture_peak_change'])
                punct.append(dict(level=level,step=step,time_M=g['time_M'],field=field,sqrt_peak_change=S,geometric_peak_change=G,
                    geometric_over_sqrt=G/S if S else '',sqrt_KO_peak=s['puncture_KO_peak'],geometric_KO_peak=g['puncture_KO_peak'],endpoint=step==74))
    a.save('t16-puncture-ratios.csv',punct)
    for name,count_field in [('t16-controls.csv','bit_mismatches'),('t16-nonlapse-identity.csv','bit_mismatches')]:
        assert all(int(r[count_field])==0 for r in rows(name)); checked(name,len(rows(name)),'all bit mismatches zero')
    for path in ROOT.rglob('*.resources.json'):
        r=json.loads(path.read_text()); assert r['peak_rss_bytes']<3_000_000_000 and r['tree_peak_bytes']<3_000_000_000
    checked('all historical resource receipts',len(list(ROOT.rglob('*.resources.json'))),'per-process and active-tree gates <3 GB, including retained failed setup attempts')
    a.save('t16-verification.csv',checks)
    print(json.dumps(dict(checks=len(checks),history_qualified=sum(r['qualified'] for r in summary),
        history_total=32*74,unqualified=len(unqualified),sampling_resolved_reductions=sum(r['reduction_resolved'] for r in summary),
        min_signal_margin=min(r['min_signal_margin_5x'] for r in summary),min_reduction_margin=min(r['min_reduction_margin_5x'] for r in summary)),indent=2))
    print('Enhanced shift steps:',[(r['level'],r['ray'],r['norm'],r['enhanced_steps']) for r in summary if r['field']=='shift'])
    print('Puncture endpoint:',[(r['level'],r['field'],r['geometric_over_sqrt']) for r in punct if r['endpoint']])
if __name__=='__main__': main()
