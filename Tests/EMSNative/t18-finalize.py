#!/usr/bin/env python3
"""Serial, detached full-window trials, with resource receipts and final marker."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t18')
PY='/Users/auroradysis/miniconda3/bin/python'
spec=importlib.util.spec_from_file_location('t18work',HERE/'t18-work.py')
work=importlib.util.module_from_spec(spec);spec.loader.exec_module(work)
replace=work.module('params',HERE/'t13-prepare.py').replace

def prepare():
    jobs=[]
    for dt in (.25,.375,.5):
        name='dt-'+str(dt)+'-window'
        old=(HERE/('t18-dt-'+str(dt)+'.txt')).read_text()
        old='\n'.join(line for line in old.splitlines() if not line.startswith('ems_ctt_data_path ='))+'\n'
        text=replace(old,dict(max_steps=3,stop_time=2,t18_stop_time=.875))
        param=HERE/('t18-'+name+'.txt');param.write_text(text)
        d=ROOT/'local'/name
        for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
        jobs.append(dict(name=name,directory=str(d),command=['env','OMP_NUM_THREADS=2',str(ROOT/'native-window.ex'),str(param)]))
    record=dict(window_M=.875,window_Rh=.875/.0060404520035922523,
        finest_native_steps=[8192,5461,4096],threads=2,process_RSS_cap_bytes=4000000000,
        physical_boundary_nearest_M=40,initial_data_after_t0='forbidden by t2 guard',
        launch_observer='pointwise native t0 subtraction on 0.00075<=r<=0.0025; every finest endpoint; '
            'Gamma Euclidean maximum, lapse absolute maximum, shift Euclidean maximum',
        launch_common_endpoint_M=.001922607421875,
        bounding_window='all valid cells, metric positivity, chi/lapse floors, puncture field maximum on fixed r<0.08',
        threshold_source='same predeclared 1e6 field, 2 shift/lapse, zero nonfinite and zero floor-crossing screen',
        prior_setup_failure='dt-0.25 never evolved: empty ems_ctt_data_path aborted ParmParse; receipt preserved',
        short_registration='superseded by this longer window before any dt evolution result',jobs=jobs)
    (HERE/'t18-window-registration.json').write_text(json.dumps(record,indent=2)+'\n')
    d=ROOT/'window';d.mkdir(exist_ok=True)
    # One worker owns the trials and classification; native failures do not skip
    # independent dt values. A memory gate terminates the series.
    plan=dict(jobs=[dict(name='window-analysis',directory=str(d/'job'),
        command=[PY,str(HERE/'t18-finalize.py'),'run'])])
    (d/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    print(d/'plan.json')

def run():
    registration=json.loads((HERE/'t18-window-registration.json').read_text())
    jobs=registration['jobs'];receipts=[]
    for job in jobs:
        d=Path(job['directory'])
        if not (d/'done.exit').exists():
            with (d/'run.log').open('wb') as log:
                rc=subprocess.run([PY,str(HERE/'t18-run.py'),'--measure',job['name'],
                    '--directory',str(d),'--',*job['command']],stdout=log,stderr=subprocess.STDOUT).returncode
            tmp=d/'done.exit.tmp';tmp.write_text(str(rc)+'\n');tmp.replace(d/'done.exit')
        receipt=json.loads((d/(job['name']+'.resources.json')).read_text())
        receipts.append(receipt)
        if receipt['gate_reason'] or receipt['peak_rss_bytes']>=4000000000:
            raise RuntimeError('memory/output gate; remaining trials not started')
    subprocess.run([PY,str(HERE/'t18-analyze.py')],check=True)
    with (HERE/'t18-dt-levels.csv').open() as f:
        levels=[r for r in csv.DictReader(f) if r['case'].endswith('-window')]
    failures=[dict(case=r['case'],level=r['level']) for r in levels if r['screen_pass']!='True']
    passed=all(r['returncode']==0 for r in receipts) and len(levels)==39 and not failures
    results=dict(status='COMPLETED' if passed else 'COMPLETED_WITH_FAILED_TRIALS',
        per_level_screen_failures=failures,
        process_peak_RSS_bytes=max(r['peak_rss_bytes'] for r in receipts),trials=[dict(case=r['process'],returncode=r['returncode'],
        wall_seconds=r['wall_seconds'],peak_RSS_bytes=r['peak_rss_bytes']) for r in receipts],
        scope='native boundedness screen and launch dt sensitivity, not a long-term merger admission')
    (HERE/'t18-window-results.json').write_text(json.dumps(results,indent=2)+'\n')
    marker=ROOT/'window/scientific.done.exit';tmp=marker.with_suffix('.tmp')
    tmp.write_text('0\n' if results['status']=='COMPLETED' else '1\n');tmp.replace(marker)
    return 0

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    else:sys.exit(run())
