#!/usr/bin/env python3
"""Receipt audit and continuation after the T17 output gate; no completed evolution repeats."""
import sys
sys.dont_write_bytecode=True
import csv, importlib.util, json, math, re
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('work',HERE/'t17-work.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
ROOT=w.ROOT
def read_rows(path):
    return list(csv.DictReader(path.open()))
def checked(d,allowed=(0,)):
    r=json.loads((d/(d.name+'.resources.json')).read_text())
    assert r['returncode'] in allowed and r['child_measurement']['returncode']==r['returncode']
    assert not r['gate_reason'] and r['peak_rss_bytes']<6e9
    assert int((d/'done.exit').read_text())==(r['returncode']%256)
    return r
def audit():
    gate=json.loads((ROOT/'pipeline/job/local-evidence.resources.json').read_text())
    assert gate['returncode']==-15 and gate['gate_reason']=='output ceiling 16000000000 bytes'
    assert int((ROOT/'pipeline/done.exit').read_text())==125
    # wait4's child result is complete, although the outer gate interrupted the
    # watchdog's final accounting. Recover only accounting, never finder state.
    d=ROOT/'finder/fine-step0-n48';target=d/(d.name+'.resources.json')
    if not target.exists():
        child=json.loads((d/(d.name+'.child.json')).read_text());assert child['returncode']==1
        surface=read_rows(d/'surfaces.csv');last=[x for x in surface if int(x['stage'])==0]
        assert len(last)==2 and all(x['status']=='UPDATE_CAP' and int(x['updates'])==2000 for x in last)
        assert 'time: sysctl kern.clockrate' in (d/'run.log').read_text()
        wall=float(re.search(r'([\d.]+) real',(d/(d.name+'.time')).read_text()).group(1))
        recovered=dict(process=d.name,returncode=1,child_measurement=child,
            peak_rss_bytes=child['peak_rss_bytes'],rss_GB=child['peak_rss_bytes']/1e9,
            cap_GB=6,gate_reason='',wall_seconds=wall,tree_peak_bytes=0,
            tree_measurement='not retained; outer gate interrupted final watchdog receipt',
            recovered_accounting=True,measurement='own completed wait4 child receipt and time -l; no probe repeated')
        target.write_text(json.dumps(recovered,indent=2)+'\n');(d/'done.exit').write_text('1\n')
    data=[]
    for category in ('census','evolution','finder','timing'):
        for d in sorted((ROOT/category).iterdir()):
            if not d.is_dir():continue
            r=checked(d,(0,-11) if category=='census' else (0,1) if category=='finder' else (0,))
            status='COMPLETED'
            detail=''
            if category=='census':
                if r['returncode']==-11:
                    status='NATIVE_SIGSEGV';detail='fine native attempt failed; registered dense retry used'
                else:assert 'GRChombo finished.' in (d/'run.log').read_text()
            if category=='evolution':
                log=(d/'run.log').read_text();assert 'T13 clean native stop' in log and 'GRChombo finished.' in log
                assert 'Read EMSTRUMPET' not in log
                expected=.00197601318359375 if d.name.startswith('fine') else .001922607421875
                assert float(read_rows(d/'t13-stop.csv')[0]['actual_time_M'])==expected
                floor=[x for path in d.glob('t13-floors-*.csv') for x in read_rows(path)];assert floor
                assert all(int(x[k])==0 for x in floor for k in ('chi_activations','lapse_activations','nonfinite'))
                detail='native stop, static guard/restart and zero launch floors verified'
            if category=='finder':
                last=read_rows(d/'surfaces.csv')[-2:];assert len(last)==2
                status='FINDER_'+last[0]['status'];detail=' '.join(x['expansion_squared'] for x in last)
                assert all(math.isfinite(float(x['expansion_squared'])) for x in last)
            if category=='timing':
                rate=read_rows(d/'t17-coarse-rates.csv')
                assert [float(x['time_M']) for x in rate]==[.4375,.875]
                assert all(float(x['seconds'])>0 and float(x['dt0_M'])==.4375 for x in rate)
                assert 'GRChombo finished.' in (d/'run.log').read_text()
                detail=' '.join(x['seconds'] for x in rate)+' s; includes unintended capture overhead'
            data.append(dict(category=category,case=d.name,status=status,actual_child_returncode=r['returncode'],
                peak_RSS_bytes=r['peak_rss_bytes'],wall_seconds=r['wall_seconds'],
                recovered_accounting=r.get('recovered_accounting',False),action='KEEP; no repeat',detail=detail))
    for case in ('fine-step0-n96','fine-step1-n48','fine-step1-n96'):
        assert not (ROOT/'finder'/case/'done.exit').exists()
        data.append(dict(category='finder',case=case,status='UNSTARTED',actual_child_returncode='',
            peak_RSS_bytes='',wall_seconds='',recovered_accounting=False,action='RUN',detail='registered probe unchanged'))
    assert not (ROOT/'timing/fine-single').exists()
    data.append(dict(category='timing',case='fine-single',status='UNSTARTED',actual_child_returncode='',
        peak_RSS_bytes='',wall_seconds='',recovered_accounting=False,action='RUN',detail='fine reference if no rung qualifies; rates and box metadata only'))
    for name in ('t17-initial-identity.csv','t17-dense-identity.csv'):
        values=read_rows(HERE/name);assert values and all(int(x['bit_mismatches'])==0 for x in values)
    endpoints=read_rows(HERE/'t17-launch-endpoint.csv');required=[x for x in endpoints if x['field']=='Gamma']
    assert len(required)==16 and all(x['qualified']=='True' and float(x['upper_5x'])<1 for x in required)
    w.save('t17-receipt-audit.csv',data)
    print('PASS receipt audit: keep all completed census/launch/binary timing and five completed finder cases.')
def continue_work():
    for item in json.loads((HERE/'t17-inputs.json').read_text()):
        assert w.digest(Path(item['path']))==item['sha256']
    # horizons() skips cases with exact matching parameters and retained marker.
    # Nothing calls census(), launch(), controls() or timed_binary() here.
    checked(ROOT/'finder/fine-step0-n48',(1,))
    w.horizons('fine')
    analysis=w.module('analysis',HERE/'t17-analyze.py');choice=analysis.resolution()
    assert choice!='coarse','Completed coarse probes failed; coarse admission would be inconsistent'
    w.timed_single(json.loads((ROOT/'choices.json').read_text()),choice or 'fine')
    checked(ROOT/'timing/fine-single')
    analysis.finish()
    print('T17 unfinished evidence complete; inspect admission verdict. No merger run.')
def rate_check():
    text=w.p.replace(w.base(1,56,64),dict(max_steps=2,stop_time=1,checkpoint_interval=1,plot_interval=1,
        ems_use_geometric_initial_lapse='true',verbosity=0))
    d=ROOT/'rate-smoke/identity';assert w.run_case(d.name,d,ROOT/'launch-rate.ex',text)==0
    count,bits=w.p.compare(ROOT/'controls/geometric-capture/chk/EMS_000002.2d.hdf5',d/'chk/EMS_000002.2d.hdf5')
    assert count>0 and bits==0 and not list(d.glob('t13-*.xz'))
    d=ROOT/'rate-smoke/output-free';text=w.p.replace(text,dict(checkpoint_interval=-1,plot_interval=-1))
    assert w.run_case(d.name,d,ROOT/'launch-rate.ex',text)==0
    assert len(read_rows(d/'t17-coarse-rates.csv'))==2
    assert not list(d.rglob('*.hdf5')) and not list(d.glob('t13-*.xz'))
    w.save('t17-rate-check.csv',[dict(check='disabled capture guard preserves native final fields',Float64_values=count,bit_mismatches=bits,
        output_free_steps=2,field_output_files=0,result='PASS')])
    print('PASS disabled capture bit identity and two output-free native step timers.')
if __name__=='__main__':
    {'audit':audit,'continue':continue_work,'rate-check':rate_check}[sys.argv[1]]()
