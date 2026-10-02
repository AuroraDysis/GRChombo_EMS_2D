#!/usr/bin/env python3
"""Final receipt/stage/pair audit; no initializer, evolution or finder invocation."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, json, math, re, statistics
from datetime import datetime, timezone
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('interim',HERE/'t17-finder-interim.py')
i=importlib.util.module_from_spec(spec);spec.loader.exec_module(i)
ROOT=i.ROOT; DIAG=i.DIAG
SCHEDULE=(1e-7,1e-10,1e-12)

def pipeline_receipt(d,label):
    r=json.loads((d/(label+'.resources.json')).read_text())
    assert int((d/'done.exit').read_text())==r['returncode']==r['child_measurement']['returncode']==0
    assert not r['gate_reason'] and r['peak_rss_bytes']<6e9 and r['tree_peak_bytes']<6.5e9
    return r

def extrapolate(d,n,known48):
    p=[x for x in i.read(d/'progress.csv') if int(x['stage'])==1]
    raw=i.read(d/'surfaces.csv'); end=next(x for x in raw if int(x['stage'])==1)
    first=next(x for x in raw if int(x['stage'])==0)
    period=(float(end['seconds'])-float(first['seconds']))/int(end['updates'])
    result=[]
    for fraction in (.5,.25,.1):
        tail=p[int(len(p)*(1-fraction)):]
        fit=statistics.linear_regression([int(x['update']) for x in tail],
                                        [math.log(float(x['expansion_squared'])) for x in tail])
        assert fit.slope<0
        forecast=float(end['seconds'])+(math.log(1e-12)-math.log(float(end['expansion_squared'])))/fit.slope*period
        result.append(dict(N=n,fit_stage=1,tail_fraction=fraction,log_residual_slope_per_update=fit.slope,
            observed_updates=int(end['updates']),observed_seconds=float(end['seconds']),
            observed_expansion_squared=float(end['expansion_squared']),seconds_per_update=period,
            projected_finish_seconds=forecast,local_pair_seconds=known48+forecast if n==96 else '',
            actual_finish_seconds=known48 if n==48 else '',
            relative_prediction_error=(forecast-known48)/known48 if n==48 else '',
            measured_strict_pair=False,scope='local frozen unboosted coarse; constant late log slope is an unverified extrapolation'))
    return result

def main():
    for x in json.loads((DIAG/'registration.json').read_text())['pins']:
        p=Path(x['path']); assert p.stat().st_size==x['bytes']
        with p.open('rb') as f: assert hashlib.file_digest(f,'sha256').hexdigest()==x['sha256'],p
    pipes=[pipeline_receipt(DIAG/'pipeline','queue'),pipeline_receipt(DIAG/'pipeline/job','finder-diagnosis')]
    controls=[]
    for d in sorted((DIAG/'initialization').iterdir()):
        r=json.loads((d/(d.name+'.resources.json')).read_text())
        assert int((d/'done.exit').read_text())==r['returncode']==r['child_measurement']['returncode']==0
        assert not r['gate_reason'] and r['peak_rss_bytes']<6e9 and r['tree_peak_bytes']<6.5e9
        assert 'GRChombo finished.' in (d/'run.log').read_text()
        controls.append(dict(case=d.name,peak_RSS_bytes=r['peak_rss_bytes'],wall_seconds=r['wall_seconds'],result='verified initialization'))
    assert len(controls)==4
    dirs=[(d,'original') for d in sorted((ROOT/'finder').iterdir()) if d.is_dir()]
    dirs += [(d,'diagnosis') for d in sorted((DIAG/'probes').iterdir()) if d.is_dir()]
    assert len(dirs)==66 and sum(origin=='diagnosis' for _,origin in dirs)==58
    stages=[]; summary=[]; evidence=[]; by_pair={}
    for index,(d,origin) in enumerate(dirs,1):
        assert (d/'done.exit').exists(),d
        final=i.inspect(d,origin) # Own child, marker, gates, full finite trace and final shape witnesses.
        p=i.parameters(d/'params.txt'); raw=i.read(d/'surfaces.csv'); count=int(p['RH_num_horizons'])
        assert all(float(x['time']) in (0.,.4375) for x in raw)
        ident=('O' if origin=='original' else 'D')+f'{index:02d}'
        seed=re.search(r'-r([\d.]+)-o(-?[\d.]+)-wall',d.name)
        radius,offset=(float(seed[1]),float(seed[2])) if seed else (1.,0.)
        passed=all(x['strict_PASS'] for x in final)
        for hole in range(count):
            attempted=[x for x in raw if int(x['search_index'])==hole and int(x['stage'])>=0]
            assert attempted
            assert [int(x['stage']) for x in attempted]==list(range(len(attempted)))
            assert len(attempted)==3 if passed else len(attempted)<=3
            for x in attempted:
                stage=int(x['stage']); assert float(x['threshold'])==SCHEDULE[stage]
                assert (x['status']=='FOUND')==(float(x['expansion_squared'])<=SCHEDULE[stage])
            if passed: assert len(attempted)==3 and all(x['status']=='FOUND' for x in attempted)
        for x in raw:
            stages.append(dict(probe_id=ident,origin=origin,case=d.name,rung=final[0]['rung'],kind=final[0]['kind'],
                radius_factor=radius,center_offset_Rh=offset,chase=p['RH_chase_speeds'].split()[0],quota=p['RH_time_step_freq'].split()[0],
                cap=int(p['offline_max_updates']),floor_window=int(p['offline_floor_window']),wall_limit=float(p['offline_seconds']),
                strict_probe_PASS=passed,wall_seconds=final[0]['wall_seconds'],peak_RSS_bytes=final[0]['peak_RSS_bytes'],
                child_returncode=final[0]['child_returncode'],directory=str(d),**x,
                expansion_RMS=math.sqrt(float(x['expansion_squared']))))
        for x in final: summary.append(dict(probe_id=ident,**x))
        r=json.loads((d/(d.name+'.resources.json')).read_text())
        evidence.append(dict(probe_id=ident,case=d.name,origin=origin,rung=final[0]['rung'],kind=final[0]['kind'],
            time_M=final[0]['time'],N_theta=final[0]['N'],chase=final[0]['chase'],quota=final[0]['quota'],
            radius_factor=radius,offset_Rh=offset,cap=final[0]['cap'],floor=final[0]['floor'],wall_limit=final[0]['wall_limit'],
            E0=max(float(x['expansion_squared']) for x in raw if int(x['stage'])==0),
            E1=max((float(x['expansion_squared']) for x in raw if int(x['stage'])==1),default='NOT_ATTEMPTED'),
            E2=max((float(x['expansion_squared']) for x in raw if int(x['stage'])==2),default='NOT_ATTEMPTED'),
            status='/'.join(sorted({x['status'] for x in final})),last_stage=final[0]['stage'],
            last_stage_updates=final[0]['updates'],total_updates=sum(int(x['updates']) for x in raw if int(x['search_index'])==0 and int(x['stage'])>=0),
            finder_seconds=max(float(x['seconds']) for x in raw),wall_seconds=r['wall_seconds'],
            peak_RSS_bytes=r['peak_rss_bytes'],child_returncode=r['returncode'],gate_reason=r['gate_reason'],strict_PASS=passed,
            recovered_accounting=r.get('recovered_accounting',False),tree_peak_bytes=r['tree_peak_bytes'],directory=str(d)))
        key=(origin,final[0]['rung'],final[0]['kind'],final[0]['time'],final[0]['chase'],final[0]['quota'],
             final[0]['cap'],final[0]['floor'],final[0]['wall_limit'],radius,offset)
        assert final[0]['N'] not in by_pair.setdefault(key,{})
        by_pair[key][final[0]['N']]=final
    pairs=[]
    for key,ns in by_pair.items():
        for hole in range(len(next(iter(ns.values())))):
            a=next((x for x in ns.get(48,[]) if x['hole']==hole),None)
            b=next((x for x in ns.get(96,[]) if x['hole']==hole),None)
            strict=bool(a and b and a['strict_PASS'] and b['strict_PASS'])
            da=abs(a['A']-b['A'])/abs(b['A']) if a and b else ''
            dq=abs(a['Q']-b['Q'])/abs(b['Q']) if a and b else ''
            pairs.append(dict(origin=key[0],rung=key[1],kind=key[2],time_M=key[3],chase=key[4],quota=key[5],cap=key[6],floor=key[7],
                wall_limit=key[8],radius_factor=key[9],offset_Rh=key[10],hole=hole,
                N48_case=a['case'] if a else 'NOT_RUN',N96_case=b['case'] if b else 'NOT_RUN',
                N48_squared=a['expansion_squared'] if a else '',N96_squared=b['expansion_squared'] if b else '',
                both_strict_residual_PASS=strict,relative_A_of_stopped_outputs=da,relative_Q_of_stopped_outputs=dq,
                angular_qualified=bool(strict and da<=1e-3 and dq<=1e-3),
                warning='' if strict else 'stopped A/Q comparisons cannot qualify horizons'))
    assert not any(x['kind'] in ('boosted','binary') and x['angular_qualified'] for x in pairs)
    passes=[x for x in evidence if x['strict_PASS']]
    assert len(passes)==1 and passes[0]['kind']=='unboosted' and passes[0]['N_theta']==48
    for name,data in [('t17-finder-final-stages.csv',stages),('t17-finder-final-probes.csv',evidence),
                      ('t17-finder-final-residuals.csv',summary),('t17-finder-final-pairs.csv',pairs),
                      ('t17-finder-final-controls.csv',controls)]: i.write(name,data)
    known=passes[0]['finder_seconds']; d48=Path(passes[0]['directory'])
    d96=DIAG/'probes/coarse-unboosted-t0-n96-s1-q1-c100000-f100001-r1-o0-wall235'
    cost=extrapolate(d48,48,known)+extrapolate(d96,96,known)
    i.write('t17-finder-final-cost.csv',cost)
    result=dict(verified_UTC=datetime.now(timezone.utc).isoformat(),diagnosis_probes=58,original_probes=8,
        initialization_controls=4,all_receipts_verified=True,pipeline_receipts=[dict(process=x['process'],returncode=x['returncode'],
            peak_RSS_bytes=x['peak_rss_bytes'],wall_seconds=x['wall_seconds']) for x in pipes],
        max_probe_RSS_bytes=max(x['peak_RSS_bytes'] for x in evidence),
        actual_probe_returncodes={str(k):sum(x['child_returncode']==k for x in evidence) for k in (0,1,2)},
        strict_probe_cases=[x['case'] for x in passes],strict_boosted_pairs=0,strict_binary_pairs=0,
        no_grid_admitted=True,policy_unchanged=True,new_probes_started=False,
        measured_local_unboosted_N48_seconds=known,measured_strict_N48_N96_pair_available=False,
        local_pair_forecast_seconds=[x['local_pair_seconds'] for x in cost if x['N']==96],
        model_warning='late constant log-slope extrapolation only; no HPC scaling or rigorous bound',
        total_T17_output_bytes=sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    assert result['total_T17_output_bytes']<=40e9
    (HERE/'t17-finder-final.json').write_text(json.dumps(result,indent=2)+'\n')
    table=['| ID | Data/rung/t/N | Chase/quota | Seed R/Rh, offset/Rh | Cap, floor, seconds | E0 | E1 | E2 | Stop, updates | Finder / process seconds | RSS GB |',
           '|---|---|---|---|---|---:|---:|---:|---|---:|---:|']
    def number(v): return '—' if isinstance(v,str) else f'{v:.6g}'
    for x in evidence:
        table.append(f"| {x['probe_id']} | {x['kind']}/{x['rung']}/{x['time_M']:g}/{x['N_theta']} | {x['chase']}/{x['quota']} | {x['radius_factor']:g}, {x['offset_Rh']:g} | {x['cap']}, {x['floor']}, {x['wall_limit']:g} | {number(x['E0'])} | {number(x['E1'])} | {number(x['E2'])} | {x['status']}, {x['last_stage_updates']} | {x['finder_seconds']:.2f} / {x['wall_seconds']:.2f} | {x['peak_RSS_bytes']/1e9:.4f} |")
    report=HERE/'t17-finder-final.md'; text=report.read_text(); marker='<!-- ALL-PROBES -->';assert marker in text
    report.write_text(text[:text.index(marker)]+marker+'\n\n'+'\n'.join(table)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
