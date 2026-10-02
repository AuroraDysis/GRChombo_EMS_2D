#!/usr/bin/env python3
"""Read completed receipts only; never calls the running diagnosis worker."""
import sys
sys.dont_write_bytecode=True
import csv, json, math, statistics
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t17')
DIAG=ROOT/'finder-diagnosis'

def parameters(path):
    return {s.split('=',1)[0].strip():s.split('=',1)[1].split('#',1)[0].strip()
            for s in path.read_text().splitlines() if '=' in s and not s.lstrip().startswith('#')}

def read(path):
    with path.open() as f: return list(csv.DictReader(f))

def write(name,data):
    assert data
    with (HERE/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=data[0]);writer.writeheader();writer.writerows(data)

def inspect(d,origin):
    r=json.loads((d/(d.name+'.resources.json')).read_text())
    assert r['returncode'] in (0,1,2) and r['returncode']==r['child_measurement']['returncode']
    assert int((d/'done.exit').read_text())==r['returncode']%256
    assert not r['gate_reason'] and r['peak_rss_bytes']<6e9
    assert r['tree_peak_bytes']<6.5e9
    assert r['tree_peak_bytes']>0 or r.get('recovered_accounting',False)
    p=parameters(d/'params.txt'); surfaces=read(d/'surfaces.csv');count=int(p['RH_num_horizons'])
    assert surfaces and all(math.isfinite(float(x['expansion_squared'])) for x in surfaces)
    final=surfaces[-count:]; strict=[x for x in surfaces if int(x['stage'])==2]
    passed=r['returncode']==0 and len(strict)==count and all(x['status']=='FOUND' and float(x['expansion_squared'])<=1e-12 for x in strict)
    assert (r['returncode']==0)==passed
    # Streaming witness of minima and late behaviour; memory independent of update count.
    progress={}
    with (d/'progress.csv').open() as f:
        for x in csv.DictReader(f):
            key=int(x['stage']),int(x['search_index']); value=float(x['expansion_squared'])
            assert math.isfinite(value)
            v=progress.setdefault(key,dict(minimum=math.inf,minimum_update=0,tail=deque(maxlen=128)))
            if value<v['minimum']: v['minimum'],v['minimum_update']=value,int(x['update'])
            v['tail'].append(value)
    group=surfaces[-count:]
    data=[]
    for x in group:
        k=int(x['stage']),int(x['search_index']);v=progress.get(k)
        seed=next(y for y in surfaces if int(y['stage'])==-1 and y['search_index']==x['search_index'])
        shape=[float(y) for y in (d/f"shape-{x['search_index']}-{x['stage']}.dat").read_text().split()[2:]]
        tail=list(v['tail']) if v else [float(x['expansion_squared'])]
        data.append(dict(origin=origin,case=d.name,rung=d.name.split('-')[0],kind=d.name.split('-')[1] if origin=='diagnosis' else 'binary',
            time=float(x['time']),N=int(x['N_theta']),hole=int(x['search_index']),chase=p['RH_chase_speeds'].split()[0],
            quota=p['RH_time_step_freq'].split()[0],cap=int(p['offline_max_updates']),floor=int(p['offline_floor_window']),
            wall_limit=float(p['offline_seconds']),seed_radius=float(p['RH_initial_radii'].split()[0]),
            seed_centre=float(p['RH_initial_centre'].split()[int(x['search_index'])]),strict_PASS=passed,
            status=x['status'],stage=int(x['stage']),updates=int(x['updates']),expansion_squared=float(x['expansion_squared']),
            expansion_RMS=math.sqrt(float(x['expansion_squared'])),minimum_in_final_stage=v['minimum'] if v else '',
            minimum_update=v['minimum_update'] if v else '',tail_min=min(tail),tail_max=max(tail),
            tail_relative_span=(max(tail)-min(tail))/max(tail) if max(tail) else 0,
            final_centre=float(x['centre']),centre_change=float(x['centre'])-float(seed['centre']),
            A=float(x['A']),Q=float(x['Q']),radius_min=min(shape),radius_max=max(shape),
            hits_author_radial_bound=any(y<=.0001 or y>=10 for y in shape),
            wall_seconds=r['wall_seconds'],finder_seconds=float(x['seconds']),peak_RSS_bytes=r['peak_rss_bytes'],
            gate_reason=r['gate_reason'],child_returncode=r['returncode'],directory=str(d)))
    return data

def main():
    now=datetime.now(timezone.utc)
    completed=sorted(d for d in (DIAG/'probes').iterdir() if (d/'done.exit').exists())
    original=sorted(d for d in (ROOT/'finder').iterdir() if (d/'done.exit').exists())
    data=[x for d in completed for x in inspect(d,'diagnosis')]
    data += [x for d in original for x in inspect(d,'original')]
    write('t17-finder-interim.csv',data)
    pending=[]
    for d in sorted((DIAG/'probes').iterdir()):
        if not (d/'done.exit').exists():
            pending.append(dict(case=d.name,start_UTC=datetime.fromtimestamp((d/'params.txt').stat().st_mtime,timezone.utc).isoformat(),
                                elapsed_seconds=now.timestamp()-(d/'params.txt').stat().st_mtime))
    registration=json.loads((DIAG/'registration.json').read_text())
    assert all(Path(x['path']).stat().st_size==x['bytes'] for x in registration['pins'])
    result=dict(snapshot_UTC=now.isoformat(),completed_diagnosis_probes=len(completed),completed_original_probes=len(original),
        active=pending,strict_pass_cases=sorted({x['case'] for x in data if x['strict_PASS']}),
        maximum_completed_RSS_bytes=max(x['peak_RSS_bytes'] for x in data),
        pipeline_done=(DIAG/'pipeline/done.exit').exists(),sources=[str(d) for d in completed])
    (HERE/'t17-finder-interim.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2),flush=True)
    for kind in ('unboosted','boosted','binary'):
        for speed,quota in (('1','1'),('32','1'),('128','1'),('512','1'),('128','16')):
            sub=[x for x in data if x['origin']=='diagnosis' and x['kind']==kind and x['rung']=='coarse' and
                 x['chase']==speed and x['quota']==quota and x['cap']==100000 and x['time']==0]
            for n in (48,96):
                rr=[x for x in sub if x['N']==n]
                if not rr:continue
                print(kind,speed,quota,n,[(x['case'],x['status'],x['stage'],x['expansion_squared'],
                    x['minimum_in_final_stage'],x['tail_relative_span'],x['finder_seconds']) for x in rr if x['hole']==0],flush=True)
if __name__=='__main__':main()
