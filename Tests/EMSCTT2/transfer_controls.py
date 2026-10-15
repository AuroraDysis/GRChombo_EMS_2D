"""T2: unchanged t=0 layouts/masks, legacy and production point transfer."""
from pathlib import Path
import hashlib,json,os,subprocess,resource,time,csv,math
R=Path(__file__).resolve().parent
OUT=R/'results/point-transfer';OUT.mkdir(exist_ok=True)
ENV={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}

def authenticate():
    old=json.loads((R/'results/authentication.json').read_text())
    auth=[]
    for name in ('a0.9-e7.trumpet','n24-final.ctt'):
        row=next(r for r in old if Path(r['path']).name==name)
        path=R/'results/input-snapshots'/name
        b=path.read_bytes();h=hashlib.sha256(b).hexdigest()
        assert h==row['sha256'],name
        if name.endswith('.ctt'):
            payload,trailer=b.rsplit(b'SHA256 ',1)
            assert hashlib.sha256(payload).hexdigest()==trailer.decode().strip()
        auth.append(dict(path=str(path),sha256=h,bytes=len(b),T1_sha256=row['sha256']))
    for path in [R/'EMSCTT2Hierarchy.cpp',R/'transfer_controls.py',next(R.glob('EMSCTT2Hierarchy*.ex')),R.parents[1]/'Source/GRChomboCore/PointAMRTransfer.hpp',R.parents[1]/'Source/GRChomboCore/SimulationParametersBase.hpp',R.parents[1]/'Source/GRChomboCore/GRAMRLevel.cpp']+list(R.glob('params-*.txt')):
        b=path.read_bytes();auth.append(dict(path=str(path),sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)))
    (OUT/'authentication.json').write_text(json.dumps(auth,indent=2)+'\n')

def run():
    authenticate();exe=str(next(R.glob('EMSCTT2Hierarchy*.ex')))
    records=[]
    for domain in ('global','near'):
        stem='' if domain=='global' else 'near-'
        for mode,ghosts,analytic in [('legacy',5,False),('legacy3',3,False),('point',3,False),('analytic',3,True),('analytic5',5,True)]:
            for n in (80,160,320):
                args=[exe,str(R/f'params-{stem}{n}.txt'),
                      'amr_transfer='+('point' if mode in ('point','analytic') else 'legacy'),
                      f'num_ghosts={ghosts}',f'ctt2_test_analytic_ghosts={int(analytic)}',
                      'ems_data_path='+str(R/'results/input-snapshots/a0.9-e7.trumpet'),
                      'ems_ctt_data_path='+str(R/'results/input-snapshots/n24-final.ctt')]
                start=time.monotonic();name=f'{domain}-{mode}-{n}'
                with (OUT/(name+'.log')).open('w') as log:
                    p=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,timeout=240,env=ENV)
                rec=dict(name=name,command=args,seconds=time.monotonic()-start,exit=p.returncode,peak_child_rss_upper_bound=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
                records.append(rec);(OUT/'cost.json').write_text(json.dumps(records,indent=2)+'\n')
                assert p.returncode==0,rec
                assert rec['peak_child_rss_upper_bound']<6_000_000_000,rec
                print(name,rec['seconds'],flush=True)

def report():
    rows=[];layouts={};counts={}
    for domain in ('global','near'):
        for mode in ('legacy','legacy3','point','analytic','analytic5'):
            prev={}
            for n in (80,160,320):
                text=(OUT/f'{domain}-{mode}-{n}.log').read_text()
                assert 'seconds=' in text
                layout=sorted(l for l in text.splitlines() if l.startswith(('BOX,','LEVEL,')))
                key=(domain,n)
                if key in layouts:assert layouts[key]==layout,(domain,mode,n,'changed layout')
                else:layouts[key]=layout
                census={}
                for line in text.splitlines():
                    if not line.startswith('NORM,'):continue
                    a=line.split(',');region=a[1];count=int(a[2]);floors=int(a[3]);staging=float(a[4]);v=list(map(float,a[5:]))
                    assert count>0 and floors==0 and staging<5e-13
                    assert len(v)==8 and all(math.isfinite(x) for x in v)
                    census[region]=count
                    for j,q in enumerate(('H','M','GE','GB')):
                        for statistic,offset in [('RMS',0),('max',1)]:
                            value=v[2*j+offset];pk=(region,q,statistic)
                            before=prev.get(pk);order=math.log2(before/value) if before and value else None
                            rows.append(dict(domain=domain,mode=mode,N1=n,region=region,diagnostic=q,statistic=statistic,cells=count,staging=staging,value=value,order=order))
                            prev[pk]=value
                if key in counts:assert counts[key]==census
                else:counts[key]=census
    assert len(rows)==960,len(rows)
    baseline=[]
    for domain in ('global','near'):
        for n in (80,160,320):
            stem='' if domain=='global' else 'near-'
            old={l.split(',')[1]:list(map(float,l.split(',')[5:])) for l in (R/f'results/hierarchy-{stem}{n}.log').read_text().splitlines() if l.startswith('NORM,')}
            for region,values in old.items():
                current=[r for r in rows if r['domain']==domain and r['mode']=='legacy' and r['N1']==n and r['region']==region]
                diffs=[abs(r['value']-values[2*('H','M','GE','GB').index(r['diagnostic'])+(r['statistic']=='max')]) for r in current]
                baseline.append(dict(domain=domain,N1=n,region=region,max_abs_difference=max(diffs)))
    (OUT/'T1-comparison.json').write_text(json.dumps(baseline,indent=2)+'\n')
    assert not subprocess.check_output(['git','diff','98778e5','--','Source','Examples']).strip(), 'production source changed'
    costs=json.loads((OUT/'cost.json').read_text())
    assert len(costs)==30 and len({r['name'] for r in costs})==30
    assert all(r['exit']==0 and r['peak_child_rss_upper_bound']<6_000_000_000 for r in costs)
    point=[r for r in rows if r['mode']=='point' and r['diagnostic']!='GB' and r['N1']!=80]
    assert all(r['order']>0 for r in point), 'nondecreasing point-transfer norm requires localization'
    assert all(r['value']==0 and r['order'] is None for r in rows if r['diagnostic']=='GB')
    (OUT/'decision.json').write_text(json.dumps(dict(production_source_unchanged_from_T1=True,all_point_RMS_and_maxima_decrease=True,positive_pair_orders=len(point),magnetic_Gauss_exactly_zero=True,minimum_point_pair_order=min(r['order'] for r in point),max_abs_T1_legacy_replay_difference=max(r['max_abs_difference'] for r in baseline)),indent=2)+'\n')
    with (OUT/'norms.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (OUT/'layout-check.json').write_text(json.dumps({'identical_boxes_and_valid_cell_counts':True,'cases':len(layouts)},indent=2)+'\n')
    lines=[]
    for domain in ('global','near'):
        for region in ('inner','pulse','outer','pulse_axis'):
            lines += [f'## {domain}: {region}','', '| Mode | N1 | Diagnostic | RMS | order | max | order |','|---|---:|---|---:|---:|---:|---:|']
            for mode in ('legacy','legacy3','point','analytic','analytic5'):
                for n in (80,160,320):
                    for q in ('H','M','GE','GB'):
                        r=[x for x in rows if (x['domain'],x['region'],x['mode'],x['N1'],x['diagnostic'])==(domain,region,mode,n,q)]
                        assert len(r)==2
                        fmt=lambda x:f'{x:.6e}' if x is not None else '—'
                        lines.append(f'| {mode} | {n} | {q} | '+ ' | '.join(fmt(x) for x in [r[0]['value'],r[0]['order'],r[1]['value'],r[1]['order']])+' |')
            lines.append('')
    (OUT/'tables.md').write_text('\n'.join(lines).rstrip()+'\n')
    for domain in ('global','near'):
        for mode in ('legacy','legacy3','point','analytic','analytic5'):
            print(domain,mode,[(x['diagnostic'],x['value'],x['order']) for x in rows if x['domain']==domain and x['mode']==mode and x['N1']==320 and x['region']=='outer' and x['statistic']=='RMS'])

if __name__=='__main__':
    import sys
    if len(sys.argv)==1:run()
    report()
