"""Production point-transfer t=0 controls, including the four omitted faces.
T1/T2 legacy evidence remains frozen in results/point-transfer and the archive.
"""
from pathlib import Path
import hashlib, json, os, subprocess, resource, time, csv, math
from evidence import R, RAW, OUT, ARCHIVE, reference
ENV = {**os.environ, 'OMP_NUM_THREADS':'1', 'OPENBLAS_NUM_THREADS':'1'}
DATA = R/'results/input-snapshots'

def authenticate():
    old = json.loads((R/'results/authentication.json').read_text())
    rows = []
    for name in ('a0.9-e7.trumpet','n24-final.ctt'):
        path = DATA/name
        b = path.read_bytes(); h = hashlib.sha256(b).hexdigest()
        assert h == next(r['sha256'] for r in old if Path(r['path']).name == name)
        if name.endswith('.ctt'):
            payload, trailer = b.rsplit(b'SHA256 ',1)
            assert hashlib.sha256(payload).hexdigest() == trailer.decode().strip()
        rows.append(dict(path=str(path), sha256=h, bytes=len(b)))
    manifest = ARCHIVE/'manifest.json'
    receipt = json.loads((R/'results/raw-archive.json').read_text())
    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == receipt['manifest_sha256']
    for row in json.loads(manifest.read_text()):
        assert hashlib.sha256((ARCHIVE/row['path']).read_bytes()).hexdigest() == row['sha256']
    for path in [R/'EMSCTT2Hierarchy.cpp', R/'transfer_controls.py', next(R.glob('EMSCTT2Hierarchy*.ex')),
                 R.parents[1]/'Examples/EMS/SimulationParameters.hpp',
                 R.parents[1]/'Source/GRChomboCore/PointAMRTransfer.hpp',
                 R.parents[1]/'Source/InitialConditions/EMSBH/1D_SOL/EMSCTTPulse_read.hpp',
                 R.parents[1]/'Source/InitialConditions/EMSBH/EMSBH_trumpet_read.impl.hpp',
                 R/'EMSCTT2Readback.cpp', next(R.glob('EMSCTT2Readback*.ex')), R/'parser.py', R/'evidence.py'] + list(R.glob('params-*.txt')):
        b=path.read_bytes(); rows.append(dict(path=str(path),sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)))
    (OUT/'authentication.json').write_text(json.dumps(rows,indent=2)+'\n')

def startup():
    exe = str(next(R.glob('EMSCTT2Hierarchy*.ex')))
    base = [exe,str(R/'params-80.txt'),'ctt2_test_check_params_only=1',
            'ems_data_path='+str(DATA/'a0.9-e7.trumpet'),
            'ems_ctt_data_path='+str(DATA/'n24-final.ctt')]
    rows=[]
    cases=[('legacy3',['amr_transfer=legacy'],False),
           ('legacy5',['amr_transfer=legacy','num_ghosts=5'],False),
           ('point5',['num_ghosts=5'],False),('point2',['num_ghosts=2'],False),
           ('point3',[],True),
           ('no-companion-legacy',['amr_transfer=legacy'],True)]
    # Omission must select the unchanged legacy default, then refuse EMSCTT 2.
    params = RAW/'missing-transfer.txt'
    params.write_text('\n'.join(l for l in (R/'params-80.txt').read_text().splitlines() if not l.startswith('amr_transfer'))+'\n')
    cases.append(('default-legacy',[],False))
    plain_params=RAW/'no-companion.txt'
    plain_params.write_text('\n'.join(l for l in (R/'params-80.txt').read_text().splitlines() if not l.startswith('ems_ctt_data_path'))+'\n')
    crlf = RAW/'CRLF.ctt'
    b=(DATA/'n24-final.ctt').read_bytes();payload=b[:b.rindex(b'SHA256 ')].replace(b'EMSCTT 2\n',b'EMSCTT 2\r\n',1)
    crlf.write_bytes(payload+b'SHA256 '+hashlib.sha256(payload).hexdigest().encode()+b'\n')
    cases.append(('CRLF',['ems_ctt_data_path='+str(crlf)],False))
    for name, args, accepted in cases:
        command=base.copy()
        if name=='default-legacy': command[1]=str(params)
        if name=='no-companion-legacy':
            command[1]=str(plain_params)
            command=[a for a in command if not a.startswith('ems_ctt_data_path=')]
        command+=args
        p=subprocess.run(command,env=ENV,capture_output=True,text=True,timeout=10)
        output=p.stdout+p.stderr
        (RAW/f'startup-{name}.log').write_text(output)
        assert (p.returncode==0)==accepted,(name,p.returncode,output)
        assert 'EMSTRUMPET profile:' not in output and 'LEVEL,' not in output
        expected='startup parameters accepted' if accepted else ('version line contains CR; use LF line endings' if name=='CRLF' else 'num_ghosts = 2 is invalid' if name=='point2' else 'EMSCTT 2 requires amr_transfer = point and num_ghosts = 3')
        assert expected in output,(name,output)
        rows.append(dict(name=name,exit=p.returncode,expected=expected,initialization_reached=False))
    # Exercise the dispatcher's CRLF diagnostic independently of startup.
    exe=str(next(R.glob('EMSCTT2Readback*.ex')))
    p=subprocess.run([exe,str(DATA/'a0.9-e7.trumpet'),str(crlf)],env=ENV,capture_output=True,text=True,timeout=10)
    assert p.returncode!=0 and 'version line contains CR; use LF line endings' in p.stdout+p.stderr
    for name in ('zero','n24'):
        p=subprocess.run([exe,str(DATA/'a0.9-e7.trumpet'),str(DATA/f'{name}-final.ctt'),'--adm-consistency'],input=reference('points.txt').read_text(),env=ENV,capture_output=True,text=True,timeout=10)
        assert p.returncode==0,(name,p.stderr)
        assert sum(len(l.split())==29 for l in p.stdout.splitlines())==146
    (OUT/'startup.json').write_text(json.dumps(dict(cases=rows,CRLF_dispatch_refused=True,ADM_cache_consistent_points_per_companion=146,ADM_cache_tolerance=5e-13),indent=2)+'\n')

def run():
    authenticate(); startup()
    exe=str(next(R.glob('EMSCTT2Hierarchy*.ex')))
    # Select the supplied panel faces in the two previously uncovered gaps.
    faces=[float(l.split()[3]) for l in (DATA/'n24-final.ctt').read_text().splitlines()
           if l.startswith('PANEL log ') and (8<float(l.split()[3])<25 or 55<float(l.split()[3])<68)]
    assert len(faces)==4,faces
    (OUT/'masks.json').write_text(json.dumps(dict(radii_M=faces,half_width_M=2.,source=str(DATA/'n24-final.ctt')),indent=2)+'\n')
    costs=[]
    for domain in ('global','near'):
        stem='' if domain=='global' else 'near-'
        for mode in ('point','analytic'):
            for n in (80,160,320):
                args=[exe,str(R/f'params-{stem}{n}.txt'),'ctt2_test_analytic_ghosts='+str(int(mode=='analytic')),
                      'ems_data_path='+str(DATA/'a0.9-e7.trumpet'),'ems_ctt_data_path='+str(DATA/'n24-final.ctt')]
                if domain=='global':args+=['ctt2_test_interfaces='+' '.join(map(str,faces)),'ctt2_test_interface_half_width=2']
                name=f'{domain}-{mode}-{n}'; start=time.monotonic()
                with (RAW/(name+'.log')).open('w') as log:
                    p=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,timeout=240,env=ENV)
                rec=dict(name=name,command=args,seconds=time.monotonic()-start,exit=p.returncode,peak_child_rss_upper_bound=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
                costs.append(rec);(OUT/'cost-transfer.json').write_text(json.dumps(costs,indent=2)+'\n')
                assert p.returncode==0 and rec['peak_child_rss_upper_bound']<6_000_000_000,rec
                print(name,rec['seconds'],flush=True)

def report():
    rows=[];layouts={};counts={};replay=[]
    old=list(csv.DictReader((R/'results/point-transfer/norms.csv').open()))
    for domain in ('global','near'):
        for mode in ('point','analytic'):
            previous={}
            for n in (80,160,320):
                text=(RAW/f'{domain}-{mode}-{n}.log').read_text()
                assert 'CONFIG,point,3,' in text and 'seconds=' in text
                layout=sorted(l for l in text.splitlines() if l.startswith(('BOX,','LEVEL,')))
                key=(domain,n)
                if key in layouts:assert layout==layouts[key]
                layouts[key]=layout
                historical=sorted(l for l in (ARCHIVE/f'point-transfer/{domain}-{mode}-{n}.log').read_text().splitlines() if l.startswith(('BOX,','LEVEL,')))
                assert layout==historical
                census={}
                for l in text.splitlines():
                    if not l.startswith('NORM,'):continue
                    a=l.split(',');region=a[1];count=int(a[2]);floor=int(a[3]);staging=float(a[4]);v=list(map(float,a[5:]))
                    assert count>0 and floor==0 and staging<5e-13 and len(v)==8 and all(math.isfinite(x) for x in v)
                    census[region]=count
                    for j,q in enumerate(('H','M','GE','GB')):
                        for statistic,offset in [('RMS',0),('max',1)]:
                            value=v[2*j+offset];pk=(region,q,statistic);before=previous.get(pk)
                            order=math.log2(before/value) if before and value else None
                            row=dict(domain=domain,mode=mode,N1=n,region=region,diagnostic=q,statistic=statistic,cells=count,staging=staging,value=value,order=order)
                            rows.append(row);previous[pk]=value
                            if not region.startswith('interface'):
                                prior=next(r for r in old if all(str(row[k])==r[k] for k in ('domain','mode','N1','region','diagnostic','statistic')))
                                err=abs(value-float(prior['value']))/(1+abs(float(prior['value'])))
                                assert err<=5e-13;replay.append(err)
                assert len(census)==(8 if domain=='global' else 4)
                if key in counts:assert counts[key]==census
                counts[key]=census
    with (OUT/'norms.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    lines=[]
    for domain in ('global','near'):
        for region in dict.fromkeys(r['region'] for r in rows if r['domain']==domain):
            lines += [f'## {domain}: {region}','','| Mode | N1 | Diagnostic | RMS | order | max | order |','|---|---:|---|---:|---:|---:|---:|']
            for mode in ('point','analytic'):
                for n in (80,160,320):
                    for q in ('H','M','GE','GB'):
                        pair=[r for r in rows if (r['domain'],r['region'],r['mode'],r['N1'],r['diagnostic'])==(domain,region,mode,n,q)]
                        fmt=lambda x:f'{x:.6e}' if x is not None else '—'
                        lines.append(f'| {mode} | {n} | {q} | '+' | '.join(fmt(x) for x in [pair[0]['value'],pair[0]['order'],pair[1]['value'],pair[1]['order']])+' |')
            lines.append('')
    (OUT/'tables.md').write_text('\n'.join(lines))
    interface_diffs=[]
    for row in rows:
        if row['mode']=='point' and row['region'].startswith('interface'):
            analytic=next(a for a in rows if a['mode']=='analytic' and all(a[k]==row[k] for k in ('domain','N1','region','diagnostic','statistic')))
            diff=abs(row['value']-analytic['value'])/(1+abs(analytic['value']))
            assert diff<=5e-13
            interface_diffs.append(diff)
    point=[r for r in rows if r['mode']=='point' and r['diagnostic']!='GB' and r['N1']!=80]
    assert all(r['value']==0 for r in rows if r['diagnostic']=='GB')
    (OUT/'decision.json').write_text(json.dumps(dict(identical_layouts_to_T2=True,interface_point_vs_analytic_max_scaled_difference=max(interface_diffs),original_mask_replay_max_scaled_difference=max(replay),replay_allowance=5e-13,nonpositive_point_orders=[r for r in point if r['order']<=0],minimum_point_pair_order=min(r['order'] for r in point),magnetic_Gauss_zero=True),indent=2)+'\n')

if __name__=='__main__':
    import sys
    if len(sys.argv)==1:run()
    report()
