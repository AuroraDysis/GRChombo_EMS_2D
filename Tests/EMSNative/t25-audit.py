#!/usr/bin/env python3
"""Seal T25 receipts, source hashes, angular comparisons and bracket products.

All big files are hashed by streaming. This never evaluates a static profile;
the single-hole header benchmark is provenance, not an input to the map.
"""
import csv,hashlib,importlib.util,json,math,subprocess
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t25')
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return list(csv.DictReader(Path(p).open()))
def write(p,rows):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
    spec=importlib.util.spec_from_file_location('barriers',HERE/'t25-brackets.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
    summary=[];angular=[]
    for kind in ['single','binary','late']:
        maps={n:read(HERE/f't25-map-{kind}-N{n}.csv') for n in [96,192]}
        for n,rs in maps.items():
            assert len(rs)==len({r['id'] for r in rs})
            bs=b.brackets(rs);bar=b.pointwise_barriers(rs)
            b.write(HERE/f't25-brackets-{kind}-N{n}.csv',bs)
            b.write(HERE/f't25-barriers-{kind}-N{n}.csv',bar)
            counts=Counter(r['status'] for r in rs)
            summary.append(dict(kind=kind,time=rs[0]['time'],N=n,surfaces=len(rs),resolved=counts['RESOLVED'],unresolved=counts['UNRESOLVED'],
                invalid=counts['INVALID_GEOMETRY'],adjacent_pointwise=sum(r['status'].startswith('POINTWISE') for r in bs),
                average_only=sum(r['status']=='AVERAGE_ONLY_WEAKER' for r in bs),pointwise_barrier_families=len(bar)))
        other={r['id']:r for r in maps[192]}
        for x in maps[96]:
            y=other[x['id']];assert x['status']==y['status']
            z=dict(kind=kind,id=x['id'],status=x['status'])
            for c in ['A','Q','theta_min','theta_max','theta_mean','theta_rms','negative_area_fraction']:
                z[c+'_N96']=x[c];z[c+'_N192']=y[c];z[c+'_spread']=abs(float(x[c])-float(y[c]))
            angular.append(z)
    write(HERE/'t25-summary.csv',summary);write(HERE/'t25-angular.csv',angular)
    receipts=[]
    for kind in ['single','binary']:
        d=ROOT/'initialization'/kind; r=json.loads((d/'initialization.resources.json').read_text())
        assert r['returncode']==0 and not r['gate_reason'] and int((d/'done.exit').read_text())==0
        assert 'GRChombo finished.' in (d/'run.log').read_text()
        assert 'max_steps = 0' in (d/'params.txt').read_text() and 'stop_time = 0' in (d/'params.txt').read_text()
        receipts.append(dict(case='initialization-'+kind,**{c:r[c] for c in ['returncode','wall_seconds','peak_rss_bytes','cap_GB','gate_reason']}))
    for d in sorted((ROOT/'maps').glob('*-N*')):
        r=json.loads((d/'native-map.resources.json').read_text())
        assert r['returncode']==0 and not r['gate_reason'] and int((d/'done.exit').read_text())==0
        assert 'T25_EXPANSION_MAP_COMPLETE; no advances' in (d/'run.log').read_text()
        assert r['peak_rss_bytes']<3e9
        assert 'ems_data_path = /T25-FORBIDDEN-static' in (d/'params.txt').read_text()
        receipts.append(dict(case=d.name,**{c:r[c] for c in ['returncode','wall_seconds','peak_rss_bytes','cap_GB','gate_reason']}))
    stages=[]
    for d in sorted((ROOT/'finders').glob('*-N*')):
        if not (d/'finder.resources.json').exists():continue
        r=json.loads((d/'finder.resources.json').read_text());rs=read(d/'finder.csv')
        assert r['returncode'] in [0,1] and not r['gate_reason'] and int((d/'done.exit').read_text())==r['returncode']
        assert 'T25_FINDER_COMPLETE' in (d/'run.log').read_text() and r['peak_rss_bytes']<3e9
        if r['returncode']==0:assert len(rs)==3 and all(s['status']=='FOUND' and float(s['expansion_squared'])<=float(s['threshold']) for s in rs)
        else:assert rs[-1]['status'] in ['TIME_CAP','FLOOR','UPDATE_CAP']
        receipts.append(dict(case=d.name,**{c:r[c] for c in ['returncode','wall_seconds','peak_rss_bytes','cap_GB','gate_reason']}))
        stages += [dict(case=d.name,**s) for s in rs]
    write(HERE/'t25-receipts.csv',receipts);write(HERE/'t25-finder-stages.csv',stages)
    base=Path('/private/tmp/ems-t24/source/Source/RHFinder')
    author=[]
    for name in ['RHSurf.hpp','RHUnion.hpp','RHBandedMatrixInverter.hpp']:
        current=sha(base/name)
        content=subprocess.check_output(['git','show','d507438:Source/RHFinder/'+name])
        assert current==hashlib.sha256(content).hexdigest()
        author.append(dict(file=name,sha256=current,unchanged_author=True))
    single=read(HERE/'t25-root-single.csv')[0]
    assert .00604<float(single['a'])<.00606
    assert float(single['finder_theta_max_difference'])<1e-9
    sources=[Path('/Users/auroradysis/Workspace/EMS/artifacts/echo-binary/a0.9-e8.trumpet'),
        Path('/Users/auroradysis/Workspace/EMS/.data/echo-binary/production-outer16-a20-r32-a10/n32-r6.ctt')]
    cps=[ROOT/'initialization'/k/'chk/EMS_000000.2d.hdf5' for k in ['single','binary']]
    cps.append(Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall/chk/EMS_000152.2d.hdf5'))
    inputs=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size,role='initialization only' if p in sources else 'offline map geometry') for p in sources+cps]
    assert inputs[1]['sha256']=='46253cb8d3bf6255b149b44bed7d4899cd97c697fe677a9846649bfbb2bdc4ee'
    header=sources[0].read_text().splitlines(); area=float(next(s.split('=',1)[1] for s in header if s.startswith('# A_H=')))
    result=dict(status='READY-EXCEPT',scope='T25 completed; binary finder remains unqualified; no early positive-time checkpoints collected',
        all_probes_verified=True,author_finder=author,inputs=inputs,initial_profile_header_area=area,
        card_area_benchmark=.2218,benchmark_discrepancy=(area/.2218-1),
        native_sign='K_ij=-(partial_t gamma_ij-Lie_beta gamma_ij)/(2 alpha); theta_plus=Div s+Kss-K',
        interpolation='unchanged native AMRInterpolator<Lagrange<4>>, including native derivative queries and production evolution ghost fill',
        source_archive='d507438',no_dirty_source_compiled=True,static_inputs_only_in_initializers=True,
        peak_map_rss_bytes=max(r['peak_rss_bytes'] for r in receipts if '-N' in r['case']),
        raw_angles_directory=str(ROOT/'maps'),
        hierarchy='low rung L0-12; h0=1.75 M_i, h12=.00042724609375 M_i; initial radius about 14 cells')
    (HERE/'t25-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],summary=summary,peak=result['peak_map_rss_bytes'])))
if __name__=='__main__':main()
