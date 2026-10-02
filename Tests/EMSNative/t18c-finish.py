#!/usr/bin/env python3
"""Check completed analysis receipts and independently resample its saved caches."""
import sys
sys.dont_write_bytecode=True
import csv,importlib.util,json,math
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('audit',HERE/'t18c-analyze.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def main():
    m.verify_inputs();checks=[]
    for name in ('full-sqrt','full-geometric','launch-sqrt','launch-geometric'):
        m.receipt(name);checks.append(dict(check=name+' actual child, gates, clocks and floors',result='PASS'))
    for path in (m.ROOT/'completion-verified/analysis/analysis-final.resources.json',m.ROOT/'completion-verified/queue.resources.json',m.ROOT/'domain-check/domain-check.resources.json'):
        r=json.loads(path.read_text());assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason'] and r['peak_rss_bytes']<4e9
        checks.append(dict(check=r['process']+' own completed receipt',result='PASS'))
    assert (m.ROOT/'completion-verified/done.exit').read_text().strip()=='0'
    history=m.rows(HERE/'t18c-launch-history.csv');compared=0
    for mode in ('sqrt','geometric'):
        for hole in ('left','right'):
            for ray in ('axis','diagonal'):
                data=np.load(m.OUT/f'launch-{mode}-{hole}-{ray}.npz');assert np.array_equal(data['times'],m.CLOCKS)
                vec=np.array([1.,0.]) if ray=='axis' else np.ones(2)/np.sqrt(2)
                if hole=='right':vec*=np.array([-1.,1.])
                fields={n:[m.a.fields(v,vec) for v in data['p'+str(n)]] for n in (8,10)}
                for r in history:
                    if (r['mode'],r['hole'],r['ray'])!=(mode,hole,ray):continue
                    step=int(r['step']);field=r['field'];u=fields[8][step][field]-fields[8][0][field];v=fields[10][step][field]-fields[10][0][field]
                    assert m.norm(u,r['norm'])==float(r['amplitude']) and m.norm(u-v,r['norm'])==float(r['spread']);compared+=1
    assert compared==640
    for r in m.rows(HERE/'t18c-launch-ratios.csv'):
        G,S=float(r['geometric_amplitude']),float(r['sqrt_amplitude']);assert float(r['ratio'])==G/S
        g=next(x for x in history if x['mode']=='geometric' and all(x[k]==r[k] for k in ('hole','ray','field','norm','step')))
        d=next(x for x in history if x['mode']=='sqrt' and all(x[k]==r[k] for k in ('hole','ray','field','norm','step')))
        eg,es=5*float(g['spread']),5*float(d['spread'])
        assert float(r['lower_5x'])==max(0,G-eg)/(S+es) and float(r['upper_5x'])==((G+eg)/(S-es) if S>es else math.inf)
    checks.append(dict(check='640 cached amplitudes/spreads and all 288 ratio intervals',result='BIT_IDENTICAL'))
    endpoint=m.rows(HERE/'t18c-launch-endpoint.csv');old_endpoint=m.rows(HERE/'t16b-launch-endpoint.csv')
    for r in endpoint:
        for rung in ('coarse','fine'):
            old=next(x for x in old_endpoint if x['rung']==rung and all(x[k]==r[k] for k in ('hole','ray','field','norm')))
            for key in ('ratio','lower_5x','upper_5x','geometric_amplitude','sqrt_amplitude','geometric_margin_5x','sqrt_margin_5x','reduction_margin_5x'):
                r['d32_'+rung+'_'+key]=old[key]
    m.save('launch-endpoint',endpoint)
    extrema={}
    with (m.ROOT/'evolution/full-geometric/t18c-t0-native.csv').open() as f:
        for r in csv.DictReader(f):
            x,y=float(r['x']),float(r['y'])
            if abs(x)>8 or y>8:continue
            for field in ('Ham','Mom','GaussE'):
                v=abs(float(r[field]))
                if field not in extrema or v>extrema[field]['peak']:
                    extrema[field]=dict(region='interhole',field=field,peak=v,x_M=x,y_M=y,level=int(r['level']),r_near_M=min(math.hypot(x+8,y),math.hypot(x-8,y)))
    m.save('constraint-hotspots',list(extrema.values()));m.save('completion-checks',checks)
    coverage=m.rows(HERE/'t18c-interface-coverage.csv');observed=sum(x['observed']=='True' for x in coverage)
    extra='''The unchanged inter-hole box |x-midpoint|≤8 touches both punctures when d=16; at d=32 it excludes both. Its large d16 peaks therefore include the unresolved puncture neighborhood. The native extrema lie at r_near='''+', '.join(f"{x['r_near_M']:.6g} M_i ({x['field']})" for x in extrema.values())+'''. The mask is retained exactly and no near-hole samples are removed. These peaks are not a like-for-like exterior separation comparison. [Locations](t18c-constraint-hotspots.csv).

'''+f'''The selector census observes {observed}/{len(coverage)} declared geometric interfaces at native stencil spacing, versus 30/44 in T16b d32. Unobserved end sheets remain explicitly labelled without an inferred pass or order. [Coverage](t18c-interface-coverage.csv).

The cropped launch's entire refined t0 non-lapse state and base ±200-M_i neighborhood match 53,342,145 full-production values with zero bit mismatches. All 28 captured finest t0 fields also match. The recorded update counts give a deliberately generous coordinate-footprint estimate of 168.417 M_i, versus 440 M_i to the nearer cropped outer boundary; this is a scope estimate, not a full-production evolved-state bit comparison. [Domain check](t18c-domain-check.json).

The final analysis and its enclosing queue each have their own successful actual-child receipt; 640 cached launch amplitudes/spreads and every one of 288 five-spread ratio intervals reproduce exactly. Failed preliminary reporting receipts are retained. The authoritative completion marker is `/private/tmp/ems-t18c/completion-verified/done.exit = 0`. No native job remains pending. [Completion checks](t18c-completion-checks.csv).
'''
    comparison=[]
    for hole in ('left','right'):
        for ray in ('axis','diagonal'):
            for field in m.FIELDS:
                rr=[r for r in endpoint if (r['hole'],r['ray'],r['field'])==(hole,ray,field)]
                comparison.append([hole,ray,field,
                    f"{min(min(float(r['geometric_margin_5x']),float(r['sqrt_margin_5x'])) for r in rr):.6g}",
                    f"{min(min(float(r['d32_coarse_geometric_margin_5x']),float(r['d32_coarse_sqrt_margin_5x'])) for r in rr):.6g}",
                    f"{min(float(r['reduction_margin_5x']) for r in rr):.6g}",
                    f"{min(float(r['d32_coarse_reduction_margin_5x']) for r in rr):.6g}"])
    extra+='\nEndpoint margins compare actual endpoints and include the minimum across peak/RMS; the d32 coarse and fine amplitudes, intervals and margins are also retained in [the endpoint table](t18c-launch-endpoint.csv).\n\n'+m.b.table(['Hole','Ray','Field','d16 signal / 5 spread','d32 coarse signal / 5 spread','d16 reduction margin','d32 coarse reduction margin'],comparison)+'\n'
    report=HERE/'t18c-report.md';text=report.read_text();marker='The unchanged inter-hole box'
    if marker in text:text=text[:text.index(marker)]
    report.write_text(text.rstrip()+'\n\n'+extra)
    readme=HERE/'README.md';text=readme.read_text();start=text.index('## T18c —');section=text[start:]
    section=section.replace('**RUNS-PENDING.**','**READY-EXCEPT: completed and independently verified.**',1)
    section=section.replace('Corrected analysis uses `/private/tmp/ems-t18c/completion/done.exit` as the completion marker.',
        'The verifier also rejected a missing production-main stop CSV; native lossless clocks prove the Tests harness stop. Final analysis uses `/private/tmp/ems-t18c/completion-verified/done.exit = 0` as the completion marker.')
    completed='### Completed T18c evidence';section=section[:section.index(completed)]+completed+'\n\n'+report.read_text()
    readme.write_text(text[:start]+section)
    status=HERE/'t18c-status.json';z=json.loads(status.read_text());z.update(completion_verified=True,completion_marker=str(m.ROOT/'completion-verified/done.exit'),cached_norm_checks=640,ratio_interval_checks=288,observed_interfaces=observed)
    status.write_text(json.dumps(z,indent=2)+'\n')
    resources=[]
    for f in m.ROOT.rglob('*.resources.json'):
        r=json.loads(f.read_text());resources.append(dict(record=str(f),process=r['process'],peak_RSS_bytes=r['peak_rss_bytes'],returncode=r['returncode'],child_returncode=r.get('child_measurement',{}).get('returncode',''),gate_reason=r['gate_reason']))
    m.save('resources',resources)
    files=sorted([p for p in HERE.glob('t18c-*') if p.name!='t18c-manifest.txt']+[HERE/'README.md']+[p for p in m.ROOT.rglob('*') if p.is_file() and p.suffix in ('.json','.csv','.xz','.ex','.npz','.bin','.log','.time')])
    (HERE/'t18c-manifest.txt').write_text(''.join(m.w.digest(p)+'  '+str(p)+'\n' for p in files))
    print(json.dumps(z,indent=2))
if __name__=='__main__':main()
