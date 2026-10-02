#!/usr/bin/env python3
"""Verify each native receipt, then apply the unchanged T16b norm/margin rule."""
import sys
sys.dont_write_bytecode=True
import csv, json, math
from pathlib import Path
import numpy as np
import importlib.util
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t18c');OUT=ROOT/'analysis';OUT.mkdir(exist_ok=True)
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
w=module('t18cwork',HERE/'t18c-work.py')
b=module('t18cbaseline',HERE/'t16b-analyze.py');a=b.a
t17=module('t18csampling',HERE.parents[1].parent/'wt-native-t4/Tests/EMSNative/t17-analyze.py')
REG=json.loads((HERE/'t18c-registration.json').read_text());CLOCKS=np.array(REG['native_clocks_M'])
FIELDS=('Gamma','metric_Gamma','shift','lapse');NORMS=('peak','RMS')
def rows(p):return list(csv.DictReader(p.open()))
def save(name,data):
    assert data,name
    with (HERE/('t18c-'+name+'.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=data[0]);writer.writeheader();writer.writerows(data)
def norm(v,kind):return float(np.max(abs(v))) if kind=='peak' else a.rms(v)
def receipt(name):
    d=ROOT/'evolution'/name;r=json.loads((d/(name+'.resources.json')).read_text())
    assert (d/'done.exit').read_text().strip()=='0'
    assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason']
    assert r['peak_rss_bytes']<4e9 and r['tree_peak_bytes']<4.5e9
    log=(d/'run.log').read_text();assert 'T16b native initial audit complete' in log
    if name.startswith('launch'):
        assert 'T18c launch finished.' in log
        pos=log.index('T18c evolution starts after completed t=0 audit')
        assert 'Read EMSTRUMPET' not in log[pos:] and 'EMSTRUMPET radial solution ready' not in log[pos:]
        # The Tests main catches the same native stop exception; its lossless
        # endpoint frames, rather than the production main's CSV, prove clocks.
        native_clocks=sorted({fr['meta'][0] for fr in a.decoder.frames(d/'t13-t7-stage-L12.xz') if fr['phase']==50 and len(fr['cells'])})
        assert np.array_equal(native_clocks,CLOCKS),(name,native_clocks)
        floors=[x for f in d.glob('t13-floors-*.csv') for x in rows(f)]
        assert floors and all(int(x[k])==0 for x in floors for k in ('chi_activations','lapse_activations','nonfinite'))
    else:assert 'T18c audit-only finished.' in log
    return r
def verify_inputs():
    for filename in ('t18c-inputs.json','t18c-executable-inputs.json'):
        for r in json.loads((HERE/filename).read_text()):
            assert w.digest(Path(r['path']))==r['sha256'],r['path']
    assert w.digest(w.DATA)==REG['companion_sha256']
def initial():
    ranges=[];sources=[];limits=[];constraints=[];interfaces=[];audits=[]
    panels={int(r['panel']):{k:(r[k] if k in ('kind','boundary') else float(r[k])) for k in r} for r in rows(HERE/'t18c-panels.csv')}
    old_constraints={kind:rows(HERE/f't16b-{kind}.csv') for kind in ('interface-constraints','constraints')}
    for scope in ('full','launch'):
        for mode in ('sqrt','geometric'):
            name=scope+'-'+mode;d=ROOT/'evolution'/name;r=receipt(name);rr=rows(d/'t18c-ranges.csv')
            assert max(int(x['level']) for x in rr)==12
            expected=[list(map(int,line.split())) for line in (HERE/f't18c-{scope}-boxes.txt').read_text().splitlines()]
            actual=[(int(x['level']),*[int(x[k]) for k in ('x0','y0','x1','y1')]) for x in rr]
            assert sorted(actual)==sorted(map(tuple,expected))
            for key in ('a_left','a_right','alpha_bin','native_lapse'):
                lo=min(float(x[key+'_min']) for x in rr);hi=max(float(x[key+'_max']) for x in rr)
                assert 0<lo<=hi<=1
                ranges.append(dict(run=name,quantity=key,minimum=lo,maximum=hi,
                    d32_coarse_min=next(float(x['minimum']) for x in rows(HERE/'t16b-ranges.csv') if x['run']=='coarse-'+mode and x['quantity']==key),
                    d32_coarse_max=next(float(x['maximum']) for x in rows(HERE/'t16b-ranges.csv') if x['run']=='coarse-'+mode and x['quantity']==key)))
            for hole in ('left','right'):
                raw=[x for x in rows(d/'t18c-source.csv') if x['hole']==hole];assert len(raw)==4
                assert all(int(x['Base_bit_mismatches'])==0 and float(x['driver_cancel'])==0 for x in raw)
                for field in ('physical_Gamma_n','KO_Gamma_n','total_Gamma_n','physical_K','KO_K','lapse_rhs'):
                    v=np.array([float(x[field]) for x in raw]);old=next(x for x in rows(HERE/'t16b-source-summary.csv') if x['run']=='coarse-'+mode and x['hole']==hole and x['field']==field)
                    sources.append(dict(run=name,hole=hole,field=field,peak=float(max(abs(v))),RMS=float(np.sqrt(np.mean(v*v))),d32_coarse_peak=float(old['peak']),d32_coarse_RMS=float(old['RMS'])))
            for x in rows(d/'t18c-limits.csv'):
                old=next(z for z in rows(HERE/'t16b-puncture-limits.csv') if z['run']=='coarse-'+mode and z['hole']==x['hole'] and float(z['R'])==float(x['R']))
                limits.append(dict(run=name,**x,d32_alpha_bin_over_a_own=old['alpha_bin_over_a_own']))
            accum={}
            def add(region,subset,v,weight,crossing):
                key=(region,subset)
                if key not in accum:accum[key]=[0,0.,np.zeros(5),np.zeros(5),0]
                z=accum[key];z[0]+=1;z[1]+=weight;z[2]=np.maximum(z[2],abs(v));z[3]+=weight*v*v;z[4]+=crossing
            with (d/'t18c-t0-native.csv').open() as f:
                for x in csv.DictReader(f):
                    X,Y=float(x['x']),float(x['y']);v=np.array([float(x[k]) for k in ('Ham','Mom','GaussE','GaussB','C_Gamma')]);weight=float(x['weight'])
                    assert np.isfinite(v).all() and weight>0
                    rl=math.hypot(X+8,Y);rrad=math.hypot(X-8,Y);rg=math.hypot(X,Y);seam=int(x['amr_seam']);regions=[]
                    if .01<=rl<=.02:regions.append('collar-left')
                    if .01<=rrad<=.02:regions.append('collar-right')
                    if abs(X)<=8 and 0<=Y<=8:regions.append('interhole')
                    if 80<=rg<=128:regions.append('far')
                    if min(rl,rrad)>=.01 and rg<=128 and not x['sheets']:regions.append('elsewhere-off-sheet')
                    for region in regions:
                        add(region,'all',v,weight,0)
                        if not seam:add(region,'AMR-face-excluded',v,weight,0)
                    for pair in x['sheets'].split(';') if x['sheets'] else []:
                        add('sheet-'+pair,'all',v,weight,int(x['native_crossing']))
                        if not seam:add('sheet-'+pair,'AMR-face-excluded',v,weight,int(x['native_crossing']))
            for (region,subset),z in accum.items():
                for k,field in enumerate(('Ham','Mom','GaussE','GaussB','C_Gamma')):
                    old=next((r for r in old_constraints['interface-constraints' if region.startswith('sheet') else 'constraints'] if r['run']=='coarse-'+mode and r['region']==region and r['subset']==subset and r['field']==field),None)
                    rec=dict(run=name,region=region,subset=subset,field=field,cells=z[0],volume=z[1],peak=z[2][k],RMS=math.sqrt(z[3][k]/z[1]),native_crossing_cells=z[4],
                        d32_coarse_peak=float(old['peak']) if old else '',d32_coarse_RMS=float(old['RMS']) if old else '')
                    if region.startswith('sheet'):
                        p,q=map(int,region[6:].split(':'));rec['interface_kind']=b.adjacent(panels[p],panels[q]);interfaces.append(rec)
                    else:constraints.append(rec)
            audits.append(dict(run=name,native_stop_M=CLOCKS[-1] if scope=='launch' else 0.,
                nonlapse_values=sum(int(x['nonlapse_values']) for x in rr),nonlapse_bit_mismatches=0 if mode=='geometric' else 'reference',
                native_RHS_bit_mismatches=0,driver_cancellation=0,selector_checks=sum(int(x['selector_checks']) for x in rr),
                geometric_native_max_error=max(float(x['geometric_native_max_error']) for x in rr),peak_RSS_bytes=r['peak_rss_bytes'],wall_seconds=r['wall_seconds']))
    for scope in ('full','launch'):
        assert w.digest(ROOT/'evolution'/(scope+'-sqrt')/'t18c-t0-native.csv')==w.digest(ROOT/'evolution'/(scope+'-geometric')/'t18c-t0-native.csv')
    save('ranges',ranges);save('sources',sources);save('limits',limits);save('constraints',constraints);save('interfaces',interfaces);save('audit',audits)
    expected={f'{p}:{q}':b.adjacent(panels[p],panels[q]) for p in panels for q in panels if p<q and b.adjacent(panels[p],panels[q])!='multiple-interface corner transition'}
    observed={x['region'][6:] for x in interfaces if x['run']=='full-geometric' and x['native_crossing_cells']>0}
    save('interface-coverage',[dict(interface=k,kind=v,observed=k in observed,status='native stencil observed' if k in observed else 'unobserved; no inferred pass/order') for k,v in expected.items()])
    return ranges,sources,limits,constraints,interfaces,audits
def capture(name,hole):
    h=REG['h_finest_M'];full=name.startswith('full');mid=2240. if full else 448.;center=mid+(-8. if hole=='left' else 8.)
    vectors={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/np.sqrt(2)}
    if hole=='right':vectors={k:v*np.array([-1.,1.]) for k,v in vectors.items()}
    cells=t17.select(h,center,vectors);times=[];group=[];previous=None
    payload=OUT/(name+'-'+hole+'-input.bin');output=OUT/(name+'-'+hole+'-q.bin')
    with payload.open('wb') as stream:
        def consume(frames):
            lo,u=a.dense(frames);rec=np.empty((len(cells),3+49*28));rec[:,0]=h;rec[:,1]=(cells[:,1]+.5)*h;rec[:,2]=(cells[:,0]+.5)*h-center
            for j in range(-3,4):
                for i in range(-3,4):
                    v=u[cells[:,1]+j-lo[1],cells[:,0]+i-lo[0]];assert np.isfinite(v).all()
                    rec[:,3+(j+3)*7+i+3::49]=v
            rec.tofile(stream);times.append(frames[0]['meta'][0])
        for fr in a.decoder.frames(ROOT/'evolution'/name/'t13-t7-stage-L12.xz'):
            if fr['phase']!=50 or not len(fr['cells']):continue
            if (np.mean((fr['cells'][:,0]+.5)*h)<mid)!=(hole=='left'):continue
            time=fr['meta'][0]
            if previous is not None and time!=previous:consume(group);group=[]
            group.append(fr);previous=time
        if group:consume(group)
    assert np.array_equal(times,[0.] if full else CLOCKS),(name,hole,times)
    w.measured('replay-'+name+'-'+hole,OUT,[ROOT/'replay.ex','--state',payload,output])
    (OUT/(name+'-'+hole+'-input.sha256')).write_text(w.digest(payload)+'\n');payload.unlink()
    q=np.memmap(output,mode='r',dtype='f8',shape=(len(times),len(cells),180));assert np.isfinite(q).all() and np.count_nonzero(q[:,:,144])==0
    np.savez(OUT/(name+'-'+hole+'-meta.npz'),cells=cells,times=times)
    profiles={}
    if not full:
        for ray,vec in vectors.items():
            profiles[ray]={n:np.array([t17.sample(cells,frame,vec,n,h,center) for frame in q]) for n in (8,10)}
            np.savez_compressed(OUT/(name+'-'+hole+'-'+ray+'.npz'),times=times,r=a.R,p8=profiles[ray][8],p10=profiles[ray][10])
    return q,cells,vectors,profiles
def launch():
    history=[];ratios=[];puncture=[];controls=[]
    for mode in ('sqrt','geometric'):
        for hole in ('left','right'):
            full,fc,_,_=capture('full-'+mode,hole);q,cells,vectors,p=capture('launch-'+mode,hole)
            assert np.array_equal(fc,cells+np.array([1024*4096,0]))
            x=np.ascontiguousarray(full[0,:,:28]);y=np.ascontiguousarray(q[0,:,:28]);bits=int(np.count_nonzero(x.view('u8')!=y.view('u8')))
            assert bits==0,(mode,hole,'translated launch t0 differs',bits)
            controls.append(dict(mode=mode,hole=hole,cells=len(cells),values=x.size,bit_mismatches=bits))
            center=448+(-8 if hole=='left' else 8);c=math.floor(center/REG['h_finest_M']-.5)
            mask=(cells[:,0]>=c)&(cells[:,0]<=c+1)&(cells[:,1]<2);assert mask.sum()==4
            for step,time in enumerate(CLOCKS):
                for field,index in [('Gamma1',11),('Gamma2',12),('shift1',14),('shift2',15),('lapse',13)]:
                    puncture.append(dict(mode=mode,hole=hole,step=step,time_M=time,field=field,peak_change=float(max(abs(q[step,mask,index]-q[0,mask,index]))),KO_peak=float(max(abs(q[step,mask,56+index])))))
            for ray,vec in vectors.items():
                fs={n:[a.fields(frame,vec) for frame in p[ray][n]] for n in (8,10)}
                for step,time in enumerate(CLOCKS):
                    for field in FIELDS:
                        u=fs[8][step][field]-fs[8][0][field];control=fs[10][step][field]-fs[10][0][field]
                        for kind in NORMS:
                            amp=norm(u,kind);spread=norm(u-control,kind)
                            history.append(dict(mode=mode,hole=hole,ray=ray,field=field,norm=kind,step=step,time_M=time,amplitude=amp,spread=spread,margin_5x=amp/(5*spread) if spread else math.inf if amp else 0))
    lookup={(r['mode'],r['hole'],r['ray'],r['field'],r['norm'],r['step']):r for r in history}
    baseline=rows(HERE/'t16b-launch-endpoint.csv')
    for g in history:
        if g['mode']!='geometric' or not g['step']:continue
        s=lookup['sqrt',g['hole'],g['ray'],g['field'],g['norm'],g['step']];G,S=g['amplitude'],s['amplitude'];eg,es=5*g['spread'],5*s['spread']
        old=next(x for x in baseline if x['rung']=='coarse' and all(x[k]==g[k] for k in ('hole','ray','field','norm')))
        ratios.append(dict(hole=g['hole'],ray=g['ray'],field=g['field'],norm=g['norm'],step=g['step'],time_M=g['time_M'],geometric_amplitude=G,sqrt_amplitude=S,
            ratio=G/S if S else '',lower_5x=max(0,G-eg)/(S+es) if S+es else '',upper_5x=(G+eg)/(S-es) if S>es else math.inf,
            geometric_margin_5x=g['margin_5x'],sqrt_margin_5x=s['margin_5x'],reduction_margin_5x=(S-G)/(eg+es) if eg+es else math.inf if S>G else -math.inf,
            qualified=G>eg and S>es,endpoint=g['step']==9,d32_coarse_endpoint_ratio=float(old['ratio'])))
    save('launch-history',history);save('launch-ratios',ratios);save('launch-endpoint',[r for r in ratios if r['endpoint']]);save('native-puncture',puncture);save('translation-control',controls)
    summary=[]
    for hole in ('left','right'):
        for ray in ('axis','diagonal'):
            for field in FIELDS:
                for kind in NORMS:
                    rs=[r for r in ratios if (r['hole'],r['ray'],r['field'],r['norm'])==(hole,ray,field,kind)]
                    summary.append(dict(hole=hole,ray=ray,field=field,norm=kind,qualified=sum(r['qualified'] for r in rs),
                        uncertainty_dominated=sum(not r['qualified'] for r in rs),resolved_reductions=sum(r['qualified'] and r['upper_5x']<1 for r in rs),enhanced_steps=' '.join(str(r['step']) for r in rs if isinstance(r['ratio'],float) and r['ratio']>1)))
    save('history-summary',summary)
    return ratios,summary
def report(initial_result,ratios,summary):
    ranges,sources,limits,constraints,interfaces,audits=initial_result;end=[r for r in ratios if r['endpoint']]
    gamma_pass=all(r['qualified'] and r['upper_5x']<1 for r in end if r['field']=='Gamma')
    source_pass=all(next(r['peak'] for r in sources if r['run']=='full-geometric' and r['hole']==hole and r['field']=='total_Gamma_n')<next(r['peak'] for r in sources if r['run']=='full-sqrt' and r['hole']==hole and r['field']=='total_Gamma_n') for hole in ('left','right'))
    admitted=gamma_pass and source_pass
    status=dict(status='READY-EXCEPT',diagnostic_lapse_rule_pass=admitted,geometric_endpoint_Gamma_rule=gamma_pass,native_total_Gamma_source_rule=source_pass,
        recommendation='ems_use_geometric_initial_lapse=true' if admitted else 'sqrt(chi); geometric lapse fails or is unresolved under T16b margin rule',
        diagnostic_low_merger_admission=admitted,convergence_admission=False,
        limits='Only LOW rung and nine positive native clocks; no T16b two-rung admission, long-time binary stability or companion certification.',
        companion='documented diagnostic, D2(logpsi) 2.8922936531e-5 > 1e-9',peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in audits))
    (HERE/'t18c-status.json').write_text(json.dumps(status,indent=2)+'\n')
    fmt=lambda v:f'{float(v):.6g}' if v!='' else 'unresolved'
    table=b.table
    text=['READY-EXCEPT',
        'Native t=0 audit and launch receipts pass the execution checks. '+status['recommendation']+'. Diagnostic LOW merger lapse rule: '+('PASS' if admitted else 'FAIL / unresolved')+'. Convergence admission remains false.',
        'LOW d16 uses h=1.75/4096, R_h/h=14.1381, box size24 and dt_multiplier=0.5. Every production cell and ghost enters the strict t0 bounds and bit comparison. The unused-array initialization fixture and translated-launch t0 checks are bit-identical. The maximum measured native-leg RSS is '+fmt(status['peak_RSS_bytes']/1e9)+' GB.',
        table(['Run','Quantity','Min','Max','d32 coarse min','d32 coarse max'],[[r['run'],r['quantity'],fmt(r['minimum']),fmt(r['maximum']),fmt(r['d32_coarse_min']),fmt(r['d32_coarse_max'])] for r in ranges if r['run'].startswith('full')]),
        table(['Hole','R','alpha_bin/a_own d16','d32'],[[r['hole'],r['R'],fmt(r['alpha_bin_over_a_own']),fmt(r['d32_alpha_bin_over_a_own'])] for r in limits if r['run']=='full-geometric']),
        table(['Run','Hole','Source','Peak','RMS','d32 peak','d32 RMS'],[[r['run'],r['hole'],r['field'],fmt(r['peak']),fmt(r['RMS']),fmt(r['d32_coarse_peak']),fmt(r['d32_coarse_RMS'])] for r in sources if r['run'].startswith('full') and r['field'] in ('physical_Gamma_n','KO_Gamma_n','total_Gamma_n')]),
        table(['Region','Constraint','d16 peak','d16 RMS','d32 peak','d32 RMS'],[[r['region'],r['field'],fmt(r['peak']),fmt(r['RMS']),fmt(r['d32_coarse_peak']),fmt(r['d32_coarse_RMS'])] for r in constraints if r['run']=='full-geometric' and r['subset']=='all']),
        'Full sheet and AMR-face-excluded norms and panel crossing counts are retained in [interfaces](t18c-interfaces.csv). Sheet masks retain the same selector/stencil definition; panel surfaces move with the d16 companion, and the LOW reference spacing differs. These are not a fixed-separation convergence pair.',
        table(['Hole','Ray','Field','Norm','Ratio','five-spread interval','signal margins G / S','reduction margin','d32 endpoint ratio'],[[r['hole'],r['ray'],r['field'],r['norm'],fmt(r['ratio']),fmt(r['lower_5x'])+' / '+fmt(r['upper_5x']),fmt(r['geometric_margin_5x'])+' / '+fmt(r['sqrt_margin_5x']),fmt(r['reduction_margin_5x']),fmt(r['d32_coarse_endpoint_ratio'])] for r in end]),
        'All nine positive native clocks are retained in [launch ratios](t18c-launch-ratios.csv); early enhancements and unresolved entries remain in [history summary](t18c-history-summary.csv). Margins are T16b five-spread sampling sensitivity estimates, not rigorous error bounds. Native four-cell changes and KO peaks are in [puncture corroboration](t18c-native-puncture.csv).',
        'The d16 endpoint is 0.001922607421875 M_i, 3.077% before d32 T16b 0.001983642578125. d32 T16b has R_h/h49.483/98.967 and dt=h/4; separation, grid phase, spatial resolution, dt and endpoint all differ. The table compares the actual endpoints without time alignment; it cannot attribute differences solely to d. No second-rung or same-grid temporal control was requested here.',
        'The native constraints are identical between lapse variants. The documented d16 D2(logpsi) jump is 2.8922936531e-5 versus d32 about8.89e-6; both exceed the unchanged1e-9 gate. This initial-lapse choice cannot repair that defect. '+('The requested LOW diagnostic lapse criterion passes at both holes; the companion may enter an explicitly diagnostic low-resolution merger with geometric lapse.' if admitted else 'The requested LOW geometric lapse criterion fails or is unresolved; this audit does not admit the proposed geometric-lapse diagnostic merger.')+' The strict T16b two-rung/convergence criteria are not met.',
        'Reproduce with t18c-work.py generate/census/prepare/build/smoke, then the measured detached pipeline plan, then t18c-analyze.py. Every native child has its own receipt; zero done.exit is execution evidence only. No SSH or commit.']
    report='\n\n'.join(text)+'\n';(HERE/'t18c-report.md').write_text(report)
    readme=HERE/'README.md';original=readme.read_text();marker='### Completed T18c evidence'
    if marker in original:original=original[:original.index(marker)]
    readme.write_text(original.rstrip()+'\n\n'+marker+'\n\n'+report)
    resources=[]
    for f in ROOT.rglob('*.resources.json'):
        r=json.loads(f.read_text());resources.append(dict(record=str(f),process=r['process'],peak_RSS_bytes=r['peak_rss_bytes'],returncode=r['returncode'],child_returncode=r.get('child_measurement',{}).get('returncode',''),gate_reason=r['gate_reason']))
    save('resources',resources)
    files=sorted([p for p in HERE.glob('t18c-*') if p.name!='t18c-manifest.txt']+[HERE/'README.md']+[p for p in ROOT.rglob('*') if p.is_file() and p.suffix in ('.json','.csv','.xz','.ex','.npz')])
    (HERE/'t18c-manifest.txt').write_text(''.join(w.digest(p)+'  '+str(p)+'\n' for p in files))
    print(json.dumps(status,indent=2))
if __name__=='__main__':
    verify_inputs()
    report(initial(),*launch())
