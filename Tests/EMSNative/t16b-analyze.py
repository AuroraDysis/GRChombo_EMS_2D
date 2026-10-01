#!/usr/bin/env python3
"""Stream native binary receipts and replay captured current states, never a seed."""
import sys
sys.dont_write_bytecode=True
import csv,hashlib,importlib.util,json,math,os,re
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t16b');OUT=ROOT/'analysis';OUT.mkdir(exist_ok=True)
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
w=module('work',HERE/'t16b-work.py');old=module('old_analysis',HERE/'t16-analyze.py');a=old.a
NAMES=['coarse-sqrt','coarse-geometric','fine-sqrt','fine-geometric']
FIELDS=('Gamma','metric_Gamma','shift','lapse');NORMS=('peak','RMS')
def rows(path):return list(csv.DictReader(path.open()))
def save(name,data):w.save('t16b-'+name+'.csv',data)
def hashfile(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()
def fn(v,norm):return float(np.max(abs(v))) if norm=='peak' else a.rms(v)
def receipt(name):
    d=ROOT/'evolution'/name;r=json.loads((d/(name+'.resources.json')).read_text())
    assert (d/'done.exit').read_text().strip()=='0' and r['returncode']==r['child_measurement']['returncode']==0
    assert not r['gate_reason'] and r['peak_rss_bytes']<6e9 and r['tree_peak_bytes']<6.5e9
    return r
def adjacent(p,q):
    if p['center']==q['center'] and p['stretch']==q['stretch']:
        overlap=min(p['mu_hi'],q['mu_hi'])>max(p['mu_lo'],q['mu_lo'])
        radial=any(abs(x-y)<=32*np.finfo(float).eps*max(abs(x),abs(y),1e-20) for x,y in [(p['hi'],q['lo']),(q['hi'],p['lo'])])
        if overlap and radial and p['kind']!='inverse' and q['kind']!='inverse':return 'radial'
        if p['lo']==q['lo'] and p['hi']==q['hi'] and (p['mu_hi']==q['mu_lo'] or q['mu_hi']==p['mu_lo']):return 'angular'
    if p['kind']==q['kind']=='bridge' and p['boundary']==q['boundary']=='plane' and p['center']!=q['center']:return 'bridge-plane'
    for b,c in [(p,q),(q,p)]:
        if b['kind']=='bridge' and b['boundary']=='sphere' and c['kind']=='inverse' and ((b['center']<0 and c['mu_hi']==0) or (b['center']>0 and c['mu_lo']==0)):return 'outer-sphere'
    return 'multiple-interface corner transition'
def initial():
    ranges=[];sources=[];constraints=[];interfaces=[];audits=[];limits=[]
    panels={int(r['panel']):{k:(r[k] if k in ('kind','boundary') else float(r[k])) for k in r} for r in rows(HERE/'t16b-panels.csv')}
    expected={f'{p}:{q}':adjacent(panels[p],panels[q]) for p in panels for q in panels if p<q and adjacent(panels[p],panels[q])!='multiple-interface corner transition'}
    for name in NAMES:
        d=ROOT/'evolution'/name;r=receipt(name);rung,mode=name.split('-');reg=next(x for x in rows(HERE/'t16b-registration.csv') if x['run']==name)
        log=(d/'run.log').read_text();assert 'T16b native initial audit complete' in log and 'GRChombo finished.' in log
        pos=log.index('T16b evolution starts after completed t=0 audit')
        assert 'Read EMSTRUMPET' not in log[pos:]
        stop=rows(d/'t13-stop.csv')[0];assert float(stop['actual_time_M'])==float(reg['native_stop_M'])
        rr=rows(d/'t16b-ranges.csv');assert max(int(x['level']) for x in rr)==13
        for key in ('a_left','a_right','alpha_bin','native_lapse'):
            lo=min(float(x[key+'_min']) for x in rr);hi=max(float(x[key+'_max']) for x in rr)
            assert 0<lo<=hi<=1
            ranges.append(dict(run=name,quantity=key,minimum=lo,maximum=hi))
        # Filled coarse/fine and physical ghosts can differ from a fresh setter
        # evaluation; record that interpolation discrepancy rather than hiding it.
        assert all(math.isfinite(float(x['geometric_native_max_error'])) for x in rr)
        for l in range(1,14):
            level=[x for x in rr if int(x['level'])==l];h=float(level[0]['h']);face=128/2**l
            assert max((int(x['y1'])+1)*h for x in level)==face
            for center in (320,352):
                cover=[x for x in level if int(x['x0'])*h<=center<(int(x['x1'])+1)*h]
                assert cover
                if l>=4:
                    local=[x for x in level if abs((int(x['x0'])+int(x['x1'])+1)*h/2-center)<face]
                    assert min(int(x['x0'])*h for x in local)==center-face
                    assert max((int(x['x1'])+1)*h for x in local)==center+face
        for hole in ('left','right'):
            ss=[x for x in rows(d/'t16b-source.csv') if x['hole']==hole];assert len(ss)==4
            assert all(int(x['Base_bit_mismatches'])==0 and float(x['driver_cancel'])==0 for x in ss)
            for key in ('physical_Gamma_n','KO_Gamma_n','total_Gamma_n','physical_K','KO_K','lapse_rhs'):
                values=[float(x[key]) for x in ss]
                sources.append(dict(run=name,hole=hole,field=key,peak=max(map(abs,values)),RMS=math.sqrt(sum(v*v for v in values)/4)))
        for x in rows(d/'t16b-limits.csv'):limits.append(dict(run=name,**x))
        accum={}
        def add(region,subset,v,weight,crossing=0):
            key=(region,subset)
            if key not in accum:accum[key]=[0,0.,np.zeros(5),np.zeros(5),0]
            z=accum[key];z[0]+=1;z[1]+=weight;z[2]=np.maximum(z[2],abs(v));z[3]+=weight*v*v;z[4]+=crossing
        with (d/'t16b-t0-native.csv').open() as f:
            for x in csv.DictReader(f):
                X=float(x['x']);Y=float(x['y']);v=np.array([float(x[k]) for k in ('Ham','Mom','GaussE','GaussB','C_Gamma')]);weight=float(x['weight'])
                assert np.isfinite(v).all() and weight>0
                rl=math.hypot(X+16,Y);rradius=math.hypot(X-16,Y);rg=math.hypot(X,Y);seam=int(x['amr_seam'])
                regions=[]
                if .01<=rl<=.02:regions.append('collar-left')
                if .01<=rradius<=.02:regions.append('collar-right')
                if abs(X)<=8 and 0<=Y<=8:regions.append('interhole')
                if 80<=rg<=128:regions.append('far')
                if min(rl,rradius)>=.01 and rg<=128 and not x['sheets']:regions.append('elsewhere-off-sheet')
                for region in regions:
                    add(region,'all',v,weight)
                    if not seam:add(region,'AMR-face-excluded',v,weight)
                for pair in x['sheets'].split(';') if x['sheets'] else []:
                    add('sheet-'+pair,'all',v,weight,int(x['native_crossing']))
                    if not seam:add('sheet-'+pair,'AMR-face-excluded',v,weight,int(x['native_crossing']))
        for (region,subset),z in accum.items():
            for k,field in enumerate(('Ham','Mom','GaussE','GaussB','C_Gamma')):
                rec=dict(run=name,rung=rung,mode=mode,region=region,subset=subset,field=field,cells=z[0],volume=z[1],peak=z[2][k],RMS=math.sqrt(z[3][k]/z[1]),native_crossing_cells=z[4])
                if region.startswith('sheet-'):
                    pair=region[6:];p1,p2=map(int,pair.split(':'));rec['interface_kind']=adjacent(panels[p1],panels[p2]);interfaces.append(rec)
                else:constraints.append(rec)
        floor=[x for f in d.glob('t13-floors-*.csv') for x in rows(f)]
        assert floor and all(int(x[k])==0 for x in floor for k in ('chi_activations','lapse_activations','nonfinite'))
        audits.append(dict(run=name,native_stop_M=stop['actual_time_M'],initial_nonlapse_values=sum(int(x['nonlapse_values']) for x in rows(d/'t16b-ranges.csv')),
            nonlapse_bit_mismatches=0 if mode=='geometric' else 'reference',chi_floors=0,lapse_floors=0,nonfinite=0,static_reader_after_advance=0,
            peak_RSS_bytes=r['peak_rss_bytes'],wall_seconds=r['wall_seconds'],selector_checks=sum(int(x['selector_checks']) for x in rows(d/'t16b-ranges.csv'))))
    save('ranges',ranges);save('source-summary',sources);save('puncture-limits',limits);save('constraints',constraints);save('interface-constraints',interfaces);save('run-audits',audits)
    differences=[]
    for g in [*constraints,*interfaces]:
        if g['mode']!='geometric':continue
        s=next(x for x in [*constraints,*interfaces] if x['rung']==g['rung'] and x['mode']=='sqrt' and x['region']==g['region'] and x['subset']==g['subset'] and x['field']==g['field'])
        assert g['cells']==s['cells'] and g['peak']==s['peak'] and g['RMS']==s['RMS']
    for s in interfaces:
        if s['run']!='coarse-geometric':continue
        fine=next((x for x in interfaces if x['run']=='fine-geometric' and x['region']==s['region'] and x['subset']==s['subset'] and x['field']==s['field']),None)
        for norm in NORMS:
            differences.append(dict(interface=s['region'][6:],kind=s['interface_kind'],subset=s['subset'],field=s['field'],norm=norm,
                coarse=s[norm],fine=fine[norm] if fine else '',fine_over_coarse=fine[norm]/s[norm] if fine and s[norm]>0 else '',
                coarse_cells=s['cells'],fine_cells=fine['cells'] if fine else 0,coarse_native_crossings=s['native_crossing_cells'],
                fine_native_crossings=fine['native_crossing_cells'] if fine else 0,status='measured' if fine and s[norm]>0 else 'zero or absent; no order'))
    save('interface-growth',differences)
    observed={x['region'][6:] for x in interfaces if x['native_crossing_cells']>0}
    save('interface-coverage',[dict(interface=k,kind=v,observed=k in observed) for k,v in expected.items()])
    return ranges,sources,constraints,differences,audits
def extract(name,hole):
    reg=next(r for r in rows(HERE/'t16b-registration.csv') if r['run']==name);h=float(reg['h_M']);center=320 if hole=='left' else 352
    vectors={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/np.sqrt(2)}
    if hole=='right':vectors={k:v*np.array([-1.,1.]) for k,v in vectors.items()}
    a.NV=vectors;cells=old.select(h);shift=int((center-336)/h);actual=cells.copy();actual[:,0]+=shift
    inp=OUT/(name+'-'+hole+'-input.bin');out=OUT/(name+'-'+hole+'-q.bin');times=[];group=[];previous=None
    with inp.open('wb') as stream:
        def consume(frames):
            lo,u=a.dense(frames);rec=np.empty((len(cells),3+49*28));rec[:,0]=h;rec[:,1]=(actual[:,1]+.5)*h;rec[:,2]=(actual[:,0]+.5)*h-center
            for j in range(-3,4):
                for i in range(-3,4):
                    v=u[actual[:,1]+j-lo[1],actual[:,0]+i-lo[0]];assert np.isfinite(v).all()
                    rec[:,3+(j+3)*7+i+3::49]=v
            rec.tofile(stream);times.append(frames[0]['meta'][0])
        for frame in a.decoder.frames(ROOT/'evolution'/name/'t13-t7-stage-L13.xz'):
            if frame['phase']!=50 or len(frame['cells'])==0:continue
            belongs=float(np.mean((frame['cells'][:,0]+.5)*h))<336
            if belongs!=(hole=='left'):continue
            time=frame['meta'][0]
            if previous is not None and time!=previous:consume(group);group=[]
            group.append(frame);previous=time
        if group:consume(group)
    assert np.array_equal(times,np.arange(66)/8192/4),(name,hole,times)
    w.measured('replay-'+name+'-'+hole,OUT,[ROOT/'replay.ex','--state',inp,out])
    (OUT/(name+'-'+hole+'-input.sha256')).write_text(hashfile(inp)+' '+str(inp)+'\n');inp.unlink()
    q=np.memmap(out,dtype='f8',mode='r',shape=(66,len(cells),180));assert np.isfinite(q).all() and np.count_nonzero(q[:,:,144])==0
    a.H=h;profiles={}
    np.savez(OUT/(name+'-'+hole+'-meta.npz'),cells=actual,times=times)
    for ray,vec in vectors.items():
        profiles[ray]={n:np.array([a.sample(cells,frame,vec,n) for frame in q]) for n in (8,10)}
        np.savez_compressed(OUT/(name+'-'+hole+'-'+ray+'.npz'),times=times,r=a.R,p8=profiles[ray][8],p10=profiles[ray][10])
    return profiles,vectors,times
def launch():
    history=[];ratios=[];puncture=[]
    for name in NAMES:
        reg=next(r for r in rows(HERE/'t16b-registration.csv') if r['run']==name);h=float(reg['h_M']);rung,mode=name.split('-')
        for hole in ('left','right'):
            p,vectors,times=extract(name,hole)
            meta=np.load(OUT/(name+'-'+hole+'-meta.npz'));cells=meta['cells'];q=np.memmap(OUT/(name+'-'+hole+'-q.bin'),mode='r',dtype='f8',shape=(66,len(cells),180));c=int((320 if hole=='left' else 352)/h)
            mask=(cells[:,0]>=c-1)&(cells[:,0]<=c)&(cells[:,1]<=1)
            for step,time in enumerate(times):
                for field,index in [('Gamma1',11),('Gamma2',12),('shift1',14),('shift2',15),('lapse',13)]:
                    puncture.append(dict(run=name,hole=hole,step=step,time_M=time,field=field,peak_change=float(np.max(abs(q[step,mask,index]-q[0,mask,index]))),KO_peak=float(np.max(abs(q[step,mask,56+index])))))
            for ray,vec in vectors.items():
                fs={n:[a.fields(frame,vec) for frame in p[ray][n]] for n in (8,10)}
                for step,time in enumerate(times):
                    for field in FIELDS:
                        u=fs[8][step][field]-fs[8][0][field];control=fs[10][step][field]-fs[10][0][field]
                        for norm in NORMS:
                            amp=fn(u,norm);spread=fn(u-control,norm)
                            history.append(dict(run=name,rung=rung,mode=mode,hole=hole,ray=ray,field=field,norm=norm,step=step,time_M=time,
                                amplitude=amp,spread=spread,margin_5x=amp/(5*spread) if spread else (math.inf if amp else 0)))
    lookup={(r['run'],r['hole'],r['ray'],r['field'],r['norm'],r['step']):r for r in history}
    for g in history:
        if g['mode']!='geometric' or g['step']==0:continue
        s=lookup[g['rung']+'-sqrt',g['hole'],g['ray'],g['field'],g['norm'],g['step']];G=g['amplitude'];S=s['amplitude'];eg=5*g['spread'];es=5*s['spread']
        ratios.append(dict(rung=g['rung'],hole=g['hole'],ray=g['ray'],field=g['field'],norm=g['norm'],step=g['step'],time_M=g['time_M'],
            geometric_amplitude=G,sqrt_amplitude=S,ratio=G/S if S else '',lower_5x=max(0,G-eg)/(S+es) if S+es else '',
            upper_5x=(G+eg)/(S-es) if S>es else math.inf,geometric_margin_5x=g['margin_5x'],sqrt_margin_5x=s['margin_5x'],
            reduction_margin_5x=(S-G)/(eg+es) if eg+es else '',qualified=G>eg and S>es,endpoint=g['step']==65))
    summary=[]
    for rung in ('coarse','fine'):
        for hole in ('left','right'):
            for ray in ('axis','diagonal'):
                for field in FIELDS:
                    for norm in NORMS:
                        rr=[r for r in ratios if (r['rung'],r['hole'],r['ray'],r['field'],r['norm'])==(rung,hole,ray,field,norm)]
                        summary.append(dict(rung=rung,hole=hole,ray=ray,field=field,norm=norm,positive_clocks=65,qualified=sum(r['qualified'] for r in rr),
                            uncertainty_dominated=sum(not r['qualified'] for r in rr),central_reductions=sum(r['ratio']<1 for r in rr),
                            resolved_reductions=sum(r['upper_5x']<1 for r in rr),ratio_min=min(r['ratio'] for r in rr),ratio_max=max(r['ratio'] for r in rr),
                            enhanced_steps=' '.join(str(r['step']) for r in rr if r['ratio']>1),min_signal_margin_5x=min(min(r['geometric_margin_5x'],r['sqrt_margin_5x']) for r in rr)))
    save('launch-history',history);save('launch-ratios',ratios);save('launch-endpoint',[r for r in ratios if r['endpoint']]);save('history-summary',summary);save('native-puncture',puncture)
    return ratios,summary
def table(header,data):return '\n'.join(['| '+' | '.join(header)+' |','| '+' | '.join(['---']*len(header))+' |']+['| '+' | '.join(map(str,r))+' |' for r in data])
def fmt(v):return f'{float(v):.5g}' if v!='' else 'undefined'
def verify():
    """Check completed receipts and recompute the registered cached comparisons."""
    checks=[]
    for d in [ROOT/'pipeline',ROOT/'pipeline/job',ROOT/'evolution']:
        assert (d/'done.exit').read_text().strip()=='0'
    for f in sorted(ROOT.rglob('*.resources.json')):
        r=json.loads(f.read_text());assert r['peak_rss_bytes']<6e9 and r['tree_peak_bytes']<6.5e9
        if 'repeat' not in f.name:
            assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason'],f
    for r in json.loads((HERE/'t16b-inputs.json').read_text()):
        assert hashfile(Path(r['path']))==r['sha256']
    results=initial()
    history=rows(HERE/'t16b-launch-history.csv')
    lookup={(r['run'],r['hole'],r['ray'],r['field'],r['norm'],int(r['step'])):r for r in history}
    assert len(history)==8448 and len(lookup)==len(history)
    initial_bits={};profiles_checked=0
    for name in NAMES:
        d=ROOT/'evolution'/name;reg=next(r for r in rows(HERE/'t16b-registration.csv') if r['run']==name);h=float(reg['h_M'])
        source=rows(d/'t16b-source.csv')
        residual=max(abs(sum(float(r[f'Gamma_term_{k}']) for k in range(12))-float(r['physical_Gamma_n']))/
            max(1.,sum(abs(float(r[f'Gamma_term_{k}'])) for k in range(12))) for r in source)
        assert residual<16*np.finfo(float).eps
        assert all(float(r['Gamma_term_12'])==float(r['KO_Gamma_n']) for r in source)
        audit=re.findall(r'Geometric lapse 4x4 audit: samples=(\d+) max_relative=([\deE.+-]+)',(d/'run.log').read_text())
        assert audit and all(int(n)==24 and float(e)<16*np.finfo(float).eps for n,e in audit)
        for hole in ('left','right'):
            center=320 if hole=='left' else 352;meta=np.load(OUT/(name+'-'+hole+'-meta.npz'));cells=meta['cells'];times=meta['times']
            assert np.array_equal(times,np.arange(66)/8192/4)
            path=OUT/(name+'-'+hole+'-q.bin');assert path.stat().st_size==66*len(cells)*180*8
            q=np.memmap(path,dtype='f8',mode='r',shape=(66,len(cells),180))
            assert np.isfinite(q).all() and np.count_nonzero(q[:,:,144])==0
            keep=[k for k in range(28) if k!=13];initial_bits[name,hole]=(cells.copy(),q[0,:,keep].copy().view('u8'))
            local=cells.copy();local[:,0]-=round((center-336)/h);a.H=h
            vectors={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/np.sqrt(2)}
            if hole=='right':vectors={k:v*np.array([-1.,1.]) for k,v in vectors.items()}
            for ray,vec in vectors.items():
                z=np.load(OUT/(name+'-'+hole+'-'+ray+'.npz'))
                assert np.array_equal(z['times'],times) and np.array_equal(z['r'],a.R)
                p={n:z[f'p{n}'] for n in (8,10)}
                for n in (8,10):
                    assert np.isfinite(p[n]).all()
                    for step in (0,1,65):assert np.array_equal(a.sample(local,q[step],vec,n),p[n][step])
                fs={n:[a.fields(frame,vec) for frame in p[n]] for n in (8,10)}
                for step in range(66):
                    for field in FIELDS:
                        u=fs[8][step][field]-fs[8][0][field];control=fs[10][step][field]-fs[10][0][field]
                        for norm in NORMS:
                            r=lookup[name,hole,ray,field,norm,step]
                            assert fn(u,norm)==float(r['amplitude']) and fn(u-control,norm)==float(r['spread'])
                profiles_checked+=1
        checks.append(dict(run=name,result='PASS',common_clocks_each_hole=66,native_Base_bit_mismatches=0,
            native_Gamma_term_relative_residual=residual,geometric_4x4_max_relative=max(float(e) for n,e in audit),
            constraints_between_options='bit-identical',floors=0))
    for rung in ('coarse','fine'):
        assert hashfile(ROOT/'evolution'/(rung+'-sqrt')/'t16b-t0-native.csv')==hashfile(ROOT/'evolution'/(rung+'-geometric')/'t16b-t0-native.csv')
        for hole in ('left','right'):
            s=initial_bits[rung+'-sqrt',hole];g=initial_bits[rung+'-geometric',hole]
            assert np.array_equal(s[0],g[0]) and np.array_equal(s[1],g[1])
    ratios=rows(HERE/'t16b-launch-ratios.csv');assert len(ratios)==4160
    for r in ratios:
        for k in ('qualified','endpoint'):r[k]=r[k]=='True'
        for k in set(r)-{'rung','hole','ray','field','norm','qualified','endpoint'}:r[k]=float(r[k])
        g=lookup[r['rung']+'-geometric',r['hole'],r['ray'],r['field'],r['norm'],int(r['step'])]
        s=lookup[r['rung']+'-sqrt',r['hole'],r['ray'],r['field'],r['norm'],int(r['step'])]
        G=float(g['amplitude']);S=float(s['amplitude']);eg=5*float(g['spread']);es=5*float(s['spread'])
        assert r['ratio']==G/S and r['lower_5x']==max(0,G-eg)/(S+es) and r['upper_5x']==(G+eg)/(S-es)
        assert r['qualified']==(G>eg and S>es) and r['endpoint']==(r['step']==65)
    summary=rows(HERE/'t16b-history-summary.csv')
    for r in summary:
        rr=[v for v in ratios if all(v[k]==r[k] for k in ('rung','hole','ray','field','norm'))]
        for k in ('qualified','uncertainty_dominated'):r[k]=int(r[k])
        assert len(rr)==65 and r['qualified']==sum(v['qualified'] for v in rr)
        assert r['uncertainty_dominated']==sum(not v['qualified'] for v in rr)
        assert int(r['resolved_reductions'])==sum(v['upper_5x']<1 for v in rr)
    save('completion-checks',checks);hotspots();report(*results,ratios,summary)
    print('Verified raw t=0 audits, 16 profile caches, 8448 history amplitudes and 4160 ratio intervals.',flush=True)
def hotspots():
    """Locate extrema on the already registered sheets; masks are unchanged."""
    expected={(r['run'],r['region'],r['subset'],r['field']):r for r in rows(HERE/'t16b-interface-constraints.csv')}
    found=[]
    for name in ('coarse-geometric','fine-geometric'):
        stats={}
        with (ROOT/'evolution'/name/'t16b-t0-native.csv').open() as f:
            for r in csv.DictReader(f):
                if not r['sheets']:continue
                weight=float(r['weight'])
                for pair in r['sheets'].split(';'):
                    for subset in ('all','AMR-face-excluded'):
                        if subset!='all' and r['amr_seam']=='1':continue
                        for field in ('Ham','Mom','GaussE','GaussB','C_Gamma'):
                            key=(name,'sheet-'+pair,subset,field);v=abs(float(r[field]))
                            if key not in stats:stats[key]=[0,0.,0.,-1.,None]
                            z=stats[key];z[0]+=1;z[1]+=weight;z[2]+=weight*v*v
                            if v>z[3]:z[3]=v;z[4]=r
        for key,z in stats.items():
            e=expected[key];r=z[4];rms=math.sqrt(z[2]/z[1])
            assert z[0]==int(e['cells']) and z[3]==float(e['peak']) and math.isclose(rms,float(e['RMS']),rel_tol=2e-13,abs_tol=0.)
            X=float(r['x']);Y=float(r['y']);near=min(math.hypot(X-16,Y),math.hypot(X+16,Y));h=float(r['h'])
            found.append(dict(run=name,interface=key[1][6:],subset=key[2],field=key[3],peak=z[3],RMS=rms,
                x_M=X,y_M=Y,h_M=h,r_nearest_hole_M=near,puncture_within_stencil=near<=3*math.sqrt(2)*h,
                outer_boundary_within_stencil=min(336-abs(X),336-Y)<=3*h,native_crossing=r['native_crossing'],cells=z[0]))
    save('interface-hotspots',found)
def refresh_report():
    """Refresh prose from verified artifacts without repeating native extraction."""
    spots=rows(HERE/'t16b-interface-hotspots.csv')
    for r in spots:
        r['outer_boundary_within_stencil']=min(336-abs(float(r['x_M'])),336-float(r['y_M']))<=3*float(r['h_M'])
        if r['interface']=='21:22' and r['field']=='Ham' and r['subset']=='AMR-face-excluded':assert r['outer_boundary_within_stencil']
        if r['interface'] in ('5:10','10:15','15:20') and r['field']=='Ham' and r['subset']=='AMR-face-excluded':
            assert not r['outer_boundary_within_stencil'] and r['puncture_within_stencil']=='False'
    save('interface-hotspots',spots)
    source=rows(HERE/'t16b-source-summary.csv')
    for r in source:
        for k in ('peak','RMS'):r[k]=float(r[k])
    ratios=rows(HERE/'t16b-launch-ratios.csv')
    for r in ratios:
        for k in ('qualified','endpoint'):r[k]=r[k]=='True'
        for k in set(r)-{'rung','hole','ray','field','norm','qualified','endpoint'}:r[k]=float(r[k])
    summary=rows(HERE/'t16b-history-summary.csv')
    for r in summary:
        for k in ('qualified','uncertainty_dominated'):r[k]=int(r[k])
    report(rows(HERE/'t16b-ranges.csv'),source,rows(HERE/'t16b-constraints.csv'),rows(HERE/'t16b-interface-growth.csv'),
        rows(HERE/'t16b-run-audits.csv'),ratios,summary)
def report(ranges,sources,constraints,growth,audits,ratios,summary):
    endpoint=[r for r in ratios if r['endpoint']];data=[]
    for rung in ('coarse','fine'):
        for hole in ('left','right'):
            for ray in ('axis','diagonal'):
                for field in FIELDS:
                    rr=[next(r for r in endpoint if (r['rung'],r['hole'],r['ray'],r['field'],r['norm'])==(rung,hole,ray,field,n)) for n in NORMS]
                    data.append([rung,hole,ray,field]+[fmt(r['ratio']) for r in rr]+[fmt(min(min(r['geometric_margin_5x'],r['sqrt_margin_5x']) for r in rr)),fmt(min(r['reduction_margin_5x'] for r in rr))])
    source_data=[]
    for rung in ('coarse','fine'):
        for hole in ('left','right'):
            vals=[next(r for r in sources if r['run']==rung+'-'+mode and r['hole']==hole and r['field']==field) for mode in ('sqrt','geometric') for field in ('physical_Gamma_n','KO_Gamma_n','total_Gamma_n')]
            source_data.append([rung,hole]+[fmt(v['peak']) for v in vals])
    constraint_data=[]
    for region in ('collar-left','collar-right','interhole','far','elsewhere-off-sheet'):
        for field in ('Ham','Mom','GaussE','GaussB','C_Gamma'):
            rr=[next(r for r in constraints if r['run']==rung+'-geometric' and r['region']==region and r['field']==field and r['subset']=='all') for rung in ('coarse','fine')]
            constraint_data.append([region,field]+[fmt(r[n]) for r in rr for n in NORMS])
    interface_data=[];interface_other_data=[]
    for pair in sorted(set(r['interface'] for r in growth if r['kind']!='multiple-interface corner transition'),key=lambda s:tuple(map(int,s.split(':')))):
        rr=[r for r in growth if r['interface']==pair and r['field']=='Ham']
        vals=[]
        for subset in ('all','AMR-face-excluded'):
            for norm in NORMS:
                r=next((r for r in rr if r['subset']==subset and r['norm']==norm),None)
                vals.append(' / '.join(fmt(r[k]) for k in ('coarse','fine','fine_over_coarse')) if r else 'absent')
        interface_data.append([pair,rr[0]['kind']]+vals)
        vals=[]
        for field in ('Mom','GaussE'):
            for norm in NORMS:
                r=next((r for r in growth if r['interface']==pair and r['subset']=='AMR-face-excluded' and r['field']==field and r['norm']==norm),None)
                vals.append(' / '.join(fmt(r[k]) for k in ('coarse','fine','fine_over_coarse')) if r else 'absent')
        interface_other_data.append([pair]+vals)
    passes=all(r['qualified'] and r['upper_5x']<1 for r in endpoint if r['field']=='Gamma')
    passes &= all(next(r for r in sources if r['run']==rung+'-geometric' and r['hole']==hole and r['field']=='total_Gamma_n')['peak']<
        next(r for r in sources if r['run']==rung+'-sqrt' and r['hole']==hole and r['field']=='total_Gamma_n')['peak'] for rung in ('coarse','fine') for hole in ('left','right'))
    recommendation='ems_use_geometric_initial_lapse=true' if passes else 'ems_use_geometric_initial_lapse=false (default sqrt(chi))'
    text=['### Completed diagnostic binary evidence',
        '**READY-EXCEPT: native binary audit and launch comparison complete; the companion remains documented diagnostic data, not convergence-certified.** Every leg has its own actual-child receipt, clean native stop, clear gates and finite bit-matched native replay. Geometric runs compare all native non-lapse values including ghosts directly against the default stream before advancing; zero mismatches. Native t=0 constraints also match exactly between lapse variants. [Run receipts/floors](t16b-run-audits.csv), [ranges](t16b-ranges.csv), [limits](t16b-puncture-limits.csv), [coarse default raw source terms](/private/tmp/ems-t16b/evolution/coarse-sqrt/t16b-source.csv). Filled-ghost differences between the native lapse and a fresh geometric setter evaluation remain recorded in each run’s raw ranges table.',
        'The resumed verification checks each native leg and enclosing pipeline against its own actual-child receipt, including the sandbox timer fallback; it recomputes all fixed native t=0 norms and cached launch amplitudes/spreads and checks the ratio intervals. The pre-refresh manifest has 253 entries with zero hash mismatches. [Completion checks](t16b-completion-checks.csv) record the native 13-term source attribution and independent resampling at t=0, the first clock and the endpoint. The closed-form/4x4 inversion route agrees to a printed maximum relative error 8.88178e-16 at 24 samples for every invocation. At both punctures alpha_bin/a_own is 0.9998919207 at R=0.01, 0.9999999983 at R=1e-4, and rounds to exactly 1 by R=1e-7; the smallest tested R is 1e-10. This is a sampled puncture-limit check, not a binary source-cancellation identity.',
        table(['Run','Quantity','Minimum','Maximum'],[[r['run'],r['quantity'],fmt(r['minimum']),fmt(r['maximum'])] for r in ranges]),
        'Native t=0 longitudinal Gamma source at the four innermost cells, peak norms. Gauge RHS remains the native gauge and no binary translation identity is imposed. Physical RHS, KO and their actual sum are separate; all Base RHS bit comparisons and initialized-driver cancellations are zero.',
        table(['Rung','Hole','sqrt physical','sqrt KO','sqrt total','geometric physical','geometric KO','geometric total'],source_data),
        'Constraints on the final CTT composite state use native derivatives and the registered fixed masks. GaussB zero rows retain zero, without invented orders.',
        table(['Region','Constraint','Coarse peak','Coarse RMS','Fine peak','Fine RMS'],constraint_data),
        'Endpoint ratios are geometric/sqrt(chi) at 0.001983642578125 M_i. Sampling sensitivity intervals use five times |I8-I10|: max(0,G-eG)/(S+eS) to (G+eG)/(S-eS), unbounded if the denominator fails. Reduction margin is (S-G)/(eG+eS). These are declared sampling sensitivity estimates, not rigorous interpolation-error bounds. Full intervals and every positive clock remain in [ratios](t16b-launch-ratios.csv).',
        table(['Rung','Hole','Ray','Field','Peak ratio','RMS ratio','Min signal margin 5x','Min reduction margin 5x'],data),
        'Across the complete history, '+str(sum(r['qualified'] for r in summary))+'/'+str(len(summary)*65)+' entries qualify; '+str(sum(r['uncertainty_dominated'] for r in summary))+' are uncertainty-dominated. Early enhanced-step lists and peak/RMS extrema for each hole/ray/field/rung are retained in [history summary](t16b-history-summary.csv), together with [native puncture changes and KO peaks](t16b-native-puncture.csv). No early row is dropped. At t=0 disturbances are zero and ratios undefined. No same-grid dt/2 temporal control was registered; changing h and dt together cannot bound temporal error.',
        'Evolved Gamma, metric Gamma and lapse have resolved reductions at all 65 positive clocks on both rays, holes and rungs. Shift has a resolved early enhancement: axis peak at steps 1–14 and RMS at steps 1–16, diagonal peak at steps 1–9 and RMS at steps 1–11, the same bands on both rungs and holes. The largest shift ratio is 8.0543. Thus endpoint reduction is not uniform suppression of every gauge disturbance during the launch.',
        'Interface sheets use the same coarse-h physical mask on both rungs. Each cell below is coarse / fine / fine-to-coarse ratio for Ham. The complete [interface-growth table](t16b-interface-growth.csv) includes Mom, GaussE, GaussB and C_Gamma, actual stencil crossings and counts. [Coverage](t16b-interface-coverage.csv) lists geometric interfaces not crossed by an actual native stencil. Non-adjacent panel-pair corner transitions are kept separately in the complete data and are not labelled a single physical interface. Full sheets and AMR-face-excluded subsets distinguish the companion seams from AMR transfers; absence on a subset stays labelled.',
        table(['Panels','Interface','Full peak','Full RMS','AMR-face-excluded peak','AMR-face-excluded RMS'],interface_data),
        'Mom and GaussE below use the same AMR-face-excluded sheets; each entry is again coarse / fine / ratio. GaussB and C_Gamma vanish exactly on every t=0 sheet in both options; their ratios/orders are undefined.',
        table(['Panels','Mom peak','Mom RMS','GaussE peak','GaussE RMS'],interface_other_data),
        'The diagnostic resolves native exterior sheet noncontraction: after excluding AMR faces, Ham peak/RMS grows by 3.5675/2.7453 on the left bridge angular seam (5:10), 3.5667/2.7454 on its right counterpart (15:20), and 3.9795/2.6965 on the bridge-plane seam (10:15). The corresponding extrema are approximately 14–63 M_i from the nearer hole and far from the physical outer boundary. The outer sphere seams (5:21 and 20:22) grow by about 2.79/1.51. These are native t=0 values, so the lapse choice cannot remove them. GaussE on those sheets decreases by about 1/16; Ham does not follow fourth-order reduction. This is consistent with the documented interface regularity failure, while two grids do not establish its asymptotic rate or exclusive cause. The compact angular seam (21:22) grows by 2.0021/1.4163, but its peak lies in the last y cell at the physical boundary, y=335.5/335.75 M_i; it is not independent evidence for a companion-only defect. Near R=1e-4 the sheet Mom peak grows from 5.8753e4 to 4.7798e5 (8.1355x); its stencil includes the puncture, so this large inner value also cannot be attributed solely to a panel interface. [Hotspot locations](t16b-interface-hotspots.csv) retain both flags without changing any mask. The bridge-sheet growth persists away from both punctures, AMR faces and outer boundaries. Thirty of the 44 geometric interface pairs are observed; the 14 unobserved pairs are the unresolved end panels below the native puncture spacing. The known companion cannot enter a design-order convergence campaign on this evidence.',
        '**Low-resolution diagnostic merger choice: `'+recommendation+'`.** This choice requires resolved evolved-Gamma endpoint reduction on both rays at both holes/rungs and a smaller native total puncture Gamma source at both holes/rungs; the data above determine it. It is a lapse recommendation for diagnostic evolution only. The known second-derivative interface mismatch remains a failure of convergence admission even if these coarse native norms decrease. Do not use this companion for a design-order merger convergence campaign until its interface regularity/constraint defect is repaired and certified. No merger run is launched here.']
    path=HERE/'README.md';before=path.read_text();marker='### Completed diagnostic binary evidence';start=before.index('## T16b —');head=before[:start];section=before[start:]
    if marker in section:section=section[:section.index(marker)]
    path.write_text(head+section.rstrip()+'\n\n'+'\n\n'.join(text)+'\n')
    (HERE/'t16b-status.json').write_text(json.dumps(dict(status='READY-EXCEPT',recommendation=recommendation,companion='documented diagnostic only; known D2 mismatch',
        qualified=sum(r['qualified'] for r in summary),uncertainty_dominated=sum(r['uncertainty_dominated'] for r in summary),
        resolved_reductions=sum(int(r['resolved_reductions']) for r in summary),convergence_admission=False,
        peak_RSS_bytes=max(int(r['peak_RSS_bytes']) for r in audits)),indent=2)+'\n')
def manifest():
    resource=[]
    for f in ROOT.rglob('*.resources.json'):
        r=json.loads(f.read_text());resource.append(dict(process=r['process'],peak_RSS_bytes=r['peak_rss_bytes'],
            peak_footprint_bytes=max((v[1] for v in r.get('per_pid_peaks',{}).values()),default=0),
            tree_peak_bytes=r['tree_peak_bytes'],returncode=r['returncode'],actual_child_returncode=r.get('child_measurement',{}).get('returncode',''),gate_reason=r['gate_reason'],record=str(f)))
    save('resources',resource)
    targets=list(HERE.glob('t16b-*'))+[HERE/'README.md',HERE/'t13-run.py',HERE/'t16-native.hpp',HERE/'t13-analyze.py',HERE/'t16-analyze.py',HERE/'t7-check.py']
    targets+=list(ROOT.rglob('*'))+[w.DATA,w.PROFILE,w.DOC/'params-initial-data.txt',w.DOC/'audit.md']
    active=os.environ.get('T13_CURRENT_MEASURE','')
    lines=['T16b baseline 31176be; no commit; companion diagnostic only; OMP=2 BLAS=1; 6 GB process cap; peak RSS '+str(max(r['peak_RSS_bytes'] for r in resource))+
        '; sampled tree peak '+str(max(r['tree_peak_bytes'] for r in resource))+'; active manifest writer receipts/log excluded until completion: '+active,'SHA256 bytes path']
    for f in sorted(set(targets)):
        # Finished run logs are pinned; the active writer cannot pin itself.
        live=active and f.name.startswith(active+'.') and ROOT in f.parents
        if f.is_file() and f!=HERE/'t16b-manifest.txt' and not live:lines.append(f'{hashfile(f)} {f.stat().st_size} {f}')
    (HERE/'t16b-manifest.txt').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':
    if len(sys.argv)==1:
        results=initial();ratios,summary=launch();report(*results,ratios,summary)
    elif sys.argv[1]=='verify':verify()
    elif sys.argv[1]=='report':refresh_report()
    manifest()
