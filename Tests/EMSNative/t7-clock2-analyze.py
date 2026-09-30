#!/usr/bin/env python3
"""Frozen T7 A norm comparison; current-field recordings only, one level per call."""
import ast, csv, json, math, resource, sys, time
from collections import Counter, defaultdict
from pathlib import Path
import importlib.util
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t7a', HERE / 't7a-analyze.py')
a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
RUN = Path('/private/tmp/ems-t7-clock2/evolution/clock-2')
TMP = Path('/private/tmp/ems-t7-clock2-analysis'); TMP.mkdir(exist_ok=True)

def read(p):
    with Path(p).open() as f: return list(csv.DictReader(f))

def save(p, rows):
    assert rows, p
    with Path(p).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def packet(lev, t, h, iv, v):
    # Same physical ROI and uncovered-cell predicate as common_packet_peaks.
    x = (iv[:, 0] + .5)*h - 256; y = (iv[:, 1] + .5)*h
    e = np.maximum(abs(x), y); valid = e > a.FACES.get(lev+1, -1) + 1e-10
    rows = []
    for face in (3., 4.):
        w = 1/24 if face == 3 else 1/12
        for ray in a.RAYS:
            tube = (x > 0) & (y < w) if ray == 'axis_plus' else abs(x) < w if ray == 'equator' else (abs(x-face) < w) & (abs(y-face) < w)
            use = valid & tube & (abs(e-face) < w)
            if not use.any(): continue
            for field, k in (('Ham', 28), ('Theta', 10)):
                rows.append(dict(time_M=t, face_M=face, ray=ray, level=lev, field=field,
                                 cells=int(use.sum()), peak_abs=float(np.max(abs(v[use, k])))))
    return rows

def snapshots(lev, directory=RUN):
    start = time.monotonic(); path = directory / f't7-snapshot-L{lev}.xz'
    metrics, peaks, audit, boxes, group, ts = [], [], [], set(), [], []
    current = None
    def flush():
        m = group[0]['meta']; t, h = m[:2]
        iv = np.concatenate([f['cells'] for f in group]); v = np.concatenate([f['values'] for f in group])
        assert len(np.unique(iv, axis=0)) == len(iv)
        metrics.extend(a.metrics(directory.name, t, h, lev, iv, v)); peaks.extend(packet(lev,t,h,iv,v)); ts.append(t)
        audit.append(dict(level=lev,time_M=t,cells=len(iv),chi_min=float(v[:,0].min()),lapse_min=float(v[:,13].min()),
                          chi_floor_cells=int(np.count_nonzero(v[:,0]<=1e-12)),lapse_floor_cells=int(np.count_nonzero(v[:,13]<=1e-12)),
                          GaussB_nonzero=int(np.count_nonzero(v[:,32])),finite=True))
        for f in group:
            x0,y0,x1,y1=f['valid']; boxes.add((x0*h-256,(x1+1)*h-256,y0*h,(y1+1)*h))
    for f in a.check.frames(path):
        assert f['phase'] == 50 and f['level'] == lev and f['values'].shape[1] == 33
        t = f['meta'][0]
        if current is not None and abs(t-current)>1e-10: flush(); group=[]
        current=t; group.append(f)
    flush(); assert len(ts)==73 and np.max(abs(np.array(ts)-a.TIMES))<1e-9
    if directory != RUN: return peaks
    expected = {tuple(float(r[k]) for k in ('xlo_M','xhi_M','ylo_M','yhi_M')) for r in read(HERE/'t7-boxes.csv') if r['case']=='clock' and int(r['level'])==lev}
    assert len(boxes)==len(expected)
    assert all(any(np.max(abs(np.array(b)-e))<1e-10 for e in expected) for b in boxes)
    for name, rows in (('metrics',metrics),('packets',peaks),('snapshots',audit)):
        save(TMP/f'{name}-L{lev}.csv',rows)
    save(TMP/f'snapshot-input-L{lev}.csv',[dict(path=str(path),bytes=path.stat().st_size,sha256=a.digest(path),frames=73*len(boxes),
         finite='PASS',columns=33,box_unions_and_seams='IDENTICAL_TO_CLOCK',seconds=time.monotonic()-start,
         peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)])
    print('SNAPSHOTS',lev,len(metrics),'seconds',time.monotonic()-start,flush=True)

def stages(lev):
    start=time.monotonic(); path=RUN/f't7-stage-L{lev}.xz'; counts=Counter(); previous={}; parts={}
    floor_count=0; floor_change=0.; stage_error=0.; rhs_error=0.; chi=math.inf; lapse=math.inf
    for f in a.check.frames(path):
        p,src,m,v=f['phase'],f['source'],f['meta'],f['values'];counts[p]+=1
        assert f['level']==lev and v.shape[1]==28
        if p in (2,3,20,21,22,40):
            assert f['stage'] in (0,1,2,3)
            stage_error=max(stage_error,abs(m[3]-m[0]-(0,.5,.5,1)[f['stage']]*m[2]))
        if p in (2,3,4,5,6,7,10,11,40):
            chi=min(chi,float(v[:,0].min()));lapse=min(lapse,float(v[:,13].min()))
        if p in (2,4,6):previous[(src,p)]=f
        if p in (3,5,7):
            before=previous[(src,p-1)]; assert np.array_equal(before['cells'],f['cells'])
            b=before['values'][:,[0,13]]; z=v[:,[0,13]]
            floor_count+=int(np.count_nonzero(b<1e-12));floor_change=max(floor_change,float(np.max(abs(z-b))))
        if p in (20,21):parts[(src,p)]=f
        if p==22:
            b,c=parts[(src,20)],parts[(src,21)]
            assert np.array_equal(b['cells'],c['cells']) and np.array_equal(b['cells'],f['cells'])
            rhs_error=max(rhs_error,float(np.max(abs(v-b['values']-c['values'])/np.maximum(1.,abs(b['values'])+abs(c['values'])))))
    old=read(a.artifact(f't7a-stage-audit-clock-L{lev}.csv'))[0]
    # The same recorder operations occur twice as often at half dt.
    assert sum(counts.values())==2*int(old['frames'])
    assert dict(counts)=={int(k):2*v for k,v in ast.literal_eval(old['phase_counts']).items()}
    assert stage_error<1e-12 and rhs_error<16*np.finfo(float).eps
    save(TMP/f'stage-input-L{lev}.csv',[dict(path=str(path),bytes=path.stat().st_size,sha256=a.digest(path),frames=sum(counts.values()),
        finite='PASS',columns=28,phase_counts=json.dumps(dict(sorted(counts.items()))),stage_time_error=stage_error,
        RHS_addition_error=rhs_error,chi_min_captured=chi,lapse_min_captured=lapse,pre_projection_below_floor=floor_count,
        projection_chi_lapse_change_max=floor_change,seconds=time.monotonic()-start,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)])
    print('STAGES',lev,sum(counts.values()),'seconds',time.monotonic()-start,flush=True)

def comparison(n0,n1,n2,s):
    first=n1-n0;second=n2-n1;spatial=s-n1
    ratio=abs(second/spatial) if spatial else (0. if second==0 else math.inf)
    verdict='KILL' if ratio>=.5 else 'PASS' if abs(second/n1)<=.05 and ratio<.2 else 'IN_BETWEEN'
    return dict(T6=n0,clock=n1,clock2=n2,space=s,clock2_over_clock_minus_1=second/n1,
         space_over_clock_minus_1=spatial/n1,signed_contrast_ratio=second/spatial if spatial else math.nan,
         absolute_contrast_ratio=ratio,clock_over_T6_minus_1=first/n0,clock2_over_T6_minus_1=(n2-n0)/n0,
         first_temporal_change=first,second_temporal_change=second,
         temporal_signed_ratio=second/first if first else math.nan,
         temporal_absolute_ratio=abs(second/first) if first else math.nan,verdict=verdict)

def finish():
    native=[r for lev in (4,5,6) for r in read(TMP/f'metrics-L{lev}.csv')];groups=defaultdict(list)
    for r in native:groups[(a.snapshot_time(float(r['time_M'])),float(r['face_M']),r['ray'],r['band'],r['field'])].append(r)
    composite=[]; lookup={}
    for key,rr in sorted(groups.items()):
        w=sum(float(r['weight']) for r in rr);ss=sum(float(r['square_sum']) for r in rr)
        row=dict(time_M=key[0],face_M=key[1],ray=key[2],band=key[3],field=key[4],
                 cells=sum(int(r['cells']) for r in rr),volume_weight=w,rms=math.sqrt(ss/w))
        composite.append(row);lookup[key]=row
    save(HERE/'t7-clock2-composite-metrics.csv',composite)
    hist=[r for lev in (4,5,6) for r in read(TMP/f'packets-L{lev}.csv')]
    save(HERE/'t7-clock2-packet-history.csv',hist)
    rows=[]
    for r in read(HERE/'t7a-common-packet-ratios.csv'):
        face=float(r['face_M']);lev=int(r['level']);ray=r['ray'];field=r['field']
        if lev!=int(9-face) or ray not in ('axis_plus','diagonal'):continue
        p=max((q for q in hist if float(q['face_M'])==face and int(q['level'])==lev and q['ray']==ray and q['field']==field),key=lambda q:float(q['peak_abs']))
        row=dict(kind='packet',face_M=face,ray=ray,field=field,time_M=float(p['time_M']),level=lev,cells=int(p['cells']),
                 T6_peak_time_M=float(r['T6_peak_time_M']),clock_peak_time_M=float(r['clock_peak_time_M']),
                 clock2_peak_time_M=float(p['time_M']),space_peak_time_M=float(r['space_peak_time_M']),
                 condition='RISING_ENDPOINT_AT_6M' if face==4 and ray=='diagonal' else 'INDEPENDENT_FULL_INTERVAL_PEAK')
        row.update(comparison(*(float(r[k]) for k in ('T6_peak','clock_peak')),float(p['peak_abs']),float(r['space_peak'])));rows.append(row)
    # Correct the preceding T7 A exact-float time grouping, without changing ROIs.
    baseline_groups=defaultdict(list)
    for r in csv.DictReader((HERE/'t7a-native-metrics.csv').open()):
        t=a.snapshot_time(float(r['time_M']))
        if t in (4.,6.) and r['ray']=='all_face':
            baseline_groups[(r['case'],t,float(r['face_M']),r['band'],r['field'])].append(r)
    fragments=defaultdict(list)
    for r in csv.DictReader((HERE/'t7a-composite-metrics.csv').open()):
        t=a.snapshot_time(float(r['time_M']))
        if t in (4.,6.) and r['ray']=='all_face':
            fragments[(r['case'],t,float(r['face_M']),r['band'],r['field'])].append(r)
    corrected={}
    for key,rr in baseline_groups.items():
        assert len({int(r['level']) for r in rr})==len(rr)
        w=sum(float(r['weight']) for r in rr);ss=sum(float(r['square_sum']) for r in rr)
        corrected[key]=dict(rms=math.sqrt(ss/w),cells=sum(int(r['cells']) for r in rr),
            volume_weight=w,levels=';'.join(r['level'] for r in rr),
            raw_times_M=';'.join(r['time_M'] for r in rr))
    supplement=[];corrections=[]
    for r in read(HERE/'t7a-ratios.csv'):
        face=float(r['face_M']);t=float(r['time_M']);field=r['field'];band=r['band'];q=lookup[(t,face,'all_face',band,field)]
        primary=(field=='Ham' and band!='packet') or (band=='packet' and field in ('Ham','Theta') and (face,t) in ((3.,4.),(4.,6.)))
        row=dict(kind='collar' if band=='packet' else band,face_M=face,ray='all_face',field=field,time_M=t,level='COMPOSITE_4_5_6',cells=q['cells'],
                 T6_peak_time_M='',clock_peak_time_M='',clock2_peak_time_M='',space_peak_time_M='',
                 condition='BEFORE_MAIN_4M_CROSSING' if face==4 and t==4 and band!='packet' else 'MATCHED_PHYSICAL_VOLUME_ROI')
        b={case:corrected[(case,t,face,band,field)] for case in ('T6_face3','clock','space')}
        assert b['clock']['cells']==q['cells'] and abs(b['clock']['volume_weight']/q['volume_weight']-1)<1e-12
        correction=dict(face_M=face,time_M=t,band=band,field=field,primary=primary)
        for case in ('T6_face3','clock','space'):
            correction.update({f'{case}_{k}':v for k,v in b[case].items()})
            correction[f'{case}_published_rms']=float(r['T6_face3' if case=='T6_face3' else case])
            published=[q for q in fragments[(case,t,face,band,field)] if float(q['rms'])==correction[f'{case}_published_rms']]
            assert len(published)==1
            correction[f'{case}_published_cells']=int(published[0]['cells'])
        corrections.append(correction)
        row.update(comparison(b['T6_face3']['rms'],b['clock']['rms'],q['rms'],b['space']['rms']));supplement.append({**row,'primary':primary})
        if primary:rows.append(row)
    assert len(rows)==20
    save(HERE/'t7-clock2-norm-comparison.csv',rows);save(HERE/'t7-clock2-supplemental-norms.csv',supplement)
    save(HERE/'t7-clock2-time-grouping-correction.csv',corrections)
    for r in rows:print(r['kind'],r['face_M'],r['ray'],r['field'],round(r['time_M'],6),*(f'{r[k]:.7g}' for k in ('clock2_over_clock_minus_1','space_over_clock_minus_1','absolute_contrast_ratio','temporal_absolute_ratio')),r['verdict'])

def ray_wakes():
    sampled={}
    for lev in (4,5,6):
        group=[];current=None
        def flush():
            t=group[0]['meta'][0]
            if min(abs(t-4),abs(t-6))>1e-9:return
            iv=np.concatenate([f['cells'] for f in group]);v=np.concatenate([f['values'] for f in group]);h=group[0]['meta'][1]
            for n,mode in ((6,'values'),(8,'values8')):sampled[(lev,round(t),mode)]=a.line_sample(iv,v,h,lev,n)
        for f in a.check.frames(RUN/f't7-snapshot-L{lev}.xz'):
            t=f['meta'][0]
            if current is not None and abs(t-current)>1e-10:flush();group=[]
            current=t;group.append(f)
        flush()
    rows=[]
    for r in read(HERE/'t7a-ray-wake-ratios.csv'):
        t=round(float(r['time_M']));ri=a.RAYS.index(r['ray']);mode=r['method'];k={'Ham':28,'Theta':10,'GaussE':31}[r['field']]
        z=np.full(len(a.R),np.nan)
        for lev in (4,5,6):
            mask=(a.R*max(a.NV[ri])<a.FACES[lev])&(a.R*max(a.NV[ri])>=a.FACES.get(lev+1,0))
            z[mask]=sampled[(lev,t,mode)][ri,k,mask]
        use=(a.R>float(r['r_min_M']))&(a.R<float(r['r_max_M']));assert np.isfinite(z[use]).all()
        n2=float(np.sqrt(np.mean(z[use]**2)))
        rows.append({**{k:r[k] for k in ('time_M','ray','field','method','r_min_M','r_max_M','sampling_condition')},
                     **comparison(*(float(r[k]) for k in ('T6_rms','clock_rms')),n2,float(r['space_rms']))})
    save(HERE/'t7-clock2-ray-wake-comparison.csv',rows)

def audit():
    def params(path):
        return {k.strip():v.strip() for k,v in (line.split('=',1) for line in Path(path).read_text().splitlines() if '=' in line)}
    old=Path('/private/tmp/ems-t7/evolution/clock');p=params(RUN/'params.txt');p0=params(old/'params.txt')
    assert {k for k in p if p[k]!=p0.get(k)}=={'dt_multiplier','max_steps'} and p.keys()==p0.keys()
    assert p['sigma']=='1' and p['amr_transfer']=='point' and p['nan_check']=='1' and p['t2_guard_initial_data_after_t0']=='true'
    assert set(p['regrid_interval'].split())=={'0'} and p['min_chi']==p['min_lapse']=='1e-12'
    assert (RUN/'params.txt').read_bytes()==(HERE/'params/t7-ref-mid-clock2.txt').read_bytes()
    assert (RUN/'done.exit').read_text().strip()=='0'
    assert a.digest('/private/tmp/ems-t7/evolution.ex')=='15e9e81bb100c3bfb508637cadd2f82c0df79814f8f8f73c22b607b5199aebdf'
    text=(RUN/'run.log').read_text();baseline=(old/'run.log').read_text()
    assert 'GRChombo finished.' in text
    assert not any(s in text.lower() for s in ('nancheck in','mayday::error','segmentation fault','nonfinite','t7 output gate:'))
    defaults=[dict(message=line) for line in text.splitlines() if 'not found' in line]
    assert Counter(r['message'] for r in defaults)==Counter(line for line in baseline.splitlines() if 'not found' in line)
    save(HERE/'t7-clock2-default-parameters.csv',defaults)
    sources=[]
    for r in read(HERE/'t7a-source-controls.csv'):
        assert a.digest(HERE.parents[1]/r['path'])==r['sha256']
        sources.append({**r,'status':'UNCHANGED_IN_CLOCK2_ANALYSIS'})
    save(HERE/'t7-clock2-source-controls.csv',sources)
    inputs=[r for lev in (4,5,6) for kind in ('stage','snapshot') for r in read(TMP/f'{kind}-input-L{lev}.csv')]
    save(HERE/'t7-clock2-input-audit.csv',[{k:r[k] for k in ('path','bytes','sha256','frames','finite','columns','seconds','peak_rss_bytes')} for r in inputs])
    s=[r for lev in (4,5,6) for r in read(TMP/f'snapshots-L{lev}.csv')];st=[read(TMP/f'stage-input-L{lev}.csv')[0] for lev in (4,5,6)]
    for lev,r in zip((4,5,6),st):
        oldcounts=ast.literal_eval(read(a.artifact(f't7a-stage-audit-clock-L{lev}.csv'))[0]['phase_counts'])
        assert json.loads(r['phase_counts'])=={str(k):2*v for k,v in oldcounts.items()}
    save(HERE/'t7-clock2-stage-audit.csv',st)
    save(HERE/'t7-clock2-snapshot-audit.csv',s)
    output=sum(f.stat().st_size for f in RUN.iterdir() if f.is_file());assert output<6_000_000_000
    row=dict(case='clock-2',exit=0,sigma=1,amr_transfer='point',fixed_hierarchy=True,static_guard=True,
         stage_frames=sum(int(r['frames']) for r in st),snapshots_per_level=73,missing_evolved_components=0,
         nonfinite='NONE_IN_ALL_DECODED_FRAMES;GLOBAL_RUNTIME_NAN_CHECK_ON',chi_min_snapshot=min(float(r['chi_min']) for r in s),
         lapse_min_snapshot=min(float(r['lapse_min']) for r in s),floor_cells_snapshot=sum(int(r['chi_floor_cells'])+int(r['lapse_floor_cells']) for r in s),
         floor_hits_stage=sum(int(r['pre_projection_below_floor']) for r in st),projection_chi_lapse_change_max=max(float(r['projection_chi_lapse_change_max']) for r in st),
         GaussB_nonzero=sum(int(r['GaussB_nonzero']) for r in s),stage_time_error_max=max(float(r['stage_time_error']) for r in st),
         RHS_addition_error_max=max(float(r['RHS_addition_error']) for r in st),wall_seconds=float((RUN/'wall_seconds').read_text()),
         evolution_peak_rss_bytes=int((RUN/'peak_rss_bytes').read_text()),output_bytes=output,
         analysis_peak_rss_bytes=max(int(r['peak_rss_bytes']) for r in inputs),default_messages=len(defaults),
         default_messages_vs_clock='IDENTICAL',box_unions_and_seams='IDENTICAL_L4_L5_L6;IDENTICAL_FIXED_PARAMETERS_L0_TO_L6',
         min_chi_file_actual_variable='DIAGNOSTIC_MOD_F_NOT_EVOLVED_CHI',params_sha256=a.digest(RUN/'params.txt'),log_sha256=a.digest(RUN/'run.log'))
    save(HERE/'t7-clock2-run-audit.csv',[row]);print(json.dumps(row,indent=2))

def figures():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scienceplots
    plt.style.use(['science','no-latex']);plt.rcParams.update({'font.size':9})
    rows=read(HERE/'t7-clock2-norm-comparison.csv');labels=[]
    for r in rows:
        kind=r['kind'];label=f"{float(r['face_M']):g} M {r['ray'].replace('axis_plus','axis').replace('all_face','')} {r['field']} "
        label+=f"peak" if kind=='packet' else f"{kind.replace('wake_','')} t={float(r['time_M']):g}"
        labels.append(label.strip())
    y=np.arange(len(rows));colors=['#b2182b' if r['verdict']=='KILL' else '#2166ac' for r in rows]
    fig,ax=plt.subplots(1,3,figsize=(12,8),sharey=True,layout='constrained')
    columns=('clock2_over_clock_minus_1','absolute_contrast_ratio','temporal_absolute_ratio')
    for col,key in enumerate(columns):
        z=np.array([float(r[key]) for r in rows])*(100 if col==0 else 1)
        ax[col].scatter(z,y,c=colors,s=22);ax[col].grid(axis='both',alpha=.2)
    ax[0].axvline(0,color='.3',lw=.5)
    ax[0].axvline(5,color='.4',ls='--');ax[0].axvline(-5,color='.4',ls='--');ax[0].set_xlim(-5.5,5.5)
    ax[0].set_xlabel('100 (clock-2/clock - 1), %')
    ax[1].axvline(.2,color='.4',ls='--');ax[1].axvline(.5,color='#b2182b',ls='--');ax[1].set_xscale('log');ax[1].set_xlim(.001,1)
    ax[1].set_xlabel('|clock-2 - clock| / |space - clock|')
    ax[2].axvline(1/16,color='.4',ls='--');ax[2].set_xscale('log');ax[2].set_xlim(.04,10)
    ax[2].set_xlabel('|clock-2 - clock| / |clock - T6|')
    ax[0].set_yticks(y,labels);ax[0].invert_yaxis()
    result='KILL' if any(r['verdict']=='KILL' for r in rows) else 'PASS' if all(r['verdict']=='PASS' for r in rows) else 'IN BETWEEN'
    fig.suptitle(f'Registered clock-2 norm test: {result}, after matching level timestamps\nDashed temporal line = 1/16. 4 M corner peaks end at 6 M.')
    for ext in ('png','pdf'):fig.savefig(HERE/f'figures/t7-clock2-norms.{ext}',dpi=220)
    plt.close(fig)

def check():
    assert a.snapshot_time(4.-2e-14)==a.snapshot_time(4.+2e-14)==4.
    try:a.snapshot_time(4.01)
    except AssertionError:pass
    else:raise AssertionError('off-ladder time accepted')
    np.testing.assert_allclose(comparison(1.,1.1,1.10625,.5)['temporal_absolute_ratio'],1/16)
    assert comparison(1.,1.,1.01,.5)['verdict']=='PASS'
    assert comparison(1.,1.,1.15,.5)['verdict']=='IN_BETWEEN'
    assert comparison(1.,1.,1.25,.5)['verdict']=='KILL'
    assert comparison(1.,1.,1.,1.)['verdict']=='PASS'
    assert comparison(1.,1.,1.01,1.)['verdict']=='KILL'
    old=read(HERE/'t7a-common-packet-ratios.csv')
    for lev in (5,6):
        rows=snapshots(lev,Path('/private/tmp/ems-t7/evolution/clock'))
        for r in old:
            if int(r['level'])!=lev:continue
            peak=max((q for q in rows if q['face_M']==float(r['face_M']) and q['ray']==r['ray'] and q['field']==r['field']),key=lambda q:q['peak_abs'])
            assert peak['peak_abs']==float(r['clock_peak']) and peak['time_M']==float(r['clock_peak_time_M'])
    for r in read(HERE/'t7-clock2-time-grouping-correction.csv'):
        assert len(r['clock_levels'].split(';'))==len(set(r['clock_levels'].split(';')))
        if r['primary']=='True':assert int(r['clock_cells'])>0
        if r['primary']=='True' and r['band']=='packet':
            assert int(r['clock_cells'])>int(r['clock_published_cells'])
    for r in read(HERE/'t7-clock2-norm-comparison.csv'):
        q=comparison(*(float(r[k]) for k in ('T6','clock','clock2','space')))
        assert q['verdict']==r['verdict'] and q['absolute_contrast_ratio']==float(r['absolute_contrast_ratio'])
    print('PASS: frozen clock packet extraction reproduced exactly; strict threshold and zero-contrast checks')

if __name__=='__main__':
    task=sys.argv[1]
    if task in ('snapshots','stages'):globals()[task](int(sys.argv[2]))
    else:globals()[task]()
