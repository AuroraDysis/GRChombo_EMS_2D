#!/usr/bin/env python3
"""Bounded analysis of current native arrays; no static fields are evaluated."""
import sys
sys.dont_write_bytecode=True
import csv, json, math, os
from pathlib import Path
import h5py
import numpy as np
import importlib.util
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t17');OUT=ROOT/'analysis';OUT.mkdir(exist_ok=True)
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
w=module('work',HERE/'t17-work.py');a=module('native_analysis',HERE/'t13-analyze.py')
FIELDS=('Gamma','metric_Gamma','shift','lapse');NORMS=('peak','RMS');CLOCK=.875/2048/4
def rows(path):return list(csv.DictReader(path.open()))
def save(name,data):
    if data:w.save('t17-'+name+'.csv',data)
def receipt(directory):
    r=json.loads((directory/(directory.name+'.resources.json')).read_text())
    assert (directory/'done.exit').read_text().strip()=='0'
    assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason']
    assert r['peak_rss_bytes']<6e9 and r['tree_peak_bytes']<6.5e9
    return r
def blocks(file):
    for l in range(int(file.attrs['num_levels'])):
        g=file[f'level_{l}'];h=float(g.attrs['dx']);ghost=g['data_attributes'].attrs['ghost'];gx,gy=map(int,ghost)
        offsets=g['data:offsets=0'][:]
        for k,raw in enumerate(g['boxes'][:]):
            box=tuple(int(raw[n]) for n in ('lo_i','lo_j','hi_i','hi_j'));nx=box[2]-box[0]+1+2*gx;ny=box[3]-box[1]+1+2*gy
            u=g['data:datatype=0'][int(offsets[k]):int(offsets[k])+28*nx*ny].reshape(28,ny,nx)
            yield l,k,box,h,u
def census():
    choices=json.loads((ROOT/'choices.json').read_text());data=[];faces=[];identity=[];coverage=[]
    for rung,choice in choices.items():
        geo=Path(choice['geometric_checkpoint']);sqrt=Path(choice['sqrt_checkpoint']);maximum=12 if rung=='coarse' else 13
        with h5py.File(geo) as f,h5py.File(sqrt) as s:
            assert int(f.attrs['num_levels'])==maximum+1 and float(f.attrs['time'])==0
            count=bits=0
            for (l,k,b,h,u),(ll,kk,bb,hh,v) in zip(blocks(f),blocks(s),strict=True):
                assert (l,k,b,h)==(ll,kk,bb,hh);assert np.isfinite(u).all() and np.isfinite(v).all()
                nonlapse=[i for i in range(28) if i!=13];count+=u[nonlapse].size
                bits+=int(np.count_nonzero(u[nonlapse].view('u8')!=v[nonlapse].view('u8')))
                chi,alpha=u[0],u[13];d=u[1]*u[3]-u[2]*u[2]
                assert (d>0).all() and (u[1]>0).all() and (u[4]>0).all()
                inverse=np.maximum(1/u[4],(u[1]+u[3]+np.hypot(u[1]-u[3],2*u[2]))/(2*d))
                beta=np.hypot(u[14],u[15]);vg=beta+np.sqrt(inverse)
                vl=beta+alpha*np.sqrt(chi*inverse);va=beta+np.sqrt(1.8*alpha*chi*inverse)
                cf=int(np.count_nonzero(chi<1e-12));af=int(np.count_nonzero(alpha<1e-12));assert cf==af==0
                assert 0<float(alpha.min())<=float(alpha.max())<=1
                data.append(dict(rung=rung,level=l,box=k,x0=b[0],y0=b[1],x1=b[2],y1=b[3],h_M=h,
                    valid_cells=(b[2]-b[0]+1)*(b[3]-b[1]+1),ghost_cells=u[0].size-(b[2]-b[0]+1)*(b[3]-b[1]+1),
                    chi_min=float(chi.min()),lapse_min=float(alpha.min()),lapse_max=float(alpha.max()),
                    chi_floors=cf,lapse_floors=af,nonfinite=0,max_beta=float(beta.max()),
                    max_light_speed=float(vl.max()),max_lapse_speed=float(va.max()),max_shift_speed=float(vg.max())))
            assert count>0 and bits==0;identity.append(dict(rung=rung,nonlapse_values=count,bit_mismatches=bits))
        for l in range(maximum+1):
            group=[r for r in data if r['rung']==rung and r['level']==l];h=group[0]['h_M']
            for hole,center in [('left',2224.),('right',2256.)]:
                axis=sorted((r['x0']*h,(r['x1']+1)*h) for r in group if r['y0']==0);merged=[]
                for lo,hi in axis:
                    if merged and lo<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],hi)
                    else:merged.append([lo,hi])
                interval=next((q for q in merged if q[0]<=center<q[1]),None);assert interval is not None
                vertical=max((r['y1']+1)*h for r in group if r['x0']*h<=center<(r['x1']+1)*h)
                faces.append(dict(rung=rung,hole=hole,level=l,h_M=h,x_left_M=interval[0],x_right_M=interval[1],
                    left_face_M=center-interval[0],right_face_M=interval[1]-center,y_face_M=vertical,
                    merged_with_other_hole=interval[0]<=2224. and interval[1]>2256.))
                if l==maximum:
                    speed=max(r['max_shift_speed'] for r in group)
                    for ray,normal in [('axis',1.),('diagonal',math.sqrt(2))]:
                        gap=min(center-interval[0],interval[1]-center,vertical)-.0025/normal-7.5*h
                        coverage.append(dict(rung=rung,hole=hole,ray=ray,gap_M=gap,speed_initial=speed,
                            travel_time_bound_M=gap/speed,margin_after_launch_M=gap/speed-.002))
                        assert gap/speed>.002
    assert max(r[k] for r in data for k in ('max_light_speed','max_lapse_speed','max_shift_speed'))<2
    save('census',data);save('faces',faces);save('initial-identity',identity);save('window-margins',coverage)
def select(h,center,vectors):
    cells=set()
    for vec in vectors.values():
        z=np.c_[(center+a.R*vec[0])/h-.5,a.R*vec[1]/h-.5]
        for start in np.floor(z).astype(int)-4:
            for j in range(10):
                y=start[1]+j;y=y if y>=0 else -y-1
                for i in range(10):cells.add((start[0]+i,y))
    c=math.floor(center/h-.5)
    for j in range(6):
        for i in range(c-5,c+7):cells.add((i,j))
    return np.array(sorted(cells,key=lambda x:(x[1],x[0])),int)
def sample(cells,q,vec,n,h,center):
    lookup={tuple(iv):i for i,iv in enumerate(cells)}
    z=np.c_[(center+a.R*vec[0])/h-.5,a.R*vec[1]/h-.5];start=np.floor(z).astype(int)-n//2+1
    wx=a.weights(z[:,0],start[:,0],n);wy=a.weights(z[:,1],start[:,1],n);answer=None
    for j in range(n):
        y=start[:,1]+j;iy=np.where(y<0,-y-1,y)
        for i in range(n):
            indices=np.array([lookup[(x,yy)] for x,yy in zip(start[:,0]+i,iy)])
            v=np.array(q[np.ix_(indices,a.SELECT)]);v[np.ix_(np.flatnonzero(y<0),a.ODD_C)]*=-1
            if answer is None:anchor=v.copy();answer=anchor.copy()
            answer+=(v-anchor)*(wx[i]*wy[j])[:,None]
    return answer
def extract(name,hole):
    rung,mode=name.split('-');level=12 if rung=='coarse' else 13;h=1.75/2**level;center=2224. if hole=='left' else 2256.
    vectors={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/math.sqrt(2)}
    if hole=='right':vectors={k:v*np.array([-1.,1.]) for k,v in vectors.items()}
    cells=select(h,center,vectors);payload=OUT/(name+'-'+hole+'-input.bin');output=OUT/(name+'-'+hole+'-q.bin')
    times=[];group=[];previous=None
    with payload.open('wb') as stream:
        def consume(frames):
            lo,u=a.dense(frames);rec=np.empty((len(cells),3+49*28));rec[:,0]=h
            rec[:,1]=(cells[:,1]+.5)*h;rec[:,2]=(cells[:,0]+.5)*h-center
            for j in range(-3,4):
                for i in range(-3,4):
                    v=u[cells[:,1]+j-lo[1],cells[:,0]+i-lo[0]];assert np.isfinite(v).all()
                    rec[:,3+(j+3)*7+i+3::49]=v
            rec.tofile(stream);times.append(frames[0]['meta'][0])
        for frame in a.decoder.frames(ROOT/'evolution'/name/f't13-t7-stage-L{level}.xz'):
            if frame['phase']!=50 or not len(frame['cells']):continue
            if (np.mean((frame['cells'][:,0]+.5)*h)<2240)!=(hole=='left'):continue
            t=frame['meta'][0]
            if abs(t/CLOCK-round(t/CLOCK))>1e-10 or t>18*CLOCK:continue
            if previous is not None and t!=previous:consume(group);group=[]
            group.append(frame);previous=t
        if group:consume(group)
    assert np.array_equal(times,np.arange(19)*CLOCK),(name,hole,times)
    assert w.measured('replay-'+name+'-'+hole,OUT/('replay-'+name+'-'+hole),[ROOT/'replay.ex','--state',payload,output])==0
    (OUT/(name+'-'+hole+'-input.sha256')).write_text(w.digest(payload)+' '+str(payload)+'\n');payload.unlink()
    q=np.memmap(output,mode='r',dtype='f8',shape=(19,len(cells),180));assert np.isfinite(q).all() and np.count_nonzero(q[:,:,144])==0
    np.savez(OUT/(name+'-'+hole+'-meta.npz'),cells=cells,times=times)
    profiles={}
    for ray,vec in vectors.items():
        profiles[ray]={n:np.array([sample(cells,frame,vec,n,h,center) for frame in q]) for n in (8,10)}
        np.savez_compressed(OUT/(name+'-'+hole+'-'+ray+'.npz'),r=a.R,times=times,p8=profiles[ray][8],p10=profiles[ray][10])
    return profiles,vectors,times,q,cells,h,center
def norm(v,kind):return float(np.max(abs(v))) if kind=='peak' else a.rms(v)
def launch():
    history=[];ratios=[];sources=[];audits=[];native=[]
    for rung in ('coarse','fine'):
        for mode in ('sqrt','geometric'):
            name=rung+'-'+mode;d=ROOT/'evolution'/name;r=receipt(d);log=(d/'run.log').read_text()
            assert 'T17 native t=0 audit complete' in log and 'T13 clean native stop' in log and 'GRChombo finished.' in log
            assert 'Read EMSTRUMPET' not in log and 'EMSTRUMPET radial solution ready' not in log
            expected=int(.002/(1.75/2**(12 if rung=='coarse' else 13)/4))*(1.75/2**(12 if rung=='coarse' else 13)/4)
            assert float(rows(d/'t13-stop.csv')[0]['actual_time_M'])==expected
            floors=[x for path in d.glob('t13-floors-*.csv') for x in rows(path)];assert floors
            totals={k:sum(int(x[k]) for x in floors) for k in ('chi_activations','lapse_activations','nonfinite')}
            assert all(v==0 for v in totals.values())
            audits.append(dict(run=name,native_stop_M=expected,chi_floors=totals['chi_activations'],lapse_floors=totals['lapse_activations'],
                nonfinite=totals['nonfinite'],static_reader_calls=0,peak_RSS_bytes=r['peak_rss_bytes']))
            for hole in ('left','right'):
                ss=[x for x in rows(d/'t17-source.csv') if x['hole']==hole];assert len(ss)==4
                assert all(int(x['Base_bit_mismatches'])==0 for x in ss)
                for field in ('physical_Gamma_n','KO_Gamma_n','total_Gamma_n'):
                    values=np.array([float(x[field]) for x in ss])
                    sources.append(dict(rung=rung,mode=mode,hole=hole,field=field,peak=float(max(abs(values))),RMS=float(np.sqrt(np.mean(values**2)))))
                profiles,vectors,times,q,cells,h,center=extract(name,hole)
                c=math.floor(center/h-.5);mask=(cells[:,0]>=c)&(cells[:,0]<=c+1)&(cells[:,1]<2);assert sum(mask)==4
                for step,time in enumerate(times):
                    for field,index in [('Gamma1',11),('Gamma2',12),('shift1',14),('shift2',15),('lapse',13)]:
                        native.append(dict(run=name,hole=hole,step=step,time_M=time,field=field,
                            peak_change=float(np.max(abs(q[step,mask,index]-q[0,mask,index]))),KO_peak=float(np.max(abs(q[step,mask,56+index])))))
                for ray,vec in vectors.items():
                    fs={n:[a.fields(p,vec) for p in profiles[ray][n]] for n in (8,10)}
                    for step,time in enumerate(times):
                        for field in FIELDS:
                            u=fs[8][step][field]-fs[8][0][field];control=fs[10][step][field]-fs[10][0][field]
                            for kind in NORMS:
                                amp=norm(u,kind);spread=norm(u-control,kind)
                                history.append(dict(run=name,rung=rung,mode=mode,hole=hole,ray=ray,field=field,norm=kind,
                                    step=step,time_M=time,amplitude=amp,spread=spread,margin_5x=amp/(5*spread) if spread else (math.inf if amp else 0)))
    lookup={(x['run'],x['hole'],x['ray'],x['field'],x['norm'],x['step']):x for x in history}
    for g in history:
        if g['mode']!='geometric' or g['step']==0:continue
        s=lookup[(g['rung']+'-sqrt',g['hole'],g['ray'],g['field'],g['norm'],g['step'])]
        G,S,eg,es=g['amplitude'],s['amplitude'],5*g['spread'],5*s['spread']
        ratios.append(dict(rung=g['rung'],hole=g['hole'],ray=g['ray'],field=g['field'],norm=g['norm'],step=g['step'],time_M=g['time_M'],
            geometric=G,sqrt=S,ratio=G/S if S else math.nan,lower_5x=max(0,G-eg)/(S+es) if S+es else math.nan,
            upper_5x=(G+eg)/(S-es) if S>es else math.inf,geometric_margin_5x=g['margin_5x'],sqrt_margin_5x=s['margin_5x'],
            reduction_margin_5x=(S-G)/(eg+es) if eg+es else math.inf,qualified=G>eg and S>es,endpoint=g['step']==18))
    summary=[]
    for rung in ('coarse','fine'):
        for hole in ('left','right'):
            for ray in ('axis','diagonal'):
                for field in FIELDS:
                    for kind in NORMS:
                        rr=[x for x in ratios if (x['rung'],x['hole'],x['ray'],x['field'],x['norm'])==(rung,hole,ray,field,kind)]
                        summary.append(dict(rung=rung,hole=hole,ray=ray,field=field,norm=kind,positive_clocks=18,
                            qualified=sum(x['qualified'] for x in rr),uncertainty_dominated=sum(not x['qualified'] for x in rr),
                            resolved_reductions=sum(x['qualified'] and x['upper_5x']<1 for x in rr),
                            enhanced_steps=' '.join(str(x['step']) for x in rr if x['ratio']>1),
                            ratio_min=min(x['ratio'] for x in rr),ratio_max=max(x['ratio'] for x in rr)))
    save('launch-history',history);save('launch-ratios',ratios);save('launch-endpoint',[x for x in ratios if x['endpoint']])
    save('history-summary',summary);save('unqualified',[x for x in ratios if not x['qualified']]);save('source-summary',sources)
    save('run-audits',audits);save('native-puncture',native)
def resolution():
    endpoint=rows(HERE/'t17-launch-endpoint.csv');sources=rows(HERE/'t17-source-summary.csv');finder=[];verdict=[]
    for rung in ('coarse','fine'):
        angular=True
        for step in (0,1):
            pairs={}
            for n in (48,96):
                d=ROOT/'finder'/f'{rung}-step{step}-n{n}';r=json.loads((d/(d.name+'.resources.json')).read_text())
                raw=rows(d/'surfaces.csv');last=[x for x in raw if int(x['stage'])==2]
                passed=(r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason'] and len(last)==2
                    and all(x['status']=='FOUND' and math.isfinite(float(x['expansion_squared'])) and float(x['expansion_squared'])<=1e-12 for x in last))
                angular &= passed
                for x in raw:
                    if int(x['stage'])<0:continue
                    finder.append({**x,'rung':rung,'step':step,'N_theta':n,'expansion_RMS':math.sqrt(float(x['expansion_squared'])),
                        'actual_child_returncode':r['returncode'],'peak_RSS_bytes':r['peak_rss_bytes']})
                if passed:pairs[n]={int(x['search_index']):x for x in last}
            if len(pairs)==2:
                for hole in (0,1):
                    for key in ('A','Q'):
                        small=float(pairs[48][hole][key]);large=float(pairs[96][hole][key]);angular &= abs(small-large)/abs(large)<=1e-3
        required=[x for x in endpoint if x['rung']==rung and x['field']=='Gamma']
        launch_pass=len(required)==8 and all(x['qualified']=='True' and float(x['upper_5x'])<1 for x in required)
        for hole in ('left','right'):
            g=next(x for x in sources if x['rung']==rung and x['mode']=='geometric' and x['hole']==hole and x['field']=='total_Gamma_n')
            s=next(x for x in sources if x['rung']==rung and x['mode']=='sqrt' and x['hole']==hole and x['field']=='total_Gamma_n')
            launch_pass &= float(g['peak'])<float(s['peak'])
        verdict.append(dict(rung=rung,launch_PASS=launch_pass,horizon_PASS=bool(angular),eligible=bool(launch_pass and angular)))
    save('finder',finder)
    spatial=[]
    if all(v['horizon_PASS'] for v in verdict):
        for step in (0,1):
            for hole in (0,1):
                pair=[next(x for x in finder if x['rung']==rung and x['step']==step and x['N_theta']==96 and int(x['search_index'])==hole and int(x['stage'])==2) for rung in ('coarse','fine')]
                Ac,Qc,Af,Qf=(float(pair[0]['A']),float(pair[0]['Q']),float(pair[1]['A']),float(pair[1]['Q']))
                ec=Qc/math.sqrt(Ac/(4*math.pi));ef=Qf/math.sqrt(Af/(4*math.pi))
                spatial.append(dict(step=step,hole=hole,A_coarse=Ac,A_fine=Af,Q_coarse=Qc,Q_fine=Qf,e_coarse=ec,e_fine=ef,
                    relative_A=abs(Ac-Af)/abs(Af),relative_Q=abs(Qc-Qf)/abs(Qf),delta_e=abs(ec-ef)))
        if any(x['relative_A']>.01 or x['relative_Q']>.01 or x['delta_e']>.1 for x in spatial):
            verdict[0]['eligible']=False;verdict[0]['spatial_reference_PASS']=False
        else:verdict[0]['spatial_reference_PASS']=True
    else:verdict[0]['spatial_reference_PASS']='unavailable';verdict[0]['eligible']=False
    verdict[1]['spatial_reference_PASS']='fine rung is the spatial reference; no third rung/order claimed'
    save('horizon-grid-differences',spatial)
    if not spatial:
        with (HERE/'t17-horizon-grid-differences.csv').open('w',newline='') as f:
            csv.DictWriter(f,fieldnames=['step','hole','A_coarse','A_fine','Q_coarse','Q_fine',
                'e_coarse','e_fine','relative_A','relative_Q','delta_e']).writeheader()
    save('resolution',verdict)
    recommendation=next((x['rung'] for x in verdict if x['eligible']),None)
    (HERE/'t17-resolution.json').write_text(json.dumps(dict(recommendation=recommendation,verdicts=verdict),indent=2)+'\n')
    return recommendation
def speed_history():
    data=[]
    for rung in ('coarse','fine'):
        d=ROOT/'timing'/(rung+'-binary');receipt(d)
        for step in (1,2):
            cp=d/f'chk/EMS_{step:06}.2d.hdf5'
            with h5py.File(cp) as f:
                t=float(f.attrs['time']);assert t==step*.4375
                values=dict(max_beta=0.,max_light_speed=0.,max_lapse_speed=0.,max_shift_speed=0.,
                    endpoint_chi_floor_values=0,endpoint_lapse_floor_values=0,chi_min=math.inf,lapse_min=math.inf)
                for l,k,b,h,u in blocks(f):
                    assert float(f[f'level_{l}'].attrs['time'])==t and np.isfinite(u).all()
                    chi,alpha=u[0],u[13];det=u[1]*u[3]-u[2]*u[2]
                    assert (chi>0).all() and (alpha>0).all() and (det>0).all() and (u[4]>0).all()
                    inv=np.maximum(1/u[4],(u[1]+u[3]+np.hypot(u[1]-u[3],2*u[2]))/(2*det))
                    beta=np.hypot(u[14],u[15])
                    for key,z in [('max_beta',beta),('max_light_speed',beta+alpha*np.sqrt(chi*inv)),
                        ('max_lapse_speed',beta+np.sqrt(1.8*alpha*chi*inv)),('max_shift_speed',beta+np.sqrt(inv))]:
                        values[key]=max(values[key],float(z.max()))
                    values['chi_min']=min(values['chi_min'],float(chi.min()));values['lapse_min']=min(values['lapse_min'],float(alpha.min()))
                    values['endpoint_chi_floor_values']+=int(np.count_nonzero(chi<=1e-12))
                    values['endpoint_lapse_floor_values']+=int(np.count_nonzero(alpha<=1e-12))
                values['registered_V2_pass']=max(values[k] for k in ('max_light_speed','max_lapse_speed','max_shift_speed'))<=2
                data.append(dict(rung=rung,step=step,time_M=t,**values))
    save('speed-history',data)
def work(file):
    with h5py.File(file) as f:
        return sum(sum((int(b['hi_i'])-int(b['lo_i'])+1)*(int(b['hi_j'])-int(b['lo_j'])+1) for b in f[f'level_{l}']['boxes'][:])*2**l for l in range(int(f.attrs['num_levels'])))
def finish():
    for item in json.loads((HERE/'t17-inputs.json').read_text()):
        assert w.digest(Path(item['path']))==item['sha256'],f'Input changed: {item["path"]}'
    recommendation=resolution();speed_history();rung=recommendation or 'fine';binary=ROOT/'timing'/(rung+'-binary');single=ROOT/'timing'/(rung+'-single')
    br=rows(binary/'t17-coarse-rates.csv');sr=rows(single/'t17-coarse-rates.csv');assert len(br)==len(sr)==2
    assert [float(x['time_M']) for x in br]==[.4375,.875] and [float(x['time_M']) for x in sr]==[.4375,.875]
    rb=float(np.median([float(x['seconds']) for x in br]));rs=float(np.median([float(x['seconds']) for x in sr]));ratio=rb/rs
    wr=work(Path('/Users/auroradysis/Workspace/EMS/.data/exp-0020/E-low/plt/EMS_Plot_000000.2d.hdf5'))
    # The unfinished single-hole rate control retains box metadata, not fields.
    ws=sum(int(x['valid_cells'])*2**int(x['level']) for x in rows(single/'t17-boxes.csv'));factor=ws/wr
    assert not list((single/'chk').glob('*.hdf5')) and not list((single/'plt').glob('*.hdf5'))
    assert not list(single.glob('t13-*.xz'))
    local=[]
    for candidate in ('coarse','fine'):
        d=ROOT/'timing'/(candidate+'-binary');r=receipt(d)
        samples=rows(d/'t17-coarse-rates.csv');assert len(samples)==2
        for step,sample in enumerate(samples,1):
            local.append(dict(rung=candidate,configuration='binary',step=step,time_M=sample['time_M'],
                seconds=sample['seconds'],capture_bytes=sum(f.stat().st_size for f in d.glob('t13-*.xz')),
                peak_RSS_bytes=r['peak_rss_bytes'],capture_in_timer=True))
    for step,sample in enumerate(sr,1):
        local.append(dict(rung=rung,configuration='single',step=step,time_M=sample['time_M'],
            seconds=sample['seconds'],capture_bytes=0,peak_RSS_bytes=receipt(single)['peak_rss_bytes'],capture_in_timer=False))
    save('local-rates',local)
    cluster_13=256*ratio*factor;cluster_15=1030*ratio*factor/3.9974387120380532
    budgets=json.loads((HERE/'t17-budget.json').read_text())['domain_budgets'];cost=[]
    for budget in budgets:
        steps=budget['coarse_steps'];model=(cluster_13+cluster_15)/2;hours=steps*model/3600
        safe_steps=int(10*3600/(2*model))
        segment_steps=8*(safe_steps//8) if safe_steps>=8 else safe_steps
        cost.append(dict(rung=rung,admitted=recommendation is not None,post_common_Mf=budget['post_common_Mf'],steps=steps,
            local_binary_s=rb,local_single_s=rs,binary_single_ratio=ratio,single_cell_step_factor=factor,
            timing_basis='completed binary rates include unintended disabled capture; single has no capture/plots/checkpoints',
            clean_binary_single_ratio_available=False,
            cluster_from_13_s=cluster_13,cluster_from_15_s=cluster_15,model_s_per_step=model,
            model_evolution_node_hours=hours,planning_evolution_low_hours=hours*.5,planning_evolution_high_hours=hours*2,
            strict_finder_pairs=47 if budget['post_common_Mf']==100 else 'add continued 5Mf common probes',
            strict_finder_low_hours=47*2006/3600 if budget['post_common_Mf']==100 else '',
            strict_finder_high_hours=47*3719/3600 if budget['post_common_Mf']==100 else '',
            segment_steps=segment_steps,checkpoint_interval=8 if safe_steps>=8 else 1,
            segment_cost_PASS=segment_steps>0,max_segment_hours=12))
    save('cost',cost)
    if recommendation:
        param=HERE/'t17-merger-draft.txt';text=param.read_text();level=12 if rung=='coarse' else 13
        text=w.p.replace(text,dict(max_level=level,regrid_interval=' '.join('4' if l<6 else '1' if l==6 else '0' for l in range(level)),
            num_mass_extraction_radii=level,mass_extraction_levels=' '.join(str(l) for l in range(1,level+1)),
            mass_extraction_radii=w.p.parameters(w.study_base(level))['mass_extraction_radii']))
        choices=json.loads((ROOT/'choices.json').read_text());text=w.p.replace(text,dict(t14_dense_initial_tags=str(choices[rung]['dense']).lower()))
        text=w.p.replace(text,dict(max_steps=cost[0]['segment_steps'],checkpoint_interval=cost[0]['checkpoint_interval']))
        text=text.replace('# Fine level is a provisional placeholder until t17-resolution.json qualifies a rung.',
            f'# Local t17-resolution.json selects {rung}, level {level}; loss diagnostics still require the listed hooks.')
        param.write_text(text)
    section=['### Completed local evidence',
        f"Least expensive qualifying candidate: **{recommendation or 'NONE'}**. See [resolution](t17-resolution.csv), [horizon residuals](t17-finder.csv), [grid differences](t17-horizon-grid-differences.csv), [endpoint ratios](t17-launch-endpoint.csv), [all histories](t17-launch-ratios.csv), and [history counts](t17-history-summary.csv). A capped/failed horizon probe is not convergence.",
        f"The local median binary/single coarse-step rates at {rung} are {rb:.6g}/{rs:.6g} s, raw ratio {ratio:.6g}. Completed binary timers include unintended disabled-capture writes; the unfinished single control writes no field output. The storage bias is unmeasured, so this is not a clean compute ratio and the cluster model is a storage-inclusive planning estimate. Completed evolutions were preserved without repetition. The single-hole cell-step work is {factor:.6g} times exp-0020's actual E-low hierarchy; normalized 13/15-level anchors estimate {cluster_13:.6g}/{cluster_15:.6g} s per binary coarse step. A factor 0.5–2 planning range covers unmeasured machine/grid evolution and coarsest regrid overhead, but is not a measured error bar or a bound on the capture bias; read [cost](t17-cost.csv).",
        'The companion remains diagnostic-only. These launch/horizon controls do not repair its C1 interfaces or establish long-duration loss accuracy. The production build and write-only wrappers have separate default-path and dense identity checks. The parameter file remains a draft: the loss protocol, t=0 large-radius mass/charge diagnostic and the bounded common-search/checkpoint orchestration must be in place before any merger submission. No merger or cluster operation was run.']
    decision=json.loads((HERE/'t17-resolution.json').read_text())
    table=['| Rung | R_h/h | Launch | N48/N96 horizons | Admission |',
        '|---|---:|---|---|---|']
    for v in decision['verdicts']:
        table.append(f"| {v['rung']} | {w.RH/(1.75/2**(12 if v['rung']=='coarse' else 13)):.6f} | {'PASS' if v['launch_PASS'] else 'FAIL'} | {'PASS' if v['horizon_PASS'] else 'FAIL'} | {'eligible' if v['eligible'] else 'not admitted'} |")
    section.append('\n'.join(table))
    if recommendation is None:
        section.append('**No registered rung is admitted for the merger.** The draft fine grid is a cost reference only. Do not infer horizon convergence from area/charge values of failed searches, and do not tighten the registered stopping controls to manufacture admission.')
        param=HERE/'t17-merger-draft.txt';text=param.read_text()
        text=text.replace('# Fine level is a provisional placeholder until t17-resolution.json qualifies a rung.',
            '# NO RUNG ADMITTED: fine level is a cost reference only; registered horizon tests failed.')
        param.write_text(text)
    path=HERE/'README.md';text=path.read_text();marker='### Completed local evidence'
    start=text.index('## T17 —');head=text[:start];tail=text[start:]
    if marker in tail:tail=tail[:tail.index(marker)]
    path.write_text(head+tail.rstrip()+'\n\n'+'\n\n'.join(section)+'\n')
    manifest()
def manifest():
    resource=[]
    for f in sorted(ROOT.rglob('*.resources.json')):
        r=json.loads(f.read_text());assert r['peak_rss_bytes']<6e9
        resource.append(dict(process=r['process'],peak_RSS_bytes=r['peak_rss_bytes'],tree_peak_bytes=r['tree_peak_bytes'],
            actual_child_returncode=r.get('child_measurement',{}).get('returncode',''),returncode=r['returncode'],gate_reason=r['gate_reason'],record=str(f)))
    save('resources',resource)
    targets=list(HERE.glob('t17-*'))+[HERE/'README.md',HERE/'t13-run.py',HERE/'t13-analyze.py',HERE/'t16-native.hpp',HERE/'T14DenseTags.hpp',w.PROFILE,w.CTT]
    targets += [w.REPO/'Examples/EMS'/f for f in ('Main_EMSBH2DBH.cpp','EMSBH2DLevel.cpp','SimulationParameters.hpp')]
    targets += [w.REPO/'Source/BlackHoles/PunctureTracker.cpp']
    diagnosis=ROOT/'finder-diagnosis'
    diagnosis_complete=(diagnosis/'pipeline/done.exit').exists()
    targets += [f for f in ROOT.rglob('*') if f.is_file() and f.name not in ('launcher.log','launcher.pid')
        and ('pipeline' not in f.parts or f.name.endswith(('.resources.json','.child.json','.time')) or f.name=='done.exit')
        and (diagnosis_complete or diagnosis not in f.parents)]
    if (diagnosis/'registration.json').exists():targets.append(diagnosis/'registration.json')
    if not diagnosis_complete:
        targets=[f for f in targets if f.name not in ('t17-finder-diagnosis.csv','t17-finder-confirmation.csv',
            't17-finder-diagnosis-status.json','t17-finder-diagnosis-angular.csv','t17-finder-diagnosis-spatial.csv')]
    for extra in (Path('/private/tmp/ems-t17-interim'),Path('/private/tmp/ems-t17-final-audit')):
        if extra.exists():targets += [f for f in extra.iterdir() if f.is_file() and f.name.endswith(('.resources.json','.child.json','.time'))]
    text=['# T17 no commit; SHA256 bytes path.','\n# Peak RSS and actual child/gate receipts: t17-resources.csv.',
        '# Completed original/continuation pipeline receipts included; launcher files omitted.',
        '# Completed finder-diagnosis outputs included.' if diagnosis_complete else
        '# Active finder-diagnosis outputs excluded until its done.exit; registration and source are included.']
    for f in sorted(set(targets)):
        if f.is_file():text.append(f'{w.digest(f)} {f.stat().st_size} {f}')
    (HERE/'COMMIT-MANIFEST-T17.txt').write_text('\n'.join(text)+'\n')
def self_check():
    # Off-grid centres are essential: translating the old 336-centred integer
    # nodes would lose the native subcell phase on this binary grid.
    vectors={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/math.sqrt(2)}
    count=0
    for h in (1.75/4096,1.75/8192):
        for center in (2224.,2256.):
            cells=select(h,center,vectors);x=(cells[:,0]+.5)*h-center;y=(cells[:,1]+.5)*h
            q=np.zeros((len(cells),180));q[:,0]=.1+x*x+y*y;q[:,12]=x*y
            for vec in vectors.values():
                for n in (8,10):
                    sampled=sample(cells,q,vec,n,h,center)
                    assert np.max(abs(sampled[:,0]-(.1+a.R*a.R)))<2e-14
                    assert np.max(abs(sampled[:,12]-(a.R*a.R*vec[0]*vec[1])))<2e-15
                    count+=1
    save('analysis-check',[dict(check='I8/I10 off-grid quadratic and odd-y reflection, both actual puncture phases and spacings',cases=count,result='PASS')])
if __name__=='__main__':
    {'census':census,'launch':launch,'resolution':resolution,'finish':finish,'manifest':manifest,'self-check':self_check}[sys.argv[1]]()
