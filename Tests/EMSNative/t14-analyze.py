#!/usr/bin/env python3
"""T14 registered current-state screen; reuse frozen T13 native operators."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, json, math, os, subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t14')
OLD=Path('/private/tmp/ems-t13');OUT=ROOT/'analysis';OUT.mkdir(exist_ok=True)
PY='/Users/auroradysis/miniconda3/bin/python'
LEGS={'max12':(12,1),'max13':(13,2),'max14':(14,4),'max14-half':(14,8)}
FIELDS=('Gamma','metric_Gamma','shift','lapse')

def module(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

a=module('t13analysis','t13-analyze.py')
BASE_FRAMES=a.decoder.frames

def configure(leg):
    level,stride=LEGS[leg];a.H=(7/12)/2**level;a.DT=a.H/4
    a.CENTRE=round(336/a.H);a.TMP=ROOT;a.OUT=OUT;a.LEGS=(leg,)
    a.SELECT=list(range(28))+list(range(145,152));a.ODD_C=a.ODD_FIELDS+[29,32]
    return level,stride

def save(name,rows,merge=False):
    assert rows,name
    path=HERE/name
    if merge and path.exists():
        old=list(csv.DictReader(path.open()));runs={x['run'] for x in rows}
        rows=[x for x in old if x['run'] not in runs]+rows
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def hash_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def native(label,mode,payload,output):
    subprocess.run([PY,str(HERE/'t14-run.py'),'--measure',label,'--directory',str(OUT/'resources'),
        '--',str(OLD/'replay.ex'),mode,str(payload),str(output)],check=True)

def cache(leg,ray):
    return (OLD/'analysis'/('maximal-'+ray+'-profiles.npz') if leg=='max12'
            else OUT/(leg+'-'+ray+'-profiles.npz'))

def extract(leg):
    level,stride=configure(leg);directory=ROOT/'evolution'/leg
    assert (directory/'done.exit').read_text().strip()=='0'
    cells=a.selected();times=[];steps=[];group=[];previous=None;all_times=[];speed=0.
    payload=OUT/(leg+'-input.bin');output=OUT/(leg+'-q.bin')
    with payload.open('wb') as stream:
        def consume(group):
            nonlocal speed
            t=float(group[0]['meta'][0]);dt=a.DT/(2 if leg.endswith('half') else 1)
            step=round(t/dt);assert abs(t-step*dt)<2e-18
            all_times.append(t)
            # Verify the speed on every recorder cell, at every native snapshot.
            for fr in group:
                q=fr['values'];det=q[:,1]*q[:,3]-q[:,2]**2
                assert np.all(det>0) and np.all(q[:,4]>0)
                inv=(q[:,1]+q[:,3]+np.sqrt((q[:,1]-q[:,3])**2+4*q[:,2]**2))/(2*det)
                speed=max(speed,float(np.max(np.sqrt(np.maximum(inv,1/q[:,4]))+np.hypot(q[:,14],q[:,15]))))
            if not (step<=4 or step%stride==0 or step==int(.002/dt)):return
            lo,state=a.dense(group)
            if step==0:np.savez(OUT/(leg+'-initial.npz'),lo=lo,a=state)
            record=np.empty((len(cells),3+49*28));record[:,0]=a.H
            record[:,1]=(cells[:,1]+.5)*a.H;record[:,2]=(cells[:,0]+.5)*a.H-336
            for j in range(-3,4):
                for i in range(-3,4):
                    values=state[cells[:,1]+j-lo[1],cells[:,0]+i-lo[0]]
                    assert np.isfinite(values).all(),(leg,step,i,j)
                    record[:,3+(j+3)*7+i+3::49]=values
            record.tofile(stream);times.append(t);steps.append(step)
        for fr in BASE_FRAMES(directory/f't13-t7-stage-L{level}.xz'):
            if fr['phase']!=50:continue
            t=float(fr['meta'][0])
            if previous is not None and t!=previous:consume(group);group=[]
            group.append(fr);previous=t
        if group:consume(group)
    native('snapshots-'+leg,'--state',payload,output)
    q=np.memmap(output,mode='r',dtype='f8',shape=(len(times),len(cells),180))
    assert np.isfinite(q).all() and np.max(q[:,:,144])==0 and np.max(abs(q[:,:,151]))==0
    assert speed<1.1
    np.savez(OUT/(leg+'-meta.npz'),cells=cells,times=times,steps=steps,
        all_times=all_times,speed=speed)
    for ray,vec in a.NV.items():
        samples={n:np.array([a.sample(cells,frame,vec,n) for frame in q]) for n in (6,8)}
        np.savez_compressed(cache(leg,ray),times=times,steps=steps,r=a.R,p6=samples[6],p8=samples[8])
    (OUT/(leg+'-input.sha256')).write_text(hash_file(payload)+' '+str(payload)+'\n')
    payload.unlink() # Own regenerable stencil payload; exact hash and native output retained.
    print(leg,'all snapshots',len(all_times),'selected',len(times),'cells',len(cells),'speed',speed,flush=True)

def stages(leg):
    level,stride=configure(leg)
    check=module('t13closure','t13-check.py');check.TMP=ROOT
    replay=ROOT/'replay.ex'
    if not replay.exists():replay.symlink_to(OLD/'replay.ex')
    dt=a.DT/(2 if leg.endswith('half') else 1)
    def early(path):
        path=Path(path)
        if path.name=='t13-t7-stage-L12.xz':path=path.with_name(f't13-t7-stage-L{level}.xz')
        for frame in BASE_FRAMES(path):
            if frame['phase']==50 and frame['meta'][0]>=4*dt-1e-18:break
            yield frame
    check.mod.frames=early
    result=check.closure(ROOT/'evolution'/leg,level)
    save('t14-stage-closure.csv',[dict(run=leg,**result)],merge=True)
    closure_hashes(leg)
    report=module('t13report','t13-report.py');report.a=a;report.TMP=ROOT;report.OUT=OUT
    a.decoder.frames=early
    a.save=lambda name,rows:save(name.replace('t13-','t14-'),rows,merge=True)
    a.timed=lambda label,cmd:subprocess.run([PY,str(HERE/'t14-run.py'),'--measure',label,
        '--directory',str(OUT/'resources'),'--',*cmd],check=True)
    report.attribution();report.predictor();puncture(leg)
    print('Stage replay and attribution complete:',leg,flush=True)

def norm(values,kind):
    return float(np.max(abs(values))) if kind=='peak' else a.rms(values)

def closure_hashes(leg):
    # Recover each overwritten native output exactly from its lossless saved q.
    for path in sorted((OUT/(leg+'-stages')).glob('frame-*.npz')):
        with np.load(path) as z:
            q=z['q'];assert q.shape[1]==180
            digest=hashlib.sha256(q.tobytes()).hexdigest()
        path.with_name(path.stem+'-output.sha256').write_text(digest+'  native closure output reconstructed from lossless q\n')

def puncture(leg):
    configure(leg);dt=a.DT/(2 if leg.endswith('half') else 1)
    directory=OLD/'analysis/maximal-stages' if leg=='max12' else OUT/(leg+'-stages')
    rows=[];deltas={};names='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split()
    for path in sorted(directory.glob('frame-*.npz')):
        with np.load(path) as z:
            cells=z['cells'];stage=int(z['stage']);step=round(float(z['meta'][0])/dt)+1
            native_q=None
            for ox,oy in ((0,0),(1,0),(0,1),(1,1)):
                ix=np.flatnonzero((cells[:,0]==a.CENTRE+ox)&(cells[:,1]==oy))
                if not len(ix):continue
                assert len(ix)==1
                if native_q is None:native_q=z['q']
                v=native_q[ix[0]];vec=np.array([ox+.5,oy+.5]);vec/=np.linalg.norm(vec)
                weight=dt*(1/6,1/3,1/3,1/6)[stage]
                physical=float(v[39:41]@vec);adv=float(v[95:97]@vec)
                parts=dict(geometric=physical-adv,advection=adv,KO=float(v[67:69]@vec),total=float(v[163:165]@vec))
                for group,value in parts.items():
                    key=(step,ox,oy,group);deltas[key]=deltas.get(key,0.)+weight*value
                components=[('Gamma_n',parts['geometric'],adv,parts['KO'],parts['total'])]
                components += [(name,float(v[28+i]-v[84+i]),float(v[84+i]),float(v[56+i]),float(v[152+i])) for i,name in enumerate(names)]
                for name,geometry,advection,ko,total in components:
                    rows.append(dict(run=leg,step=step,stage=stage,stage_time_M=float(z['meta'][3]),
                        x_cell=ox,y_cell=oy,R_over_h=math.hypot(ox+.5,oy+.5),variable=name,
                        geometric_RHS=geometry,advection_RHS=advection,KO_RHS=ko,total_RHS=total,
                        direct_KO_over_total_abs=abs(ko/total) if total else float('nan')))
    save('t14-puncture-stages.csv',rows,merge=True)
    budgets=[]
    for step in range(1,5):
        for ox,oy in ((0,0),(1,0),(0,1),(1,1)):
            parts={group:deltas[step,ox,oy,group] for group in ('geometric','advection','KO','total')}
            budgets.append(dict(run=leg,step=step,x_cell=ox,y_cell=oy,**parts,
                direct_KO_over_total_abs=abs(parts['KO']/parts['total']) if parts['total'] else float('nan')))
    save('t14-puncture-budgets.csv',budgets,merge=True)

def arrays(leg,ray):
    z=np.load(cache(leg,ray));vec=a.NV[ray]
    result={n:{key:np.array([a.fields(row,vec)[key] for row in z[f'p{n}']])
        for key in a.fields(z['p6'][0],vec)} for n in (6,8)}
    return z['times'],result

def audit():
    rows=[]
    for leg in LEGS:
        if leg=='max12':
            old=next(x for x in csv.DictReader((HERE/'t13-run-audit.csv').open()) if x['run']=='maximal')
            directory=OLD/'evolution/maximal';resource=json.loads((directory/'maximal.resources.json').read_text())
            count=57;speed=float(old['current_shift_speed_proxy_max']);endpoint=float(old['end_time_M'])
        else:
            directory=ROOT/'evolution'/leg;resource=json.loads((directory/(leg+'.resources.json')).read_text())
            z=np.load(OUT/(leg+'-meta.npz'));count=len(z['all_times']);speed=float(z['speed']);endpoint=float(z['all_times'][-1])
            pp=(HERE/'params'/f't14-E-mid-{leg}.txt').read_text()
            for literal in ('sigma = 1\n','amr_transfer = point','t2_guard_initial_data_after_t0 = true',
                            'ems_use_maximal_initial_lapse = true'):assert literal in pp
        assert (directory/'done.exit').read_text().strip()=='0' and resource['returncode']==0 and not resource['gate_reason']
        floor=[x for path in directory.glob('t13-floors-*.csv') for x in csv.DictReader(path.open())]
        log=(directory/'run.log').read_text();start=log.index('GRAMRLevel::advance level 0 at time 0')
        row=dict(run=leg,native_snapshots=count,finest_steps=count-1,actual_endpoint_M=endpoint,
            chi_activations=sum(int(x['chi_activations']) for x in floor),
            lapse_activations=sum(int(x['lapse_activations']) for x in floor),
            nonfinite=sum(int(x['nonfinite']) for x in floor),
            chi_min=min(float(x['chi_min']) for x in floor),lapse_min=min(float(x['lapse_min']) for x in floor),
            reader_messages_after_first_advance=log[start:].count('Read EMSTRUMPET'),
            current_shift_speed_max=speed,peak_rss_bytes=resource['peak_rss_bytes'],wall_s=resource['wall_seconds'])
        assert row['chi_activations']==row['lapse_activations']==row['nonfinite']==row['reader_messages_after_first_advance']==0
        assert speed<1.1 and resource['peak_rss_bytes']<3e9
        assert count=={'max12':57,'max13':113,'max14':225,'max14-half':450}[leg]
        rows.append(row)
    save('t14-evolved-audit.csv',rows)
    return rows

def lobe(r,values):
    index=int(np.argmax(abs(values)));height=abs(values[index]);sign=np.sign(values[index])
    if height==0:return dict(width_cells=0.,width_M=0.,above_half_cells=0,clipped=True,peak_R_M=r[index])
    low=high=index
    while low>0 and sign*values[low-1]>=height/2:low-=1
    while high<len(r)-1 and sign*values[high+1]>=height/2:high+=1
    clipped=low==0 or high==len(r)-1
    def crossing(i,j):
        fraction=(height/2-sign*values[i])/(sign*(values[j]-values[i]))
        return i+fraction*(j-i),r[i]+fraction*(r[j]-r[i])
    lc,lr=(float(low),r[low]) if low==0 else crossing(low,low-1)
    hc,hr=(float(high),r[high]) if high==len(r)-1 else crossing(high,high+1)
    return dict(width_cells=hc-lc,width_M=hr-lr,above_half_cells=high-low+1,clipped=clipped,peak_R_M=r[index])

def screen():
    audit();history=[];profiles=[];constraints=[];constraint_history=[];native_rows=[];screen_rows=[];temporal=[];summaries=[]
    for ray,vec in a.NV.items():
        data={leg:arrays(leg,ray) for leg in LEGS}
        clocks=data['max12'][0];indices={}
        for leg,(times,fields) in data.items():
            indices[leg]=np.array([np.argmin(abs(times-t)) for t in clocks])
            defect=float(np.max(abs(times[indices[leg]]-clocks)))
            assert defect<2e-18,(leg,ray,defect)
            for ti,t in enumerate(times):
                gamma=fields[6]['Gamma'][ti]-fields[6]['Gamma'][0]
                spread=gamma-(fields[8]['Gamma'][ti]-fields[8]['Gamma'][0])
                history.append(dict(run=leg,ray=ray,time_M=t,Gamma_peak=norm(gamma,'peak'),Gamma_RMS=norm(gamma,'RMS'),
                    interpolation_peak=norm(spread,'peak'),interpolation_RMS=norm(spread,'RMS'),
                    probe_Gamma=float(np.interp(.0015,a.R,gamma))))
            for ci in (0,1,14,28,42,56):
                ti=indices[leg][ci]
                for j,r in enumerate(a.R):
                    profiles.append(dict(run=leg,ray=ray,time_M=clocks[ci],R_M=r,
                        **{k:float(fields[6][k][ti,j]) for k in FIELDS},
                        Gamma_disturbance=float(fields[6]['Gamma'][ti,j]-fields[6]['Gamma'][0,j]),
                        Gamma_spread=float(fields[6]['Gamma'][ti,j]-fields[8]['Gamma'][ti,j])))
            for field in ('C_Gamma','Ham','Mom','GaussE'):
                for kind in ('peak','RMS'):
                    values=[norm(v,kind) for v in fields[6][field]]
                    errors=[norm(x-y,kind) for x,y in zip(fields[6][field],fields[8][field])]
                    ti=indices[leg][-1]
                    constraints.append(dict(run=leg,ray=ray,field=field,norm=kind,initial=values[0],
                        endpoint=values[ti],history_sup=max(values),endpoint_interpolation=errors[ti],
                        history_sup_interpolation=max(errors)))
                    constraint_history += [dict(run=leg,ray=ray,field=field,norm=kind,time_M=t,
                        value=value,interpolation=error) for t,value,error in zip(times,values,errors)]
            configure(leg)
            directory=OLD/'analysis' if leg=='max12' else OUT
            tag='maximal' if leg=='max12' else leg
            meta=np.load(directory/(tag+'-meta.npz'));cells=meta['cells'];times=meta['times']
            q=np.memmap(directory/(tag+'-q.bin'),mode='r',dtype='f8',shape=(len(times),len(cells),180))
            x=(cells[:,0]+.5)*a.H-336;y=(cells[:,1]+.5)*a.H;r=np.hypot(x,y)
            mask=(cells[:,1]==0)&(x>0) if ray=='axis' else cells[:,0]-a.CENTRE==cells[:,1]
            mask &= (r>=a.R[0])&(r<=a.R[-1]);ix=np.flatnonzero(mask);ix=ix[np.argsort(r[ix])]
            assert len(ix)>1
            ti=int(np.argmin(abs(times-clocks[-1])))
            values=(q[ti,ix,11:13]-q[0,ix,11:13])@vec
            native_rows.append(dict(run=leg,ray=ray,native_cells=len(ix),Gamma_peak=float(max(abs(values))),
                Gamma_cell_RMS=float(np.sqrt(np.mean(values**2))),**lobe(r[ix],values),
                ray_cell_spacing_M=a.H*(1 if ray=='axis' else math.sqrt(2)),
                condition='first cartoon row y=h/2 axis proxy' if ray=='axis' else 'native diagonal centres'))
        for field in FIELDS:
            v={leg:data[leg][1][6][field][indices[leg]] for leg in LEGS}
            e={leg:v[leg]-data[leg][1][8][field][indices[leg]] for leg in LEGS}
            baseline_half=np.load(OLD/'analysis'/('maximal-half-'+ray+'-profiles.npz'))
            baseline_half_field=np.array([a.fields(row,vec)[field] for row in baseline_half['p6']])[::2]
            assert baseline_half_field.shape==v['max12'].shape
            for kind in ('peak','RMS'):
                collected=[]
                # Empirical zero-field floor from the captured numerical t=0
                # Gamma/metric-Gamma, not a static reader or a fitted background.
                floors={leg:norm(v[leg][0],kind) if field in ('Gamma','metric_Gamma') else 0. for leg in LEGS}
                floor0=floors['max12']+floors['max13'];floor1=floors['max13']+floors['max14']
                for ci,t in enumerate(clocks):
                    d0=norm(v['max12'][ci]-v['max13'][ci],kind);d1=norm(v['max13'][ci]-v['max14'][ci],kind)
                    b0=norm(e['max12'][ci],kind)+norm(e['max13'][ci],kind)
                    b1=norm(e['max13'][ci],kind)+norm(e['max14'][ci],kind)
                    te=norm(v['max14'][ci]-v['max14-half'][ci],kind)
                    old_te=norm(v['max12'][ci]-baseline_half_field[ci],kind)
                    floor_dominated=d0<=5*floor0 or d1<=5*floor1
                    significant=d0>5*b0 and d1>5*b1 and d0>0 and d1>0 and not floor_dominated
                    temporal_ok=te<=.2*min(d0,d1)
                    rho=d1/d0 if d0 and not floor_dominated else float('nan')
                    status=('numerical-floor' if floor_dominated else 'uncertainty-dominated' if not significant else 'temporal-unqualified' if not temporal_ok
                        else 'contracting' if rho<=.8 else 'inconclusive' if rho<1 else 'noncontracting')
                    row=dict(ray=ray,field=field,norm=kind,time_M=t,D12_13=d0,D13_14=d1,rho_D=rho,
                        rho_lower=max(0.,d1-b1)/(d0+b0) if d0+b0 and not floor_dominated else float('nan'),
                        rho_upper=(d1+b1)/(d0-b0) if d0>b0 and not floor_dominated else float('nan'),
                        numerical_floor12_13=floor0,numerical_floor13_14=floor1,
                        interpolation12_13=b0,interpolation13_14=b1,
                        margin12_13=d0/(5*b0) if b0 else float('inf'),margin13_14=d1/(5*b1) if b1 else float('inf'),
                        temporal_difference=te,temporal_over_D12_13=te/d0 if d0 else float('nan'),
                        temporal_over_D13_14=te/d1 if d1 else float('nan'),max12_temporal_difference=old_te,
                        max12_temporal_over_D12_13=old_te/d0 if d0 else float('nan'),status=status)
                    collected.append(row);screen_rows.append(row)
                    temporal.append({k:row[k] for k in ('ray','field','norm','time_M','temporal_difference',
                        'temporal_over_D12_13','temporal_over_D13_14','max12_temporal_difference','max12_temporal_over_D12_13')})
                nonzero=collected[1:];qualified=[r for r in nonzero if r['status'] not in ('numerical-floor','uncertainty-dominated','temporal-unqualified')]
                end=collected[-1]
                summaries.append(dict(ray=ray,field=field,norm=kind,endpoint_rho=end['rho_D'],endpoint_status=end['status'],
                    endpoint_margin12_13=end['margin12_13'],endpoint_margin13_14=end['margin13_14'],
                    endpoint_temporal_over_D12_13=end['temporal_over_D12_13'],endpoint_temporal_over_D13_14=end['temporal_over_D13_14'],
                    history_D12_13_sup=max(r['D12_13'] for r in nonzero),history_D13_14_sup=max(r['D13_14'] for r in nonzero),
                    history_sup_ratio=max(r['D13_14'] for r in nonzero)/max(r['D12_13'] for r in nonzero),
                    history_max_rho=max(r['rho_D'] for r in nonzero),
                    qualified_history_max_rho=max((r['rho_D'] for r in qualified),default=float('nan')),
                    qualified_times=len(qualified),contracting_times=sum(r['status']=='contracting' for r in nonzero),
                    noncontracting_times=sum(r['status']=='noncontracting' for r in nonzero),
                    inconclusive_times=sum(r['status']=='inconclusive' for r in nonzero),
                    uncertainty_dominated_times=sum(r['status']=='uncertainty-dominated' for r in nonzero),
                    numerical_floor_times=sum(r['status']=='numerical-floor' for r in nonzero),
                    temporal_unqualified_times=sum(r['status']=='temporal-unqualified' for r in nonzero),
                    history_min_margin12_13=min(r['margin12_13'] for r in nonzero),
                    history_min_margin13_14=min(r['margin13_14'] for r in nonzero),
                    history_max_temporal_fraction=max(max(r['temporal_over_D12_13'],r['temporal_over_D13_14']) for r in nonzero),
                    maximum_clock_defect_M=max(float(np.max(abs(data[leg][0][indices[leg]]-clocks))) for leg in LEGS)))
    for name,rows in [('t14-history.csv',history),('t14-profiles.csv',profiles),('t14-constraints.csv',constraints),
                      ('t14-constraint-history.csv',constraint_history),
                      ('t14-native-lobes.csv',native_rows),('t14-screen-history.csv',screen_rows),
                      ('t14-temporal.csv',temporal),('t14-screen.csv',summaries)]:save(name,rows)
    print(json.dumps(summaries,indent=2),flush=True)

if __name__=='__main__':
    command=sys.argv[1]
    if command=='extract':extract(sys.argv[2])
    elif command=='stages':stages(sys.argv[2])
    elif command=='screen':screen()
    else:raise ValueError(command)
