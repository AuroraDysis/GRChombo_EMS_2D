#!/usr/bin/env python3
"""One predeclared I8/I10 rescreen of stored current states; no evolution."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, json, math, os, subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t14d');ROOT.mkdir(exist_ok=True)
PY='/Users/auroradysis/miniconda3/bin/python'
spec=importlib.util.spec_from_file_location('t13',HERE/'t13-analyze.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
OLD=Path('/private/tmp/ems-t13');T14=Path('/private/tmp/ems-t14')
LEGS={'max12':(12,1,OLD/'evolution/maximal'),
      'max13':(13,2,T14/'evolution/max13'),
      'max14':(14,4,T14/'evolution/max14'),
      'max14-half':(14,8,T14/'evolution/max14-half')}
FIELDS=('Gamma','metric_Gamma','shift','lapse')
CONSTRAINTS=('C_Gamma','Ham','Mom','GaussE')

def read(name):return list(csv.DictReader((HERE/name).open()))
def save(name,rows):
    assert rows,name
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def clocks():
    return np.array([float(x['time_M']) for x in read('t14-screen-history.csv')
        if x['field']=='Gamma' and x['ray']=='axis' and x['norm']=='peak'])
def configure(leg):
    level,stride,directory=LEGS[leg]
    a.H=(7/12)/2**level;a.DT=a.H/4;a.CENTRE=round(336/a.H)
    return level,stride,directory
def cache(leg,ray):return ROOT/f't14d-{leg}-{ray}-profiles.npz'
def norm(v,kind):return float(np.max(abs(v))) if kind=='peak' else a.rms(v)

def register():
    inputs=[OLD/'replay.ex',HERE/'t13-analyze.py',HERE/'t7-check.py',HERE/'t14-screen-history.csv']
    inputs += [d/f't13-t7-stage-L{level}.xz' for level,stride,d in LEGS.values()]
    missing=[str(p) for p in inputs if not p.is_file()]
    if missing:raise FileNotFoundError('Required stored inputs missing: '+', '.join(missing))
    assert len(clocks())==57
    for level,stride,d in LEGS.values():assert (d/'done.exit').read_text().strip()=='0'
    record=dict(baseline='ba95217',central='I8',spread='abs(I8-I10)',
        unchanged=['runs','57 clocks','W','peak/radial RMS','5x significance','temporal <=0.2 of both spatial differences',
                   'empirical initial Gamma floor','status classes','registered decision table'],
        W_M=[float(a.R[0]),float(a.R[-1])],samples=len(a.R),times_M=clocks().tolist(),
        history_start_rule='first positive common clock when all four fields, both rays and both norms exceed 5x interpolation',
        Stage_D_exclusion_R_M=.005,operator_variants_remaining=0,
        inputs=[dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for p in inputs])
    (HERE/'t14d-registration.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Registered one I8/I10 variant; all required stored inputs present.',flush=True)

def selected():
    cells=set()
    for vec in a.NV.values():
        xy=a.R[:,None]*vec
        start=np.floor(np.c_[(336+xy[:,0])/a.H-.5,xy[:,1]/a.H-.5]).astype(int)-4
        for sx,sy in start:
            for j in range(10):
                for i in range(10):cells.add((sx+i,sy+j if sy+j>=0 else -sy-j-1))
    return np.array(sorted(cells,key=lambda p:(p[1],p[0])),int)

def extract(leg):
    level,stride,directory=configure(leg);cells=selected();wanted=clocks()
    payload=ROOT/f't14d-{leg}-input.bin';output=ROOT/f't14d-{leg}-q.bin'
    times=[];group=[];previous=None;coverage=0;all_snapshots=0
    with payload.open('wb') as stream:
        def consume(group):
            nonlocal coverage,all_snapshots
            t=float(group[0]['meta'][0]);all_snapshots+=1
            ci=int(np.argmin(abs(wanted-t)))
            if abs(wanted[ci]-t)>2e-18:return
            assert ci==len(times),(leg,ci,len(times))
            lo,state=a.dense(group)
            record=np.empty((len(cells),3+49*28));record[:,0]=a.H
            record[:,1]=(cells[:,1]+.5)*a.H;record[:,2]=(cells[:,0]+.5)*a.H-336
            for j in range(-3,4):
                for i in range(-3,4):
                    ix=cells[:,0]+i-lo[0];iy=cells[:,1]+j-lo[1]
                    assert np.all((ix>=0)&(ix<state.shape[1])&(iy>=0)&(iy<state.shape[0])),('Stored recorder support missing',leg,t,i,j)
                    values=state[iy,ix]
                    assert np.isfinite(values).all(),('Stored recorder cells missing',leg,t,i,j)
                    record[:,3+(j+3)*7+i+3::49]=values;coverage+=len(cells)
            record.tofile(stream);times.append(t)
        for fr in a.decoder.frames(directory/f't13-t7-stage-L{level}.xz'):
            if fr['phase']!=50:continue
            t=float(fr['meta'][0])
            if previous is not None and t!=previous:consume(group);group=[]
            group.append(fr);previous=t
        if group:consume(group)
    assert len(times)==57 and max(abs(np.array(times)-wanted))<2e-18
    subprocess.run([PY,str(HERE/'t14-run.py'),'--measure','t14d-native-'+leg,
        '--directory',str(ROOT/'resources'),'--',str(OLD/'replay.ex'),'--state',str(payload),str(output)],check=True)
    q=np.memmap(output,mode='r',dtype='f8',shape=(57,len(cells),180))
    assert np.isfinite(q).all() and np.max(abs(q[:,:,144]))==0 and np.max(abs(q[:,:,151]))==0
    identity=[]
    for ray,vec in a.NV.items():
        samples={n:np.array([a.sample(cells,frame,vec,n) for frame in q]) for n in (8,10)}
        np.savez_compressed(cache(leg,ray),times=times,r=a.R,p8=samples[8],p10=samples[10])
        original=(OLD/'analysis'/f'maximal-{ray}-profiles.npz' if leg=='max12'
                  else T14/'analysis'/f'{leg}-{ray}-profiles.npz')
        with np.load(original) as z:
            indices=np.array([np.argmin(abs(z['times']-t)) for t in wanted])
            old=z['p8'][indices];mismatches=int(np.count_nonzero(old.view('u8')!=samples[8].view('u8')))
            assert mismatches==0,(leg,ray,mismatches)
        identity.append(dict(run=leg,ray=ray,I8_bit_mismatches=mismatches,Float64_values=old.size))
    save(f't14d-identity-{leg}.csv',identity)
    meta=dict(run=leg,native_cells=len(cells),common_snapshots=len(times),all_native_snapshots=all_snapshots,
        native_stencil_cells_checked=coverage,maximum_clock_defect_M=float(max(abs(np.array(times)-wanted))),
        input_bytes=payload.stat().st_size,input_sha256=digest(payload),q_bytes=output.stat().st_size,q_sha256=digest(output),
        support_buffer_cells=8,static_reader_calls=0)
    (ROOT/f't14d-{leg}-extraction.json').write_text(json.dumps(meta,indent=2)+'\n')
    payload.unlink() # Own regenerable payload; exact hash and native q retained.
    print(json.dumps(meta),flush=True)

def arrays(leg,ray):
    with np.load(cache(leg,ray)) as z:
        return {n:{key:np.array([a.fields(row,a.NV[ray])[key] for row in z[f'p{n}']])
            for key in (*FIELDS,*CONSTRAINTS)} for n in (8,10)}

def bounds(d0,d1,b0,b1,factor,floor=False):
    if floor:return float('nan'),float('nan')
    low=max(0.,d1-factor*b1)/(d0+factor*b0) if d0+factor*b0 else float('nan')
    high=(d1+factor*b1)/(d0-factor*b0) if d0>factor*b0 else float('inf')
    return low,high

def screen():
    registration=json.loads((HERE/'t14d-registration.json').read_text())
    assert registration['central']=='I8' and registration['spread']=='abs(I8-I10)'
    times=clocks();history=[];summaries=[];constraint_rows=[];orders=[];order_summaries=[];amplitudes=[]
    for ray,vec in a.NV.items():
        data={leg:arrays(leg,ray) for leg in LEGS}
        for field in FIELDS:
            v={leg:data[leg][8][field] for leg in LEGS}
            e={leg:v[leg]-data[leg][10][field] for leg in LEGS}
            for kind in ('peak','RMS'):
                floor={leg:norm(v[leg][0],kind) if field in ('Gamma','metric_Gamma') else 0. for leg in LEGS}
                f0=floor['max12']+floor['max13'];f1=floor['max13']+floor['max14'];rows=[]
                for step,t in enumerate(times):
                    d0=norm(v['max12'][step]-v['max13'][step],kind);d1=norm(v['max13'][step]-v['max14'][step],kind)
                    b0=norm(e['max12'][step],kind)+norm(e['max13'][step],kind)
                    b1=norm(e['max13'][step],kind)+norm(e['max14'][step],kind)
                    te=norm(v['max14'][step]-v['max14-half'][step],kind)
                    floored=d0<=5*f0 or d1<=5*f1
                    significant=d0>5*b0 and d1>5*b1 and d0>0 and d1>0 and not floored
                    temporal=te<=.2*min(d0,d1)
                    rho=d1/d0 if d0 and not floored else float('nan')
                    status=('numerical-floor' if floored else 'uncertainty-dominated' if not significant
                            else 'temporal-unqualified' if not temporal else 'contracting' if rho<=.8
                            else 'inconclusive' if rho<1 else 'noncontracting')
                    r1=bounds(d0,d1,b0,b1,1,floored);r5=bounds(d0,d1,b0,b1,5,floored)
                    failed=[pair for pair,d,b in (('12,13',d0,b0),('13,14',d1,b1)) if d<=5*b or d==0]
                    ff=[pair for pair,d,f in (('12,13',d0,f0),('13,14',d1,f1)) if d<=5*f]
                    row=dict(ray=ray,field=field,norm=kind,time_M=t,step=step,max13_step=2*step,max14_step=4*step,control_step=8*step,
                        D12_13=d0,D13_14=d1,rho_D=rho,rho_lower_1x=r1[0],rho_upper_1x=r1[1],rho_lower_5x=r5[0],rho_upper_5x=r5[1],
                        interpolation12_13=b0,interpolation13_14=b1,numerical_floor12_13=f0,numerical_floor13_14=f1,
                        margin12_13=d0/(5*b0) if b0 else float('inf'),margin13_14=d1/(5*b1) if b1 else float('inf'),
                        temporal_difference=te,temporal_over_D12_13=te/d0 if d0 else float('nan'),
                        temporal_over_D13_14=te/d1 if d1 else float('nan'),significance_failed_pairs=';'.join(failed),
                        floor_failed_pairs=';'.join(ff),status=status)
                    rows.append(row);history.append(row)
                positive=rows[1:];qualified=[x for x in positive if x['status'] in ('contracting','inconclusive','noncontracting')];end=rows[-1]
                summaries.append(dict(ray=ray,field=field,norm=kind,
                    endpoint_rho=end['rho_D'],endpoint_status=end['status'],
                    endpoint_margin12_13=end['margin12_13'],endpoint_margin13_14=end['margin13_14'],
                    endpoint_rho_lower_5x=end['rho_lower_5x'],endpoint_rho_upper_5x=end['rho_upper_5x'],
                    qualified_times=len(qualified),contracting_times=sum(x['status']=='contracting' for x in positive),
                    uncertainty_dominated_times=sum(x['status']=='uncertainty-dominated' for x in positive),
                    numerical_floor_times=sum(x['status']=='numerical-floor' for x in positive),
                    temporal_unqualified_times=sum(x['status']=='temporal-unqualified' for x in positive),
                    inconclusive_times=sum(x['status']=='inconclusive' for x in positive),noncontracting_times=sum(x['status']=='noncontracting' for x in positive),
                    initial_status=rows[0]['status'],qualified_history_max_rho=max((x['rho_D'] for x in qualified),default=float('nan')),
                    history_sup_ratio=max(x['D13_14'] for x in positive)/max(x['D12_13'] for x in positive),
                    history_max_temporal_fraction=max(max(x['temporal_over_D12_13'],x['temporal_over_D13_14']) for x in positive)))
        for field in CONSTRAINTS:
            for kind in ('peak','RMS'):
                norms={leg:[norm(x,kind) for x in data[leg][8][field]] for leg in LEGS}
                spreads={leg:[norm(x-y,kind) for x,y in zip(data[leg][8][field],data[leg][10][field])] for leg in LEGS}
                for leg in LEGS:
                    constraint_rows += [dict(run=leg,ray=ray,field=field,norm=kind,time_M=t,step=step,value=norms[leg][step],interpolation=spreads[leg][step])
                        for step,t in enumerate(times)]
                for coarse,fine in (('max12','max13'),('max13','max14')):
                    pair=coarse[3:]+','+fine[3:];rows=[]
                    for step,t in enumerate(times):
                        c=norms[coarse][step];f=norms[fine][step];ec=spreads[coarse][step];ef=spreads[fine][step]
                        p=math.log2(c/f) if c>0 and f>0 else float('nan')
                        qualified=c>5*ec and f>5*ef and c>0 and f>0
                        lo=math.log2((c-5*ec)/(f+5*ef)) if qualified else float('nan')
                        hi=math.log2((c+5*ec)/(f-5*ef)) if qualified else float('nan')
                        row=dict(ray=ray,field=field,norm=kind,pair=pair,time_M=t,step=step,coarse_norm=c,fine_norm=f,
                            coarse_interpolation=ec,fine_interpolation=ef,p=p,p_lower_5x=lo,p_upper_5x=hi,
                            qualified=qualified,design_p=4.,p_minus_design=p-4)
                        rows.append(row);orders.append(row)
                    cs=max(norms[coarse]);fs=max(norms[fine]);positive=rows[1:];qp=[x for x in positive if x['qualified']]
                    order_summaries.append(dict(ray=ray,field=field,norm=kind,pair=pair,initial_p=rows[0]['p'],
                        endpoint_p=rows[-1]['p'],endpoint_qualified=rows[-1]['qualified'],endpoint_p_lower_5x=rows[-1]['p_lower_5x'],
                        endpoint_p_upper_5x=rows[-1]['p_upper_5x'],history_sup_p=math.log2(cs/fs) if cs>0 and fs>0 else float('nan'),
                        history_min_p=min(x['p'] for x in positive),history_max_p=max(x['p'] for x in positive),
                        history_qualified_times=len(qp),qualified_history_min_p=min((x['p'] for x in qp),default=float('nan')),
                        qualified_history_max_p=max((x['p'] for x in qp),default=float('nan')),design_p=4.))
        for leg in LEGS:
            v=data[leg][8]['Gamma'];w=data[leg][10]['Gamma']
            for step,t in enumerate(times):
                g=v[step]-v[0];spread=g-(w[step]-w[0])
                amplitudes.append(dict(run=leg,ray=ray,time_M=t,step=step,peak=norm(g,'peak'),RMS=norm(g,'RMS'),
                    interpolation_peak=norm(spread,'peak'),interpolation_RMS=norm(spread,'RMS')))
    unqualified=[x for x in history if x['status'] in ('uncertainty-dominated','numerical-floor','temporal-unqualified')]
    positive_bad=[x for x in unqualified if x['step']>0]
    noncontracting=[x for x in history if x['step']>0 and x['status']=='noncontracting']
    marginal=[x for x in history if x['step']>0 and x['status']=='inconclusive']
    verdict='insufficient' if noncontracting else 'inconclusive' if positive_bad or marginal else 'pass'
    starts=[step for step in range(1,57) if all(x['margin12_13']>1 and x['margin13_14']>1 for x in history if x['step']==step)]
    result=dict(verdict=verdict,endpoint_pass=all(x['endpoint_status']=='contracting' for x in summaries),
        positive_unqualified=len(positive_bad),initial_unqualified_rows=sum(x['step']==0 for x in unqualified),
        initial_floor_rows=sum(x['step']==0 and x['status']=='numerical-floor' for x in unqualified),
        initial_uncertainty_rows=sum(x['step']==0 and x['status']=='uncertainty-dominated' for x in unqualified),
        noncontracting=len(noncontracting),marginal=len(marginal),
        first_all_significant_step=starts[0] if starts else None,first_all_significant_time_M=float(times[starts[0]]) if starts else None,
        later_sampling_exceptions=[step for step in range(starts[0],57) if step not in starts] if starts else [],
        max_temporal_fraction=max(x['history_max_temporal_fraction'] for x in summaries),
        max_qualified_rho=max(x['qualified_history_max_rho'] for x in summaries),operator_variants_remaining=0)
    for name,rows in [('screen-history',history),('screen',summaries),('unqualified',unqualified),('constraint-history',constraint_rows),
                      ('constraint-orders-history',orders),('constraint-orders',order_summaries),('amplitude-history',amplitudes)]:
        save('t14d-'+name+'.csv',rows)
    (HERE/'t14d-verdict.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    if sys.argv[1]=='register':register()
    elif sys.argv[1]=='extract':extract(sys.argv[2])
    elif sys.argv[1]=='screen':screen()
    else:raise ValueError(sys.argv[1])
