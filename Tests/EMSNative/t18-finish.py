#!/usr/bin/env python3
"""Close measured timestep evidence and update the low-rung performance budget."""
import csv,hashlib,importlib.util,json,math,sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t18')
def read(name):return list(csv.DictReader((HERE/name).open()))
def save(name,rows):
    assert rows
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def main():
    windows=read('t18-window-audit.csv');controls=read('t18-launch-full-grid-control.csv')
    assert len(windows)==3 and all(r['all_13_level_screens_pass']=='True' for r in windows)
    assert len(controls)==4 and all(int(r['differing_bits'])==0 for r in controls)
    errors=read('t18-temporal-spatial.csv');history=read('t18-launch-constraints.csv')
    for dt in (.375,.5):
        rows=[r for r in errors if float(r['dt_multiplier'])==dt]
        assert len(rows)==528
        assert all(r['verdict']=='resolved subdominant' for r in rows if r['field'] not in ('GaussB','Lambda'))
        assert all(float(r['temporal_difference'])==float(r['spatial_difference'])==0 for r in rows if r['field'] in ('GaussB','Lambda'))
    decision=dict(status='READY-EXCEPT',dt_multiplier=.5,calibration_rung='coarse',
        finest_level=12,finest_spacing_M=1.75/2**12,largest_tested=True,
        reason='All 39 registered window screens pass; binary temporal errors are resolved below the T17 rung difference for every nonzero channel at all common launch clocks, rays, holes and norms.',
        full_grid_native_values=sum(int(r['values']) for r in controls),full_grid_differing_bits=0,
        unresolved_zero_channels=['GaussB','Lambda'],
        limitations=['0.875 Mi isolated finite window and 0.001922607421875 Mi binary launch; no long-merger theorem',
            'CFL margins are frozen scalar-wave proxies',
            'Floor observations cover registered RK inputs/endpoints; final pre-clamp activation counts were not instrumented',
            'Cluster matrix remains a prediction; T17 horizon/spatial admission remains unavailable'],
        unchanged=['gauge','KO','equations','point transfers','puncture','Float64','RHFinder'])
    (HERE/'t18-dt-recommendation.json').write_text(json.dumps(decision,indent=2)+'\n')
    index={(r['run'],r['hole'],r['ray'],r['field'],r['measure'],r['time_M'],r['norm']):r for r in history}
    comparisons=[]
    for r in history:
        if r['run'] not in ('coarse-dt0.375','coarse-dt0.5') or float(r['time_M'])==0:continue
        base=index[('coarse-dt0.25',r['hole'],r['ray'],r['field'],r['measure'],r['time_M'],r['norm'])]
        value,ref=float(r['amplitude']),float(base['amplitude'])
        comparisons.append({**r,'dt025_amplitude':ref,'norm_difference_from_dt025':value-ref,
            'fractional_norm_difference_from_dt025':value/ref-1 if ref else 'undefined'})
    save('t18-launch-comparison.csv',comparisons)
    # The original finite windows provide Euclidean volume maxima; the binary
    # captures separately provide radial peak/RMS profiles and native constraints.
    common={};clocks=[.000640869140625,.00128173828125,.001922607421875]
    for dt in (.25,.375,.5):
        data=list(csv.DictReader((ROOT/'local'/('dt-'+str(dt)+'-window')/'t18-launch-history.csv').open()))
        for t in clocks:
            common[dt,t]=next(r for r in data if abs(float(r['time_M'])-t)<1e-14)
    volume=[]
    for dt in (.25,.375,.5):
        for t in clocks:
            row=dict(dt_multiplier=dt,time_M=t)
            for k in ('delta_Gamma_W_max','delta_lapse_W_max','delta_shift_W_max'):
                value=float(common[dt,t][k]);base=float(common[.25,t][k])
                row[k]=value;row[k+'_difference_from_dt025']=value-base
                row[k+'_fractional_change_from_dt025']=value/base-1 if base else 'undefined'
            volume.append(row)
    save('t18-launch-common.csv',volume)
    # No new model algebra: regenerate the existing witnessed fit with selected
    # numerical metadata. A timestep changes the physical clock per coarse step.
    spec=importlib.util.spec_from_file_location('model',HERE/'t18-model.py')
    model=importlib.util.module_from_spec(spec);spec.loader.exec_module(model)
    model.model()
    matrix=read('t18-cluster-matrix.csv')
    best={rung:min((r for r in matrix if r['rung']==rung),key=lambda r:float(r['predicted_s_per_step'])) for rung in ('coarse','fine')}
    wall=[]
    for rung in ('coarse','fine'):
        for dt in sorted({.25,decision['dt_multiplier']}):
            rate=float(best[rung]['predicted_s_per_step'])
            for d in (32,12,16):
                infall=306*(d/32)**1.5;duration=infall+100*2.002174604
                steps=math.ceil(duration/(1.75*dt));evolution=steps*rate/3600
                wall.append(dict(rung=rung,dt_multiplier=dt,separation_Mi=d,coarse_steps=steps,
                    s_per_step=rate,evolution_hours=evolution,routine_finder_hours=.01*evolution,
                    cold_qualification_pairs=4,cold_qualification_low_hours=4*5522/3600,
                    cold_qualification_high_hours=4*6164/3600,
                    combined_low_hours=1.01*evolution+4*5522/3600,
                    combined_high_hours=1.01*evolution+4*6164/3600,
                    evolution_half_to_double_combined_low_hours=.505*evolution+4*5522/3600,
                    evolution_half_to_double_combined_high_hours=2.02*evolution+4*6164/3600,
                    recommended_dt_scope=rung=='coarse' and dt==.5,
                    limitation='predicted layout; extra cold searches, I/O, restart and extra qualification excluded; d12/16 require new initial data'))
    save('t18-wallclock.csv',wall)
    # The worker receives exactly the low-resolution production hierarchy.
    spec=importlib.util.spec_from_file_location('work',HERE/'t18-launch-work.py')
    work=importlib.util.module_from_spec(spec);spec.loader.exec_module(work)
    base=(work.w.T17/'t17-merger-draft.txt').read_text()
    radii=next(line.split('=',1)[1].split() for line in base.splitlines() if line.startswith('mass_extraction_radii ='))[:12]
    assert len(radii)==12
    for size in (16,24,32,48,128):
        text=work.replace(base,dict(max_level=12,max_box_size=size,dt_multiplier=.5,
            regrid_interval='4 4 4 4 4 4 1 0 0 0 0 0 0',
            max_steps=5,stop_time=4.375,checkpoint_interval=5,plot_interval=-1,
            RH_activate='false',activate_extraction=0,activate_rs_extraction=0,activate_em_extraction=0,
            activate_mq_extraction=0,activate_gw_extraction=0,activate_mass_extraction=0,
            num_mass_extraction_radii=12,mass_extraction_levels=' '.join(map(str,range(1,13))),
            mass_extraction_radii=' '.join(radii),
            t13_launch_stop_time=0,t18_monitor='false'))
        (HERE/('t18-cluster-b'+str(size)+'.txt')).write_text(text)
    result=json.loads((HERE/'t18-window-results.json').read_text())
    result.update(status='VERIFIED_FINITE_WINDOWS',per_level_screens_verified=39,
        verification='t18-window-audit.csv',launch_comparison='t18-launch-comparison.csv',
        temporal_spatial_comparison='t18-temporal-spatial.csv',recommendation=decision)
    (HERE/'t18-window-results.json').write_text(json.dumps(result,indent=2)+'\n')
    receipts=[]
    for p in sorted(ROOT.rglob('*.resources.json')):
        r=json.loads(p.read_text())
        receipts.append(dict(process=r['process'],receipt=str(p),wall_seconds=r['wall_seconds'],
            peak_RSS_bytes=r['peak_rss_bytes'],returncode=r['returncode'],gate_reason=r['gate_reason']))
    save('t18-resources.csv',receipts)
    inputs=json.loads((HERE/'t18-inputs.json').read_text());seen={r['path'] for r in inputs}
    for p in [ROOT/'launch.ex',ROOT/'launch-replay.ex',HERE/'t18-launch.cpp',HERE/'t16-replay.cpp',HERE/'t18-launch-registration.json']:
        record=dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
        if str(p) not in seen:inputs.append(record)
        else:inputs=[record if r['path']==str(p) else r for r in inputs]
    (HERE/'t18-inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    print('Recommended coarse dt=.5; measured maximum process RSS',max(int(r['peak_RSS_bytes']) for r in receipts))
    for r in wall:
        if r['recommended_dt_scope']:print(r)
if __name__=='__main__':main()
