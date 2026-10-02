#!/usr/bin/env python3
"""All-field, coordinate-matched bit checks and native timestep summaries."""
import csv
import json
from pathlib import Path
import sys
import h5py
import numpy as np
import hashlib

HERE = Path(__file__).resolve().parent
ROOT = Path('/private/tmp/ems-t18/local')

def save(name, rows):
    if not rows:
        return
    with (HERE/name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

def read_level(group, components):
    boxes = [[int(b[k]) for k in ('lo_i','lo_j','hi_i','hi_j')] for b in group['boxes'][:]]
    offsets = group['data:offsets=0'][:]
    coords, values = [], []
    for k, (x0, y0, x1, y1) in enumerate(boxes):
        nx, ny = x1-x0+1, y1-y0+1
        chunk = group['data:datatype=0'][int(offsets[k]):int(offsets[k+1])]
        # HDF5 checkpoints/plots retain three ghost cells in these controls.
        ghost = 3 if chunk.size == components*(nx+6)*(ny+6) else 0
        assert chunk.size == components*(nx+2*ghost)*(ny+2*ghost)
        a = chunk.reshape(components, ny+2*ghost, nx+2*ghost)
        a = a[:, ghost:ghost+ny, ghost:ghost+nx]
        x, y = np.meshgrid(np.arange(x0,x1+1), np.arange(y0,y1+1))
        coords.append(np.column_stack((x.ravel(), y.ravel())))
        values.append(a.reshape(components,-1).T)
    coords, values = np.concatenate(coords), np.concatenate(values)
    order = np.lexsort((coords[:,0], coords[:,1]))
    coords, values = coords[order], values[order]
    assert len(np.unique(coords, axis=0)) == len(coords)
    return coords, np.ascontiguousarray(values)

def neutrality():
    reference = ROOT/'neutral-b128-t1'
    rows = []
    for candidate in sorted(ROOT.glob('neutral-*')):
        if not (candidate/'done.exit').exists() or (candidate/'done.exit').read_text().strip() != '0':
            continue
        for kind, name in [('t0','chk/EMS_000000.2d.hdf5'), ('step3','chk/EMS_000003.2d.hdf5')]:
            with h5py.File(reference/name) as base, h5py.File(candidate/name) as other:
                assert int(base.attrs['num_levels']) == int(other.attrs['num_levels']) == 13
                assert float(base.attrs['time']) == float(other.attrs['time'])
                names = [str(base.attrs['component_'+str(k)]) for k in range(int(base.attrs['num_components']))]
                assert int(other.attrs['num_components']) == len(names) == 28
                for level in range(13):
                    x,a = read_level(base[f'level_{level}'],28)
                    y,b = read_level(other[f'level_{level}'],28)
                    assert np.array_equal(x,y), (candidate,kind,level,'coverage changed')
                    assert np.isfinite(a).all() and np.isfinite(b).all(), (candidate,kind,level,'nonfinite field')
                    diff = a.view('u8') != b.view('u8')
                    for k, field in enumerate(names):
                        rows.append(dict(case=candidate.name, clock=kind, level=level,
                            field=field, values=len(a), differing_bits=int(diff[:,k].sum()),
                            max_absolute_difference=float(np.max(np.abs(a[:,k]-b[:,k]))),
                            compared_domain='every valid evolved cell, including covered coarse cells'))
    save('t18-neutrality.csv',rows)
    for name in sorted({r['case'] for r in rows}):
        rs = [r for r in rows if r['case']==name]
        print(name, sum(r['values'] for r in rs), sum(r['differing_bits'] for r in rs))

def timestep():
    rows = []
    for directory in sorted(ROOT.glob('dt-*')):
        complete = (directory/'done.exit').exists() and (directory/'done.exit').read_text().strip()=='0'
        for path in sorted(directory.glob('t18-monitor-l*.csv')):
            row = next(csv.DictReader(path.open()))
            row.update(case=directory.name, completed=complete,
                requested_finest_window_M=.875 if directory.name.endswith('-window') else .1025390625,
                scope='RK stage inputs and endpoints; levels have different local clocks before a partial coarse stop')
            speed = max(float(row[k]) for k in ('max_light_speed','max_lapse_speed','max_shift_speed'))
            c = float(row['dt'])/float(row['h'])*speed
            row['max_characteristic_CFL'] = c
            row['frozen_scalar_wave_RK4_margin'] = np.sqrt(3)/2/c-1 if c else 0
            row['margin_scope'] = '2D fourth-order scalar-wave proxy; not full Cartoon/AMR nonlinear stability'
            row['screen_pass'] = (int(row['nonfinite_or_bad_metric'])==0 and
                int(row['chi_floor_crossings'])==int(row['lapse_floor_crossings'])==0 and
                float(row['puncture_field_max'])<1e6 and float(row['max_beta'])<2 and
                float(row['lapse_max'])<2 and float(row['chi_min'])>0 and float(row['lapse_min'])>0)
            rows.append(row)
    save('t18-dt-levels.csv',rows)

def launch():
    endpoint=.001922607421875
    rows=[]
    for directory in sorted(ROOT.glob('dt-*-window')):
        path=directory/'t18-launch-history.csv'
        if not path.exists():continue
        data=list(csv.DictReader(path.open()))
        early=[r for r in data if float(r['time_M'])<=endpoint+1e-14]
        end=[r for r in data if abs(float(r['time_M'])-endpoint)<1e-14]
        if not end:continue
        row=dict(case=directory.name,common_endpoint_M=endpoint,early_endpoints=len(early))
        for key in ('delta_Gamma_W_max','delta_lapse_W_max','delta_shift_W_max'):
            row[key+'_endpoint']=float(end[0][key])
            row[key+'_launch_peak']=max(float(r[key]) for r in early)
        rows.append(row)
    if rows:
        base=next((r for r in rows if r['case']=='dt-0.25-window'),None)
        if base:
            for row in rows:
                for key in ('delta_Gamma_W_max','delta_lapse_W_max','delta_shift_W_max'):
                    value=key+'_endpoint'
                    row[key+'_fractional_change_from_dt025']=row[value]/base[value]-1 if base[value] else 0
        save('t18-launch-dt.csv',rows)

def audit_windows():
    registration=json.loads((HERE/'t18-window-registration.json').read_text())
    pinned=json.loads((HERE/'t18-inputs.json').read_text())
    binary=next(r for r in pinned if r['path'].endswith('/native-window.ex'))
    assert hashlib.sha256(Path(binary['path']).read_bytes()).hexdigest()==binary['sha256']
    levels=list(csv.DictReader((HERE/'t18-dt-levels.csv').open()))
    audits=[]
    for job,steps in zip(registration['jobs'],registration['finest_native_steps'],strict=True):
        directory=Path(job['directory']);name=job['name']
        receipt=json.loads((directory/(name+'.resources.json')).read_text())
        assert (directory/'done.exit').read_text().strip()=='0'
        assert receipt['returncode']==receipt['child_measurement']['returncode']==0
        assert not receipt['gate_reason'] and receipt['peak_rss_bytes']<4e9
        assert json.loads(receipt['command'])==job['command']
        rows=sorted((r for r in levels if r['case']==name),key=lambda r:int(r['level']))
        assert len(rows)==13 and [int(r['level']) for r in rows]==list(range(13))
        assert all(r['completed']=='True' and r['screen_pass']=='True' for r in rows)
        history=list(csv.DictReader((directory/'t18-launch-history.csv').open()))
        finest=rows[-1];dt=float(finest['dt']);stop=float((directory/'t18-stop.txt').read_text())
        assert len(history)==steps and int(finest['samples'])==1+5*steps
        assert stop==float(finest['time_M'])==float(history[-1]['time_M'])==steps*dt
        assert stop<=registration['window_M']<stop+dt+1e-14
        assert all(abs(float(r['time_M'])-(i+1)*dt)<1e-14 for i,r in enumerate(history))
        assert all(float(r['puncture_field_max'])<1e6 for r in history)
        audits.append(dict(case=name,measured_returncode=receipt['returncode'],gate_reason=receipt['gate_reason'],
            finest_steps=steps,finest_stop_M=stop,all_13_level_screens_pass=True,
            observed_nonfinite_or_bad_metric=sum(int(r['nonfinite_or_bad_metric']) for r in rows),
            observed_chi_floor_crossings=sum(int(r['chi_floor_crossings']) for r in rows),
            observed_lapse_floor_crossings=sum(int(r['lapse_floor_crossings']) for r in rows),
            chi_min=min(float(r['chi_min']) for r in rows),lapse_min=min(float(r['lapse_min']) for r in rows),
            lapse_max=max(float(r['lapse_max']) for r in rows),max_beta=max(float(r['max_beta']) for r in rows),
            max_light_speed=max(float(r['max_light_speed']) for r in rows),
            max_lapse_speed=max(float(r['max_lapse_speed']) for r in rows),
            max_shift_speed=max(float(r['max_shift_speed']) for r in rows),
            max_characteristic_CFL=max(float(r['max_characteristic_CFL']) for r in rows),
            minimum_scalar_proxy_margin=min(float(r['frozen_scalar_wave_RK4_margin']) for r in rows),
            puncture_field_max=max(float(r['puncture_field_max']) for r in rows),
            wall_seconds=receipt['wall_seconds'],peak_RSS_bytes=receipt['peak_rss_bytes'],
            floor_scope='registered RK inputs and endpoints; final pre-clamp activation count not instrumented',
            receipt=str(directory/(name+'.resources.json'))))
    save('t18-window-audit.csv',audits)
    print('Verified independent receipts, native histories and all 39 level screens')

if __name__ == '__main__':
    neutrality()
    timestep()
    launch()
    # The already detached worker loads this analyzer after all three trials.
    # An unsuccessful level screen must prevent a successful scientific marker.
    windows=[ROOT/('dt-'+str(dt)+'-window') for dt in (.25,.375,.5)]
    if all((p/'done.exit').exists() for p in windows):
        audit_windows()
        with (HERE/'t18-dt-levels.csv').open() as f:
            levels=[r for r in csv.DictReader(f) if r['case'].endswith('-window')]
        passed=(all((p/'done.exit').read_text().strip()=='0' for p in windows) and
                len(levels)==39 and all(r['screen_pass']=='True' for r in levels))
        if not passed:
            result=dict(status='COMPLETED_WITH_FAILED_TRIALS',
                failures=[dict(case=r['case'],level=r['level']) for r in levels if r['screen_pass']!='True'],
                exits={p.name:(p/'done.exit').read_text().strip() for p in windows})
            (HERE/'t18-window-results.json').write_text(json.dumps(result,indent=2)+'\n')
            marker=ROOT.parent/'window/scientific.done.exit'
            temporary=marker.with_suffix('.tmp');temporary.write_text('1\n');temporary.replace(marker)
            raise SystemExit(1)
