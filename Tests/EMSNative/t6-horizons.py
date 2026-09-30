#!/usr/bin/env python3
"""One bounded frozen-checkpoint find per invocation; no evolution/static input."""
import argparse
import csv
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import h5py

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'EMSRHFinder'))
from offline import history, params

INPUT = Path('/Users/auroradysis/Workspace/EMS/.data/exp-0019')
PARAMS = Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0019/submissions/exp-0019')
OUT = Path('/private/tmp/ems-t6/horizons')


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(name, rows):
    with (HERE/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)


def checkpoint_info(path):
    with h5py.File(path,'r') as f:
        levels = [f[f'level_{i}'] for i in range(int(f.attrs['num_levels']))]
        t = float(levels[0].attrs['time'])
        times = np.array([float(g.attrs['time']) for g in levels])
        spread = float(np.max(np.abs(times-t)))
        assert np.isfinite(times).all() and t >= 0
        assert spread <= 64*np.finfo(float).eps*max(1.,abs(t)), (path,times)
        return dict(time=t,level_time_spread_M=spread)


def run(scale, checkpoint, n):
    root = INPUT/scale
    cp = root/'chk'/checkpoint
    info = checkpoint_info(cp)
    parameter = PARAMS/f'params-{scale}.txt'
    surf, shape = root/'rh_surf_0.dat', root/'rh_f0.dat'
    pairs = {float(r[0]): r for r in history(shape)}
    rows = history(surf)
    if info['time'] == 0:
        seed = rows[0]  # Explicit numerical seed only; fields always come from t=0.
    else:
        match = [r for r in rows if float(r[0]) == float(format(info['time'], '.8e'))]
        assert match, (scale, checkpoint, info['time'])
        seed = match[-1]
    f = pairs[float(seed[0])]
    assert len(seed)==20 and len(f)==97 and seed[-1] in ('found','close','far')
    assert np.isfinite(np.array(f, float)).all() and min(map(float,f[1:])) > 0
    d = OUT/scale/cp.stem/f'n{n}'
    d.mkdir(parents=True)
    for name, row in [('rh_surf_0.dat',seed), ('rh_f0.dat',f)]:
        (d/name).write_text(' '.join([format(info['time'],'.17g'), *row[1:]])+'\n')
    p = params(parameter)
    for key in ('output_path','hdf5_subpath','data_subpath','pout_subpath'):
        p.pop(key,None)
    p.update(restart_file=str(cp), checkpoint_interval=-1, plot_interval=-1, verbosity=0,
             offline_points=n, offline_max_updates=3000, offline_seconds=115,
             RH_initial_radii=repr(float(np.mean(np.array(f[1:],float)))), RH_initial_centre=seed[2],
             ems_data_path='/private/tmp/ems-t6/ABSENT-STATIC-INPUT', t6_checkpoint_diagnostics='true')
    (d/'params.txt').write_text(''.join(f'{k} = {v}\n' for k,v in p.items()))
    exe = next((HERE.parent/'EMSRHFinder').glob('EMSRHCheckpoint2d.*.ex')).resolve()
    inputs = [cp,parameter,surf,shape,exe,Path('Source/RHFinder/RHSurf.hpp'),Path('Source/RHFinder/RHUnion.hpp')]
    before = {str(p.resolve()):sha(p) for p in inputs}
    start = time.monotonic()
    with (d/'run.log').open('w') as log:
        try:
            result = subprocess.run([str(exe),'params.txt'],cwd=d,stdout=log,stderr=subprocess.STDOUT,
                                    env=dict(os.environ,OMP_NUM_THREADS='4'),timeout=120).returncode
        except subprocess.TimeoutExpired:
            result = 124
    state = dict(scale=scale,checkpoint=str(cp),time_M=info['time'],N_theta=n,
                 level_time_spread_M=info['level_time_spread_M'],
                 seed_time_M=float(seed[0]),seed_mode=seed[-1],seed_source='EARLIEST_NUMERICAL_SHAPE' if info['time']==0 else 'SAME_TIME_NUMERICAL_SHAPE',
                 exit=result,wall_seconds=time.monotonic()-start,
                 peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,sha256=before)
    (d/'status.json').write_text(json.dumps(state,indent=2)+'\n')
    assert before == {str(p.resolve()):sha(p) for p in inputs}
    print(json.dumps({k:v for k,v in state.items() if k!='sha256'}),flush=True)
    if (d/'surfaces.csv').exists():
        print((d/'surfaces.csv').read_text().splitlines()[-1],flush=True)
    return result


def report():
    histories, finals, stages, spheres, audits = [], [], [], [], []
    inputs = []
    for scale in ('E-low','E-mid'):
        for v in history(INPUT/scale/'rh_surf_0.dat'):
            err, mean = float(v[9]),float(v[7])
            assert err >= 0 and math.sqrt(err)+1e-8 >= abs(mean)
            histories.append(dict(scale=scale,time_M=float(v[0]),A=float(v[4]),Q=float(v[17]),
                                  theta_plus_mean=mean,theta_plus_mean_square=err,
                                  theta_plus_rms=math.sqrt(err),mode=v[-1]))
        for cp in sorted((INPUT/scale/'chk').glob('*.hdf5')):
            measurements = []
            for n in (48,96):
                d = OUT/scale/cp.stem/f'n{n}'
                status = json.loads((d/'status.json').read_text())
                status.setdefault('level_time_spread_M', checkpoint_info(cp)['level_time_spread_M'])
                assert status['exit'] == 0
                audits.append({k:v for k,v in status.items() if k!='sha256'})
                for path,digest in status['sha256'].items():
                    inputs.append(dict(scale=scale,checkpoint=cp.name,N_theta=n,
                                       input_path=path,sha256=digest))
                data = list(csv.DictReader((d/'surfaces.csv').open()))
                for r in data:
                    stages.append(dict(scale=scale,checkpoint=cp.name,
                                       **r,theta_plus_rms=math.sqrt(float(r['expansion_squared']))))
                last = data[-1]
                measurements.append(last)
                for r in csv.DictReader((d/'fixed-spheres.csv').open()):
                    spheres.append(dict(scale=scale,checkpoint=cp.name,**r))
            a,b = measurements
            strict = (a['status']==b['status']=='FOUND' and a['stage']==b['stage']=='2')
            assert strict and max(float(a['expansion_squared']),float(b['expansion_squared'])) <= 1e-12
            directory = OUT/scale/cp.stem/'n96'
            all96 = list(csv.DictReader((directory/'surfaces.csv').open()))
            previous = next((r for r in all96 if r['stage']=='1'),None)
            seed = all96[0]
            finals.append(dict(scale=scale,checkpoint=cp.name,time_M=float(b['time']),status='QUALIFIED_FROZEN_GRID' if strict else 'INCOMPLETE',
                               A48=float(a['A']),A96=float(b['A']),Q48=float(a['Q']),Q96=float(b['Q']),
                               angular_delta_A=abs(float(a['A'])-float(b['A'])),
                               angular_delta_Q=abs(float(a['Q'])-float(b['Q'])),
                               stopping_delta_A=abs(float(previous['A'])-float(b['A'])) if previous else math.nan,
                               stopping_delta_Q=abs(float(previous['Q'])-float(b['Q'])) if previous else math.nan,
                               replay_A=float(seed['A']),replay_Q=float(seed['Q']),
                               replay_theta_plus_rms=math.sqrt(float(seed['expansion_squared'])),
                               theta_plus_rms48=math.sqrt(float(a['expansion_squared'])),
                               theta_plus_rms96=math.sqrt(float(b['expansion_squared'])),
                               theta_minus96=float(b['theta_minus'])))
    save('t6-finder-history.csv',histories);save('t6-finder-stages.csv',stages)
    save('t6-finder-final.csv',finals);save('t6-fixed-spheres.csv',spheres);save('t6-finder-runs.csv',audits)
    save('t6-finder-input-audit.csv',inputs)
    drifts = []
    for scale in ('E-low','E-mid'):
        h = sorted([r for r in finals if r['scale']==scale],key=lambda r:r['time_M'])
        for n in (48,96):
            akey,qkey='A'+str(n),'Q'+str(n)
            for r in h:
                drifts.append(dict(scale=scale,time_M=r['time_M'],N_theta=n,surface='HORIZON',radius_M=math.nan,
                                   A=r[akey],Q=r[qkey],baseline_time_M=0.,
                                   relative_A_from_t0=(r[akey]-h[0][akey])/h[0][akey],
                                   relative_Q_from_t0=(r[qkey]-h[0][qkey])/h[0][qkey],
                                   relative_Q_from_first_late=(r[qkey]-h[1][qkey])/h[1][qkey] if r['time_M']>=h[1]['time_M'] else math.nan))
            for radius in (.02,.05,.1):
                q = sorted([r for r in spheres if r['scale']==scale and int(r['N_theta'])==n and float(r['radius'])==radius],key=lambda r:float(r['time']))
                for r in q:
                    drifts.append(dict(scale=scale,time_M=float(r['time']),N_theta=n,surface='FIXED_SPHERE',radius_M=radius,
                                       A=float(r['A']),Q=float(r['Q']),baseline_time_M=0.,
                                       relative_A_from_t0=(float(r['A'])-float(q[0]['A']))/float(q[0]['A']),
                                       relative_Q_from_t0=(float(r['Q'])-float(q[0]['Q']))/float(q[0]['Q']),
                                       relative_Q_from_first_late=(float(r['Q'])-float(q[1]['Q']))/float(q[1]['Q']) if float(r['time'])>=float(q[1]['time']) else math.nan))
    save('t6-charge-drift.csv',drifts)
    native = []
    for scale in ('E-low','E-mid'):
        data = [r for r in histories if r['scale']==scale]
        times = [r['time_M'] for r in data]
        for start,end in ((2.,10.),(data[0]['time_M'],data[-1]['time_M'])):
            q0,q1 = np.interp([start,end],times,[r['Q'] for r in data])
            native.append(dict(scale=scale,start_time_M=start,end_time_M=end,
                               Q_start=q0,Q_end=q1,relative_Q_change=(q1-q0)/q0,
                               qualification='UNQUALIFIED_HISTORY_LINEAR_INTERPOLATION'))
    save('t6-native-drift.csv',native)
    print('Finder history, residual/angle sensitivity and numerical-baseline drift exported.')


if __name__=='__main__':
    ap = argparse.ArgumentParser();ap.add_argument('action',choices=('run','report'))
    ap.add_argument('scale',nargs='?');ap.add_argument('checkpoint',nargs='?');ap.add_argument('n',nargs='?',type=int)
    args=ap.parse_args()
    if args.action=='run':
        assert args.scale in ('E-low','E-mid') and args.n in (48,96)
        sys.exit(run(args.scale,args.checkpoint,args.n))
    report()
