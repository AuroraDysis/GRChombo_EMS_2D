"""Serial, bounded t=0 finds. Each invocation is independently capped at 240 s."""
from pathlib import Path
import json
import math
import hashlib
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXE = next(HERE.glob('EMSRHSurfaceTest2d.*.ex'))
EXE_SHA256 = hashlib.sha256(EXE.read_bytes()).hexdigest()


def parameters(case, n, resolution, threshold='1e-7', duplicates=False):
    p = {}
    for line in (ROOT/'Examples/EMS/params-trumpet-evolve-rh.txt').read_text().splitlines():
        if '=' in line and not line.startswith('#'):
            k, v = line.split('=', 1)
            p[k.strip()] = v.strip()
    rn = case == 'rn'
    p.update(N1=str(4*resolution), N2=str(2*resolution), L='4', center='2 0', star_centre='2 0',
             max_steps='0', stop_time='0', max_box_size='64', RH_num_points=str(n),
             RH_time_step_freq='400', RH_initial_radii='0.575' if rn else '0.64',
             RH_initial_centre='2', ems_test_max_updates='1000', verbosity='0',
             ems_test_expansion_squared=threshold,
             boosted='true' if case in ('plus', 'minus') else 'false',
             boost_rapidity='-0.20273255' if case == 'minus' else '0.20273255' if case == 'plus' else '0')
    if rn:
        p['ems_data_path'] = str(HERE/'fixtures/rn.trumpet')
    if duplicates:
        p['RH_num_horizons'] = '2'
        for key in tuple(p):
            if key.startswith('RH_') and key not in ('RH_num_horizons', 'RH_activate'):
                p[key] += ' ' + p[key]
    p.update(plot_prefix='plt/EMS_', plot_interval='1', num_plot_vars='13',
             plot_vars='chi h11 h12 h22 hww K A11 A12 A22 Aww phi Ex Ey')
    return p


def run(root, case, n, resolution, threshold='1e-7', duplicates=False, **changes):
    path = root/f'{case}-n{n}-dx{resolution}'
    path.mkdir(parents=True, exist_ok=True)
    (path/'plt').mkdir(exist_ok=True)
    p = parameters(case, n, resolution, threshold, duplicates)
    p.update(changes)
    text='\n'.join(f'{k} = {v}' for k, v in p.items())+'\n'
    if (path/'status.json').exists() and (path/'params.txt').read_text()==text:
        prior=json.loads((path/'status.json').read_text())
        if prior['exit']==0 and prior.get('exe_sha256')==EXE_SHA256:
            print(json.dumps(dict(prior,reused=True)),flush=True)
            return 0
    (path/'params.txt').write_text(text)
    start = time.monotonic()
    with (path/'run.log').open('w') as log:
        try:
            status = subprocess.run([str(EXE), 'params.txt'], cwd=path, stdout=log, stderr=subprocess.STDOUT,
                                    env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=240).returncode
        except subprocess.TimeoutExpired:
            status = 124
    if 'ems_test_initial_shape' in p and 'EMSRH_TEST_INITIAL_SHAPE ' not in (path/'run.log').read_text():
        status = 2
    record = dict(case=case, N_theta=n, resolution=resolution, threshold=threshold,
                  seconds=time.monotonic()-start, exit=status, exe_sha256=EXE_SHA256)
    (path/'status.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record), flush=True)
    return status


def killing_seed(root, case, n, resolution):
    """Known stationary control: X_rest=cosh(eta)*x, r_rest=R_h.

    This is an initial guess for the unchanged finder, not a replacement finder.
    Keep these measurements separate from the unseeded native search gates.
    """
    root.mkdir(parents=True,exist_ok=True)
    radius=0.5701934158668871 if case=='rn' else 0.63593977642346233
    gamma=math.cosh(0.20273255) if case in ('plus','minus') else 1.
    shape=[radius/math.sqrt(gamma**2*math.cos((j+0.5)*math.pi/n)**2+
                           math.sin((j+0.5)*math.pi/n)**2) for j in range(n)]
    path=root/f'{case}-n{n}.shape'
    path.write_text('0 '+' '.join(format(f,'.17g') for f in shape)+'\n')
    return run(root,case,n,resolution,ems_test_initial_shape=str(path.resolve()))


def verify(root, source):
    """Replay saved final geometry through the author's public interpolation/expansion."""
    records=[]
    for directory in sorted(source.glob('*-n*-dx*')):
        status=json.loads((directory/'status.json').read_text())
        if status['exit'] != 0: raise RuntimeError(f'failed source run: {directory}')
        shape=directory/'rh_f0.dat'
        rows=[line.split() for line in (directory/'rh_surf_0.dat').read_text().splitlines()
              if line.strip() and not line.lstrip().startswith('#')]
        if rows[-1][-1]!='found': raise RuntimeError('source surface not found')
        changes={k.strip():v.split('#')[0].strip()
                 for line in (directory/'params.txt').read_text().splitlines()
                 if '=' in line and not line.startswith('#') for k,v in [line.split('=',1)]}
        if changes['RH_num_horizons']!='1': raise ValueError('single-surface replay only')
        changes.update(ems_test_initial_shape=str(shape),RH_initial_centre=rows[-1][2],ems_test_verify_only='true')
        code=run(root,status['case'],status['N_theta'],status['resolution'],**changes)
        log=(root/directory.name/'run.log').read_text().splitlines()
        values=[list(map(float,line.split()[1:])) for line in log if line.startswith('EMSRH_FRESH ')]
        if len(values)!=1: raise RuntimeError(f'missing fresh expansion: {directory}')
        A,Q,plus,minus,error=values[0]
        records.append(dict(source=str(directory),shape_sha256=hashlib.sha256(shape.read_bytes()).hexdigest(),
                            A=A,Q=Q,theta_plus=plus,theta_minus=minus,expansion_squared=error,
                            native_expansion_squared=float(rows[-1][9]),exit=code,exe_sha256=EXE_SHA256))
    (root/'fresh.json').write_text(json.dumps(records,indent=2)+'\n')
    return int(any(r['exit'] for r in records))


if __name__ == '__main__':
    root = Path(sys.argv[1]).resolve()
    mode = sys.argv[2]
    status = 0
    if mode == 'verify':
        status=verify(root,Path(sys.argv[3]).resolve())
    elif mode == 'pilot':
        for case, n, dx in [('plus',48,32), ('minus',48,32), ('rn',48,32), ('plus',192,64)]:
            status |= run(root,case,n,dx)
    elif mode in ('matrix', 'tight', 'killing'):
        for case in ('static','plus','minus','rn'):
            for n in (48,96,192):
                for dx in ((32,64,128) if mode=='killing' else (32,64)):
                    if mode=='killing':
                        status |= killing_seed(root,case,n,dx)
                        continue
                    changes={}
                    if mode=='tight':
                        native=Path(sys.argv[3]).resolve()/f'{case}-n{n}-dx{dx}'
                        last=[line.split() for line in (native/'rh_surf_0.dat').read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')][-1]
                        if last[-1]!='found': raise RuntimeError('seed surface not found')
                        changes=dict(ems_test_initial_shape=str(native/'rh_f0.dat'), RH_initial_centre=last[2], RH_chase_speeds='0.25')
                    status |= run(root,case,n,dx, '1e-12' if mode=='tight' else '1e-7', **changes)
    else:
        raise ValueError(mode)
    (root/'complete.json').write_text(json.dumps({'exit':status})+'\n')
    sys.exit(status)
