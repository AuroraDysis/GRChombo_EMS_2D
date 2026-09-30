#!/usr/bin/env python3
"""Frozen T6 finder criteria, two threads; missing rows use preceding saved shape."""
import argparse, importlib.util, json, os, resource, subprocess, time
from pathlib import Path
import numpy as np
spec = importlib.util.spec_from_file_location('t6_horizons', Path(__file__).with_name('t6-horizons.py'))
t6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t6)
OUT = Path('/private/tmp/ems-t7e/horizons')
def run(checkpoint,n):
    root=t6.INPUT/'E-high'; cp=root/'chk'/checkpoint
    info=t6.checkpoint_info(cp); target=info['time']
    surf,shape=root/'rh_surf_0.dat',root/'rh_f0.dat'
    rows=t6.history(surf); pairs={float(r[0]):r for r in t6.history(shape)}
    exact=[r for r in rows if float(r[0])==float(format(target,'.8e'))]
    seed=exact[-1] if exact else (rows[0] if target==0 else max((r for r in rows if float(r[0])<target),key=lambda r:float(r[0])))
    f=pairs[float(seed[0])]
    assert len(seed)==20 and len(f)==97 and np.isfinite(np.array(f,float)).all() and min(map(float,f[1:]))>0
    d=OUT/'E-high'/cp.stem/f'n{n}';d.mkdir(parents=True)
    for name,row in [('rh_surf_0.dat',seed),('rh_f0.dat',f)]:
        (d/name).write_text(' '.join([format(target,'.17g'),*row[1:]])+'\n')
    parameter=t6.PARAMS/'params-E-high.txt';p=t6.params(parameter)
    for key in ('output_path','hdf5_subpath','data_subpath','pout_subpath'):p.pop(key,None)
    p.update(restart_file=str(cp),checkpoint_interval=-1,plot_interval=-1,verbosity=0,
        offline_points=n,offline_max_updates=3000,offline_seconds=115,
        RH_initial_radii=repr(float(np.mean(np.array(f[1:],float)))),RH_initial_centre=seed[2],
        ems_data_path='/private/tmp/ems-t6/ABSENT-STATIC-INPUT',t6_checkpoint_diagnostics='true')
    (d/'params.txt').write_text(''.join(f'{k} = {v}\n' for k,v in p.items()))
    exe=next((t6.HERE.parent/'EMSRHFinder').glob('EMSRHCheckpoint2d.*.ex')).resolve()
    inputs=[cp,parameter,surf,shape,exe,Path('Source/RHFinder/RHSurf.hpp'),Path('Source/RHFinder/RHUnion.hpp')]
    before={str(p.resolve()):t6.sha(p) for p in inputs};start=time.monotonic()
    with (d/'run.log').open('w') as log:
        try:code=subprocess.run([str(exe),'params.txt'],cwd=d,stdout=log,stderr=subprocess.STDOUT,
            env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'),timeout=120).returncode
        except subprocess.TimeoutExpired:code=124
    state=dict(scale='E-high',checkpoint=str(cp),time_M=target,N_theta=n,
        level_time_spread_M=info['level_time_spread_M'],seed_time_M=float(seed[0]),seed_mode=seed[-1],
        seed_source='SAME_TIME_NUMERICAL_SHAPE' if exact else ('EARLIEST_NUMERICAL_SHAPE' if target==0 else 'PRECEDING_NUMERICAL_SHAPE'),
        exit=code,wall_seconds=time.monotonic()-start,peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,sha256=before)
    (d/'status.json').write_text(json.dumps(state,indent=2)+'\n')
    assert before=={str(p.resolve()):t6.sha(p) for p in inputs}
    print(json.dumps({k:v for k,v in state.items() if k!='sha256'}),flush=True)
    if (d/'surfaces.csv').exists():print((d/'surfaces.csv').read_text().splitlines()[-1],flush=True)
    return code
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('checkpoint');ap.add_argument('n',type=int,choices=(48,96));a=ap.parse_args()
    raise SystemExit(run(a.checkpoint,a.n))
