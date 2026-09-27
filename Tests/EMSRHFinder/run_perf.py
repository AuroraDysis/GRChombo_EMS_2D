"""Bounded local identity/timing runs; no MPI, SSH or build side effects.

python3 run_perf.py OUTPUT BASELINE_TREE {identity,restart,boxes,weyl,final}
Build baseline and candidate completely first. Raw outputs stay in OUTPUT.
Each process has a 240-second cap; timing comparisons run sequentially.
"""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
EMS = Path('/Users/auroradysis/Workspace/EMS')
OUT = Path(sys.argv[1]).resolve()
BASE = Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True, exist_ok=True)


def exe(root, directory, name):
    return next((root/directory).glob(name+'2d.*.ex'))


def params(file):
    return {k.strip(): v.split('#')[0].strip()
            for line in file.read_text().splitlines()
            if '=' in line and not line.startswith('#')
            for k, v in [line.split('=', 1)]}


def pilot():
    p = params(ROOT/'Examples/EMS/params-radiation-pilot.txt')
    p.update(ems_data_path=EMS/'artifacts/reference-alpha20-qfile.trumpet',
             ems_ctt_data_path=EMS/'.data/binary-ctt-t6/n28-r6.ctt',
             verbosity=0, plot_interval=-1, checkpoint_interval=1)
    return p


def run(name, binary, p, threads, *extra):
    path = OUT/name
    path.mkdir(parents=True, exist_ok=True)
    for directory in ('chk', 'plt', 'data/extraction'):
        (path/directory).mkdir(parents=True, exist_ok=True)
    (path/'params.txt').write_text(''.join(f'{k} = {v}\n' for k, v in p.items()))
    start = time.monotonic()
    with (path/'run.log').open('w') as log:
        try:
            status = subprocess.run([str(binary), 'params.txt', *extra], cwd=path,
                stdout=log, stderr=subprocess.STDOUT, timeout=240,
                env=dict(os.environ, OMP_NUM_THREADS=str(threads), CH_TIMER='1')).returncode
        except subprocess.TimeoutExpired:
            status = 124
    result = dict(name=name, threads=threads, seconds=time.monotonic()-start,
                  exit=status, exe_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
    (path/'status.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    if status:
        raise RuntimeError(f'{name} failed: {status}; see {path}/run.log')


def reduced():
    p = pilot()
    p.update(N1=192, N2=96, L=96, center='48 0', star_centre='48 0',
             max_level=3, regrid_interval='2 2 2', max_box_size=32,
             extraction_center='48 0', extraction_radii='20 24 28',
             extraction_levels='1 1 1', num_points_theta=65,
             ems_radiation_wave_level=1, ems_radiation_wave_radius=30,
             mq_extraction_center='48 0 0', mq_num_extraction_radii=2,
             mq_extraction_radii='20 28', mq_extraction_levels='0 0',
             mq_num_points_theta=32, RH_initial_centre='32 64 48',
             max_steps=3, stop_time=.375)
    return p


def identity():
    p = params(ROOT/'Examples/EMS/params-trumpet-evolve-rh.txt')
    p.update(N1=128, N2=64, L=4, center='2 0', star_centre='2 0',
             RH_initial_centre=2, RH_initial_radii=.64, RH_num_points=96,
             RH_time_step_freq=400, RH_chase_speeds=.125,
             ems_rh_expansion_threshold='1e-12', ems_test_max_updates=1000,
             verbosity=0, max_box_size=32, plot_interval=1)
    for label, root in (('baseline', BASE), ('candidate', ROOT)):
        for threads in (1, 8):
            run(f'{label}-rh-{threads}', exe(root, 'Tests/EMSRHFinder', 'EMSRHSurfaceTest'), p, threads)
            run(f'{label}-pilot-{threads}', exe(root, 'Examples/EMS', 'Main_EMSBH2DBH'), reduced(), threads)


def restart():
    source = OUT/'candidate-pilot-8'
    for threads in (1, 8):
        name = f'restart-{threads}'
        shutil.copytree(source, OUT/name)
        p = params(source/'params.txt')
        p['restart_file'] = 'chk/EMS_000002.2d.hdf5'
        run(name, exe(ROOT, 'Examples/EMS', 'Main_EMSBH2DBH'), p, threads)
        log = (OUT/name/'run.log').read_text()
        assert log.count('restored centre and 96 radii from t = 0.25') == 2, log[-2000:]
        assert 'surface 2 restored' not in log  # still dormant before t=70
        for i in (0, 1):
            for filename in (f'rh_surf_{i}.dat', f'rh_f{i}.dat'):
                def rows(path):
                    return [line for line in path.read_text().splitlines()
                            if line.strip() and not line.lstrip().startswith('#')]
                before, after = rows(source/filename), rows(OUT/name/filename)
                prefix = [line for line in before if float(line.split()[0]) <= .25]
                assert after[:len(prefix)] == prefix
                times = [float(line.split()[0]) for line in after]
                assert times == sorted(set(times)), (filename, times)
                assert times[-1] == .375


def boxes():
    # Initialization only: no evolution or finder chase. The full production
    # tagging/radiation floors generate the actual level-by-level box layouts.
    for size in (64, 32, 16):
        p = pilot()
        p.update(max_steps=0, stop_time=0, RH_activate=0, max_box_size=size)
        run(f'boxes-{size}', exe(ROOT, 'Examples/EMS', 'Main_EMSBH2DBH'), p, 8)


def weyl():
    p = reduced()
    p.update(max_box_size=16, RH_activate=0, write_extraction=1)
    for label, root, threads in (('baseline', BASE, 1), ('candidate', ROOT, 8)):
        run(f'weyl-{label}', exe(root, 'Examples/EMS', 'Main_EMSBH2DBH'), p, threads)


def final():
    p = reduced()
    p.update(max_box_size=16, write_extraction=1)
    for label, root, threads in (('baseline', BASE, 1), ('candidate', ROOT, 1),
                                  ('candidate', ROOT, 8)):
        run(f'final-{label}-{threads}', exe(root, 'Examples/EMS', 'Main_EMSBH2DBH'), p, threads)


{'identity': identity, 'restart': restart, 'boxes': boxes, 'weyl': weyl,
 'final': final}[sys.argv[3]]()
