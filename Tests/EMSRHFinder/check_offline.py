"""Reproduce the bounded q=.7 local checkpoint pilot and frozen controls."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

from postprocess import params, self_check
from offline import seeds

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def run(exe, p, directory, cap=240):
    directory.mkdir(parents=True)
    for name in ('chk','plt','data/extraction'):
        (directory/name).mkdir(parents=True, exist_ok=True)
    (directory/'params.txt').write_text(''.join(f'{k} = {v}\n' for k,v in p.items()))
    start = time.monotonic()
    with (directory/'run.log').open('w') as log:
        result = subprocess.run([str(exe),'params.txt'], cwd=directory,
            env=dict(os.environ, OMP_NUM_THREADS='8'), stdout=log, stderr=subprocess.STDOUT, timeout=cap)
    (directory/'status.json').write_text(json.dumps(dict(exit=result.returncode, seconds=time.monotonic()-start)))
    result.check_returncode()
    print(directory, time.monotonic()-start, flush=True)


def prepare(out, ems):
    p = params(ROOT/'Examples/EMS/params-radiation-pilot.txt')
    p.update(ems_data_path=ems/'artifacts/reference-alpha20-q0p7.trumpet',
             ems_ctt_data_path=ems/'.data/binary-ctt-q0p7/n28-r6.ctt',
             bh_charge=.7, boost_rapidity='0.2027325540540822',
             N1=192, N2=96, L=96, center='48 0', star_centre='48 0',
             max_level=3, regrid_interval='2 2 2', max_box_size=16,
             extraction_center='48 0', extraction_radii='20 24 28',
             extraction_levels='1 1 1', num_points_theta=65,
             ems_radiation_wave_level=1, ems_radiation_wave_radius=30,
             mq_extraction_center='48 0 0', mq_num_extraction_radii=2,
             mq_extraction_radii='20 28', mq_extraction_levels='0 0',
             mq_num_points_theta=32, RH_initial_centre='32 64 48',
             RH_initial_radii='0.63611642452210804 0.63611642452210804 4',
             RH_time_step_freq='400 400 400', RH_chase_speeds='1 1 1',
             ems_rh_expansion_threshold='1e-7', checkpoint_interval=1,
             plot_interval=-1, max_steps=3, stop_time=.375, verbosity=0)
    exe = next((ROOT/'Examples/EMS').glob('Main_EMSBH2DBH2d.*.ex'))
    run(exe, p, out/'binary')
    q = params(ROOT/'Examples/EMS/params-trumpet-evolve-rh.txt')
    q.update(ems_data_path=p['ems_data_path'], bh_charge=.7, N1=128, N2=64,
             L=4, center='2 0', star_centre='2 0', max_steps=0, stop_time=0,
             max_box_size=32, RH_initial_centre=2, RH_initial_radii=.64,
             RH_num_points=96, RH_time_step_freq=400, ems_test_max_updates=1000,
             ems_test_checkpoint=1, checkpoint_interval=1, verbosity=0)
    run(next(HERE.glob('EMSRHSurfaceTest2d.*.ex')), q, out/'static32')


def check(out):
    self_check()
    # Exact time matching; a missing pair must not silently use an old surface.
    p = params(out/'binary/params.txt')
    saved, skipped = seeds(p, out/'binary', .375)
    assert set(saved) == {0,1} and skipped == {2:'DORMANT'}
    try:
        seeds(p, out/'binary', .374)
    except ValueError:
        pass
    else:
        raise AssertionError('accepted a missing checkpoint-time seed')
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        for i, pair in saved.items():
            for name, row in zip((f'rh_surf_{i}.dat', f'rh_f{i}.dat'), pair):
                (directory/name).write_text(' '.join(row)+'\n')
        good = (directory/'rh_f0.dat').read_text()
        for bad in (good.replace(good.split()[1], '-1', 1),
                    good.replace(good.split()[0], '.25', 1),
                    good.replace(good.split()[1], 'nan', 1)):
            (directory/'rh_f0.dat').write_text(bad)
            try:
                seeds(p, directory, .375)
            except ValueError:
                pass
            else:
                raise AssertionError('accepted malformed history')
        (directory/'rh_f0.dat').write_text(good)
        dead = saved[0][0].copy(); dead[-1] = 'dead'
        (directory/'rh_surf_0.dat').write_text(' '.join(dead)+'\n')
        assert seeds(p, directory, .375)[1] == {0:'DEAD',2:'DORMANT'}
    print('offline input/postprocessor self-checks passed')


def results(directory):
    with (directory/'masses.csv').open() as f:
        rows = list(csv.DictReader(f))
    final = [r for r in rows if r['stage'] == '2']
    assert final and all(r['status'] == 'FOUND' for r in final)
    assert all(float(r['expansion_squared']) <= 1e-12 and float(r['theta_minus']) < 0 for r in final)
    assert all(r['equilibrium_status'] == 'BRANCH_UNVERIFIED_STATIC_PROJECTION' for r in rows)
    assert all(float(r['delta_M']) > 0 and float(r['phi_rms']) >= 0 for r in rows)
    report = []
    for index in sorted({r['search_index'] for r in final}):
        values = sorted((r for r in final if r['search_index'] == index), key=lambda r: int(r['N_theta']))
        assert [int(r['N_theta']) for r in values] == [48,96,192]
        a = [float(r['A']) for r in values]
        order = math.log2(abs((a[0]-a[1])/(a[1]-a[2])))
        assert 1.8 < order < 2.2, order  # midpoint area rule: expected order two
        native = next(r for r in rows if r['search_index'] == index and r['stage'] == '-2')
        replay = next(r for r in rows if r['search_index'] == index and r['stage'] == '-1')
        assert abs(float(replay['A'])/float(native['A'])-1) < 1e-7
        assert abs(float(replay['Q'])/float(native['Q'])-1) < 1e-7
        if math.isfinite(float(values[0]['expected_M'])):
            assert all(float(r['mass_error']) <= float(r['delta_M']) for r in values)
            assert all(abs(float(r['phi_mean'])-.052274528945421062) < 1e-7 for r in values)
        report.append(dict(search_index=index, area_order=order,
                           native={k:native[k] for k in ('A','Q','M_eq')},
                           strict=[{k:r[k] for k in ('N_theta','A','Q','M_eq','delta_M','phi_mean','phi_rms','expansion_squared')} for r in values]))
    (directory/'checks.json').write_text(json.dumps(dict(passed=True, surfaces=report), indent=2)+'\n')
    print(json.dumps(report, indent=2))


def replay(directory):
    """Full-precision final shapes survive the production checkpoint reader."""
    checks = []
    for source in sorted(directory.glob('n[0-9]*')):
        p = params(source/'params.txt')
        n = int(source.name[1:])
        p['RH_num_points'] = ' '.join([str(n)]*int(p['RH_num_horizons']))
        with (source/'surfaces.csv').open() as f:
            final = {int(r['search_index']):r for r in csv.DictReader(f) if r['stage'] == '2'}
        assert final
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp)
            for i, r in final.items():
                shape = (source/f'shape-{i}-2.dat').read_text().split()
                d = (source/f'rh_surf_{i}.dat').read_text().split()
                d[0], d[2], d[4], d[9], d[17], d[-1] = shape[0],shape[1],r['A'],r['expansion_squared'],r['Q'],'found'
                (destination/f'rh_surf_{i}.dat').write_text(' '.join(d)+'\n')
                (destination/f'rh_f{i}.dat').write_text(' '.join([shape[0],*shape[2:]])+'\n')
            # Exercise the largest valid driver cap without doing a long chase.
            p.update(offline_max_updates=1,offline_seconds=1795)
            (destination/'params.txt').write_text(''.join(f'{k} = {v}\n' for k,v in p.items()))
            with (destination/'run.log').open('w') as f:
                code = subprocess.run([str(next(HERE.glob('EMSRHCheckpoint2d.*.ex'))),'params.txt'],
                        cwd=destination, env=dict(os.environ,OMP_NUM_THREADS='8'), stdout=f, stderr=f, timeout=30).returncode
            assert code in (0,1), (destination/'run.log').read_text()
            assert 'GRAMRLevel::advance' not in (destination/'run.log').read_text()
            with (destination/'surfaces.csv').open() as f:
                restored = [r for r in csv.DictReader(f) if r['stage'] == '-1']
            assert len(restored) == len(final)
            for r in restored:
                expected = final[int(r['search_index'])]
                for key in ('A','Q','phi_mean','phi_rms','theta_minus'):
                    assert abs(float(r[key])-float(expected[key])) < 1e-12, (key,r,expected)
                assert abs(float(r['expansion_squared'])-float(expected['expansion_squared'])) < 1e-20
        checks.append(dict(N_theta=n, surfaces=len(final)))
    assert checks
    (directory/'replay-check.json').write_text(json.dumps(dict(passed=True, checks=checks,
        absolute_tolerance=1e-12, residual_absolute_tolerance=1e-20), indent=2)+'\n')
    print('full-precision checkpoint replay passed')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode', choices=('prepare','check','results','replay'))
    ap.add_argument('output', type=Path)
    ap.add_argument('--ems-root', type=Path, default=Path('/Users/auroradysis/Workspace/EMS'))
    a = ap.parse_args()
    if a.mode == 'prepare': prepare(a.output.resolve(), a.ems_root.resolve())
    elif a.mode == 'check': check(a.output.resolve())
    elif a.mode == 'results': results(a.output.resolve())
    else: replay(a.output.resolve())
