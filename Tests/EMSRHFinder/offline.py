"""Batch frozen checkpoint finds. Run --help; all outputs go to a fresh directory."""
import argparse
import csv
import hashlib
import json
import math
import os
import resource
import sys
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np
from postprocess import coincident, params

HERE = Path(__file__).resolve().parent


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def history(path):
    return [line.split() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def seeds(p, directory, t):
    """Exact checkpoint time at the native writer's nine significant digits."""
    saved, skipped = {}, {}
    starts = list(map(float, p.get('RH_start_times', '0 '*int(p['RH_num_horizons'])).split()))
    counts = list(map(int, p['RH_num_points'].split()))
    target = float(format(t, '.8e'))
    for i, (start, count) in enumerate(zip(starts, counts, strict=True)):
        if t < start:
            skipped[i] = 'DORMANT'
            continue
        paths = [directory/f'rh_surf_{i}.dat', directory/f'rh_f{i}.dat']
        rows = [history(path) for path in paths]
        matches = [[r for r in rs if float(r[0]) == target] for rs in rows]
        # A previously dead surface can have no new row; preserve it, never revive it.
        if not all(matches):
            before = [[r for r in rs if float(r[0]) <= target] for rs in rows]
            if all(before) and before[0][-1][-1] == 'dead':
                matches = [[rs[-1]] for rs in before]
            else:
                raise ValueError(f'surface {i}: missing checkpoint-time rh_surf/rh_f pair at {t}')
        d, f = (rs[-1] for rs in matches)
        if (len(d) != 20 or len(f) != count+1 or int(d[1]) != i or d[0] != f[0]
                or d[-1] not in ('found','close','medium','far','dead','dormnt')
                or not all(math.isfinite(float(x)) for x in d[:-1]+f)
                or any(float(x) <= 0 for x in f[1:])):
            raise ValueError(f'surface {i}: invalid or mismatched saved geometry')
        saved[i] = (d, f)
        if d[-1] == 'dead':
            skipped[i] = 'DEAD'
    if len(starts) != int(p['RH_num_horizons']):
        raise ValueError('RH parameter counts differ')
    return saved, skipped


def checkpoint_info(path):
    with h5py.File(path, 'r') as f:
        levels = sorted((k for k in f if k.startswith('level_')), key=lambda k: int(k[6:]))
        t = float(f['level_0'].attrs['time'])
        if not math.isfinite(t) or t < 0 or any(float(f[k].attrs['time']) != t for k in levels):
            raise ValueError('checkpoint times must be finite, nonnegative and synchronized')
        return dict(time=t, levels=len(levels), dx=float(f['level_0'].attrs['dx']),
                    field_bytes=sum(f[k]['data:datatype=0'].size*f[k]['data:datatype=0'].dtype.itemsize for k in levels))


def collect(runs, p, saved, args, output):
    rows = []
    for directory in runs:
        with (directory/'surfaces.csv').open() as f:
            for raw in csv.DictReader(f):
                r = {k: (v if k == 'status' else float(v)) for k, v in raw.items()}
                for k in ('N_theta', 'search_index', 'stage', 'updates'):
                    r[k] = int(r[k])
                if r['stage'] == -1 and directory != runs[0]:
                    continue
                shape = np.loadtxt(directory/f"shape-{r['search_index']}-{r['stage']}.dat", ndmin=1)
                if (len(shape) != r['N_theta']+2 or not np.isfinite(shape).all()
                        or shape[0] != r['time'] or shape[1] != r['centre'] or (shape[2:] <= 0).any()):
                    raise ValueError('full-precision shape does not match its measurement')
                r['_shape'] = (shape[2:], shape[1])
                rows.append(r)
                if r['stage'] == -1:
                    d = saved[r['search_index']][0]
                    native = dict(r, stage=-2, status='INRUN_'+d[-1].upper(), A=float(d[4]),
                                  Q=float(d[17]), expansion_squared=float(d[9]))
                    rows.append(native)
    for r in rows:
        n, stage, index = (r[k] for k in ('N_theta', 'stage', 'search_index'))
        prior = next((x for x in rows if x['N_theta'] == n and x['search_index'] == index
                      and x['stage'] == stage-1), None) if stage > 0 else None
        angular = sorted((x for x in rows if x['search_index'] == index and x['stage'] == stage
                          and x['status'] == 'FOUND'),
                         key=lambda x: x['N_theta'])
        neighbor = next((x for x in angular if x['N_theta'] == 2*n), None)
        neighbor = neighbor or next((x for x in angular if x['N_theta']*2 == n), None)
        bias = (math.pi/(2*n))/math.sin(math.pi/(2*n))-1
        for key in ('A', 'Q'):
            stop = abs(r[key]-prior[key]) if prior else 0.
            quad = abs(r[key])*bias
            if neighbor:
                quad = max(quad, abs(r[key]-neighbor[key])*(4/3 if neighbor['N_theta'] > n else 1/3))
            r['stopping_delta_'+key] = stop
            r['angular_delta_'+key] = quad
            r['measurement_delta_'+key] = stop+quad
        r.update(ems_alpha=float(p['ems_alpha']), f0=float(p.get('ems_f0', 0)),
                 f1=float(p.get('ems_f1', 0)), f2=float(p.get('ems_f2', 0)),
                 phi_inf=float(p.get('ems_radiation_phi_inf', 0)), branch='scalarized-positive-n0', duplicate_count=1,
                 duplicate_A_Q=f"{r['A']:.17g}:{r['Q']:.17g}", spread_A=0., spread_Q=0.,
                 expected_M=args.expected_mass, expected_A=args.expected_area, expected_Q=args.expected_charge,
                 error_scope=('ANGULAR_STOPPING_TABLE_ONLY_GRID_UNMEASURED' if stage > 0 and r['status'] == 'FOUND'
                              else 'DIAGNOSTIC_ONLY'))
    unique = []
    for r in rows:
        same = next((x for x in unique if (x['N_theta'],x['stage'],x['status']) == (r['N_theta'],r['stage'],r['status'])
                     and coincident(x['_shape'], r['_shape'])), None)
        if same:
            same['duplicate_count'] += 1
            same['duplicate_A_Q'] += ';'+r['duplicate_A_Q']
            values = [list(map(float, pair.split(':'))) for pair in same['duplicate_A_Q'].split(';')]
            same['spread_A'] = max(v[0] for v in values)-min(v[0] for v in values)
            same['spread_Q'] = max(v[1] for v in values)-min(v[1] for v in values)
        else:
            unique.append(r)
    if not unique:
        return
    for r in unique:
        del r['_shape']
    with (output/'measurements.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=unique[0])
        writer.writeheader(); writer.writerows(unique)
    with (output/'julia.log').open('w') as log:
        subprocess.run(['julia', '--project=test', str(HERE/'postprocess.jl'),
                        str(output/'measurements.csv'), str(output/'masses.csv')],
                       cwd=args.ems_root, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('parameters', type=Path)
    ap.add_argument('history', type=Path, help='directory containing original rh_surf/rh_f histories')
    ap.add_argument('output', type=Path, help='new output directory (never reused)')
    ap.add_argument('checkpoints', type=Path, nargs='+')
    ap.add_argument('--ems-root', type=Path, help='unchanged EMS.jl tree, required for mass projection')
    ap.add_argument('--executable', type=Path)
    ap.add_argument('--points', type=int, nargs='+', default=[96,192])
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--cap', type=int, default=240, help='seconds per checkpoint/angular resolution; max 1800')
    ap.add_argument('--dry-run', action='store_true', help='validate seeds and print resource/command plan only')
    ap.add_argument('--expected-mass', type=float, default=float('nan'), help='isolated control only')
    ap.add_argument('--expected-area', type=float, default=float('nan'))
    ap.add_argument('--expected-charge', type=float, default=float('nan'))
    ap.add_argument('--sphere-seed', action='store_true', help='seed one puncture-centred sphere at KS2 r_h')
    ap.add_argument('--skip-mass', action='store_true', help='measure surfaces without the horizon-family projection')
    args = ap.parse_args()
    if (not 10 <= args.cap <= 1800 or args.threads < 1 or len(set(args.points)) != len(args.points)
            or any(n < 8 or n > 4096 for n in args.points)):
        ap.error('invalid cap, threads or angular counts')
    if args.sphere_seed and not args.skip_mass:
        ap.error('--sphere-seed requires --skip-mass')
    if not args.skip_mass and args.ems_root is None:
        ap.error('--ems-root is required unless --skip-mass is set')
    args.ems_root = args.ems_root.resolve() if args.ems_root else HERE
    args.output = args.output.resolve()
    p = params(args.parameters)
    if args.sphere_seed:
        ks2 = Path(p['ems_data_path'])
        metadata = {k: v for line in ks2.read_text().splitlines() if line.startswith('# ') and '=' in line
                    for k, v in [line[2:].split('=', 1)]}
        radius = float(metadata['r_h'])
        centre = float(p['star_centre'].split()[0])
        p.update(RH_num_horizons='1', RH_initial_radii=format(radius, '.17g'),
                 RH_initial_centre=format(centre, '.17g'), RH_num_points='96',
                 RH_level='0', RH_time_step_freq='400', RH_chase_speeds='0.125',
                 RH_start_times='0', RH_newton_crit='0')
    binaries = list(HERE.glob('EMSRHCheckpoint2d.*.ex'))
    if not args.executable and len(binaries) != 1:
        ap.error('supply --executable when there is not exactly one build')
    exe = args.executable or binaries[0]
    exe = exe.resolve()
    if args.output.exists():
        raise ValueError('output directory already exists; use a fresh path')
    plans = []
    for cp in args.checkpoints:
        cp = cp.resolve()
        info = checkpoint_info(cp)
        saved, skipped = ({}, {}) if args.sphere_seed else seeds(p, args.history, info['time'])
        plans.append((cp, info, saved, skipped))
    print(json.dumps(dict(checkpoints=[dict(path=str(cp), **info, skipped=skipped) for cp, info, _, skipped in plans],
                          mpi_ranks=1, threads=args.threads, points=args.points, cap_per_case=args.cap,
                          max_batch_seconds=len(plans)*len(args.points)*args.cap,
                          memory_estimate_bytes=2*1024**3+6*max(i['field_bytes'] for _,i,_,_ in plans),
                          executable=str(exe)), indent=2), flush=True)
    if args.dry_run:
        return
    args.output.mkdir(parents=True)
    failed = False
    batch_start = time.monotonic()
    for number, (cp, info, saved, skipped) in enumerate(plans):
        output = args.output/f'{number:04d}-{cp.stem}'
        output.mkdir()
        sources = [cp, args.parameters.resolve(), exe, HERE/'offline.py']
        if args.sphere_seed:
            sources.append(Path(p['ems_data_path']).resolve())
        if not args.skip_mass:
            sources += [HERE/'postprocess.py', HERE/'postprocess.jl',
                        args.ems_root/'artifacts/horizon-family/production/alpha20-positive.hfamily']
        sources += [args.history/f for i in saved for f in (f'rh_surf_{i}.dat', f'rh_f{i}.dat')]
        before = {str(s.resolve()):digest(s) for s in sources}
        runs, statuses = [], []
        for n in args.points:
            directory = output/f'n{n}'
            directory.mkdir()
            q = dict(p)
            for key in ('output_path','hdf5_subpath','data_subpath','pout_subpath'):
                q.pop(key, None)
            q.update(restart_file=str(cp), checkpoint_interval=-1, plot_interval=-1,
                     offline_points=n, offline_seconds=args.cap-5, verbosity=0,
                     ignore_checkpoint_name_mismatch=0)
            if args.sphere_seed:
                q['offline_sphere_seed'] = '1'
            if args.skip_mass:
                q['offline_skip_mass'] = '1'
            (directory/'params.txt').write_text(''.join(f'{k} = {v}\n' for k,v in q.items()))
            for i, pair in saved.items():
                for name, row in zip((f'rh_surf_{i}.dat', f'rh_f{i}.dat'), pair):
                    (directory/name).write_text(' '.join([format(info['time'], '.17g'), *row[1:]])+'\n')
            start = time.monotonic()
            with (directory/'run.log').open('w') as log:
                try:
                    code = subprocess.run([str(exe), 'params.txt'], cwd=directory, stdout=log,
                            stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS=str(args.threads)),
                            timeout=args.cap).returncode
                except subprocess.TimeoutExpired:
                    code = 124
            rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            result = dict(N_theta=n, exit=code, seconds=time.monotonic()-start,
                          peak_child_rss_bytes=int(rss*(1 if sys.platform == 'darwin' else 1024)),
                          stage=f'checkpoint_{number}/N_{n}',
                          elapsed_seconds=time.monotonic()-batch_start,
                          items_done=number*len(args.points)+len(statuses)+1,
                          items_total=len(plans)*len(args.points))
            statuses.append(result)
            (output/'status.json').write_text(json.dumps(statuses, indent=2)+'\n')
            print(json.dumps(result), flush=True)
            if code not in (0,1):
                raise RuntimeError(f'checkpoint find failed: {directory}; exit {code}')
            failed |= code != 0
            runs.append(directory)
        if not args.skip_mass:
            collect(runs, p, saved, args, output)
        if before != {str(s.resolve()):digest(s) for s in sources}:
            raise RuntimeError('an input changed during the offline run')
        (output/'inputs.json').write_text(json.dumps(dict(info=info, skipped=skipped, sha256=before,
                                                        threads=args.threads), indent=2)+'\n')
    raise SystemExit(int(failed))


if __name__ == '__main__':
    main()
