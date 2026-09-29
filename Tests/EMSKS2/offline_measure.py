"""Measure frozen reference-stationary KS2 checkpoints; run --help for usage."""
import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


def rows(path):
    with path.open(newline='') as file:
        return list(csv.DictReader(file))


def sinc_correction(n):
    half = math.pi / (2*n)
    return math.sin(half) / half


def collect(work, output):
    records = []
    for checkpoint in sorted(work.iterdir()):
        if not checkpoint.is_dir():
            continue
        by_n, exterior_by_key = {}, {}
        for directory in sorted(checkpoint.glob('n*')):
            surfaces = rows(directory/'surfaces.csv')
            final = next(r for r in surfaces if int(r['stage']) == 2)
            n = int(final['N_theta'])
            correction = sinc_correction(n)
            row = dict(final, kind='horizon', A_corrected=float(final['A'])*correction,
                       Q_corrected=float(final['Q'])*correction,
                       area_delta_corrected=float(final['area_delta'])*correction,
                       charge_delta_corrected=float(final['charge_delta'])*correction)
            row['R_corrected'] = math.sqrt(row['A_corrected']/(4*math.pi))
            by_n[n] = row
            records.append(row)
            for exterior in rows(directory/'exterior.csv'):
                exterior.update(kind='exterior_'+exterior['kind'],
                                Q_corrected=float(exterior['Q_S'])*correction)
                for key in ('m', 'phi', 'lapse', 'K'):
                    exterior[key+'_minus_static'] = (float(exterior[key+'_mean'])
                                                     - float(exterior[key+'_static']))
                exterior['Q_minus_static'] = exterior['Q_corrected']-float(exterior['Q_static'])
                if exterior['kind'] == 'exterior_areal':
                    exterior['areal_root_rel_tol'] = 1e-12
                exterior_by_key[(n, exterior['kind'], exterior['target_over_rh'])] = exterior
                records.append(exterior)
        for n, row in by_n.items():
            neighbor = by_n.get(2*n) or by_n.get(n//2)
            if neighbor:
                factor = 4/3 if 2*n in by_n else 1/3
                for key in ('A_corrected', 'Q_corrected'):
                    row['angular_delta_'+key] = factor*abs(row[key]-neighbor[key])
            else:
                row['angular_delta_A_corrected'] = math.nan
                row['angular_delta_Q_corrected'] = math.nan
        for (n, kind, radius), row in exterior_by_key.items():
            neighbor = exterior_by_key.get((2*n, kind, radius)) or exterior_by_key.get((n//2, kind, radius))
            row['angular_delta_Q_corrected'] = (abs(row['Q_corrected']-neighbor['Q_corrected'])
                                                  * (4/3 if (2*n, kind, radius) in exterior_by_key else 1/3)
                                                  if neighbor else math.nan)
    baseline = {(r['N_theta'], r['target_over_rh']): float(r['R_mean'])
                for r in records if r['kind'] == 'exterior_coordinate' and float(r['time']) == 0.}
    for row in records:
        if row['kind'] == 'exterior_coordinate':
            key = row['N_theta'], row['target_over_rh']
            if key in baseline:
                row['coordinate_R_drift_rel'] = float(row['R_mean'])/baseline[key]-1
    fields = list(dict.fromkeys(k for row in records for k in row))
    with output.open('w', newline='') as file:
        writer = csv.DictWriter(file, fields)
        writer.writeheader()
        writer.writerows(records)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path, help='run directory with params.txt and chk/')
    parser.add_argument('output', type=Path, help='new CSV path; work stays beside it')
    parser.add_argument('--checkpoints', type=Path, nargs='*', help='default: all chk/*.hdf5')
    parser.add_argument('--points', type=int, nargs='+', default=[48, 96, 192])
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--cap', type=int, default=240, help='seconds per checkpoint and angular count')
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    run, output = args.run.resolve(), args.output.resolve()
    checkpoints = args.checkpoints or sorted((run/'chk').glob('*.hdf5'))
    if not checkpoints:
        parser.error('no checkpoints found')
    if not (run/'params.txt').is_file() or output.exists() or output.with_suffix('.work').exists():
        parser.error('params.txt missing or output already exists')
    command = [sys.executable, str(Path(__file__).resolve().parents[1]/'EMSRHFinder'/'offline.py'),
               str(run/'params.txt'), str(run), str(output.with_suffix('.work')),
               *(str(p.resolve()) for p in checkpoints), '--points', *(str(n) for n in args.points),
               '--threads', str(args.threads), '--cap', str(args.cap), '--sphere-seed', '--skip-mass']
    if args.executable:
        command += ['--executable', str(args.executable.resolve())]
    if args.dry_run:
        command.append('--dry-run')
    result = subprocess.run(command, check=False)
    if not args.dry_run:
        collect(output.with_suffix('.work'), output)
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
