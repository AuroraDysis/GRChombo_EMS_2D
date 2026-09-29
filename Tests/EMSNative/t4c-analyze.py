#!/usr/bin/env python3
"""Completed T4 reference runs: native masks and shared physical class support.

Reads evolved fields, run logs and T3's published CSV only; never reads initial
data. The 3/2 cell-centred grids have no coincident centres. Class comparisons
use the intersection of their cell footprints, with native values (no field
interpolation), and the original T2 cylindrical norm convention.
"""
import csv
import hashlib
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t4_analysis', HERE / 't4-analyze.py')
t4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t4)
SCALES = t4.SCALES
CLASSES = t4.CLASSES
MASKS = {'horizon': (t4.t2.RH['R'], 2*t4.t2.RH['R']),
         'near_hole': (.1, 2.), 'ring': (2., 4.), 'far': (4., 8.)}
FIELDS = {'Ham': ('Ham',), 'Mom': ('Mom1', 'Mom2'), 'GaussE': ('GaussE',)}


def save(name, rows):
    with (HERE / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def names(f):
    return [f.attrs[f'component_{i}'].decode()
            for i in range(int(f.attrs['num_components']))]


def boxes(g):
    return np.array([[int(b[k]) for k in ('lo_i', 'lo_j', 'hi_i', 'hi_j')]
                     for b in g['boxes'][:]])


def dense(path):
    """Native T2 classes assembled on rectangular level unions; checked below."""
    levels = []
    with h5py.File(path) as f:
        components = names(f)
        domain = f['level_0'].attrs['prob_domain']
        for lev in range(int(f.attrs['num_levels'])):
            g = f[f'level_{lev}']
            b = boxes(g)
            lo, hi = b[:, :2].min(axis=0), b[:, 2:].max(axis=0)
            nx, ny = hi - lo + 1
            X, Y = np.meshgrid(np.arange(lo[0], hi[0]+1), np.arange(lo[1], hi[1]+1))
            data = np.empty((len(components), ny, nx))
            seen = np.zeros((ny, nx), int)
            seam = np.zeros((ny, nx), bool)
            offsets = g['data:offsets=0'][:]
            assert not any(g['data_attributes'].attrs['outputGhost'])
            for bi, (x0, y0, x1, y1) in enumerate(b):
                ys, xs = slice(y0-lo[1], y1-lo[1]+1), slice(x0-lo[0], x1-lo[0]+1)
                data[:, ys, xs] = g['data:datatype=0'][int(offsets[bi]):int(offsets[bi+1])].reshape(
                    (len(components), y1-y0+1, x1-x0+1))
                seen[ys, xs] += 1
                seam[ys, xs] = ((X[ys, xs]-x0 < 2) | (x1-X[ys, xs] < 2) |
                               (Y[ys, xs]-y0 < 2) | (y1-Y[ys, xs] < 2))
            assert np.all(seen == 1), 'classifier requires a tiled rectangular union'
            valid = np.ones((ny, nx), bool)
            if lev+1 < int(f.attrs['num_levels']):
                for a, c, d, e in boxes(f[f'level_{lev+1}']):
                    valid &= ~((X >= a//2) & (X <= d//2) & (Y >= c//2) & (Y <= e//2))
            ex = ((X-lo[0] < 2) | (hi[0]-X < 2)) & (lev > 0)
            ey = ((Y-lo[1] < 2) | (hi[1]-Y < 2)) & (lev > 0) & (Y >= 2)
            flags = {'all': np.ones((ny, nx), bool), 'convex_corner': ex & ey,
                     'patch_edge': (ex ^ ey) & ~((Y < 2) & ex),
                     'axis_patch_junction': (Y < 2) & ex,
                     'cartoon_axis': (Y < 2) & ~ex,
                     'outer_boundary': ((X < 2) | (X >= (int(domain['hi_i'])+1)*2**lev-2) |
                                        (Y >= (int(domain['hi_j'])+1)*2**lev-2)) & (lev == 0)}
            flags['box_seam'] = seam & ~(ex | ey) & (Y >= 2) & ~flags['outer_boundary']
            flags['interior'] = ~np.logical_or.reduce([flags[k] for k in CLASSES if k not in ('all', 'interior')])
            h = float(g.attrs['dx'])
            rho = np.hypot((X+.5)*h-256., (Y+.5)*h)
            masks = {k: valid & (rho >= a) & (rho <= c) for k, (a, c) in MASKS.items()}
            fields = {k: sum(data[components.index(c)]**2 for c in cs) for k, cs in FIELDS.items()}
            levels.append(dict(h=h, lo=lo, hi=hi, boxes=b, masks=masks, flags=flags,
                               fields=fields, y=(Y+.5)*h))
    return levels


def audit_file(scale, path, kind):
    row = dict(scale=scale, kind=kind, path=str(path), time_M=0.,
               components=0, missing_components='', nonfinite_stored=0,
               nonfinite_valid=0, min_chi=math.inf, min_lapse=math.inf,
               min_chi_stored=math.inf, min_lapse_stored=math.inf,
               chi_stored_at_or_below_floor=0, lapse_stored_at_or_below_floor=0,
               chi_below_floor=0, chi_at_floor=0, lapse_below_floor=0,
               lapse_at_floor=0, chi_below_constraint_regularizer=0,
               max_abs_GaussB='')
    with h5py.File(path) as f:
        cs = names(f)
        expected = ('chi', 'lapse', 'phi', 'Theta', *[c for v in FIELDS.values() for c in v], 'GaussB', 'Qscalar') if kind == 'plot' else (
            'chi', 'h11', 'h12', 'h22', 'hww', 'K', 'A11', 'A12', 'A22', 'Aww',
            'Theta', 'Gamma1', 'Gamma2', 'lapse', 'shift1', 'shift2', 'B1', 'B2',
            'phi', 'Pi', 'Lambda', 'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Xi')
        row.update(time_M=float(f['level_0'].attrs['time']), components=len(cs),
                   missing_components=';'.join(c for c in expected if c not in cs))
        assert not row['missing_components']
        assert len(cs) == len(expected)
        if kind == 'plot':
            row['max_abs_GaussB'] = 0.
        for lev in range(int(f.attrs['num_levels'])):
            g = f[f'level_{lev}']
            ghost = tuple(g['data_attributes'].attrs['outputGhost'])
            offsets = g['data:offsets=0'][:]
            assert g['data:datatype=0'].dtype == np.dtype('float64')
            for bi, (x0, y0, x1, y1) in enumerate(boxes(g)):
                nx, ny = x1-x0+1, y1-y0+1
                v = g['data:datatype=0'][int(offsets[bi]):int(offsets[bi+1])].reshape(
                    (len(cs), ny+2*ghost[1], nx+2*ghost[0]))
                row['nonfinite_stored'] += int(np.count_nonzero(~np.isfinite(v)))
                for c in ('chi', 'lapse'):
                    a = v[cs.index(c)]
                    row['min_'+c+'_stored'] = min(row['min_'+c+'_stored'], float(a.min()))
                    row[c+'_stored_at_or_below_floor'] += int(np.count_nonzero(a <= 1e-12))
                v = v[:, ghost[1]:ghost[1]+ny, ghost[0]:ghost[0]+nx]
                row['nonfinite_valid'] += int(np.count_nonzero(~np.isfinite(v)))
                for c in ('chi', 'lapse'):
                    a = v[cs.index(c)]
                    row['min_'+c] = min(row['min_'+c], float(a.min()))
                    row[c+'_below_floor'] += int(np.count_nonzero(a < 1e-12))
                    row[c+'_at_floor'] += int(np.count_nonzero(a == 1e-12))
                row['chi_below_constraint_regularizer'] += int(np.count_nonzero(v[cs.index('chi')] < 1e-6))
                if kind == 'plot':
                    row['max_abs_GaussB'] = max(row['max_abs_GaussB'], float(np.abs(v[cs.index('GaussB')]).max()))
    return row


def main():
    root, legacy_root, legacy_csv = map(Path, sys.argv[1:4])
    published = {(float(r['time_M']), r['mask'], r['cell_class'], r['constraint']): r
                 for r in csv.DictReader(legacy_csv.open()) if r['restriction'] == 'off'}
    audits, logs, native, shared_rows, support = [], [], {}, [], []
    native_rows = []
    for scale in SCALES:
        run = root / 'evolution' / f'ref-{scale}'
        assert (run / 'done.exit').read_text().strip() == '0'
        log_paths = [run/'run.log'] + sorted(run.glob('**/pout*'))
        log = '\n'.join(p.read_text(errors='replace') for p in log_paths if p.is_file())
        params = (HERE / 'params' / f't4-ref-{scale}-2p5M.txt').read_text()
        parameter = dict(re.findall(r'^\s*(\w+)\s*=\s*([^#\n]+)', params, re.M))
        assert float(parameter['sigma']) == 1.
        assert parameter['amr_transfer'].strip() == 'point'
        assert float(parameter['min_chi']) == float(parameter['min_lapse']) == 1e-12
        defaults = sorted(set(re.findall(r'Parameter: (\w+) not found', log)))
        bad = [s for s in log.splitlines() if re.search(
            r'\b(?:nan|inf|infinity|nonfinite|fatal|error|warning)\b|unknown variable|variable.*not found', s, re.I)]
        first_advance = next(i for i, s in enumerate(log.splitlines()) if 'GRAMRLevel::advance' in s)
        reads = [i for i, s in enumerate(log.splitlines()) if 'Read EMSTRUMPET' in s]
        rh_files = sorted(str(p) for p in run.rglob('*') if p.is_file() and re.match(r'rh_(?:surf|f)\d*', p.name))
        logs.append(dict(scale=scale, exit=0, log_files=';'.join(map(str, log_paths)),
                         finished='GRChombo finished.' in log, amr_transfer='point', sigma=parameter['sigma'].strip(),
                         min_chi_floor=parameter['min_chi'].strip(), min_lapse_floor=parameter['min_lapse'].strip(),
                         nan_check=parameter['nan_check'].strip(), defaulted_parameters=';'.join(defaults),
                         missing_required_variables=0, error_nonfinite_lines=';'.join(bad),
                         nonfinite_dat_tokens=sum(len(re.findall(r'\b(?:nan|inf|infinity)\b', p.read_text(), re.I))
                                                  for p in run.glob('*.dat')),
                         initial_reader_messages=len(reads), reader_messages_after_first_advance=sum(i > first_advance for i in reads),
                         RH_activate=parameter['RH_activate'].strip(), RH_outputs=';'.join(rh_files),
                         horizon_area='', horizon_charge='', horizon_status='DISABLED_NO_OUTPUT',
                         wall_seconds=float((run/'wall_seconds').read_text()),
                         conservative_peak_rss_bytes=int((run/'peak_rss_bytes').read_text()),
                         legacy_csv=str(legacy_csv), legacy_csv_sha256=hashlib.sha256(legacy_csv.read_bytes()).hexdigest()))
        assert not bad and logs[-1]['finished'] and not logs[-1]['reader_messages_after_first_advance']
        assert not logs[-1]['nonfinite_dat_tokens']
        paths = sorted((run/'plt').glob('*.hdf5'))
        assert len(paths) == 6
        audits.extend(audit_file(scale, p, 'plot') for p in paths)
        audits.append(audit_file(scale, sorted((run/'chk').glob('*.hdf5'))[-1], 'final_checkpoint'))
    for time, steps in ((1., (8, 12, 18)), (2.5, (20, 30, 45))):
        grids = {}
        for scheme in ('legacy', 'point'):
            grids[scheme] = []
            for scale, step in zip(SCALES, steps):
                run = root/'evolution'/f'ref-{scale}' if scheme == 'point' else legacy_root/f't3-ref-match-{scale}-off'
                # Original T3 low/mid were attached and have no exit markers.
                if scheme == 'point' or scale == 'high':
                    assert (run/'done.exit').read_text().strip() == '0'
                path = run/'plt'/f'EMS_Plot_{step:06d}.2d.hdf5'
                with h5py.File(path) as f:
                    assert abs(float(f['level_0'].attrs['time'])-time) < 1e-8
                grids[scheme].append(dense(path))
                measured = t4.t2.measure(run.name, path, 'R', 256., CLASSES)
                for _, _, mask, cl, field, n, rms, maximum in measured:
                    native[scheme, scale, time, mask, cl, field] = (n, rms)
                    native_rows.append(dict(scheme=scheme, scale=scale, time_M=time,
                                            mask=mask, cell_class=cl, constraint=field,
                                            cells=n, cylindrical_rms=rms, maximum=maximum, path=str(path)))
                    # Independent assembly must reproduce original class and norm.
                    ss = w = count = 0
                    for g in grids[scheme][-1]:
                        selected = g['masks'][mask] & g['flags'][cl]
                        ss += float(np.sum(g['y'][selected]*g['fields'][field][selected]))
                        w += float(np.sum(g['y'][selected]))
                        count += int(np.count_nonzero(selected))
                    assert count == n
                    assert math.isclose(math.sqrt(ss/w), rms, rel_tol=2e-13, abs_tol=1e-30)
                    if scheme == 'legacy' and (time, mask, cl, field) in published:
                        old = published[time, mask, cl, field]
                        assert int(old[scale+'_cells']) == n
                        assert math.isclose(float(old[scale+'_rms']), rms, rel_tol=2e-13)
        accum = defaultdict(lambda: [0., 0., 0])
        for lev in range(len(grids['point'][0])):
            gs = [grids['point'][i][lev] for i in range(3)]
            for i, g in enumerate(gs):
                old = grids['legacy'][i][lev]
                assert np.array_equal(old['boxes'], g['boxes']) and old['h'] == g['h']
                assert np.allclose(g['lo']*g['h'], gs[0]['lo']*gs[0]['h'], rtol=0, atol=1e-12)
                assert np.allclose((g['hi']+1)*g['h'], (gs[0]['hi']+1)*gs[0]['h'], rtol=0, atol=1e-12)
            for mask in MASKS:
                for cl in CLASSES[1:]:
                    chosen = [g['masks'][mask] & g['flags'][cl] for g in gs]
                    if not all(a.any() for a in chosen):
                        continue
                    repeated = [a.repeat(r, axis=0).repeat(r, axis=1) for a, r in zip(chosen, (9, 6, 4))]
                    assert repeated[0].shape == repeated[1].shape == repeated[2].shape
                    common = repeated[0] & repeated[1] & repeated[2]
                    if not common.any():
                        continue
                    delta = gs[0]['h']/9
                    y = (gs[0]['lo'][1]*9+np.arange(common.shape[0])+.5)*delta
                    weights = common*y[:, None]
                    census = []
                    for i, (g, r, scale) in enumerate(zip(gs, (9, 6, 4), SCALES)):
                        ny, nx = chosen[i].shape
                        native_weight = weights.reshape(ny, r, nx, r).sum(axis=(1, 3))
                        count = int(np.count_nonzero(native_weight))
                        census.append(count)
                        assert np.all((native_weight > 0) <= chosen[i])
                        # Full-footprint reduction preserves the T2 per-cell y weights.
                        overlay_y = np.broadcast_to(y[:, None], common.shape)
                        assert np.allclose(overlay_y.reshape(ny, r, nx, r).sum(axis=(1, 3)),
                                           g['y']*r*r, rtol=2e-14, atol=1e-14)
                        for scheme in grids:
                            source = grids[scheme][i][lev]
                            for field in FIELDS:
                                a = accum[scheme, scale, mask, cl, field]
                                a[0] += float(np.sum(native_weight*source['fields'][field]))
                                a[1] += float(native_weight.sum())
                                a[2] += count
                    support.append(dict(time_M=time, mask=mask, cell_class=cl, level=lev,
                                        overlay_dx_M=delta, overlay_cells=int(common.sum()),
                                        cylindrical_physical_measure_M3=float(weights.sum()*delta*delta),
                                        low_cells=census[0], mid_cells=census[1], high_cells=census[2]))
        for (scheme, scale, mask, cl, field), (ss, w, n) in sorted(accum.items()):
            shared_rows.append(dict(scheme=scheme, scale=scale, time_M=time, mask=mask,
                                    cell_class=cl, constraint=field, cells=n,
                                    cylindrical_rms=math.sqrt(ss/w), overlay_weight=w))
    lookup = {(r['scheme'], r['scale'], r['time_M'], r['mask'], r['cell_class'], r['constraint']):
              (r['cells'], r['cylindrical_rms']) for r in shared_rows}
    orders = []
    for time in (1., 2.5):
        for mask in MASKS:
            for cl in CLASSES:
                for field in FIELDS:
                    row = dict(time_M=time, mask=mask, cell_class=cl, constraint=field,
                               sampling='native_T2' if cl == 'all' else 'common_footprints_same_level')
                    old = published.get((time, mask, cl, field))
                    for k in ('order_low_mid', 'order_mid_high', 'low_cells', 'mid_cells', 'high_cells'):
                        row['legacy_published_'+k] = old[k] if old else ''
                    for scheme in ('legacy', 'point'):
                        table = native if cl == 'all' else lookup
                        values = [table.get((scheme, s, time, mask, cl, field)) for s in SCALES]
                        for s, v in zip(SCALES, values):
                            row[scheme+'_'+s+'_rms'] = v[1] if v else ''
                            row[scheme+'_'+s+'_cells'] = v[0] if v else 0
                        p = [math.log(values[i][1]/values[i+1][1])/math.log(1.5)
                             if values[i] and values[i+1] and min(values[i][1], values[i+1][1]) > 0 else ''
                             for i in range(2)]
                        row[scheme+'_order_low_mid'], row[scheme+'_order_mid_high'] = p
                        row[scheme+'_status'] = ('ABSENT_COMMON' if not all(values) else
                                                'ZERO' if not any(v[1] for v in values) else
                                                'PASS' if all(isinstance(q, float) and q >= 3 for q in p) else 'FAIL')
                    orders.append(row)
    save('t4c-run-audit.csv', logs)
    save('t4c-field-checks.csv', audits)
    save('t4c-native-norms.csv', native_rows)
    save('t4c-common-norms.csv', shared_rows)
    save('t4c-common-support.csv', support)
    save('t4c-orders.csv', orders)
    print('Verified native T2 norms/counts and complete supplied T3 CSV; matched legacy/point layouts.')
    print('All output files finite:', all(not r['nonfinite_stored'] for r in audits))
    print('Orders:', len(orders), 'point pass:', sum(r['point_status'] == 'PASS' for r in orders),
          'absent common:', sum(r['point_status'] == 'ABSENT_COMMON' for r in orders))
    for r in orders:
        if r['cell_class'] == 'all':
            print(r['time_M'], r['mask'], r['constraint'],
                  r['legacy_order_low_mid'], r['legacy_order_mid_high'],
                  r['point_order_low_mid'], r['point_order_mid_high'])


if __name__ == '__main__':
    main()
