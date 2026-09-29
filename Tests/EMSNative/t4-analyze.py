#!/usr/bin/env python3
"""T4: reuse registered T2/T3 masks/classes; compare HDF5 value bits."""
import csv
import importlib.util
import math
from pathlib import Path
import sys

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t2', HERE / 't2-localize.py')
t2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t2)
CLASSES = ('all', 'patch_edge', 'convex_corner', 'axis_patch_junction',
           'cartoon_axis', 'box_seam', 'outer_boundary', 'interior')
SCALES = ('low', 'mid', 'high')


def save(name, header, rows):
    with (HERE / name).open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def value_equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.dtype.names:
        return all(value_equal(a[k], b[k]) for k in a.dtype.names)
    if a.dtype.kind in 'OU':
        return bool(np.array_equal(a, b))
    return a.tobytes() == b.tobytes()


def compare(a, b):
    datasets, attributes, failures = 0, 0, []
    with h5py.File(a) as f, h5py.File(b) as g:
        fn, gn = [], []
        f.visit(fn.append)
        g.visit(gn.append)
        if fn != gn:
            failures.append('object_paths')
        for name in [''] + fn:
            if name not in g:
                continue
            x, y = f[name or '/'], g[name or '/']
            if set(x.attrs) != set(y.attrs):
                failures.append(name + ':attribute_names')
            for key in x.attrs:
                attributes += 1
                if key not in y.attrs or not value_equal(x.attrs[key], y.attrs[key]):
                    failures.append(name + ':' + key)
            if isinstance(x, h5py.Dataset):
                datasets += 1
                if not value_equal(x[:], y[:]):
                    failures.append(name)
    return datasets, attributes, failures


def identities(root):
    rows = []
    for member in ('E', 'ref'):
        for test, a, b, subpath in (
                ('default_t0', 'baseline', 'off', 'plt/EMS_Plot_000000.2d.hdf5'),
                ('default_registered_t0', 'registered-t0-baseline', 'registered-t0-off', 'plt/EMS_Plot_000000.2d.hdf5'),
                ('default_two_steps', 'baseline', 'off', 'chk/EMS_000002.2d.hdf5'),
                ('default_regrid_two_steps', 'regrid-baseline', 'regrid-off', 'chk/EMS_000002.2d.hdf5'),
                ('point_restart_four_steps', 'point-full', 'point-restart', 'chk/EMS_000004.2d.hdf5')):
            left, right = (root / 'control' / f'{member}-{x}' / subpath for x in (a, b))
            if not left.exists() or not right.exists():
                rows.append((member, test, '', '', 'PENDING', str(left), str(right), ''))
                continue
            n, m, failure = compare(left, right)
            rows.append((member, test, n, m, 'FAIL' if failure else 'PASS',
                         str(left), str(right), ';'.join(failure)))
    save('t4-bit-identity.csv', ('member', 'test', 'datasets', 'attributes', 'status',
                               'left', 'right', 'differences'), rows)
    print('identity:', [(x[0], x[1], x[4]) for x in rows])


def constraints(root, evolved=False):
    phase = 'evolved' if evolved else 't0'
    rows, layouts, checks = [], [], []
    for member in (('ref',) if evolved else ('E', 'ref')):
        for scale, final_step in zip(SCALES, (20, 30, 45)):
            directory = root / ('evolution' if evolved else 't0') / f'{member}-{scale}'
            marker = directory / 'done.exit'
            if not marker.exists() or marker.read_text().strip() != '0':
                continue
            steps = (0, final_step) if evolved else (0,)
            for step in steps:
                path = directory / 'plt' / f'EMS_Plot_{step:06d}.2d.hdf5'
                with h5py.File(path) as f:
                    time = float(f['level_0'].attrs['time'])
                    names = [f.attrs[f'component_{i}'].decode()
                             for i in range(int(f.attrs['num_components']))]
                    nonfinite, bmax = 0, 0.
                    for level in range(int(f.attrs['num_levels'])):
                        g = f[f'level_{level}']
                        offsets = g['data:offsets=0'][:]
                        for bi, box in enumerate(g['boxes'][:]):
                            cells = ((int(box['hi_i'])-int(box['lo_i'])+1) *
                                     (int(box['hi_j'])-int(box['lo_j'])+1))
                            data = g['data:datatype=0'][int(offsets[bi]):int(offsets[bi+1])]
                            nonfinite += int(np.count_nonzero(~np.isfinite(data)))
                            b = data[names.index('GaussB')*cells:(names.index('GaussB')+1)*cells]
                            bmax = max(bmax, float(np.max(np.abs(b))))
                        if step == 0:
                            boxes = g['boxes'][:]
                            h = float(g.attrs['dx'])
                            center = 224. if member == 'E' else 256.
                            layouts.append((phase, member, scale, level, h, len(boxes),
                                min(int(b['lo_i']) for b in boxes)*h-center,
                                (max(int(b['hi_i']) for b in boxes)+1)*h-center,
                                min(int(b['lo_j']) for b in boxes)*h,
                                (max(int(b['hi_j']) for b in boxes)+1)*h))
                    checks.append((member, scale, time, nonfinite, bmax))
                for _, _, mask, cl, field, n, rms, maximum in t2.measure(
                        directory.name, path, 'E' if member == 'E' else 'R',
                        224. if member == 'E' else 256., CLASSES,
                        {'cavity_core': (.055, .075), 'outer_ring_core': (.85, 1.),
                         'far_core': (6.5, 8.)} if member == 'E' else None):
                    rows.append((member, scale, time, mask, cl, field, n, rms, maximum))
    save(f't4-{phase}-norms.csv', ('member', 'scale', 'time_M', 'mask', 'cell_class',
                                 'constraint', 'cells', 'cylindrical_rms', 'maximum'), rows)
    save(f't4-{phase}-layout.csv', ('phase', 'member', 'scale', 'level', 'dx', 'boxes',
                                  'xlo_M', 'xhi_M', 'ylo_M', 'yhi_M'), layouts)
    save(f't4-{phase}-field-checks.csv', ('member', 'scale', 'time_M',
                                       'nonfinite_saved_plot_values', 'max_abs_GaussB'), checks)
    lookup = {(m, s, round(t, 8), mask, cl, field): (n, rms)
              for m, s, t, mask, cl, field, n, rms, _ in rows}
    orders = []
    for m, t, mask, cl, field in sorted({(m, round(t, 8), mask, cl, field)
                                        for m, _, t, mask, cl, field, _, _, _ in rows}):
        values = [lookup.get((m, s, t, mask, cl, field)) for s in SCALES]
        rms = [v[1] if v else '' for v in values]
        order = [math.log(rms[i]/rms[i+1])/math.log(1.5)
                 if values[i] and values[i+1] and rms[i] > 0 and rms[i+1] > 0
                 else '' for i in range(2)]
        complete = all((root / ('evolution' if evolved else 't0') / f'{m}-{scale}' /
                        'done.exit').exists() and
                       (root / ('evolution' if evolved else 't0') / f'{m}-{scale}' /
                        'done.exit').read_text().strip() == '0' for scale in SCALES)
        status = ('ABSENT_ON_GRID' if not all(values) and complete else
                  'PENDING' if not all(values) else
                  'ZERO' if not any(rms) else
                  'ROUNDOFF' if max(rms) < 1e-13 else
                  'PASS' if all(p >= 3.5 for p in order) else 'FAIL')
        orders.append((m, t, mask, cl, field, *rms, *order, status,
                       *(v[0] if v else '' for v in values)))
    save(f't4-{phase}-orders.csv', ('member', 'time_M', 'mask', 'cell_class', 'constraint',
                                  'low_rms', 'mid_rms', 'high_rms', 'order_low_mid',
                                  'order_mid_high', 'status', 'low_cells', 'mid_cells',
                                  'high_cells'), orders)
    print(phase, 'orders:', len(orders), 'failures:', sum(r[10] == 'FAIL' for r in orders))


def wave_line(path, component):
    x, u = [], []
    with h5py.File(path) as f:
        for level in range(int(f.attrs['num_levels'])):
            g = f[f'level_{level}']
            fine = f[f'level_{level+1}/boxes'][:] if level + 1 < int(f.attrs['num_levels']) else []
            h = float(g.attrs['dx'])
            offsets = g['data:offsets=0'][:]
            for bi, box in enumerate(g['boxes'][:]):
                lo, hi = int(box['lo_i']), int(box['hi_i'])
                yl, yh = int(box['lo_j']), int(box['hi_j'])
                if yl != 0:
                    continue
                xx = np.arange(lo, hi+1)
                valid = np.ones(len(xx), bool)
                for b in fine:
                    if int(b['lo_j']) == 0:
                        valid &= ~((xx >= int(b['lo_i'])//2) & (xx <= int(b['hi_i'])//2))
                data = g['data:datatype=0'][int(offsets[bi]):int(offsets[bi+1])]
                values = data.reshape((2, yh-yl+1, hi-lo+1))[component, 0]
                x.extend(((xx[valid]+.5)*h).tolist())
                u.extend(values[valid].tolist())
    indices = np.argsort(x)
    return np.asarray(x)[indices], np.asarray(u)[indices]


def sample(x, u, targets):
    # Independent six-point interpolation on the composite NONUNIFORM nodes.
    result = []
    for target in targets:
        first = max(0, min(len(x)-6, int(np.searchsorted(x, target))-3))
        nodes, vals = x[first:first+6], u[first:first+6]
        value = 0.
        for i in range(6):
            w = 1.
            for j in range(6):
                if i != j:
                    w *= (target-nodes[j])/(nodes[i]-nodes[j])
            value += w*vals[i]
        result.append(value)
    return np.asarray(result)


def waves(root):
    rows = []
    targets = np.linspace(1., 4., 481)
    for mode in ('legacy', 'point'):
        for component, field in enumerate(('u', 'v')):
            values = []
            for scale, step in zip(SCALES, (48, 96, 192)):
                directory = root / 'wave' / f'{scale}-{mode}'
                marker = directory / 'done.exit'
                if not marker.exists() or marker.read_text().strip() != '0':
                    break
                path = directory / 'plt' / f'EMS_Plot_{step:06d}.2d.hdf5'
                values.append(sample(*wave_line(path, component), targets))
            if len(values) != 3:
                rows.append((mode, field, .75, '', '', '', '', '', '', 'PENDING'))
                continue
            differences = [values[i]-values[i+1] for i in range(2)]
            l2 = [float(np.sqrt(np.mean(d*d))) for d in differences]
            linf = [float(np.max(np.abs(d))) for d in differences]
            p2, pi = math.log2(l2[0]/l2[1]), math.log2(linf[0]/linf[1])
            rows.append((mode, field, .75, *l2, p2, *linf, pi,
                         'PASS' if p2 >= 3.5 else 'FAIL'))
    save('t4-wave-orders.csv', ('scheme', 'field', 'time', 'low_mid_rms', 'mid_high_rms',
                               'self_order_rms', 'low_mid_max', 'mid_high_max',
                               'self_order_max', 'status'), rows)
    print('wave:', [(x[0], x[1], x[5], x[-1]) for x in rows])


if __name__ == '__main__':
    root = Path(sys.argv[1])
    identities(root)
    constraints(root)
    constraints(root, evolved=True)
    waves(root)
