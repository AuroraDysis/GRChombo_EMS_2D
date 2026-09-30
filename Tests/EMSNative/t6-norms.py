#!/usr/bin/env python3
"""Coordinate-volume and legacy cell-weighted norms on frozen reference plots."""
import csv
import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t5', HERE/'t5-analyze.py')
t5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t5)
ROOT = Path('/Users/auroradysis/Workspace/EMS/.data/exp-0019')
TIMES = (1., 2.5, 5., 7.5, 10.)
MASKS = {'horizon': (.63593977642346233, 1.2718795528469247),
         'near_hole': (.1, 2.), 'ring': (2., 4.), 'far': (4., 8.),
         'wake_2p5_3': (2.5, 3.), 'wake_4p5_5p5': (4.5, 5.5)}


def norm(c, w):
    assert len(c) and np.all(w > 0) and np.all(np.isfinite(c))
    return np.sqrt(np.sum(w*c*c)/np.sum(w))


def main():
    # Constant data and split-cell volumes give the same coordinate-volume norm.
    assert norm(np.ones(4)*3, np.array([1., 2., 3., 4.])) == 3
    assert norm(np.array([1., 2.]), np.array([4., 4.])) == norm(
        np.array([1., 2., 2., 2., 2.]), np.array([4., 1., 1., 1., 1.]))
    rows = []
    for scale in t5.SCALES:
        plots = sorted((ROOT/f'ref-{scale}'/'plt').glob('*.hdf5'))
        published = {(round(float(v['time_M']), 8), v['mask'], v['constraint']): v for v in
                     csv.DictReader((ROOT/f'ref-{scale}'/'constraint-rms.csv').open())
                     if v['cell_class'] == 'all'}
        for time in TIMES:
            tt, levels = t5.read(plots[int(round(2*time))])
            assert abs(tt-time) < 1e-8
            assert all(np.count_nonzero(g['a'][t5.COMPONENTS.index('GaussB')]) == 0 for g in levels)
            xy, v, lev, _ = t5.nodes(levels)
            h = np.array([g['h'] for g in levels])[lev]
            legacy_w, volume_w = xy[:, 1], 2*np.pi*xy[:, 1]*h*h
            r = np.hypot(xy[:, 0], xy[:, 1])
            fields = t5.constrained(v)
            for mask, (lo, hi) in MASKS.items():
                use = (r >= lo) & (r <= hi)
                for bulk in (False, True) if mask in ('ring', 'far', 'wake_2p5_3', 'wake_4p5_5p5') else (False,):
                    sel = use & t5.bulk_points(xy, lev) if bulk else use
                    name = 'bulk_'+mask if bulk else mask
                    for ci, constraint in enumerate(t5.CONSTRAINTS):
                        for weighting, w in (('cell_weighted', legacy_w), ('coordinate_volume', volume_w)):
                            rms = norm(fields[ci, sel], w[sel])
                            if weighting == 'cell_weighted' and not bulk:
                                old = published[time, mask, constraint] if (time, mask, constraint) in published else None
                                if old:
                                    assert int(old['cells']) == int(sel.sum())
                                    assert np.isclose(rms, float(old['cylindrical_rms']), rtol=2e-12, atol=1e-30)
                            rows.append(dict(time_M=time, scale=scale, mask=name, constraint=constraint,
                                             weighting=weighting, cells=int(sel.sum()), weight_sum=float(w[sel].sum()),
                                             rms=rms, h_min_M=h[sel].min(), h_max_M=h[sel].max()))
    t5.save('t6-norms.csv', rows)
    lookup = {(v['time_M'], v['mask'], v['constraint'], v['weighting'], v['scale']): v for v in rows}
    orders = []
    for key in sorted({(v['time_M'], v['mask'], v['constraint'], v['weighting']) for v in rows}):
        values = [lookup[*key, s] for s in t5.SCALES]
        rms = np.array([v['rms'] for v in values])
        p = np.log(rms[:-1]/rms[1:])/np.log(1.5)
        orders.append(dict(time_M=key[0], mask=key[1], constraint=key[2], weighting=key[3],
                           **{s+'_rms': rms[i] for i,s in enumerate(t5.SCALES)},
                           order_low_mid=p[0], order_mid_high=p[1],
                           **{s+'_cells': values[i]['cells'] for i,s in enumerate(t5.SCALES)}))
    t5.save('t6-orders.csv', orders)
    print('PASS: constant/split-cell norm checks, legacy reductions; volume norms exported.')


if __name__ == '__main__':
    main()
