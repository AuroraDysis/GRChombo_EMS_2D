#!/usr/bin/env python3
"""T5: read exp-0019 only; compare current resolutions at common coordinates."""
import argparse
import csv
import hashlib
import time
import resource
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
SCALES = ('low', 'mid', 'high')
VARIABLES = ('chi', 'lapse', 'phi', 'Theta', 'Qscalar')
CONSTRAINTS = ('Ham', 'Mom', 'GaussE')
COMPONENTS = ('chi', 'lapse', 'phi', 'Theta', 'Ham', 'Mom1', 'Mom2', 'GaussE', 'GaussB', 'Qscalar')
EVOLVED_COMPONENTS = ('chi', 'h11', 'h12', 'h22', 'hww', 'K', 'A11', 'A12', 'A22', 'Aww',
                      'Theta', 'Gamma1', 'Gamma2', 'lapse', 'shift1', 'shift2', 'B1', 'B2',
                      'phi', 'Pi', 'Lambda', 'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Xi')
EDGES = np.arange(0., 20.0001, .125)
RADII = (EDGES[:-1]+EDGES[1:])/2
LINE_RADII = np.arange(.00390625, 20., .0078125)
OUTER_RADII = np.arange(200.5, 256., 1.)
GLOBAL_EDGES = np.arange(0., 364.0001, 2.)
GLOBAL_RADII = (GLOBAL_EDGES[:-1]+GLOBAL_EDGES[1:])/2
DIRECTIONS = ('all', 'angle_0_30', 'angle_30_60', 'angle_60_90',
              'angle_90_120', 'angle_120_150', 'angle_150_180',
              'axis_plus', 'axis_minus', 'equator', 'diagonal')
Q4 = 1.5**4


def save(name, rows):
    with (HERE/name).open('w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def read(path):
    levels = []
    with h5py.File(path, 'r') as f:
        cs = tuple(f.attrs[f'component_{i}'].decode() for i in range(int(f.attrs['num_components'])))
        assert cs == COMPONENTS
        for lev in range(int(f.attrs['num_levels'])):
            g = f[f'level_{lev}']
            b = np.array([[int(v[k]) for k in ('lo_i', 'lo_j', 'hi_i', 'hi_j')] for v in g['boxes'][:]])
            lo, hi = b[:, :2].min(axis=0), b[:, 2:].max(axis=0)
            nx, ny = hi-lo+1
            a = np.empty((len(cs), ny, nx))
            seen = np.zeros((ny, nx), int)
            off = g['data:offsets=0'][:]
            assert not any(g['data_attributes'].attrs['outputGhost'])
            assert g['data:datatype=0'].dtype == np.dtype('float64')
            for j, (x0, y0, x1, y1) in enumerate(b):
                ys, xs = slice(y0-lo[1], y1-lo[1]+1), slice(x0-lo[0], x1-lo[0]+1)
                a[:, ys, xs] = g['data:datatype=0'][int(off[j]):int(off[j+1])].reshape(len(cs), y1-y0+1, x1-x0+1)
                seen[ys, xs] += 1
            assert np.all(seen == 1) and np.all(np.isfinite(a))
            h = float(g.attrs['dx'])
            levels.append(dict(h=h, lo=lo, hi=hi, a=a, boxes=b,
                               bounds=(lo[0]*h-256., (hi[0]+1)*h-256., lo[1]*h, (hi[1]+1)*h)))
        t = float(f['level_0'].attrs['time'])
    return t, levels


def nodes(levels, maximum=20.):
    """Composite native points and an interior subset excluding two-cell features."""
    coords, values, levs, interiors = [], [], [], []
    for lev, g in enumerate(levels):
        h, lo, hi = g['h'], g['lo'], g['hi']
        X, Y = np.meshgrid(np.arange(lo[0], hi[0]+1), np.arange(lo[1], hi[1]+1))
        valid = np.ones(X.shape, bool)
        if lev+1 < len(levels):
            for x0, y0, x1, y1 in levels[lev+1]['boxes']:
                valid &= ~((X >= x0//2) & (X <= x1//2) & (Y >= y0//2) & (Y <= y1//2))
        x, y = (X+.5)*h-256., (Y+.5)*h
        valid &= np.hypot(x, y) <= maximum
        seam = np.zeros(X.shape, bool)
        for x0, y0, x1, y1 in g['boxes']:
            in_box = (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)
            seam |= in_box & ((X-x0 < 2) | (x1-X < 2) | (Y-y0 < 2) | (y1-Y < 2))
        interior = ~seam & (Y >= 2)
        coords.append(np.stack((x[valid], y[valid]), axis=1))
        values.append(g['a'][:, valid])
        levs.append(np.full(valid.sum(), lev))
        interiors.append(interior[valid])
    return np.concatenate(coords), np.concatenate(values, axis=1), np.concatenate(levs), np.concatenate(interiors)


def weights(q, first, n):
    out = np.ones((len(q), n))
    for i in range(n):
        for j in range(n):
            if i != j:
                out[:, i] *= (q-first-j)/(i-j)
    return out


def sample(levels, xy, n):
    """Tensor point interpolation; parity at y=0; one-sided at patch boundaries."""
    chosen = np.full(len(xy), -1, int)
    for lev, g in enumerate(levels):
        a, b, c, d = g['bounds']
        chosen[(xy[:, 0] >= a) & (xy[:, 0] <= b) & (xy[:, 1] >= c) & (xy[:, 1] <= d)] = lev
    assert np.all(chosen >= 0)
    out = np.empty((len(COMPONENTS), len(xy)))
    for lev, g in enumerate(levels):
        mask = chosen == lev
        if not mask.any():
            continue
        h, lo, a = g['h'], g['lo'], g['a']
        qx = (xy[mask, 0]+256.)/h-.5-lo[0]
        qy = xy[mask, 1]/h-.5-lo[1]
        startx = np.clip(np.floor(qx).astype(int)-(n//2-1), 0, a.shape[2]-n)
        # y=0 parity ghosts belong to current arrays; all unions meet the axis.
        assert lo[1] == 0
        starty = np.minimum(np.floor(qy).astype(int)-(n//2-1), a.shape[1]-n)
        ix = startx[:, None]+np.arange(n)
        rawy = starty[:, None]+np.arange(n)
        iy = np.where(rawy < 0, -rawy-1, rawy)
        wx, wy = weights(qx, startx, n), weights(qy, starty, n)
        anchor = a[:, iy[:, n//2], ix[:, n//2]]
        result = anchor.copy()
        for j in range(n):
            parity = np.where(rawy[:, j] < 0, -1., 1.)
            for i in range(n):
                v = a[:, iy[:, j], ix[:, i]].copy()
                v[COMPONENTS.index('Mom2')] *= parity
                result += (v-anchor)*(wy[:, j]*wx[:, i])[None, :]
        out[:, mask] = result
    return out, chosen


def bins(xy, edges=EDGES):
    r = np.hypot(xy[:, 0], xy[:, 1])
    radial = np.minimum(np.searchsorted(edges, r, side='right')-1, len(edges)-2)
    angular = np.minimum((np.arctan2(xy[:, 1], xy[:, 0])*6/np.pi).astype(int), 5)
    return radial, angular


def rms_bins(a, xy, subset=None, edges=EDGES):
    rb, ab = bins(xy, edges)
    nb = len(edges)-1
    subset = np.ones(len(xy), bool) if subset is None else subset
    result = np.full((7, a.shape[0], nb), np.nan)
    counts = np.zeros((7, nb), int)
    for d in range(7):
        use = subset & ((ab == d-1) if d else True)
        w = np.bincount(rb[use], weights=xy[use, 1], minlength=nb)
        counts[d] = np.bincount(rb[use], minlength=nb)
        for k in range(a.shape[0]):
            ss = np.bincount(rb[use], weights=xy[use, 1]*a[k, use]**2, minlength=nb)
            np.sqrt(np.divide(ss, w, out=np.full_like(ss, np.nan), where=w > 0), out=result[d, k])
    return result, counts


def constrained(v):
    return np.stack((v[4], np.hypot(v[5], v[6]), v[7]))


def bulk_points(xy, levs):
    """Fixed physical exclusion: four low-grid coarse cells from either face."""
    x, y = np.abs(xy[:, 0]), xy[:, 1]
    selected = y > 4*(2./2**levs)
    for lev, a in ((1, 64.), (2, 32.), (3, 16.), (4, 8.), (5, 4.), (6, 2.)):
        inside = (x <= a) & (y <= a)
        distance = np.where(inside, np.minimum(a-x, a-y), np.hypot(np.maximum(x-a, 0), np.maximum(y-a, 0)))
        selected &= distance > 4*(2./2**(lev-1))
    return selected


def selfcheck():
    for n in (6, 8):
        q = np.array([-.5, .125, .25, .375, .625, .75, .875, 3.])
        first = np.floor(q).astype(int)-(n//2-1)
        w = weights(q, first, n)
        for power in range(n):
            assert np.allclose(np.sum(w*(first[:, None]+np.arange(n))**power, axis=1), q**power, rtol=2e-11, atol=2e-11)
    assert Q4 == 5.0625
    xy = np.array([[-256+.1, 0.], [-256+.3, .012], [-256+.6, .35], [-256+.9, .9]])
    x, y = np.meshgrid((np.arange(32)+.5)/32, (np.arange(32)+.5)/32)
    g = dict(h=1/32, lo=np.array([0, 0]), a=np.zeros((10, 32, 32)), bounds=(-256., -255., 0., 1.))
    for n in (6, 8):
        for a in range(n):
            for b in range(n):
                g['a'][:] = 0.
                c = 0 if b % 2 == 0 else 6
                g['a'][c] = x**a*y**b
                v, _ = sample([g], xy, n)
                expected = (xy[:, 0]+256.)**a*xy[:, 1]**b
                assert np.allclose(v[c], expected, rtol=3e-13, atol=3e-13), (n, a, b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path)
    parser.add_argument('--times', default='all')
    parser.add_argument('--output', default='t5')
    parser.add_argument('--cache', type=Path, default=Path('/private/tmp/ems-t5/t5-maps.npz'))
    args = parser.parse_args()
    selfcheck()
    requested = np.arange(0., 10.0001, .5) if args.times == 'all' else np.array([float(t) for t in args.times.split(',')])
    paths = []
    for scale in SCALES:
        ps = sorted((args.root/f'ref-{scale}'/'plt').glob('*.hdf5'))
        assert len(ps) == 21
        paths.append(ps)
    vi = [COMPONENTS.index(v) for v in VARIABLES]
    data, audit, layout, summaries = defaultdict(list), [], [], []
    geometry = {}
    started = time.monotonic()
    for t in requested:
        frame = int(round(t*2))
        grids = []
        for i, scale in enumerate(SCALES):
            p = paths[i][frame]
            tt, g = read(p)
            assert abs(tt-t) < 1e-8
            if scale not in geometry:
                geometry[scale] = [(v['h'], v['bounds'], v['boxes'].copy()) for v in g]
            assert len(g) == 7
            for v, (h, bounds, b) in zip(g, geometry[scale]):
                assert v['h'] == h and v['bounds'] == bounds and np.array_equal(v['boxes'], b)
            grids.append(g)
            nx, nv, _, _ = nodes(g, maximum=364.)
            audit.append(dict(scale=scale, time_M=t, path=str(p), bytes=p.stat().st_size,
                              sha256=hashlib.sha256(p.read_bytes()).hexdigest(), components=';'.join(COMPONENTS),
                              fixed_geometry=True, nonfinite=0,
                              min_chi_all_valid=min(float(v['a'][0].min()) for v in g),
                              min_lapse_all_valid=min(float(v['a'][1].min()) for v in g),
                              min_chi_uncovered=float(nv[0].min()), min_lapse_uncovered=float(nv[1].min()),
                              chi_at_or_below_floor_uncovered=int((nv[0] <= 1e-12).sum()),
                              lapse_at_or_below_floor_uncovered=int((nv[1] <= 1e-12).sum()),
                              max_abs_GaussB=max(float(np.abs(v['a'][8]).max()) for v in g)))
            if frame == 0:
                for lev, v in enumerate(g):
                    a, b, c, d = v['bounds']
                    layout.append(dict(scale=scale, level=lev, dx_M=v['h'], boxes=len(v['boxes']),
                                       xlo_M=a, xhi_M=b, ylo_M=c, yhi_M=d,
                                       axis_boundary_radius_M=b, diagonal_boundary_radius_M=np.hypot(b, d)))
        xy, low, levs, interior = nodes(grids[0])
        mid, lm = sample(grids[1], xy, 6)
        high, lh = sample(grids[2], xy, 6)
        assert np.array_equal(lm, levs) and np.array_equal(lh, levs)
        mid8, _ = sample(grids[1], xy, 8)
        high8, _ = sample(grids[2], xy, 8)
        dlm, dmh = low[vi]-mid[vi], mid[vi]-high[vi]
        d8lm, d8mh = low[vi]-mid8[vi], mid8[vi]-high8[vi]
        native = []
        for g in grids:
            nx, nv, _, _ = nodes(g)
            nr, _ = rms_bins(constrained(nv), nx)
            native.append(nr)
        common = [rms_bins(constrained(v), xy)[0] for v in (low, mid, high)]
        cinterior = [rms_bins(constrained(v), xy, interior)[0] for v in (low, mid, high)]
        diff = [rms_bins(v, xy)[0] for v in (dlm, dmh)]
        diff8 = [rms_bins(v, xy)[0] for v in (d8lm, d8mh)]
        mismatch = rms_bins(dlm-Q4*dmh, xy)[0]
        sensitivity = [rms_bins(v, xy)[0] for v in (dlm-d8lm, dmh-d8mh)]
        counts = rms_bins(constrained(low), xy)[1]
        line_xy = np.concatenate((np.stack((LINE_RADII, np.zeros_like(LINE_RADII)), axis=1),
                                  np.stack((-LINE_RADII, np.zeros_like(LINE_RADII)), axis=1),
                                  np.stack((np.zeros_like(LINE_RADII), LINE_RADII), axis=1),
                                  np.stack((LINE_RADII*np.cos(np.pi/4), LINE_RADII*np.sin(np.pi/4)), axis=1)))
        lines, lines8 = [], []
        for g in grids:
            v, _ = sample(g, line_xy, 6)
            v8, _ = sample(g, line_xy, 8)
            lines.append(v.reshape(10, 4, len(LINE_RADII)).transpose(1, 0, 2))
            lines8.append(v8.reshape(10, 4, len(LINE_RADII)).transpose(1, 0, 2))
        outer_xy = np.concatenate((np.stack((OUTER_RADII, np.zeros_like(OUTER_RADII)), axis=1),
                                   np.stack((-OUTER_RADII, np.zeros_like(OUTER_RADII)), axis=1),
                                   np.stack((np.zeros_like(OUTER_RADII), OUTER_RADII), axis=1)))
        outer = [sample(g, outer_xy, 6)[0].reshape(10, 3, len(OUTER_RADII)).transpose(1, 0, 2) for g in grids]
        gx, gv, _, _ = nodes(grids[0], maximum=364.)
        gm, _ = sample(grids[1], gx, 6)
        gh, _ = sample(grids[2], gx, 6)
        global_common = [rms_bins(constrained(v), gx, edges=GLOBAL_EDGES)[0] for v in (gv, gm, gh)]
        global_diff = [rms_bins(v, gx, edges=GLOBAL_EDGES)[0] for v in (gv[vi]-gm[vi], gm[vi]-gh[vi])]
        global_native = []
        for g in grids:
            px, pv, _, _ = nodes(g, maximum=364.)
            global_native.append(rms_bins(constrained(pv), px, edges=GLOBAL_EDGES)[0])
        keys = dict(time=t, common=np.stack(common), native=np.stack(native), interior=np.stack(cinterior),
                    diff=np.stack(diff), diff8=np.stack(diff8), mismatch=mismatch,
                    sensitivity=np.stack(sensitivity), counts=counts,
                    lines=np.stack(lines), lines8=np.stack(lines8), outer=np.stack(outer),
                    global_common=np.stack(global_common), global_native=np.stack(global_native),
                    global_diff=np.stack(global_diff))
        for k, v in keys.items():
            data[k].append(v)
        # Original masks, using only native fields; compare the supplied reductions.
        for si, (scale, g) in enumerate(zip(SCALES, grids)):
            nx, nv, native_levels, _ = nodes(g)
            r = np.hypot(nx[:, 0], nx[:, 1])
            for mask, a, b in (('ring', 2., 4.), ('far', 4., 8.), ('near_hole', .1, 2.),
                               ('puncture_core', 0., .1), ('inner_0p1_0p25', .1, .25),
                               ('inner_0p25_0p5', .25, .5), ('inner_0p5_1', .5, 1.),
                               ('wake_2p5_3', 2.5, 3.), ('wake_3_3p5', 3., 3.5),
                               ('wake_4p5_5p5', 4.5, 5.5), ('wake_6p5_7p5', 6.5, 7.5)):
                sel = (r >= a) & (r <= b)
                vals = constrained(nv)
                for ci, name in enumerate(CONSTRAINTS):
                    rms = np.sqrt(np.sum(nx[sel, 1]*vals[ci, sel]**2)/np.sum(nx[sel, 1]))
                    summaries.append(dict(time_M=t, scale=scale, mask=mask, constraint=name,
                                          cells=int(sel.sum()), cylindrical_rms=rms))
                    if mask in ('ring', 'far', 'wake_2p5_3', 'wake_4p5_5p5'):
                        strict = sel & bulk_points(nx, native_levels)
                        assert strict.any()
                        rms = np.sqrt(np.sum(nx[strict, 1]*vals[ci, strict]**2)/np.sum(nx[strict, 1]))
                        summaries.append(dict(time_M=t, scale=scale, mask='bulk_'+mask, constraint=name,
                                              cells=int(strict.sum()), cylindrical_rms=rms))
        print(f't={t:g} read/sample complete; elapsed {time.monotonic()-started:.2f}s', flush=True)
    # Validate original shell norms against the independent on-cluster reduction.
    for si, scale in enumerate(SCALES):
        published = {(round(float(v['time_M']), 8), v['mask'], v['constraint']): v for v in
                     csv.DictReader((args.root/f'ref-{scale}'/'constraint-rms.csv').open()) if v['cell_class'] == 'all'}
        for v in summaries:
            key = (v['time_M'], v['mask'], v['constraint'])
            if v['scale'] == scale and key in published:
                old = published[key]
                assert int(old['cells']) == v['cells']
                assert np.isclose(float(old['cylindrical_rms']), v['cylindrical_rms'], rtol=2e-12, atol=1e-30)
    args.cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.cache, radii=RADII, global_radii=GLOBAL_RADII, line_radii=LINE_RADII,
                        outer_radii=OUTER_RADII, directions=DIRECTIONS,
                        variables=VARIABLES, constraints=CONSTRAINTS,
                        **{k: np.array(v) for k, v in data.items()})
    save(args.output+'-input-audit.csv', audit)
    save(args.output+'-variable-availability.csv', [
        dict(variable=v, role='EVOLVED' if v in EVOLVED_COMPONENTS else 'DIAGNOSTIC',
             status='SAVED' if v in COMPONENTS else 'UNAVAILABLE',
             saved_times_M='0:0.5:10' if v in COMPONENTS else '',
             basis='63 plot component inventories; pulled data contain no checkpoint or full-field file')
        for v in EVOLVED_COMPONENTS+tuple(c for c in COMPONENTS if c not in EVOLVED_COMPONENTS)])
    if layout:
        save(args.output+'-layout.csv', layout)
    save(args.output+'-mask-norms.csv', summaries)
    save(args.output+'-analysis-performance.csv', [dict(wall_seconds=time.monotonic()-started,
                                                      peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)])
    print('DONE elapsed', time.monotonic()-started, flush=True)


if __name__ == '__main__':
    main()
