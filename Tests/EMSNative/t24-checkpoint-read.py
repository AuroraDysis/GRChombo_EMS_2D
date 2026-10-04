#!/usr/bin/env python3
"""Read numerical checkpoint fields and rank logs; never advance or read static data."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import resource
import time

import h5py
import numpy as np

KEYS = ('lo_i', 'lo_j', 'hi_i', 'hi_j')


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(p, rows):
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('evidence', type=Path)
    ap.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument('--cell', type=int, nargs=2, default=(1310826, 27))
    ap.add_argument('--cell-level', type=int, default=10)
    ap.add_argument('--levels', type=int, nargs=2, default=(8, 12))
    ap.add_argument('--centre', type=float, nargs=2, default=(2240., 0.))
    ap.add_argument('--hole-radius', type=float, default=.075)
    ap.add_argument('--cell-radius', type=int, default=6)
    ap.add_argument('--snapshot-region', type=float, nargs=4, default=(-.4, 0., .4, .15))
    args = ap.parse_args()
    start = time.monotonic()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    verified = {}
    for attempt in ('ev1', 'ev2'):
        manifest = args.evidence / attempt / ('stall-evidence-' + attempt[-1] + '.sha256')
        count = 0
        for line in manifest.read_text().splitlines():
            digest, name = line.split(None, 1)
            p = manifest.parent / name.lstrip('*')
            assert sha(p) == digest, p
            count += 1
        verified[attempt] = dict(files=count, manifest_sha256=sha(manifest))

    checkpoint = args.evidence / 'chk/EMS_000152.2d.hdf5'
    track = args.evidence / 'ev1/runs/exp-0024/merger/punctures.dat'
    tracks = np.loadtxt(track)
    row = tracks[np.argmin(abs(tracks[:, 0] - 133.))]
    assert row[0] == 133.
    holes = [tuple(row[1:3]), tuple(row[3:5])]
    extrema, neighbourhood, nearest, negatives, coverage_rows, storage = [], [], [], [], [], []
    cached = {}
    with h5py.File(checkpoint, 'r') as f:
        names = [f.attrs[f'component_{i}'].decode() for i in range(int(f.attrs['num_components']))]
        nlevels = int(f.attrs['num_levels'])
        layouts = [[tuple(int(b[k]) for k in KEYS) for b in f[f'level_{l}/boxes'][:]]
                   for l in range(nlevels)]
        href = float(f[f'level_{args.cell_level}'].attrs['dx'])
        bad = tuple((i + .5) * href for i in args.cell)
        mirror = (2 * args.centre[0] - bad[0], bad[1])
        points = dict(failing_cell=bad, mirrored_cell=mirror, puncture1=holes[0], puncture2=holes[1])
        desired = 'chi lapse K h11 h12 h22 hww A11 A12 A22 Aww shift1 shift2 Theta'.split()
        for level in range(nlevels):
            g = f[f'level_{level}']
            h = float(g.attrs['dx'])
            offsets = g['data:offsets=0'][:]
            gx, gy = map(int, g['data_attributes'].attrs['outputGhost'])
            children = [(a // 2, b // 2, c // 2, d // 2) for a, b, c, d in layouts[level + 1]] if level + 1 < nlevels else []
            accum = {}
            def add(region, mask, a, xx, yy, box):
                if not mask.any():
                    return
                for name in desired:
                    v = a[names.index(name)][mask]
                    key = (region, name)
                    z = accum.setdefault(key, dict(n=0, sumsq=0., lo=float('inf'), hi=-float('inf'),
                         nonfinite=0, nonpositive=0, floor=0, below_floor=0, minimum_cell=''))
                    z['n'] += len(v)
                    z['nonfinite'] += int((~np.isfinite(v)).sum())
                    z['sumsq'] += float(np.sum(v * v))
                    z['nonpositive'] += int((v <= 0).sum())
                    z['floor'] += int((v == 1e-12).sum())
                    z['below_floor'] += int((v < 1e-12).sum())
                    lo, hi = float(np.min(v)), float(np.max(v))
                    if lo < z['lo']:
                        ij = np.argwhere(mask)[int(np.argmin(v))]
                        z['lo'] = lo
                        z['minimum_cell'] = f'{box[0]+int(ij[1])},{box[1]+int(ij[0])}'
                    z['hi'] = max(z['hi'], hi)

            snapshot_cells = snapshot_ghost_bound = 0
            for bi, box in enumerate(layouts[level]):
                x0, y0, x1, y1 = box
                nx, ny = x1 - x0 + 1, y1 - y0 + 1
                raw = g['data:datatype=0'][int(offsets[bi]):int(offsets[bi + 1])]
                full = raw.reshape(len(names), ny + 2 * gy, nx + 2 * gx)
                a = full[:, gy:gy + ny, gx:gx + nx].copy()
                if args.levels[0] <= level <= args.levels[1]:
                    cached.setdefault(level, []).append((box, a))
                xx, yy = np.meshgrid((np.arange(x0, x1 + 1) + .5) * h,
                                     (np.arange(y0, y1 + 1) + .5) * h)
                covered = np.zeros((ny, nx), dtype=bool)
                for cx0, cy0, cx1, cy1 in children:
                    if cx1 < x0 or cx0 > x1 or cy1 < y0 or cy0 > y1:
                        continue
                    covered[max(0, cy0-y0):min(ny, cy1-y0+1), max(0, cx0-x0):min(nx, cx1-x0+1)] = True
                add('all_valid', np.ones((ny, nx), dtype=bool), a, xx, yy, box)
                add('active_valid', ~covered, a, xx, yy, box)
                add('covered_valid', covered, a, xx, yy, box)
                if args.levels[0] <= level <= args.levels[1]:
                    ghost = np.ones(full.shape[1:], dtype=bool)
                    ghost[gy:gy+ny,gx:gx+nx] = False
                    gxx, gyy = np.meshgrid((np.arange(x0-gx,x1+gx+1)+.5)*h,
                                          (np.arange(y0-gy,y1+gy+1)+.5)*h)
                    add('saved_ghosts_with_duplicates',ghost,full,gxx,gyy,(x0-gx,y0-gy,x1+gx,y1+gy))
                for name in ('chi', 'lapse'):
                    v = a[names.index(name)]
                    for iy, ix in np.argwhere(v <= 0)[:8]:
                        negatives.append(dict(level=level, field=name, i=x0+int(ix), j=y0+int(iy),
                            x_M=float(xx[iy, ix]-args.centre[0]), y_M=float(yy[iy, ix]-args.centre[1]),
                            value=float(v[iy, ix]), covered_by_child=bool(covered[iy, ix])))
                if not args.levels[0] <= level <= args.levels[1]:
                    continue
                sx0, sy0, sx1, sy1 = args.snapshot_region
                snap = (xx-args.centre[0] >= sx0) & (xx-args.centre[0] <= sx1) & (yy-args.centre[1] >= sy0) & (yy-args.centre[1] <= sy1)
                snapshot_cells += int(snap.sum())
                if snap.any():
                    # Conservative per-box bound: save the entire intersecting FAB, ghosts included.
                    snapshot_ghost_bound += (nx + 2*gx) * (ny + 2*gy)
                for label, (px, py) in points.items():
                    mask = ((xx-px)**2+(yy-py)**2 <= args.hole_radius**2) if label.startswith('puncture') else ((abs(xx-px) <= args.cell_radius*href) & (abs(yy-py) <= args.cell_radius*href))
                    add(label, mask, a, xx, yy, box)
                    i, j = int(np.floor(px / h)), int(np.floor(py / h))
                    if x0 <= i <= x1 and y0 <= j <= y1:
                        v = a[:, j-y0, i-x0]
                        nearest.append(dict(level=level, region=label, i=i, j=j, box=bi,
                            x_M=(i+.5)*h-args.centre[0], y_M=(j+.5)*h-args.centre[1],
                            covered_by_child=bool(covered[j-y0, i-x0]), **{n:float(v[names.index(n)]) for n in names}))
                        for dj in range(-2, 3):
                            for di in range(-2, 3):
                                qi, qj = i+di, j+dj
                                if x0 <= qi <= x1 and y0 <= qj <= y1:
                                    q = a[:, qj-y0, qi-x0]
                                    neighbourhood.append(dict(level=level, region=label, i=qi, j=qj, box=bi,
                                        x_M=(qi+.5)*h-args.centre[0], y_M=(qj+.5)*h-args.centre[1],
                                        covered_by_child=bool(covered[qj-y0, qi-x0]),
                                        **{n:float(q[names.index(n)]) for n in desired}))
            for (region, name), z in accum.items():
                extrema.append(dict(level=level, region=region, field=name, valid_cells=z['n'],
                    minimum=z['lo'], maximum=z['hi'], rms=np.sqrt(z['sumsq']/z['n']),
                    nonfinite=z['nonfinite'], nonpositive=z['nonpositive'], floor_cells=z['floor'],
                    below_floor_cells=z['below_floor'], minimum_cell=z['minimum_cell']))
            if args.levels[0] <= level <= args.levels[1]:
                storage.append(dict(level=level, h_M=h, roi_valid_cells=snapshot_cells,
                    intersecting_full_FAB_cells_with_ghosts=snapshot_ghost_bound,
                    valid_field_bytes=snapshot_cells*len(names)*8, FAB_field_bytes=snapshot_ghost_bound*len(names)*8))
                for label, (px, py) in points.items():
                    i, j = int(np.floor(px/h)), int(np.floor(py/h))
                    bb = [bi for bi, (x0,y0,x1,y1) in enumerate(layouts[level]) if x0 <= i <= x1 and y0 <= j <= y1]
                    coverage_rows.append(dict(level=level, region=label, i=i,j=j,covered=bool(bb),box=bb[0] if bb else '',h_M=h))

    def value(level, i, j, field):
        for (x0,y0,x1,y1), a in cached[level]:
            if x0 <= i <= x1 and y0 <= j <= y1:
                return float(a[names.index(field),j-y0,i-x0])
        raise KeyError((level,i,j,field))

    # Replay the existing restriction arithmetic, including its anchor and loop order.
    w = []
    for i in range(6):
        wi = 1.
        for j in range(6):
            if j != i:
                wi *= (.5 - (-2+j)) / (i-j)
        w.append(wi)
    restriction = []
    witnesses = [(args.cell_level,*args.cell,n) for n in desired]
    witnesses += [(r['level'],*map(int,r['minimum_cell'].split(',')),r['field']) for r in extrema
                   if r['region']=='all_valid' and r['field'] in ('chi','lapse') and r['minimum'] < 0
                   and r['level'] in cached and r['level']+1 in cached]
    for level, i, j, field in witnesses:
        try:
            base = value(level+1,2*i,2*j,field)
            correction = 0.
            values = []
            for y in range(6):
                for x in range(6):
                    v = value(level+1,2*i+x-2,2*j+y-2,field)
                    values.append(v)
                    correction += w[x] * w[y] * (v-base)
            predicted = base + correction
            actual = value(level,i,j,field)
            restriction.append(dict(level=level, i=i,j=j,field=field,actual=actual,
                point_restriction=predicted,difference=predicted-actual,
                bit_identical=bool(np.float64(predicted).view(np.uint64)==np.float64(actual).view(np.uint64)),
                stencil_min=min(values),stencil_max=max(values),stencil='36 globally valid child cells'))
        except KeyError:
            restriction.append(dict(level=level,i=i,j=j,field=field,actual=value(level,i,j,field),
                point_restriction='',difference='',bit_identical='',stencil_min='',stencil_max='',
                stencil='not wholly in saved globally valid child cells'))

    witnesses, box_counts, transitions = [], [], []
    for attempt in ('ev1','ev2'):
        for p in sorted((args.evidence/attempt).rglob('pout.*')):
            text = p.read_text(errors='replace')
            if 'Values have become nan.' in text:
                start_dump = text.index('NaNCheck in specific Advance:')
                dump = text[start_dump:]
                witnesses.append(dict(attempt=attempt,path=str(p),dump=dump,
                    dump_sha256=hashlib.sha256(dump.encode()).hexdigest()))
            if p.name == 'pout.0' and 'engine-step000156' in str(p):
                last_time = ''
                previous = {}
                for line in text.splitlines():
                    m = re.search(r'GRAMRLevel::regrid level (\d+) at time (\S+)',line)
                    if m:
                        last_time = m.group(2)
                    m = re.search(r'level\s+(\d+)\s+.*(?:total|boxes)',line)
                    # Preserve full log line instead of assuming a progress-line format.
                    if 'level 10' in line and ('240' in line or '256' in line):
                        box_counts.append(dict(attempt=attempt,last_regrid_time=last_time,line=line))
                    m = re.search(r'GRAMRLevel::advance level (\d+) at time (\S+).*Boxes on this rank: (\d+) / (\d+)',line)
                    if m and args.levels[0] <= int(m.group(1)) <= args.levels[1]:
                        l, n = int(m.group(1)), int(m.group(4))
                        if previous.get(l) != n:
                            transitions.append(dict(attempt=attempt,level=l,printed_time=m.group(2),
                                old_boxes=previous.get(l,''),new_boxes=n,line=line))
                        previous[l] = n
    save(out/'t24-checkpoint-extrema.csv',extrema)
    save(out/'t24-checkpoint-nearest.csv',nearest)
    save(out/'t24-checkpoint-neighbourhood.csv',neighbourhood)
    if negatives:
        save(out/'t24-negative-valid-examples.csv',negatives)
    save(out/'t24-point-restriction.csv',restriction)
    save(out/'t24-cell-coverage.csv',coverage_rows)
    save(out/'t24-recording-size.csv',storage)
    if box_counts:
        save(out/'t24-level10-count-history.csv',box_counts)
    if transitions:
        save(out/'t24-box-count-transitions.csv',transitions)
    receipt = dict(returncode=0, elapsed_s=time.monotonic()-start,
        peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        arguments=vars(args) | dict(evidence=str(args.evidence),output=str(out)),
        input_manifests=verified, checkpoint_sha256=sha(checkpoint),
        witness_count=len(witnesses), witnesses=witnesses,
        restriction_weights=w, hole_coordinates_t133=holes,
        neighbourhood_rows=len(neighbourhood), negative_examples=len(negatives),
        negative_examples_scope='first eight per FAB and field; counts are in extrema CSV',
        snapshot_bytes_with_whole_intersecting_FABs=sum(r['FAB_field_bytes'] for r in storage),
        static_reference_used=False, advances=0)
    (out/'t24-offline-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('witnesses','arguments')},indent=2),flush=True)


if __name__ == '__main__':
    main()
