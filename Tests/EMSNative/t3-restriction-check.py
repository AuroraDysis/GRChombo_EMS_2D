#!/usr/bin/env python3
"""Compare every covered coarse plot cell with arithmetic restriction of fine cells."""
import csv
import sys
from pathlib import Path

import h5py
import numpy as np


def level_arrays(group, ncomp):
    boxes = group['boxes'][:]
    offsets = group['data:offsets=0'][:]
    data = group['data:datatype=0']
    ghost = group['data_attributes'].attrs['outputGhost']
    gx, gy = int(ghost['intvecti']), int(ghost['intvectj'])
    for bi, box in enumerate(boxes):
        x0, y0, x1, y1 = (int(box[k]) for k in ('lo_i', 'lo_j', 'hi_i', 'hi_j'))
        values = np.asarray(data[int(offsets[bi]):int(offsets[bi + 1])])
        values = values.reshape((-1, y1-y0+1+2*gy, x1-x0+1+2*gx))[:ncomp]
        yield (x0, y0, x1, y1), values[:,gy:gy+y1-y0+1,gx:gx+x1-x0+1]


def check(path, label):
    rows = []
    with h5py.File(path) as f:
        ncomp = min(28, int(f.attrs['num_components']))
        names = [f.attrs[f'component_{i}'].decode() for i in range(ncomp)]
        selected = range(28) if ncomp == 28 else [i for i, name in enumerate(names)
                                                if name in ('chi', 'lapse', 'phi', 'Theta')]
        for level in range(int(f.attrs['num_levels']) - 1):
            coarse = list(level_arrays(f[f'level_{level}'], ncomp))
            fine = list(level_arrays(f[f'level_{level+1}'], ncomp))
            fx0 = min(b[0] for b, _ in fine) // 2
            fx1 = max(b[2] for b, _ in fine) // 2
            fy1 = max(b[3] for b, _ in fine) // 2
            errors = [[] for _ in names]
            near_errors = [[] for _ in names]
            for (x0, y0, x1, y1), values in fine:
                assert x0 % 2 == y0 % 2 == 0 and x1 % 2 == y1 % 2 == 1
                assert (x1-x0+1) % 2 == (y1-y0+1) % 2 == 0
                avg = values.reshape(ncomp, (y1-y0+1)//2, 2, (x1-x0+1)//2, 2).mean(axis=(2, 4))
                cx0, cy0, cx1, cy1 = x0//2, y0//2, x1//2, y1//2
                X, Y = np.meshgrid(np.arange(cx0,cx1+1), np.arange(cy0,cy1+1))
                near = (X-fx0 < 2) | (fx1-X < 2) | (fy1-Y < 2)
                matched = 0
                for (bx0, by0, bx1, by1), cv in coarse:
                    ix0, iy0, ix1, iy1 = max(cx0,bx0), max(cy0,by0), min(cx1,bx1), min(cy1,by1)
                    if ix0 > ix1 or iy0 > iy1: continue
                    matched += (ix1-ix0+1) * (iy1-iy0+1)
                    a = avg[:,iy0-cy0:iy1-cy0+1,ix0-cx0:ix1-cx0+1]
                    c = cv[:,iy0-by0:iy1-by0+1,ix0-bx0:ix1-bx0+1]
                    d = np.abs(c-a)
                    nm = near[iy0-cy0:iy1-cy0+1,ix0-cx0:ix1-cx0+1]
                    for k in selected:
                        errors[k].append(d[k].ravel())
                        near_errors[k].append(d[k][nm])
                assert matched == avg.shape[1] * avg.shape[2]
            for k in selected:
                name = names[k]
                d = np.concatenate(errors[k]); nd = np.concatenate(near_errors[k])
                rows.append((label, level, name, d.size, nd.size, float(d.max()), float(np.sqrt(np.mean(d*d))), float(nd.max()), float(np.sqrt(np.mean(nd*nd)))))
    return rows


if __name__ == '__main__':
    output = Path(sys.argv[1])
    specs = (arg.split('=', 1) for arg in sys.argv[2:])
    with output.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(('run','coarse_level','field','covered_cells','near_interface_cells','max_abs_difference','rms_difference','near_max_abs_difference','near_rms_difference'))
        for label, path in specs:
            w.writerows(check(path, label))
