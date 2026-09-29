"""Compare narrow and full-path fine coverage on common fine cells."""
import csv
import math
import sys
from pathlib import Path

import h5py
import numpy as np

L = 36.54961070210609
FIELDS = ("lapse", "K", "Gamma1", "chi")


def cells(path, level, field):
    with h5py.File(path) as f:
        names = [f.attrs[f"component_{i}"].decode()
                 for i in range(int(f.attrs["num_components"]))]
        component = names.index(field)
        group = f[f"level_{level}"]
        values = group["data:datatype=0"][:]
        offsets = group["data:offsets=0"][:]
        result = {}
        for box, start, end in zip(group["boxes"][:], offsets[:-1], offsets[1:]):
            lo_x, lo_y, hi_x, hi_y = map(int, box)
            nx, ny = hi_x - lo_x + 7, hi_y - lo_y + 7
            size = nx * ny
            chunk = values[int(start) + component * size:int(start) + (component + 1) * size]
            assert len(chunk) == size and int(end - start) == size * len(names)
            arr = chunk.reshape(ny, nx)
            for j in range(lo_y, hi_y + 1):
                for i in range(lo_x, hi_x + 1):
                    result[i, j] = float(arr[j - lo_y + 3, i - lo_x + 3])
        return float(f.attrs["time"]), result


def measure(a, b, N, mask):
    dx = L / (N * 2)
    errors = []
    for (i, j), value in a.items():
        x, y = (i + .5) * dx - L / 2, (j + .5) * dx
        if mask == "incoming" and not (1.5 <= x <= 2.5 and y <= .7):
            continue
        if mask == "interface" and not (2.5 <= x <= 3.0 and y <= .7):
            continue
        if (i, j) in b:
            errors.append(value - b[i, j])
    assert errors, (N, mask)
    return len(errors), math.sqrt(sum(e * e for e in errors) / len(errors)), max(map(abs, errors))


def main(root, output):
    rows = []
    for N in (256, 384, 512):
        for mode in ("raw", "relative"):
            for step in (0, 60 * N // 256):
                paths = [Path(root) / f"N{N}" / mode / coverage / "chk" /
                         f"EMS_{step:06d}.2d.hdf5" for coverage in ("narrow", "wide")]
                for field in FIELDS:
                    ta, a = cells(paths[0], 1, field)
                    tb, b = cells(paths[1], 1, field)
                    assert abs(ta - tb) < 1e-12
                    for mask in ("incoming", "interface"):
                        count, rms, maximum = measure(a, b, N, mask)
                        rows.append((N, mode, step, ta, field, mask, count, rms, maximum))
    with open(output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("N", "transfer", "step", "time", "field", "mask", "count", "RMS", "Linf"))
        writer.writerows(rows)
    for mode in ("raw", "relative"):
        for field in FIELDS:
            chosen = [r for r in rows if r[1] == mode and r[4] == field and
                      r[5] == "incoming" and r[2] > 0]
            orders = [math.log(a[7] / b[7]) / math.log(b[0] / a[0])
                      for a, b in zip(chosen, chosen[1:])]
            print(mode, field, "RMS", *(f"{r[7]:.6g}" for r in chosen),
                  "orders", *(f"{p:.3f}" for p in orders))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
