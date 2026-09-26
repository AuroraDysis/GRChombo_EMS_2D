#!/usr/bin/env python3
"""Analyze the C3 Chombo plots. Run with: uv run --with h5py --with numpy python analyze.py RUN_ROOT."""

import csv
import math
import sys
from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np

OUT = Path(__file__).parent
FIELDS = ("chi", "h11", "K", "A11", "lapse", "shift1", "phi", "Pi", "Ex")
RESIDUALS = ("Ham", "Mom1", "Mom2", "GaussE", "GaussB", "Theta")
TIMES = (0, 0.0625, 0.125, 0.25)


@lru_cache(None)
def plot(path, level=0):
    with h5py.File(path) as f:
        names = [f.attrs[f"component_{i}"].decode() for i in range(int(f.attrs["num_components"]))]
        g = f[f"level_{level}"]
        dx, time = float(g.attrs["dx"]), float(g.attrs["time"])
        domain = g.attrs["prob_domain"]
        nx, ny = int(domain["hi_i"]) + 1, int(domain["hi_j"]) + 1
        a = np.full((len(names), ny, nx), np.nan)
        covered = np.zeros((ny, nx), dtype=bool)
        flat = g["data:datatype=0"][:]
        offsets = g["data:offsets=0"][:]
        for box, lo, hi in zip(g["boxes"][:], offsets[:-1], offsets[1:]):
            x0, y0 = int(box["lo_i"]), int(box["lo_j"])
            x1, y1 = int(box["hi_i"]) + 1, int(box["hi_j"]) + 1
            assert hi - lo == len(names) * (x1 - x0) * (y1 - y0)
            a[:, y0:y1, x0:x1] = flat[lo:hi].reshape(len(names), y1 - y0, x1 - x0)
            covered[y0:y1, x0:x1] = True
        assert np.isfinite(a[:, covered]).all(), (path, level, "nonfinite plotted field")
        return {name: a[i] for i, name in enumerate(names)}, covered, dx, time


def file_at(root, case, t):
    files = sorted((root / case / "plt").glob("*.hdf5"))
    for path in files:
        if abs(plot(path)[3] - t) < 1e-12:
            return path
    raise FileNotFoundError((case, t, files))


def mask(covered, dx, centre, inner=0.75, outer=1.5, stretch=1.0):
    ny, nx = covered.shape
    x = (np.arange(nx) + 0.5) * dx - centre
    y = (np.arange(ny) + 0.5) * dx
    r = np.hypot(stretch * x[None, :], y[:, None])
    return covered & (r >= inner) & (r <= outer)


def norms(a):
    return math.sqrt(float(np.mean(a * a))), float(np.max(np.abs(a)))


def write(name, rows):
    assert rows
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def interp_six(a, dx, x, y):
    """Tensor-product degree-five Lagrange interpolation to physical cell centres."""
    def nodes_weights(z, n):
        position = z / dx - 0.5
        start = np.clip(np.floor(position).astype(int) - 2, 0, n - 6)
        nodes = start[:, None] + np.arange(6)
        weights = np.ones_like(nodes, dtype=float)
        for k in range(6):
            for j in range(6):
                if j != k:
                    weights[:, k] *= (position - nodes[:, j]) / (nodes[:, k] - nodes[:, j])
        return nodes, weights
    xi, wx = nodes_weights(x, a.shape[1])
    yi, wy = nodes_weights(y, a.shape[0])
    assert np.isfinite(a[yi[:, :, None], xi[:, None, :]]).all()
    return np.einsum("pa,pb,pab->p", wy, wx, a[yi[:, :, None], xi[:, None, :]])


def single_constraints(root):
    rows = []
    previous = {}
    for case, n, eta in (("n16", 16, 0), ("n32", 32, 0), ("n64", 64, 0),
                         ("boosted", 32, 0.20273255)):
        for t in TIMES:
            data, covered, dx, _ = plot(file_at(root, case, t))
            region = mask(covered, dx, 4, stretch=math.cosh(eta))
            row = {"case": case, "N": n, "eta": eta, "time": t, "count": int(region.sum()),
                   "chi_below_floor_all": int((data["chi"][covered] < 1e-12).sum()),
                   "lapse_below_floor_all": int((data["lapse"][covered] < 1e-12).sum()),
                   "chi_below_floor_mask": int((data["chi"][region] < 1e-12).sum()),
                   "lapse_below_floor_mask": int((data["lapse"][region] < 1e-12).sum()),
                   "min_chi": float(data["chi"][covered].min()),
                   "min_lapse": float(data["lapse"][covered].min())}
            for q in RESIDUALS:
                l2, linf = norms(data[q][region])
                row[q + "_L2"], row[q + "_Linf"] = l2, linf
                prior = previous.get((n // 2, t, q)) if case != "boosted" else None
                row[q + "_pL2"] = math.log2(prior[0] / l2) if prior and prior[0] and l2 else ""
                row[q + "_pLinf"] = math.log2(prior[1] / linf) if prior and prior[1] and linf else ""
                if case != "boosted":
                    previous[n, t, q] = l2, linf
            rows.append(row)
    write("single-constraints.csv", rows)


def field_convergence(root):
    rows = []
    for t in TIMES:
        p16 = plot(file_at(root, "n16", t))
        p32 = plot(file_at(root, "n32", t))
        p64 = plot(file_at(root, "n64", t))
        region = mask(p16[1], p16[2], 4)
        jj, ii = np.nonzero(region)
        x, y = (ii + 0.5) * p16[2], (jj + 0.5) * p16[2]
        for field in FIELDS:
            coarse = p16[0][field][jj, ii]
            medium = interp_six(p32[0][field], p32[2], x, y)
            fine = interp_six(p64[0][field], p64[2], x, y)
            a2, ai = norms(coarse - medium)
            b2, bi = norms(medium - fine)
            rows.append({"time": t, "field": field, "count": len(ii),
                         "difference_16_32_L2": a2, "difference_32_64_L2": b2,
                         "Q_L2": a2 / b2 if b2 else "",
                         "p_L2": math.log2(a2 / b2) if a2 and b2 else "",
                         "difference_16_32_Linf": ai, "difference_32_64_Linf": bi,
                         "Q_Linf": ai / bi if bi else "",
                         "p_Linf": math.log2(ai / bi) if ai and bi else ""})
    write("field-convergence.csv", rows)


def common_constraints(root):
    rows = []
    for t in TIMES:
        plots = [plot(file_at(root, case, t)) for case in ("n16", "n32", "n64")]
        jj, ii = np.nonzero(mask(plots[0][1], plots[0][2], 4))
        x, y = (ii + 0.5) * plots[0][2], (jj + 0.5) * plots[0][2]
        for q in RESIDUALS:
            samples = [plots[0][0][q][jj, ii]]
            samples += [interp_six(p[0][q], p[2], x, y) for p in plots[1:]]
            values = [norms(s) for s in samples]
            row = {"time": t, "quantity": q, "count": len(ii)}
            for n, value in zip((16, 32, 64), values):
                row[f"N{n}_L2"], row[f"N{n}_Linf"] = value
            for k, pair in enumerate(((16, 32), (32, 64))):
                for j, name in enumerate(("L2", "Linf")):
                    old, new = values[k][j], values[k + 1][j]
                    row[f"p{pair[0]}_{pair[1]}_{name}"] = math.log2(old / new) if old and new else ""
            rows.append(row)
    write("common-constraints.csv", rows)


def binary_field_convergence(root):
    rows = []
    for t in (0, 0.25):
        plots = [plot(file_at(root, case, t)) for case in ("binary", "binary-n32", "binary-n64")]
        _, covered, dx, _ = plots[0]
        ny, nx = covered.shape
        xgrid = (np.arange(nx) + 0.5) * dx - 20
        ygrid = (np.arange(ny) + 0.5) * dx
        radius = np.minimum(np.hypot(xgrid[None, :] - 16, ygrid[:, None]),
                            np.hypot(xgrid[None, :] + 16, ygrid[:, None]))
        jj, ii = np.nonzero(covered & (radius >= 0.75) & (radius <= 1.5))
        x, y = (ii + 0.5) * dx, (jj + 0.5) * dx
        for field in FIELDS:
            samples = [plots[0][0][field][jj, ii]]
            samples += [interp_six(p[0][field], p[2], x, y) for p in plots[1:]]
            a2, ai = norms(samples[0] - samples[1])
            b2, bi = norms(samples[1] - samples[2])
            rows.append({"time": t, "field": field, "count": len(ii),
                         "difference_16_32_L2": a2, "difference_32_64_L2": b2,
                         "Q_L2": a2 / b2 if b2 else "",
                         "p_L2": math.log2(a2 / b2) if a2 and b2 else "",
                         "difference_16_32_Linf": ai, "difference_32_64_Linf": bi,
                         "Q_Linf": ai / bi if bi else "",
                         "p_Linf": math.log2(ai / bi) if ai and bi else ""})
    write("binary-field-convergence.csv", rows)


def binary_and_amr(root):
    binary_rows, symmetry_rows, amr_rows = [], [], []
    odd = {"h12", "A12", "Gamma1", "shift1", "B1", "By", "Bz", "Ex", "Mom1", "Lambda"}
    for t in (0, 0.25):
        data, covered, dx, _ = plot(file_at(root, "binary", t))
        ny, nx = covered.shape
        x = (np.arange(nx) + 0.5) * dx - 20
        y = (np.arange(ny) + 0.5) * dx
        radius = np.minimum(np.hypot(x[None, :] - 16, y[:, None]),
                            np.hypot(x[None, :] + 16, y[:, None]))
        shell = covered & (radius >= 0.75) & (radius <= 1.5)
        row = {"time": t, "shell_count": int(shell.sum()),
               "chi_below_floor_all": int((data["chi"][covered] < 1e-12).sum()),
               "lapse_below_floor_all": int((data["lapse"][covered] < 1e-12).sum()),
               "chi_below_floor_shell": int((data["chi"][shell] < 1e-12).sum()),
               "lapse_below_floor_shell": int((data["lapse"][shell] < 1e-12).sum())}
        for q in RESIDUALS:
            row[q + "_L2"], row[q + "_Linf"] = norms(data[q][shell])
        binary_rows.append(row)
        for q, a in data.items():
            sign = -1 if q in odd else 1
            difference = np.abs(a - sign * a[:, ::-1])
            symmetry_rows.append({"time": t, "field": q, "parity": sign,
                                  "max_abs": float(np.max(difference[covered])),
                                  "max_normalized": float(np.max(difference[covered]) /
                                                          (1 + np.max(np.abs(a[covered]))))})
    for path in sorted((root / "amr" / "plt").glob("*.hdf5")):
        with h5py.File(path) as f:
            levels = sorted(int(k.split("_")[1]) for k in f if k.startswith("level_"))
        for level in levels:
            data, covered, dx, t = plot(path, level)
            shell = mask(covered, dx, 8)
            outside = mask(covered, dx, 8, inner=0.25, outer=math.inf)
            row = {"time": t, "level": level, "dx": dx, "cells": int(covered.sum()),
                   "shell_cells": int(shell.sum()),
                   "chi_below_floor_all": int((data["chi"][covered] < 1e-12).sum()),
                   "lapse_below_floor_all": int((data["lapse"][covered] < 1e-12).sum()),
                   "chi_below_floor_outside_R025": int((data["chi"][outside] < 1e-12).sum()),
                   "lapse_below_floor_outside_R025": int((data["lapse"][outside] < 1e-12).sum())}
            for q in RESIDUALS:
                row[q + "_L2"], row[q + "_Linf"] = norms(data[q][shell]) if shell.any() else ("", "")
            amr_rows.append(row)
    write("binary-constraints.csv", binary_rows)
    write("binary-symmetry.csv", symmetry_rows)
    write("amr-levels.csv", amr_rows)


def snapshot_status(root):
    rows = []
    for case, centre, binary in (("n16", 4, False), ("n32", 4, False),
                                 ("n64", 4, False), ("boosted", 4, False),
                                 ("binary", 20, True), ("binary-n32", 20, True),
                                 ("binary-n64", 20, True), ("amr", 8, False)):
        for path in sorted((root / case / "plt").glob("*.hdf5")):
            with h5py.File(path) as f:
                levels = sorted(int(k.split("_")[1]) for k in f if k.startswith("level_"))
            for level in levels:
                data, covered, dx, t = plot(path, level)
                ny, nx = covered.shape
                x = (np.arange(nx) + 0.5) * dx - centre
                y = (np.arange(ny) + 0.5) * dx
                if binary:
                    r = np.minimum(np.hypot(x[None, :] - 16, y[:, None]),
                                   np.hypot(x[None, :] + 16, y[:, None]))
                else:
                    r = np.hypot(x[None, :], y[:, None])
                outside = covered & (r > 0.25)
                rows.append({"case": case, "time": t, "level": level, "dx": dx,
                             "cells": int(covered.sum()),
                             "chi_below_floor": int((data["chi"][covered] < 1e-12).sum()),
                             "lapse_below_floor": int((data["lapse"][covered] < 1e-12).sum()),
                             "chi_below_floor_outside_R025": int((data["chi"][outside] < 1e-12).sum()),
                             "lapse_below_floor_outside_R025": int((data["lapse"][outside] < 1e-12).sum()),
                             "min_chi": float(data["chi"][covered].min()),
                             "min_lapse": float(data["lapse"][covered].min())})
    write("snapshot-status.csv", rows)


def horizon(root):
    path = root / "rh" / "rh_surf_0.dat"
    if not path.exists():
        return
    rows = []
    for line in path.read_text().splitlines()[1:]:
        parts = line.split()
        if not parts:
            continue
        area = float(parts[4])
        rows.append({"time": float(parts[0]), "coordinate_radius": float(parts[3]),
                     "area": area, "areal_radius": math.sqrt(area / (4 * math.pi)),
                     "expected_areal_radius": 1.7505898524629477,
                     "theta_plus_mean": float(parts[7]), "solver_error": float(parts[9]),
                     "mode": parts[-1]})
    write("rh-results.csv", rows)


def fixed_strip(log):
    lines = Path(log).read_text().splitlines()
    rows, previous = [], {}
    for k in range(0, len(lines), 4):
        header = lines[k].split()
        n = int(header[0].split("=")[1])
        eta = f'{float(header[1].split("=")[1]):.8f}'.rstrip('0').rstrip('.')
        names = lines[k + 1].split(',')
        for line in lines[k + 2:k + 4]:
            row = {"N": n, "eta": eta, **dict(zip(names, line.split(',')))}
            if row["region"] != "fixed_axis_strip":
                continue
            prior = previous.get(eta)
            for q in ("H", "Mx", "My", "GaussE", "GaussB"):
                for norm in ("L2", "Linf"):
                    col = q + "_" + norm
                    row[q + "_p" + norm] = (math.log2(float(prior[col]) / float(row[col]))
                                                if prior and float(prior[col]) and float(row[col]) else "")
            rows.append(row)
            previous[eta] = row
    assert len(rows) == 9
    write("fixed-axis-strip.csv", rows)


if __name__ == "__main__":
    grid = (np.arange(20) + 0.5) / 16
    polynomial = grid[None, :] ** 5 + grid[:, None] ** 5
    test_x, test_y = np.array([0.11, 0.47, 0.93]), np.array([0.07, 0.55, 1.01])
    assert np.max(np.abs(interp_six(polynomial, 1 / 16, test_x, test_y) -
                         (test_x ** 5 + test_y ** 5))) < 1e-12
    root = Path(sys.argv[1])
    single_constraints(root)
    field_convergence(root)
    common_constraints(root)
    binary_and_amr(root)
    binary_field_convergence(root)
    snapshot_status(root)
    horizon(root)
    if len(sys.argv) > 2:
        fixed_strip(sys.argv[2])
    print("wrote evolution CSV files")
