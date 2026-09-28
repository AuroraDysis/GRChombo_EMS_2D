"""Audit the reference smoke's finest boxes over the full 5 r_h disk."""
import csv
from pathlib import Path

import h5py
import numpy as np

root = Path(__file__).resolve().parent
mass = 0.5710876672204076
with (root / "s2-coverage.csv").open("w", newline="") as target:
    writer = csv.writer(target)
    writer.writerow(("M_over_h", "finest_boxes", "finest_dx_over_M",
                     "cells_inside_5rh", "uncovered_cells"))
    for n, step in ((32, 48), (64, 32)):
        checkpoint = root / f"s2-reference-{n}-stage{step}" / "chk" / (
            f"EMS_{step:06}.2d.hdf5")
        with h5py.File(checkpoint) as handle:
            grid = handle["level_2"]
            boxes = grid["boxes"][:]
            dx = float(grid.attrs["dx"])
            domain = handle["level_0"].attrs["prob_domain"]
            nx, ny = 4 * (int(domain["hi_i"]) + 1), 4 * (int(domain["hi_j"]) + 1)
            covered = np.zeros((ny, nx), dtype=bool)
            for box in boxes:
                covered[box["lo_j"]:box["hi_j"] + 1,
                        box["lo_i"]:box["hi_i"] + 1] = True
            x = (np.arange(nx) + .5) * dx - 16 * mass
            y = (np.arange(ny) + .5) * dx
            required = x[None, :] ** 2 + y[:, None] ** 2 <= 25
            missing = int(np.count_nonzero(required & ~covered))
            writer.writerow((n, len(boxes), dx / mass,
                             int(np.count_nonzero(required)), missing))
            assert missing == 0
