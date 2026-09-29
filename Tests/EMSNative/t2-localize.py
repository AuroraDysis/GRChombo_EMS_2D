#!/usr/bin/env python3
"""Classify t=0 composite constraint cells in existing Chombo plot files."""
import csv
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).parent
RINGS = defaultdict(list)
for row in csv.DictReader((ROOT / "rings.csv").open()):
    RINGS[row["member"]].append(float(row["isotropic_R_M"]))
RH = {"E": 0.0060286129392158804, "B": 0.0014189019194714414,
      "R": 0.63593977642346233}


def measure(name, path, member, center, classes, extra_masks=None):
    if member == "R":
        masks = {"horizon": (RH[member], 2 * RH[member]),
                 "near_hole": (0.1, 2.), "ring": (2., 4.), "far": (4., 8.)}
    else:
        masks = {
            "horizon": (RH[member], 2 * RH[member]),
            "inside_inner_ring": (2 * RH[member], RINGS[member][0]),
            "cavity": (RINGS[member][0], RINGS[member][1]),
            "between_rings": (RINGS[member][0], RINGS[member][2]),
            "outer_ring": (.8 * RINGS[member][2], 1.2 * RINGS[member][2]),
            "far": (4., 8.),
            "outer_boundary_shell": (0., 1.),
        }
    masks.update(extra_masks or {})
    acc = defaultdict(lambda: [0., 0., 0, 0.])
    with h5py.File(path) as f:
        component_names = [f.attrs[f"component_{i}"].decode() for i in range(int(f.attrs["num_components"]))]
        fields = {"Ham": (component_names.index("Ham"),),
                  "Mom": (component_names.index("Mom1"), component_names.index("Mom2")),
                  "GaussE": (component_names.index("GaussE"),)}
        nl = int(f.attrs["num_levels"])
        base_domain = f["level_0"].attrs["prob_domain"]
        domain_length_x = (int(base_domain["hi_i"]) + 1) * float(f["level_0"].attrs["dx"])
        domain_length_y = (int(base_domain["hi_j"]) + 1) * float(f["level_0"].attrs["dx"])
        for level in range(nl):
            g = f[f"level_{level}"]
            h = float(g.attrs["dx"])
            boxes = g["boxes"][:]
            bounds = np.array([[int(b[k]) for k in ("lo_i", "lo_j", "hi_i", "hi_j")] for b in boxes])
            xmin, ymin = bounds[:, :2].min(axis=0)
            xmax, ymax = bounds[:, 2:].max(axis=0)
            rectangle = int(np.sum((bounds[:, 2] - bounds[:, 0] + 1) * (bounds[:, 3] - bounds[:, 1] + 1))) == (xmax - xmin + 1) * (ymax - ymin + 1)
            if not rectangle:
                raise ValueError(f"{name} level {level} union has re-entrant corners; classifier needs generalization")
            finer = f[f"level_{level + 1}/boxes"][:] if level + 1 < nl else []
            offsets = g["data:offsets=0"][:]
            data = g["data:datatype=0"]
            nx_domain = (int(base_domain["hi_i"]) + 1) * (2 ** level)
            ny_domain = (int(base_domain["hi_j"]) + 1) * (2 ** level)
            for bi, (x0, y0, x1, y1) in enumerate(bounds):
                xx = np.arange(x0, x1 + 1)
                yy = np.arange(y0, y1 + 1)
                X, Y = np.meshgrid(xx, yy)
                rho = np.hypot((X + .5) * h - center, (Y + .5) * h)
                outer_distance = np.minimum.reduce(((X + .5) * h,
                    domain_length_x - (X + .5) * h,
                    domain_length_y - (Y + .5) * h))
                valid = np.ones(X.shape, bool)
                for fb in finer:
                    a, c, d, e = [int(fb[k]) for k in ("lo_i", "lo_j", "hi_i", "hi_j")]
                    valid &= ~((X >= a // 2) & (X <= d // 2) & (Y >= c // 2) & (Y <= e // 2))
                v = np.asarray(data[int(offsets[bi]):int(offsets[bi + 1])]).reshape((len(component_names), len(yy), len(xx)))
                edge_x = ((X - xmin < 2) | (xmax - X < 2)) & (level > 0)
                edge_y = ((Y - ymin < 2) | (ymax - Y < 2)) & (level > 0) & (Y >= 2)
                flags = {
                    "all": np.ones(X.shape, bool),
                    "convex_corner": edge_x & edge_y,
                    "patch_edge": (edge_x ^ edge_y) & ~((Y < 2) & edge_x),
                    "axis_patch_junction": (Y < 2) & edge_x,
                    "cartoon_axis": (Y < 2) & ~edge_x,
                    "box_seam": ((X - x0 < 2) | (x1 - X < 2) | (Y - y0 < 2) | (y1 - Y < 2)) & ~(edge_x | edge_y) & (Y >= 2),
                    "outer_boundary": ((X < 2) | (X >= nx_domain - 2) |
                                       (Y >= ny_domain - 2)) & (level == 0),
                }
                flags["box_seam"] &= ~flags["outer_boundary"]
                flags["interior"] = ~(flags["convex_corner"] | flags["patch_edge"] | flags["axis_patch_junction"] | flags["cartoon_axis"] | flags["box_seam"] | flags["outer_boundary"])
                for mask, (lo, hi) in masks.items():
                    in_mask = valid & (outer_distance >= lo) & (outer_distance <= hi) & ((Y + .5) * h > 1.) if mask == "outer_boundary_shell" else valid & (rho >= lo) & (rho <= hi)
                    if not in_mask.any():
                        continue
                    for cell_class in classes:
                        chosen = in_mask & flags[cell_class]
                        if not chosen.any():
                            continue
                        weight = np.broadcast_to((yy + .5)[:, None] * h, X.shape)[chosen]
                        for field, components in fields.items():
                            sq = sum(v[c][chosen] ** 2 for c in components)
                            a = acc[mask, cell_class, field]
                            a[0] += float(np.sum(weight * sq))
                            a[1] += float(np.sum(weight))
                            a[2] += int(np.count_nonzero(chosen))
                            a[3] = max(a[3], float(np.sqrt(sq.max())))
    return [[name, member, mask, cell_class, field, n, math.sqrt(ss / w), maximum]
            for (mask, cell_class, field), (ss, w, n, maximum) in sorted(acc.items())]


if __name__ == "__main__":
    rows = []
    classes = ("all", "convex_corner", "patch_edge", "axis_patch_junction", "cartoon_axis", "box_seam", "outer_boundary", "interior")
    for spec in sys.argv[1:]:
        name, member, center, path = spec.split(":", 3)
        rows += measure(name, path, member, float(center), classes)
    with Path(os.environ.get("T2_LOCALIZE_OUTPUT", ROOT / "t2-localization.csv")).open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["run", "member", "mask", "cell_class", "constraint", "cells", "cylindrical_rms", "maximum"])
        writer.writerows(rows)
