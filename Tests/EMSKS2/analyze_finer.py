#!/usr/bin/env python3
"""Combine the unchanged three-level ladder with the authorized finer grids."""

import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
LADDERS = {"B": (512, 1024, 2048, 4096),
           "E": (512, 1024, 2048, 4096),
           "reference": (32, 64, 128, 256)}
MASKS = ("join", "KS_collar")
FIELDS = ("H", "Mx", "My", "GaussE", "GaussB", "Gamma_error")


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, columns, rows):
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns)
        writer.writeheader()
        writer.writerows(rows)


def order(a, b):
    return math.log2(a / b) if a > 0 and b > 0 else ""


old = {(r["member"], int(r["N"]), r["mask"], r["residual"]): r
       for r in read_csv(HERE / "grid-results.csv")}
results, orders, profiles, peaks, overlays, concentration = [], [], [], [], [], []
max_rerun_error = 0.0
for member, ns in LADDERS.items():
    fine_rows = {(r["mask"], r["residual"]): r for r in
                 read_csv(HERE / f"grid-{member}-{ns[-1]}.out")
                 if r["member"] == member}
    for n in ns[1:-1]:
        rerun = {(r["mask"], r["residual"]): r for r in
                 read_csv(HERE / f"grid-{member}-{n}.out")
                 if r["member"] == member}
        for key, row in rerun.items():
            if key[0] not in MASKS:
                continue
            saved = old[(member, n, *key)]
            for field in ("L2", "Linf"):
                x, y = float(row[field]), float(saved[field])
                error = abs(x - y) / abs(y) if y else abs(x)
                max_rerun_error = max(max_rerun_error, error)
                assert error < 1e-10, (member, n, key, field, error)
    for mask in MASKS:
        for field in FIELDS:
            rows = [old[(member, n, mask, field)] for n in ns[:-1]]
            rows.append(fine_rows[(mask, field)])
            assert all(int(r["chi_floor"]) == 0 and
                       int(r["lapse_floor"]) == 0 for r in rows)
            for n, row in zip(ns, rows):
                results.append({**row, "min_lapse": 1e-9 if n == ns[-1] else 1e-8,
                                "min_chi": 1e-9 if n == ns[-1] else 1e-8})
            l2 = [float(r["L2"]) for r in rows]
            linf = [float(r["Linf"]) for r in rows]
            orders.append({"member": member, "mask": mask, "residual": field,
                           **{f"N{i}": n for i, n in enumerate(ns)},
                           **{f"L2_{i}": v for i, v in enumerate(l2)},
                           **{f"Linf_{i}": v for i, v in enumerate(linf)},
                           **{f"L2_order_{i}{i+1}": order(l2[i], l2[i+1])
                              for i in range(3)},
                           **{f"Linf_order_{i}{i+1}": order(linf[i], linf[i+1])
                              for i in range(3)}})
    by_n = {n: read_csv(HERE / f"profile-{member}-{n}.csv") for n in ns[1:]}
    for n, rows in by_n.items():
        assert len(rows) == 40 and all(int(r["count"]) > 0 for r in rows)
        join_count = int((fine_rows if n == ns[-1] else old)
                         [("join", "H") if n == ns[-1] else
                          (member, n, "join", "H")]["count"])
        assert sum(int(r["count"]) for r in rows) == join_count
        profiles.extend(rows)
        for field, column in (("H", "H_rms_h4"), ("M", "M_rms_h4")):
            peak = max(rows, key=lambda r: float(r[column]))
            peaks.append({"member": member, "N": n, "field": field,
                          "bin": peak["bin"], "rho": peak["rho_mid"],
                          "r": peak["r_mean"], "z": peak["z_mean"],
                          "rms_h4": peak[column],
                          "max_h4": peak[f"{field}_max_h4"],
                          "relative_rms": peak[f"{field}_ratio_rms"]})
    for a, b in zip(ns[1:-1], ns[2:]):
        for field in ("H", "M"):
            column = f"{field}_rms_h4"
            x = [float(r[column]) for r in by_n[a]]
            y = [float(r[column]) for r in by_n[b]]
            overlays.append({"member": member, "field": field,
                             "N_coarse": a, "N_fine": b,
                             "relative_curve_L2": math.sqrt(
                                 sum((u-v)**2 for u, v in zip(x, y)) /
                                 sum(v*v for v in y))})
    for field in ("H", "M"):
        rows = by_n[ns[-1]]
        power = [int(r["count"]) * float(r[f"{field}_rms_h4"])**2
                 for r in rows]
        concentration.append({"member": member, "N": ns[-1],
                              "field": field, "bins_37_39_L2_power_fraction":
                              sum(power[37:]) / sum(power)})

write_csv("grid-results-four.csv", list(results[0]), results)
write_csv("grid-orders-four.csv", list(orders[0]), orders)
write_csv("join-profiles.csv", list(profiles[0]), profiles)
write_csv("join-localization.csv", list(peaks[0]), peaks)
write_csv("join-overlay.csv", list(overlays[0]), overlays)
write_csv("join-concentration.csv", list(concentration[0]), concentration)
stage_lines, resources = [], []
for member, ns in LADDERS.items():
    for n in ns[1:]:
        lines = (HERE / f"grid-{member}-{n}.log").read_text().splitlines()
        stage_lines.extend(lines)
        stage = {line.split()[0].split("=")[1]:
                 dict(part.split("=", 1) for part in line.split()[1:])
                 for line in lines}
        resources.append({"member": member, "N": n,
                          "min_lapse": 1e-9 if n == ns[-1] else 1e-8,
                          "min_chi": 1e-9 if n == ns[-1] else 1e-8,
                          "min_alpha_margin": stage["setter"]["min_alpha_margin"],
                          "min_chi_margin": stage["setter"]["min_chi_margin"],
                          "seconds": stage["report"]["elapsed"],
                          "peak_RSS_bytes": stage["report"]["peak_RSS_bytes"]})
(HERE / "grid-finer-stages.log").write_text("\n".join(stage_lines) + "\n")
write_csv("grid-finer-resources.csv", list(resources[0]), resources)
print(f"rows: results={len(results)} orders={len(orders)} "
      f"profiles={len(profiles)} peaks={len(peaks)} overlays={len(overlays)}")
print(f"maximum old-grid rerun relative norm difference={max_rerun_error:.3e}")
join_gate = all(float(r["L2_order_23"]) >= 3.7 for r in orders
                if r["mask"] == "join" and r["residual"] != "GaussB")
overlay_gate = all(
    next(r["relative_curve_L2"] for r in overlays if
         r["member"] == member and r["field"] == field and
         r["N_coarse"] == ns[2]) <
    next(r["relative_curve_L2"] for r in overlays if
         r["member"] == member and r["field"] == field and
         r["N_coarse"] == ns[1])
    for member, ns in LADDERS.items() for field in ("H", "M"))
print(f"join_L2_gate={join_gate} overlay_tightens={overlay_gate} "
      f"outcome={'PRE-ASYMPTOTIC' if join_gate and overlay_gate else 'UNCLASSIFIED'}")
