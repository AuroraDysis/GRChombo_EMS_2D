#!/usr/bin/env python3
"""Collect the production RHS ladder without subtracting a reference RHS."""

import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
LADDERS = {"B": (512, 1024, 2048, 4096),
           "E": (512, 1024, 2048, 4096),
           "reference": (32, 64, 128, 256)}


def order(coarse, fine):
    return math.log2(coarse / fine) if coarse > 0 and fine > 0 else ""


def write_csv(name, rows, fields):
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    results = []
    for member, levels in LADDERS.items():
        for n in levels:
            with (HERE / f"rhs-{member}-{n}.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            assert len(rows) == 2 * 4 * 14
            for row in rows:
                assert row["member"] == member and int(row["N"]) == n
                assert int(row["chi_floor"]) == int(row["lapse_floor"]) == 0
                assert min(float(row["min_alpha_margin"]),
                           float(row["min_chi_margin"])) >= 100
                results.append({**row, "source": "production_RHS"})
                if row["group"] == "beta":
                    assert float(row["L2"]) == float(row["Linf"]) == 0
                    results.append({**row, "group": "B_gauge",
                                    "source": "exact_zero_assertion"})
    fields = list(results[0])
    write_csv("rhs-results.csv", results, fields)
    grouped = {}
    for row in results:
        grouped.setdefault((row["member"], row["KO"], row["mask"],
                            row["group"]), {})[int(row["N"])] = row
    orders = []
    for (member, ko, mask, group), by_n in grouped.items():
        ns = LADDERS[member]
        for c, f in zip(ns[:-1], ns[1:]):
            a, b = by_n[c], by_n[f]
            orders.append({"member": member, "KO": ko, "mask": mask,
                           "group": group, "N_coarse": c, "N_fine": f,
                           "L2_coarse": a["L2"], "L2_fine": b["L2"],
                           "L2_order": order(float(a["L2"]), float(b["L2"])),
                           "Linf_coarse": a["Linf"], "Linf_fine": b["Linf"],
                           "Linf_order": order(float(a["Linf"]), float(b["Linf"]))})
    write_csv("rhs-orders.csv", orders, list(orders[0]))
    failures = [r for r in orders if r["N_fine"] == LADDERS[r["member"]][-1]
                and r["L2_order"] != "" and r["L2_order"] < 3.5]
    print(f"results={len(results)} orders={len(orders)} finest_L2_below_3.5="
          f"{len(failures)} outcome=NOT-STATIONARY")
    for member in LADDERS:
        join = [r for r in orders if r["member"] == member and
                r["KO"] == "1" and r["mask"] == "join" and
                r["N_fine"] == LADDERS[member][-1] and r["L2_order"] != ""]
        print(member, "join_KO_on_min_order", min(r["L2_order"] for r in join))


if __name__ == "__main__":
    main()
