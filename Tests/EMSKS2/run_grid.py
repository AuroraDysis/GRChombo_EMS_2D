#!/usr/bin/env python3
"""Run the registered EMSKS 2 grid ladder and write its four result files."""

import csv
import math
import os
from pathlib import Path
import subprocess
import sys


def main(profile_root):
    here = Path(__file__).resolve().parent
    exe = here / "EMSKS2GridConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex"
    rows, audits, stages = [], [], []
    for member, levels in (("B", (512, 1024, 2048)),
                           ("E", (512, 1024, 2048)),
                           ("reference", (32, 64, 128))):
        for n in levels:
            run = subprocess.run(
                [str(exe), member, str(n), str(profile_root / f"{member}.ks2")],
                env={**os.environ, "OMP_NUM_THREADS": "1"},
                text=True, capture_output=True, timeout=240, check=True)
            stages.extend(run.stderr.splitlines())
            print(run.stderr, file=sys.stderr, end="", flush=True)
            for line in run.stdout.splitlines():
                if line.startswith("audit,"):
                    z = line.split(",")
                    audits.append([z[1], z[2], *z[4::2]])
                elif line.startswith("member,"):
                    header = line
                else:
                    rows.append(line)
    (here / "grid-results.csv").write_text(header + "\n" + "\n".join(rows) + "\n")
    (here / "grid-stages.log").write_text("\n".join(stages) + "\n")
    with (here / "puncture-audit.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(("member", "N", "rho", "alpha", "alpha_over_rho_nu",
                         "chi", "chi_power_ratio", "E2", "E2_limit_ratio"))
        writer.writerows(audits)
    by_key = {}
    for row in csv.DictReader((here / "grid-results.csv").open()):
        by_key.setdefault((row["member"], row["mask"], row["residual"]), {})[
            int(row["N"])] = row
    with (here / "grid-orders.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(("member", "mask", "residual", "N_coarse", "N_middle",
                         "N_fine", "L2_coarse", "L2_middle", "L2_fine",
                         "L2_order_coarse_middle", "L2_order_middle_fine",
                         "Linf_coarse", "Linf_middle", "Linf_fine",
                         "Linf_order_coarse_middle", "Linf_order_middle_fine"))
        for key, by_n in by_key.items():
            ns = sorted(by_n)
            values = [[float(by_n[n][field]) for n in ns]
                      for field in ("L2", "Linf")]
            def order(a, b):
                return math.log2(a / b) if a > 0 and b > 0 else ""
            writer.writerow([*key, *ns, *values[0],
                             order(values[0][0], values[0][1]),
                             order(values[0][1], values[0][2]), *values[1],
                             order(values[1][0], values[1][1]),
                             order(values[1][1], values[1][2])])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_grid.py path/to/t3")
    main(Path(sys.argv[1]))
