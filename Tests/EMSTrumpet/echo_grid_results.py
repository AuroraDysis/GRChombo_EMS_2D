"""Summarize B grid logs; inputs are the N=128,256,512 stdout files."""
import csv
import math
import re
import sys

fields = ("H", "Mx", "My", "GaussE", "GaussB")
runs = {}
for path in sys.argv[1:]:
    lines = open(path).read().splitlines()
    n = int(re.search(r"\bN=(\d+)", next(x for x in lines if x.startswith("N="))).group(1))
    header = next(i for i, x in enumerate(lines) if x.startswith("region,count,"))
    row = next(x for x in csv.DictReader(lines[header:]) if x["region"] == "full")
    assert n not in runs and int(row["chi_floor"]) == int(row["lapse_floor"]) == 0
    runs[n] = row
assert sorted(runs) == [128, 256, 512]
names = ["N", "dx", "cells", "count", "field_error"]
for field in fields:
    names += [f"{field}_L2", f"{field}_pL2", f"{field}_Linf", f"{field}_pLinf"]
writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(names)
for n in sorted(runs):
    row, prev = runs[n], runs.get(n // 2)
    out = [n, 1 / n, n * n // 2, row["count"], row["field_error"]]
    for field in fields:
        for norm in ("L2", "Linf"):
            key = f"{field}_{norm}"
            value = float(row[key])
            old = float(prev[key]) if prev else 0
            out += [row[key], f"{math.log2(old / value):.8g}" if old > 0 and value > 0 else ""]
    writer.writerow(out)
