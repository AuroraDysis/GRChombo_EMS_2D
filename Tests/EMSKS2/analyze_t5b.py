"""Assemble the registered T5b comparisons from raw run diagnostics."""
import csv
import hashlib
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
JL = Path("/Users/auroradysis/Workspace/EMS")


def read(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def write(name, rows):
    rows = list(rows)
    assert rows
    with (ROOT / name).open("w", newline="") as target:
        writer = csv.DictWriter(target, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(name, len(rows))


def at(rows, mask, time):
    found = [row for row in rows if row["mask"] == mask and
             abs(float(row["time_over_M"]) - time) < 1e-8]
    assert len(found) == 1, (mask, time, len(found))
    return found[0]


s1 = {n: read(ROOT / f"s1-B-{n}-pilot/diagnostics.csv") for n in (256, 512)}
rhs512 = read(ROOT / "rhs-B-512.csv")
for mask in ("join", "KS_collar", "exterior", "far"):
    expected_count = next(int(row["count"]) for row in rhs512 if
                          row["KO"] == "1" and row["group"] == "K" and
                          row["mask"] == mask)
    assert int(s1[512][("join", "KS_collar", "exterior", "far").index(mask)]["count"]) == expected_count
write("s1-time-series.csv", (
    {"M_over_h": n, **row} for n, rows in s1.items() for row in rows))
orders = []
for time in (1 / 256, 1 / 128):
    for mask in ("join", "KS_collar", "exterior", "far"):
        coarse, fine = (at(s1[n], mask, time) for n in (256, 512))
        for field in ("H_L2", "M_L2", "GaussE_L2"):
            a, b = float(coarse[field]), float(fine[field])
            orders.append({"time_over_M": time, "mask": mask, "field": field,
                           "M256": a, "M512": b,
                           "first_pair_order": math.log2(a / b) if a and b else ""})
write("s1-first-pair-orders.csv", orders)

sources = (("harmonic", 32, "s2-reference-32-stage48"),
           ("harmonic", 64, "s2-reference-64-stage32"),
           ("onepluslog", 32, "s3-reference-32-onepluslog-stage48"))
series = {(gauge, n): read(ROOT / directory / "diagnostics.csv")
          for gauge, n, directory in sources}
for rows in (*s1.values(), *series.values()):
    assert all(row["nonfinite"] == "0" and row["floor_count"] == "0"
               for row in rows)
write("s2-s3-time-series.csv", (
    {"reference_f": gauge, "M_over_h": n, **row}
    for (gauge, n), rows in series.items() for row in rows))
ratios = []
for step in range(1, 17):
    time = step / 32
    for mask in ("join", "KS_collar", "exterior", "far"):
        a = at(series["harmonic", 32], mask, time)
        b = at(series["harmonic", 64], mask, time)
        for field in ("H_L2", "M_L2", "GaussE_L2", "max_alpha_drift"):
            av, bv = float(a[field]), float(b[field])
            ratios.append({"time_over_M": time, "mask": mask, "field": field,
                           "M32": av, "M64": bv,
                           "M32_over_M64": av / bv if bv else ""})
write("s2-common-ratios.csv", ratios)
snapshots = []
for gauge, n, directory in sources:
    rows = read(ROOT / directory / "diagnostics.csv.radial.csv")
    for time in ((.5, 1.5) if n == 32 else (.5,)):
        selected = [row for row in rows if
                    abs(float(row["time_over_M"]) - time) < 1e-8]
        assert len(selected) == 60
        snapshots.extend({"reference_f": gauge, "M_over_h": n, **row}
                         for row in selected)
write("s2-s3-radial-snapshots.csv", snapshots)
write("s1-radial-common.csv", (
    {"M_over_h": n, **row}
    for n in (256, 512)
    for row in read(ROOT / f"s1-B-{n}-pilot/diagnostics.csv.radial.csv")
    if abs(float(row["time_over_M"]) - 1 / 128) < 1e-8))

s3 = []
for time in (.5, 1., 1.5):
    for mask, field in (("exterior", "H_L2"), ("join", "max_alpha_drift")):
        h = float(at(series["harmonic", 32], mask, time)[field])
        o = float(at(series["onepluslog", 32], mask, time)[field])
        s3.append({"time_over_M": time, "mask": mask, "field": field,
                   "harmonic": h, "onepluslog": o,
                   "onepluslog_over_harmonic": o / h})
write("s3-comparison.csv", s3)

run_dirs = [path.parent for path in ROOT.glob("g4-controls/*/*/run.log")]
run_dirs += [path.parent for path in ROOT.glob("g4-restart/*/run.log")]
run_dirs += [path.parent for path in ROOT.glob("s[123]-*/run.log")]
run_dirs += [ROOT / "g4-ks2-experimental", ROOT / "g4-ks2-t4"]
defaults = set()
pattern = re.compile(r"Parameter: (\S+) not found.*default value = (.*)\.")
for directory in run_dirs:
    log = directory / "run.log"
    if not log.exists():
        continue
    for line in log.read_text(errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            defaults.add((str(directory.relative_to(ROOT)), *match.groups()))
write("parameter-defaults.csv", (
    {"run": run, "parameter": parameter, "default_value": value}
    for run, parameter, value in sorted(defaults)))

costs = []
for path in sorted(ROOT.glob("s[123]-*cost.log")):
    for line in path.read_text().splitlines():
        match = re.search(r"stage=(\S+) elapsed=([\d.]+) peak_RSS_bytes=(\d+) .*exit=(\d+)", line)
        if match:
            stage, wall, rss, code = match.groups()
            costs.append({"stage": stage, "wall_seconds": wall,
                          "peak_RSS_bytes": rss, "exit": code})
write("t5b-cost.csv", costs)
total_ram = 25769803776
free_fraction = .52  # memory_pressure -Q immediately after the M/512 pilot
projected = 4 * 2732638208  # measured M/512 peak; same boxes at half h
write("s1-resource-decision.csv", [{
    "M512_two_coarse_wall_seconds": 76.290,
    "M512_two_coarse_evolution_seconds": .546247 * 60,
    "M1024_projected_two_coarse_seconds":
        4 * (76.290 - .546247 * 60) + 8 * (.546247 * 60),
    "M1024_projected_three_coarse_seconds":
        4 * (76.290 - .546247 * 60) + 12 * (.546247 * 60),
    "M512_peak_RSS_bytes": 2732638208,
    "M1024_projected_RSS_bytes": projected,
    "total_RAM_bytes": total_ram,
    "free_fraction_before_dispatch": free_fraction,
    "max_dispatch_RSS_with_20pct_reserve_bytes":
        int((free_fraction - .2) * total_ram),
    "decision": "BLOCKED: projected peak exceeds free memory after reserve",
}])

expected = {
    JL / "docs/ks-puncture-format.md": "82947f2ceaadc9330040898c89c17c7e1c9c86f98c6ea8d93d0c675dfacb5914",
    JL / "artifacts/echo-evolution/t3/B.ks2": "6fe2e42892c3842a78d15ed8b1cdaccacf75d545aad8917085b8049a55a17d01",
    JL / "artifacts/echo-evolution/t3/E.ks2": "b2faa29e1c81cf0c871e6a4a6c8eab39d0ab01f26f0a6c0675534b8895768912",
    JL / "artifacts/echo-evolution/t3/reference.ks2": "fb46ceecd7565c5afb5b04b4ace4dd0790a832844d3f57a28333700230870533",
}
fingerprints = []
for path, supplied in expected.items():
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == supplied, path
    fingerprints.append({"path": str(path), "sha256": digest, "status": "verified"})
for name in ("smoke-B-harmonic.csv", "smoke-cost.csv", "g4-evolution-hashes.csv",
             "g4-restart-hashes.csv", "g4-ks2-init-hashes.csv", "s1-time-series.csv",
             "s2-s3-time-series.csv", "s2-s3-radial-snapshots.csv",
             "s2-common-ratios.csv", "s3-comparison.csv",
             "s1-resource-decision.csv", "s2-coverage.csv"):
    path = ROOT / name
    fingerprints.append({"path": str(path.relative_to(ROOT)),
                         "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "status": "output"})
write("t5b-fingerprints.csv", fingerprints)
