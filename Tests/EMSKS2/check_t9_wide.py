"""Check one emitted wide profile for complete, finite composite-AMR bins."""
import csv
import math
import sys

rows = list(csv.DictReader(open(sys.argv[1], newline="")))
assert len(rows) == 80
assert [int(row["bin"]) for row in rows] == list(range(80))
assert math.isclose(float(rows[0]["rho_lo_over_rh"]), 1.0)
assert float(rows[-1]["rho_hi_over_rh"]) >= float(sys.argv[2])
for previous, row in zip(rows, rows[1:]):
    assert math.isclose(float(previous["rho_hi_over_rh"]),
                        float(row["rho_lo_over_rh"]), rel_tol=1e-14)
for row in rows:
    assert int(row["count"]) > 0
    assert int(row["finest_level"]) >= 0
    assert all(math.isfinite(float(row[key])) for key in
               ("H_RMS", "M_RMS", "GaussE_RMS", "K_minus_Kstar_RMS",
                "alpha_drift_RMS"))
print(f"wide profile: {len(rows)} finite populated bins through {sys.argv[2]} r_h")
