"""Check the independent writer cadence on the six exp-0012 schedules."""
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FILES = ("ref-M32", "ref-M64", "ref-M128", "ref-M64-k1", "B-M512", "B-M1024")


def parameters(name):
    values = {}
    for line in (ROOT / f"Examples/EMS/params-cluster-{name}.txt").read_text().splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def crosses(time, dt, interval):
    eps = 1e-12 * max(interval, dt)
    return math.floor((time + eps) / interval) > math.floor(
        (time - dt + eps) / interval
    )


def main():
    writer = csv.writer(sys.stdout)
    writer.writerow(("case", "finest_dt", "final_time_over_M", "stop_time_over_M",
                     "radial_boundaries_reached", "new_captured", "new_missing",
                     "old_dropped_nominal", "epsilon_window_old_drops",
                     "missing_if_one_coarse_step_short"))
    predicted = {"ref-M128": (8.,), "B-M512": (3., 3.5, 4., 4.5, 5., 9.5),
                 "B-M1024": (1.5, 2., 2.5, 5.)}
    for name in FILES:
        p = parameters(name)
        mass = float(p["bh_mass"])
        levels = int(p["max_level"])
        ratio = 2 ** levels
        coarse_dt = float(p["dt_multiplier"]) * float(p["L"]) / int(p["N1"])
        dt = coarse_dt / ratio
        summary = float(p["reference_diagnostics_interval"])
        radial = float(p["reference_radial_interval"])
        stop = float(p["stop_time"])
        steps = min(int(p["max_steps"]), math.ceil(stop / coarse_dt))
        reached = set()
        old = set()
        for coarse in range(steps):
            for fine in range(1, ratio + 1):
                time = coarse * coarse_dt + fine * dt
                if crosses(time, dt, radial):
                    boundary = math.floor((time + 1e-12 * max(radial, dt)) / radial)
                    reached.add(boundary)
                    if crosses(time, dt, summary):
                        old.add(boundary)
        final = steps * coarse_dt
        expected = set(range(1, math.floor((min(final, stop) +
                           1e-12 * radial) / radial) + 1))
        assert reached == expected, (name, sorted(expected - reached))
        vulnerable = 0
        for over_mass in predicted.get(name, ()):
            boundary = round(over_mass * mass / radial) * radial
            radial_eps = 1e-12 * max(radial, dt)
            summary_eps = 1e-12 * max(summary, dt)
            near = boundary - (radial_eps + summary_eps) / 2
            assert crosses(near, dt, radial) and not crosses(near, dt, summary)
            vulnerable += 1
        short_final = math.floor(((steps - 1) * coarse_dt +
                                  1e-12 * radial) / radial)
        missing_if_short = len(expected) - short_final
        assert missing_if_short >= 0
        writer.writerow((name, f"{dt:.17g}", f"{final / mass:.17g}",
                         f"{stop / mass:.17g}", len(expected), len(reached), 0,
                         len(reached - old), vulnerable, missing_if_short))


if __name__ == "__main__":
    main()
