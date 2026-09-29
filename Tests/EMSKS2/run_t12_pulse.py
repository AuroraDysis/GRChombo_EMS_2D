"""Generate and run the T12 two-level pulse comparisons (one process at a time)."""
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

root, executable = Path(sys.argv[1]), sys.argv[2]
base = (Path(__file__).with_name("params-t12-static-raw.txt")
        .read_text().splitlines())
for N in (256, 384, 512):
    for transfer in ("raw", "relative"):
        for coverage, radius in (("narrow", 2), ("wide", 5)):
            directory = root / f"N{N}" / transfer / coverage
            (directory / "chk").mkdir(parents=True, exist_ok=True)
            steps = 60 * N // 256
            changes = {
                "N1": str(N), "N2": str(N // 2),
                "max_level": "1", "regrid_interval": "0",
                "num_mass_extraction_radii": "1",
                "mass_extraction_levels": "1",
                "mass_extraction_radii": str(radius),
                "max_steps": str(steps), "stop_time": "10",
                "checkpoint_interval": str(steps), "dt_multiplier": "0.175",
                "reference_transfer": transfer,
            }
            lines = [f'{line.split("=", 1)[0].strip()} = {changes[line.split("=", 1)[0].strip()]}'
                     if "=" in line and line.split("=", 1)[0].strip() in changes
                     else line for line in base
                     if not line.startswith("reference_transfer_probe_path")]
            lines += ["reference_transfer_pulse_amplitude = 0.0001",
                      "reference_transfer_pulse_x = 2.3",
                      "reference_transfer_pulse_width = 0.35"]
            (directory / "params.txt").write_text("\n".join(lines) + "\n")
            start = time.perf_counter()
            with (directory / "run.log").open("w") as log:
                result = subprocess.run([executable, "params.txt"], cwd=directory,
                                        env={**os.environ, "OMP_NUM_THREADS": "1"},
                                        stdout=log, stderr=subprocess.STDOUT,
                                        timeout=300)
            print(f"stage=N{N}-{transfer}-{coverage} "
                  f"elapsed={time.perf_counter() - start:.3f}s "
                  f"peak_RSS_bytes_so_far={resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss} "
                  f"exit={result.returncode} items={steps}/{steps}", flush=True)
            if result.returncode:
                raise SystemExit(result.returncode)
