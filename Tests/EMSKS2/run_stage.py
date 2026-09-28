"""Run one serial evolution stage with a hard cap and record wall/RSS."""
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

directory, name, done, total = sys.argv[1:5]
executable = Path(__file__).resolve().parents[2] / "Examples/EMS" / (
    "Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex")
start = time.monotonic()
with (Path(directory) / "run.log").open("w") as log:
    result = subprocess.run(
        ["/opt/homebrew/bin/gtimeout", "-k", "10s", "600s", str(executable),
         "params.txt"], cwd=directory, stdout=log, stderr=subprocess.STDOUT,
        env={**os.environ, "OMP_NUM_THREADS": "1"}, check=False)
elapsed = time.monotonic() - start
print(f"stage={name} elapsed={elapsed:.3f} "
      f"peak_RSS_bytes={resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss} "
      f"items={done}/{total} exit={result.returncode}", flush=True)
if result.returncode:
    raise SystemExit(result.returncode)
