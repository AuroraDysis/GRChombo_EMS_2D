"""Run the four matched 34f2ef0/current evolution controls sequentially."""
from pathlib import Path
import argparse
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
CONTROLS = ROOT / "Tests/EMSKS2/g4-controls"
EXE = "Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex"
CASES = ("trumpet-reference", "trumpet-echo-B", "EMSCTT", "legacy_dat")

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=CONTROLS)
parser.add_argument("--baseline", type=Path,
                    default=Path("/private/tmp/emsks2-t5b-baseline/Examples/EMS") / EXE)
parser.add_argument("--current", type=Path, default=ROOT / "Examples/EMS" / EXE)
args = parser.parse_args()
executables = {"baseline": args.baseline, "current": args.current}

for index, (case, branch) in enumerate(
    ((case, branch) for case in CASES for branch in executables), 1
):
    directory = args.output / case / branch
    if args.output != CONTROLS:
        (directory / "chk").mkdir(parents=True, exist_ok=True)
        shutil.copy2(CONTROLS / case / branch / "params.txt", directory / "params.txt")
    checkpoint = directory / "chk/EMS_000004.2d.hdf5"
    start = time.monotonic()
    with (directory / "run.log").open("w") as log:
        result = subprocess.run(
            ["/usr/bin/time", "-p", "/opt/homebrew/bin/gtimeout", "-k", "10s",
             "600s", "env", "OMP_NUM_THREADS=1", str(executables[branch]),
             "params.txt"],
            cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=False,
        )
    elapsed = time.monotonic() - start
    print(f"stage={case}/{branch} elapsed={elapsed:.3f} "
          f"items={index}/8 exit={result.returncode}", flush=True)
    if result.returncode or not checkpoint.exists():
        raise SystemExit(f"{case}/{branch}: missing checkpoint or failed run")
