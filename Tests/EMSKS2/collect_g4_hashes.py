"""Collect the unchanged per-component G4 checkpoint comparisons."""
import contextlib
import csv
import io
import sys
from pathlib import Path

from compare_state import compare


root = Path(sys.argv[1])
writer = csv.writer(sys.stdout)
for step in (0, 4):
    for case in ("trumpet-reference", "trumpet-echo-B", "EMSCTT", "legacy_dat"):
        name = f"EMS_{step:06d}.2d.hdf5"
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            compare(root / case / "baseline" / "chk" / name,
                    root / case / "current" / "chk" / name, case)
        rows = csv.reader(io.StringIO(buffer.getvalue()))
        if step == 0 and case == "trumpet-reference":
            writer.writerow(next(rows))
        else:
            next(rows)
        writer.writerows(rows)
