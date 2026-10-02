#!/usr/bin/env python3
"""T18 uses the existing measured runner; local jobs are serial and bounded."""
import os
import runpy
from pathlib import Path

os.environ.update(T13_OUTPUT_ROOT='/private/tmp/ems-t18',
                  T13_PROCESS_CAP_BYTES='4000000000',
                  T13_TREE_CAP_BYTES='4500000000',
                  T13_DISK_CAP_BYTES='12000000000')
runpy.run_path(str(Path(__file__).with_name('t13-run.py')), run_name='__main__')
