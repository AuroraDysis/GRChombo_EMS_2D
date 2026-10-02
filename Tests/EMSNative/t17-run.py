#!/usr/bin/env python3
"""Reuse the serial native watchdog with T17's resource limits."""
import os, runpy
from pathlib import Path
os.environ.update(T13_OUTPUT_ROOT='/private/tmp/ems-t17',
                  T13_PROCESS_CAP_BYTES='6000000000',
                  T13_TREE_CAP_BYTES='6500000000',
                  T13_DISK_CAP_BYTES='40000000000')
runpy.run_path(str(Path(__file__).with_name('t13-run.py')),run_name='__main__')
