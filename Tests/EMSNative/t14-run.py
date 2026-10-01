#!/usr/bin/env python3
"""T14 uses the existing measured serial runner, with stricter gates."""
import csv, json, os, runpy, sys
from pathlib import Path

disk_cap=int(os.environ.get('T14_DISK_CAP_BYTES','5700000000'))
assert 0<disk_cap<=8000000000, 'T14 output ceiling must be <=8 GB'
os.environ.update(T13_OUTPUT_ROOT=os.environ.get('T14_OUTPUT_ROOT','/private/tmp/ems-t14'),
                  T13_PROCESS_CAP_BYTES='3000000000',
                  T13_TREE_CAP_BYTES='3000000000',
                  T13_DISK_CAP_BYTES=str(disk_cap))
if '--detach' in sys.argv:
    root=Path('/private/tmp/ems-t14')
    assert json.loads((root/'census/verified.json').read_text())['old_boxes_identical']
    rows=list(csv.DictReader(Path(__file__).with_name('t14-controls.csv').open()))
    assert len(rows)==16 and all(int(a['bit_mismatches'])==0 for a in rows)
runpy.run_path(str(Path(__file__).with_name('t13-run.py')),run_name='__main__')
