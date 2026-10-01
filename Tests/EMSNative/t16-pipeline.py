#!/usr/bin/env python3
"""Serial detached launch and bounded current-state analysis, with done markers."""
import sys
sys.dont_write_bytecode=True
import subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t16');PY='/Users/auroradysis/miniconda3/bin/python'
def measured(label,directory,command):
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/(label+'.log')).open('wb') as log:
        subprocess.run([PY,str(HERE/'t16-run.py'),'--measure',label,'--directory',str(directory),'--',*command],
            stdout=log,stderr=subprocess.STDOUT,check=True)
measured('evolution-worker',ROOT/'evolution',[PY,str(HERE/'t16-run.py'),'--worker',str(ROOT/'evolution/plan.json')])
measured('analysis',ROOT/'analysis',[PY,str(HERE/'t16-analyze.py'),'analyze'])
measured('report',ROOT/'report',[PY,str(HERE/'t16-report.py')])
# Pin the report as well as the scientific tables produced before it.
measured('manifest',ROOT/'manifest',[PY,str(HERE/'t16-analyze.py'),'resources'])
print('T16 isolated pipeline finished; binary awaits certified tranche-2a companion.',flush=True)
