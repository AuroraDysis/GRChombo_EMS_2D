"""Existing t=0 grid regression matrix, serial, each executable capped at 240 s."""
from pathlib import Path
import json
import os
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'results/regressions'
EMS=Path('/Users/auroradysis/Workspace/EMS')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    executable=lambda folder,name: str(next((ROOT/'Tests'/folder).glob(name+'2d.*.ex')))
    trumpet=executable('EMSTrumpet','EMSTrumpetGridConvergence')
    binary=executable('EMSTrumpet','EMSTrumpetBinaryGridConvergence')
    ctt=executable('EMSCTT','EMSCTTGridConvergence')
    profile=str(ROOT/'Tests/EMSTrumpet/fixtures/reference.trumpet')
    companion=str(EMS/'.data/binary-ctt-t6/n28-r6.ctt')
    jobs=[]
    for eta in ('0','0.20273255','-0.20273255'):
        for n in (32,64,128):
            jobs.append((f'single-{eta}-{n}',[trumpet,str(n),eta,profile]))
    for n in (128,256,512):
        jobs.append((f'echo-B-{n}',[trumpet,str(n),'0',str(ROOT/'Tests/EMSTrumpet/fixtures-echo/B/reference.trumpet'),'echo']))
    for mode in ('density','old'):
        for n in (32,64,128):
            jobs.append((f'binary-{mode}-{n}',[binary,str(n),'0.20273255',profile,mode]))
    for mode in ('density','ctt'):
        for refinement in ('0.5','1','2','4'):
            jobs.append((f'{mode}-{refinement}',[ctt,refinement,profile,companion if mode=='ctt' else 'density']))
    records=[]
    for name,cmd in jobs:
        started=time.monotonic()
        with (OUT/(name+'.log')).open('w') as log:
            try:
                status=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,
                                      env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=240).returncode
            except subprocess.TimeoutExpired:
                status=124
        record=dict(name=name,exit=status,seconds=time.monotonic()-started)
        records.append(record);print(json.dumps(record),flush=True)
        (OUT/'status.json').write_text(json.dumps(records,indent=2)+'\n')
    if any(r['exit'] for r in records): raise SystemExit(1)
    (OUT/'complete.json').write_text(json.dumps(dict(count=len(records),passed=True))+'\n')


if __name__=='__main__': main()
