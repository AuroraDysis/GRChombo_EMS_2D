"""Syntax-check the actual production headers using the selected local compiler."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import time
import sys
from config import LOG_ROOT, CHOMBO_HOME, CXX

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'Tests/EMSNP/results-t2'
env=dict(os.environ,CHOMBO_HOME=CHOMBO_HOME)
dry=subprocess.run(['make','-n','-W','EMSBH2DLevel.cpp','-C','Examples/EMS','all','DIM=2','CXX='+CXX],cwd=ROOT,
                   env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
logs=LOG_ROOT/OUT.name;logs.mkdir(parents=True,exist_ok=True)
(logs/'production-dry-build.log').write_text(dry.stdout)
command=next(line for line in dry.stdout.splitlines()
             if line.startswith(CXX+' ') and '-c EMSBH2DLevel.cpp' in line)
args=shlex.split(command)
index=args.index('-o');del args[index:index+2];args.remove('-c')
args[1:1]=['-fsyntax-only']
start=time.monotonic()
result=subprocess.run(args,cwd=ROOT/'Examples/EMS',text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
report={'command':shlex.join(args),'seconds':time.monotonic()-start,
        'exit':result.returncode,'output':result.stdout,'workaround':None}
(OUT/'production-syntax.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(result.returncode)
