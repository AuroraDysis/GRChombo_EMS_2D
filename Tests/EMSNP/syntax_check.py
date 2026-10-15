"""Production syntax check with the operator's temporary access workaround."""
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import time
import sys

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'Tests/EMSNP/results-t2'
env=dict(os.environ,CHOMBO_HOME='/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
dry=subprocess.run(['make','-n','-C','Examples/EMS','all','DIM=2'],cwd=ROOT,
                   env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
(OUT/'production-dry-build.log').write_text(dry.stdout)
command=next(line for line in dry.stdout.splitlines()
             if line.startswith('/opt/homebrew/bin/g++-16 ') and '-c EMSBH2DLevel.cpp' in line)
args=shlex.split(command)
index=args.index('-o');del args[index:index+2];args.remove('-c')
with tempfile.TemporaryDirectory(prefix='ems-np-headers-',dir='/private/tmp') as directory:
    for source in (ROOT/'Source/GRChomboCore').glob('*.hpp'):shutil.copy2(source,directory)
    header=Path(directory)/'GRAMRLevel.hpp';lines=header.read_text().splitlines(keepends=True)
    assert lines[50].strip()=='private:'
    lines[50]=lines[50].replace('private:','protected:');header.write_text(''.join(lines))
    args[1:1]=['-I'+directory,'-fsyntax-only']
    start=time.monotonic()
    result=subprocess.run(args,cwd=ROOT/'Examples/EMS',text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    report={'command':shlex.join(args),'seconds':time.monotonic()-start,
            'exit':result.returncode,'output':result.stdout,'workaround':'temporary GRAMRLevel.hpp line 51 private -> protected'}
(OUT/'production-syntax.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(result.returncode)
