#!/usr/bin/env python3
"""Serial T7 queue with atomic markers and the registered shared disk ceiling."""
import json,os,resource,subprocess,sys,time,tempfile,shutil
from pathlib import Path

def size(root):
    return sum(p.stat().st_size for p in root.rglob('*') if p.is_file())

def run(plan):
    jobs=json.loads(Path(plan).read_text());root=Path(jobs[0]['directory']).parent
    for job in jobs:
        d=Path(job['directory']);d.mkdir(exist_ok=True)
        assert job['memory_bound_GB']<12
        if (d/'done.exit').exists():
            if (d/'done.exit').read_text().strip()!='0':return 1
            continue
        start=time.monotonic();reason='';code=0
        with (d/'run.log').open('w') as log:
            child=subprocess.Popen(job['command'],cwd=d,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                                   env={**os.environ,'OMP_NUM_THREADS':'4','CHOMBO_HOME':'/Users/auroradysis/Workspace/EMS-deps/Chombo/lib'})
            while child.poll() is None:
                if size(root)>5_700_000_000:
                    reason='registered collective evolution-output ceiling 5.7 GB';child.terminate()
                    try:child.wait(timeout=10)
                    except subprocess.TimeoutExpired:child.kill();child.wait()
                    break
                time.sleep(10)
            code=child.wait()
            if reason:log.write('\nT7 OUTPUT GATE: '+reason+'\n');code=125
        (d/'wall_seconds').write_text(f'{time.monotonic()-start:.6f}\n')
        (d/'peak_rss_bytes').write_text(f'{resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}\n')
        p=d/'done.exit.tmp';p.write_text(f'{code}\n');p.replace(d/'done.exit')
        if code:return code
    return 0

def selfcheck():
    with tempfile.TemporaryDirectory(prefix='t7-queue-check-',dir='/private/tmp') as name:
        root=Path(name);d=root/'case';d.mkdir();p=root/'plan.json'
        p.write_text(json.dumps([dict(directory=str(d),memory_bound_GB=1,command=[sys.executable,'-c','pass'])]))
        assert run(p)==0 and (d/'done.exit').read_text().strip()=='0' and not (d/'done.exit.tmp').exists()
        assert run(p)==0
        (d/'done.exit').unlink()
        p.write_text(json.dumps([dict(directory=str(d),memory_bound_GB=1,command=[sys.executable,'-c','import time; time.sleep(2)'])]))
        globals()['size']=lambda _:5_700_000_001
        assert run(p)==125 and (d/'done.exit').read_text().strip()=='125'
    print('PASS: atomic completion, idempotent success, shared output gate')

if __name__=='__main__':
    if sys.argv[1]=='--selfcheck':selfcheck();sys.exit(0)
    assert shutil.which('xz'), 'T7 diagnostic output requires installed xz'
    plan=str(Path(sys.argv[-1]).resolve())
    if sys.argv[1]=='--detach':
        root=Path(plan).parent
        with (root/'launcher.log').open('w') as log:
            child=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),plan],stdin=subprocess.DEVNULL,
                                   stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
        (root/'launcher.pid').write_text(f'{child.pid}\n');print(child.pid)
    else:sys.exit(run(plan))
