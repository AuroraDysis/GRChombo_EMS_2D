"""Sequential restart replay with peak resident-memory accounting.

The sandbox denies ps/sysctl. Each command records its peak RSS with getrusage;
the final peak must also be inspected against the operator's 6 GB limit.
"""
import json
from pathlib import Path
import resource
import subprocess
import time

out=Path(__file__).resolve().parent
start=time.monotonic()
with (out/'run-final.log').open('w') as log:
    process=subprocess.Popen(['zsh','Tests/EMSNP/run_t3.sh'],
        stdout=log,stderr=subprocess.STDOUT)
    status=process.wait()
    log.write(f'\nRUN_EXIT={status}\n')
(out/'resume-memory.json').write_text(json.dumps({
    'elapsed_seconds':time.monotonic()-start,'exit':status,
    'usage':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
    'live_rss_monitor':'unavailable: sandbox denies ps/sysctl; inspect MEASURE logs'},indent=2)+'\n')
raise SystemExit(status)
