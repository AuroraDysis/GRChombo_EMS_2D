#!/usr/bin/env python3
"""Small real-child check: isolated disk gate and RSS after gate termination."""
import json, os, subprocess, sys, tempfile
from pathlib import Path

here=Path(__file__).resolve().parent
root=Path('/private/tmp/ems-t14/runner-check');root.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(dir=root) as temporary:
    directory=Path(temporary)
    env=dict(os.environ,T14_OUTPUT_ROOT=temporary,T14_DISK_CAP_BYTES='1024')
    child="from pathlib import Path; import time; a=bytearray(48000000); Path('payload').write_bytes(b'x'*2048); time.sleep(30)"
    result=subprocess.run([sys.executable,str(here/'t14-run.py'),'--measure','gate-check',
        '--directory',temporary,'--',sys.executable,'-c',child],env=env,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    record=json.loads((directory/'gate-check.resources.json').read_text())
    assert result.returncode==125, result.stdout.decode()
    assert record['gate_reason']=='output ceiling 1024 bytes'
    assert record['cap_GB']==3 and record['peak_rss_bytes']>=48000000
    assert record['peak_rss_bytes']>=max(v[0] for v in record['per_pid_peaks'].values())
    (root/'verified.json').write_text(json.dumps(record,indent=2)+'\n')
print('PASS: isolated disk gate, fixed 3 GB cap, sampled RSS retained after SIGTERM')
