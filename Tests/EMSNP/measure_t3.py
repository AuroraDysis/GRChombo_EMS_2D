"""Time a local control and report Darwin child peak RSS without sysctl/ps."""
import resource
import subprocess
import sys
import time

start=time.monotonic()
status=subprocess.run(sys.argv[1:]).returncode
peak=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
print(f'MEASURE seconds={time.monotonic()-start:.6f} peak_rss_bytes={peak}',file=sys.stderr)
if peak>6_000_000_000:
    print('FAIL: operator per-process memory budget exceeded',file=sys.stderr)
    status=1
raise SystemExit(status)
