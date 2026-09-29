#!/usr/bin/env python3
"""Serial detached queue; each process owns its log and atomic exit marker."""
import json
import os
import resource
from pathlib import Path
import subprocess
import sys
import time


def run(plan):
    for job in json.loads(Path(plan).read_text()):
        directory = Path(job['directory'])
        directory.mkdir(parents=True, exist_ok=True)
        for name in ('chk', 'plt'):
            (directory / name).mkdir(exist_ok=True)
        if (directory / 'done.exit').exists():
            if (directory / 'done.exit').read_text().strip() != '0':
                return 1
            continue
        # Jobs are serial. Registered per-job allocation bounds must be <12 GB.
        assert job['memory_bound_GB'] < 12
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(job.get('threads', 4))
        env['CHOMBO_HOME'] = '/Users/auroradysis/Workspace/EMS-deps/Chombo/lib'
        start = time.time()
        with (directory / 'run.log').open('w') as log:
            result = subprocess.run(job['command'],
                                    cwd=directory, env=env, stdin=subprocess.DEVNULL,
                                    stdout=log, stderr=subprocess.STDOUT)
        (directory / 'wall_seconds').write_text(f'{time.time()-start:.6f}\n')
        # macOS reports bytes. In a serial queue this is a conservative bound
        # on each child: getrusage retains the largest previous child's peak.
        (directory / 'peak_rss_bytes').write_text(
            f'{resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}\n')
        marker = directory / 'done.exit.tmp'
        marker.write_text(f'{result.returncode}\n')
        marker.replace(directory / 'done.exit')
        if result.returncode:
            return result.returncode
    return 0


if __name__ == '__main__':
    plan = str(Path(sys.argv[-1]).resolve())
    if sys.argv[1] == '--detach':
        root = Path(plan).parent
        with (root / 'launcher.log').open('w') as log:
            process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), plan],
                                       stdin=subprocess.DEVNULL, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True,
                                       close_fds=True)
        (root / 'launcher.pid').write_text(f'{process.pid}\n')
        print(process.pid)
    else:
        sys.exit(run(plan))
