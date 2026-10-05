"""Local frozen-checkpoint Newton/flow comparison; no evolution or remote writes."""
import argparse, hashlib, json, os, shutil, subprocess, time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('solver', choices=['flow', 'newton'])
p.add_argument('steps', nargs='*', type=int, default=[4, 34, 40, 54, 72])
p.add_argument('--executable', required=True)
p.add_argument('--label')
p.add_argument('--set', action='append', default=[])
p.add_argument('--output', default='/private/tmp/rh-newton-run')
p.add_argument('--data', default='/Users/auroradysis/Workspace/EMS/.data/exp-0027')
a = p.parse_args()
root = Path(a.output)
root.mkdir(parents=True, exist_ok=True)
data = Path(a.data)
for step in a.steps:
    chk = data/'perf'/f'EMS_{step:06d}.2d.hdf5'
    if not chk.exists():
        print(f'{step}: checkpoint absent/incomplete', flush=True)
        continue
    receipt = json.loads((data/'collect/X/readouts'/f'{step:06d}'/'readout.receipt.json').read_text())
    digest = hashlib.file_digest(chk.open('rb'), 'sha256').hexdigest()
    assert digest == receipt['checkpoint_sha256'], (step, digest, receipt['checkpoint_sha256'])
    attempt = receipt['horizon']['attempts'][0]
    seed = data/'collect/X/readouts'/f'{step:06d}'/'horizon'/Path(attempt['receipt']).parent.name
    directory = root/f'{step:06d}-{a.label or a.solver}'
    directory.mkdir()
    for name in ('params.txt', 'rh_surf_0.dat', 'rh_f0.dat'):
        shutil.copy2(seed/name, directory/name)
    seed_sha256 = {name: hashlib.file_digest((seed/name).open('rb'), 'sha256').hexdigest()
                   for name in ('params.txt', 'rh_surf_0.dat', 'rh_f0.dat')}
    params = directory/'params.txt'
    lines = params.read_text().splitlines()
    assert sum(x.startswith('restart_file = ') for x in lines) == 1
    lines = [f'restart_file = {chk}' if x.startswith('restart_file = ') else x for x in lines]
    if a.solver == 'newton':
        lines += ['offline_solver = newton', 'offline_check_jacobian = true']
    lines += a.set
    params.write_text('\n'.join(lines)+'\n')
    exe = Path(a.executable).resolve()
    start = time.monotonic()
    with (directory/'run.log').open('w') as log:
        r = subprocess.run([str(exe), 'params.txt'], cwd=directory,
                           env=dict(os.environ, OMP_NUM_THREADS='1', CH_TIMER='1'),
                           stdout=log, stderr=subprocess.STDOUT)
    result = dict(step=step, solver=a.solver, executable=str(exe), seed=str(seed),
                  checkpoint_sha256=digest, elapsed_s=time.monotonic()-start,
                  returncode=r.returncode,
                  executable_sha256=hashlib.file_digest(exe.open('rb'), 'sha256').hexdigest())
    result['seed_sha256'] = seed_sha256
    assert seed_sha256 == {name: hashlib.file_digest((seed/name).open('rb'), 'sha256').hexdigest()
                           for name in seed_sha256}
    (directory/'run.receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
