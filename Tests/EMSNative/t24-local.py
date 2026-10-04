#!/usr/bin/env python3
"""Bounded local d507438 restart investigation; no remote or production edits."""
import argparse
import csv
from collections import defaultdict
import difflib
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = Path('/private/tmp/ems-t24')
INPUT = Path('/Users/auroradysis/Workspace/EMS/.data/exp-0024/merger-stall')
RUN = INPUT/'ev1/runs/exp-0024/merger'
SUBMISSION = Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0024/submissions/exp-0024')
MPI = Path('/Users/auroradysis/.julia/artifacts/2e0bbdf5bae18755b0ebba03bc02d6fb568c05bc')
GCC = '/opt/homebrew/bin/g++-16'
PY = sys.executable

def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(path)

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def mpi_environment():
    hw = []
    for p in Path('/Users/auroradysis/.julia/artifacts').glob('*/lib/libhwloc.15.dylib'):
        if 'arm64' in subprocess.check_output(['file', str(p)], text=True):
            hw.append(str(p.parent))
    return dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1',
                MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
                TMPDIR='/private/tmp', MPIR_CVAR_NOLOCAL='1',
                DYLD_LIBRARY_PATH=':'.join([str(MPI/'lib'), *hw]))

def mpi_probe():
    ROOT.mkdir(exist_ok=True)
    env = mpi_environment()
    cmd = [GCC, '-O2', '-fno-fast-math', '-ffp-contract=off', '-fopenmp',
           '-I'+str(MPI/'include'), str(HERE/'t24-mpi-probe.cpp'),
           '-L'+str(MPI/'lib'), '-Wl,-rpath,'+str(MPI/'lib'),
           '-lmpi', '-o', str(ROOT/'mpi-probe.ex')]
    subprocess.run(cmd, env=env, check=True)
    launch = [str(MPI/'bin/mpiexec.hydra'), '-launcher', 'fork', '-n', '2',
              str(ROOT/'mpi-probe.ex')]
    try:
        p = subprocess.run(launch, env=env, capture_output=True, text=True, timeout=15)
        q = dict(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr,
                 command=launch, mpi=str(MPI),
                 env={k:env[k] for k in ('DYLD_LIBRARY_PATH', 'OMP_NUM_THREADS', 'MPIR_CVAR_NOLOCAL',
                                      'MPIR_CVAR_NUM_CLIQUES','MPIR_CVAR_CH4_SHM') if k in env})
    except subprocess.TimeoutExpired as e:
        q = dict(returncode=124, stdout=str(e.stdout), stderr=str(e.stderr), command=launch)
    old = ROOT/'mpi-probe.json'
    if old.exists():
        old.rename(ROOT/f'mpi-probe-previous-{time.time_ns()}.json')
    write(old, q)
    print(json.dumps(q, indent=2), flush=True)
    return q['returncode']

def debugger_probe():
    d = ROOT/'debugger-probe'
    rc = measured('lldb', ['/usr/bin/lldb', '--batch', '-o',
                  'settings set target.disable-aslr false', '-o', 'run',
                  '-o', 'thread backtrace all', '--',
                  '/private/tmp/ems-t23/builds/tracking/production.ex',
                  '/private/tmp/ems-t23/runs/tracking-production/params.txt',
                  'just_check_params=1'], d)
    text = (d/'run.log').read_text()
    write(d/'debugger.receipt.json', dict(returncode=rc, available=rc==0,
           diagnostic=text[-5000:], scope='sandbox launch availability, no advances'))
    # An unavailable debugger is a recorded limitation, not a failed evolution.
    print('T24_DEBUGGER_AVAILABLE', rc==0, flush=True)

def sparse_probe():
    import h5py
    d = ROOT/'sparse-probe'
    d.mkdir(exist_ok=True)
    audit = load_module('t24_boxreader', SUBMISSION/'smoke-audit.py')
    lines = []
    with h5py.File(INPUT/'chk/EMS_000152.2d.hdf5') as f:
        for level in range(1,13):
            for b in f[f'level_{level}/boxes'][:]:
                lines.append(' '.join(map(str,[level,*[int(b[k]) for k in audit.BOX_KEYS]]))+'\n')
    (d/'boxes.txt').write_text(''.join(lines))
    ch = Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
    spec = json.loads(Path('/private/tmp/ems-t23/builds/tracking/build-spec.json').read_text())
    flags = spec['flags'].copy()
    flags = [f for f in flags if f!='-O3']
    flags += ['-O1','-g','-fsanitize=bounds','-fno-sanitize-recover=bounds',
              '-fno-fast-math','-ffp-contract=off','-fno-omit-frame-pointer']
    cmds = [flags+['-c',str(ch/'src/BoxTools/TreeIntVectSet.cpp'),'-o',str(d/'TreeIntVectSet.o')],
            flags+[str(HERE/'t24-sparse-probe.cpp'),str(d/'TreeIntVectSet.o'),
                   *spec['libs'],'-o',str(d/'sparse-probe.ex')]]
    for i, command in enumerate(cmds):
        rc = measured('compile-'+str(i),command,d/('compile-'+str(i)))
        assert rc==0, (i,rc)
    rc = measured('sparse-native',[d/'sparse-probe.ex',d/'boxes.txt'],d/'native',
                  dict(UBSAN_OPTIONS='print_stacktrace=1'))
    text = (d/'native/run.log').read_text()
    write(HERE/'t24-sparse-probe.json',dict(returncode=rc,
        status='ARRAY_BOUNDS_DEFECT_OBSERVED' if 'out of bounds' in text else
               ('PROBE_COMPLETED' if rc==0 else 'OTHER_PROBE_FAILURE'),
        diagnostic=text[-7000:], commands=cmds,
        interpretation='a bounds defect is not yet causal proof of the step155 stall',
        native_receipt=str(d/'native/sparse-native.resources.json')))
    print('T24_SPARSE_PROBE_ENDED',rc,flush=True)

def replace_params(text, changes):
    for k, v in changes.items():
        pat = r'^'+re.escape(k)+r'\s*=.*$'
        line = f'{k} = {v}'
        text = (re.sub(pat, line, text, flags=re.M) if re.search(pat, text, re.M)
                else text+'\n'+line+'\n')
    return text

def input_audit():
    import h5py
    import numpy as np
    manifest = INPUT/'ev1/stall-evidence-1.sha256'
    verified = []
    for line in manifest.read_text().splitlines():
        digest, name = line.split(None, 1)
        p = INPUT/'ev1'/name.lstrip('*')
        assert p.is_file() and sha(p) == digest, p
        verified.append(dict(path=str(p), sha256=digest, bytes=p.stat().st_size))
    checkpoint = INPUT/'chk/EMS_000152.2d.hdf5'
    seal = json.loads(checkpoint.with_suffix('.sealed.json').read_text())
    assert seal['complete'] and checkpoint.stat().st_size == seal['bytes']
    audit = load_module('t24_published_audit', SUBMISSION/'smoke-audit.py')
    actual = audit.state(checkpoint, 152)
    assert actual == seal['endpoint'], (actual, seal['endpoint'])
    rows = []
    layouts = []
    with h5py.File(checkpoint) as f:
        for level in range(int(f.attrs['num_levels'])):
            g = f[f'level_{level}']
            boxes = [tuple(int(b[k]) for k in audit.BOX_KEYS) for b in g['boxes'][:]]
            domain = tuple(int(g.attrs['prob_domain'][k]) for k in audit.BOX_KEYS)
            layouts.append((boxes, domain))
            for i, b in enumerate(boxes):
                rows.append(dict(level=level, time_M=float(g.attrs['time']),
                                 h_M=float(g.attrs['dx']), box=i,
                                 **dict(zip(audit.BOX_KEYS, b))))
    save_csv(HERE/'t24-checkpoint-boxes.csv', rows)
    nesting = []
    for level in range(1, len(layouts)):
        parent, domain = layouts[level-1]
        children, _ = layouts[level]
        intervals = defaultdict(list)
        for x0, y0, x1, y1 in parent:
            for y in range(y0, y1+1):
                intervals[y].append((x0, x1))
        for y in intervals:
            intervals[y].sort()
        for pad in (0, 4, 8):
            missing = 0
            example = None
            for x0, y0, x1, y1 in children:
                x0, y0, x1, y1 = x0//2-pad, y0//2-pad, x1//2+pad, y1//2+pad
                x0, y0 = max(x0, domain[0]), max(y0, domain[1])
                x1, y1 = min(x1, domain[2]), min(y1, domain[3])
                for y in range(y0, y1+1):
                    cursor = x0
                    for a, b in intervals.get(y, []):
                        if b < cursor: continue
                        if a > x1: break
                        if a > cursor:
                            missing += min(a-1, x1)-cursor+1
                            if example is None: example = (cursor, y)
                        cursor = max(cursor, b+1)
                        if cursor > x1: break
                    if cursor <= x1:
                        missing += x1-cursor+1
                        if example is None: example = (cursor, y)
            nesting.append(dict(child_level=level, parent_buffer_cells=pad,
                                missing_parent_cells_with_multiplicity=missing,
                                first_missing=str(example),
                                scope='saved step152; clipped physical boundaries'))
    save_csv(HERE/'t24-checkpoint-nesting.csv', nesting)
    assert all(r['missing_parent_cells_with_multiplicity']==0 for r in nesting if r['parent_buffer_cells']==0)
    lb = RUN/'segments/engine-step000156/interrupted-1791066562058016119/LB.txt'
    blocks = re.split(r'(?=\s*LoadBalance\s+sizes:)', lb.read_text())
    counts = []
    for i, block in enumerate(blocks):
        n = re.search(r'a_boxes=(\d+)', block)
        if not n: continue
        rankrows = re.findall(r'rank(\d+)\s+load=(\d+)\s+boxes\((\d+)\)', block)
        counts.append(dict(block=i, boxes=int(n[1]), ranks_recorded=len(rankrows),
                           min_boxes_per_recorded_rank=min(int(r[2]) for r in rankrows) if rankrows else '',
                           max_boxes_per_recorded_rank=max(int(r[2]) for r in rankrows) if rankrows else '',
                           complete_224_rank_counts=len(rankrows)==224))
    save_csv(HERE/'t24-load-balance-counts.csv', counts)
    punctures = np.loadtxt(RUN/'punctures.dat')
    last = [dict(time_M=float(row[0]), x1=float(row[1]), y1=float(row[2]),
                 x2=float(row[3]), y2=float(row[4]), separation_M=float(row[3]-row[1]))
            for row in punctures if row[0]>=133.0]
    save_csv(HERE/'t24-last-punctures.csv', last)
    result = dict(status='INPUTS_VERIFIED', checkpoint=str(checkpoint),
                  checkpoint_bytes=checkpoint.stat().st_size,
                  checkpoint_whole_file_sha256=sha(checkpoint),
                  seal_endpoint=actual, verified_small_files=len(verified),
                  small_manifest_sha256=sha(manifest),
                  published_audit_sha256=sha(SUBMISSION/'smoke-audit.py'),
                  nesting=nesting, load_balance_last_blocks=counts[-13:],
                  post_stall_box_coordinates_available=False,
                  warning='LB reports counts/loads, not coordinates; checkpoint-layout.txt is Lustre striping')
    write(HERE/'t24-input-audit.json', result)
    print('T24_INPUT_AUDIT_COMPLETE', len(verified), actual['valid_field_sha256'], flush=True)

def save_csv(path, rows):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def prepare():
    ROOT.mkdir(exist_ok=True)
    src = ROOT/'source'
    commit = subprocess.check_output(['git', 'rev-parse', 'd507438'], cwd=REPO, text=True).strip()
    if not src.exists():
        src.mkdir()
        archive = ROOT/'source.tar'
        with archive.open('wb') as f:
            subprocess.run(['git', 'archive', commit], cwd=REPO, stdout=f, check=True)
        subprocess.run(['tar', '-xf', str(archive), '-C', str(src)], check=True)
    original = src/'Source/GRChomboCore/GRAMRLevel.cpp'
    text = original.read_text()
    needle = '\n}\n\n// initialize grid\n'
    assert text.count(needle)==1
    instrumented = '#include "t24-observer.hpp"\n'+text.replace(
        needle, '\n    t24_record_layout(m_grids, m_level, m_time, m_dx);'+needle)
    observer = ROOT/'GRAMRLevel-observer.cpp'
    observer.write_text(instrumented)
    (HERE/'t24-observer.patch').write_text(''.join(difflib.unified_diff(
        text.splitlines(True), instrumented.splitlines(True),
        fromfile='a/Source/GRChomboCore/GRAMRLevel.cpp',
        tofile='b/Source/GRChomboCore/GRAMRLevel.cpp')))
    original_ch = Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
    ch = ROOT/'Chombo/lib'
    if not ch.exists():
        ch.mkdir(parents=True)
        ignore = shutil.ignore_patterns('.git', 'o', 'd', '*.o', '*.d', '*.a', '*.ex', '*.dSYM')
        for name in ('src', 'mk', 'util'):
            shutil.copytree(original_ch/name, ch/name, ignore=ignore)
        shutil.copy2(original_ch/'GNUmakefile', ch/'GNUmakefile')
        local = (original_ch/'mk/Make.defs.local').read_text()
        strict = '-O3 -g -fno-fast-math -ffp-contract=off -fno-omit-frame-pointer'
        local = re.sub(r'^cxxoptflags\s*=.*$', 'cxxoptflags = '+strict, local, flags=re.M)
        local = re.sub(r'^foptflags\s*=.*$', 'foptflags = '+strict, local, flags=re.M)
        (ch/'mk/Make.defs.local').write_text(local)
    common = dict(static_read_policy='t=0 only; restart paths deliberately absent',
                  fork_commit=commit, compiler=GCC,
                  production_source_modified=False, instrumentation='private observer build only',
                  Float64=True, fast_math=False, FMA_contraction='off',
                  rank_count=1, OpenMP_threads=2,
                  layout_scope='original physical grid; local load balancing differs from 224 ranks',
                  MPI_limit='installed MPICH fails POSIX SHM bootstrap inside sandbox',
                  cap_bytes=6000000000, preventive_gate_bytes=5700000000)
    write(HERE/'t24-build-registration.json', common)
    jobs = []
    def job(name, directory, command):
        d = ROOT/directory
        d.mkdir(parents=True, exist_ok=True)
        jobs.append(dict(name=name, directory=str(d), command=list(map(str, command))))
    job('inputs', 'inputs', [PY, HERE/'t24-local.py', 'input_audit'])
    job('debugger-probe', 'debugger-probe-job', [PY, HERE/'t24-local.py', 'debugger_probe'])
    job('strict-build', 'build', [PY, HERE/'t24-local.py', 'build'])
    job('observer-regression', 'regression', [PY, HERE/'t24-local.py', 'regression'])
    job('restart152', 'restart152', [PY, HERE/'t24-local.py', 'restart'])
    pipeline = ROOT/'pipeline'
    pipeline.mkdir(exist_ok=True)
    write(pipeline/'plan.json', dict(jobs=jobs))
    print('T24_PLAN_READY', pipeline/'plan.json', flush=True)

def build():
    ch = ROOT/'Chombo/lib'
    d = ROOT/'build'
    d.mkdir(exist_ok=True)
    env = dict(os.environ, TMPDIR='/private/tmp', OMP_NUM_THREADS='2')
    commands = []
    old_commands = []
    if (d/'commands.json').exists():
        old_commands = json.loads((d/'commands.json').read_text())
        (d/'commands.json').rename(d/f'commands-before-{time.time_ns()}.json')
    def call(cmd, cwd=d):
        commands.append(dict(command=list(map(str,cmd)), directory=str(cwd)))
        write(d/'commands.json', commands)
        subprocess.run(list(map(str, cmd)), cwd=cwd, env=env, check=True)
    call(['make', '-j2', 'BaseTools', 'BoxTools', 'AMRTools', 'AMRTimeDependent'], ch)
    src = ROOT/'source'
    hdf = Path('/Users/auroradysis/Workspace/EMS-deps/hdf5-1.14')
    flags = [GCC, '-O3', '-g', '-fno-fast-math', '-ffp-contract=off',
             '-fno-omit-frame-pointer', '-fno-lto', '-std=c++17', '-fopenmp']
    flags += ['-D'+s for s in ('CH_SPACEDIM=2','CH_Darwin','NDEBUG','CH_USE_COMPLEX',
              'CH_USE_64','CH_USE_DOUBLE','CH_USE_HDF5','H5_USE_16_API','CH_USE_LAPACK',
              'CH_FORT_UNDERSCORE','CH_LANG_CC')]
    flags += ['-I'+str(src/'Examples/EMS'), '-I'+str(HERE), '-I'+str(hdf/'include')]
    flags += ['-I'+str(p) for p in (src/'Source').rglob('*') if p.is_dir()]
    flags += ['-I'+str(ch/'src'/s) for s in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
    libs = ['-L'+str(ch)]+['-l'+s+'2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
                                        for s in ('amrtimedependent','amrtools','boxtools','basetools')]
    libs += ['-L'+str(hdf/'lib'), '-Wl,-rpath,'+str(hdf/'lib'), '-lhdf5', '-lz',
             '-L/opt/homebrew/lib/gcc/16', '-lgfortran', '-lm', '-lgomp', '-framework', 'Accelerate']
    objects = []
    for stem in ('GRAMRLevel','GRAMR','GRLevelData','BoundaryConditions','SmallDataIO',
                 'PunctureTracker','PETScCommunicator','EMSBH2DLevel'):
        matches = list(src.rglob(stem+'.cpp'))
        assert len(matches)==1, matches
        obj = d/(stem+'.o')
        cmd = list(map(str,flags+['-c',matches[0],'-o',obj]))
        # Only unchanged, fully compiled units from this strict-source build.
        if obj.exists() and any(r['command']==cmd for r in old_commands):
            commands.append(dict(command=cmd, directory=str(d),
                                 reuse='unchanged d507438 source and identical strict flags',
                                 object_sha256=sha(obj)))
            write(d/'commands.json', commands)
        else:
            call(cmd)
        objects.append(obj)
    call(flags+['-c',ROOT/'GRAMRLevel-observer.cpp','-o',d/'GRAMRLevel-observer.o'])
    call(flags+['-c',HERE/'t24-stack.cpp','-o',d/'stack.o'])
    main = src/'Examples/EMS/Main_EMSBH2DBH.cpp'
    call(flags+[main,*objects,*libs,'-o',d/'production.ex'])
    observer_objects = [d/'GRAMRLevel-observer.o' if p.name=='GRAMRLevel.o' else p for p in objects]
    call(flags+[main,*observer_objects,d/'stack.o',*libs,'-o',d/'observer.ex'])
    call(flags+[HERE/'t24-stack-selftest.cpp',d/'stack.o','-o',d/'stack-selftest.ex'])
    write(d/'build-spec.json', dict(flags=flags, libraries=libs,
          compiler=subprocess.check_output([GCC,'--version'],text=True).splitlines()[0],
          executable_sha256={n:sha(d/(n+'.ex')) for n in ('production','observer','stack-selftest')},
          chombo_library_sha256={p.name:sha(p) for p in ch.glob('*.a')},
          original_GRAMRLevel_sha256=sha(src/'Source/GRChomboCore/GRAMRLevel.cpp'),
          observer_cpp_sha256=sha(ROOT/'GRAMRLevel-observer.cpp')))
    print('T24_STRICT_BUILD_COMPLETE', flush=True)

def measured(name, command, directory, extra_env=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TMPDIR='/private/tmp', T13_OUTPUT_ROOT=str(ROOT),
               T13_PROCESS_CAP_BYTES='5700000000', T13_TREE_CAP_BYTES='7500000000',
               T13_DISK_CAP_BYTES='16000000000')
    env.update(extra_env or {})
    with (directory/'run.log').open('wb') as log:
        p = subprocess.run([PY, str(HERE/'t13-run.py'), '--measure', name,
                            '--directory', str(directory), '--', *map(str, command)],
                           stdout=log, stderr=subprocess.STDOUT, env=env)
    (directory/'done.exit').write_text(str(p.returncode)+'\n')
    return p.returncode

def regression():
    compare = load_module('t24_compare', HERE/'t23-compare.py')
    base = (Path('/private/tmp/ems-t23/runs/tracking-production/params.txt')).read_text()
    # Same source parameters as the prior ordinary production-main fixture;
    # enable a regrid to exercise the new read-only layout observer.
    base = replace_params(base, dict(regrid_interval='1', max_steps=1))
    summary = []
    cases = [('plain','production',False), ('observer-off','observer',False),
             ('observer-on','observer',True), ('plain-regrid','production',False),
             ('observer-off-regrid','observer',False), ('observer-on-regrid','observer',True)]
    for name, executable, on in cases:
        d = ROOT/'regression'/name
        for sub in ('chk','plt'): (d/sub).mkdir(parents=True, exist_ok=True)
        params = replace_params(base,dict(max_steps=2)) if name.endswith('-regrid') else base
        env = dict(T24_STACK_PREFIX=str(d/'stacks'), T24_LAYOUT_FILE=str(d/'layouts.csv')) if on else {}
        if (d/'done.exit').exists() and (d/'done.exit').read_text().strip()=='0':
            assert (d/'params.txt').read_text()==params
            q = json.loads((d/(name+'.resources.json')).read_text())
            assert q['returncode']==0 and not q['gate_reason'] and q['peak_rss_bytes']<6e9
            rc = 0
        else:
            (d/'params.txt').write_text(params)
            rc = measured(name, [ROOT/'build'/(executable+'.ex'), d/'params.txt'], d, env)
        assert rc==0, (name, rc)
        text = (d/'run.log').read_text()
        assert 'GRChombo finished.' in text
        if on:
            assert len(list(d.glob('stacks.thread-*.txt')))==2
            if name.endswith('-regrid'):
                assert (d/'layouts.csv').is_file() and (d/'layouts.csv').stat().st_size>0
        else:
            assert not list(d.glob('stacks.thread-*.txt')) and not (d/'layouts.csv').exists()
    for candidate in ('observer-off', 'observer-on','observer-off-regrid','observer-on-regrid'):
        reference = 'plain-regrid' if candidate.endswith('-regrid') else 'plain'
        for step in ((0,1,2) if candidate.endswith('-regrid') else (0,1)):
            for prefix in ('chk/EMS_', 'plt/EMS_Plot_'):
                file = f'{prefix}{step:06d}.2d.hdf5'
                rows = compare.comparisons(ROOT/'regression'/reference/file,
                                           ROOT/'regression'/candidate/file,
                                           'plain vs '+candidate)
                evolved = set('chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split())
                defined = [r for r in rows if r['field'] in evolved or r['region'] in ('valid','same_level_ghost')]
                assert all(r['bit_mismatches']==0 for r in defined), defined
                summary.append(dict(candidate=candidate, file=file,
                    defined_values=sum(r['values'] for r in defined), bit_mismatches=0,
                    undefined_plot_diagnostic_ghosts_excluded=True))
    d = ROOT/'regression/stack-selftest'
    d.mkdir(exist_ok=True)
    rc = measured('stack-selftest', [ROOT/'build/stack-selftest.ex'], d,
                  dict(T24_STACK_PREFIX=str(d/'stacks')))
    assert rc==0
    assert len(list(d.glob('stacks.thread-*.txt')))==2
    for p in d.glob('stacks.thread-*.txt'):
        assert 'T24_SIGNAL_STACK_BEGIN' in p.read_text() and 'T24_SIGNAL_STACK_END' in p.read_text()
    save_csv(HERE/'t24-observer-regression.csv', summary)
    write(HERE/'t24-observer-qualification.json', dict(status='PASS', comparisons=summary,
          thread_stack_selftest='both two OpenMP threads dumped',
          production_source_changed=False, native_RHS_or_evolution_changed=False))
    print('T24_OBSERVER_REGRESSION_COMPLETE', flush=True)

def restart():
    assert json.loads((HERE/'t24-observer-qualification.json').read_text())['status']=='PASS'
    assert json.loads((HERE/'t24-input-audit.json').read_text())['status']=='INPUTS_VERIFIED'
    d = Path(os.environ.get('T24_RESTART_DIRECTORY',str(ROOT/'restart152')))
    rss_only = os.environ.get('T24_RSS_GATE','false')=='true'
    for sub in ('chk','plt'): (d/sub).mkdir(parents=True, exist_ok=True)
    # The cluster driver copies its numerical histories into a writable run
    # directory. The tracker reads the line at t=133 and removes later lines.
    shutil.copy2(RUN/'punctures.dat', d/'punctures.dat')
    absent = d/'ABSENT-TRUMPET-AND-COMPANION'
    assert not absent.exists()
    changes = dict(restart_file=INPUT/'chk/EMS_000152.2d.hdf5',
                   ems_data_path=absent, ems_ctt_data_path=absent,
                   max_steps=156, stop_time='1e100',
                   chk_prefix=d/'chk/EMS_', plot_prefix=d/'plt/EMS_Plot_',
                   data_subpath=d)
    text = replace_params((RUN/'runtime/params-production.txt').read_text(), changes)
    (d/'params.txt').write_text(text)
    write(d/'overrides.json', {k:str(v) for k,v in changes.items()})
    (d/'layouts.csv').write_text('phase,level,time_M,h_M,box,rank,lo_i,lo_j,hi_i,hi_j\n')
    env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', TMPDIR='/private/tmp',
               T24_STACK_PREFIX=str(d/'stacks'), T24_LAYOUT_FILE=str(d/'layouts.csv'))
    command = [str(ROOT/'build/observer.ex'), str(d/'params.txt')]
    start = time.monotonic()
    snapshots = []
    reason = ''
    last_progress = start
    signatures = None
    next_dump = 300
    monitor = load_module('t24_monitor', HERE/'t13-run.py')
    peaks = {}
    tree_peak = 0
    with (d/'native.log').open('wb') as log:
        p = subprocess.Popen(command, cwd=d, env=env, stdout=log,
                             stderr=subprocess.STDOUT, start_new_session=True)
        (d/'native.pid').write_text(str(p.pid)+'\n')
        while p.poll() is None:
            now = time.monotonic()
            usage = monitor.tree_usage(p.pid)
            current_tree = sum(max(v) for v in usage.values())
            tree_peak = max(tree_peak, current_tree)
            for pid, values in usage.items():
                old = peaks.get(pid, (0,0))
                peaks[pid] = tuple(max(a,b) for a,b in zip(old,values))
            if any((v[0] if rss_only else max(v))>5.7e9 for v in usage.values()):
                reason = 'PREVENTIVE_MEMORY_GATE_5.7GB_BEFORE_6GB_CAP'
            if current_tree>7.5e9:
                reason = 'AGGREGATE_FOOTPRINT_GATE_7.5GB'
            files = [d/'native.log', *d.glob('pout.*')]
            current = tuple((str(f),f.stat().st_size) for f in files if f.is_file())
            if current != signatures:
                signatures = current
                last_progress = now
                next_dump = 300
            quiet = now-last_progress
            if quiet >= next_dump:
                os.kill(p.pid, signal.SIGUSR1)
                snapshots.append(dict(elapsed_s=now-start, quiet_s=quiet, signal='SIGUSR1'))
                write(d/'watchdog.json', dict(snapshots=snapshots, status='STACK_REQUESTED'))
                next_dump += 300
                if quiet >= 900:
                    reason = 'QUIET_WATCHDOG_AFTER_THREE_STACK_REQUESTS_NOT_A_NUMERICAL_VERDICT'
            if now-start>4*3600:
                reason = 'LOCAL_WALL_CAP_4H_NOT_A_NUMERICAL_VERDICT'
            if reason:
                os.killpg(p.pid, signal.SIGTERM)
                try: p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                break
            time.sleep(.25)
        rc = p.wait()
    elapsed = time.monotonic()-start
    receipt = dict(native_returncode=rc, gate_reason=reason, elapsed_s=elapsed,
                   command=command, per_pid_peak_rss_and_footprint_bytes=peaks,
                   tree_peak_bytes=tree_peak, watchdog_requests=snapshots,
                   complete_to_step156=False, no_static_files_exist=not absent.exists(),
                   rank_count=1, OpenMP_threads=2,
                   per_process_cap_kind='RSS' if rss_only else 'RSS and footprint',
                   per_process_RSS_cap_bytes=6000000000,
                   preventive_RSS_gate_bytes=5700000000,
                   native_wait4_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                   watchdog_process_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if rc==0:
        assert 'GRChombo finished.' in (d/'native.log').read_text()
        audit = load_module('t24_endpoint', SUBMISSION/'smoke-audit.py')
        receipt['endpoint'] = audit.state(d/'chk/EMS_000156.2d.hdf5',156)
        receipt['complete_to_step156'] = True
    write(d/'native.receipt.json', receipt)
    print('T24_NATIVE_RESTART_ENDED', json.dumps(receipt), flush=True)
    return 125 if reason else rc

def replay_worker():
    """Direct owner: native RSS gate remains enforced; no outer footprint gate."""
    d = ROOT/'restart152-rss'
    os.environ['T24_RESTART_DIRECTORY'] = str(d)
    os.environ['T24_RSS_GATE'] = 'true'
    try:
        rc = restart()
    except BaseException as exc:
        write(d/'worker-error.json',dict(type=type(exc).__name__,message=str(exc)))
        rc = 1
        raise
    finally:
        tmp = d/'done.exit.tmp'
        tmp.write_text(str(rc)+'\n')
        tmp.replace(d/'done.exit')
    return rc

def detach_replay():
    d = ROOT/'restart152-rss'
    d.mkdir(exist_ok=False)
    with (d/'run.log').open('wb') as log:
        p = subprocess.Popen([PY,str(HERE/'t24-local.py'),'replay_worker'],
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
            env=dict(os.environ,TMPDIR='/private/tmp'))
    (d/'launcher.pid').write_text(str(p.pid)+'\n')
    write(d/'launch.json',dict(pid=p.pid, mode='RSS-capped replay',
        native_RSS_cap_bytes=6000000000, preventive_RSS_gate_bytes=5700000000,
        aggregate_footprint_gate_bytes=7500000000, OpenMP_threads=2,
        source='same qualified observer executable and same physical parameters',
        done_marker=str(d/'done.exit')))
    print('T24_DETACHED_RSS_REPLAY_PID',p.pid,'MARKER',d/'done.exit',flush=True)

def packet():
    """Refresh the provisional packet without treating launcher exits as passes."""
    if (HERE/'t24-offline-receipt.json').exists():
        # The completed read-only tranche supersedes the historical replay packet.
        return subprocess.call([PY, str(HERE/'t24-seal.py')])
    audit = json.loads((HERE/'t24-input-audit.json').read_text())
    qualification_path = HERE/'t24-observer-qualification.json'
    qualified = json.loads(qualification_path.read_text()) if qualification_path.exists() else None
    active = ROOT/'restart152-rss' if (ROOT/'restart152-rss/launch.json').exists() else ROOT/'continuation-3'
    native_directory = ROOT/'restart152-rss' if active.name=='restart152-rss' else ROOT/'restart152'
    native_path = native_directory/'native.receipt.json'
    native = json.loads(native_path.read_text()) if native_path.exists() else None
    resources = []
    for p in ROOT.rglob('*.resources.json'):
        q = json.loads(p.read_text())
        resources.append(dict(receipt=str(p), process=q['process'],
                              returncode=q['returncode'], peak_RSS_bytes=q['peak_rss_bytes'],
                              wall_seconds=q['wall_seconds'], gate_reason=q['gate_reason']))
    save_csv(HERE/'t24-resources.csv', resources)
    peak = max(r['peak_RSS_bytes'] for r in resources)
    sha_cp = audit['checkpoint_whole_file_sha256']
    status = dict(status='READY-EXCEPT', cause='UNESTABLISHED', minimal_fix='NONE_APPLIED',
                  checkpoint_verified=True, observer_qualification=qualified,
                  restart_receipt=native, maximum_completed_job_peak_RSS_bytes=peak,
                  original_pipeline_marker=str(ROOT/'pipeline/done.exit'),
                  continuation_marker=str(active/'done.exit'),
                  native_job_marker=str(native_directory/'done.exit'))
    write(HERE/'t24-status.json', status)
    report = f'''# T24: deterministic binary merger stall — source audit and local replay

**READY-EXCEPT. The causal stall site and minimal fix are not established.** The
step-152 numerical input is verified and the bounded replay is registered. This
packet records completed evidence separately from the pending restart. No
production source, equations, gauge, KO, transfer rule, or refinement footprint
has been changed. In particular, no speculative fix is proposed as qualified.

The pinned source is `d50743867e48bcb466a60ca16869b07fe1546c27`, archived under
`/private/tmp/ems-t24/source/`; the dirty working tree is not used to build it.
The source line numbers below refer to that archive. The current T22/T23 edits
are preserved. The archived static setter is initial-data-only, and the native
restart has both static input paths set to a deliberately absent local path.
The guard remains enabled.

The supplied checkpoint is **1,500,693,480 bytes**, step **152**, time **133 M_i**.
Its 28 valid evolved fields, including covered coarse cells, reproduce the
published partition-independent seal digest
`7146c763475ad64b73a62f19ab2e17f5e091a3da6891698ccda8a18c03b389cc` exactly.
The independently recorded whole-file SHA-256 is `{sha_cp}`. All **1,855**
files in `ev1/stall-evidence-1.sha256` match. The published audit's own checksum
is recorded in [t24-input-audit.json](t24-input-audit.json). This audit's own
return code is 0 and peak RSS is **896,090,112 bytes**. The sandboxed time
wrapper's `sysctl` error is separate from the native/analysis exit.

## What the source establishes

`Examples/EMS/EMSBH2DLevel.cpp:430–480` first computes the ordinary criterion,
removing the midpoint's forced extraction radii for a binary, then takes the
cellwise maximum of two centred criteria. The finite centre loop starts at
line 455; there is no division by separation or search-until-disjoint loop.
The native radius forcing is `r < 1.2*r_extraction` in
`Source/TaggingCriteria/EMSExtractionTaggingCriterion.hpp`. Coincident finite
centres make the forced-tag union idempotent; this statement does not prove
the sparse-set or distributed regrid implementations safe.

The last recorded separations are **0.4519610**, **0.4149546**, and **0.3813584 M_i**
at **133**, **133.875**, and **134.75 M_i**, respectively, as independently read
in [t24-last-punctures.csv](t24-last-punctures.csv). The prescribed forcing
radius for creating level 10 is **0.21875 M_i** around each centre. Those
nominal circles overlap once separation falls below **0.4375 M_i**. Overlap is
a geometric inference from the parameters and recorded centres, not a measured
post-stall box-layout audit.

`Examples/EMS/EMSBH2DLevel.cpp:486–488` executes tracking on level 6, after its
children complete; it writes punctures at coarse clocks. However,
`Examples/EMS/Main_EMSBH2DBH.cpp:44–45` fixes the interpolation minimum level to
`max_level-1`, hence **11**, rather than to the tracking execution level.
`Source/BlackHoles/PunctureTracker.cpp:145–179` advances a fixed-size list with
a trapezoidal shift estimate; `:188–249` refreshes the interpolator, fills
shift ghosts from level 11, and submits both coordinates collectively. No
coincidence transition or automatic label merge exists. Duplicate queries are
legal in the query interface. The implementation's Lagrange stencil growth in
`Source/AMRInterpolator/Lagrange.impl.hpp:277` either accepts another stencil
point or disables a growth direction; there is no separation-dependent
unbounded search in this function. MPI collectives in
`Source/AMRInterpolator/MPIContext.impl.hpp:49,77,96` remain possible wait sites
if ranks enter different collective sequences; no such divergence is proved.

`Source/GRChomboCore/GRAMRLevel.cpp:387–470` installs the new box list, load
balances it, reconstructs exchange and coarse-fine support, interpolates new
cells, restores old overlaps, fills boundaries and redefines auxiliary storage.
The nonzero-rank `regrid` progress message is **at the end** (line 470).
Advance progress at `:179–181` is printed before the advance's copy/RK work,
and normally only on rank 0. Therefore the other ranks' last `regrid` messages
do not locate their current stalls inside regrid. They can already be in an
unprinted advance or collective. Rank 0's last advance message likewise does
not identify an RK stage. `evalRHS` exchanges evolved stencils and performs
point time interpolation before evaluating the native RHS and KO.

The point-transfer ghost-set and interpolation loops in
`Source/GRChomboCore/PointAMRTransfer.hpp:57–118` are finite in their valid
inputs. A low-level data-structure defect is still possible: the external
Chombo `TreeIntVectSet.H:314–319` uses fixed **24-entry**, OpenMP thread-private
traversal scratch arrays. A narrow bounds-instrumented probe of grow,
subtract-own-box and iteration completed for **2,704** saved level-1–12 boxes
without a bounds report. See [t24-sparse-probe.json](t24-sparse-probe.json) and
[t24-sparse-probe.cpp](t24-sparse-probe.cpp). This does not cover neighbour
subtraction, tag unions/growth, BR mesh refinement, the new layout or corrupted
inputs, and does not eliminate the mechanism.

Another source-level nontermination condition exists in external Chombo
`AMR.cpp:1210–1219`: its `while(stepsLeft>0)` does not consume `stepsLeft` if
the current level ceases to have a finer level during a subcycle. That condition
is not established here: the last reported hierarchy still has levels 11 and
12 with 256 boxes. It is a candidate to check in stacks, not this report's cause.

## Hierarchy and load balance

[t24-checkpoint-boxes.csv](t24-checkpoint-boxes.csv) records the actual step-152
layout. [t24-checkpoint-nesting.csv](t24-checkpoint-nesting.csv) checks every
coarsened child rectangle against the parent union, with **0**, **4**, and **8**
parent-cell buffers, clipped only at the physical domain. All **36** tests have
zero missing parent cells. This proves the stated coverage at **step 152**;
it does not certify the layout created at **135.174 M_i**.

The supplied `checkpoint-layout.txt` is a Lustre striping report. The native
`LB.txt` has loads and box counts, not box coordinates or IDs, because external
Chombo `LoadBalance.cpp:443` sets `printEveryBox=false`. Complete 240-box records
include **all 224 ranks**, each with **1–2 boxes**, so the reduction from 256
to 240 is not evidence of an empty rank at that level. The final 256-box record
is truncated after 203 rank records. These distinctions are tabulated in
[t24-load-balance-counts.csv](t24-load-balance-counts.csv). The decrease can be
consistent with changing geometry/partitioning as footprints overlap, but its
actual geometric cause and post-regrid proper nesting require the missing
coordinates. The private observer writes those coordinates after every local
regrid.

## Bounded local replay and current receipts

[t24-local.py](t24-local.py) registers, builds and runs the evidence queue.
The initial input audit and debugger availability probe completed. The first
strict build stopped on a missing `LayoutIterator.H` include in the Tests-only
observer. That include is repaired. The next link exposed that the Chombo top-level
AMRTimeDependent target had not built its dependency libraries; the build now names
all four libraries explicitly. Completed unchanged native objects are reused.
Both failed receipts and markers remain under `/private/tmp/ems-t24/pipeline/`,
`continuation/` and `build/`; all repaired build and regression jobs completed.
The first regression qualification also stopped because its one-step fixture
does not reach the first scheduled regrid. Those completed native cases are
retained and reused; additional two-step fixtures now exercise and compare the
actual layout callback. No bitwise check is weakened. All one-step and two-step
observer cases and the two-thread stack self-test now pass.

LLDB fails before launching the target (`process exited with status -1 (no such
process)`); its own receipt is `debugger-probe/debugger.receipt.json`. The
installed local MPICH also fails during POSIX shared-memory bootstrap before
the probe collective, despite tested network-only/clique controls. These are
environment limitations, not merger failures. The fallback is **one process,
two OpenMP threads**, with the original physical hierarchy and strict floating
flags, rather than the cluster's 224-rank decomposition. A serial success would
not disprove a distributed layout/communication defect. The local compiler is
GCC 16 on arm64, not the cluster's compiler/architecture; this difference is
explicitly registered in [t24-build-registration.json](t24-build-registration.json).

All four Chombo libraries and all native units are compiled with Float64,
`-fno-fast-math -ffp-contract=off`; the diagnostic build also has debug symbols
and frame pointers. [t24-stack.cpp](t24-stack.cpp) is a default-off linked
observer; it pre-opens per-thread files and pre-warms the native unwinder.
SIGUSR1 requests self-backtraces from both registered OpenMP threads. This is
best-effort signal-time unwinding, not a debugger backtrace. A two-thread stack
self-test and a one-step production-main regression compare the plain binary
against observers omitted/off/on before the production-sized restart is
allowed. All defined evolved values and valid/same-level-supported diagnostic
values must be bit identical at t=0 and one step; the known undefined saved
plot diagnostic halos are excluded exactly as in T23. No pass is claimed until
the per-job native receipts and comparison CSV exist.

The restart copies only the small numerical puncture history into a writable
private directory. The sealed HDF5 checkpoint is read directly from its
collected absolute path. Both static file paths are absent. The only parameter
overrides are restart/output paths and the target step **156**; the production
checkpoint interval **4**, disabled plots, physical diagnostics, hierarchy,
time steps and regrid intervals are retained. The target spans the stalled
time and ends at **136.5 M_i**. It does not purport to reach coincidence.

The hard allowed cap is **6 GB peak RSS per process**. The first serial attempt
was stopped at **21.75 s**, at level-0 advance from t=133: measured RSS was
**3.746 GB**, physical footprint **5.803 GB**. The conservative combined gate
thus stopped below the requested RSS cap; this is not the cluster stall. The
partial attempt and all receipts remain under `/private/tmp/ems-t24/restart152/`.
The new, separately detached attempt uses a preventive **5.7 GB RSS** gate
and a **7.5 GB aggregate footprint** gate. It does not pass through an outer
runner that would reapply the combined per-process footprint gate. Its direct
owner writes native wait4 and sampled RSS/footprint receipts and a done marker.
The physical parameters, executable and loaded state are unchanged.
The original queue enforced a **16 GB** output ceiling. The direct retry keeps
plots disabled and writes only its final scheduled checkpoint plus logs/layouts;
its storage budget remains **16 GB** and does not have that outer disk gate.
The watchdog requests stacks after **300**, **600** and **900 seconds** with no
log growth, then terminates that diagnostic attempt. A separate **4-hour**
wall cap prevents an unbounded local replay. Either gate is reported as a
controlled stop, never as a numerical failure. Completed measured resources
are listed in [t24-resources.csv](t24-resources.csv); the largest currently
recorded per-process RSS is **{peak/1e9:.3f} GB**.

Active done marker: `{active/'done.exit'}`.
Native restart receipt: `{native_directory/'native.receipt.json'}`.
Observer qualification: {qualified['status'] if qualified else '**pending**'}.
Native restart result: {native['gate_reason'] or ('completed to step 156' if native['complete_to_step156'] else 'native error') if native else '**pending**'}.
The outer launcher exit alone is never acceptance.

## Fix scope and merger policy

There is **no qualified fix** yet. A fix must be selected from captured stacks
and layout evidence, then reproduce unchanged evolved bits through steps
152–154 and advance through the affected time and a declared coincidence test.
Those comparisons are outstanding. If the cause is in a gauge, equation, KO
or transfer rule, work stops for the controller as ordered. If a tagging or
centre-transition fix changes the refinement footprint, that change requires
a new run review; it must not be presented as an identity-preserving restart.

For exp-0025, the two new hooks are bypassed when `ems_binary_refinement=false`
and `ems_track_punctures=false`; the two-centre refinement and tracker paths
then cannot themselves execute. Generic Chombo regrid, sparse sets, exchange
and point-transfer paths are shared. Until the cause is identified, single-hole
immunity is not established.

The following design items are **inferred requirements**, not measured cures.
Keep the existing union and fixed labels through overlap until a numerical
merger/coincidence transition is specified and audited. If one common centre
is to replace two refinement centres, predeclare its definition and transition
condition and review the resulting footprint; separation alone is not a
qualified horizon decision. Use current numerical centres to seed the common
finder, keep the inherited strict N48/N96 qualification before horizon-based
losses/clock start, and treat individual warm tracking according to the existing
APPROXIMATE/UNRESOLVED policy. Add explicit layout/phase evidence at the
transition, so per-rank quiet progress is not misread as an exact call site.
No new merger threshold, finder tolerance or refinement radius is invented here.

The manifest [t24-manifest.txt](t24-manifest.txt) seals the current local packet;
it records pending evidence as pending. Refresh with `t24-local.py packet`
after checking each own receipt. No SSH, cluster operation or commit was made.
'''
    (HERE/'t24-merger-stall.md').write_text(report)
    readme = HERE/'README.md'
    title = '## T24 — binary merger stall, d507438'
    text = readme.read_text()
    if title in text:
        start = text.index(title)
        end = text.find('\n## ',start+len(title))
        text = text[:start].rstrip()+('' if end<0 else '\n\n'+text[end+1:])
    text += '\n\n'+title+'\n\n'+('**READY-EXCEPT; cause and fix unestablished.** '
          '[Source audit and bounded replay](t24-merger-stall.md), '
          '[driver](t24-local.py), [input audit](t24-input-audit.json), '
          '[nesting](t24-checkpoint-nesting.csv), [LB counts](t24-load-balance-counts.csv), '
          '[status](t24-status.json), [manifest](t24-manifest.txt). '
          'All 1,855 small-file hashes and the step-152 numerical seal match; '
          '36 saved parent-support checks pass. Two-centre tagging is a finite union; '
          'the quiet-rank progress policy does not locate the stall. LLDB and local MPI '
          'cannot launch usable debugger/multi-rank probes in this sandbox. The strict '
          'serial/two-thread observer regression and unchanged-grid restart to step 156 '
          'are in the bounded evidence package. Observer regression now passes; the '
          'first serial restart was stopped by a conservative footprint gate at '
          '3.746 GB RSS. An RSS-capped retry is detached. The first observer build required '
          'one missing include and explicit Chombo dependency targets; failed receipts are retained. Markers: '
          '`/private/tmp/ems-t24/restart152-rss/done.exit`. No production patch or refinement footprint '
          'change is made; pre-stall fix identity and coincidence evolution remain pending.\n')
    readme.write_text(text)
    paths = [p for p in HERE.glob('t24-*') if p.is_file() and p.name!='t24-manifest.txt']
    paths += [HERE/'README.md',HERE/'t13-run.py',HERE/'t23-compare.py',
              ROOT/'pipeline/plan.json',ROOT/'continuation/plan.json',ROOT/'continuation-2/plan.json',
              ROOT/'continuation-3/plan.json',ROOT/'restart152-rss/launch.json',
              ROOT/'source/Examples/EMS/EMSBH2DLevel.cpp',
              ROOT/'source/Examples/EMS/Main_EMSBH2DBH.cpp',
              ROOT/'source/Source/GRChomboCore/GRAMRLevel.cpp',
              ROOT/'source/Source/BlackHoles/PunctureTracker.cpp',
              ROOT/'source/Source/GRChomboCore/PointAMRTransfer.hpp']
    (HERE/'t24-manifest.txt').write_text('# T24 provisional source audit / bounded local replay; no commit\n'
          '# Cause/fix UNESTABLISHED; running receipts are not sealed as completed\n'
          '# SHA256 bytes absolute_path\n'+''.join(f'{sha(p)}  {p.stat().st_size}  {p}\n' for p in sorted(set(paths))))
    print('T24_PACKET_REFRESHED',status['status'],flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['mpi_probe','debugger_probe','sparse_probe','input_audit','prepare','build','regression','restart','packet','replay_worker','detach_replay'])
    args = parser.parse_args()
    if args.mode in ('restart','detach_replay') and (HERE/'t24-offline-receipt.json').exists():
        parser.error('Full local replay is retired after its footprint-gate stop; use the T24 diagnostic restart design.')
    sys.exit(globals()[args.mode]())
