#!/usr/bin/env python3
"""Reproducible T18 local setup, using existing build and resource helpers."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = Path('/private/tmp/ems-t18')
PY = '/Users/auroradysis/miniconda3/bin/python'
T17 = REPO.parent/'wt-native-t4/Tests/EMSNative'

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

def measured(name, command):
    directory = ROOT/name
    directory.mkdir(parents=True, exist_ok=True)
    with (directory/'run.log').open('w') as log:
        subprocess.run([PY, str(HERE/'t18-run.py'), '--measure', name,
                        '--directory', str(directory), '--', *map(str, command)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)

def build(which):
    ch = Path('/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
    if which == 'census':
        command = ['/opt/homebrew/bin/g++-16', '-O3', '-std=c++17', '-fopenmp']
        command += ['-D'+s for s in ('CH_SPACEDIM=2', 'CH_Darwin', 'NDEBUG',
                                     'CH_USE_64', 'CH_USE_DOUBLE', 'CH_LANG_CC')]
        command += ['-I'+str(ch/'src'/s) for s in ('BoxTools', 'BaseTools')]
        command += [str(HERE/'t18-census.cpp'), '-L'+str(ch)]
        command += ['-l'+s+'2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC'
                    for s in ('boxtools', 'basetools')]
        command += ['-o', str(ROOT/'census.ex')]
    else:
        builder = module('t18builder', HERE/'t14-prepare.py')
        builder.ROOT = ROOT
        builder.build('t18-native.cpp', 'native-window.ex' if which=='window' else 'native.ex')
        return
    (ROOT/'census-build-command.json').write_text(json.dumps(command, indent=2)+'\n')
    subprocess.run(command, check=True)

def prepare():
    replace = module('t18params', HERE/'t13-prepare.py').replace
    base = (T17/'t17-coarse-geometric.txt').read_text()
    jobs = []
    # Full-depth, narrow fixed AMR patches: identical coverage across splits.
    # Three synchronized coarse steps exercise all 13 levels and point transfers.
    cases = [('neutral-b128-t1', 128, 1), ('neutral-b16-t1', 16, 1),
             ('neutral-b24-t1', 24, 1), ('neutral-b32-t1', 32, 1),
             ('neutral-b48-t1', 48, 1), ('neutral-b16-t2', 16, 2),
             ('neutral-b16-t4', 16, 4)]
    for name, size, threads in cases:
        directory = ROOT/'local'/name
        for sub in ('chk', 'plt'):
            (directory/sub).mkdir(parents=True, exist_ok=True)
        text = replace(base, dict(N1=128, N2=64, L=224, center='112 0',
            star_centre='112 0', mass_extraction_center='112 0',
            max_box_size=size, max_steps=3, stop_time=2,
            checkpoint_interval=3, plot_interval=1, t18_fixed_hierarchy='true',
            ems_binary_refinement='false', ems_track_punctures='false',
            t18_monitor='false', t18_stop_time=0))
        parameter = HERE/(name.replace('neutral', 't18-neutral')+'.txt')
        parameter.write_text(text)
        jobs.append(dict(name=name, directory=str(directory),
            command=['env', 'OMP_NUM_THREADS='+str(threads), str(ROOT/'native.ex'), str(parameter)]))
    # Same h0, finest h, phase and per-hole native tag radii as the coarse draft.
    for dt in (.25, .375, .5):
        name = 'dt-'+str(dt)
        directory = ROOT/'local'/name
        for sub in ('chk', 'plt'):
            (directory/sub).mkdir(parents=True, exist_ok=True)
        text = replace(base, dict(N1=64, N2=32, L=112, center='40 0',
            star_centre='40 0', mass_extraction_center='40 0',
            binary='false', ems_binary_refinement='false', ems_ctt_data_path='',
            max_box_size=128, max_steps=1, stop_time=1, dt_multiplier=dt,
            checkpoint_interval=-1, plot_interval=-1,
            t18_fixed_hierarchy='false', t18_monitor='true',
            t18_stop_time=.1025390625,
            mass_extraction_radii=' '.join(format((224/2**i)/1.2, '.17g') for i in range(1,13))))
        text = '\n'.join(line for line in text.splitlines()
                         if not line.startswith('ems_ctt_data_path ='))+'\n'
        parameter = HERE/('t18-'+name+'.txt')
        parameter.write_text(text)
        jobs.append(dict(name=name, directory=str(directory),
            command=['env', 'OMP_NUM_THREADS=2', str(ROOT/'native.ex'), str(parameter)]))
    plan = ROOT/'local/plan.json'
    plan.write_text(json.dumps(dict(jobs=jobs), indent=2)+'\n')
    registration = dict(baseline='77c5b6f', process_RSS_cap_bytes=4000000000,
        maximum_threads=4, neutrality_steps=3, neutrality_levels=list(range(13)),
        dt_window_M=.1025390625, dt_window_Rh=16.975394601194597,
        bound_before_results='All fields finite; no chi/lapse floor crossings; '
            'positive metric; puncture max |Gamma|, |K|, |Theta|, |Pi| < 1e6; '
            'max |shift| < 2, max lapse < 2. These are a blow-up screen, '
            'not a long-time or temporal-convergence proof.',
        jobs=jobs)
    (HERE/'t18-registration.json').write_text(json.dumps(registration, indent=2)+'\n')
    print(plan)

if __name__ == '__main__':
    ROOT.mkdir(exist_ok=True)
    if sys.argv[1] == 'build':
        build(sys.argv[2])
    elif sys.argv[1] == 'census':
        measured('census', [ROOT/'census.ex', HERE/'t18-census'])
    elif sys.argv[1] == 'prepare':
        prepare()
    else:
        raise SystemExit('use build census|native, census, or prepare')
