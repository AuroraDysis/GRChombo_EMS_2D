#!/usr/bin/env python3
"""Detached serial T25 continuation; each phase and native probe has a receipt.

The atomic overall done.exit records workflow success, never finder success.
The full time table preserves capped/weak/unresolved results explicitly.
At most 18 physical searches, 235 s native cap each; no finder on midpoint.
"""
import csv,importlib.util,json,resource,subprocess,sys,time,traceback
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('early',HERE/'t25-early.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
PY='/Users/auroradysis/miniconda3/bin/python'
def report():
    # A completed post-flow map has the authoritative final report. Preserve
    # its provenance rather than replacing it with the preliminary seed table.
    if (HERE/'t25-early-final-audit.json').exists():
        e.load_module('t25-early-finalize').main()
        return
    rows=e.read(HERE/'t25-early-time-table.csv')
    complete=(e.ROOT/'finders-complete.json').exists()
    text=['## Time-resolved exp-0024 map, t = 0 to 28 M_i',
          '', '**Completed offline reduction.**' if complete else '**PENDING bounded offline analysis; no evolution.**', '',
          'The time-zero state is regenerated with exp-0024\'s d=16, rapidity '
          '0.05778205303580913 and its sealed companion (SHA-256 '
          '`107370ae00d8b69dc3122dc023e34bbff23b90afc8401bb233743c1182e6083c`). '
          'It replaces the d=32 T17 control only in this time series. The original '
          'control and its results above remain separate. No static or CTT file '
          'is accessible in any positive-time map or search.', '',
          'The nested radial input starts at three finest cells and ends at '
          '0.05 M_i, with 4% initial sphere spacing (8% for the broader fallback). Sign boundaries are bisected up to '
          'four times, to 0.5% local spacing. The irreducible angular zero band '
          'can leave a wider pointwise barrier; that width is reported, not '
          'declared artificially small. Where a centred sphere lacks a barrier, '
          'the map includes offsets −4, −2, 0, 2, 4 finest cells and axial/cylindrical '
          'aspect ratios 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 2. Surfaces below three local cells '
          'are excluded. N96 and N192 are retained, and selected endpoint signs '
          'are independently checked at N96.', '',
          'One best pointwise negative-to-positive pair is chosen per hole and '
          'time, using width and then endpoint angular RMS after the shape '
          'refinement, without an area target. In its absence, one mean-only negative-to-positive pair can '
          'seed the bounded search through the explicit Tests-only parameter '
          '`map_find_allow_average=true` (default false). This follows the '
          'continuation\'s request to search from each bracket; it does not '
          'promote a mean sign flip into a pointwise barrier. If no such pair '
          'exists, no finder runs. All selected searches use unchanged N48, '
          'chase=quota=1, 235 s, update cap 100000, inert floor window 100001, '
          'and stages 1e-7/1e-10/1e-12.', '',
          'The linear root is a trial surface obtained from area-weighted mean '
          'expansion at the two endpoints. Even for a pointwise barrier this '
          'single radius is not an angularly solved MOTS. The CSV reports its '
          'signed range/RMS and its A, coupling-weighted Q, lapse, χ, K and '
          'K−2Θ. Horizon claims and gauge-invariant loss measurements require '
          'the actual finder to qualify; a TIME_CAP remains unqualified.', '',
          '| t / M_i | hole | bracket | a / M_i | min extent / h | A (trial) | Q (trial) | lapse | finder squared residual | status |',
          '|---:|---:|:---|---:|---:|---:|---:|---:|---:|:---|']
    def fmt(x):return f'{float(x):.7g}' if x!='' else '—'
    for r in rows:text.append('| '+' | '.join([fmt(r['time']),r['hole'],r['bracket_status'],fmt(r['r_cyl']),fmt(r['cells']),fmt(r['A_trial']),fmt(r['Q_trial']),fmt(r['lapse_mean']),fmt(r['finder_squared']),r['finder_status']])+' |')
    text+=['', 'The complete machine-readable table is '
        '[t25-early-time-table.csv](t25-early-time-table.csv). Per-checkpoint '
        '`t25-early-step*.csv` packet files contain centred spheres and the '
        'selected families; all mapped families and every angular row remain '
        'under `/private/tmp/ems-t25/early/step*/`. Bracket/barrier CSVs retain '
        'all families. No bulk checkpoint is copied to the worktree. The '
        'midpoint sanity map at t=28 is '
        '[t25-early-midpoint-t28.csv](t25-early-midpoint-t28.csv); it uses '
        'enclosing spheres and prolate surfaces with axial extents derived '
        'from the measured separation and has no finder.', '',
        'Every native receipt must have no gate reason, the expected child '
        'return code, finite stage data and the native completion marker. '
        'Workflow exit zero is not a finder qualification. '
        '[t25-early-receipts.csv](t25-early-receipts.csv) gives the measured '
        'per-process RSS and cost. The initial-data process has a 6 GB cap; '
        'maps retain the tighter 3 GB cap; finders have the task\'s 6 GB cap. All use one thread. '
        '[t25-early-regression.csv](t25-early-regression.csv) records 55,944 '
        'bit-identical existing values and both weaker-seed gate checks.', '']
    if complete:
        strong=[r for r in rows if r['bracket_status']=='POINTWISE_NEGATIVE_TO_POSITIVE']
        weak=[r for r in rows if r['bracket_status']=='AVERAGE_ONLY_WEAKER']
        absent=[r for r in rows if r['bracket_status']=='NO_RESOLVED_BRACKET']
        passed=[r for r in rows if r['finder_status']=='FOUND']
        text+=['### Scope of the time-history statement', '',
            f'The 18 hole/time rows contain {len(strong)} sampled pointwise barriers, '
            f'{len(weak)} weaker mean-only crossings and {len(absent)} rows without '
            f'a resolved bracket in the sampled families. {len(passed)} N48 '
            'searches pass all three stages. This count is not a new angular '
            'N48/N96 strict-finder pair. The question of coordinate squeezing '
            'must be read only from supported roots; below-three-cell rows '
            'and absent brackets provide no measured horizon radius.', '',
            'The table lists A and Q of the interpolated trial surface, '
            'not a certified horizon when its finder is capped. Therefore '
            'those rows alone cannot establish preservation or loss of '
            'the progenitor A=0.2245 and Q=1.0693. The finder A and Q columns '
            'are also trial reductions unless all three stages pass. Absence '
            'of a bracket in this finite sphere/spheroid family does not '
            'prove absence of an apparent horizon.', '']
    p=HERE/'t25-expansion-map.md';old=p.read_text();heading='## Time-resolved exp-0024 map, t = 0 to 28 M_i'
    old=old.split(heading)[0].rstrip()
    p.write_text(old+'\n\n'+'\n'.join(text)+'\n')
    p=HERE/'README.md';old=p.read_text();heading='### T25 early-checkpoint continuation'
    old=old.split(heading)[0].rstrip()
    p.write_text(old+'\n\n'+heading+'\n\n'+
        ('The bounded offline pipeline is complete. ' if complete else 'The bounded offline pipeline is pending. ')+
        'The matching d=16 time-zero state uses the exact sealed exp-0024 inputs; '
        'the earlier d=32 control remains separate. The native map adds area-weighted '
        'K and K−2Θ without changing 55,944 existing values. '
        'The pipeline SHA-256 checks all 8 collected checkpoints and processes them serially, '
        'one at a time. See [t25-early-time-table.csv](t25-early-time-table.csv), '
        'per-checkpoint t25-early-step*.csv, '
        '[t25-early-receipts.csv](t25-early-receipts.csv) and the new section in '
        '[t25-expansion-map.md](t25-expansion-map.md). Mean-only crossings remain '
        'weaker, unresolved extents are excluded, and one bounded N48 search '
        'runs per selected bracket with the unchanged 235 s cap. '
        'Midpoint t=28 is map-only. Completion marker: '
        '`/private/tmp/ems-t25/early/pipeline/done.exit`. No Source/Examples '
        'edit, evolution, SSH or commit is part of this continuation.\n')
def main():
    d=e.ROOT/'pipeline';d.mkdir(parents=True,exist_ok=True);start=time.monotonic();rc=1
    try:
        assert json.loads((HERE/'t25-early-regression.json').read_text())['status']=='PASS'
        import hashlib
        assert hashlib.sha256(e.EXE.read_bytes()).hexdigest()==json.loads((HERE/'t25-early-regression.json').read_text())['executable_sha256']
        # The matching t=0 regeneration may already be running independently.
        # Never launch a duplicate checkpoint process while it owns the inputs.
        job=e.ROOT/'initialization-job'
        if job.exists():
            deadline=time.monotonic()+3000;interval=15
            while not (job/'done.exit').exists():
                if time.monotonic()>deadline:raise RuntimeError('t0 dependency exceeded 3000 s')
                print('WAIT_T0_DEPENDENCY',str(job/'done.exit'),flush=True)
                time.sleep(interval);interval=min(60,interval*2)
            if int((job/'done.exit').read_text()):raise RuntimeError('t0 dependency failed')
        for phase in ('initialize','maps','finders'):
            pd=d/phase;pd.mkdir(exist_ok=True)
            with (pd/'run.log').open('wb') as log:
                result=subprocess.run([PY,str(HERE/'t25-early.py'),phase],stdout=log,stderr=subprocess.STDOUT)
            (pd/'done.exit.tmp').write_text(str(result.returncode)+'\n');(pd/'done.exit.tmp').replace(pd/'done.exit')
            print('PHASE_DONE',phase,result.returncode,flush=True)
            if result.returncode:raise RuntimeError('failed phase '+phase)
            e.summarize();report()
        assert len(e.read(HERE/'t25-early-time-table.csv'))==18
        rc=0
    except Exception:
        traceback.print_exc()
    finally:
        e.atomic_json(d/'receipt.json',dict(returncode=rc,wall_seconds=time.monotonic()-start,
            peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            no_advances=True,serial=True,physical_search_seconds_cap=235))
        (d/'done.exit.tmp').write_text(str(rc)+'\n');(d/'done.exit.tmp').replace(d/'done.exit')
        # Seal only after all writers in the serial pipeline have returned.
        subprocess.run([PY,str(HERE/'t25-manifest.py')],check=True)
    return rc
if __name__=='__main__':sys.exit(main())
