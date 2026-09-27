#!/usr/bin/env python3
"""Summarize the fixture and fixed-mask t=0 runs; no fitted convergence floor."""
import csv
import io
import math
from pathlib import Path

root = Path(__file__).parent
results = root / 'results'
regions = ('left', 'right', 'middle', 'outer', 'left_axis', 'right_axis', 'middle_axis')
fields = ('H', 'Mx', 'My', 'M', 'GaussE', 'GaussB')
levels = ('0.5', '1', '2', '4')
rows = []
times = []
for mode in ('ctt', 'density'):
    for level in levels:
        path = results / f'{mode}-{level}.log'
        assert (results / f'{mode}-{level}.exit').read_text().strip() == '0', path
        lines = path.read_text().splitlines()
        header = next(line for line in lines if line.startswith('mode,level,'))
        data = [line for line in lines if line.startswith(('CTT,', 'density,'))]
        parsed = list(csv.DictReader(io.StringIO('\n'.join([header] + data))))
        assert {r['region'] for r in parsed} == set(regions), path
        for row in parsed:
            assert int(row['count']) > 0 and int(row['chi_floor']) == int(row['lapse_floor']) == 0
            assert float(row['field_error']) < 5e-13
            for name in fields:
                for norm in ('L2', 'Linf'):
                    assert math.isfinite(float(row[f'{name}_{norm}']))
        rows += parsed
        times.append((mode, level, next(line[8:] for line in lines if line.startswith('seconds='))))
lookup = {(r['mode'], r['region'], float(r['level'])): r for r in rows}
for r in rows:
    previous = lookup.get((r['mode'], r['region'], float(r['level']) / 2))
    for name in fields:
        for norm in ('L2', 'Linf'):
            key = f'{name}_{norm}'
            value = float(r[key])
            r[key + '_order'] = '' if not previous or value == 0 or float(previous[key]) == 0 else math.log2(float(previous[key]) / value)
with (root / 'convergence.csv').open('w', newline='') as out:
    writer = csv.DictWriter(out, list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

baseline = (results / 'baseline-fixtures.log').read_bytes()
current = (results / 'current-fixtures.log').read_bytes()
assert baseline == current, 'no-companion fixture outputs changed'
assert baseline.count(b'all-output fingerprint=') == 13

parity = (results / 'parity.log').read_text().splitlines()
assert any(line.startswith('PASS seconds=') for line in parity)
parity_rows = list(csv.DictReader(io.StringIO('\n'.join(line for line in parity if line.startswith(('table,', 'physical,', 'ccz4,'))))))
assert all(float(row['max_scaled_error']) <= 5e-13 for row in parity_rows)

out = ['# Measured C++ handoff results', '',
       'Decision: fixture parity and no-companion bit identity pass. The CTT constraints show fourth-order decay before the declared numerical floors; the density seed retains its nonzero continuum residual. No evolution or HPC claim is made.', '',
       'Errors use `abs(C++ − Julia)/(1 + abs(Julia)) ≤ 5e-13`, the existing density-binary fixture tolerance. All 42 points pass in each table.', '',
       '| Table | Group | Max scaled error | Worst point | Component |',
       '|---|---|---:|---:|---|']
for r in parity_rows:
    out.append(f"| {r['table']} | {r['group']} | {float(r['max_scaled_error']):.4e} | {r['point']} | {r['field']} |")
out += ['', 'Each entry below is unweighted cell RMS (L2), with `log2(previous/current)` in parentheses after the first row. '
        '`M = hypot(Mx, My)`. GaussB is exactly zero at every level; its order is undefined. '
        'Full-precision L2 and L∞ values, both momentum components, orders, counts, staging errors and the reflection diagnostic are in [convergence.csv](convergence.csv).', '',
        'The core spacings are 1/16, 1/32, 1/64, 1/128 M; the outer spacings are 8, 4, 2, 1 M. '
        'The 1/16 level was added to expose the asymptotic inter-hole H range before the finest-grid numerical floor. All masks and fields are fixed across the final exterior ladder.']
for region in regions:
    out += ['', f'## {region}', '', '| Mode | h/M | Cells | H L2 (p) | Mx L2 (p) | My L2 (p) | M L2 (p) | GaussE L2 (p) | GaussB L2 |',
            '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for level in levels:
        for mode in ('CTT', 'density'):
            r = lookup[mode, region, float(level)]
            values = []
            for name in fields:
                p = r[f'{name}_L2_order']
                values.append(f"{float(r[f'{name}_L2']):.4e}" + (f' ({p:.3f})' if p != '' else ''))
            out.append(f"| {mode} | {float(r['dx']):.7g} | {r['count']} | " + ' | '.join(values) + ' |')
out += ['', '## Inter-hole Hamiltonian reflection diagnostic', '',
        'The symmetric head-on configuration has an even Hamiltonian. The antisymmetric RMS is computed as '
        '`RMS((H(x,y) − H(−x,y))/2)` using the same grid. It measures symmetry-breaking evaluation noise; it is not subtracted from any reported norm.', '',
        '| h/M | H RMS | Antisymmetric RMS |', '|---:|---:|---:|']
for level in levels:
    r = lookup['CTT', 'middle', float(level)]
    out.append(f"| {r['dx']} | {float(r['H_L2']):.4e} | {float(r['H_antisym_L2']):.4e} |")
out += ['', 'Inter-hole H has orders 3.997 and 3.995 on h=1/16, 1/32, 1/64. At h=1/128 its RMS is 5.67e-11 and the antisymmetric RMS is 2.38e-11. The latter grows approximately fourfold per halving of h, consistent with Float64 field-evaluation noise amplified by second differences. The raw finest order 1.802 is retained; no residual is subtracted or threshold relaxed.', '', 'The collars are 0.75 ≤ R_rest ≤ 1.5 M (R_h=0.63593977642346233 M in the supplied profile). Inter-hole and outer masks and the fixed axis strips are defined in [README.md](README.md). The earlier 0.5-inner-radius collar control is archived separately. This ladder does not resolve the companion’s separate scaled D¹C interface plateau; that declared limitation is retained.', '', '## Timings', '', '| Mode | Refinement | Program seconds |', '|---|---:|---:|']
for mode, level, seconds in times:
    out.append(f'| {mode} | {level} | {float(seconds):.3f} |')
out += ['', '## Bit identity', '', 'Baseline and current fixture stdout compare byte-for-byte equal. The 13 fingerprints cover every computed numeric value in the reference, density-binary and B/E echo fixture tables. Frozen baseline was built from e1ee1b8 with only fingerprint printing added.', '',
        '```text']
for line in baseline.decode().splitlines():
    if 'all-output fingerprint=' in line or line.startswith('Tests/EMSTrumpet/'):
        out.append(line)
out += ['```', '']
(root / 'results.md').write_text('\n'.join(out))
print(f'{len(rows)} mask rows, {len(parity_rows)} parity groups; 13 unchanged bit fingerprints; report written')
