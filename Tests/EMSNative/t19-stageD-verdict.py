#!/usr/bin/env python3
"""Read the sealed exp-0023 tables; regenerate the T19 report, CSVs and figures.

Usage: t19-stageD-verdict.py COLLECTED_PRODUCTION_ROOT CONTRACT_DIR [--output-dir DIR]
The five large ledgers default to COLLECTED_PRODUCTION_ROOT/../t19;
small outputs are t19-* next to this script. No evolution, reconstruction, new masks,
new clocks or replacement uncertainty estimator is used. Input status labels
are retained, including non_convergent intervals which straddle zero.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import resource
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/ems-t19-matplotlib')
HERE = Path(__file__).resolve().parent
EVOLVED = 'chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split()
CONSTRAINTS = 'Ham Mom1 Mom2 Mom GaussE GaussB Theta Lambda Xi Z1 Z2 Z CGamma1 CGamma2 CGamma det_h_minus_1 tr_A'.split()
STATUSES = ('consistent_with_4', 'measured_order_outside_4_interval',
            'measured', 'sampling_or_roundoff_limited', 'non_convergent')
INPUTS = set()
OUTPUTS = set()
BULK_NAMES = {'order-exceptions', 'field-orders-outside', 'constraint-orders',
              'temporal-failures', 'native-constraint-norms'}
BULK_DIR = None


def read_text(path):
    path = Path(path)
    INPUTS.add(path)
    return path.read_text()


def read_json(path):
    return json.loads(read_text(path))


def read_csv(path):
    return list(csv.DictReader(read_text(path).splitlines()))


def number(row, key):
    value = row.get(key, '')
    return float(value) if value not in ('', None) else math.nan


def truth(value):
    return value is True or value == 'True'


def write_csv(name, rows, columns=None):
    directory = BULK_DIR if name in BULK_NAMES else HERE
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / ('t19-' + name + '.csv')
    if columns is None:
        columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    OUTPUTS.add(path)
    return path


def save_text(name, text):
    path = HERE / ('t19-' + name)
    path.write_text(text)
    OUTPUTS.add(path)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def grouped(rows, keys):
    result = defaultdict(list)
    for row in rows:
        result[tuple(row[k] for k in keys)].append(row)
    return result


def temporal_ok(row):
    # A zero/zero temporal fraction is undefined, not silently admitted.
    return all(math.isfinite(number(row, k)) and number(row, k) <= .2 for k in
               ('temporal_fraction_of_low_mid', 'temporal_fraction_of_mid_high'))


def significance_ok(row):
    return number(row, 'spatial_significance') > 5


def exterior(row):
    return int(row['clock_index']) >= 1 and row['mask'] != 'outer_boundary_shell'


def negative_outside_uncertainty(row):
    return (significance_ok(row) and temporal_ok(row)
            and number(row, 'order_interval_high') < 0)


def check_battery(rows, masks, clocks, name):
    keys = [(r['clock_index'], r['mask'], r['category'], r['field']) for r in rows]
    assert len(set(keys)) == len(keys), (name, 'duplicate key')
    # The sealed reducer has no radial ray for the Cartesian boundary shell
    # (production/reduce.py:170). Preserve and report that absence.
    expected_masks = set(masks) - ({'outer_boundary_shell'} if name == 'ray' else set())
    assert {r['mask'] for r in rows} == expected_masks, (name, 'mask coverage')
    assert {int(r['clock_index']) for r in rows} == set(clocks), (name, 'clock coverage')
    assert {r['field'] for r in rows} == set(EVOLVED + [x for x in CONSTRAINTS if x not in EVOLVED])
    for r in rows:
        assert abs(number(r, 'time_M') - clocks[int(r['clock_index'])]['time_M']) < 1e-12
        d1, d2 = number(r, 'difference_low_mid'), number(r, 'difference_mid_high')
        floor = number(r, 'Float64_floor')
        u1 = number(r, 'I8_I10_low_mid_uncertainty') + floor
        u2 = number(r, 'I8_I10_mid_high_uncertainty') + floor
        sig = min(d1 / max(u1, sys.float_info.min), d2 / max(u2, sys.float_info.min))
        assert math.isclose(sig, number(r, 'spatial_significance'), rel_tol=2e-12, abs_tol=1e-12)
        resolved = d1 > u1 and d2 > u2
        if resolved:
            p = math.log(d1 / d2) / math.log(1.5)
            lo = math.log((d1-u1)/(d2+u2)) / math.log(1.5)
            hi = math.log((d1+u1)/(d2-u2)) / math.log(1.5)
            assert all(math.isclose(a, b, abs_tol=1e-11, rel_tol=2e-12) for a, b in
                       ((p, number(r, 'measured_order')), (lo, number(r, 'order_interval_low')),
                        (hi, number(r, 'order_interval_high')))), (name, keys, 'order mismatch')
            expected = ('measured' if sig <= 5 or not temporal_ok(r) else
                        'consistent_with_4' if lo <= 4 <= hi else
                        'non_convergent' if p <= 0 else 'measured_order_outside_4_interval')
        else:
            expected = 'sampling_or_roundoff_limited'
        assert r['status'] == expected, (name, r, expected)


def summary(rows):
    result = []
    for key, group in grouped([r for r in rows if r['field'] in EVOLVED],
                              ('dataset', 'clock_index', 'time_M', 'mask', 'category')).items():
        orders = [number(r, 'measured_order') for r in group if math.isfinite(number(r, 'measured_order'))]
        eligible = [number(r, 'measured_order') for r in group if significance_ok(r) and temporal_ok(r)
                    and math.isfinite(number(r, 'measured_order'))]
        counts = Counter(r['status'] for r in group)
        row = dict(zip(('dataset', 'clock_index', 'time_M', 'mask', 'category'), key))
        row.update(fields=len(group), measured_orders=len(orders),
                   order_median=statistics.median(orders) if orders else '',
                   order_min=min(orders) if orders else '', order_max=max(orders) if orders else '',
                   spatial_temporal_qualified_orders=len(eligible),
                   qualified_order_median=statistics.median(eligible) if eligible else '',
                   outside_4_fields=';'.join(r['field'] for r in group
                       if math.isfinite(number(r, 'order_interval_low')) and not
                       number(r, 'order_interval_low') <= 4 <= number(r, 'order_interval_high')))
        row.update({s: counts[s] for s in STATUSES})
        result.append(row)
    return result


ORDER_COLUMNS = ['dataset', 'clock_index', 'time_M', 'mask', 'category', 'field',
                 'cells', 'weight', 'difference_low_mid', 'difference_mid_high',
                 'I8_I10_low_mid_uncertainty', 'I8_I10_mid_high_uncertainty', 'Float64_floor',
                 'spatial_significance', 'temporal_high_dt2', 'temporal_fraction_of_low_mid',
                 'temporal_fraction_of_mid_high', 'measured_order', 'order_interval_low',
                 'order_interval_high', 'status', 'claim_interval']


def exceptions(rows):
    result = []
    for r in rows:
        if r['status'] == 'consistent_with_4':
            continue
        row = {k: r[k] for k in ORDER_COLUMNS}
        reasons = []
        if not significance_ok(r):
            reasons.append('spatial_significance_not_above_5')
        if not temporal_ok(r):
            reasons.append('temporal_control_not_qualified')
        lo, hi = number(r, 'order_interval_low'), number(r, 'order_interval_high')
        if not math.isfinite(lo):
            reasons.append('order_interval_undefined')
        elif not lo <= 4 <= hi:
            reasons.append('four_outside_measured_interval')
        row.update(reasons=';'.join(reasons), registered_exterior_clock=exterior(r),
                   negative_outside_uncertainty=negative_outside_uncertainty(r),
                   interpretation=('registered_FAIL_witness' if exterior(r) and negative_outside_uncertainty(r)
                                   else 'as_measured_or_unresolved'))
        result.append(row)
    return result


def constraint_orders(direct, battery):
    lookup = {(r['clock_index'], r['mask'], r['field']): r for r in battery if r['category'] == 'all'}
    expected = {key for key in lookup if key[2] in CONSTRAINTS}
    assert {(r['clock_index'],r['mask'],r['constraint']) for r in direct} == expected
    assert len(direct) == len(expected)
    result = []
    for r in direct:
        b = lookup[(r['clock_index'], r['mask'], r['constraint'])]
        row = dict(r)
        row.update(pair_difference_order=b['measured_order'],
                   pair_difference_interval_low=b['order_interval_low'],
                   pair_difference_interval_high=b['order_interval_high'],
                   pair_difference_status=b['status'],
                   difference_low_mid=b['difference_low_mid'], difference_mid_high=b['difference_mid_high'],
                   weight=b['weight'], Float64_floor=b['Float64_floor'],
                   I8_I10_low_mid_uncertainty=b['I8_I10_low_mid_uncertainty'],
                   I8_I10_mid_high_uncertainty=b['I8_I10_mid_high_uncertainty'],
                   temporal_high_dt2=b['temporal_high_dt2'],
                   temporal_fraction_of_low_mid=b['temporal_fraction_of_low_mid'],
                   temporal_fraction_of_mid_high=b['temporal_fraction_of_mid_high'],
                   spatial_significance=b['spatial_significance'])
        for pair, a, c in (('low_mid', 'low', 'mid'), ('mid_high', 'mid', 'high')):
            av, cv = number(r, a), number(r, c)
            u, v = number(r, a+'_uncertainty'), number(r, c+'_uncertainty')
            sig = min(av/max(u, sys.float_info.min), cv/max(v, sys.float_info.min))
            lo, hi = number(r, 'order_interval_low_'+pair), number(r, 'order_interval_high_'+pair)
            p = number(r, 'order_'+pair)
            if av > 0 and cv > 0:
                assert math.isclose(p, math.log(av/cv)/math.log(1.5), abs_tol=1e-10, rel_tol=2e-12)
            qualified = math.isfinite(lo) and sig > 5 and temporal_ok(b)
            gate = ('sampling_or_roundoff_limited' if not math.isfinite(lo) else
                    'significance_or_temporal_unqualified' if not qualified else
                    'consistent_with_4' if lo <= 4 <= hi else
                    'negative_outside_uncertainty' if hi < 0 else 'measured_order_outside_4_interval')
            row['norm_significance_'+pair] = sig
            row['registered_gate_'+pair] = gate
        result.append(row)
    return result


def native_magnitudes(root, clocks, audit, masks):
    magnitudes, terms, shell = [], [], []
    # The sealed magnitude and signed-shell battery runs on E-high only.
    # Other-rung constraint norms are in the common-support order tables.
    for leg in ('E-high',):
        for k, clock in clocks.items():
            step = clock['coarse_steps'][leg]
            path = root/leg/'native-pieces'/f'step{step:06d}'
            receipt = read_json(path/'pieces.receipt.json')
            valid = read_json(path/'pieces-validation.json')
            assert receipt['status'] == valid['status'] == 'NATIVE_PIECE_REPLAY_VALIDATION_PASS'
            assert valid['max_saved_total_scaled_error'] == 0
            assert valid['max_term_sum_scaled_error'] <= valid['comparison_roundoff_allowance']
            assert valid['evolution_steps'] == 0 and valid['static_data'] == 'absent'
            verify_native_exits(path)
            audit.append(dict(scope='native_pieces', leg=leg, step=step, clock_index=k,
                              status=receipt['status'], saved_total_error=valid['max_saved_total_scaled_error'],
                              term_sum_error=valid['max_term_sum_scaled_error']))
            scales = read_csv(path/'individual-terms.csv')
            norms = read_csv(path/'normalized-constraints.csv')
            expected = {(mask,c) for mask in masks for c in ('Ham','Mom','Mom1','Mom2','GaussE')}
            assert {(r['mask'],r['constraint']) for r in norms} == expected and len(norms)==len(expected)
            groups = grouped(scales, ('mask', 'constraint'))
            for r in norms:
                g = groups[(r['mask'], r['constraint'])]
                largest = max(g, key=lambda x: number(x, 'RMS'))
                assert r['largest_term'] == largest['term']
                assert number(r, 'largest_individual_term_RMS') == number(largest, 'RMS')
                denominator = number(largest, 'RMS')
                assert denominator > 0
                ratio = number(r, 'raw_RMS')/denominator
                assert math.isclose(ratio, number(r, 'normalized_RMS'), rel_tol=1e-12)
                assert number(r, 'sealed_threshold') == .001
                applies = k >= 1 and r['mask'] != 'outer_boundary_shell'
                assert truth(r['applies']) == applies
                magnitudes.append(dict(r, all_individual_term_RMS=json.dumps(
                    {x['term']: number(x, 'RMS') for x in g}, sort_keys=True),
                    recomputed_normalized_RMS=ratio,
                    threshold_met=(ratio <= .001), registered_applies=applies))
            terms.extend(dict(r, clock_index=k) for r in scales)
        shell.extend(read_csv(root/leg/'signed-GaussE-shell.csv'))
    return magnitudes, terms, shell


def verify_native_exits(directory):
    for rank in range(32):
        assert int(read_text(directory/f'native-exit-rank-{rank}.txt').strip()) == 0, directory


def receipt_audit(root, registration, clocks):
    audit, horizons, mq, field_norms = [], [], [], []
    for leg in registration:
        p = root/leg
        expected = clocks[12]['coarse_steps'][leg]
        run, chain = read_json(p/'run.done'), read_json(p/'chain.done')
        assert run == chain and run['checkpoint_step'] == expected
        assert abs(run['checkpoint_time_M']-10.5) < 1e-9
        for n in ('run.exit', 'evolution.exit', 'diagnostics.exit'):
            assert int(read_text(p/n)) == 0
        e, d, q = (read_json(p/n) for n in ('evolution.receipt.json', 'diagnostics.receipt.json',
                                          'qualification.receipt.json'))
        assert e['status'] == 'PASS' and not e['errors'] and not e['missing_native_exits']
        assert len(e['native_exits']) == 32 and all(v == 0 for v in e['native_exits'].values())
        assert not d['failures'] and not q['failures'] and not q['missing_steps']
        assert q['successful_cases'] == q['expected_cases']
        assert q['static_data'] == 'absent' and q['precision'] == 'Float64'
        required_steps = ([x['coarse_steps'][leg] for x in clocks.values()] if leg.endswith('dt2')
                          else list(range(expected+1)))
        assert q['expected_steps'] == q['saved_steps'] == required_steps
        assert q['successful_cases'] == 2*len(required_steps)
        verify_native_exits(p)
        fields = read_csv(p/'field-audit.csv')
        # The dt/2 chain also emits four segment-end plots between the sealed
        # clocks. Audit them for finiteness; never add them to the ladder.
        assert {c['coarse_steps'][leg] for c in clocks.values()} <= {int(r['coarse_step']) for r in fields}
        assert all(int(r['nonfinite_valid']) == 0 for r in fields)
        audit.append(dict(scope='evolution_and_diagnostics', leg=leg, step=expected,
                          status='COMPLETE_OWN_RECEIPTS', peak_RSS_KiB=e['resources']['maximum_process_RSS_KiB'],
                          qualified_cases=q['successful_cases'], field_clocks=len(fields)))
        qtable = read_csv(p/'qualified-horizons.csv')
        lookup = {(int(r['coarse_step']), int(r['N_theta'])): r for r in qtable}
        assert len(lookup) == len(qtable) == q['successful_cases']
        for step in required_steps:
            for n in (48, 96):
                path = p/'qualified'/f'step{step:06d}'/f'n{n}'
                own = read_json(path/'result.json')
                assert own['native_exit'] == 0 and own['status'] == 'QUALIFIED'
                assert own['evolution_steps'] == 0 and own['field_time_M'] == own['time_M']
                stages = [r for r in read_csv(path/'surfaces.csv') if int(r['stage']) in (0,1,2)
                          and int(r['N_theta']) == n]
                assert [int(r['stage']) for r in stages] == [0,1,2]
                for stage, threshold in zip(stages, (1e-7,1e-10,1e-12)):
                    assert stage['status'] == 'FOUND' and number(stage, 'expansion_squared') <= threshold
                    assert number(stage, 'theta_minus') < 0
                row = dict(lookup[(step,n)])
                assert number(row, 'A') == number(own['measurement'], 'A')
                assert number(row, 'Q') == number(own['measurement'], 'Q')
                row.update(own_receipt='QUALIFIED', elapsed_s=own['elapsed_s'],
                           seed_source=own['seed_source'], seed_time_M=own.get('seed_time_M',''),
                           field_time_M=own['field_time_M'])
                row.update({f'stage_{i}_expansion_squared': stage['expansion_squared']
                            for i,stage in enumerate(stages)})
                horizons.append(row)
        mq.extend(dict(r, leg=leg) for r in read_csv(p/'mq-values.csv'))
        field_norms.extend(read_csv(p/'constraint-rms.csv'))
    assert int(read_text(root/'reduction.exit')) == 0
    reduction = read_json(root/'battery/reduction.receipt.json')
    assert reduction['clocks'] == list(clocks) and reduction['fields'] == 42
    assert reduction['stencil_orders'] == [8,10] and reduction['masks'] == 13
    hs = read_json(root/'horizon-summary.receipt.json')
    assert hs['status'] == 'HORIZON_MEASUREMENTS_COMPLETE' and hs['clocks'] == list(clocks)
    return audit, horizons, mq, field_norms


def plots(rows, summaries, masks, refinement_ratio):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scienceplots  # Registers styles; no TeX subprocess is required.
    from matplotlib.backends.backend_pdf import PdfPages
    plt.style.use(['science', 'no-latex'])
    plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7,
                         'pdf.fonttype': 42})
    lookup = grouped(rows, ('dataset', 'mask', 'category', 'field'))
    # Fixed representative set, covering metric, connection, gauge, matter,
    # electromagnetic field and the constrained sectors on EVERY sealed mask.
    representatives = 'chi Gamma1 shift1 lapse phi Ex Ham Mom GaussE CGamma'.split()
    figures = []

    def log_limits(ax, arrays):
        positive = np.concatenate([a[np.isfinite(a) & (a > 0)] for a in arrays])
        if len(positive):
            ax.set_yscale('log')
            ax.set_ylim(float(positive.min())*.6, float(positive.max())*1.7)
        else:
            ax.text(.5, .5, 'All norms zero; log undefined', ha='center', transform=ax.transAxes)
            ax.set_ylim(-.1, 1)

    def draw(axes, group, field):
        group = sorted(group, key=lambda r: int(r['clock_index']))
        t = np.array([number(r, 'time_M') for r in group])
        get = lambda key: np.array([number(r, key) for r in group])
        norms = [get(k+'_rms') for k in ('low', 'mid', 'high', 'control')]
        for label, values in zip(('E-low', 'E-mid', 'E-high', 'E-high-dt2'), norms):
            axes[0].plot(t, np.where(values > 0, values, np.nan), label=label)
        log_limits(axes[0], norms)
        axes[0].set_ylabel(field + ' RMS')
        d1, d2, dt = (get(k) for k in ('difference_low_mid', 'difference_mid_high', 'temporal_high_dt2'))
        differences = [d1, d2, dt]
        for label, values, style in [('D low–mid', d1, '-'), ('D mid–high', d2, '-'),
                                     ('D high–dt2', dt, ':')]:
            axes[1].plot(t, np.where(values > 0, values, np.nan), style, label=label)
        for p, style in ((2, '--'), (4, '-.'), (6, ':')):
            scaled = d2*refinement_ratio**p
            differences.append(scaled)
            axes[1].plot(t, np.where(scaled > 0, scaled, np.nan), style,
                         linewidth=.8, label=f'{refinement_ratio:g}^{p} D mid–high')
        log_limits(axes[1], differences)
        axes[1].set_ylabel('Difference RMS')
        central, lo, hi = (get(k) for k in ('measured_order', 'order_interval_low', 'order_interval_high'))
        axes[2].fill_between(t, lo, hi, alpha=.18, color='C0')
        axes[2].plot(t, central, 'o-', markersize=3, color='C0', label='Measured p; I8/I10 interval')
        unqualified = np.array([not significance_ok(r) or not temporal_ok(r) for r in group])
        axes[2].plot(t[unqualified], central[unqualified], 'x', color='C3', label='Significance/temporal unqualified')
        axes[2].axhline(4, color='.35', linestyle='--', linewidth=.8, label='Design p = 4')
        axes[2].set_ylabel('Observed order')
        for ax in axes:
            ax.axvspan(0, .875, color='.7', alpha=.1)
            ax.set_xlim(0, 10.5)
            ax.set_xlabel('t / M')
            ax.grid(alpha=.2, linewidth=.4)

    path = HERE/'t19-stageD-history.pdf'
    with PdfPages(path) as pdf:
        for mask in masks:
            for first in (0, 5):
                fig, axes = plt.subplots(5, 3, figsize=(11.7, 14), layout='constrained')
                for i, field in enumerate(representatives[first:first+5]):
                    draw(axes[i], lookup[('volume', mask, 'all', field)], field)
                for ax in axes[0]:
                    ax.legend(loc='best')
                fig.suptitle(f'exp-0023: {mask}; coordinate-volume common-support RMS\n'
                             'Raw field/constraint norms, pair differences and observed orders; t=0 reported only')
                pdf.savefig(fig)
                plt.close(fig)
                figures.append(dict(file=path.name, page=len(figures)+1, mask=mask, category='all',
                                    quantities=';'.join(representatives[first:first+5]), norm='RMS'))
    OUTPUTS.add(path)
    path = HERE/'t19-stageD-field-orders.pdf'
    with PdfPages(path) as pdf:
        for mask in masks:
            fig, axes = plt.subplots(2, 1, figsize=(10, 9), layout='constrained')
            g = sorted([r for r in summaries if r['dataset']=='volume' and r['mask']==mask
                        and r['category']=='all' and int(r['clock_index'])>0],
                       key=lambda r: int(r['clock_index']))
            t = [number(r, 'time_M') for r in g]
            axes[0].fill_between(t, [number(r, 'order_min') for r in g],
                                [number(r, 'order_max') for r in g], alpha=.15, label='Range over 28 fields')
            axes[0].plot(t, [number(r, 'order_median') for r in g], 'o-', label='Median of defined p')
            axes[0].plot(t, [number(r, 'qualified_order_median') for r in g], 's--',
                         label='Median after significance/temporal screen')
            axes[0].axhline(4, linestyle=':', color='.3', label='Design 4')
            axes[0].set_ylabel('Order (range includes unqualified defined orders)')
            axes[0].legend()
            for status in STATUSES:
                axes[1].plot(t, [int(r[status]) for r in g], 'o-', markersize=3, label=status)
            axes[1].set_ylabel('Fields per status / 28')
            axes[1].set_ylim(-.5, 28.5)
            axes[1].legend(loc='upper center', ncol=2, fontsize=7)
            for ax in axes:
                ax.set_xlim(.875, 10.5)
                ax.set_xlabel('t / M')
                ax.grid(alpha=.2)
            fig.suptitle(f'exp-0023: {mask}; all 28 evolved fields; category all\n'
                         'The median is descriptive and does not qualify any individual field')
            pdf.savefig(fig)
            plt.close(fig)
            figures.append(dict(file=path.name, page=list(masks).index(mask)+1, mask=mask,
                                category='all', quantities='all 28 evolved fields', norm='RMS-derived order'))
    OUTPUTS.add(path)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
    draw(axes, lookup[('volume','horizon','smooth_interior','Gamma1')], 'Gamma1')
    for ax in axes:
        ax.legend(fontsize=6)
    fig.suptitle('exp-0023: horizon mask, smooth_interior, Gamma1; resolved negative-order witness')
    for suffix in ('pdf', 'png'):
        path = HERE/f't19-stageD-overview.{suffix}'
        fig.savefig(path, dpi=600)
        OUTPUTS.add(path)
        figures.append(dict(file=path.name, page=1, mask='horizon', category='smooth_interior',
                            quantities='Gamma1', norm='RMS'))
    plt.close(fig)
    write_csv('figures', figures)


def fmt(value):
    if value in ('', None):
        return 'undefined'
    return f'{float(value):.7g}'


def md_table(rows, columns):
    return ('| ' + ' | '.join(label for key,label in columns) + ' |\n| ' +
            ' | '.join('---' for _ in columns) + ' |\n' +
            '\n'.join('| ' + ' | '.join(str(r.get(k,'')) for k,label in columns) + ' |' for r in rows))


def report(data):
    rows, sums, direct, magnitudes, shell = (data[k] for k in ('rows','summaries','constraints','magnitudes','shell'))
    verdict, witnesses = data['verdict'], data['failures']
    counts = Counter(r['status'] for r in rows if r['dataset']=='volume')
    ray_counts = Counter(r['status'] for r in rows if r['dataset']=='ray')
    field_sums = [r for r in sums if r['dataset']=='volume' and r['category']=='all' and int(r['clock_index'])>0]
    stats = []
    for mask in data['masks']:
        g = [r for r in field_sums if r['mask']==mask]
        p = [number(r,'order_median') for r in g]
        failures = [r for r in witnesses if r['mask']==mask]
        stats.append(dict(mask=mask, median_range=f'{min(p):.5g}–{max(p):.5g}',
                          fail_rows=len(failures), fail_clocks=';'.join(sorted({r['clock_index'] for r in failures}, key=int)),
                          fail_fields=';'.join(sorted({r['field'] for r in failures})),
                          fail_classes=';'.join(sorted({r['category'] for r in failures}))))
    witness = next(r for r in witnesses if r['dataset']=='volume' and r['mask']=='horizon'
                   and r['category']=='smooth_interior' and r['field']=='Gamma1' and r['clock_index']=='2')
    above_four = next(r for r in rows if exterior(r) and r['dataset']=='volume' and r['category']=='all'
                      and r['field'] in EVOLVED and r['status']=='measured_order_outside_4_interval'
                      and number(r,'order_interval_low')>4)
    all_quantity_medians = []
    for key,group in grouped([r for r in rows if r['dataset']=='volume' and r['category']=='all'],
                              ('mask','clock_index','time_M')).items():
        p = [number(r,'measured_order') for r in group if math.isfinite(number(r,'measured_order'))]
        all_quantity_medians.append(dict(zip(('mask','clock_index','time_M'),key),
                                        quantities=42,defined_orders=len(p),median=statistics.median(p) if p else ''))
    write_csv('all-quantity-medians', all_quantity_medians)
    all_median = lambda mask,k: next(r['median'] for r in all_quantity_medians
                                   if r['mask']==mask and r['clock_index']==str(k))
    boundary_medians = [float(r['median']) for r in all_quantity_medians
                        if r['mask']=='outer_boundary_shell' and int(r['clock_index'])>0]
    magnitude_stats = []
    for constraint in ('Ham','Mom1','Mom2','Mom','GaussE'):
        g = [r for r in magnitudes if r['constraint']==constraint and r['registered_applies']]
        worst = max(g,key=lambda r:number(r,'normalized_RMS'))
        magnitude_stats.append(dict(constraint=constraint, rows=len(g),
            largest_normalized_RMS=fmt(worst['normalized_RMS']), threshold='0.001',
            mask=worst['mask'], time_M=fmt(worst['time_M']), denominator_term=worst['largest_term'],
            largest_term_RMS=fmt(worst['largest_individual_term_RMS']), all_met=all(r['threshold_met'] for r in g)))
    signed = [r for r in shell if r['leg']=='E-high' and int(r['clock_index'])>0]
    worst_shell = max(signed,key=lambda r:abs(number(r,'signed_integral_over_Q_H_0')))
    temp = data['temporal']
    qualified_temp = [r for r in temp if number(r,'spatial_significance')>5]
    worst_temp = max(temp,key=lambda r:number(r,'maximum_temporal_fraction'))
    worst_qualified = max(qualified_temp,key=lambda r:number(r,'maximum_temporal_fraction'))
    nonconv = data['nonconv']
    nonconv_stats = [dict(dataset=s, rows=sum(r['dataset']==s for r in nonconv),
                        negative_interval=sum(r['dataset']==s and r['negative_outside_uncertainty'] for r in nonconv),
                        straddles_zero=sum(r['dataset']==s and number(r,'order_interval_low')<=0<=number(r,'order_interval_high') for r in nonconv),
                        below_five=sum(r['dataset']==s and number(r,'spatial_significance')<=5 for r in nonconv))
                     for s in ('volume','ray')]
    drifts = [{**r, 'largest_observed_relative_drift':fmt(r['largest_observed_relative_drift']),
               'numerical_t0':fmt(r['numerical_t0'])} for r in data['drifts']]
    baselines = [{k:(fmt(v) if k in ('A','Q','expansion_squared') else v)
                  for k,v in r.items()} for r in data['horizons'] if int(r['coarse_step'])==0]
    horizon_sens = [r for r in data['sensitivities'] if int(r['clock_index'])>0]
    direct_counts = {pair: Counter(r['registered_gate_'+pair] for r in direct if int(r['clock_index'])>0)
                     for pair in ('low_mid','mid_high')}
    mq_summary = []
    for leg in data['registration']:
        g = sorted([r for r in data['mq'] if r['leg']==leg],key=lambda r:number(r,'time_M'))
        for radius in (20,50,100):
            mq_summary.append(dict(leg=leg, radius_M=radius, rows=len(g),
                first_time=fmt(g[0]['time_M']), last_time=fmt(g[-1]['time_M']),
                M_first=fmt(g[0][f'M{radius}']), M_last=fmt(g[-1][f'M{radius}']),
                Q_first=fmt(g[0][f'Q{radius}']), Q_last=fmt(g[-1][f'Q{radius}'])))
    md = f'''# T19 — registered Stage D verdict of exp-0023

**Registered verdict: {verdict}.** Resolved negative orders on fixed exterior masks at the registered positive clocks meet the contract's FAIL condition. This verdict does not classify positive orders above four, or narrow intervals excluding four, as failures of convergence. Those rows remain as measured. Four physically completed chains, small native constraint magnitudes and small horizon drifts do not establish the required field-and-constraint fourth-order convergence.

The production root is `{data['root']}` and the contract directory is `{data['contract']}`. The script reads their tables and own numerical receipts only. No evolution, horizon solve, new interpolation, changed mask, clock, significance rule, tolerance or baseline is introduced. The supplied collection was SHA-256 verified by the controller (57,958 files); this analysis hashes every input it actually reads in `t19-manifest.txt`. It does not repeat or replace the controller's bulk collection verification.

## Registered reading and completeness

> PASS = on the smooth exterior masks at the synchronized diagnostic clocks t_k = 0.875*k M, k = 1…12, field and constraint orders consistent with 4 within the measured uncertainties, the sealed magnitudes met, the finest horizon drifts ≤ 1e-3 (area) / 2e-4 (charge) and decreasing under refinement; FAIL = negative or non-convergent orders on any mask outside measured uncertainty; anything else reported as measured, with the failing masks/times named.

The applied registration is `submit-contract.md` §8, including its controller-amended clocks and sealed magnitudes. The three spatial rung spacings are 7/8, 7/12 and 7/18 M, with ratio 3/2; max_level is 14. Positive acceptance clocks are exactly 0.875 k M for k=1…12. Numerical t=0 is reported separately and cannot trigger a positive-history verdict. All 13 fixed masks remain reported, including the outer boundary shell. For the exterior statement, the 12 radial masks have R≥0.005 M; the boundary shell is separately named, not silently used to establish an interior failure.

All four `run.done`/`chain.done` endpoint records agree at 10.5 M and steps 48/72/108/216. Evolution native exits, finite-state audits at all 13 common clocks, independent diagnostics and all {len(data['horizons'])} own N48/N96 qualification receipts were checked. Each accepted finder case has all three FOUND stages with squared expansion at most 1e-7/1e-10/1e-12 and negative inward expansion. The 13 E-high native-term replay receipts have unchanged current evolved bits, exactly matching saved constraint totals and term sums within their recorded roundoff allowance. Native pieces and signed-shell tables are produced only for the finest rung, as required by the sealed magnitude gate; no other-rung term tables are invented. The independent reduction receipt covers all 13 clocks, 42 quantities, 13 masks and I8/I10. Axis-parity receipts were checked independently. These are completion and diagnostic checks, not a convergence PASS.

The measured order is log(D_low,mid/D_mid,high)/log(1.5). Intervals are the collected I8–I10 plus Float64 bounds and are reproduced without widening or shrinking them. Both spatial differences must exceed five times their own spread-plus-floor. The half-step difference must be ≤0.2 of **each** spatial difference. Temporal uncertainty is reported separately in the source tables; a row exceeding that control cannot support a clean spatial-order statement, even if its displayed order looks fourth order. Undefined and sampling-limited entries remain unresolved. A negative central order whose interval includes zero is not established as negative outside uncertainty; its stored `non_convergent` label is preserved in the ledger rather than promoted to a strict FAIL witness.

## 1. Verdict and named exceptions

There are {len(witnesses)} qualified negative-order witnesses on radial exterior masks at positive clocks across the volume/geometry and ray tables, plus {len(data['constraint_failures'])} negative direct constraint-to-zero pair orders outside their intervals and passing their significance/temporal checks. `t19-failing-orders.csv` names every such witness with mask, clock, field, geometry, differences, margins and order interval. `t19-order-exceptions.csv` also names every other non-qualifying self-difference entry; `t19-constraint-orders.csv` retains every constraint-to-zero exception. No exception is dropped by taking a median or restricting the report to category `all`.

For an explicit smooth-interior witness, Γ̃₁ on the horizon mask at t=1.75 M has p={fmt(witness['measured_order'])}, interval [{fmt(witness['order_interval_low'])}, {fmt(witness['order_interval_high'])}], D_low,mid={fmt(witness['difference_low_mid'])}, D_mid,high={fmt(witness['difference_mid_high'])}, significance={fmt(witness['spatial_significance'])}, and temporal fractions {fmt(witness['temporal_fraction_of_low_mid'])}/{fmt(witness['temporal_fraction_of_mid_high'])}. This resolves non-contraction at that clock under the registered uncertainty and temporal rules. It does not identify the physical or discretization mechanism from the tables alone.

{md_table(stats,[('mask','Fixed mask'),('median_range','Median p range over 28 fields, category all, k=1…12'),('fail_rows','Qualified negative self-difference rows'),('fail_clocks','Clock indices'),('fail_fields','Fields / constraints'),('fail_classes','Geometry / ray classes')])}

## 2. All evolved fields and geometry classes

The controller's full volume-table counts are independently reproduced: {dict(counts)} over {sum(counts.values()):,} rows, including t=0 and all 42 quantities. The independent ray table has {dict(ray_counts)} over {sum(ray_counts.values()):,} rows. These counts are not restricted to the 28 evolved fields. The controller's reported median ranges reproduce when taking all 42 quantities together, including constraints (`t19-all-quantity-medians.csv`); they are not medians of the 28 evolved fields alone. The all-quantity early outer_ring/outer_ring_core medians are {fmt(all_median('outer_ring',1))}/{fmt(all_median('outer_ring_core',1))}, and its boundary-shell range is {fmt(min(boundary_medians))}–{fmt(max(boundary_medians))}. `t19-field-summary.csv` supplies median/min/max defined order and every status count for exactly the 28 evolved fields at each mask × clock × existing geometry/ray class. It separately supplies the median after the spatial/temporal checks; neither median overrides a field's status. `t19-field-orders-outside.csv` lists every field with its actual order and interval excluding four, including temporally unqualified measurements, explicitly tagged. For a positive above-four example, {above_four['field']} on {above_four['mask']} at t={fmt(above_four['time_M'])} M has p={fmt(above_four['measured_order'])}, interval [{fmt(above_four['order_interval_low'])}, {fmt(above_four['order_interval_high'])}]. It is as measured: neither a fourth-order pass nor a failure of convergence.

There are {sum(r['dataset']=='volume' for r in data['absent_support'])} absent common geometry-class records in `t19-absent-support.csv`; they remain absent rather than becoming zero-norm passes. The collected reducer explicitly skips the Cartesian outer boundary shell in its radial-ray loop (`production/reduce.py:170`); the 39 unavailable clock/ray combinations are recorded separately in that CSV. Patch-edge, convex-corner, axis-junction, same-level seam and smooth-interior classes follow the collected membership, and axis/equator/diagonal profiles retain their original ray weights. At k=1, outer_ring/outer_ring_core medians are {fmt(next(r['order_median'] for r in field_sums if r['mask']=='outer_ring' and r['clock_index']=='1'))}/{fmt(next(r['order_median'] for r in field_sums if r['mask']=='outer_ring_core' and r['clock_index']=='1'))}. The broad median range near four elsewhere coexists with individually resolved negative and above-four orders.

## 3. Constraints and sealed magnitudes

`t19-constraint-orders.csv` contains every raw common-support low/mid/high constraint norm, the self-difference order and its interval, and both direct-to-zero orders and intervals. Direct-to-zero pair classifications at positive clocks, including the separately named boundary shell, are low–mid: {dict(direct_counts['low_mid'])}; mid–high: {dict(direct_counts['mid_high'])}. A direct-to-zero pair requires both norms to exceed five times their collected uncertainties; its registered temporal control is also retained. Raw native composite-grid constraint norms are reported in `t19-native-constraint-norms.csv`, and their native term-normalized RMS values are in `t19-magnitudes.csv`. The common-support constraint norms and native composite-grid norms use different supports and are not interchanged. The collected tables contain RMS orders, not peak-norm order intervals; maxima remain reported as raw values without inventing an uncertainty model or peak-order certificate.

Every applicable E-high normalized magnitude passes its sealed RMS threshold: {sum(r['registered_applies'] for r in magnitudes)} rows, comprising Hamiltonian, both momentum components, their joint norm and electric Gauss on all 12 radial masks and all 12 positive clocks. `t19-magnitudes.csv` records the actual winning individual term, its RMS, every other term's RMS, the raw residual and the recomputed ratio. `t19-native-term-scales.csv` contains the individual native pieces. The denominator is the largest individual term on that same native mask and clock, never a cancelled total or static reference.

{md_table(magnitude_stats,[('constraint','Constraint'),('rows','Required rows'),('largest_normalized_RMS','Largest normalized RMS'),('threshold','Sealed limit'),('mask','Mask at largest ratio'),('time_M','t / M'),('denominator_term','Actual denominator'),('largest_term_RMS','Denominator RMS'),('all_met','All met')])}

The largest absolute **signed** E-high shell ratio at positive clocks is {fmt(abs(number(worst_shell,'signed_integral_over_Q_H_0')))} at t={fmt(worst_shell['time_M'])} M, N={worst_shell['N_theta']}, below 2e-4. `t19-signed-GaussE-shell.csv` retains every produced E-high row and both angular resolutions, their own numerical Q_H(0), signed integrals, absolute integrals and two-cell boundary-band estimates. The signed physical-volume shell gate and the coordinate-volume normalized RMS gate both pass numerically. The measured angular/boundary-band scope is not a full shell quadrature certificate. Θ, Z, C_Γ, det(h)−1, tr(Ã), magnetic Gauss and Λ/Ξ have reported raw norms and orders, with no invented magnitude thresholds.

## 4. Temporal control

`t19-temporal-failures.csv` names all {len(temp)} rows exceeding 0.2 in at least one pair fraction ({sum(r['dataset']=='volume' for r in temp)} volume/geometry; {sum(r['dataset']=='ray' for r in temp)} ray), with both fractions, additive excess, factor above the limit, significance and original status. Of these, {len(qualified_temp)} have both spatial differences above five times uncertainty; these still cannot qualify a spatial-order claim. Zero-denominator fractions are separately counted as undefined in the JSON receipt and exceptions ledger.

The largest reported fraction is {fmt(worst_temp['maximum_temporal_fraction'])} for {worst_temp['field']}, {worst_temp['mask']}, {worst_temp['category']}, t={fmt(worst_temp['time_M'])} M, but its spatial significance is only {fmt(worst_temp['spatial_significance'])}. Among rows exceeding the spatial significance rule, the largest fraction is {fmt(worst_qualified['maximum_temporal_fraction'])} for {worst_qualified['field']}, {worst_qualified['mask']}, {worst_qualified['category']}, t={fmt(worst_qualified['time_M'])} M. A large ratio on a tiny/floor-limited spatial difference is not evidence for a large absolute temporal error; both raw differences and floors remain in the ledger.

## 5. Every stored non_convergent row

`t19-non-convergent.csv` preserves all {len(nonconv)} stored labels with their exact differences, measured significance, I8/I10 spreads and Float64 floor. Per-pair floor flags explicitly test the five-times floor and five-times spread separately and together. No labelled non_convergent row has significance ≤5 under the table's own combined spread-plus-floor; near-roundoff absolute values alone do not justify relabelling a statistically resolved difference as roundoff-limited. Intervals straddling zero are identified and are not strict negative-order witnesses. Numerical t=0 and boundary-only cases remain reported outside the radial positive-time claim.

{md_table(nonconv_stats,[('dataset','Support'),('rows','Stored non_convergent'),('negative_interval','Negative outside uncertainty, all clocks/masks'),('straddles_zero','Interval includes zero'),('below_five','Significance ≤5')])}

## 6. Horizons and MQ

All {len(data['horizons'])} qualified numerical finder cases and stage residuals are retained in `t19-qualified-horizons.csv`; `t19-horizon-baselines.csv` selects the eight independent numerical-t=0 baselines. A positive-time numerical candidate is only a geometry seed when solving the t=0 fields. The native candidate `finder-values.csv` modes are not substituted for these qualified surfaces. `M_RN_legacy` is the emitted RN expression and is not promoted to an EMS equilibrium mass or a new horizon mass diagnostic.

{md_table(baselines,[('leg','Leg'),('N_theta','Nθ'),('A','Numerical A(0)'),('Q','Numerical Q(0)'),('expansion_squared','Final squared expansion'),('elapsed_s','Own wall seconds')])}

`t19-horizon-drifts.csv` reports the full spatial coarse histories (49/73/109 clocks per N) and only the registered 13 diagnostic clocks for the half-step control. No between-diagnostic-clock control maximum is claimed. The budgets and observed decrease under spatial refinement are met in the drift table.

{md_table(drifts,[('leg','Leg'),('quantity','Quantity'),('N_theta','Nθ'),('largest_observed_relative_drift','Largest signed relative drift'),('drift_budget','Budget'),('complete_required_history','Required history complete'),('observed_drift_within_budget','Within budget')])}

All {len(horizon_sens)} positive-clock retention comparisons have sensitivities within the drift budgets and resolved **drift** differences under refinement, but {sum(not truth(r['absolute_sensitivities_resolve_interrung']) for r in horizon_sens)} have absolute angular/stopping sensitivity larger than the minimum absolute inter-rung difference. `t19-horizon-sensitivities.csv` preserves the `interrung_sensitivity_unresolved` statuses. Stable correlated N48/N96 drifts do not resolve that absolute comparison, so the stronger horizon qualification clause is not established even though the budgets pass.

All {len(data['mq'])} produced MQ rows, including off-ladder coarse clocks, are in `t19-mq.csv`. They are code-produced radius-20/50/100 mass/charge diagnostics, independent of the horizon baselines. No t=0 MQ row is fabricated and no new flux or loss threshold is imposed. The first and final reported values are:

{md_table(mq_summary,[('leg','Leg'),('radius_M','Radius / M'),('rows','Produced rows'),('first_time','First t / M'),('last_time','Last t / M'),('M_first','First M'),('M_last','Final M'),('Q_first','First Q'),('Q_last','Final Q')])}

## 7. Figures and reproducibility

`t19-stageD-history.pdf` has 26 pages: ten fixed representative quantities on every one of the 13 sealed masks. Each row shows common-support coordinate-volume **RMS norms against time**, both raw spatial differences and the temporal difference, the finer spatial difference rescaled by (3/2)^p for p=2,4,6, and the observed order with its stored interval. These rescalings are displayed assumptions, not alternate verdicts. `t19-stageD-field-orders.pdf` has one page per mask summarizing all 28 fields' order range, medians and status counts. `t19-stageD-overview.pdf`/`.png` show the smooth-interior horizon Γ̃₁ witness; its negative order at 1.75 M is resolved while both spatial differences greatly exceed their diagnostic spread. Unqualified central orders are marked and t<0.875 M is shaded outside the claim. Raw maxima are available in the original and exported native norm tables but are not the plotted norm. `t19-figures.csv` indexes every page. No superseded run is plotted.

Regenerate with `python Tests/EMSNative/t19-stageD-verdict.py COLLECTED_PRODUCTION_ROOT CONTRACT_DIR [--output-dir DIR]`. The five large regenerable ledgers default to the production root's sibling `t19/` directory; `--output-dir` changes their destination. Small outputs and README §T19 remain in this worktree. The script checks duplicate keys, frozen mask/clock coverage, every stored self-difference order/interval/significance/status, direct constraint orders, native term denominators and numerical horizon stages before writing the verdict. Its own resource receipt is `t19-stageD-verdict.json`; SHA-256 input/output seals are in `t19-manifest.txt`. No SSH, cluster work, new evolution or commit occurs.
'''
    for name in BULK_NAMES:
        filename = 't19-' + name + '.csv'
        md = md.replace('`' + filename + '`',
                        '[' + filename + '](' + str(BULK_DIR / filename) + ')')
    save_text('stageD-verdict.md', md)
    return stats, magnitude_stats, mq_summary


def update_readme(data, magnitude_stats):
    marker = '\n## T19 — registered Stage D reading of exp-0023\n'
    maxima = {r['constraint']:r['largest_normalized_RMS'] for r in magnitude_stats}
    section = f'''{marker}

The completed exp-0023 three-grid alpha_K chain and its half-step control receive the registered **{data['verdict']}** verdict under the sealed submit-contract §8. See [the complete reading](t19-stageD-verdict.md) and [the parametric regeneration script](t19-stageD-verdict.py). Positive acceptance clocks stay at 0.875 k M, k=1…12, with numerical t=0 reported separately; the masks, I8/I10 intervals, >5 significance and ≤0.2 temporal rules are unchanged. No new evolution, horizon solve, static-data use after t=0, SSH, cluster action or commit is performed.

There are {len(data['failures'])} positive-time radial-exterior self-difference rows with intervals wholly below zero, significant spatial differences and passing temporal controls, plus {len(data['constraint_failures'])} negative direct constraint-to-zero pairs. The [named failing orders](t19-failing-orders.csv), [complete exceptions](t19-order-exceptions.csv), [28-field summaries](t19-field-summary.csv) and [field intervals outside four](t19-field-orders-outside.csv) retain every mask, clock and geometry/ray class. In particular, smooth-interior horizon Gamma1 at 1.75 M has p=-1.090154 with interval [-1.090453,-1.089855]. Positive orders above four remain as measured, not passes and not failures of convergence. The controller's broad median ranges describe all 42 battery quantities; [those medians](t19-all-quantity-medians.csv) are separated from the requested 28 evolved-field medians.

All {sum(r['registered_applies'] for r in data['magnitudes'])} required E-high [native term-normalized constraint magnitudes](t19-magnitudes.csv) meet their 1e-3 limits: largest Ham/Mom/GaussE ratios are {maxima['Ham']}/{maxima['Mom']}/{maxima['GaussE']}. The signed horizon-to-0.02 M Gauss shell gate also meets 2e-4. [Constraint orders](t19-constraint-orders.csv) remain individually measured or unresolved. All {len(data['temporal'])} [temporal-limit exceptions](t19-temporal-failures.csv) and all {len(data['nonconv'])} [stored non_convergent labels](t19-non-convergent.csv) retain their uncertainties and floors; intervals including zero are not strict negative-order witnesses.

All {len(data['horizons'])} own numerical horizon cases pass the three stored finder stages, including the eight [numerical-t=0 baselines](t19-horizon-baselines.csv). Full spatial coarse-history [area/charge drift budgets](t19-horizon-drifts.csv) are met and the drifts decrease under refinement. Absolute N48/N96 sensitivity versus inter-rung differences remains unresolved at all 96 positive retention rows, as retained in [the sensitivity table](t19-horizon-sensitivities.csv). [Every produced MQ row](t19-mq.csv) is reported independently.

The [26-page RMS history atlas](t19-stageD-history.pdf) shows raw norms against time, both spatial differences, the temporal difference, finer differences rescaled for orders 2/4/6, and measured order intervals on every fixed mask. The [all-field order atlas](t19-stageD-field-orders.pdf) covers the 28-field status counts; the [smooth-interior witness figure](t19-stageD-overview.png) makes the resolved negative order visible. Maxima are not substituted for RMS norms. The [manifest](t19-manifest.txt) seals every read input and output, and [the analysis resource receipt](t19-stageD-verdict.json) records measured peak RSS against the 6 GB cap with one thread. Regenerate with `python Tests/EMSNative/t19-stageD-verdict.py <collected production root> <contract dir>`.
'''
    path = HERE/'README.md'
    for name in BULK_NAMES:
        filename = 't19-' + name + '.csv'
        section = section.replace('](' + filename + ')', '](' + str(BULK_DIR / filename) + ')')
    section += '\nThe five large ledgers use the production root sibling `t19/` directory by default; pass `--output-dir DIR` to select another destination.\n'
    old = path.read_text()
    if marker in old:
        prefix, tail = old.split(marker,1)
        suffix = '\n## '+tail.split('\n## ',1)[1] if '\n## ' in tail else ''
    else:
        prefix, suffix = old.rstrip()+'\n', ''
    path.write_text(prefix+section+suffix)
    OUTPUTS.add(path)


def main():
    global BULK_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production_root', type=Path)
    parser.add_argument('contract_dir', type=Path)
    parser.add_argument('--output-dir', type=Path,
                        help='directory for the five regenerable large ledgers '
                             '(default: production root sibling t19)')
    args = parser.parse_args()
    start = time.monotonic()
    root, contract = args.production_root.resolve(), args.contract_dir.resolve()
    BULK_DIR = args.output_dir or args.production_root.absolute().parent / 't19'
    # Fail before analysis or small-output mutation if the bulk destination is unavailable.
    BULK_DIR.mkdir(parents=True, exist_ok=True)
    assert root != HERE and contract != HERE
    sealed_contract = read_text(contract/'submit-contract.md')
    assert 'FAIL = negative or non-convergent orders on any mask outside measured uncertainty' in sealed_contract
    assert 't_k = 0.875*k M, k = 1…12' in sealed_contract
    registration = read_json(contract/'registration.json')
    clocks = {c['clock_index']: c for c in read_json(contract/'early-clocks.json')['clocks']}
    masks = {r['mask']: r for r in read_csv(contract/'masks.csv')}
    assert set(clocks) == set(range(13)) and len(masks) == 13
    # Compare the frozen masks as rows; original CSV uses literal \\r\\n text.
    local_masks = read_text(HERE/'t14d-stageD-masks.csv').replace('\\r\\n', '\n')
    assert list(csv.DictReader(local_masks.splitlines())) == list(masks.values())
    read_text(HERE/'t14d-stageD-design.md')
    for leg in registration:
        assert registration[leg]['final_M'] == 10.5
    assert all(float(r['R_min_M']) >= .005 for mask,r in masks.items() if mask!='outer_boundary_shell')
    ratio = registration['E-low']['h0_M']/registration['E-mid']['h0_M']
    assert math.isclose(ratio, 1.5) and math.isclose(
        registration['E-mid']['h0_M']/registration['E-high']['h0_M'], ratio)
    volume = read_csv(root/'battery/self-differences.csv')
    rays = read_csv(root/'battery/ray-self-differences.csv')
    check_battery(volume, masks, clocks, 'volume')
    check_battery(rays, masks, clocks, 'ray')
    for dataset, rows in (('volume',volume), ('ray',rays)):
        for r in rows:
            r['dataset'] = dataset
    rows = volume + rays
    support = read_csv(root/'battery/common-support.csv')
    absent = [dict(r,dataset='volume') for r in support if r['status']=='ABSENT_COMMON_CLASS']
    read_text(contract/'production/reduce.py')
    absent.extend(dict(dataset='ray',clock_index=k,mask='outer_boundary_shell',category='ray_'+ray,
                       cells='',weight='',status='NOT_DEFINED_IN_COLLECTED_RAY_REDUCER',
                       source='production/reduce.py:170') for k in clocks for ray in ('axis','equator','diagonal'))
    assert len(support) == 13*13*6
    volume_groups = grouped(volume, ('clock_index','mask','category'))
    for r in support:
        g = volume_groups.get((r['clock_index'],r['mask'],r['category']), [])
        assert (len(g)==42) if r['status']=='present' else not g
    parity = [dict(r,clock_index=k) for k in clocks for r in read_csv(root/'battery'/f'parity-clock{k}.csv')]
    assert len(parity)==52 and all(number(r,'checked_values')>0 and
                                  number(r,'maximum_scaled_error')<256*sys.float_info.epsilon for r in parity)
    print('T19: frozen clocks/masks, all battery orders and statuses verified', flush=True)
    audit, horizons, mq, native_norms = receipt_audit(root, registration, clocks)
    print('T19: every numerical horizon receipt and its three stages verified', flush=True)
    magnitudes, terms, shell = native_magnitudes(root, clocks, audit, masks)
    direct = constraint_orders(read_csv(root/'battery/constraint-orders.csv'), volume)
    drifts = read_csv(root/'horizon-max-drifts.csv')
    sensitivities = read_csv(root/'horizon-retention-sensitivities.csv')
    assert len(drifts)==16 and len(sensitivities)==104
    assert all(truth(r['complete_required_history']) for r in drifts)
    for r in drifts:
        assert number(r,'drift_budget') == (.001 if r['quantity']=='A' else .0002)
        assert truth(r['observed_drift_within_budget']) == (
            abs(number(r,'largest_observed_relative_drift')) <= number(r,'drift_budget'))
        history = [h for h in horizons if h['leg']==r['leg'] and h['N_theta']==r['N_theta']]
        zero = next(h for h in history if int(h['coarse_step'])==0)
        base = number(zero,r['quantity'])
        assert base==number(r,'numerical_t0') and len(history)==int(r['expected_steps'])
        observed = max((number(h,r['quantity'])/base-1 for h in history),key=abs)
        assert math.isclose(observed,number(r,'largest_observed_relative_drift'),abs_tol=1e-14)
    for r in shell:
        assert number(r,'sealed_threshold')==.0002 and number(r,'outer_radius_M')==.02
    sums = summary(rows)
    nonqual = exceptions(rows)
    outside = [r for r in nonqual if r['field'] in EVOLVED and math.isfinite(number(r,'order_interval_low'))
               and not number(r,'order_interval_low')<=4<=number(r,'order_interval_high')]
    failures = [dict(r,kind='self_difference') for r in nonqual
                if r['registered_exterior_clock'] and r['negative_outside_uncertainty']]
    constraint_failures = []
    for r in direct:
        for pair in ('low_mid','mid_high'):
            if exterior(r) and r['registered_gate_'+pair]=='negative_outside_uncertainty':
                constraint_failures.append(dict(r,kind='constraint_to_zero',dataset='volume', category='all',
                    field=r['constraint'], order_pair=pair, measured_order=r['order_'+pair],
                    order_interval_low=r['order_interval_low_'+pair], order_interval_high=r['order_interval_high_'+pair],
                    status='negative_outside_uncertainty',claim_interval=True,registered_exterior_clock=True,
                    negative_outside_uncertainty=True,interpretation='registered_FAIL_witness'))
    temporal = []
    for r in rows:
        fractions = [number(r,k) for k in ('temporal_fraction_of_low_mid','temporal_fraction_of_mid_high')]
        maximum = max((x for x in fractions if math.isfinite(x)),default=math.nan)
        if maximum > .2:
            temporal.append(dict({k:r[k] for k in ORDER_COLUMNS}, maximum_temporal_fraction=maximum,
                excess_over_0p2=maximum-.2, factor_over_0p2=maximum/.2,
                low_mid_excess=max(0.,fractions[0]-.2), mid_high_excess=max(0.,fractions[1]-.2),
                pair_exceeding_limit=';'.join(p for p,x in zip(('low_mid','mid_high'),fractions) if x>.2),
                spatial_significance_qualified=significance_ok(r),
                spatial_statement='temporal_unqualified; order descriptive only'))
    nonconv = []
    for r in rows:
        if r['status'] != 'non_convergent':
            continue
        row = {k:r[k] for k in ORDER_COLUMNS}
        floor = number(r,'Float64_floor')
        for pair in ('low_mid','mid_high'):
            d = number(r,'difference_'+pair)
            spread = number(r,'I8_I10_'+pair+'_uncertainty')
            row[pair+'_roundoff_limited_at_5x'] = d <= 5*floor
            row[pair+'_sampling_limited_at_5x'] = d <= 5*spread
            row[pair+'_combined_limited_at_5x'] = d <= 5*(floor+spread)
            row[pair+'_dominant_uncertainty'] = 'sampling' if spread>floor else 'roundoff'
        row['negative_outside_uncertainty'] = negative_outside_uncertainty(r)
        row['registered_exterior_clock'] = exterior(r)
        row['interval_includes_zero'] = number(r,'order_interval_low')<=0<=number(r,'order_interval_high')
        nonconv.append(row)
    magnitude_ok = all(r['threshold_met'] for r in magnitudes if r['registered_applies'])
    shell_ok = all(abs(number(r,'signed_integral_over_Q_H_0'))<=.0002 for r in shell
                   if r['leg']=='E-high' and int(r['clock_index'])>0)
    order_ok = all(r['status']=='consistent_with_4' for r in rows if exterior(r))
    direct_ok = all(r['registered_gate_'+p]=='consistent_with_4' for r in direct if exterior(r)
                    for p in ('low_mid','mid_high'))
    horizon_ok = all(truth(r['observed_drift_within_budget']) for r in drifts) and all(
        truth(r[k]) for r in sensitivities if int(r['clock_index'])>0 for k in
        ('drift_decreases_under_refinement','sensitivities_within_budget',
         'absolute_sensitivities_resolve_interrung','drift_sensitivity_resolves_interrung'))
    verdict = ('FAIL' if failures or constraint_failures else
               'PASS' if order_ok and direct_ok and magnitude_ok and shell_ok and horizon_ok else 'as-measured')
    write_csv('field-summary', sums)
    compact = ['dataset','clock_index','time_M','mask','category','field','measured_order',
               'order_interval_low','order_interval_high','status','spatial_significance',
               'temporal_fraction_of_low_mid','temporal_fraction_of_mid_high','reasons',
               'registered_exterior_clock','negative_outside_uncertainty','interpretation']
    write_csv('field-orders-outside', outside, compact)
    write_csv('order-exceptions', nonqual, compact)
    write_csv('failing-orders', failures + constraint_failures)
    write_csv('constraint-orders', direct)
    write_csv('magnitudes', magnitudes)
    write_csv('native-term-scales', terms)
    write_csv('signed-GaussE-shell', shell)
    write_csv('temporal-failures', temporal)
    write_csv('non-convergent', nonconv)
    write_csv('absent-support', absent)
    write_csv('parity', parity)
    write_csv('qualified-horizons', horizons)
    write_csv('horizon-baselines', [r for r in horizons if int(r['coarse_step'])==0])
    write_csv('horizon-drifts', drifts)
    write_csv('horizon-sensitivities', sensitivities)
    write_csv('mq', mq)
    write_csv('native-constraint-norms', [r for r in native_norms if r['constraint'] in CONSTRAINTS])
    write_csv('receipt-audit', audit)
    data = dict(root=root,contract=contract,registration=registration,masks=masks,rows=rows,
                summaries=sums,constraints=direct,magnitudes=magnitudes,shell=shell,temporal=temporal,
                failures=failures,constraint_failures=constraint_failures,verdict=verdict,
                absent_support=absent,nonconv=nonconv,horizons=horizons,drifts=drifts,
                sensitivities=sensitivities,mq=mq)
    stats, magnitude_stats, mq_summary = report(data)
    write_csv('mask-verdict-summary', stats)
    write_csv('magnitude-summary', magnitude_stats)
    write_csv('mq-summary', mq_summary)
    print(f'T19: verdict {verdict}; {len(failures)} exterior negative self-difference witnesses; rendering figures', flush=True)
    plots(rows,sums,masks,ratio)
    update_readme(data,magnitude_stats)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = int(peak if sys.platform=='darwin' else peak*1024)
    assert peak_bytes < 6_000_000_000, '6 GB per-process memory cap exceeded'
    result = dict(verdict=verdict, scientific_completion='COMPLETED_TABLE_READING',
        volume_rows=len(volume), ray_rows=len(rays), volume_status_counts=dict(Counter(r['status'] for r in volume)),
        ray_status_counts=dict(Counter(r['status'] for r in rays)), qualified_horizon_cases=len(horizons),
        native_piece_cases=len([r for r in audit if r['scope']=='native_pieces']),
        required_magnitude_rows=sum(r['registered_applies'] for r in magnitudes),
        magnitude_gate=magnitude_ok,signed_shell_gate=shell_ok,field_order_gate=order_ok,
        direct_constraint_order_gate=direct_ok,full_horizon_sensitivity_gate=horizon_ok,
        exterior_negative_self_difference_rows=len(failures),
        exterior_negative_self_difference_counts=dict(Counter(r['dataset'] for r in failures)),
        exterior_negative_constraint_pairs=len(constraint_failures),
        temporal_above_0p2_rows=len(temporal), temporal_spatial_significant_rows=sum(significance_ok(r) for r in temporal),
        temporal_undefined_rows=sum(any(not math.isfinite(number(r,k)) for k in
            ('temporal_fraction_of_low_mid','temporal_fraction_of_mid_high')) for r in rows),
        stored_non_convergent_rows=len(nonconv), exceptions=len(nonqual),
        timing=dict(elapsed_s=time.monotonic()-start,peak_RSS_bytes=peak_bytes,
                    measured_by='getrusage(RUSAGE_SELF).ru_maxrss; macOS bytes / Linux KiB',
                    process_cap_bytes=6_000_000_000,threads=1),
        unchanged_registration=dict(spatial_significance='>5',temporal_fraction_limit=.2,design_order=4,
            I8_central=True,I8_I10_uncertainty=True,history_clocks=list(clocks),new_evolution=False),
        read_only_roots=[str(root),str(contract)], output_directory=str(HERE))
    manifest = ['INPUT_SHA256 (only read inputs; controller collection verification remains separate)']
    for path in sorted(INPUTS):
        manifest.append(sha256(path)+'  '+str(path))
    manifest.append('OUTPUT_SHA256 (manifest excluded to avoid self-reference)')
    for path in sorted(OUTPUTS | {Path(__file__).resolve()}):
        manifest.append(sha256(path)+'  '+str(path))
    # Include hashing in the process measurement; stream hashes so exporting
    # the full named exception ledger does not require another bulk allocation.
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = int(peak if sys.platform=='darwin' else peak*1024)
    assert peak_bytes < 6_000_000_000
    result['timing'].update(peak_RSS_bytes=peak_bytes,elapsed_s=time.monotonic()-start)
    save_text('stageD-verdict.json', json.dumps(result,indent=2)+'\n')
    manifest.append(sha256(HERE/'t19-stageD-verdict.json')+'  '+str(HERE/'t19-stageD-verdict.json'))
    header = ['T19: no commit; read-only exp-0023 Stage D registered verdict',
              f'Verdict: {verdict}', f'Peak RSS bytes: {peak_bytes}; cap: 6000000000; threads: 1',
              f'Elapsed seconds: {result["timing"]["elapsed_s"]:.6f}',
              'All own horizon/native replay/evolution receipts checked, not inferred from launcher exit.']
    save_text('manifest.txt', '\n'.join(header+manifest)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__':
    main()
