#!/usr/bin/env python3
"""T5 CSVs, measured ridge tracks and figures from t5-analyze's cache."""
import csv
import os
import time
import resource
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path.home()/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scienceplots
import numpy as np
from scipy.signal import find_peaks, peak_widths

HERE = Path(__file__).resolve().parent
CACHE = Path('/private/tmp/ems-t5/t5-maps.npz')
with np.load(CACHE) as archive:
    F = dict(archive)
T, R, L = F['time'], F['radii'], F['line_radii']
SCALES = ('low', 'mid', 'high')
RAYS = ('axis_plus', 'axis_minus', 'equator', 'diagonal')
FIELDS = ('chi', 'lapse', 'phi', 'Theta', 'Qscalar')
INDEX = (0, 1, 2, 3, 9)
Q4 = 5.0625
LINES, LINES8 = F['lines'], F['lines8']
DIFF, DIFF8 = F['diff'], F['diff8']
COMMON, NATIVE, INTERIOR = F['common'], F['native'], F['interior']


def order(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where((a > 0) & (b > 0), np.log(a/b)/np.log(1.5), np.nan)


def ratio(a, b):
    return np.divide(a, b, out=np.full(np.broadcast_shapes(np.shape(a), np.shape(b)), np.nan), where=b != 0)


def save(name, rows):
    with (HERE/name).open('w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def level(r, ray):
    extent = r/np.sqrt(2) if ray == 'diagonal' else r
    return max([0]+[lev for lev, bound in ((1, 64), (2, 32), (3, 16), (4, 8), (5, 4), (6, 2)) if extent <= bound])


def widths(y, k):
    _, _, left, right = peak_widths(y, [k], rel_height=.5)
    return float((right[0]-left[0])*(L[1]-L[0]))


def fit(name, ray, method, times, radii):
    times, radii = np.array(times), np.array(radii)
    good = np.isfinite(radii)
    times, radii = times[good], radii[good]
    if len(times) < 3:
        return None
    (v, intercept), cov = np.polyfit(times, radii, 1, cov=True)
    resid = radii-v*times-intercept
    return dict(branch=name, ray=ray, method=method, time_start_M=times[0], time_end_M=times[-1],
                samples=len(times), speed_dr_dt=v, fit_standard_error=np.sqrt(cov[0, 0]),
                intercept_M=intercept, radial_residual_rms_M=np.sqrt(np.mean(resid**2)),
                track_status='AMBIGUOUS' if name.startswith('inward') and np.any(np.diff(radii) > .015625) else 'TRACKED')


def csvs():
    constraints, fields, lines, outer = [], [], [], []
    directions = F['directions'][:7]
    for ti, t in enumerate(T):
        for d, direction in enumerate(directions):
            for ri, r in enumerate(R):
                if not F['counts'][ti, d, ri]:
                    continue
                for c, name in enumerate(('Ham', 'Mom', 'GaussE')):
                    row = dict(time_M=t, direction=direction, rlo_M=r-.0625, rhi_M=r+.0625,
                               common_centres=int(F['counts'][ti, d, ri]), constraint=name)
                    for typ, a in (('common', COMMON), ('native', NATIVE), ('interior', INTERIOR)):
                        v = a[ti, :, d, c, ri]
                        row.update({typ+'_'+s+'_rms': v[i] for i, s in enumerate(SCALES)})
                        row[typ+'_order_low_mid'], row[typ+'_order_mid_high'] = order(v[:-1], v[1:])
                    constraints.append(row)
                for c, name in enumerate(FIELDS):
                    a, b = DIFF[ti, :, d, c, ri]
                    a8, b8 = DIFF8[ti, :, d, c, ri]
                    fields.append(dict(time_M=t, direction=direction, rlo_M=r-.0625, rhi_M=r+.0625,
                                       common_centres=int(F['counts'][ti, d, ri]), field=name,
                                       d_low_mid_rms=a, d_mid_high_rms=b, scaled_mid_high_rms=Q4*b,
                                       self_order=order(a, b), self_order_P8=order(a8, b8),
                                       fourth_residual_rms=F['mismatch'][ti, d, c, ri],
                                       fourth_residual_relative=ratio(F['mismatch'][ti, d, c, ri], a),
                                       P6_P8_low_mid_relative=ratio(F['sensitivity'][ti, 0, d, c, ri], a),
                                       P6_P8_mid_high_relative=ratio(F['sensitivity'][ti, 1, d, c, ri], b)))
        for d, ray in enumerate(RAYS):
            for ri in range(4, len(L), 8):
                row = dict(time_M=t, ray=ray, radius_M=L[ri])
                for si, scale in enumerate(SCALES):
                    for c, name in enumerate(('chi', 'lapse', 'phi', 'Theta', 'Ham', 'Mom1', 'Mom2', 'GaussE', 'GaussB', 'Qscalar')):
                        row[scale+'_'+name] = LINES[ti, si, d, c, ri]
                lines.append(row)
        for d, ray in enumerate(RAYS[:3]):
            for ri, r in enumerate(F['outer_radii']):
                row = dict(time_M=t, ray=ray, radius_M=r)
                for si, scale in enumerate(SCALES):
                    for c, name in enumerate(('chi', 'lapse', 'phi', 'Theta', 'Ham', 'Mom1', 'Mom2', 'GaussE', 'GaussB', 'Qscalar')):
                        row[scale+'_'+name] = F['outer'][ti, si, d, c, ri]
                outer.append(row)
    save('t5-constraint-profiles.csv', constraints)
    save('t5-field-differences.csv', fields)
    save('t5-line-profiles.csv', lines)
    save('t5-outer-profiles.csv', outer)
    global_rows = []
    for ti, t in enumerate(T):
        for d, direction in enumerate(F['directions'][:7]):
            for ri, r in enumerate(F['global_radii']):
                if not np.isfinite(F['global_native'][ti, 0, d, 0, ri]):
                    continue
                row = dict(time_M=t, direction=direction, rlo_M=r-1, rhi_M=r+1)
                for c, name in enumerate(('Ham', 'Mom', 'GaussE')):
                    for si, scale in enumerate(SCALES):
                        row[scale+'_'+name+'_native_rms'] = F['global_native'][ti, si, d, c, ri]
                for c, name in enumerate(FIELDS):
                    row[name+'_low_mid_rms'], row[name+'_mid_high_rms'] = F['global_diff'][ti, :, d, c, ri]
                    row[name+'_scaled_mid_high_rms'] = Q4*row[name+'_mid_high_rms']
                    row[name+'_self_order'] = order(row[name+'_low_mid_rms'], row[name+'_mid_high_rms'])
                global_rows.append(row)
    save('t5-global-profiles.csv', global_rows)
    norms = list(csv.DictReader((HERE/'t5-mask-norms.csv').open()))
    lookup = {(float(r['time_M']), r['mask'], r['constraint'], r['scale']): r for r in norms}
    orders = []
    for t, mask, field in sorted({(float(r['time_M']), r['mask'], r['constraint']) for r in norms}):
        v = [lookup[t, mask, field, s] for s in SCALES]
        a = np.array([float(r['cylindrical_rms']) for r in v])
        orders.append(dict(time_M=t, mask=mask, constraint=field,
                           **{s+'_rms': a[i] for i, s in enumerate(SCALES)},
                           order_low_mid=order(a[0], a[1]), order_mid_high=order(a[1], a[2]),
                           **{s+'_cells': v[i]['cells'] for i, s in enumerate(SCALES)}))
    save('t5-mask-orders.csv', orders)


def tracks():
    rows, fits, crossings, pulses, packets = [], [], [], [], []
    for d, ray in enumerate(RAYS):
        for method, line in (('P6', LINES), ('P8', LINES8)):
            for pair, (s0, s1) in zip(('low_mid', 'mid_high'), ((0, 1), (1, 2))):
                previous = None
                tr, tt = [], []
                for ti, t in enumerate(T):
                    if t < 1.5:
                        continue
                    y = np.abs(line[ti, s0, d, 0]-line[ti, s1, d, 0])
                    p, props = find_peaks(y, prominence=0)
                    keep = (L[p] >= .75) & (L[p] < 15)
                    if previous is not None:
                        keep &= (L[p] > previous) & (L[p] <= previous+.75)
                    candidates = np.where(keep)[0]
                    assert len(candidates), (ray, method, pair, t)
                    pos = candidates[np.argmax(props['prominences'][candidates])]
                    k = p[pos]
                    previous = L[k]
                    h = np.array([2., 4/3, 8/9])/2**level(L[k], ray)
                    width = widths(y, k)
                    alpha, chi = LINES[ti, 2, d, (1, 0), k]
                    assert alpha >= 0 and chi > 0
                    rows.append(dict(branch='outgoing_chi_difference', ray=ray, method=method, pair=pair,
                                     time_M=t, radius_M=L[k], amplitude=y[k], fwhm_prominence_M=width,
                                     local_level=level(L[k], ray),
                                     **{s+'_dx_M': h[i] for i, s in enumerate(SCALES)},
                                     **{s+'_width_cells': width/h[i] for i, s in enumerate(SCALES)},
                                     lapse_high=alpha, chi_high=chi,
                                     gauge_speed_1p8_proxy=np.sqrt(1.8*alpha*chi),
                                     gauge_speed_2_proxy=np.sqrt(2*alpha*chi),
                                     light_speed_zero_shift_conformal_proxy=alpha*np.sqrt(chi)))
                    tr.append(L[k]);tt.append(t)
                    if method == 'P6' and pair == 'low_mid':
                        region = (L >= L[k]-.5) & (L <= L[k]+.5)
                        for comp, field in zip(INDEX, FIELDS):
                            a = np.sqrt(np.mean((line[ti, 0, d, comp, region]-line[ti, 1, d, comp, region])**2))
                            b = np.sqrt(np.mean((line[ti, 1, d, comp, region]-line[ti, 2, d, comp, region])**2))
                            a8 = np.sqrt(np.mean((LINES8[ti, 0, d, comp, region]-LINES8[ti, 1, d, comp, region])**2))
                            b8 = np.sqrt(np.mean((LINES8[ti, 1, d, comp, region]-LINES8[ti, 2, d, comp, region])**2))
                            packets.append(dict(time_M=t, ray=ray, field=field, packet_center_M=L[k], half_width_M=.5,
                                                low_mid_ray_rms=a, mid_high_ray_rms=b, self_order=order(a, b), self_order_P8=order(a8, b8)))
                for start in (1.5, 2.5, 4.5):
                    use = np.array(tt) >= start
                    fits.append(fit('outgoing_chi_'+pair, ray, method, np.array(tt)[use], np.array(tr)[use]))
                if method == 'P6' and pair == 'low_mid':
                    for a in (2., 4., 8.):
                        bound = a*np.sqrt(2) if ray == 'diagonal' else a
                        if not tr[0] <= bound <= tr[-1]:
                            continue
                        cross = float(np.interp(bound, tr, tt))
                        crossings.append(dict(ray=ray, boundary_axis_M=a, boundary_ray_M=bound,
                                              chi_difference_crossing_time_M=cross))
                        # Locate a burst by the largest half-M logarithmic growth of
                        # high |Ham| in ±0.35 M about this interface, near crossing.
                        mask = (L >= bound-.35) & (L <= bound+.35)
                        amp = np.max(np.abs(LINES[:, 2, d, 4, :][:, mask]), axis=1)
                        growth = np.log10(amp[1:]/amp[:-1])
                        candidate = np.where((T[1:] >= cross-.5) & (T[1:] <= cross+1.))[0]
                        ki = candidate[np.argmax(growth[candidate])]+1
                        ts = T[ki]
                        crossings[-1].update(ham_burst_time_M=ts, ham_growth_over_previous_output=amp[ki]/amp[ki-1],
                                             high_ham_band_peak=amp[ki])
                        duration = 1.5 if a != 4 else 2.
                        times, positions = [], []
                        for ti in np.where((T >= ts) & (T <= ts+duration))[0]:
                            y = np.abs(LINES[ti, 2, d, 4])
                            local = (L >= max(.6, bound-(T[ti]-ts)*.85-.5)) & (L <= bound+.1)
                            k = np.where(local)[0][np.argmax(y[local])]
                            times.append(T[ti]);positions.append(L[k])
                            rows.append(dict(branch=f'inward_Ham_from_{a:g}', ray=ray, method=method, pair='high',
                                             time_M=T[ti], radius_M=L[k], amplitude=y[k], fwhm_prominence_M=np.nan,
                                             local_level=level(L[k], ray),
                                             **{s+'_dx_M': np.nan for s in SCALES},
                                             **{s+'_width_cells': np.nan for s in SCALES}, lapse_high=np.nan, chi_high=np.nan,
                                             gauge_speed_1p8_proxy=np.nan, gauge_speed_2_proxy=np.nan,
                                             light_speed_zero_shift_conformal_proxy=np.nan))
                        fits.append(fit(f'inward_Ham_from_{a:g}', ray, 'P6', times, positions))
                        for si, scale in enumerate(SCALES):
                            for c, name in ((4, 'Ham'), (3, 'Theta')):
                                for method_width, source in (('P6', LINES), ('P8', LINES8)):
                                    y = np.abs(source[ki, si, d, c])
                                    p, _ = find_peaks(y)
                                    p = p[(L[p] >= bound-.5) & (L[p] <= bound+.35)]
                                    k = p[np.argmax(y[p])]
                                    w = widths(y, k)
                                    h = np.array([2., 4/3, 8/9])[si]/2**level(L[k], ray)
                                    pulses.append(dict(ray=ray, boundary_axis_M=a, boundary_ray_M=bound,
                                                       burst_time_M=ts, scale=scale, constraint=name, method=method_width,
                                                       peak_radius_M=L[k], peak_amplitude=y[k],
                                                       fwhm_prominence_M=w, local_level=level(L[k], ray), dx_M=h,
                                                       width_cells=w/h,
                                                       width_status='STENCIL_SWITCH_OR_SUBTWO_CELL' if abs(L[k]-bound) < w or w < 2*h else 'INTERPOLATED_WIDTH'))
    save('t5-front-tracks.csv', rows)
    save('t5-front-fits.csv', [r for r in fits if r])
    save('t5-interface-events.csv', crossings)
    save('t5-pulse-widths.csv', pulses)
    save('t5-packet-orders.csv', packets)


def figures():
    plt.style.use(['science', 'no-latex'])
    plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9})
    target = HERE/'figures'
    target.mkdir(exist_ok=True)

    def mapplot(ax, values, radii, title, limits, cmap='magma', ray=None):
        v = values if cmap == 'RdYlBu' else np.where(values > 0, np.log10(np.maximum(values, 1e-300)), np.nan)
        im = ax.pcolormesh(radii, T, v, shading='nearest', cmap=cmap, vmin=limits[0], vmax=limits[1], rasterized=True)
        for b in (2., 4., 8., 16.):
            if ray == 'diagonal':
                b *= np.sqrt(2)
            ax.axvline(b, color='0.6', lw=.65, ls='--')
        ax.set(xlim=(0, 12), ylim=(0, 10), xlabel='r / M', ylabel='t / M', title=title)
        fig.colorbar(im, ax=ax, shrink=.85)

    def finish(name):
        fig.savefig(target/(name+'.png'), dpi=180)
        fig.savefig(target/(name+'.pdf'), dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(12, 9), layout='constrained')
    for c, name in enumerate(('Hamiltonian', 'Momentum', 'GaussE')):
        mapplot(axes[c, 0], NATIVE[:, 2, 0, c], R, name+' high RMS; log10', (-12, -2))
        for j in range(2):
            mapplot(axes[c, j+1], order(COMMON[:, j, 0, c], COMMON[:, j+1, 0, c]), R,
                    name+(' low→mid p' if j == 0 else ' mid→high p'), (-1, 6), 'RdYlBu')
    finish('t5-constraints-spacetime')
    fig, axes = plt.subplots(3, 6, figsize=(17, 8), layout='constrained')
    for c, name in enumerate(('Ham', 'Mom', 'GaussE')):
        for d in range(6):
            mapplot(axes[c, d], order(COMMON[:, 1, d+1, c], COMMON[:, 2, d+1, c]), R,
                    f'{name}, {30*d}–{30*(d+1)} deg, p mid→high', (-1, 6), 'RdYlBu')
    finish('t5-angle-orders')
    fig, axes = plt.subplots(5, 4, figsize=(15, 14), layout='constrained')
    for c, name in enumerate(FIELDS):
        a, b = DIFF[:, 0, 0, c], DIFF[:, 1, 0, c]
        mapplot(axes[c, 0], a, R, name+' ||D low-mid||; log10', (-13, -3))
        mapplot(axes[c, 1], Q4*b, R, name+' 5.0625 ||D mid-high||', (-13, -3))
        mapplot(axes[c, 2], order(a, b), R, name+' self order', (-1, 6), 'RdYlBu')
        mapplot(axes[c, 3], ratio(F['mismatch'][:, 0, c], a), R, name+' ||Dlm-5.0625Dmh||/||Dlm||', (-3, 1))
    finish('t5-field-differences')
    fig, axes = plt.subplots(3, 4, figsize=(14, 9), layout='constrained')
    for c, comp in enumerate((4, 3, 0)):
        for d, ray in enumerate(RAYS):
            if comp in (4, 3):
                v = np.abs(LINES[:, 2, d, comp]);name = 'high |Ham|' if comp == 4 else 'high |Theta|'
            else:
                v = np.abs(LINES[:, 0, d, comp]-LINES[:, 1, d, comp]);name = '|D chi low-mid|'
            mapplot(axes[c, d], v, L, name+' '+ray, (-12, -3), ray=ray)
    finish('t5-rays-and-interfaces')
    fig, axes = plt.subplots(5, 4, figsize=(14, 14), layout='constrained')
    for c, comp in enumerate(INDEX):
        for d, ray in enumerate(RAYS):
            a = np.abs(LINES[:, 0, d, comp]-LINES[:, 1, d, comp])
            b = np.abs(LINES[:, 1, d, comp]-LINES[:, 2, d, comp])
            mapplot(axes[c, d], order(a, b), L, FIELDS[c]+' '+ray+' self p', (-1, 6), 'RdYlBu', ray)
    finish('t5-ray-field-orders')
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), layout='constrained')
    for c, comp in enumerate((1, 0)):
        for d, ray in enumerate(RAYS[:3]):
            v = np.abs(np.gradient(LINES[:, 2, d, comp], .5, axis=0))
            mapplot(axes[c, d], v, L, ('|dt lapse|' if c == 0 else '|dt chi|')+' '+ray, (-7, -1))
    finish('t5-temporal-gauge-change')
    fig, axes = plt.subplots(3, 3, figsize=(12, 8), layout='constrained')
    for c, comp in enumerate((4, 5, 7)):
        for d, ray in enumerate(RAYS[:3]):
            v = np.abs(F['outer'][:, 2, d, comp])
            mapplot(axes[c, d], v, F['outer_radii'], ('Ham', 'Mom1', 'GaussE')[c]+' '+ray+' log10', (-12, -3))
            axes[c, d].set_xlim(200, 256)
            axes[c, d].axvline(256, color='k')
    finish('t5-outer-boundary')
    track = [r for r in csv.DictReader((HERE/'t5-front-tracks.csv').open())
             if r['branch']=='outgoing_chi_difference' and r['ray']=='axis_plus' and r['method']=='P6' and r['pair']=='low_mid']
    tt = np.array([float(r['time_M']) for r in track])
    rr = np.array([float(r['radius_M']) for r in track])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    axes[0].plot(tt, rr, 'o-', label='tracked chi difference')
    for speed in (1., .87):
        axes[0].plot(tt, rr[0]+speed*(tt-tt[0]), '--', label=f'constant speed {speed:g}')
    axes[0].set(xlabel='t / M', ylabel='r / M', title='Outward error ridge');axes[0].legend()
    axes[1].plot(tt, np.gradient(rr, tt), 'o-', label='measured dr/dt')
    for key, name in (('gauge_speed_1p8_proxy', 'sqrt(1.8 alpha chi)'),
                      ('gauge_speed_2_proxy', 'sqrt(2 alpha chi), parameter comparison'),
                      ('light_speed_zero_shift_conformal_proxy', 'alpha sqrt(chi)')):
        axes[1].plot(tt, [float(r[key]) for r in track], label=name)
    axes[1].axhline(1., color='0.5', ls='--');axes[1].axhline(.87, color='0.7', ls=':')
    axes[1].set(xlabel='t / M', ylabel='coordinate speed / proxy', ylim=(0, 1.4), title='Missing shift/metric prevent mode identification')
    axes[1].legend(fontsize=8)
    finish('t5-front-speed')
    fig, axes = plt.subplots(3, 3, figsize=(12, 9), layout='constrained')
    for row, t in enumerate((2.5, 4.5, 8.5)):
        ti = int(round(2*t))
        tr = next(r for r in track if float(r['time_M']) == t)
        crest = float(tr['radius_M'])
        window = (L >= crest-.7) & (L <= crest+.7)
        for si, scale in enumerate(SCALES):
            axes[row, 0].plot(L[window], LINES[ti, si, 0, 4, window], label=scale)
            axes[row, 1].plot(L[window], LINES[ti, si, 0, 3, window], label=scale)
        axes[row, 2].plot(L[window], LINES[ti, 0, 0, 0, window]-LINES[ti, 1, 0, 0, window], label='D low-mid')
        axes[row, 2].plot(L[window], Q4*(LINES[ti, 1, 0, 0, window]-LINES[ti, 2, 0, 0, window]), label='5.0625 D mid-high')
        for col, title in enumerate(('Ham reconstructed', 'Theta reconstructed', 'signed chi differences')):
            axes[row, col].set(xlim=(crest-.7, crest+.7), xlabel='r / M', title=f't={t:g} M: {title}')
            axes[row, col].legend(fontsize=8)
    finish('t5-front-packets')
    fig, axes = plt.subplots(3, 3, figsize=(12, 9), layout='constrained')
    for c, name in enumerate(('Ham', 'Mom', 'GaussE')):
        for si, scale in enumerate(SCALES):
            mapplot(axes[c, si], F['global_native'][:, si, 0, c], F['global_radii'],
                    scale+' '+name+' global log10 RMS', (-14, -2))
            axes[c, si].set_xlim(0, 364)
            for a in (32, 64, 256):axes[c, si].axvline(a, color='0.6', lw=.65, ls='--')
    finish('t5-global-constraints')
    fig, axes = plt.subplots(5, 3, figsize=(12, 14), layout='constrained')
    for c, field in enumerate(FIELDS):
        a, b = F['global_diff'][:, 0, 0, c], F['global_diff'][:, 1, 0, c]
        for j, (v, label) in enumerate(((a, '||D low-mid||'), (Q4*b, '5.0625 ||D mid-high||'), (order(a, b), 'self order'))):
            mapplot(axes[c, j], v, F['global_radii'], field+' global '+label,
                    (-1, 6) if j == 2 else (-14, -3), 'RdYlBu' if j == 2 else 'magma')
            axes[c, j].set_xlim(0, 364)
            for boundary in (32, 64, 256):axes[c, j].axvline(boundary, color='0.6', lw=.65, ls='--')
    finish('t5-global-field-differences')


if __name__ == '__main__':
    started = time.monotonic()
    csvs()
    tracks()
    figures()
    save('t5-report-performance.csv', [dict(wall_seconds=time.monotonic()-started,
                                           peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)])
    print('T5 CSVs, tracks, fits and eleven PNG/PDF figure pairs written.')
