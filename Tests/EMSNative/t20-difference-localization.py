#!/usr/bin/env python3
"""Localize exp-0023 numerical self-differences; T19's FAIL is not reopened.

Usage: t20-difference-localization.py COLLECTED_PRODUCTION_ROOT
No reconstruction from static data, new evolution, or replacement acceptance rule.
Outputs are t20-* beside this script. Only cached I8/I10 numerical profiles are used.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import resource
import shutil
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/ems-t20-matplotlib')
import numpy as np

HERE = Path(__file__).resolve().parent
FIELDS = 'A11 A22 Aww Gamma1 Gamma2 chi K lapse phi Ex'.split()
RAYS = ('axis', 'equator', 'diagonal')
LEGS = ('E-low', 'E-mid', 'E-high')
INPUTS, OUTPUTS = set(), set()
TIMES = np.arange(13) * .875
ANCHORS = (.0055, .00608, .0075, .02, .3, 3.5)
REGIONS = {'near': (.005, .013671875), 'inner': (.013671875, 1.),
           'exterior_near': (1., 8.), 'outer_cache': (8., 20.)}
BOX = re.compile(r'^\d+: \(\((-?\d+),(-?\d+)\) \((-?\d+),(-?\d+)\)')


def digest(path):
    d = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            d.update(chunk)
    return d.hexdigest()


def rows(path):
    INPUTS.add(path)
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def output(name, data):
    p = HERE / ('t20-' + name + '.csv')
    cols = list(dict.fromkeys(k for r in data for k in r))
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(data)
    OUTPUTS.add(p)


def text_output(name, value):
    p = HERE / ('t20-' + name)
    p.write_text(value)
    OUTPUTS.add(p)


def num(r, k):
    return float(r[k]) if r.get(k) not in ('', None) else math.nan


def rms(a, weights):
    return float(np.sqrt(np.sum(a*a*weights)/np.sum(weights)))


def clock_index(t):
    return int(round(t/.875))


def layouts(root):
    """Verify initial boxes against every collected production restart layout.

    Finder pouts do not print checkpoint boxes. Read only the restart prefix of
    production pouts, never the millions of fine-step lines that follow it.
    """
    census = root.parent/'submission/phase1-evidence/smoke/census'
    faces, audits, seams = [], [], []
    initial, domains = {}, {}
    for leg in LEGS:
        c = rows(census/(leg+'-dense')/'layout-census.csv')
        b = rows(census/(leg+'-dense')/'box-census.csv')
        initial[leg] = {l: sorted(tuple(int(r[k]) for k in
            ('lo_i', 'lo_j', 'hi_i', 'hi_j')) for r in b if int(r['level'])==l)
            for l in range(15)}
        domains[leg] = c
        for r in c:
            l, h = int(r['level']), num(r, 'dx_M')
            face = num(r, 'x_max_M')
            faces.append(dict(leg=leg, level=l, h_M=h, x_face_M=face,
                y_face_M=num(r, 'y_max_M'), axis_ray_face_M=face,
                equator_ray_face_M=num(r, 'y_max_M'),
                diagonal_ray_face_M=math.sqrt(2)*min(face, num(r, 'y_max_M')),
                boxes=int(r['boxes'])))
            # Same-level box faces differ even when their union is identical.
            sx = sorted({(v*h-336.) for box in initial[leg][l]
                for v in (box[0],box[2]+1) if -face < v*h-336. < face})
            sy = sorted({v*h for box in initial[leg][l]
                for v in (box[1],box[3]+1) if 0 < v*h < num(r, 'y_max_M')})
            seams.append(dict(leg=leg, level=l,
                internal_x_faces_M=';'.join(f'{v:.12g}' for v in sx),
                internal_y_faces_M=';'.join(f'{v:.12g}' for v in sy)))
        paths = [root/leg/'pout.0.gz'] + sorted((root/leg/'segments').glob('attempt-*/pout.0.gz'))
        for p in paths:
            INPUTS.add(p)
            levels, reading, current, t, disabled = [], False, [], None, set()
            with gzip.open(p, 'rt') as f:
                for i, line in enumerate(f):
                    m = re.match(r'regrid_interval_(\d+)\s*:\s*(\d+)', line)
                    if m and int(m[2]) == 0:
                        disabled.add(int(m[1]))
                    if 'read time =' in line:
                        t = float(line.split('=')[1])
                    if 'read level domain:' in line:
                        reading, current = True, []
                        continue
                    if reading:
                        m = BOX.match(line)
                        if m:
                            current.append(tuple(map(int, m.groups())))
                        else:
                            levels.append(sorted(current))
                            reading = False
                            if len(levels) == 15:
                                break
                    if i >= 4999:  # Initial evolutions have no checkpoint-read prefix.
                        break
            matched = len(levels) == 15 and all(levels[l] == initial[leg][l] for l in range(15))
            audits.append(dict(leg=leg, source=str(p), time_M=t,
                printed_levels=len(levels), exact_census_boxes_match=matched if levels else '',
                zero_regrid_intervals=len(disabled),
                interpretation='restart_exact_layout' if levels else 'initial_run_no_read_layout'))
            if levels:
                assert matched, f'Changed boxes: {p}'
            # Parameters are printed in the production binary's checkpoint header.
            if levels:
                assert len(disabled)==15, f'Nonzero/unprinted regrid setting: {p}'
    for l in range(15):
        f = [r for r in faces if r['level']==l]
        assert max(r['x_face_M'] for r in f)-min(r['x_face_M'] for r in f)<1e-12
        assert max(r['y_face_M'] for r in f)-min(r['y_face_M'] for r in f)<1e-12
    support = rows(root/'battery/common-support.csv')
    for key in {(r['mask'],r['category']) for r in support}:
        g = [r for r in support if (r['mask'],r['category'])==key]
        assert len(g)==13 and len({(r['cells'],r['weight'],r['status']) for r in g})==1
    for leg in LEGS:
        lb = list((root/leg).glob('qualified/step*/n48/LB.txt'))
        audits.append(dict(leg=leg, source='qualified/*/n48/LB.txt',
            printed_levels='', exact_census_boxes_match='', zero_regrid_intervals='',
            files=len(lb), interpretation='rank_loads_only_no_box_coordinates' if lb else 'not_collected'))
        if lb:
            INPUTS.add(lb[0])
    output('hierarchy-faces', faces)
    output('hierarchy-seams', seams)
    output('hierarchy-audit', audits)
    return faces, audits


def characteristic_witness(cache, names):
    # Witness first: in Q[a,b,beta,lambda], det(M-lambda I) equals
    # (lambda+beta)^2-ab for M=[[-beta,a],[b,-beta]]. A zero polynomial
    # proves the characteristic equation; positive a,b give the real sqrt branch.
    import sympy as sp
    a,b,beta,lam = sp.symbols('a b beta lambda')
    residual = sp.Poly(sp.expand(sp.Matrix([[-beta,a],[b,-beta]]).det() -
        (beta**2-a*b)), a,b,beta)
    char = sp.Poly(sp.expand((sp.Matrix([[-beta,a],[b,-beta]])-lam*sp.eye(2)).det()
        - ((lam+beta)**2-a*b)), a,b,beta,lam)
    assert residual.is_zero and char.is_zero
    errors = []
    for ray in RAYS:
        n = np.array((1.,0.) if ray=='axis' else (0.,1.) if ray=='equator'
                     else (1/math.sqrt(2),1/math.sqrt(2)))
        for k in (0,4,11,12):
            v = cache[k,ray]['H8']
            for j in (0,256,1024,2047):
                h = np.array([[v[names['h11'],j],v[names['h12'],j]],
                              [v[names['h12'],j],v[names['h22'],j]]])
                assert np.linalg.eigvalsh(h).min()>0
                hn = n@np.linalg.inv(h)@n
                alpha, chi = v[names['lapse'],j],v[names['chi'],j]
                assert alpha>0 and chi>0
                bn = n@v[[names['shift1'],names['shift2']],j]
                for aa,bb in ((alpha,alpha*chi*hn),(1.8*alpha,chi*hn),(.75,4*hn/3)):
                    numeric = np.sort(np.linalg.eigvals([[-bn,aa],[bb,-bn]]))
                    formula = np.array([-bn-math.sqrt(aa*bb),-bn+math.sqrt(aa*bb)])
                    errors.append(float(np.max(np.abs(numeric-formula))))
    evidence = dict(status='PROVED', claim='2x2 characteristic polynomial only',
        ring='Q[a,b,beta,lambda]', exact_witness='det(M-lambda I)-(lambda+beta)^2+a*b = 0',
        assumptions='real beta; positive a,b for the ordered real square-root eigenvalues',
        poles='none; metric inverse requires positive definite conformal metric',
        exact_coefficients=[int(v) for v in char.coeffs()], backend='SymPy '+sp.__version__,
        numeric_falsifier='independent numpy.linalg.eigvals, Float64, no RNG',
        numerical_samples=len(errors), sample_clocks=[0,4,11,12],
        sample_radius_indices=[0,256,1024,2047], all_rays=True,
        worst_eigenvalue_absolute_residual=max(errors),
        precision='Float64 falsifier; exact symbolic polynomial closure supplies proof',
        scope='frozen principal subblocks; full coupled gauge characteristic spectrum NOT-CLOSED')
    assert max(errors)<1e-12
    text_output('cas-witness.json', json.dumps(evidence, indent=2)+'\n')


def fit_history(t, y):
    """Descriptive fits, not an instability/period certificate.

    Fit a sinusoid plus linear trend at periods with >=2 cycles; also report
    log-amplitude slopes. Their fit quality and sign changes remain explicit.
    """
    y=np.asarray(y)
    scale=max(np.max(np.abs(y)),1e-300)
    z=y/scale
    base=np.column_stack((np.ones(len(t)),t))
    sse0=np.sum((z-base@np.linalg.lstsq(base,z,rcond=None)[0])**2)
    best=(math.inf,math.nan)
    for period in np.linspace(2*(t[1]-t[0]),(t[-1]-t[0])/2,180):
        x=np.column_stack((base,np.cos(2*np.pi*t/period),np.sin(2*np.pi*t/period)))
        sse=np.sum((z-x@np.linalg.lstsq(x,z,rcond=None)[0])**2)
        if sse<best[0]: best=(sse,period)
    signs=np.sign(y)
    crossings=int(np.count_nonzero(signs[1:]*signs[:-1]<0))
    nz=np.abs(y)>scale*1e-14
    x=t[nz]; logy=np.log(np.abs(y[nz]))
    slope, intercept=np.polyfit(x,logy,1) if len(x)>2 else (math.nan,math.nan)
    denom=np.sum((logy-np.mean(logy))**2)
    r2=1-np.sum((logy-intercept-slope*x)**2)/denom if denom>0 else math.nan
    return dict(sign_changes=crossings, candidate_period_M=best[1],
        sinusoid_improvement_over_linear=1-best[0]/sse0 if sse0>1e-30 else math.nan,
        log_abs_slope_per_M=slope, log_abs_fit_R2=r2,
        interpretation='descriptive_13_clock_fit_no_unique_period_or_exponential_claim')


def cached_analysis(root):
    gauge_coefficients=[]
    for leg in LEGS:
        p=root/leg/'qualified/step000000/n48/params.txt'
        INPUTS.add(p)
        m=re.search(r'^shift_Gamma_coeff\s*=\s*(\S+)',p.read_text(),re.M)
        assert m, 'Missing shift_Gamma_coeff: '+str(p)
        gauge_coefficients.append(float(m[1]))
    assert len(set(gauge_coefficients))==1 and gauge_coefficients[0]==.75
    # Registered masks remain unchanged; extra descriptive partitions have distinct names.
    p=HERE/'t14d-stageD-masks.csv'
    INPUTS.add(p)
    registered=list(csv.DictReader(p.read_text().replace('\\r\\n','\n').splitlines()))
    regions=dict(REGIONS)
    regions.update({r['mask']:(num(r,'R_min_M'),num(r,'R_max_M')) for r in registered
                    if r['mask']!='outer_boundary_shell'})
    cache={}
    for k in range(13):
        for ray in RAYS:
            p=root/f'battery/ray-cache-clock{k}-{ray}.npz'
            INPUTS.add(p)
            with np.load(p,allow_pickle=False) as f:
                cache[k,ray]={name:f[name] for name in f.files}
    radii=cache[0,'axis']['radii']
    names={str(f):i for i,f in enumerate(cache[0,'axis']['fields'])}
    assert len(radii)==2048 and abs(radii[0]-.005)<1e-14 and abs(radii[-1]-20)<1e-12
    for c in cache.values():
        assert np.array_equal(c['radii'],radii)
        assert np.array_equal(c['fields'],cache[0,'axis']['fields'])
    characteristic_witness(cache,names)
    weights=np.gradient(radii)
    stats=defaultdict(list); fixed=[]; history=[]; speeds=[]; packet=[]
    for k in range(13):
        for ray in RAYS:
            c=cache[k,ray]
            n=np.array((1.,0.) if ray=='axis' else (0.,1.) if ray=='equator'
                       else (1/math.sqrt(2),1/math.sqrt(2)))
            face_factor=1 if ray!='diagonal' else math.sqrt(2)
            face_r=np.array([56/2**(l-1)*face_factor for l in range(1,15)])
            for field in FIELDS:
                i=names[field]
                for pair,a,b in (('low_mid','L','M'),('mid_high','M','H')):
                    d=c[a+'8'][i]-c[b+'8'][i]
                    spread=np.abs(d-(c[a+'10'][i]-c[b+'10'][i]))
                    floor=128*np.finfo(float).eps*max(1.,max(float(np.max(abs(c[v+'8'][i]))) for v in 'LMHC'))
                    for region,(lo,hi) in regions.items():
                        mask=(radii>=lo)&(radii<=hi)
                        rr,dd,ww=radii[mask],d[mask],weights[mask]
                        peak=int(np.argmax(np.abs(dd)))
                        energy=dd**2*ww; total=np.sum(energy)
                        fractions=np.cumsum(energy)/total if total>0 else np.zeros(len(energy))
                        q=[rr[min(np.searchsorted(fractions,v),len(rr)-1)] for v in (.1,.5,.9)]
                        h0=.875 if pair=='low_mid' else 7/12
                        level=np.zeros(len(rr),dtype=int)
                        for l,fr in enumerate(face_r,1):level[rr<fr]=l
                        hh=h0/2**level
                        distance=np.min(np.abs(rr[:,None]-face_r[None,:]),axis=1)
                        band=distance<=4*hh  # Descriptive four-local-cell band, no acceptance change.
                        band_weight=np.sum(ww[band])/np.sum(ww)
                        band_energy=float(np.sum(energy[band])/total) if total>0 else math.nan
                        qualified=np.abs(dd)>5*(spread[mask]+floor)
                        s=dd[qualified]; zeros=int(np.count_nonzero(s[1:]*s[:-1]<0))
                        stats[region].append(dict(clock=k,time_M=TIMES[k],ray=ray,field=field,pair=pair,
                            radial_dr_RMS=rms(dd,ww),spread_dr_RMS=rms(spread[mask],ww),
                            temporal_dr_RMS=rms(c['H8'][i,mask]-c['C8'][i,mask],ww),
                            peak_radius_M=rr[peak],signed_peak=dd[peak],
                            peak_at_selection_edge=peak in (0,len(rr)-1),
                            energy_R10_M=q[0],energy_R50_M=q[1],energy_R90_M=q[2],
                            significant_sign_changes=zeros,significant_points=int(np.sum(qualified)),
                            nearest_refinement_face_to_peak_M=face_r[np.argmin(abs(face_r-rr[peak]))],
                            face_band_energy_fraction=band_energy,face_band_dr_fraction=band_weight,
                            face_band_enrichment=band_energy/band_weight if band_weight else math.nan))
                for target in ANCHORS:
                    j=int(np.argmin(abs(radii-target)))
                    fixed.append(dict(clock=k,time_M=TIMES[k],ray=ray,field=field,
                        target_radius_M=target,radius_M=radii[j],
                        low=c['L8'][i,j],mid=c['M8'][i,j],high=c['H8'][i,j],
                        low_mid=c['L8'][i,j]-c['M8'][i,j],mid_high=c['M8'][i,j]-c['H8'][i,j],
                        temporal=c['H8'][i,j]-c['C8'][i,j],
                        spread_low_mid=abs((c['L8'][i,j]-c['M8'][i,j])-(c['L10'][i,j]-c['M10'][i,j])),
                        spread_mid_high=abs((c['M8'][i,j]-c['H8'][i,j])-(c['M10'][i,j]-c['H10'][i,j]))))
            # Compare ray-normal characteristic estimates at anchors and the travelling Γ peak.
            for target in ANCHORS+(max(.005,min(20,TIMES[k])),):
                j=int(np.argmin(abs(radii-target)));v=c['H8']
                h=np.array([[v[names['h11'],j],v[names['h12'],j]],
                            [v[names['h12'],j],v[names['h22'],j]]])
                assert np.linalg.eigvalsh(h).min()>0
                hn=n@np.linalg.inv(h)@n
                alpha,chi=v[names['lapse'],j],v[names['chi'],j]
                assert alpha>0 and chi>0
                bn=n@v[[names['shift1'],names['shift2']],j]
                speeds.append(dict(clock=k,time_M=TIMES[k],ray=ray,radius_M=radii[j],
                    beta_n=bn,light_out=-bn+alpha*math.sqrt(chi*hn),
                    light_in=-bn-alpha*math.sqrt(chi*hn),
                    lapse_subblock_out=-bn+math.sqrt(1.8*alpha*chi*hn),
                    Gamma_long_subblock_out=-bn+math.sqrt(hn),
                    Gamma_trans_subblock_out=-bn+math.sqrt(.75*hn)))
            # Follow the outward leading packet in a declared kinematic corridor R/t=0.8..1.2,
            # separately from the later interior packet. This is exploratory localization.
            if k>=1:
                mask=(radii>=max(.005,.8*TIMES[k]))&(radii<=min(20.,1.2*TIMES[k]))
                for field in FIELDS:
                    i=names[field]
                    for pair,a,b in (('low_mid','L','M'),('mid_high','M','H')):
                        d=c[a+'8'][i]-c[b+'8'][i]
                        rr=radii[mask];dd=d[mask]
                        if not len(rr):continue
                        j=int(np.argmax(np.abs(dd)))
                        packet.append(dict(clock=k,time_M=TIMES[k],ray=ray,field=field,pair=pair,
                            radius_M=rr[j],signed_value=dd[j],edge_peak=j in (0,len(rr)-1),
                            corridor_low_M=rr[0],corridor_high_M=rr[-1],
                            note='exploratory_R_over_t_0.8_to_1.2_not_unique_packet_identification'))
    for ray in RAYS:
        for field in FIELDS:
            for target in ANCHORS:
                g=[r for r in fixed if r['ray']==ray and r['field']==field and r['target_radius_M']==target]
                for key in ('low_mid','mid_high'):
                    history.append(dict(ray=ray,field=field,pair=key,radius_M=g[0]['radius_M'],
                        **fit_history(TIMES[1:],np.array([r[key] for r in g])[1:])))
    growth=[];extrema=[]
    for region,rs in stats.items():
        for ray in RAYS:
            for field in FIELDS:
                for pair in ('low_mid','mid_high'):
                    g=[r for r in rs if r['ray']==ray and r['field']==field and r['pair']==pair and r['clock']>0]
                    for window,lo,hi in (('all_positive',.875,10.5),('early',.875,5.25),('late',5.25,10.5)):
                        gg=[r for r in g if lo<=r['time_M']<=hi]
                        t=np.array([r['time_M'] for r in gg]);y=np.array([r['radial_dr_RMS'] for r in gg])
                        logy=np.log(np.maximum(y,1e-300));slope,intercept=np.polyfit(t,logy,1)
                        denom=np.sum((logy-logy.mean())**2)
                        r2=1-np.sum((logy-intercept-slope*t)**2)/denom if denom else math.nan
                        growth.append(dict(region=region,ray=ray,field=field,pair=pair,window=window,
                            start_M=lo,end_M=hi,log_RMS_slope_per_M=slope,fit_R2=r2,
                            decreasing_intervals=int(np.count_nonzero(np.diff(y)<0)),
                            note='descriptive_log_fit_not_an_exponential_instability_certificate'))
    for ray in RAYS:
        for field in FIELDS:
            for target in ANCHORS:
                g=[r for r in fixed if r['ray']==ray and r['field']==field and r['target_radius_M']==target]
                for leg,key in zip(LEGS,('low','mid','high')):
                    y=np.array([r[key] for r in g])
                    changes=np.flatnonzero(np.diff(y)[1:]*np.diff(y)[:-1]<0)+1
                    for kind in ('maximum','minimum'):
                        ix=[j for j in changes if (y[j]>y[j-1])==(kind=='maximum')]
                        for j,jj in zip(ix,ix[1:]):
                            extrema.append(dict(ray=ray,field=field,leg=leg,radius_M=g[0]['radius_M'],
                                kind=kind,first_M=TIMES[j],next_M=TIMES[jj],
                                extremum_spacing_M=TIMES[jj]-TIMES[j],
                                neighbour_bracket_spacing_low_M=max(0,TIMES[jj]-TIMES[j]-1.75),
                                neighbour_bracket_spacing_high_M=TIMES[jj]-TIMES[j]+1.75,
                                note='conditional_on_one_corresponding_extremum_per_neighbour_bracket_not_unique_frequency'))
    for region,rs in stats.items():output('localization-'+region,rs)
    output('fixed-radius',fixed)
    output('history-fits',history)
    output('characteristic-estimates',speeds)
    output('leading-packet',packet)
    # Keep each regenerable ledger small.
    output('growth-fits',growth)
    output('sampled-extrema',extrema)
    return cache,names,radii,stats,fixed,history,speeds,packet


def additional_histories(root):
    data=[]; coverage=[]
    for leg in LEGS:
        p=root/leg/'qualified-horizons.csv'
        rs=rows(p)
        g=[r for r in rs if int(r['N_theta'])==96]
        times=sorted({num(r,'time_M') for r in g})
        coverage.append(dict(leg=leg,source='qualified-horizons.csv N96',clocks=len(times),
            first_M=times[0],last_M=times[-1],cadence_M=times[1]-times[0],
            geometry='moving numerical horizon, NOT fixed-radius profiles'))
        for r in g:
            data.append(dict(leg=leg,time_M=num(r,'time_M'),quantity='phi_mean_on_numerical_horizon',
                             value=num(r,'phi_mean')))
        p=root/leg/'constraint_norms.dat'
        INPUTS.add(p)
        a=np.loadtxt(p)
        for t,h,m in a:
            data += [dict(leg=leg,time_M=t,quantity='global_Ham_L2',value=h),
                     dict(leg=leg,time_M=t,quantity='global_Mom1_Mom2_L2',value=m)]
        coverage.append(dict(leg=leg,source='constraint_norms.dat',clocks=len(a),
            first_M=a[0,0],last_M=a[-1,0],cadence_M=a[1,0]-a[0,0],
            geometry='native whole composite grid, NOT registered masks'))
        # Every coarse checkpoint has a finder candidate radius history as well.
        rs=rows(root/leg/'finder-values.csv')
        for r in rs:
            if '<r>' in r:
                data.append(dict(leg=leg,time_M=num(r,'t'),quantity='candidate_mean_radius',value=num(r,'<r>')))
    output('finer-histories',data)
    output('history-coverage',coverage)
    return data,coverage


def horizon_radii(root):
    result={}
    table=[]
    for leg,multiple in zip(LEGS,(4,6,9)):
        for k in range(13):
            p=root/leg/f'qualified/step{k*multiple:06d}/n96/shape-0-2.dat'
            INPUTS.add(p)
            a=np.loadtxt(p).reshape(-1)
            assert abs(a[0]-TIMES[k])<1e-9
            result[leg,k]=float(np.mean(a[2:]))
            table.append(dict(leg=leg,clock=k,time_M=TIMES[k],radius_min_M=np.min(a[2:]),
                              radius_mean_M=result[leg,k],radius_max_M=np.max(a[2:])))
    output('numerical-horizon-radii',table)
    return result


def witnesses(root,stats):
    # Preserve the exact T19 ledger, and check it against the source table keys.
    rs=rows(HERE/'t19-failing-orders.csv')
    source={}
    for dataset,file in (('volume','self-differences.csv'),('ray','ray-self-differences.csv')):
        for r in rows(root/'battery'/file):
            source[dataset,r['clock_index'],r['mask'],r['category'],r['field']]=r
    for r in rs:
        if r['kind']=='self_difference':
            s=source[r['dataset'],r['clock_index'],r['mask'],r['category'],r['field']]
            for key in ('difference_low_mid','difference_mid_high','measured_order','order_interval_high',
                        'spatial_significance','temporal_high_dt2'):
                assert r[key]==s[key]
            assert num(r,'order_interval_high')<0 and num(r,'spatial_significance')>5
            assert max(num(r,'temporal_fraction_of_low_mid'),num(r,'temporal_fraction_of_mid_high'))<=.2
    roundtrip=[]
    for k in range(13):
        for ray in RAYS:
            errors=[]
            for mask in stats:
                if ('ray',str(k),mask,'ray_'+ray,FIELDS[0]) not in source:continue
                for field in FIELDS:
                    s=source['ray',str(k),mask,'ray_'+ray,field]
                    for pair in ('low_mid','mid_high'):
                        a=next(r for r in stats[mask] if (r['clock'],r['ray'],r['field'],r['pair'])==(k,ray,field,pair))
                        for key,skey in (('radial_dr_RMS','difference_'+pair),
                                         ('spread_dr_RMS','I8_I10_'+pair+'_uncertainty'),
                                         ('temporal_dr_RMS','temporal_high_dt2')):
                            expected=num(s,skey);actual=a[key]
                            relative=abs(expected-actual)/max(abs(expected),1e-300)
                            assert relative<2e-12 or abs(expected-actual)<1e-28, (k,ray,mask,field,pair,skey)
                            errors.append(relative)
            roundtrip.append(dict(clock=k,ray=ray,norms_checked=len(errors),
                max_relative_roundtrip_error=max(errors),status='PASS_cache_reproduces_collected_ray_table'))
    output('cache-roundtrip',roundtrip)
    summary=[]
    for key,count in sorted(Counter((r['clock_index'],r['time_M'],r['dataset'],r['mask'],r['category'])
                                   for r in rs).items(),key=lambda v:(int(v[0][0]),v[0][2:])):
        summary.append(dict(clock=key[0],time_M=key[1],dataset=key[2],mask=key[3],category=key[4],witnesses=count))
    output('failing-witness-localization',summary)
    selected=[]
    for (dataset,k,mask,cat,field),r in source.items():
        if dataset=='volume' and cat=='all' and ((mask=='horizon' and field=='A11') or
            (mask in ('far','exterior_wake','receiving_side_3p5') and field in FIELDS)):
            selected.append(dict(dataset=dataset,**{key:r[key] for key in
                ('clock_index','time_M','mask','category','field','difference_low_mid','difference_mid_high',
                 'temporal_high_dt2','spatial_significance','measured_order','order_interval_low','order_interval_high','status')}))
    output('mask-histories',selected)
    return rs,selected


def figures(cache,names,radii,rh,finer):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import SymLogNorm
    import scienceplots
    plt.style.use(['science','no-latex'])
    plt.rcParams.update({'font.size':8,'axes.labelsize':8,'legend.fontsize':6})
    p=HERE/'t20-signed-ray-atlas.pdf'
    with PdfPages(p) as pdf:
        for k in range(13):
            for group in (FIELDS[:5],FIELDS[5:]):
                fig,axes=plt.subplots(5,3,figsize=(10,11),layout='constrained')
                for row,field in enumerate(group):
                    i=names[field]
                    for col,ray in enumerate(RAYS):
                        ax=axes[row,col];c=cache[k,ray]
                        d1,d2=c['L8'][i]-c['M8'][i],c['M8'][i]-c['H8'][i]
                        spread=np.maximum(abs(d1-(c['L10'][i]-c['M10'][i])),
                                          abs(d2-(c['M10'][i]-c['H10'][i])))
                        ax.plot(radii,d1,lw=.7,label='low - mid')
                        ax.plot(radii,d2,lw=.7,label='mid - high')
                        ax.fill_between(radii,-5*spread,5*spread,color='gray',alpha=.15,label='5 x I8/I10 spread')
                        magnitude=max(float(np.max(abs(np.r_[d1,d2]))),1e-15)
                        upper=10.**math.ceil(math.log10(magnitude))
                        threshold=max(upper*1e-4,1e-15)
                        ax.set(xscale='log',yscale='symlog',xlim=(.005,20),
                               ylabel=field if col==0 else '',title=ray if row==0 else '')
                        ax.set_yscale('symlog',linthresh=threshold)
                        ax.set_ylim(-upper,upper)
                        ax.set_yticks([-upper,0,upper] if threshold>=upper/100 else
                                      [-upper,-upper/100,0,upper/100,upper])
                        ax.axhline(0,color='gray',lw=.4)
                        factor=math.sqrt(2) if ray=='diagonal' else 1
                        for l in range(3,15):ax.axvline(56/2**(l-1)*factor,color='gray',lw=.4,alpha=.45)
                        for leg in LEGS:ax.axvline(rh[leg,k],color='black',lw=.6,alpha=.6)
                        if row==4:ax.set_xlabel('Numerical coordinate radius R/M')
                axes[0,0].legend(loc='best')
                fig.suptitle(f't = {TIMES[k]:g} M; black: numerical N96 horizons; gray: common refinement faces\n'
                             'Signed differences; cache starts at R=0.005 M; symlog linear band is for display only')
                pdf.savefig(fig);plt.close(fig)
    OUTPUTS.add(p)
    p=HERE/'t20-space-time.pdf'
    with PdfPages(p) as pdf:
        for field in FIELDS:
            fig,axes=plt.subplots(3,2,figsize=(10,9),layout='constrained')
            for row,ray in enumerate(RAYS):
                i=names[field]
                for col,(a,b) in enumerate((('L','M'),('M','H'))):
                    ax=axes[row,col]
                    v=np.array([cache[k,ray][a+'8'][i]-cache[k,ray][b+'8'][i] for k in range(13)])
                    scale=max(np.max(abs(v)),1e-15)
                    im=ax.pcolormesh(radii,TIMES,v,shading='nearest',cmap='RdBu_r',
                        norm=SymLogNorm(linthresh=scale*1e-3,vmin=-scale,vmax=scale),rasterized=True)
                    ax.set(xscale='log',xlim=(.005,20),ylim=(-.4375,10.9375),
                           ylabel=f'{ray}: t/M',title=a+' - '+b)
                    factor=math.sqrt(2) if ray=='diagonal' else 1
                    for l in range(3,15):ax.axvline(56/2**(l-1)*factor,color='gray',lw=.4,alpha=.4)
                    ax.plot([rh['E-high',k] for k in range(13)],TIMES,'k-',lw=.7)
                    fig.colorbar(im,ax=ax,label='signed '+field+' difference',
                        ticks=[-scale,-scale/100,0,scale/100,scale],format='%.0e')
                    ax.set_xlabel('R/M')
            fig.suptitle(field+' numerical self-differences: 13 stored ray clocks, no added time samples')
            pdf.savefig(fig);plt.close(fig)
    OUTPUTS.add(p)
    p=HERE/'t20-time-traces.pdf'
    with PdfPages(p) as pdf:
        for ray in RAYS:
            for target in (.0075,3.5):
                j=int(np.argmin(abs(radii-target)))
                fig,axes=plt.subplots(5,2,figsize=(10,10),layout='constrained')
                for ax,field in zip(axes.flat,FIELDS):
                    i=names[field]
                    for a,b,label in (('L','M','low - mid'),('M','H','mid - high')):
                        ax.plot(TIMES,[cache[k,ray][a+'8'][i,j]-cache[k,ray][b+'8'][i,j]
                                      for k in range(13)],'.-',lw=.8,label=label)
                    ax.set(title=field,xlabel='t/M');ax.axhline(0,color='gray',lw=.4)
                axes[0,0].legend()
                fig.suptitle(f'{ray}, fixed R={radii[j]:.7g} M; signed differences at the stored clocks')
                pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
        quantities=('phi_mean_on_numerical_horizon','global_Ham_L2','global_Mom1_Mom2_L2')
        for ax,q in zip(axes.flat,quantities):
            for leg in LEGS:
                g=sorted((r for r in finer if r['leg']==leg and r['quantity']==q),key=lambda r:r['time_M'])
                ax.plot([r['time_M'] for r in g],[r['value'] for r in g],label=leg,lw=.8)
            ax.set(title=q,xlabel='t/M');ax.legend()
        axes[1,1].axis('off')
        axes[1,1].text(0,1,'All collected finer scalar histories.\nHorizon scalar is on a moving surface.\n'
            'Whole-grid constraint norms include the puncture.\nNo finer fixed-radius field profiles were collected.',va='top')
        pdf.savefig(fig);plt.close(fig)
    OUTPUTS.add(p)
    p=HERE/'t20-exterior-packet.png'
    fig,axes=plt.subplots(2,2,figsize=(9,6),layout='constrained')
    for row,field in enumerate(('Gamma1','lapse')):
        i=names[field]
        for col,(a,b) in enumerate((('L','M'),('M','H'))):
            ax=axes[row,col]
            visible=(radii>=1)&(radii<=7)
            for k in (3,4,5,6):
                d=cache[k,'axis'][a+'8'][i]-cache[k,'axis'][b+'8'][i]
                ax.plot(radii[visible],d[visible],label=f'{TIMES[k]:g} M',lw=.8)
            ax.set(xlim=(1,7),xlabel='R/M',ylabel=field+' '+a+' - '+b);ax.legend()
            for face in (1.75,3.5,7):ax.axvline(face,color='gray',lw=.5)
    fig.savefig(p,dpi=600);OUTPUTS.add(p)
    q=p.with_suffix('.pdf');fig.savefig(q);OUTPUTS.add(q);plt.close(fig)


def housekeeping(root):
    """Move only the five requested regenerable T19 ledgers, with hash verification.

    A sandbox denial is a recorded incomplete housekeeping item, never a reason
    to delete the source or to try another spelling of the protected directory.
    """
    destination=root.parent/'t19'
    names=('order-exceptions','field-orders-outside','constraint-orders',
           'temporal-failures','native-constraint-norms')
    denial=''
    try:
        destination.mkdir(parents=True,exist_ok=True)
    except OSError as e:
        denial=str(e)
    receipts=[]
    for name in names:
        filename='t19-'+name+'.csv';src=HERE/filename;dst=destination/filename
        if src.exists():
            n=src.stat().st_size;sha=digest(src)
            assert n>1000000
            if denial:
                status='BLOCKED_filesystem_sources_retained'
            else:
                try:
                    if dst.exists():
                        assert digest(dst)==sha, 'Existing destination differs: '+str(dst)
                    else:
                        pending=dst.with_suffix('.csv.pending')
                        shutil.copyfile(src,pending)
                        assert digest(pending)==sha
                        pending.replace(dst)
                    assert digest(dst)==sha
                    src.unlink()
                    status='MOVED_SHA256_verified'
                except OSError as e:
                    denial=str(e);status='BLOCKED_filesystem_sources_retained'
        else:
            assert dst.exists(), 'Missing T19 ledger at source and destination: '+filename
            n=dst.stat().st_size;sha=digest(dst);status='ALREADY_MOVED_SHA256_recorded'
        receipts.append(dict(filename=filename,bytes=n,sha256=sha,destination=str(dst),
                             status=status,error=denial))
    output('housekeeping',receipts)
    # Keep current references valid; future generator reports use BULK_DIR paths.
    p=HERE/'t19-stageD-verdict.md'
    md=p.read_text().split('\n## T20 ledger housekeeping\n')[0]
    for r in receipts:
        if r['status'].startswith(('MOVED','ALREADY')):
            md=md.replace('`'+r['filename']+'`','['+r['filename']+']('+r['destination']+')')
    md=md.replace('COLLECTED_PRODUCTION_ROOT CONTRACT_DIR`',
                  'COLLECTED_PRODUCTION_ROOT CONTRACT_DIR [--output-dir DIR]`')
    md+='\n## T20 ledger housekeeping\n\n'
    md+=f'The five large ledgers default to `{destination}` in the generator; `--output-dir DIR` overrides that directory. '
    md+=('The move is blocked by the session filesystem sandbox; all five worktree sources and their current references are retained. '
         if denial else 'The five ledgers have been moved out of the worktree after SHA-256 verification. ')
    md+='The per-file sizes, hashes, destinations and status are in [t20-housekeeping.csv](t20-housekeeping.csv). The original T19 numerical reading and its resource receipt are unchanged.\n'
    p.write_text(md)
    return receipts,denial


def report(root,stats,history,speeds,packet,rh,rs,audits,coverage,receipts,denial):
    packet_rows=[]
    for ray,field in (('axis','Gamma1'),('equator','Gamma2'),('diagonal','Gamma1'),('diagonal','Gamma2')):
        g=[r for r in packet if r['ray']==ray and r['field']==field and r['pair']=='mid_high' and 1<=r['clock']<=6]
        t=np.array([r['time_M'] for r in g]);r=np.array([r['radius_M'] for r in g])
        v,offset=np.polyfit(t,r,1)
        residual=float(np.max(abs(r-offset-v*t)))
        packet_rows.append(dict(ray=ray,field=field,pair='mid_high',start_M=t[0],end_M=t[-1],
            crest_speed_M_per_M=v,intercept_radius_M=offset,max_crest_line_residual_M=residual,
            radii_M=';'.join(f'{x:.8g}' for x in r),
            note='crest/envelope tracking, not a phase-speed or causal-source certificate'))
    output('packet-speed',packet_rows)
    mechanisms=[
        dict(rank=1,mechanism='resolution-dependent phase/amplitude of a gauge-driven oscillatory launch',
             for_evidence='A11,K,lapse,phi share a damped ~5.25 M sampled oscillation; differences change sign; pair extrema occur at different clocks; outward Gamma packet speed ~1',
             against_or_unknown='No finer fixed-radius profiles; no distinct rung frequency is resolved; not a proof of dephasing or of pure gauge',
             cheapest_discriminator='First re-extract already retained coarse checkpoints at their own clocks, if operator permits future compute-node postprocessing. If a new run is required: one high rung to 10.5 M with native field/RHS capture at the original 13 clocks. Gauge/equations unchanged.'),
        dict(rank=2,mechanism='puncture regularity/source reaching the horizon',
             for_evidence='Late A11 error has broad support from numerical horizon outward; cache-edge peaks often point toward R<0.005; inner field response is largest',
             against_or_unknown='No collected profile below R=0.005; origin, regularity and causal path cannot be measured; numerical horizon drifts remain very small',
             cheapest_discriminator='One high-rung pilot through 10.5 M with one extra puncture level and all exterior faces unchanged; original masks/clocks/thresholds retained; add puncture RHS capture at those clocks.'),
        dict(rank=3,mechanism='refinement-face or box-seam errors/reflections',
             for_evidence='Small stationary parity-odd Gamma2 axis spikes occur near successive faces; some inner profiles have many radial sign reversals; box seams differ between rungs',
             against_or_unknown='Dominant late A11 crest/energy is broad and shared across rays rather than tied to their different face radii; 25 volume witnesses are smooth_interior; outward leading packet traverses faces',
             cheapest_discriminator='One high-rung pilot with a parametrically wider finest-level footprint, then if needed a box-size-only repeat. Same h, gauge, KO, transfers, masks and diagnostic clocks; compare whether the feature follows the moved face or seam.'),
        dict(rank=4,mechanism='coarsest rung outside the asymptotic regime',
             for_evidence='Largest raw spatial pair differences grow in the inner region; instantaneous pair orders swing widely',
             against_or_unknown='Low has about 113 horizon-radius cells and fine about 254; smooth exterior medians generally converge; both pair errors grow and cancel, so low alone is not identified',
             cheapest_discriminator='One extra globally finer rung h0=7/27, same max_level and common physical faces, through the same 10.5 M and original masks/clocks; use mid/high/new only as a separately declared causal pilot, not a replacement T19 verdict.'),
        dict(rank=5,mechanism='slow numerical instability',
             for_evidence='Inner differences are larger late than early',
             against_or_unknown='RMS falls repeatedly; fitted log slopes depend strongly on pair/window; native oscillation damps; no resolution-scaled positive eigen-growth rate is established',
             cheapest_discriminator='Independent longer high-rung continuation after the launch has settled, retaining the original registered samples exactly. Later-time growth observations would be outside the T19 registration, with no new acceptance threshold.'),
        dict(rank=6,mechanism='different or evolving refinement-region geometry',
             for_evidence='Internal box partitions differ from t=0',
             against_or_unknown='All level union faces identical; nine restart box lists match each census exactly; regrid disabled; all common-support counts/weights constant',
             cheapest_discriminator='No new geometry-mismatch run warranted. Box-seam dependence remains in rank 3, distinct from a changing refinement footprint.')]
    output('mechanisms',mechanisms)
    def metric(region,k,ray,field,pair):
        return next(r for r in stats[region] if (r['clock'],r['ray'],r['field'],r['pair'])==(k,ray,field,pair))
    horizon_table=[]
    for ray in RAYS:
        r=metric('horizon',11,ray,'A11','mid_high')
        horizon_table.append(f"| {ray} | {r['peak_radius_M']:.8g} | {r['energy_R10_M']:.8g}–{r['energy_R90_M']:.8g} | {r['radial_dr_RMS']:.7g} |")
    speed_table=[]
    for r in packet_rows:
        speed_table.append(f"| {r['ray']} {r['field']} | {r['crest_speed_M_per_M']:.5g} | {r['max_crest_line_residual_M']:.4g} |")
    selected=[s for s in speeds if s['ray']=='axis' and s['clock'] in (3,4,5,6)
              and abs(s['radius_M']-s['time_M'])<.1]
    selected=list({(s['clock'],s['radius_M']):s for s in selected}.values())
    vlight=[s['light_out'] for s in selected];vlapse=[s['lapse_subblock_out'] for s in selected]
    ps=[r['crest_speed_M_per_M'] for r in packet_rows]
    late1=metric('inner',11,'axis','Gamma1','low_mid')['peak_radius_M']
    late2=metric('inner',12,'axis','Gamma1','low_mid')['peak_radius_M']
    restarts=sum(r['interpretation']=='restart_exact_layout' for r in audits)
    md=f'''# T20 — spatial and temporal localization of exp-0023 differences

T19's registered **FAIL stands**. This analysis identifies an outward Γ̃-dominated transient and a later, broad near-horizon oscillatory spatial error. It establishes neither the puncture origin nor a unique phase-error or instability mechanism. The independently checked ledger has {len(rs)} failing witnesses, including {sum(int(r['clock_index'])==11 for r in rs)} at k=11 (9.625 M). The full clocks/masks/classes remain in [t20-failing-witness-localization.csv](t20-failing-witness-localization.csv). No static solution is read, explained against, or subtracted; only numerical rungs are compared.

## Spatial support

The [26-page signed ray atlas](t20-signed-ray-atlas.pdf) contains low−mid and mid−high for all ten requested fields, three rays and thirteen clocks. Black lines are each leg's **numerical** N96 horizon radius and gray lines are their common refinement faces. The [ten-page space–time atlas](t20-space-time.pdf) and [time traces](t20-time-traces.pdf) expose signed lobes and moving packets. Tables `t20-localization-*.csv` give dr-weighted RMS, I8/I10 sensitivity, energy quantile radii, sign changes and distance to a face, including the unchanged registered radial masks. These ray norms use the reduction's radial dr weights; they are not volume RMS. The four extra descriptive partitions are explicitly named near, inner, exterior_near and outer_cache and carry no acceptance verdict.

At 9.625 M, the A11 mid−high difference on the fixed horizon mask is broad:

| Ray | Peak R/M | 10–90% radial error-energy radii | radial RMS |
|---|---:|---:|---:|
{chr(10).join(horizon_table)}

The peak is at the mask's first sample, R≈0.006049, while the numerical high horizon is R≈{rh['E-high',11]:.8g}. It is a clipped peak, not a measured interior maximum. Its broad energy extends across the axis/equator finest face at 0.0068359375 M and the diagonal face at 0.0096674755 M. The common radial peak on all three rays is evidence against an error confined solely to a refinement face. A11, A22, Aww, K, lapse, χ, φ and Ex have broad inner lobes; Γ̃ has those lobes plus travelling exterior structure. Profiles outside the first few inner cells also have sign-reversing ripples. Small axis Γ̃2 features sit near successive faces and persist in time, distinct from the dominant broad A11 error. Components constrained by reflection parity and their I8/I10 spreads are retained rather than interpreted as a physical wave. There are 25 smooth-interior volume FAIL witnesses overall, 14 at k=11; face/corner/seam/axis effects alone are not established as the cause of every failure.

All caches start at R=0.005 M and stop at 20 M. Many inner-region peaks are clipped at that first radius. The collected data cannot resolve the puncture itself or determine whether the visible inner lobe was generated there.

## Hierarchy identity

[Faces](t20-hierarchy-faces.csv), [internal seams](t20-hierarchy-seams.csv) and [own-record checks](t20-hierarchy-audit.csv) distinguish refinement geometry from box partition. Each rung has levels 0–14, h0=7/8, 7/12, 7/18 M. Their union faces are identical: level 0 ±336 M, then 56, 28, 14, 7, 3.5, 1.75, 0.875, 0.4375, 0.21875, 0.109375, 0.0546875, 0.02734375, 0.013671875 and 0.0068359375 M. The initial census and {restarts} production checkpoint-read layouts (each rung at 2.625, 5.25 and 7.875 M) have exactly identical integer box lists within each leg. Checkpoint headers print all fifteen regrid intervals zero. Every common-support cell count, weight and status is identical at all thirteen clocks. Together these establish fixed, common refinement footprints throughout the run, rather than comparisons across different evolving level regions.

Internal partitions differ **from t=0**: refined-level boxes number 2/8/32. For example, the finest equatorial positive-x seams are absent on low, 0.00341796875 M on mid, and 0.001708984375/0.00341796875/0.005126953125 M on high. Thus the common-footprint h-refinement also changes same-level seams. E-low LB files were not collected; mid/high LB files contain rank loads, not coordinates. Qualified finder logs do not print checkpoint boxes, and bulk HDF5 layouts are absent locally. No unprinted per-clock layout is claimed as independently inspected; the all-clock statement follows from exact restart matches, disabled regridding and invariant common support.

## Time behaviour and available sampling

The controller's A11 horizon RMS history is reproduced in [mask histories](t20-mask-histories.csv): at 9.625 M D_low,mid=5.1131602e−6 and D_mid,high=2.3002972e−5, p=−3.7088414; at 10.5 M the differences are 3.0058147e−5 and 2.4772644e−5. Both pairs rose well above their early amplitudes, but repeatedly decreased too. The low−mid minimum and the mid−high peak are asynchronous. The temporal A11 difference remains ≤1.1891016e−10 in the table, so it cannot account for this reversal. All 344 failing self-difference witnesses retain their own passing temporal fractions; a blanket statement that every field/mask has clean temporal control would be false (T19 recorded 4,610 exceptions).

At fixed axis R=0.0074978675 M, the high-rung A11 sampled minima at 2.625 and 7.875 M, K maxima at 1.75 and 7 M, and lapse maxima at 0.875 and 6.125 M each give a **5.25 M extremum spacing**. Each sampled extremum is bracketed by its two neighbours; the resulting spacing bracket is **3.5–7 M**, conditional on one corresponding extremum in each bracket. No nearest-sample half-cadence assumption is used. The oscillation damps in the native fields. This is a sampled underlying-field timescale, not a certified difference-mode frequency. Signed pair differences at that radius cross zero multiple times and are not a single sinusoid. [Sampled extrema](t20-sampled-extrema.csv), [fixed-radius values](t20-fixed-radius.csv), [descriptive period fits](t20-history-fits.csv) and [growth fits](t20-growth-fits.csv) preserve their ambiguities. Positive-time log-RMS slopes depend on pair and early/late window and include decreasing intervals; no common, resolution-scaled exponential growth rate is established. Numerical t=0 roundoff is excluded from these growth fits rather than used to inflate a rate.

[Coverage](t20-history-coverage.csv) and [all collected finer histories](t20-finer-histories.csv) use the 49/73/109 numerical N96 horizon clocks (cadences 0.21875/0.14583333/0.09722222 M) and every native whole-grid constraint-norm clock. φ on the moving numerical horizon varies only by about 4.3e−8/2.0e−8/1.7e−8; small step-like changes prevent treating it as a clean fixed-radius phase trace. These records are not replaced by interpolated ray histories. No finer fixed-radius A, Γ̃, χ, K, lapse, φ or Ex profile series is present in the collected inputs, so finer phase/frequency and near-puncture causal estimates remain unavailable.

The late inner Γ̃ lobe also shifts outward: the axis low−mid crest in 0.013671875<R<1 moves from {late1:.6g} to {late2:.6g} M between 9.625 and 10.5 M, an apparent crest shift of {(late2-late1)/.875:.5g} M/M. Only two such endpoints and changing multi-lobe shapes are available here; this is not a measured characteristic speed or a link to the earlier outward packet.

## Exterior transient, 3.5–5.25 M

The Γ̃ packet is already visible at 0.875 M near R=0.8761 M, progresses to 1.7236 and 2.6162 M, then 3.4461, 4.3590 and 5.2309 M on the axis. Equatorial Γ̃2 follows the same radii; diagonal Γ̃1 and Γ̃2 follow 0.8761, 1.7306, 2.6268, 3.4601, 4.3767 and 5.2521 M to the cache accuracy. It traverses the **fixed** receiving-side mask [3.5,4.375] and far mask [4,8]; it is not stationary at a level boundary. Earliest observed radius/time is not a measured source location: t=0 ray differences are tiny, and the interval before 0.875 M and R<0.005 is not sampled.

| Tracked mid−high crest | Linear speed M/M, 0.875–5.25 M | max radius residual from fitted line M |
|---|---:|---:|
{chr(10).join(speed_table)}

The observed outward crest speeds span {min(ps):.5g}–{max(ps):.5g}. At the axis packet the **numerically evaluated** outward coordinate light speeds span {min(vlight):.4g}–{max(vlight):.4g}, frozen lapse-subblock speeds {min(vlapse):.4g}–{max(vlapse):.4g}, and longitudinal Γ-driver estimates are approximately 1. The evidence supports a Γ-driver-dominated gauge transient rather than a light-speed matter packet. The radial Γ component carries the clearest coherent crest; A/K/lapse/χ carry weaker broad features and ripples, while φ/Ex do not identify an independent crest at this speed. Crest selection can switch lobes: [the packet ledger](t20-leading-packet.csv) marks corridor-edge maxima, and no speed is assigned to those field rows. The [exterior figure](t20-exterior-packet.png) shows the signed curves, not just their peak locations.

Code basis: `Source/CCZ4/ExperimentalGauge.hpp:84–91` uses lapse coefficient 1.8 and shift Γ coefficient 0.75; `Source/Cartoon/CCZ4Cartoon.impl.hpp:400–405` supplies the shift Laplacian and longitudinal 1/3 term. The local ray-normal estimates are −β_n+α√(χ h^nn), −β_n+√(1.8 α χ h^nn), and −β_n+√((4/3)F h^nn), with F=0.75. [Characteristic estimates](t20-characteristic-estimates.csv) use the numerical high-rung state only. [The exact CAS witness](t20-cas-witness.json) proves the 2×2 characteristic polynomial and checks independent numerical eigenvalues; the **full coupled gauge spectrum is not proved**. These are frozen-principal-subblock comparisons, not a full-system characteristic certificate.

Volume Γ̃1 on far has p=4.0291 at 3.5 M, 0.14150 at 4.375 M, 1.03850 at 5.25 M and 4.01891 at 9.625 M. The dip accompanies the packet entering that mask. At 4.375 M the axis low−mid crest is at 4.4663144 M (signed 1.4942707e−9), while mid−high is at 4.3590438 M (signed −1.4369343e−9), a 0.1072707 M separation. The latter lies inside receiving_side_3p5 and the former outside its upper face. These are directly observed differences of spatial lobe positions/signs; native rung frequency error is not thereby proved. Cancellation and which part of the packet lies on the fixed mask change the norm ratio. The seven direct-to-zero constraint FAIL witnesses at 5.25 M on far/receiving_side_3p5 remain separately recorded in T19. The profiles locate a spatial transient but do not prove that all seven constraint witnesses have the same causal origin.

## Mechanism ranking and discriminating work

[t20-mechanisms.csv](t20-mechanisms.csv) ranks the six candidates and gives a cheapest separate test for each still standing. Oscillatory launch response with resolution-dependent phase/amplitude is strongest; puncture/source sensitivity and face/seam sensitivity remain plausible; coarsest-only under-resolution and a slow instability have weaker support. Different moving refinement footprints are contradicted by the records. A new diagnostic pilot must preserve the code's gauge, KO, equations, transfers and initial-data-only static use, and must leave the original registered masks, clocks and thresholds unchanged. None was run, and none can revise T19's FAIL.

## T19 bulk housekeeping and reproduction

The updated [T19 generator](t19-stageD-verdict.py) accepts `--output-dir DIR` and defaults the five >1 MB ledgers to `{root.parent/'t19'}`. [The housekeeping receipt](t20-housekeeping.csv) records every source/destination, byte size and SHA-256. {'Creation of that destination is denied by the session sandbox; the actual move is BLOCKED and all five original ledgers are retained with valid references. No protected-path alias or permission bypass was attempted.' if denial else 'All five ledgers were moved and SHA-256 verified before removal of worktree sources.'} T19 report/README references and the manifest record this state; the numerical verdict and old numerical resource receipt are unchanged.

Regenerate with `python Tests/EMSNative/t20-difference-localization.py COLLECTED_PRODUCTION_ROOT`. Missing required inputs or changed census boxes stop with a named error. The [resource receipt](t20-resources.json) and [manifest](t20-manifest.txt) record measured peak RSS, single-thread execution and every used input/output hash. No bulk is copied into the worktree, no new evolution or horizon solve is started, and no commit or remote operation is performed.
'''
    text_output('difference-localization.md',md)
    p=HERE/'README.md'
    old=p.read_text()
    marker='\n## T20 — exp-0023 difference localization\n'
    section=f'''{marker}

T19's registered **FAIL is unchanged**. [The T20 reading](t20-difference-localization.md) locates an outward Γ̃ packet with crest speed about 1 M/M and a later broad near-horizon spatial error. The [signed ray atlas](t20-signed-ray-atlas.pdf) covers all ten requested fields, axis/equator/diagonal and all thirteen clocks; [space–time maps](t20-space-time.pdf), [fixed-radius time traces](t20-time-traces.pdf), [exterior packet](t20-exterior-packet.png), numerical horizon radii and the per-mask `t20-localization-*.csv` tables preserve signs and sampling limitations. The 351 T19 witnesses include 163 at 9.625 M. All nine production restart layouts match their own census boxes; physical refinement faces coincide across rungs and regridding is disabled, while box seams differ from t=0. [Hierarchy evidence](t20-hierarchy-audit.csv) names the missing E-low LB files and the scope of the inference.

Near R=0.0075 M the native A/K/lapse oscillation has sampled extremum spacing 5.25 M (conditional neighbouring-sample spacing bracket 3.5–7 M), but pair differences are not a single sinusoid and no common exponential growth rate is established. All collected finer horizon/whole-grid scalar histories are used; fixed-radius ray histories exist only at the thirteen diagnostic clocks and no profile exists below R=0.005 M. [Mechanisms and separate discriminating pilots](t20-mechanisms.csv) remain hypotheses; no masks, clocks, uncertainty rules or thresholds are changed, and no pilot is run.

The five T19 large-ledger destinations now default to the production sibling `t19/` directory with an explicit `--output-dir DIR` override. {'The actual move is blocked by the filesystem sandbox, so all five sources and valid local references remain; ' if denial else 'The five ledgers were moved with SHA-256 verification; '}[the receipt](t20-housekeeping.csv) names every file and its status. Regenerate T20 with `python Tests/EMSNative/t20-difference-localization.py <collected production root>`. [Measured resources](t20-resources.json) and [SHA-256 manifest](t20-manifest.txt) accompany the outputs. Numerical rungs only; no static-data subtraction, evolution, horizon solve, SSH or commit.
'''
    if marker in old:
        prefix,tail=old.split(marker,1)
        rest='\n## '+tail.split('\n## ',1)[1] if '\n## ' in tail else ''
        old=prefix+section+rest
    else:old+=section
    note='\nT20 housekeeping: the five large T19 ledgers now default to the production sibling `t19/` directory (`--output-dir DIR` overrides). '
    note+=('The sandbox blocks the actual move; current local references and all original ledgers remain valid. ' if denial else 'The actual move is SHA-256 verified. ')
    note+='See [the per-file receipt](t20-housekeeping.csv).\n'
    start=old.index('\n## T19 — registered Stage D reading of exp-0023\n')
    end=old.find('\n## ',start+4)
    if end<0:end=len(old)
    prior=old[start:end]
    prior=re.sub(r'\nT20 housekeeping:.*?\n','\n',prior)
    for r in receipts:
        if r['status'].startswith(('MOVED','ALREADY')):
            prior=prior.replace(']('+r['filename']+')',']('+r['destination']+')')
    old=old[:start]+prior+note+old[end:]
    p.write_text(old)
    OUTPUTS.add(p)
    # Re-seal changed T19 path/report/script outputs without rereading its numerical analysis.
    p=HERE/'t19-manifest.txt'
    lines=p.read_text().splitlines()
    new=[];output_section=False
    moved={str(HERE/r['filename']):r['destination'] for r in receipts if r['status'].startswith(('MOVED','ALREADY'))}
    for line in lines:
        if line.startswith('OUTPUT_SHA256'):output_section=True
        m=re.match(r'^[0-9a-f]{64}  (.*)$',line)
        if output_section and m:
            target=Path(moved.get(m[1],m[1]))
            assert target.exists(), 'Missing previously sealed output: '+str(target)
            line=digest(target)+'  '+str(target)
        if not line.startswith('T20 housekeeping status:'):new.append(line)
    new.append('T20 housekeeping status: '+('BLOCKED filesystem; sources retained' if denial else 'MOVED and SHA256 verified')+'; numerical T19 reading unchanged')
    p.write_text('\n'.join(new)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production_root',type=Path)
    args=parser.parse_args()
    root=args.production_root.absolute()
    start=time.monotonic()
    for p in (Path('Source/CCZ4/ExperimentalGauge.hpp'),Path('Source/Cartoon/CCZ4Cartoon.impl.hpp'),
              Path('Examples/EMS/EMSBH2DLevel.cpp'),Path('Examples/EMS/DiagnosticVariables.hpp')):
        INPUTS.add(HERE.parents[1]/p)
    faces,audits=layouts(root)
    print('T20: census/restart layouts and constant common support checked',flush=True)
    cache,names,radii,stats,fixed,history,speeds,packet=cached_analysis(root)
    finer,coverage=additional_histories(root)
    rh=horizon_radii(root)
    rs,maskhist=witnesses(root,stats)
    print('T20: all 39 caches, numerical horizons and T19 witnesses checked; rendering',flush=True)
    figures(cache,names,radii,rh,finer)
    receipts,denial=housekeeping(root)
    report(root,stats,history,speeds,packet,rh,rs,audits,coverage,receipts,denial)
    result=dict(task='T20 localization; registered T19 FAIL unchanged',root=str(root),
        cache_count=len(cache),ray_clocks=13,ray_radius_min_M=float(radii[0]),ray_radius_max_M=float(radii[-1]),
        witnesses=len(rs),witnesses_at_k11=sum(int(r['clock_index'])==11 for r in rs),
        restart_layouts_verified=sum(r['interpretation']=='restart_exact_layout' for r in audits),
        elapsed_seconds=time.monotonic()-start,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        process_cap_bytes=6000000000,threads=1,
        missing='profiles R<0.005 M; fixed-radius fields between the 13 ray clocks; checkpoint HDF5 bulk; E-low LB',
        no_static_data=True,no_evolution=True,no_SSH=True,no_commit=True)
    result['housekeeping_status']='BLOCKED_filesystem_sources_retained' if denial else 'MOVED_verified'
    # macOS reports ru_maxrss in bytes; Linux reports KiB.
    if sys.platform!='darwin':result['peak_RSS_bytes']*=1024
    assert result['peak_RSS_bytes']<result['process_cap_bytes']
    manifest=['T20: numerical rungs only; T19 FAIL not reopened; no commit',
              f"Measured peak RSS: {result['peak_RSS_bytes']} bytes; cap 6000000000; threads 1",
              'INPUT_SHA256']
    for p in sorted(INPUTS):manifest.append(digest(p)+'  '+str(p))
    result['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
    result['elapsed_seconds']=time.monotonic()-start
    assert result['peak_RSS_bytes']<result['process_cap_bytes']
    manifest[1]=f"Measured peak RSS: {result['peak_RSS_bytes']} bytes; cap 6000000000; threads 1"
    text_output('resources.json',json.dumps(result,indent=2)+'\n')
    manifest.append('OUTPUT_SHA256 (manifest excluded)')
    for p in sorted(OUTPUTS|{Path(__file__).resolve()}):manifest.append(digest(p)+'  '+str(p))
    text_output('manifest.txt','\n'.join(manifest)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    main()
