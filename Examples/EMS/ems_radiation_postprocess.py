#!/usr/bin/env python3
"""Native EMS radiation: Float64, no filters, fixed-retarded-time extrapolation.

Input: ems_radiation.csv from Main_EMSBH2DBH. NumPy is the only dependency.
Energies start at the first retained sample, never imply missing early flux=0.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def cumulative(y, x):
    y = np.asarray(y)
    dx = np.diff(x).reshape((-1,) + (1,) * (y.ndim - 1))
    return np.concatenate([np.zeros_like(y[:1]),
                           np.cumsum((y[1:] + y[:-1]) * dx / 2, axis=0)])


def load(path):
    records = {}
    block=[]
    previous_time=-np.inf
    with Path(path).open() as stream:
        for row in csv.DictReader(stream):
            row = {k: float(v) for k, v in row.items()}
            if not all(np.isfinite(v) for v in row.values()):
                raise ValueError("nonfinite radiation record")
            if row['m'] != 0 or row['l'] not in range(7):
                raise ValueError("expected axisymmetric l=0..6 records")
            if row['l'] != len(block) or (block and
                    (row['r'],row['t']) != (block[0]['r'],block[0]['t'])):
                raise ValueError("incomplete or unordered seven-mode block")
            if not block and row['t'] < previous_time:
                # An appended restart supersedes the entire old future, even
                # when the resumed run stops earlier than its predecessor.
                records={key:value for key,value in records.items() if key[1]<row['t']}
            block.append(row)
            previous_time=row['t']
            if len(block)==7:
                for value in block:
                    records[value['r'],value['t'],int(value['l'])]=value
                block=[]
    if block:
        raise ValueError("truncated final seven-mode block")
    result = {}
    for radius in sorted({key[0] for key in records}):
        times = sorted({key[1] for key in records if key[0] == radius})
        if len(times) < 3:
            raise ValueError("need at least three times at each radius")
        rows = [[records[radius, t, l] for l in range(7)] for t in times]
        data = {k: np.array([row[0][k] for row in rows])
                for k in ('R', 'alpha2', 'F_inf', 'L_scalar_stress', 'L_EM_stress')}
        for row in rows:
            for k in data:
                if any(v[k] != row[0][k] for v in row):
                    raise ValueError("inconsistent geometry across modes (incomplete restart row)")
        for name in ('S', 'D', 'P', 'W'):
            data[name] = np.array([[v[name+'_re'] + 1j*v[name+'_im']
                                    for v in row] for row in rows])
        data['t'] = np.array(times)
        result[radius] = data
    if not result:
        raise ValueError("empty radiation input")
    return result


def prepare(data, mass, coordinate_time=False, initial_news=None):
    """Boyle-style lapse-corrected clock; mass is an independently supplied ADM.

    alpha2 is the proper-area average of -1/g^tt. This is an asymptotic clock
    construction, not an exact gauge correction of the finite-radius flux.
    """
    t, radius = data['t'], data['R']
    if mass < 0 or not np.isfinite(mass) or np.any(radius <= 2*mass):
        raise ValueError("need finite M>=0 and areal radii R>2M")
    if np.any(data['alpha2'] <= 0) or np.any(data['F_inf'] <= 0):
        raise ValueError("nonpositive lapse squared or asymptotic coupling")
    if not np.all(data['F_inf'] == data['F_inf'][0]):
        raise ValueError("F_inf changed during run")
    rate = np.ones_like(t) if coordinate_time else np.sqrt(data['alpha2']/(1-2*mass/radius))
    clock = t[0] + cumulative(rate, t)
    rstar = radius if mass == 0 else radius + 2*mass*np.log(radius/(2*mass)-1)
    u = clock-rstar
    if np.any(np.diff(u) <= 0):
        raise ValueError("retarded time is not monotone")
    udot = np.gradient(u, t, edge_order=2)
    if np.any(udot <= 0):
        raise ValueError("nonpositive retarded-time derivative")
    # S=R*(phi-phi_inf); D is R*partial_t(phi), so include measured Rdot.
    ds = (data['D'] + (np.gradient(radius,t,edge_order=2)/radius)[:,None]*data['S'])/udot[:,None]
    news = cumulative(data['W'],u)
    if initial_news is not None:
        initial_news=np.asarray(initial_news)
        if initial_news.shape!=(7,) or not np.all(np.isfinite(initial_news)):
            raise ValueError("initial news must contain seven finite complex constants")
        news += initial_news
    return dict(data, u=u, ds=ds, news=news,
                dS_fd=np.gradient(data['S'],u,axis=0,edge_order=2))


def luminosities(news, ds, p, finf):
    return np.column_stack((np.sum(abs(news[:,2:])**2,axis=1)/(16*np.pi),
                            2*np.sum(abs(ds)**2,axis=1),
                            2*finf*np.sum(abs(p[:,1:])**2,axis=1)))


def save(path, u, luminosity, extra=None):
    energy = cumulative(luminosity,u)
    columns = [u, *luminosity.T, *energy.T]
    names = 'u,L_GW,L_phi,L_EM,E_GW,E_phi,E_EM'
    if extra:
        names += ',' + ','.join(extra)
        columns.extend(extra.values())
    np.savetxt(path,np.column_stack(columns),delimiter=',',header=names,comments='')
    return energy[-1].tolist()


def extrapolate(series, order, step):
    if len(series) < order+1:
        raise ValueError(f"order {order} needs at least {order+1} radii")
    lo=max(d['u'][0] for d in series); hi=min(d['u'][-1] for d in series)
    if hi <= lo or not np.isfinite(step) or step <= 0:
        raise ValueError("no common retarded-time window or invalid --du")
    u=np.arange(lo,hi+step*0.01,step)
    u=u[u<=hi]
    if len(u)<3:
        raise ValueError("common retarded-time window too short")
    R=np.array([np.interp(u,d['u'],d['R']) for d in series])
    fits={}
    for key in ('news','ds','P'):
        values=np.array([[np.interp(u,d['u'],d[key][:,l]) for l in range(7)]
                         for d in series]).transpose(2,0,1)
        fits[key]=np.array([np.polynomial.polynomial.polyfit(1/R[:,i],v,order)[0]
                            for i,v in enumerate(values)])
    return u,luminosities(fits['news'],fits['ds'],fits['P'],series[0]['F_inf'][0])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--mass',type=float,required=True,
                        help='independent ADM mass; never M_ADM.dat or MQ mass')
    parser.add_argument('--coordinate-time',action='store_true',help='gauge comparison: t-r*')
    parser.add_argument('--radii',type=float,nargs='+',help='radius subset, coordinate radii')
    parser.add_argument('--du',type=float,help='common retarded-time sampling interval')
    parser.add_argument('--initial-news',type=Path,
                        help='JSON mapping coordinate radius to 7 [real,imag] constants')
    args=parser.parse_args()
    raw=load(args.input)
    initials=json.loads(args.initial_news.read_text()) if args.initial_news else {}
    radii=args.radii or list(raw)
    if len(set(radii))!=len(radii):
        raise ValueError('duplicate radius subset')
    series={r:prepare(raw[r],args.mass,args.coordinate_time,
                     [complex(*v) for v in initials[str(r)]] if initials else None) for r in radii}
    if len({d['F_inf'][0] for d in series.values()})!=1:
        raise ValueError('inconsistent F_inf across radii')
    args.output.mkdir(parents=True,exist_ok=True)
    summary={'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
             'mass':args.mass,'clock':'coordinate' if args.coordinate_time else 'lapse-corrected',
             'initial_news':'supplied' if initials else 'assumed zero at first sample',
             'energy_order':['GW','phi','EM'],'radii':{},'extrapolation':{},
             'limitations':'No filtering. Per-radius energies start at each first sample. '
             'Extrapolated energies cover only the common u interval; early radiation '
             'outside the spheres, omitted modes, integration constants and late tails '
             'need independent bounds. Three radii give no quadratic fit redundancy.'}
    for r,d in series.items():
        lum=luminosities(d['news'],d['ds'],d['P'],d['F_inf'][0])
        extra={'t':d['t'],'R':d['R'],'L_scalar_stress':d['L_scalar_stress'],
               'L_EM_stress':d['L_EM_stress'],
               'L_matter_stress':d['L_scalar_stress']+d['L_EM_stress'],
               'E_matter_stress':cumulative(d['L_scalar_stress']+d['L_EM_stress'],d['t']),
               'L_phi_S_derivative':2*np.sum(abs(d['dS_fd'])**2,axis=1)}
        summary['radii'][str(r)]={'u_window':[d['u'][0],d['u'][-1]],
            'energy':save(args.output/f'r{r:g}.csv',d['u'],lum,extra)}
    step=args.du or max(np.median(np.diff(d['u'])) for d in series.values())
    for order in (1,2):
        if len(series)<order+1:
            summary['extrapolation'][str(order)]={'status':'insufficient radii'}
            continue
        u,lum=extrapolate(list(series.values()),order,step)
        summary['extrapolation'][str(order)]={'u_window':[u[0],u[-1]],
            'energy':save(args.output/f'infinity_order{order}.csv',u,lum)}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
