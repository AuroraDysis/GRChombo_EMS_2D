"""Read unchanged RHFinder output and matching Chombo plots; call Julia for mass.

uv run --with h5py --with numpy python postprocess.py RUN_ROOT OUTPUT_DIR EMS_ROOT [FRESH_JSON]
Only the last row per search/time can be accepted, and it must say `found`.
The small t=0 test hierarchy is one uniform level with reflection at y=0.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np


def params(path):
    return {k.strip(): v.split('#')[0].strip()
            for line in path.read_text().splitlines() if '=' in line and not line.startswith('#')
            for k, v in [line.split('=', 1)]}


def plot(path):
    with h5py.File(path) as f:
        levels = [key for key in f if key.startswith('level_')]
        if levels != ['level_0']:
            raise ValueError('this isolated test reader requires one uniform level')
        names = [f.attrs[f'component_{i}'].decode() for i in range(int(f.attrs['num_components']))]
        g = f['level_0']
        dx, time = float(g.attrs['dx']), float(g.attrs['time'])
        domain = g.attrs['prob_domain']
        nx, ny = int(domain['hi_i'])+1, int(domain['hi_j'])+1
        a = np.full((len(names), ny, nx), np.nan)
        flat, offsets = g['data:datatype=0'][:], g['data:offsets=0'][:]
        for box, lo, hi in zip(g['boxes'][:], offsets[:-1], offsets[1:]):
            x0, y0, x1, y1 = (int(box[k]) for k in ('lo_i','lo_j','hi_i','hi_j'))
            x1 += 1; y1 += 1
            if hi-lo != len(names)*(x1-x0)*(y1-y0):
                raise ValueError('plot must omit ghost cells')
            a[:,y0:y1,x0:x1] = flat[lo:hi].reshape(len(names),y1-y0,x1-x0)
        if not np.isfinite(a).all():
            raise ValueError('nonfinite or uncovered plotted cells')
        return dict(zip(names,a)), dx, time


def interp(a, dx, x, y, odd_y=False):
    """The same four-point, degree-three stencil as Lagrange<4>, including reflection."""
    nodes, weights = [], []
    for z in (x, y):
        u = z/dx-0.5
        indices = np.floor(u).astype(int)[:,None]-1+np.arange(4)
        w = np.ones(indices.shape)
        for k in range(4):
            for j in range(4):
                if j != k:
                    w[:,k] *= (u-indices[:,j])/(indices[:,k]-indices[:,j])
        nodes.append(indices); weights.append(w)
    xi, yi = nodes
    sign = np.where(yi < 0, -1 if odd_y else 1, 1)
    yi = np.where(yi < 0, -yi-1, yi)
    if xi.min()<0 or xi.max()>=a.shape[1] or yi.max()>=a.shape[0]:
        raise ValueError('surface stencil reaches unsupported outer boundary')
    return np.einsum('pi,pj,pij->p', weights[1]*sign, weights[0], a[yi[:,:,None],xi[:,None,:]])


def scalar_surface(fields, dx, centre, shape):
    n = len(shape); dt = math.pi/n
    th = (np.arange(n)+0.5)*dt
    f = np.pad(shape,2,mode='symmetric')
    fp = (-f[4:]+8*f[3:-1]-8*f[1:-3]+f[:-4])/(12*dt)
    x, y = centre+shape*np.cos(th), shape*np.sin(th)
    v = {key:interp(fields[key],dx,x,y,key in ('h12','Ey'))
         for key in ('chi','h11','h12','h22','hww','phi','Ex','Ey')}
    tx, ty = fp*np.cos(th)-shape*np.sin(th), fp*np.sin(th)+shape*np.cos(th)
    htt = v['h11']*tx*tx+2*v['h12']*tx*ty+v['h22']*ty*ty
    dA = 2*math.pi*shape*np.sin(th)/v['chi']*np.sqrt(htt*v['hww'])*dt
    if not np.isfinite(dA).all() or (dA <= 0).any():
        raise ValueError('invalid induced surface area')
    area = float(dA.sum())
    mean = float(np.dot(dA,v['phi'])/area)
    rms = math.sqrt(float(np.dot(dA,(v['phi']-mean)**2)/area))
    return area, mean, rms


def coincident(a, b, tolerance=0.002):
    def enclosed(outer, inner):
        shape, centre = inner
        th = (np.arange(len(shape))+0.5)*math.pi/len(shape)
        x, y = centre-outer[1]+shape*np.cos(th), shape*np.sin(th)
        oth = (np.arange(len(outer[0]))+0.5)*math.pi/len(outer[0])
        radius = np.interp(np.arctan2(y,x),oth,outer[0])
        return bool((np.hypot(x,y) <= (1+tolerance)*radius).all())
    return enclosed(a,b) and enclosed(b,a)


def self_check():
    dx = 1/32
    y, x = np.mgrid[0:64,0:128]
    x=(x+0.5)*dx; y=(y+0.5)*dx
    xx=np.array([1.1,1.3,1.7]); yy=np.array([0.001,0.2,0.8])
    assert np.max(abs(interp(2+x*x+y*y,dx,xx,yy)-(2+xx*xx+yy*yy))) < 1e-13
    assert np.max(abs(interp(x*y,dx,xx,yy,True)-xx*yy)) < 1e-13
    for n in (48,96,192):
        fields={key:np.ones_like(x) for key in ('chi','h11','h22','hww','phi')}
        fields.update(h12=np.zeros_like(x), Ex=np.zeros_like(x), Ey=np.zeros_like(x))
        area,mean,rms=scalar_surface(fields,dx,2.,np.full(n,0.6))
        bias=(math.pi/(2*n))/math.sin(math.pi/(2*n))
        assert abs(area/(4*math.pi*0.6**2)-bias)<1e-13 and abs(mean-1)<1e-13 and rms<1e-13
    assert coincident((np.ones(48),0.),(np.ones(96),0.))
    assert not coincident((np.ones(48),0.),(np.ones(48),4.))
    assert not coincident((np.ones(48),0.),(np.full(48,0.9),0.))


def process(root, output, ems, fresh_file=None):
    output.mkdir(parents=True,exist_ok=True)
    rows, manifest = [], []
    fresh={r['source']:r for r in json.loads(fresh_file.read_text())} if fresh_file else {}
    for directory in sorted(root.glob('*-n*-dx*')):
        status=json.loads((directory/'status.json').read_text())
        if status['exit'] != 0:
            raise RuntimeError(f'failed find: {directory}: {status}')
        p=params(directory/'params.txt')
        if p['max_level']!='0' or p.get('star_centre')!='2 0':
            raise ValueError('unsupported isolated hierarchy/centre')
        plots=[(file,plot(file)) for file in (directory/'plt').glob('*.hdf5')]
        candidates=[]
        for surface in sorted(directory.glob('rh_surf_*.dat')):
            index=int(surface.stem.split('_')[-1])
            lines=[line.split() for line in surface.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
            shape_file=directory/f'rh_f{index}.dat'
            shapes=[list(map(float,line.split())) for line in shape_file.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
            if len(lines)!=len(shapes):
                raise ValueError('shape and diagnostic row counts differ')
            latest={float(row[0]):i for i,row in enumerate(lines)}
            for time,i in sorted(latest.items()):
                row=lines[i]
                if row[-1]!='found': continue
                if shapes[i][0]!=time: raise ValueError('shape time mismatch')
                match=[(f,data) for f,data in plots if abs(data[2]-time)<1e-12]
                if len(match)!=1: raise ValueError('missing/ambiguous matching plot')
                plot_file,(fields,dx,_) = match[0]
                shape=np.array(shapes[i][1:]); centre=float(row[2])
                area,phi,rms=scalar_surface(fields,dx,centre,shape)
                A,Q=float(row[4]),float(row[17])
                # Serialization is nine significant digits in both shape and native diagnostics.
                if abs(area/A-1)>1e-7:
                    raise ValueError(f'plot/shape quadrature does not reproduce native A: {area}, {A}')
                r=dict(case=status['case'],N_theta=len(shape),resolution=status['resolution'],time=time,
                       horizon_id=0,search_index=index,duplicate_count=1,A=A,Q=Q,
                       phi_mean=phi,phi_rms=rms,area_from_plot=area,area_readback_relative=area/A-1,
                       expansion_squared=float(row[9]),theta_minus=float(row[8]),M_RN_legacy=float(row[18]),
                       spread_A=0.,spread_Q=0.,ems_alpha=float(p['ems_alpha']),f0=float(p['ems_f0']),
                       f1=float(p['ems_f1']),f2=float(p['ems_f2']),phi_inf=0.,
                       branch='rn' if status['case']=='rn' else 'scalarized-positive-n0',
                       threshold=status['threshold'],seconds=status['seconds'])
                r.update(expansion_squared_native=r['expansion_squared'],expansion_status='NATIVE_UNVERIFIED')
                if fresh_file:
                    check=fresh[str(directory)]
                    if index!=0 or check['exit'] not in (0,1) or check['shape_sha256']!=hashlib.sha256(shape_file.read_bytes()).hexdigest():
                        raise ValueError('fresh check missing or saved geometry changed')
                    if not all(math.isfinite(check[k]) for k in ('A','Q','theta_minus','expansion_squared')):
                        raise ValueError('nonfinite fresh expansion measurement')
                    if abs(check['A']/A-1)>1e-7 or abs(check['Q']/Q-1)>1e-7:
                        raise ValueError('fresh geometry measurement differs from native output')
                    r.update(expansion_squared=check['expansion_squared'],theta_minus=check['theta_minus'],
                             expansion_status='FRESH_REPLAY' if check['exit']==0 else 'FRESH_REPLAY_ABOVE_THRESHOLD')
                candidates.append((r,(shape,centre)))
                for file in (surface,shape_file,plot_file,directory/'params.txt'):
                    manifest.append(dict(path=str(file),sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
                # Compact raw evidence retained in the fork; complete plots stay in the run directory.
                evidence=output/'accepted'/directory.name
                evidence.mkdir(parents=True,exist_ok=True)
                (evidence/surface.name).write_text(' '.join(row)+'\n')
                (evidence/shape_file.name).write_text(' '.join(map(str,shapes[i]))+'\n')
                (evidence/'params.txt').write_text((directory/'params.txt').read_text())
        unique=[]
        for row,shape in candidates:
            group=next(((r,s,rs) for r,s,rs in unique if r['time']==row['time'] and coincident(s,shape)),None)
            if group:
                r=group[0];r['duplicate_count']+=1
                group[2].append(row)
                r['spread_A']=max(v['A'] for v in group[2])-min(v['A'] for v in group[2])
                r['spread_Q']=max(v['Q'] for v in group[2])-min(v['Q'] for v in group[2])
            else:
                row['horizon_id']=sum(r['time']==row['time'] for r,_,_ in unique)
                unique.append((row,shape,[row]))
        if not unique: raise RuntimeError(f'no found rows: {directory}')
        for r,_,members in unique:
            r['duplicate_A_Q']=';'.join(f"{m['A']:.17g}:{m['Q']:.17g}" for m in members)
            rows.append(r)
    if not rows: raise RuntimeError('no runs')
    table=output/'measurements.csv'
    with table.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    (output/'inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')
    script=Path(__file__).with_name('postprocess.jl')
    subprocess.run(['julia','--project=test',str(script),str(table.resolve()),str((output/'gates.csv').resolve())],cwd=ems,check=True)
    print(f'Postprocessed {len(rows)} physical horizon records into {output}',flush=True)


if __name__=='__main__':
    self_check()
    process(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),Path(sys.argv[3]).resolve(),
            Path(sys.argv[4]).resolve() if len(sys.argv)>4 else None)
