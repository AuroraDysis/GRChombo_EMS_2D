#!/usr/bin/env python3
"""Offline E t=0 census and native Float64 sensitivity audit; no evolution reads."""
import csv
import importlib.util
import math
import os
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t2', HERE / 't2-localize.py')
t2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t2)
SCALES = ('low', 'mid', 'high')
TARGETS = ('box_seam', 'patch_edge', 'convex_corner')
MASKS = {'horizon': (t2.RH['E'], 2*t2.RH['E']),
         'inside_inner_ring': (2*t2.RH['E'], t2.RINGS['E'][0]),
         'cavity': (t2.RINGS['E'][0], t2.RINGS['E'][1]),
         'between_rings': (t2.RINGS['E'][0], t2.RINGS['E'][2]),
         'outer_ring': (.8*t2.RINGS['E'][2], 1.2*t2.RINGS['E'][2]),
         'far': (4., 8.), 'cavity_core': (.055, .075),
         'outer_ring_core': (.85, 1.), 'far_core': (6.5, 8.)}


def save(name, rows):
    with (HERE / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def boxes(g):
    return np.array([[int(b[k]) for k in ('lo_i','lo_j','hi_i','hi_j')]
                     for b in g['boxes'][:]])


def rms(values, weights):
    return math.sqrt(float(np.sum(weights*values**2)/np.sum(weights)))


def prolong8(ix, iy, coarse, origin):
    """Operation-localization probe: eight current coarse chi samples per axis."""
    t=[ix/2.-.25,iy/2.-.25]
    starts=[np.floor(q).astype(int)-3 for q in t]
    weights=[]
    for q,first in zip(t,starts):
        ww=[]
        for a in range(8):
            w=np.ones(q.shape)
            for b in range(8):
                if b!=a: w *= (q-first-b)/(a-b)
            ww.append(w)
        weights.append(ww)
    def sample(x,y):
        y=np.where(y<0,-y-1,y)  # current chi's declared reflective parity
        return coarse[y-origin[1],x-origin[0]]
    base=sample(starts[0]+3,starts[1]+3); correction=np.zeros(base.shape)
    for a in range(8):
        for b in range(8):
            correction += weights[0][a]*weights[1][b]*(sample(starts[0]+a,starts[1]+b)-base)
    return base+correction


def build_audit(root):
    source=HERE/'T4bRoundoff.cpp'; binary=root/'t4b/roundoff.ex'
    chombo=Path(os.environ.get('CHOMBO_HOME','/Users/auroradysis/Workspace/EMS-deps/Chombo/lib'))
    repo=HERE.parents[1]
    command=[os.environ.get('CXX','/opt/homebrew/bin/g++-16'),'-O3','-std=c++17','-fopenmp',
             '-DCH_SPACEDIM=2','-DCH_Darwin','-DCH_LANG_CC','-DNDEBUG','-DCH_USE_64','-DCH_USE_DOUBLE',
             '-I'+str(repo/'Examples/EMS')]
    command+=['-I'+str(repo/'Source'/d) for d in ('utils','simd','CCZ4','Matter','Cartoon','BoxUtils','GRChomboCore')]
    command+=['-I'+str(chombo/'src'/d) for d in ('AMRTimeDependent','AMRTools','BoxTools','BaseTools')]
    command+=[str(source),'-L'+str(chombo),
              '-lboxtools2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC',
              '-lbasetools2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC','-o',str(binary)]
    with (root/'t4b/build.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)


def main(root):
    build_audit(root)
    failed = [r for r in csv.DictReader((HERE/'t4-t0-orders.csv').open())
              if r['member']=='E' and r['status']=='FAIL' and r['cell_class'] in TARGETS]
    needed = set((r['mask'],r['cell_class']) for r in failed)
    raw, geometry, identity, localization = [], [], [], []
    for scale in SCALES:
        oldpath = root/'t0'/f'E-{scale}'/'plt/EMS_Plot_000000.2d.hdf5'
        newpath = root/'t4b'/f'E-{scale}'/'plt/EMS_Plot_000000.2d.hdf5'
        marker = newpath.parents[1]/'done.exit'
        assert marker.read_text().strip()=='0'
        parts = []; binary = root/'t4b'/f'{scale}-stencils.bin'
        differences = samples = ghost_differences = ghost_samples = 0
        probe=root/'t4b'/f'{scale}-p8-stencils.bin'
        coarse_chi=coarse_origin=None
        with h5py.File(oldpath) as f, h5py.File(newpath) as full, binary.open('wb') as bout, probe.open('wb') as pout:
            oldnames = [f.attrs[f'component_{i}'].decode() for i in range(int(f.attrs['num_components']))]
            names = [full.attrs[f'component_{i}'].decode() for i in range(int(full.attrs['num_components']))]
            nlev = int(f.attrs['num_levels']); dx0 = float(f['level_0'].attrs['dx'])
            for lev in range(nlev):
                g, q = f[f'level_{lev}'], full[f'level_{lev}']
                b = boxes(g); assert np.array_equal(b,boxes(q))
                h = float(g.attrs['dx']); assert h==float(q.attrs['dx'])
                assert float(q.attrs['time'])==0.
                low = b[:,:2].min(axis=0); high=b[:,2:].max(axis=0)
                fine = boxes(f[f'level_{lev+1}'])//2 if lev+1<nlev else np.empty((0,4),int)
                fo,qo=g['data:offsets=0'][:],q['data:offsets=0'][:]
                allfull=[]
                for i,(x0,y0,x1,y1) in enumerate(b):
                    nx,ny=x1-x0+1,y1-y0+1
                    a=q['data:datatype=0'][int(qo[i]):int(qo[i+1])].reshape(len(names),ny+6,nx+6)
                    old=g['data:datatype=0'][int(fo[i]):int(fo[i+1])].reshape(len(oldnames),ny,nx)
                    for k,n in enumerate(oldnames):
                        aa=a[names.index(n),3:-3,3:-3]
                        differences+=np.count_nonzero(aa.view('u8')!=old[k].view('u8')); samples+=aa.size
                    allfull.append(a)
                seamsx=sorted(set(v[0] for v in b if v[0]!=low[0]))
                seamphysicalx=[round(v*h-224.,15) for v in seamsx]
                seamphysicaly=sorted(set(round(v[1]*h,15) for v in b if v[1]!=low[1]))
                geometry.append(dict(scale=scale,level=lev,dx_M=h,boxes=len(b),
                    xmin_M=low[0]*h-224.,xmax_M=(high[0]+1)*h-224.,
                    ymax_M=(high[1]+1)*h,internal_x_M=';'.join(map(str,seamphysicalx)),
                    internal_y_M=';'.join(map(str,seamphysicaly))))
                for bi,(x0,y0,x1,y1) in enumerate(b):
                    a=allfull[bi]; X,Y=np.meshgrid(np.arange(x0,x1+1),np.arange(y0,y1+1))
                    rho=np.hypot((X+.5)*h-224.,(Y+.5)*h)
                    valid=np.ones(X.shape,bool); near=np.zeros(X.shape,bool)
                    for cx0,cy0,cx1,cy1 in fine:
                        cov=(X>=cx0)&(X<=cx1)&(Y>=cy0)&(Y<=cy1)
                        valid &= ~cov
                        near |= (X>=cx0-2)&(X<=cx1+2)&(Y>=cy0-2)&(Y<=cy1+2)&~cov
                    ex=((X-low[0]<2)|(high[0]-X<2))&(lev>0)
                    ey=((Y-low[1]<2)|(high[1]-Y<2))&(lev>0)&(Y>=2)
                    flags={'convex_corner':ex&ey,
                           'patch_edge':(ex^ey)&~((Y<2)&ex),
                           'box_seam':((X-x0<2)|(x1-X<2)|(Y-y0<2)|(y1-Y<2))&~(ex|ey)&(Y>=2)}
                    nxdom=int(f['level_0'].attrs['prob_domain']['hi_i'])+1
                    nydom=int(f['level_0'].attrs['prob_domain']['hi_j'])+1
                    outer=((X<2)|(X>=nxdom*2**lev-2)|(Y>=nydom*2**lev-2))&(lev==0)
                    flags['box_seam'] &= ~outer
                    # Every overlap is compared to the neighbouring valid owner,
                    # including the covered coarse cells, before composite masking.
                    for bj,(cx0,cy0,cx1,cy1) in enumerate(b):
                        if bj==bi: continue
                        lx,hx=max(x0-3,cx0),min(x1+3,cx1)
                        ly,hy=max(y0-3,cy0),min(y1+3,cy1)
                        if lx>hx or ly>hy: continue
                        av=a[:28,ly-y0+3:hy-y0+4,lx-x0+3:hx-x0+4]
                        bv=allfull[bj][:28,ly-cy0+3:hy-cy0+4,lx-cx0+3:hx-cx0+4]
                        ghost_differences+=np.count_nonzero(av.view('u8')!=bv.view('u8')); ghost_samples+=av.size
                    selected=np.zeros(X.shape,bool)
                    maskflags={m:(rho>=MASKS[m][0])&(rho<=MASKS[m][1]) for m,c in needed}
                    for m,c in needed: selected |= valid&maskflags[m]&flags[c]
                    jj,ii=np.where(selected)
                    if not len(ii): continue
                    records=np.empty((len(ii),3+28*25),dtype='f8'); records[:,0]=h; records[:,1]=(Y[selected]+.5)*h
                    records[:,2]=(X[selected]+.5)*h-224.
                    for z in range(-2,3):
                        for x in range(-2,3):
                            records[:,3+(z+2)*5+x+2::25]=a[:28,jj+3+z,ii+3+x].T
                    records.tofile(bout)
                    # Change only CF ghost chi, never any current valid cell,
                    # same-level ghost, axis ghost, reference or target.
                    alternate=records.copy()
                    for z in range(-2,3):
                        for x in range(-2,3):
                            sx=X[selected]+x; sy=Y[selected]+z
                            cf=((sx<low[0])|(sx>high[0])|(sy>high[1]))&(sy>=0)&(lev>0)
                            if cf.any():
                                alternate[cf,3+(z+2)*5+x+2]=prolong8(sx[cf],sy[cf],coarse_chi,coarse_origin)
                    alternate.tofile(pout)
                    p={'level':np.full(len(ii),lev),'h':np.full(len(ii),h),'y':(Y[selected]+.5)*h,
                       'r':rho[selected],'near_finer':near[selected],
                       'covered':~valid[selected],'axis':Y[selected]<2,
                       'axis_nested_stencil':Y[selected]<4,
                       'central':abs((X[selected]+.5)*h-224.)<2*h,
                       'chi':a[names.index('chi'),jj+3,ii+3],
                       'Ham':a[names.index('Ham'),jj+3,ii+3],
                       'Mom':np.hypot(a[names.index('Mom1'),jj+3,ii+3],a[names.index('Mom2'),jj+3,ii+3]),
                       'GaussE':a[names.index('GaussE'),jj+3,ii+3]}
                    for c in TARGETS: p[c]=flags[c][selected]
                    for m in maskflags: p[m]=maskflags[m][selected]
                    parts.append(p)
                # Rectangular current coarse union, including its saved halo.
                coarse_origin=low-3
                coarse_chi=np.empty((high[1]-low[1]+7,high[0]-low[0]+7))
                for a,(x0,y0,x1,y1) in zip(allfull,b):
                    coarse_chi[y0-low[1]:y1-low[1]+7,x0-low[0]:x1-low[0]+7]=a[names.index('chi')]
        identity.append(dict(scale=scale,valid_numeric_values=samples,valid_bit_differences=int(differences),
                             same_level_ghost_values=ghost_samples,same_level_ghost_bit_differences=int(ghost_differences)))
        assert differences==ghost_differences==0
        output=root/'t4b'/f'{scale}-budgets.bin'
        subprocess.run([str(root/'t4b/roundoff.ex'),str(binary),str(output)],check=True)
        pout=root/'t4b'/f'{scale}-p8-budgets.bin'
        subprocess.run([str(root/'t4b/roundoff.ex'),str(probe),str(pout)],check=True)
        data={k:np.concatenate([p[k] for p in parts]) for k in parts[0]}
        budgets=np.fromfile(output,dtype='f8').reshape(-1,13)
        p8=np.fromfile(pout,dtype='f8').reshape(-1,13)
        assert len(budgets)==len(data['y']) and np.isfinite(budgets).all()
        for row in failed:
            m,c,n=row['mask'],row['cell_class'],row['constraint']; ci=('Ham','Mom','GaussE').index(n)
            choose=data[m]&data[c]; w=data['y'][choose]; value=data[n][choose]
            error=budgets[choose,3+ci]; actual=rms(value,w)
            assert np.count_nonzero(choose)==int(row[f'{scale}_cells'])
            assert math.isclose(actual,float(row[f'{scale}_rms']),rel_tol=3e-14)
            raw.append(dict(scale=scale,mask=m,cell_class=c,constraint=n,cells=len(w),
                raw_rms=actual,min_dx_M=float(data['h'][choose].min()),max_dx_M=float(data['h'][choose].max()),
                roundoff_budget_rms=rms(error,w),residual_over_budget=actual/rms(error,w),
                initial_roundoff_budget_rms=rms(budgets[choose,10+ci],w),
                residual_over_initial_budget=actual/rms(budgets[choose,10+ci],w),
                eps_chi_over_h2_rms=rms(np.finfo(float).eps*abs(data['chi'][choose])/data['h'][choose]**2,w),
                native_replay_difference_rms=rms(budgets[choose,ci]-value,w),
                metric_ricci_rms=rms(budgets[choose,6],w),
                remaining_hamiltonian_rms=rms(budgets[choose,7],w),
                local_metric_span_max=float(budgets[choose,8].max()),
                stored_gamma_replay_difference_rms=rms(budgets[choose,9],w),
                p8_chi_ghost_hamiltonian_rms=rms(p8[choose,0],w),
                p8_chi_ghost_hamiltonian_change_rms=rms(p8[choose,0]-budgets[choose,0],w),
                near_finer_cells=int(np.count_nonzero(data['near_finer'][choose])),
                covered_cells=int(np.count_nonzero(data['covered'][choose])),axis_cells=int(np.count_nonzero(data['axis'][choose])),
                axis_nested_stencil_cells=int(np.count_nonzero(data['axis_nested_stencil'][choose])),
                central_seam_cells=int(np.count_nonzero(data['central'][choose])),
                min_r_M=float(data['r'][choose].min()),max_r_M=float(data['r'][choose].max())))
            for subset, flag in [('central',data['central']),('other_seams',~data['central']),
                                 ('near_finer',data['near_finer']),('away_finer',~data['near_finer'])]:
                ch=choose&flag
                if not ch.any(): continue
                ww=data['y'][ch]
                localization.append(dict(scale=scale,mask=m,cell_class=c,constraint=n,subset=subset,
                    cells=int(ch.sum()),rms=rms(data[n][ch],ww),roundoff_budget_rms=rms(budgets[ch,3+ci],ww)))
    save('t4b-roundoff.csv',raw);save('t4b-layout.csv',geometry)
    save('t4b-identity.csv',identity);save('t4b-localization.csv',localization)
    summary=[]
    for r in sorted(failed,key=lambda r:(r['cell_class'],r['mask'],r['constraint'])):
        v=[next(x for x in raw if x['scale']==s and all(x[k]==r[k] for k in ('mask','cell_class','constraint'))) for s in SCALES]
        row={k:r[k] for k in ('cell_class','mask','constraint','order_low_mid','order_mid_high')}
        for scale,x in zip(SCALES,v):
            for k in ('cells','raw_rms','min_dx_M','max_dx_M','roundoff_budget_rms','initial_roundoff_budget_rms',
                      'residual_over_initial_budget','near_finer_cells','covered_cells','axis_cells',
                      'axis_nested_stencil_cells','p8_chi_ghost_hamiltonian_rms'):
                row[scale+'_'+k]=x[k]
        row['finest_grid_assessment']='FLOOR_COMPATIBLE_WITH_INITIALIZATION_BUDGET' if v[2]['residual_over_initial_budget']<=1 else 'ABOVE_INITIALIZATION_BUDGET'
        summary.append(row)
    save('t4b-case-summary.csv',summary)
    print('All existing valid values and same-level ghost overlaps are bitwise equal.')
    for r in failed:
        rr=[next(x for x in raw if x['scale']==s and x['mask']==r['mask'] and x['cell_class']==r['cell_class'] and x['constraint']==r['constraint']) for s in SCALES]
        print(r['cell_class'],r['mask'],r['constraint'],'budgets',*[f"{x['roundoff_budget_rms']:.3e}" for x in rr],
              'ratios',*[f"{x['residual_over_budget']:.3g}" for x in rr])


if __name__=='__main__':
    main(Path(sys.argv[1] if len(sys.argv)>1 else '/private/tmp/ems-t4'))
