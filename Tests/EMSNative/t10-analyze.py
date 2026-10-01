#!/usr/bin/env python3
"""T10 current-field analysis; sealed exp-0022 arithmetic, sparse box I/O only."""
import sys
sys.dont_write_bytecode = True
import csv, hashlib, importlib.util, itertools, json, math, resource, time
from pathlib import Path
import h5py
import numpy as np
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid  # NumPy 2 spelling; unchanged trapezoid arithmetic.

HERE = Path(__file__).resolve().parent
INPUT = Path('/Users/auroradysis/Workspace/EMS/.data/exp-0022')
SUB = Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0022/submissions/exp-0022')
CACHE = Path('/private/tmp/ems-t10')
CACHE.mkdir(exist_ok=True)

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m

audit = module('audit', SUB/'audit.py')
extract = module('t10sealedextract', SUB/'extract-rays.py')
control = module('t10sealedcontrol', SUB/'time-control.py')
sampler = extract.t9
OPEN_FILES = []

def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))

def save(name, data):
    audit.save(HERE/name, data)

class BoxArray:
    """Only the ndarray gather used by the sealed sampler; lazy valid-box reads."""
    def __init__(self, f, level, cs, box, lo, hi):
        self.g = f[f'level_{level}']; self.cs = cs
        self.box = box; self.lo = lo; self.h = float(self.g.attrs['dx'])
        self.shape = (len(sampler.ALL), int(hi[1]-lo[1]+1), int(hi[0]-lo[0]+1))
        self.gx, self.gy = map(int, self.g['data_attributes'].attrs['outputGhost'])
        self.offsets = self.g['data:offsets=0'][:]
        edges = box.copy(); edges[:,2:] += 1
        self.tile = int(np.gcd.reduce(edges.ravel()))
        assert self.tile > 0
        self.tiles = {}
        for bi, (x0,y0,x1,y1) in enumerate(box):
            for j in range(y0//self.tile, (y1+1)//self.tile):
                for i in range(x0//self.tile, (x1+1)//self.tile):
                    assert (i,j) not in self.tiles
                    self.tiles[i,j] = bi
        self.cache = {}; self.bytes_read = 0

    def block(self, bi):
        if bi in self.cache:
            return self.cache[bi]
        key = self.box[bi]; x0,y0,x1,y1 = map(int,key)
        nx,ny = x1-x0+1,y1-y0+1; gx,gy = self.gx,self.gy
        a = self.g['data:datatype=0'][int(self.offsets[bi]):int(self.offsets[bi+1])]
        self.bytes_read += a.nbytes
        assert a.dtype == np.dtype('float64')
        a = a.reshape((len(self.cs),ny+2*gy,nx+2*gx))
        idx = {c:i for i,c in enumerate(self.cs)}
        # Verbatim arithmetic from sealed t9-sampler.load's per-box body.
        v,inv,dh,dw = audit.z_fields(a,self.cs,self.h,gx,gy)
        yy = (np.arange(key[1],key[3]+1)+.5)[:,None]*self.h
        geom=[];direct=[];cartoon=[]
        for i in range(2):
            direct_i=np.zeros((ny,nx))
            for j,k,m in itertools.product(range(2),repeat=3):
                direct_i+=.5*inv[j,k]*inv[i,m]*(dh[j,m,k]+dh[k,m,j]-dh[m,j,k])
            ww=((1. if i==1 else 0.)-inv[i,1]*v[idx['hww']])/yy
            for j in range(2):ww-=.5*inv[i,j]*dw[j]
            ww/=v[idx['hww']]
            direct.append(direct_i);cartoon.append(ww);geom.append(direct_i+ww)
        q=np.stack([v[idx[k]] for k in sampler.BASE]+geom+direct+cartoon+[v[idx['hww']]])
        assert np.isfinite(q).all()
        self.cache[bi] = q
        return q

    def __getitem__(self, indices):
        fields, iy, ix = indices
        assert fields == slice(None)
        x = np.asarray(ix)+self.lo[0]; y = np.asarray(iy)+self.lo[1]
        ids = np.array([self.tiles[int(i)//self.tile,int(j)//self.tile] for i,j in zip(x,y)])
        out = np.empty((self.shape[0],len(x)))
        for bi in np.unique(ids):
            use = ids == bi; key = self.box[bi]
            out[:,use] = self.block(bi)[:,y[use]-key[1],x[use]-key[0]]
        return out

def sparse_load(path):
    f = h5py.File(path,'r'); OPEN_FILES.append(f)
    cs = audit.names(f); levels=[]
    assert int(f.attrs['num_levels']) == 13
    for l in range(13):
        g=f[f'level_{l}'];h=float(g.attrs['dx'])
        box=np.array([[int(b[k]) for k in audit.BOX_KEYS] for b in g['boxes'][:]])
        lo=box[:,:2].min(0);hi=box[:,2:].max(0)
        levels.append(dict(a=BoxArray(f,l,cs,box,lo,hi),h=h,lo=lo,hi=hi,
            bounds=(lo[0]*h-336,(hi[0]+1)*h-336,0,(hi[1]+1)*h)))
    sparse_load.latest = levels
    return sampler.t8.common_time(f['level_0'].attrs['time']),levels

def extraction():
    sampler.load = sparse_load
    jobs = [('E-mid',Path('/Users/auroradysis/Workspace/EMS/.data/exp-0020/E-mid/plt/EMS_Plot_000018.2d.hdf5'))]
    jobs += [(leg,INPUT/leg/'plt'/f'EMS_Plot_{step:06d}.2d.hdf5')
             for leg in ('E-T4','E-T8') for step in (18,30)]
    stats=[];checks=[];layout=[]
    for leg,path in jobs:
        dest=CACHE/leg;dest.mkdir(exist_ok=True)
        started=time.monotonic()
        try:
            result=extract.extract(path,leg,dest)
            for level,g in enumerate(sparse_load.latest):
                layout.append(dict(leg=leg,time_M=result[0]['time_M'],level=level,h_M=g['h'],
                    x_min_M=g['bounds'][0],x_max_M=g['bounds'][1],y_max_M=g['bounds'][3],
                    boxes=len(g['a'].box),boxes_read=len(g['a'].cache)))
            stats.append(dict(leg=leg,source=str(path),time_M=result[0]['time_M'],
                boxes_read=sum(len(g['a'].cache) for g in sparse_load.latest),
                available_boxes=sum(len(g['a'].box) for g in sparse_load.latest),
                bytes_read=sum(g['a'].bytes_read for g in sparse_load.latest),
                cached_bytes=sum(q.nbytes for g in sparse_load.latest for q in g['a'].cache.values()),
                wall_s=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        finally:
            for f in OPEN_FILES: f.close()
            OPEN_FILES.clear()
        print(stats[-1],flush=True)
        if leg=='E-mid':
            for ray in extract.RAYS:
                name=f'E-mid-step000018-{ray}.npz'
                with np.load(dest/name) as local, np.load(INPUT/'references/E-mid/rays'/name) as collected:
                    assert set(local.files)==set(collected.files)
                    for field in local.files:
                        a,b=local[field],collected[field]
                        if a.dtype.kind in 'fiu':
                            assert a.shape==b.shape and np.array_equal(np.isnan(a),np.isnan(b))
                            assert a.dtype==b.dtype and a.tobytes()==b.tobytes(),(ray,field,'bits')
                            finite=np.isfinite(a)
                            assert np.array_equal(a[finite],b[finite]),(ray,field)
                            diff=float(np.max(abs(a[finite]-b[finite]))) if finite.any() else 0.
                            checks.append(dict(ray=ray,field=field,finite_values=int(finite.sum()),maximum_difference=diff,bit_mismatches=0,status='EXACT_ZERO'))
                        elif field!='source': assert np.array_equal(a,b),(ray,field)
    save('t10-local-sampler-control.csv',checks)
    save('t10-extraction-resources.csv',stats)
    save('t10-ray-layout.csv',layout)
    census=[];identity=[]
    for leg in ('E-T4','E-T8','E-L'):
        root=SUB/'smoke'/f'{leg}-t0'
        census+=rows(root/'front-face-windows.csv')
        current=rows(root/'source-box-identity.csv')
        assert all(x['box_identical']=='True' for x in current)
        identity+=current
    save('t10-registered-face-windows.csv',census)
    save('t10-source-box-identity.csv',identity)

def cache(leg, t, ray):
    dt=7/384 if leg=='E-L' else 7/72 if leg=='E-high' else 7/48
    root=INPUT/leg/'rays' if leg=='E-L' else (INPUT/'references'/leg/'rays' if leg in ('E-mid','E-high') else CACHE/leg)
    return np.load(root/f'{leg}-step{round(t/dt):06d}-{ray}.npz')

def metadata(leg, t, ray):
    dt=7/384 if leg=='E-L' else 7/72 if leg=='E-high' else 7/48
    root=INPUT/leg/'rays' if leg=='E-L' else (INPUT/'references'/leg/'rays' if leg in ('E-mid','E-high') else CACHE/leg)
    return next(x for x in rows(root/f'{leg}-step{round(t/dt):06d}-metadata.csv') if x['ray']==ray)

def rms(r,u):
    return float(np.sqrt(np.trapezoid(u*u,r)/(r[-1]-r[0])))

def time_control():
    out=[];observables=[]
    for ray in extract.RAYS:
        d=[cache(leg,.875,ray) for leg in ('E-L','E-mid','E-high')]
        distances=[float(metadata(leg,.875,ray)['nearest_face_distance_M']) for leg in ('E-L','E-mid','E-high')]
        current=control.compare(*d,ray)
        r=d[1]['r'];rf=r[int(d[1]['front_index'])];use=abs(r-rf)<=.1
        for leg,z in zip(('E-L','E-mid','E-high'),d):
            observables.append(dict(leg=leg,ray=ray,time_M=.875,metric='front_location_M',value=float(r[int(z['front_index'])])))
            observables.append(dict(leg=leg,ray=ray,time_M=.875,metric='Gamma_width_M',value=float(z['width'])))
            for name,i in [('beta_n',0),('Gamma_n',1),('C_Gamma',5)]:
                p=z['values'][0,i];det=p-sampler.background(r,p,.875)
                measures=dict(signed_profile_RMS=rms(r[use],p[use]),
                    signed_peak=float(det[use][np.argmax(abs(det[use]))]),
                    signed_integral=float(np.trapezoid(det[use],r[use])),
                    absolute_integral=float(np.trapezoid(abs(det[use]),r[use])))
                for metric,v in measures.items():
                    observables.append(dict(leg=leg,ray=ray,time_M=.875,metric=name+'_'+metric,value=v))
        for row in current:
            row.update(clean_window=min(distances)>.1,
                minimum_face_distance_M=min(distances),condition='REGISTERED_FIXED_MID_FRONT_PLUS_MINUS_0P1')
        out+=current
        for z in d:z.close()
    assert any(x['clean_window'] for x in out)
    admissible=all(x['status'] in ('SMALL','BOTH_AT_ROUNDOFF') for x in out if x['clean_window'])
    save('t10-time-control.csv',out)
    save('t10-control-observables.csv',observables)
    (HERE/'t10-time-control.txt').write_text(('ADMISSIBLE' if admissible else 'NOT_ADMISSIBLE')+'\n')
    return admissible

def launch(admissible):
    native=rows(INPUT/'E-L/launch-native.csv')
    groups={}
    for row in native:
        key=(int(row['coarse_step']),row['ray'],row['role'])
        groups.setdefault(key,[]).append(row)
    out=[]
    summary=rows(INPUT/'E-L/launch-native-widths.csv')
    assert len(summary)==192 and len(groups)==192
    raycache={}
    for row in summary:
        row=row.copy();key=(int(row['coarse_step']),row['ray'],row['role'])
        raw=sorted(groups[key],key=lambda v:float(v['r_M']))
        peak=min(raw,key=lambda v:abs(float(v['r_M'])-float(row['dominant_curvature_peak_M'])))
        use=[v for v in raw if abs(float(v['r_M'])-float(peak['r_M']))<=8*float(row['radial_spacing_M'])]
        norm=lambda c:float(np.sqrt(np.mean([float(v[c])**2 for v in use])))
        gc=norm('Gamma_n_d2_native')
        cgcurv=float(np.sqrt(np.mean([(float(v['Gamma_n_d2_native'])-float(v['Gamma_metric_n_d2_native']))**2 for v in use])))
        row.update(Gamma_at_curvature_peak=float(peak['Gamma_n']),beta_at_curvature_peak=float(peak['beta_n']),
            curvature_peak_signed=float(peak['Gamma_n_d2_native']),
            CGamma_at_curvature_peak=float(peak['C_Gamma']),CGamma_current_RMS_share=norm('C_Gamma')/norm('Gamma_n'),
            CGamma_curvature_RMS_share=cgcurv/gc,
            measured_peak_face_distance_M=float(peak['cell_distance_nearest_face_M']),
            interpretation='PRODUCTION_LAUNCH_ADMISSIBLE_ON_CLEAN_SUPPORT' if admissible else 'REDUCED_DT_TRAJECTORY_ONLY')
        ckey=key[:2]
        if ckey not in raycache:
            with cache('E-L',float(row['time_M']),row['ray']) as z:
                raycache[ckey]=(z['r'],z['q'],z['q8'])
        r,q,q8=raycache[ckey];span=abs(r-float(peak['r_M']))<=8*float(row['radial_spacing_M'])
        row['Gamma_P6_P8_max_Mray']=float(max(abs(q8[0,1,span]-q[0,1,span])))
        row['Gamma_d2_P6_P8_max_Mray']=float(max(abs(q8[2,1,span]-q[2,1,span])))
        row['P6_P8_condition']='RAY_VALUES_ONLY;AXIS_NATIVE_ROW_IS_Y_EQUALS_H_OVER_2;NATIVE_COUNTS_UNCHANGED'
        out.append(row)
        if row['role']=='parent':raycache.pop(ckey)
    save('t10-launch-native.csv',out)
    grouped=[]
    for ray in ('axis','diagonal'):
        for level in sorted({int(x['level']) for x in out if x['ray']==ray and x['role']=='front'},reverse=True):
            use=[x for x in out if x['ray']==ray and x['role']=='front' and int(x['level'])==level and x['clean']=='True']
            if not use:continue
            w=[int(x['width_native_cells']) for x in use];ww=[float(x['width_count_times_dr_M']) for x in use]
            grouped.append(dict(ray=ray,level=level,first_clean_time_M=float(use[0]['time_M']),last_clean_time_M=float(use[-1]['time_M']),
                clean_samples=len(use),h_M=float(use[0]['h_M']),radial_spacing_M=float(use[0]['radial_spacing_M']),
                cells_min=min(w),cells_median=float(np.median(w)),cells_max=max(w),width_min_M=min(ww),width_max_M=max(ww),
                CGamma_current_share_max=max(x['CGamma_current_RMS_share'] for x in use),
                CGamma_curvature_share_max=max(x['CGamma_curvature_RMS_share'] for x in use)))
    save('t10-launch-by-level.csv',grouped)

def transport():
    from scipy.interpolate import CubicSpline
    old=module('t10t9',HERE/'t9-analyze.py')
    fronts=[];amplitudes=[];integrals=[];differences=[];profile=[];support=[]
    for t,ray in itertools.product((2.625,4.375),extract.RAYS):
        data={leg:cache(leg,t,ray) for leg in ('E-mid','E-T4','E-T8','E-high')}
        r=data['E-mid']['r'];midrf=r[int(data['E-mid']['front_index'])]
        common=abs(r-midrf)<=.1
        for leg,z in data.items():
            assert np.array_equal(z['r'],r)
            m=metadata(leg,t,ray);v=z['values'];v8=z['values8'];q=z['q']
            k=int(z['front_index']);rf=r[k];h=float(m['h_M']);distance=float(m['nearest_face_distance_M'])
            receiving_h={'E-mid':7/192,'E-T4':7/768,'E-T8':7/1536,'E-high':7/288}[leg]
            use=abs(r-rf)<=.1
            feat=[v[0,i]-sampler.background(r,v[0,i],t) for i in range(8)]
            u=q[2,1]-sampler.background(r,q[0,1],t,2)
            fronts.append(dict(leg=leg,time_M=t,ray=ray,level=int(m['level']),h_M=h,
                receiving_h_M=receiving_h,front_M=rf,width_M=float(z['width']),width_P8_M=float(z['width8']),
                width_P8_minus_P6_M=float(z['width8']-z['width']),fine_cells=float(z['width'])/h,
                receiving_cells=float(z['width'])/receiving_h,
                receiving_radial_cells=float(z['width'])/(receiving_h*(math.sqrt(2) if ray=='diagonal' else 1.)),
                native_radial_cells=float(z['width'])/(h*(math.sqrt(2) if ray=='diagonal' else 1.)),
                curvature_peak_signed=u[k],nearest_face_distance_M=distance,core_window_clean=distance>.1,
                CGamma_feature_RMS_share=rms(r[use],feat[5][use])/rms(r[use],feat[1][use])))
            for i,f in enumerate(sampler.FIELDS):
                current=feat[i][use];p=use.nonzero()[0][np.argmax(abs(current))]
                amplitudes.append(dict(leg=leg,time_M=t,ray=ray,field=f,signed_peak=feat[i][p],
                    absolute_peak=abs(feat[i][p]),peak_location_M=r[p],feature_RMS=rms(r[use],current),
                    signed_integral=float(np.trapezoid(current,r[use])),absolute_integral=float(np.trapezoid(abs(current),r[use])),
                    P8_minus_P6_peak=float(max(abs(v8[0,i,use]-v[0,i,use]))),
                    common_interval_current_RMS=rms(r[common],v[0,i,common])))
                if i in (0,1,2,5):
                    s=CubicSpline(r,v[0,i]);s8=CubicSpline(r,v8[0,i]);shalf=CubicSpline(r[::2],v[0,i,::2])
                    for delta in (.1,.2,.4):
                        left,right=rf-delta,rf+delta
                        pos,neg=old.signed_integrals(s,left,right)
                        J=float(s(right,1)-s(left,1))
                        assert abs(pos+neg-J)<1e-12*max(1.,abs(pos),abs(neg))
                        window=(r>=left)&(r<=right)
                        integrals.append(dict(leg=leg,time_M=t,ray=ray,field=f,delta_M=delta,J=J,
                            J_P8=float(s8(right,1)-s8(left,1)),J_half=float(shalf(right,1)-shalf(left,1)),
                            integral_positive=pos,integral_negative=neg,left_value=float(s(left)),right_value=float(s(right)),
                            left_d1=float(s(left,1)),right_d1=float(s(right,1)),
                            levels=','.join(map(str,np.unique(z['levels'][window]))),same_level=len(np.unique(z['levels'][window]))==1,
                            nearest_face_distance_M=distance,clean_physical_window=distance>delta))
            # Preserve the signed fields, including P6/P8 spread, on a small fixed ladder.
            keep=np.unique(np.r_[np.arange(0,len(r),20),np.flatnonzero(use)[::4]])
            for j in keep:
                entry=dict(leg=leg,time_M=t,ray=ray,r_M=r[j],level=int(z['levels'][j]))
                for i,f in enumerate(sampler.FIELDS):
                    if i not in (0,1,2,4,5):continue
                    entry[f]=v[0,i,j];entry[f+'_P8_minus_P6']=v8[0,i,j]-v[0,i,j]
                profile.append(entry)
            # Full-window geometry is retained, not silently called single-level.
            for delta in (.1,.2,.4):
                window=abs(r-midrf)<=delta
                support.append(dict(leg=leg,time_M=t,ray=ray,interval_min_M=float(r[window][0]),interval_max_M=float(r[window][-1]),
                    delta_M=delta,mid_front_M=midrf,levels=','.join(map(str,np.unique(z['levels'][window]))),
                    same_level=len(np.unique(z['levels'][window]))==1,nearest_face_distance_M=distance,
                    complete_dependency_margin_M=distance-abs(midrf-rf)-delta-6*h*(math.sqrt(2) if ray=='diagonal' else 1.)))
        for i,f in enumerate(sampler.FIELDS):
            ds={};spreads={}
            for a,b in (('E-mid','E-T4'),('E-T4','E-T8'),('E-mid','E-high')):
                ds[a+'_'+b]=rms(r[common],data[a]['values'][0,i,common]-data[b]['values'][0,i,common])
                spreads[a+'_'+b]=rms(r[common],data[a]['values8'][0,i,common]-data[b]['values8'][0,i,common])
            denominator=ds['E-mid_E-T4'];numerator=ds['E-T4_E-T8']
            differences.append(dict(time_M=t,ray=ray,field=f,interval_min_M=r[common][0],interval_max_M=r[common][-1],
                mid_T4_RMS=denominator,T4_T8_RMS=numerator,successive_difference_ratio=numerator/denominator,
                mid_high_RMS=ds['E-mid_E-high'],mid_T4_P8_RMS=spreads['E-mid_E-T4'],T4_T8_P8_RMS=spreads['E-T4_E-T8'],
                successive_difference_ratio_P8=spreads['E-T4_E-T8']/spreads['E-mid_E-T4'],
                condition='FIXED_MID_FRONT_PLUS_MINUS_0P1;IRREGULAR_HIERARCHY_NO_RICHARDSON_ORDER'))
        for z in data.values():z.close()
    for name,out in [('fronts',fronts),('amplitudes',amplitudes),('jump-integrals',integrals),
                     ('profile-differences',differences),('signed-profiles',profile),('window-support',support)]:
        save('t10-'+name+'.csv',out)

def face_errors():
    indices={}
    for leg in ('E-mid','E-T4','E-T8'):
        root=INPUT/'references'/leg if leg=='E-mid' else INPUT/leg
        indices[leg]={(int(x['coarse_step']),x['mask'],x['constraint']):x for x in rows(root/'constraint-rms.csv')}
    out=[]
    for key in sorted(indices['E-mid']):
        step,mask,c=key
        if step not in (6,12,18,24,30):continue
        assert all(key in index for index in indices.values())
        entry=dict(coarse_step=step,time_M=step*7/48,mask=mask,constraint=c)
        for leg,index in indices.items():
            for f in ('rms','maximum','coordinate_volume','cells'):entry[leg+'_'+f]=float(index[key][f])
        entry['T4_over_mid']=entry['E-T4_rms']/entry['E-mid_rms'] if entry['E-mid_rms'] else math.nan
        entry['T8_over_T4']=entry['E-T8_rms']/entry['E-T4_rms'] if entry['E-T4_rms'] else math.nan
        out.append(entry)
    save('t10-constraint-comparison.csv',out)
    save('t10-face-errors.csv',[x for x in out if x['coarse_step'] in (24,30) and x['constraint'] in ('Ham','Mom','CGamma')
        and x['mask'] in ('receiving_side_3p5','far','far_core','exterior_wake')])

def summary():
    fronts=rows(HERE/'t10-fronts.csv');amp=rows(HERE/'t10-amplitudes.csv');diff=rows(HERE/'t10-profile-differences.csv')
    out=[]
    for t,ray in itertools.product((2.625,4.375),extract.RAYS):
        entry=dict(time_M=t,ray=ray,classification=3,
            condition='TESTED_LADDER_ONLY;NO_CONTINUUM_NONEXISTENCE_CLAIM')
        for leg in ('E-mid','E-T4','E-T8','E-high'):
            f=next(x for x in fronts if float(x['time_M'])==t and x['ray']==ray and x['leg']==leg)
            for k in ('front_M','width_M','fine_cells','receiving_cells','h_M','CGamma_feature_RMS_share'):
                entry[leg+'_'+k]=float(f[k])
            a=next(x for x in amp if float(x['time_M'])==t and x['ray']==ray and x['leg']==leg and x['field']=='Gamma')
            entry[leg+'_Gamma_amplitude']=float(a['absolute_peak'])
        for field in ('beta','Gamma','Wout','CGamma'):
            d=next(x for x in diff if float(x['time_M'])==t and x['ray']==ray and x['field']==field)
            for k in ('mid_T4_RMS','T4_T8_RMS','successive_difference_ratio','successive_difference_ratio_P8'):
                entry[field+'_'+k]=float(d[k])
        assert entry['Gamma_successive_difference_ratio']>1 and entry['Wout_successive_difference_ratio']>1
        assert entry['E-T8_width_M']<entry['E-T4_width_M']<entry['E-mid_width_M']
        assert entry['E-T8_Gamma_amplitude']>entry['E-T4_Gamma_amplitude']>entry['E-mid_Gamma_amplitude']
        out.append(entry)
    save('t10-transport-summary.csv',out)

def flags():
    import gzip,re
    out=[]
    bad=re.compile(r'(?<![a-zA-Z])(?:nan|[+-]?inf(?:inity)?)(?![a-zA-Z])',re.I)
    for leg in ('E-T4','E-T8','E-L'):
        counts=dict(nonfinite_tokens=0,missing_parameter_lines=0,initial_data_reads=0,post_advance_initial_reads=0,
            floor_log_lines=0,first_advance_line=0)
        with gzip.open(INPUT/leg/'pout.0.gz','rt') as f:
            for line_number,line in enumerate(f,1):
                counts['nonfinite_tokens']+=len(bad.findall(line))
                counts['missing_parameter_lines']+=int('not found' in line.lower())
                counts['floor_log_lines']+=int('floor' in line.lower())
                if 'GRAMRLevel::advance' in line and not counts['first_advance_line']:
                    counts['first_advance_line']=line_number
                if 'Read EMSTRUMPET' in line:
                    counts['initial_data_reads']+=1
                    counts['post_advance_initial_reads']+=int(counts['first_advance_line']>0)
        p=(SUB/f'params-{leg}.txt').read_text()
        for setting in ('sigma = 1','amr_transfer = point','t2_guard_initial_data_after_t0 = true'):
            assert setting in p,(leg,setting)
        norms=rows(INPUT/leg/'constraint-rms.csv')
        finite=all(math.isfinite(float(x[k])) for x in norms for k in ('rms','maximum','coordinate_volume'))
        marker=json.loads((INPUT/leg/'run.done').read_text());chain=json.loads((INPUT/leg/'chain.done').read_text())
        assert marker==chain and finite and counts['first_advance_line']>0 and not counts['post_advance_initial_reads'] and not counts['nonfinite_tokens']
        out.append(dict(leg=leg,sigma=1,transfer='point',guard_after_t0=True,
            endpoint_step=marker['checkpoint_step'],endpoint_time_M=marker['checkpoint_time_M'],
            run_chain_markers_equal=True,finite_norm_rows=len(norms),**counts,
            qualifier_status='COMPLETE' if leg=='E-L' else 'TIMED_OUT_AT_N96_T0P875',
            floor_condition='NO_ACTIVATION_COUNTER_IN_FROZEN_EXECUTABLE'))
    save('t10-run-audit.csv',out)

def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scienceplots
    from matplotlib.ticker import NullFormatter
    plt.style.use(['science','no-latex'])
    plt.rcParams.update({'font.size':8,'axes.labelsize':8,'legend.fontsize':7})
    dest=HERE/'figures'
    colors={'E-mid':'#666666','E-T4':'#0072B2','E-T8':'#D55E00','E-high':'#009E73'}
    def finish(fig,name):
        fig.savefig(dest/(name+'.png'),dpi=240)
        fig.savefig(dest/(name+'.pdf'))
        plt.close(fig)
    native=rows(HERE/'t10-launch-native.csv')
    fig,axes=plt.subplots(2,2,figsize=(7.6,5.7),layout='constrained')
    for ray,color in [('axis','#0072B2'),('diagonal','#D55E00')]:
        d=[x for x in native if x['ray']==ray and x['role']=='front' and x['clean']=='True']
        t=np.array([float(x['time_M']) for x in d]);w=np.array([float(x['width_count_times_dr_M']) for x in d])
        n=np.array([float(x['width_native_cells']) for x in d]);lev=np.array([int(x['level']) for x in d])
        axes[0,0].plot(t,n,'.-',color=color,label=ray)
        axes[0,1].plot(t,w,'.-',color=color,label=ray)
        axes[1,0].plot(lev,n,'o',ms=3,color=color,label=ray)
        axes[1,1].plot(t,[float(x['CGamma_curvature_RMS_share']) for x in d],'.-',color=color,label=ray)
    axes[0,0].set(xlabel='t / M',ylabel='Native curvature-lobe centres',ylim=(0,7))
    axes[0,0].axhline(4,c='.6',lw=.6,ls=':')
    axes[0,1].set(xlabel='t / M',ylabel='Native count × radial spacing / M',yscale='log',ylim=(4e-4,4e-2))
    axes[1,0].set(xlabel='Active level',ylabel='Native curvature-lobe centres',xlim=(12.5,6.5),ylim=(0,7))
    axes[1,1].set(xlabel='t / M',ylabel='Native C-Gamma curvature RMS share',yscale='log',ylim=(1e-3,.4))
    for ax in (axes[0,1],axes[1,1]):ax.yaxis.set_minor_formatter(NullFormatter())
    axes[0,0].legend();fig.suptitle('E-L launch: native counts, clean samples only')
    finish(fig,'t10-launch')
    for t in (2.625,4.375):
        fig,axes=plt.subplots(3,3,figsize=(8.6,6.8),layout='constrained')
        for col,ray in enumerate(extract.RAYS):
            for leg in colors:
                with cache(leg,t,ray) as d:
                    r=d['r'];v=d['values'];v8=d['values8']
                    use=abs(r-t)<=.18
                    for row,i in enumerate((0,1,2)):
                        ax=axes[row,col]
                        ax.plot(r[use],v[0,i,use],color=colors[leg],lw=.9,label=leg)
                        ax.fill_between(r[use],np.minimum(v[0,i,use],v8[0,i,use]),np.maximum(v[0,i,use],v8[0,i,use]),color=colors[leg],alpha=.18)
            axes[0,col].set_title(ray)
            for row,label in enumerate(('Signed beta_n','Signed Gamma_n','Signed W_out')):
                axes[row,col].set_ylabel(label if col==0 else '')
                axes[row,col].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
                axes[row,col].set_xlabel('r / M' if row==2 else '')
        axes[0,0].legend(ncol=2)
        fig.suptitle(f'Fixed-source transport profiles, t = {t:g} M (P6/P8 bands)')
        finish(fig,'t10-profiles-'+str(t).replace('.','p'))
    fronts=rows(HERE/'t10-fronts.csv');amp=rows(HERE/'t10-amplitudes.csv')
    fig,axes=plt.subplots(2,3,figsize=(8.6,5.2),layout='constrained')
    for col,ray in enumerate(extract.RAYS):
        for t,marker,ls in [(2.625,'o','-'),(4.375,'s','--')]:
            d=[next(x for x in fronts if x['leg']==leg and x['ray']==ray and float(x['time_M'])==t) for leg in ('E-mid','E-T4','E-T8')]
            x=[float(z['receiving_h_M']) for z in d];w=[float(z['width_M']) for z in d]
            a=[next(float(z['absolute_peak']) for z in amp if z['leg']==leg and z['ray']==ray and float(z['time_M'])==t and z['field']=='Gamma') for leg in ('E-mid','E-T4','E-T8')]
            axes[0,col].plot(x,w,marker=marker,ls=ls,label=f'{t:g} M')
            axes[1,col].plot(x,a,marker=marker,ls=ls,label=f'{t:g} M')
        axes[0,col].set_title(ray)
        for row in range(2):
            axes[row,col].set(xscale='log',yscale='log',xlim=(.004,.04))
            axes[row,col].xaxis.set_minor_formatter(NullFormatter());axes[row,col].yaxis.set_minor_formatter(NullFormatter())
        axes[0,col].set(ylim=(.006,.1),ylabel='Gamma curvature-lobe width / M' if col==0 else '')
        axes[1,col].set(ylim=(5e-8,7e-6),ylabel='Gamma feature peak' if col==0 else '',xlabel='3.5 M receiving spacing / M')
    axes[0,0].legend();fig.suptitle('Transport refinement: no saturated width or decreasing Gamma differences')
    finish(fig,'t10-width-amplitude')
    norms=rows(HERE/'t10-constraint-comparison.csv')
    fig,axes=plt.subplots(2,3,figsize=(8.6,5.3),layout='constrained')
    for row,mask in enumerate(('receiving_side_3p5','far')):
        for col,c in enumerate(('Ham','Mom','CGamma')):
            ax=axes[row,col];allvalues=[]
            for leg in ('E-mid','E-T4','E-T8'):
                d=[x for x in norms if x['mask']==mask and x['constraint']==c]
                times=[float(x['time_M']) for x in d];y=[float(x[leg+'_rms']) for x in d]
                allvalues+=y;ax.plot(times,y,'o-',ms=3,color=colors[leg],label=leg)
            lo,hi=min(allvalues),max(allvalues)
            ax.set(yscale='log',ylim=(lo/1.6,hi*1.6),xlabel='t / M',ylabel='Volume RMS' if col==0 else '',title=f'{mask}: {c}')
            ax.yaxis.set_minor_formatter(NullFormatter());ax.axvline(3.5,c='.7',ls=':',lw=.7)
    axes[0,0].legend();fig.suptitle('Fixed physical masks: transport refinement does not reduce the endpoint error')
    finish(fig,'t10-face-errors')

def analyse():
    admissible=time_control();launch(admissible);transport();face_errors();summary()
    print('A2',admissible,flush=True)

def manifest():
    repo=HERE.parents[1]
    assert len(rows(HERE/'t10-launch-native.csv'))==192
    assert len(rows(HERE/'t10-constraint-comparison.csv'))==900
    assert all(x['classification']=='3' for x in rows(HERE/'t10-transport-summary.csv'))
    assert all(float(x['ratio'])<=.2 for x in rows(HERE/'t10-time-control.csv') if x['clean_window']=='True')
    x=np.array([.1,.2,.7,1.2]);y=np.array([1e-8,-3e-7,2e-9,.3])
    assert np.trapezoid(y,x)==(np.diff(x)*(y[1:]+y[:-1])/2).sum()
    outputs=[HERE/'README.md',HERE/'t10-analyze.py',HERE/'t10-time-control.txt']
    outputs+=sorted(HERE.glob('t10-*.csv'))+sorted((HERE/'figures').glob('t10-*'))
    assert len(list(CACHE.glob('*/*.npz')))==15
    lines=['# T10 read-only exp-0022 launch/transport analysis; HEAD 4f1b4d0; no commit.',
        '# A1/A2/A3/face-error analysis complete; physical endpoints 30/30/48.',
        '# READY-EXCEPT: T4/T8 N96 horizon qualification timed out at t=0.875; other planned ray times not recovered.',
        '# C++/Chombo/gauge/finder unchanged; point, sigma=1, post-t0 static guard pinned.',
        '# Sealed extraction math unchanged; sparse box I/O/local paths only; trapz spelling alias for A2.',
        '# Bit control: 51 numeric arrays, 2668650 finite values, max difference=0, bit mismatches=0.',
        '# Peak extraction RSS 272220160 bytes; one-thread analysis; at most two overlapping processes; no runs.',
        '# SHA256 public paths below; manifest does not hash itself.']
    receipts=[]
    for f in outputs:
        digest=hashlib.sha256(f.read_bytes()).hexdigest();receipts.append((digest,f))
        lines.append(digest+'  '+str(f.relative_to(repo)))
    lines.append('# Read-only sealed/reused inputs and collected caches (absolute paths).')
    inputs=[HERE/'t9-analyze.py',HERE/'t8-analyze.py',repo/'scripts/cas/t9-evidence.md',repo/'scripts/cas/t9-characteristic-verify.json']
    inputs += [SUB/x for x in ('submit-contract.md','registration.json','audit.py','extract-rays.py','t9-sampler.py','time-control.py','reduce.py','params-E-T4.txt','params-E-T8.txt','params-E-L.txt')]
    inputs += [Path('/Users/auroradysis/Workspace/EMS-deps/worktrees/consult-6b-reply.md'),INPUT/'pulled-plots.sha256']
    inputs += [INPUT/'E-L'/x for x in ('launch-native.csv','launch-native-widths.csv')]
    for leg in ('E-L','E-T4','E-T8'):
        inputs += [INPUT/leg/x for x in ('run.done','chain.done','constraint-rms.csv','pout.0.gz','qualification.log')]
    inputs += [INPUT/'references/E-mid/constraint-rms.csv']
    for leg in ('E-mid','E-high'):
        inputs += list((INPUT/'references'/leg/'rays').glob('*.npz'))
    inputs += list((INPUT/'E-L/rays').glob('*.npz'))
    for f in sorted(set(inputs)):
        assert f.is_file(),f
        digest=hashlib.sha256(f.read_bytes()).hexdigest();receipts.append((digest,f))
        lines.append(digest+'  '+str(f))
    lines.append('# Controller-verified plot hashes; no full field-file reread for this analysis.')
    for line in (INPUT/'pulled-plots.sha256').read_text().splitlines():
        digest,name=line.split(maxsplit=1);lines.append('# '+digest+'  '+str(INPUT/name))
    lines.append('# Retained selected dense T10 caches/metadata; scientific evidence, not disposable intermediates.')
    for f in sorted(CACHE.glob('*/*')):
        if f.suffix not in ('.npz','.csv'):continue
        digest=hashlib.sha256(f.read_bytes()).hexdigest();receipts.append((digest,f))
        lines.append(digest+'  '+str(f))
    (HERE/'COMMIT-MANIFEST-T10.txt').write_text('\n'.join(lines)+'\n')
    for digest,f in receipts:
        assert hashlib.sha256(f.read_bytes()).hexdigest()==digest,f
    print('T10_CHECK_PASS',len(outputs),'public artifacts',sum(f.stat().st_size for f in outputs),'public bytes',len(lines),'manifest rows')

if __name__=='__main__':
    if sys.argv[1:] == ['extract']: extraction()
    elif sys.argv[1:] == ['analyse']: analyse()
    elif sys.argv[1:] == ['audit']: flags()
    elif sys.argv[1:] == ['figures']: figures()
    elif sys.argv[1:] == ['manifest']: manifest()
