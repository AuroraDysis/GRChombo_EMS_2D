#!/usr/bin/env python3
"""Gauge-labelled native frames, I8/I10 ray profiles and fixed-radius histories.

No evolution, static input, fitted phase or gauge-dependent B comparison.
The native replay executable supplies geometry and RHS; this file only samples.
"""
import argparse, csv, hashlib, json, math, resource, struct, subprocess, sys
from pathlib import Path
import numpy as np
from collections import defaultdict
sys.dont_write_bytecode=True

GEOMETRY='R_area_density Rdot_total Rdot_preKO dRdx dRdy X_total X_preKO m_MS_proxy angular_anisotropy radial_tangent_cross geometry_valid'.split()
ODD_Y={2,4,6,7}

def frames(path):
    meta=json.loads(path.with_suffix('.json').read_text())
    n=len(meta['fields']); dtype=np.dtype([('ij','i4',(2,)),('values','f8',(4*n+11,))])
    with path.open('rb') as f:
        while True:
            magic=f.read(8)
            if not magic:break
            if magic!=b'T21RHS01':raise ValueError(f'{path}: bad/truncated frame magic')
            raw=f.read(76)
            if len(raw)!=76:raise ValueError(f'{path}: truncated header')
            l,rank,stage,call,nf,count,snapshot=struct.unpack('=7I',raw[:28])
            t,start,h,dt,cx,cy=struct.unpack('=6d',raw[28:])
            if nf!=n:raise ValueError(f'{path}: metadata field mismatch')
            a=np.fromfile(f,dtype=dtype,count=count)
            if len(a)!=count:raise ValueError(f'{path}: truncated cell payload')
            if not np.isfinite(a['values'][:,:4*n]).all():
                raise ValueError(f'{path}: nonfinite native state or RHS; no coarser fallback permitted')
            yield dict(level=l,rank=rank,stage=stage,call=call,time=t,start=start,h=h,dt=dt,
                       centre=(cx,cy),snapshot=bool(snapshot),meta=meta,data=a)

def labels(meta):
    fields=[f'{s}[{meta["gauge"]}:{meta["driver_semantics"]}]' if s in ('B1','B2') else s for s in meta['fields']]
    return [f'{part}:{s}' for part in ('state','preKO','KO','total') for s in fields]+GEOMETRY

def save_csv(path,rows):
    if not rows:raise ValueError(f'empty output {path}')
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def stages(args):
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for p in sorted(Path(args.frames).glob('*.bin')):
        for f in frames(p):
            if f['snapshot']:continue
            a=f['data']['values'];meta=f['meta'];n=len(meta['fields'])
            y=(f['data']['ij'][:,1]+.5)*f['h']-f['centre'][1]
            volume=2*np.pi*np.abs(y)*f['h']**2
            def norm(v):return float(np.sqrt(np.sum(volume*v*v)/np.sum(volume)))
            for i,s in enumerate(meta['fields']):
                pre,ko,total=a[:,n+i],a[:,2*n+i],a[:,3*n+i]
                scale=np.maximum(1.,np.maximum(np.abs(pre),np.maximum(np.abs(ko),np.abs(total))))
                reconstruction=float(np.max(np.abs(pre+ko-total)/scale))
                # Sequential direction-wise KO additions can round differently from pre+KO.
                if reconstruction>128*np.finfo(float).eps:
                    raise ValueError(f'{p}: native RHS reconstruction exceeds roundoff, {s}')
                rows.append(dict(file=p.name,level=f['level'],rank=f['rank'],stage=f['stage'],call=f['call'],
                    time_M=f['time'],level_start_M=f['start'],field=labels(meta)[i][6:],gauge=meta['gauge'],
                    driver_semantics=meta['driver_semantics'],native_kernel_type=meta['native_kernel_type'],
                    cells=len(a),volume_weight_sum=float(np.sum(volume)),preKO_RMS=norm(pre),KO_RMS=norm(ko),
                    total_RMS=norm(total),KO_fraction_of_total_RMS=norm(ko)/norm(total) if norm(total) else '',
                    cancellation_amplification=(norm(pre)+norm(ko))/norm(total) if norm(total) else '',
                    reconstruction_scaled_max=reconstruction,geometry_invalid=int(np.count_nonzero(a[:,-1]!=1))))
    save_csv(out/'t21-stage-summary.csv',rows)
    groups=defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in ('level','stage','call','time_M','level_start_M','field','gauge','driver_semantics'))].append(row)
    aggregate=[]
    for key,group in groups.items():
        total_weight=sum(a['volume_weight_sum'] for a in group)
        norms={k:math.sqrt(sum(a['volume_weight_sum']*a[k]**2 for a in group)/total_weight) for k in ('preKO_RMS','KO_RMS','total_RMS')}
        aggregate.append(dict(zip(('level','stage','call','time_M','level_start_M','field','gauge','driver_semantics'),key),
            cells=sum(a['cells'] for a in group),volume_weight_sum=total_weight,**norms,
            KO_fraction_of_total_RMS=norms['KO_RMS']/norms['total_RMS'] if norms['total_RMS'] else '',
            geometry_invalid=sum(a['geometry_invalid'] for a in group),
            reconstruction_scaled_max=max(a['reconstruction_scaled_max'] for a in group)))
    save_csv(out/'t21-stage-global-summary.csv',aggregate)
    return dict(frames=len({(r['file'],r['call']) for r in rows}),rows=len(rows),gauge_labels=sorted({r['gauge'] for r in rows}))

def load_grids(directory):
    groups={};meta=None
    for p in sorted(directory.glob('*.bin')):
        for f in frames(p):
            if not f['snapshot']:continue
            if meta is None:meta=f['meta']
            if meta!=f['meta']:raise ValueError('mixed gauge, field or native kernel metadata')
            groups.setdefault(f['level'],[]).append(f)
    if not groups:raise ValueError(f'{directory}: no snapshot frames')
    result=[];clock=None
    for level,group in sorted(groups.items(),reverse=True):
        f=group[0];ij=np.concatenate([g['data']['ij'] for g in group]);values=np.concatenate([g['data']['values'] for g in group])
        if clock is None:clock=f['time']
        if any(abs(g['time']-clock)>1e-9 or g['h']!=f['h'] or g['centre']!=f['centre'] for g in group):
            raise ValueError('mixed snapshot clocks or geometry')
        unique=np.unique(ij,axis=0)
        if len(unique)!=len(ij):raise ValueError('duplicate valid native cells; use one replay attempt')
        lo=ij.min(0);hi=ij.max(0)
        grid=np.full((hi[1]-lo[1]+1,hi[0]-lo[0]+1,values.shape[1]),np.nan)
        grid[ij[:,1]-lo[1],ij[:,0]-lo[0]]=values
        n=len(meta['fields']);sign=[-1 if c in ODD_Y else 1 for c in meta['parities']]*4+[1,1,1,1,-1,1,1,1,1,-1,1]
        result.append(dict(level=level,time=f['time'],h=f['h'],centre=f['centre'],lo=lo,hi=hi,grid=grid,sign=np.array(sign),n=n))
    return clock,meta,result

def weights(x,n):
    nodes=np.arange(math.floor(x)-(n//2-1),math.floor(x)+(n//2+1),dtype=int)
    w=np.ones(n)
    for k in range(n):
        for j in range(n):
            if k!=j:w[k]*=(x-nodes[j])/(nodes[k]-nodes[j])
    return nodes,w

def sample(g,x,y,n):
    ix,wx=weights((x+g['centre'][0])/g['h']-.5,n)
    jy,wy=weights((y+g['centre'][1])/g['h']-.5,n)
    original=jy.copy();jy=np.where(jy<0,-jy-1,jy)
    if ix.min()<g['lo'][0] or ix.max()>g['hi'][0] or jy.min()<g['lo'][1] or jy.max()>g['hi'][1]:return None
    a=g['grid'][jy[:,None]-g['lo'][1],ix[None,:]-g['lo'][0]].copy()
    # The symmetry axis is y=0 in these EMS cartoon runs.
    a[original<0]*=g['sign']
    if not np.isfinite(a[:,:,:4*g['n']]).all():return None
    return np.einsum('i,j,ijv->v',wy,wx,a)

def profiles(args):
    cfg=json.loads(Path(args.config).read_text());out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    time,meta,grids=load_grids(Path(args.frames));names=labels(meta)
    radii=np.unique(np.r_[np.geomspace(cfg['ray_radius_min_M'],cfg['ray_radius_max_M'],cfg['ray_points']),cfg['fixed_radii_M']])
    dirs=cfg['rays'];values=np.full((2,len(dirs),len(radii),len(names)),np.nan)
    levels=np.full((len(dirs),len(radii)),-1,dtype='i4');steps=np.full(levels.shape,np.nan)
    for q,d in enumerate(dirs):
        for i,r in enumerate(radii):
            x,y=np.array(d['direction'])*r
            for g in grids:
                ten=sample(g,x,y,10)
                if ten is None:continue
                eight=sample(g,x,y,8)
                if eight is None:raise ValueError('I8 missing where I10 exists')
                values[:,q,i]=eight,ten;levels[q,i]=g['level'];steps[q,i]=g['h'];break
    if np.any(levels<0):raise ValueError('missing I10 support at registered radius; add replay levels, do not move radii')
    spread=np.abs(values[0]-values[1]);floor=128*np.finfo(float).eps*np.maximum(1.,np.abs(values[0]))
    fixed=np.array([int(np.flatnonzero(radii==r)[0]) for r in cfg['fixed_radii_M']])
    stem=f't21-profiles-step{args.step:06}'
    np.savez_compressed(out/(stem+'.npz'),time_M=time,step=args.step,radii_M=radii,
        rays=np.array([d['name'] for d in dirs]),fields=np.array(names),I8=values[0],I10=values[1],
        spread=spread,floor=floor,levels=levels,h_M=steps,gauge=meta['gauge'],driver_semantics=meta['driver_semantics'],
        native_level_times_M=json.dumps({g['level']:g['time'] for g in grids}))
    np.savez_compressed(out/f't21-fixed-step{args.step:06}.npz',time_M=time,step=args.step,
        radii_M=radii[fixed],rays=np.array([d['name'] for d in dirs]),fields=np.array(names),
        I8=values[0][:,fixed],I10=values[1][:,fixed],spread=spread[:,fixed],floor=floor[:,fixed],
        levels=levels[:,fixed],h_M=steps[:,fixed],gauge=meta['gauge'],driver_semantics=meta['driver_semantics'])
    # Angular spread is retained as a diagnostic, never an average or an error bound.
    geom_indices=[names.index(s) for s in GEOMETRY]
    rows=[]
    for k in geom_indices:
        finite=np.isfinite(values[0,:,:,k]).all(0)
        rows.append(dict(step=args.step,time_M=time,quantity=names[k],radii=len(radii),
            invalid_radii=int(np.count_nonzero(~finite)),
            ray_disagreement_RMS=float(np.sqrt(np.mean(np.ptp(values[0,...,k][:,finite],axis=0)**2))) if finite.any() else '',
            interpolation_spread_RMS=float(np.sqrt(np.nanmean(spread[:,:,k]**2))),
            interpretation='angular area-density proxy; exact areal/Misner-Sharp reading requires spherical symmetry'))
    save_csv(out/(stem+'-geometry.csv'),rows)
    return dict(step=args.step,time_M=time,rays=len(dirs),radii=len(radii),columns=len(names),metadata=meta,
                native_level_times_M={g['level']:g['time'] for g in grids})

def run(args):
    """Compute-node sequential replay of retained checkpoints; no launcher hidden here."""
    cfg=json.loads(Path(args.config).read_text());root=Path(args.checkpoints);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    rows=[]
    def digest(path):
        with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
    packet_hashes={s:digest(getattr(args,s)) for s in ('config','params','executable')}
    for p in sorted(root.glob(args.pattern)):
        # A completed checkpoint can be selected without loading its field arrays in Python.
        import h5py
        with h5py.File(p) as f:
            t=float(f['level_0'].attrs['time']);h=float(f['level_0'].attrs['dx'])
        step=round(t/(h*cfg['dt_multiplier']))
        if abs(t-step*h*cfg['dt_multiplier'])>1e-9:raise ValueError(f'{p}: not a native coarse clock')
        d=out/f'step{step:06}';d.mkdir(exist_ok=True)
        if (d/'done.exit').exists():
            receipt=json.loads((d/'t21-replay.receipt.json').read_text())
            if receipt['returncode']!=0 or receipt['packet_hashes']!=packet_hashes:raise ValueError(f'{d}: failed or different prior replay retained')
            rows.append(receipt);continue
        overrides=[f'restart_file={p.resolve()}','ems_data_path=/T21_STATIC_FILE_MUST_NOT_BE_READ',
                   't21_rhs_capture=true',f't21_rhs_prefix={d}/native']
        command=[*json.loads(args.launch_prefix),str(Path(args.executable).resolve()),str(Path(args.params).resolve()),*overrides]
        if (d/'native.done.exit').exists():
            receipt=json.loads((d/'t21-replay.receipt.json').read_text());rc=receipt['returncode']
            if receipt['packet_hashes']!=packet_hashes:raise ValueError(f'{d}: input packet changed; use a new replay attempt')
        else:
            with (d/'run.log').open('wb') as log:
                rc=subprocess.run(command,cwd=d,stdout=log,stderr=subprocess.STDOUT).returncode
            receipt=dict(step=step,time_M=t,checkpoint=str(p),checkpoint_sha256=digest(p),command=command,returncode=rc,packet_hashes=packet_hashes)
            (d/'t21-replay.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
            (d/'native.done.exit').write_text(str(rc)+'\n')
        if rc:
            (d/'done.exit').write_text(str(rc)+'\n');raise RuntimeError(f'{d}: native replay exit {rc}')
        receipt['profiles']=profiles(argparse.Namespace(frames=str(d),config=args.config,output=str(out/'profiles'),step=step))
        audits=list(d.glob('replay-valid-*.csv'))
        if not audits:raise ValueError('missing native replay identity receipts')
        for a in audits:
            for r in csv.DictReader(a.open()):
                if any(int(r[k]) for k in ('bit_mismatches','nonfinite','advances','static_reads')):raise ValueError(f'{a}: invalid native receipt')
        (d/'t21-replay.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        stamp=d/'done.exit.tmp';stamp.write_text('0\n');stamp.replace(d/'done.exit');rows.append(receipt)
    if not rows:raise ValueError('no retained checkpoints matched')
    (out/'t21-history.receipt.json').write_text(json.dumps(rows,indent=2)+'\n')
    return dict(checkpoints=len(rows),first_step=rows[0]['step'],last_step=rows[-1]['step'])

def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='mode',required=True)
    p=sub.add_parser('stages');p.add_argument('--frames',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('profiles');p.add_argument('--frames',required=True);p.add_argument('--config',required=True);p.add_argument('--output',required=True);p.add_argument('--step',type=int,required=True)
    p=sub.add_parser('run')
    for s in ('checkpoints','config','output','executable','params'):p.add_argument('--'+s,required=True)
    p.add_argument('--pattern',default='EMS_*.2d.hdf5');p.add_argument('--launch-prefix',default='[]',help='JSON command prefix, e.g. the frozen MPI launcher')
    a=ap.parse_args();result=globals()[a.mode](a)
    result['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
