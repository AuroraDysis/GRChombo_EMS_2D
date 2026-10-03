#!/usr/bin/env python3
"""Consume every local child's own receipt, then test evolved bits and recorder output."""
import csv, importlib.util, json, resource, sys
from pathlib import Path
import h5py
import numpy as np
import mpmath as mp
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t21')
spec=importlib.util.spec_from_file_location('t21_record',HERE/'t21-record.py');r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)

def own(directory,name):
    a=json.loads((directory/(name+'.resources.json')).read_text())
    assert (directory/'done.exit').read_text().strip()=='0',(directory,'done marker')
    assert a['returncode']==a['child_measurement']['returncode']==0 and not a['gate_reason'],a
    assert a['peak_rss_bytes']<6_000_000_000,a
    return dict(case=name,returncode=a['returncode'],gate_reason=a['gate_reason'],peak_RSS_bytes=a['peak_rss_bytes'],wall_seconds=a['wall_seconds'])

def compare(a,b):
    checked=mismatches=0
    with h5py.File(a) as f,h5py.File(b) as g:
        assert int(f.attrs['num_levels'])==int(g.attrs['num_levels'])
        assert int(f.attrs['num_components'])==int(g.attrs['num_components'])
        for i in range(int(f.attrs['num_components'])):assert f.attrs[f'component_{i}']==g.attrs[f'component_{i}']
        for l in range(int(f.attrs['num_levels'])):
            x,y=f[f'level_{l}'],g[f'level_{l}']
            assert x.attrs['dx']==y.attrs['dx'] and x.attrs['time']==y.attrs['time']
            assert np.array_equal(x['boxes'][:],y['boxes'][:])
            assert np.array_equal(x['data:offsets=0'][:],y['data:offsets=0'][:])
            for start in range(0,len(x['data:datatype=0']),1_000_000):
                xx=x['data:datatype=0'][start:start+1_000_000].view('u8')
                yy=y['data:datatype=0'][start:start+1_000_000].view('u8')
                checked+=len(xx);mismatches+=int(np.count_nonzero(xx!=yy))
    assert mismatches==0,(a,b,mismatches)
    return checked

receipts=[own(ROOT/'build','build')];bits=[]
for case in ('original','omitted','off','on'):
    receipts.append(own(ROOT/'controls'/case,case))
    assert 'GRChombo finished.' in (ROOT/'controls'/case/'run.log').read_text()
for case in ('omitted','off','on'):
    for prefix in ('chk/EMS_','plt/EMS_Plot_'):
        for step in (0,1):
            file=f'{prefix}{step:06}.2d.hdf5'
            a=ROOT/'controls/original'/file;b=ROOT/'controls'/case/file
            count=compare(a,b)
            bits.append(dict(case=case,file=file,Float64_values=count,bit_mismatches=0,
                scope='production RK4 two-level fixture; initial and one synchronized coarse step; no scientific launch claim'))
assert not list((ROOT/'controls/omitted').glob('native*'))
assert not list((ROOT/'controls/off').glob('native*'))
stage=r.stages(type('Args',(),dict(frames=str(ROOT/'controls/on'),output=str(ROOT/'qualification')))())
mp.mp.dps=50;geometry_samples=0;worst_X=mp.mpf(0)
for path in sorted((ROOT/'controls/on').glob('*.bin')):
    for frame in r.frames(path):
        names=frame['meta']['fields'];index={s:i for i,s in enumerate(names)};n=len(names)
        for cell in frame['data'][::max(1,len(frame['data'])//4)]:
            values=cell['values'];geometric=values[4*n:]
            if geometric[-1]!=1:continue
            v={s:mp.mpf(float(values[k])) for s,k in index.items()}
            spatial=mp.matrix([[v['h11'],v['h12'],0],[v['h12'],v['h22'],0],[0,0,v['hww']]])/v['chi']
            beta=mp.matrix([v['shift1'],v['shift2'],0]);metric=mp.matrix(4)
            metric[0,0]=-v['lapse']**2+(beta.T*spatial*beta)[0]
            for i in range(3):
                metric[0,i+1]=metric[i+1,0]=(spatial*beta)[i]
                for j in range(3):metric[i+1,j+1]=spatial[i,j]
            derivative=mp.matrix([mp.mpf(float(geometric[i])) for i in (1,3,4)]+[0])
            brute=(derivative.T*(metric**-1)*derivative)[0]
            spatial_part=(derivative[1:4,0].T*(spatial**-1)*derivative[1:4,0])[0]
            scale=max(mp.mpf(1),abs(spatial_part),abs(brute),abs(mp.mpf(float(geometric[5]))))
            discrepancy=abs(brute-mp.mpf(float(geometric[5])))/scale
            assert discrepancy<128*np.finfo(float).eps,('native X versus independent 4x4 inverse',discrepancy)
            worst_X=max(worst_X,discrepancy);geometry_samples+=1
assert geometry_samples>0
cfg=json.loads((HERE/'t21-recording.json').read_text());cfg.update(ray_radius_min_M=.05,ray_radius_max_M=1.,ray_points=64,fixed_radii_M=[.1,.2,.3,.5,.9],levels=[0,1])
config=ROOT/'qualification/fixture-recording.json';config.write_text(json.dumps(cfg)+'\n')
for step in (0,1):
    directory=ROOT/'replay'/f'step{step:06}';receipts.append(own(directory,f'replay-{step}'))
    for p in directory.glob('replay-valid*.csv'):
        for a in csv.DictReader(p.open()):
            assert all(int(a[k])==0 for k in ('bit_mismatches','nonfinite','advances','static_reads')),a
    result=r.profiles(type('Args',(),dict(frames=str(directory),config=str(config),output=str(ROOT/'qualification/profiles'),step=step))())
    assert result['metadata']['gauge']=='experimental'
    assert 'ExperimentalGauge' in result['metadata']['native_kernel_type']
    with np.load(ROOT/f'qualification/profiles/t21-profiles-step{step:06}.npz') as data:
        assert any(s=='state:B1[experimental:integrated_offset_decaying_0p1]' for s in data['fields'])
        assert np.isfinite(data['I8'][...,:112]).all()
        assert np.all(data['levels']>=0)
r.save_csv(HERE/'t21-recording-identity.csv',bits);r.save_csv(HERE/'t21-local-resources.csv',receipts)
summary=dict(status='RECORDING_BIT_IDENTITY_PASS',files=len(bits),checked_values=sum(a['Float64_values'] for a in bits),
    bit_mismatches=0,native_RHS_KO_capture=stage,checkpoint_replay_numerical_t0_and_positive_time=True,
    native_geometry_inversion_samples=geometry_samples,native_X_versus_50digit_inverse_scaled_max=str(worst_X),
    static_path_absent_during_replay=True,peak_RSS_bytes=max(a['peak_RSS_bytes'] for a in receipts),
    qualifier_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    scope='serial production-RK4 recording qualification only; MPI and target hierarchy preflight required on one exclusive node',
    level15_status='SUPERSEDED; no initialization or evolution performed')
(HERE/'t21-qualification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
