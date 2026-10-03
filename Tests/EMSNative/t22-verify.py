#!/usr/bin/env python3
"""Verify native receipts first, then scientific identities. No static reads."""
import argparse, csv, importlib.util, json, math, resource, subprocess, sys
from pathlib import Path
import h5py
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t22');EPS=np.finfo(float).eps
FIELDS='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split()
N=28;COL=4*N+12
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(name,rows):
    assert rows,name
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def own(directory,name,marker=True):
    a=json.loads((directory/(name+'.resources.json')).read_text())
    if marker:assert (directory/'done.exit').read_text().strip()=='0',(directory,'marker')
    assert a['returncode']==a['child_measurement']['returncode']==0 and not a['gate_reason'],a
    assert a['peak_rss_bytes']<6_000_000_000,a
    return dict(case=name,returncode=a['returncode'],gate_reason=a['gate_reason'],peak_RSS_bytes=a['peak_rss_bytes'],wall_seconds=a['wall_seconds'])
def compare(a,b,skip=()):
    count=mismatches=0
    with h5py.File(a) as f,h5py.File(b) as g:
        assert int(f.attrs['num_levels'])==int(g.attrs['num_levels'])
        for l in range(int(f.attrs['num_levels'])):
            x,y=f[f'level_{l}'],g[f'level_{l}'];assert x.attrs['dx']==y.attrs['dx'] and x.attrs['time']==y.attrs['time']
            assert np.array_equal(x['boxes'][:],y['boxes'][:])
            assert np.array_equal(x['data:offsets=0'][:],y['data:offsets=0'][:])
            if skip:
                off=x['data:offsets=0'][:];n=int(f.attrs['num_components'])
                for start,end in zip(off[:-1],off[1:]):
                    xx=x['data:datatype=0'][int(start):int(end)].reshape(n,-1)
                    yy=y['data:datatype=0'][int(start):int(end)].reshape(n,-1)
                    take=[i for i in range(n) if i not in skip]
                    count+=xx[take].size;mismatches+=int(np.count_nonzero(xx[take].copy().view('u8')!=yy[take].copy().view('u8')))
            else:
                for start in range(0,len(x['data:datatype=0']),1_000_000):
                    xx=x['data:datatype=0'][start:start+1_000_000].view('u8');yy=y['data:datatype=0'][start:start+1_000_000].view('u8')
                    count+=len(xx);mismatches+=int(np.count_nonzero(xx!=yy))
    assert mismatches==0,(a,b,mismatches)
    return dict(Float64_values=int(count),bit_mismatches=mismatches)
def formulas(q,gauge):
    assert q.ndim==2 and q.shape[1]==COL and np.isfinite(q).all()
    assert np.all(q[:,-1]==0),'native total RHS differs bitwise from orchestration'
    pre=q[:,N:2*N];ko=q[:,2*N:3*N];total=q[:,3*N:4*N]
    geo=q[:,4*N:4*N+2];matter=q[:,4*N+2:4*N+4]
    gamma=pre[:,11:13]
    defect=np.abs(gamma-(geo+matter));budget=128*EPS*np.maximum(1,np.abs(geo)+np.abs(matter))
    assert np.all(defect<=budget),('Gamma stress-tensor reconstruction',np.max(defect/budget))
    d=dict(native_samples=len(q),nonzero_momentum_samples=int(np.count_nonzero(np.linalg.norm(q[:,4*N+4:4*N+6],axis=1)>0)),
        matter_Gamma_peak=float(np.max(abs(matter))),Gamma_sum_budget_fraction=float(np.max(defect/budget)),production_total_bit_mismatches=0)
    if gauge=='moving_puncture':
        expected=np.c_[q[:,4*N+10],q[:,4*N+8:4*N+10],q[:,4*N+6:4*N+8]]
        actual=pre[:,13:18]
        residual=np.abs(actual-expected);scale=128*EPS*np.maximum(1,np.abs(actual)+np.abs(expected))
        assert np.all(residual<=scale),('MPG equations',np.max(residual/scale))
        # The post-KO driver relation is Bdot = Gamma_preKO - eta B + KO(B)
        # (with the configured advection terms); not Gamma_total + KO(B).
        assert np.all(np.abs(total[:,16:18]-(expected[:,3:5]+ko[:,16:18]))<=128*EPS*np.maximum(1,np.abs(total[:,16:18])+np.abs(expected[:,3:5])+np.abs(ko[:,16:18])))
        d.update(gauge_formula_budget_fraction=float(np.max(residual/scale)),B_equation_KO_convention='full Gamma before KO, ordinary KO(B) once')
    else:d.update(gauge_formula_budget_fraction='',B_equation_KO_convention='legacy integrated offset; not an MPG driver')
    return d
def restart_probe(exe,param,name):
    with Path('probe.log').open('wb') as f:p=subprocess.run([exe,param],stdout=f,stderr=subprocess.STDOUT)
    text=Path('probe.log').read_text();rejected=name.startswith('restart-cross')
    passed=(p.returncode!=0 and 'EMS checkpoint gauge does not match ems_gauge' in text) if rejected else (p.returncode==0 and 'GRChombo finished.' in text)
    Path('probe.receipt.json').write_text(json.dumps(dict(actual_child_returncode=p.returncode,expected='gauge mismatch refusal' if rejected else 'same-gauge guarded restart succeeds',pass_control=passed,static_path='/T22_STATIC_FILE_MUST_NOT_BE_READ'),indent=2)+'\n')
    assert passed,(name,p.returncode,text[-1000:])
def quick():
    receipts=[own(ROOT/'build','build')];bits=[]
    for case in ('baseline','omitted','explicit','capture','moving','moving-capture'):
        d=ROOT/'controls'/case;receipts.append(own(d,case));assert 'GRChombo finished.' in (d/'run.log').read_text()
    for case in ('omitted','explicit','capture'):
        for prefix in ('chk/EMS_','plt/EMS_Plot_'):
            for step in (0,1):
                file=f'{prefix}{step:06}.2d.hdf5'
                bits.append(dict(case=case,file=file,**compare(ROOT/'controls/baseline'/file,ROOT/'controls'/case/file)))
    for prefix in ('chk/EMS_','plt/EMS_Plot_'):
        for step in (0,1):
            file=f'{prefix}{step:06}.2d.hdf5'
            bits.append(dict(case='moving capture vs off',file=file,**compare(ROOT/'controls/moving'/file,ROOT/'controls/moving-capture'/file)))
    bits.append(dict(case='selected initialization; physical fields retained',file='chk/EMS_000000.2d.hdf5',**compare(ROOT/'controls/baseline/chk/EMS_000000.2d.hdf5',ROOT/'controls/moving/chk/EMS_000000.2d.hdf5',skip=(16,17))))
    rows=[]
    for source in ('baseline','src'):
        for gauge in ('experimental','moving_puncture'):
            name=source+'-'+gauge;d=ROOT/'generic'/name;receipts.append(own(d,name))
            q=np.fromfile(d/'native.bin').reshape(-1,COL)
            if source=='src':rows.append(dict(case=name,**formulas(q,gauge)))
    old=np.fromfile(ROOT/'generic/baseline-experimental/native.bin').view('u8')
    new=np.fromfile(ROOT/'generic/src-experimental/native.bin').view('u8');assert np.array_equal(old,new)
    rows.append(dict(case='generic experimental pinned vs selected',native_samples=64,nonzero_momentum_samples=64,matter_Gamma_peak='',Gamma_sum_budget_fraction='',production_total_bit_mismatches=int(np.count_nonzero(old!=new)),gauge_formula_budget_fraction='',B_equation_KO_convention='all native columns bit-identical'))
    old=np.fromfile(ROOT/'generic/baseline-moving_puncture/native.bin').reshape(-1,COL)
    new=np.fromfile(ROOT/'generic/src-moving_puncture/native.bin').reshape(-1,COL)
    assert np.all(np.linalg.norm(new[:,4*N+4:4*N+6],axis=1)>0)
    change=new[:,N+16:N+18]-old[:,N+16:N+18];matter=new[:,4*N+2:4*N+4]
    assert np.all(np.abs(change-matter)<128*EPS*np.maximum(1,np.abs(change)+abs(matter)))
    # Exactly the gauge rows can differ between old buggy and repaired MPG.
    fields=[i for i in range(N) if i not in (13,14,15,16,17)]
    assert np.array_equal(old[:,N:2*N][:,fields].copy().view('u8'),new[:,N:2*N][:,fields].copy().view('u8'))
    save('t22-generic-equations.csv',rows);save('t22-default-identity.csv',bits)
    recorder=module('record',HERE/'t21-record.py');caprows=[]
    for case,gauge in [('capture','experimental'),('moving-capture','moving_puncture')]:
        d=ROOT/'controls'/case;checked=0;worst=0.
        for path in d.glob('*.bin'):
            for fr in recorder.frames(path):
                assert fr['meta']['gauge']==gauge
                v=fr['data']['values'];checked+=len(v)
                if gauge=='moving_puncture':
                    pre=v[:,N:2*N];u=v[:,:N];expected=pre[:,11:13]-u[:,16:18]
                    fraction=np.max(abs(pre[:,16:18]-expected)/(128*EPS*np.maximum(1,abs(pre[:,16:18])+abs(expected))))
                    worst=max(worst,float(fraction));assert fraction<=1
        assert checked>0
        caprows.append(dict(case=case,gauge=gauge,native_capture_cells=checked,B_complete_Gamma_budget_fraction=worst if gauge=='moving_puncture' else '',driver_comparison='within gauge only'))
    save('t22-native-capture.csv',caprows)
    restart=[]
    for name in ('restart-same','restart-legacy','restart-cross','restart-cross-tagged','restart-same-experimental'):
        d=ROOT/'controls'/name;receipts.append(own(d,name));p=json.loads((d/'probe.receipt.json').read_text());assert p['pass_control'];restart.append(dict(case=name,**p))
    save('t22-restart-controls.csv',restart);save('t22-quick-resources.csv',receipts)
    summary=dict(status='QUICK_AUDITS_PASS',bit_comparisons=sum(a['Float64_values'] for a in bits)+old.size,bit_mismatches=0,generic_states=64,matter_source_present=True,all_six_coefficients_exercised=True,capture_labels_verified=True,restart_compatibility_verified=True,peak_RSS_bytes=max(a['peak_RSS_bytes'] for a in receipts),scope='two-level production initialization and one coarse RK4 step; generic native equations; not full E-high equilibrium or launch admission')
    (HERE/'t22-quick-qualification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
def equilibrium():
    receipts=[];rows=[];summaries=[];states=[]
    trace=json.loads((HERE/'t22-trace-qualification.json').read_text())
    assert trace['status']=='K_ROUNDOFF_THETA_EXACT_ZERO'
    trace_receipt=json.loads((ROOT/'trace-audit-job/trace-audit.resources.json').read_text())
    assert trace_receipt['returncode']==trace_receipt['child_measurement']['returncode']==0 and not trace_receipt['gate_reason']
    measured=list(csv.DictReader((HERE/'t22-KTheta-census.csv').open()))
    lapse_rates=list(csv.DictReader((HERE/'t22-initial-lapse-rates.csv').open()))
    expected=list(csv.DictReader(Path('/Users/auroradysis/Workspace/EMS/.data/exp-0023/submission/phase1-evidence/smoke/census/E-high-dense/box-census.csv').open()))
    expected_boxes={l:sorted(tuple(int(a[k]) for k in ('lo_i','lo_j','hi_i','hi_j')) for a in expected if int(a['level'])==l) for l in range(15)}
    for gauge in ('experimental','moving_puncture'):
        d=ROOT/'equilibrium'/gauge;receipts.append(own(d,'equilibrium-'+gauge));assert 'T22_NATIVE_INITIAL_AUDIT_COMPLETE' in (d/'run.log').read_text()
        for l in range(15):
            census=list(csv.DictReader((d/f't22-boxes-L{l}.csv').open()));actual=sorted(tuple(int(a[k]) for k in ('lo_i','lo_j','hi_i','hi_j')) for a in census)
            assert actual==expected_boxes[l],(gauge,l,'different hierarchy')
            ss=list(csv.DictReader((d/f't22-state-L{l}.csv').open()))[0];assert int(ss['nonfinite'])==0
            # Explicit correction of the former exact-zero assertion: the
            # setter contracts an analytically traceless tensor in Float64.
            # K is measured against 64 eps times its local contraction scale;
            # Theta still must be exact zero. The Stage D screen is untouched.
            kr=next(a for a in measured if a['gauge']==gauge and int(a['level'])==l and a['scope']=='valid')
            assert float(kr['Theta_peak'])==0 and float(kr['K_over_eps_contraction_scale_peak'])<=trace['declared_valid_bound_eps']
            assert int(ss['nonzero_KTheta_cells'])==int(kr['K_nonzero'])
            if gauge=='moving_puncture':assert int(ss['nonzero_B_cells'])==0
            states.append(dict(gauge=gauge,**ss,boxes=len(census),boxes_equal_collected_E_high=True,K_peak=kr['K_peak'],K_cell_RMS=kr['K_cell_RMS'],Theta_peak=kr['Theta_peak'],K_over_eps_contraction_scale_peak=kr['K_over_eps_contraction_scale_peak'],replaced_zero_assertion='Theta exact zero; K <=64 eps local contraction scale'))
        data=np.load(ROOT/'trace-audit'/(gauge+'-native-all-levels.npz'))['data']
        for l in range(15):
            subset=data[data[:,0]==l];meta=subset[:,:5];q=subset[:,5:]
            f=formulas(q,gauge)
            if gauge=='moving_puncture':
                assert np.all(q[:,N+14:N+16]==0),'initial pre-KO shift must vanish'
                expected=-2*q[:,13]*(q[:,5]-2*q[:,10])
                scale=abs(expected)+abs(q[:,N+13])
                assert np.all(abs(q[:,N+13]-expected)<=128*EPS*np.maximum(scale,np.finfo(float).tiny)),'measured MPG lapse formula'
                assert np.array_equal(q[:,N+16:N+18].copy().view('u8'),q[:,N+11:N+13].copy().view('u8'))
            lr=next(a for a in lapse_rates if a['gauge']==gauge and int(a['level'])==l)
            summaries.append(dict(gauge=gauge,level=l,**f,preKO_lapse_peak_all_valid=float(lr['preKO_lapse_peak']),preKO_lapse_cell_RMS_all_valid=float(lr['preKO_lapse_cell_RMS']),preKO_lapse_volume_RMS_all_valid=float(lr['preKO_lapse_volume_RMS']),preKO_lapse_peak=float(abs(q[:,N+13]).max()),preKO_lapse_RMS_sampled=float(np.sqrt(np.mean(q[:,N+13]**2))),preKO_shift_peak=float(abs(q[:,N+14:N+16]).max()),preKO_B_peak=float(abs(q[:,N+16:N+18]).max()),preKO_B_RMS_sampled=float(np.sqrt(np.mean(q[:,N+16:N+18]**2))),KO_lapse_peak=float(abs(q[:,2*N+13]).max()),KO_shift_peak=float(abs(q[:,2*N+14:2*N+16]).max()),KO_B_peak=float(abs(q[:,2*N+16:2*N+18]).max()),Gamma_geometric_peak=float(abs(q[:,4*N:4*N+2]).max()),Gamma_geometric_RMS_sampled=float(np.sqrt(np.mean(q[:,4*N:4*N+2]**2))),Gamma_matter_peak=float(abs(q[:,4*N+2:4*N+4]).max()),Gamma_matter_RMS_sampled=float(np.sqrt(np.mean(q[:,4*N+2:4*N+4]**2))),Gamma_full_preKO_peak=float(abs(q[:,N+11:N+13]).max()),Gamma_full_preKO_RMS_sampled=float(np.sqrt(np.mean(q[:,N+11:N+13]**2))),KO_Gamma_peak=float(abs(q[:,2*N+11:2*N+13]).max()),KO_Gamma_RMS_sampled=float(np.sqrt(np.mean(q[:,2*N+11:2*N+13]**2)))))
            for m,v in zip(meta,q):
                row=dict(gauge=gauge,level=int(m[0]),i=int(m[1]),j=int(m[2]),h_M=m[3],time_M=m[4])
                for name,offset in [('state',0),('preKO',N),('KO',2*N),('total',3*N)]:
                    for field in ('K','Theta','Gamma1','Gamma2','lapse','shift1','shift2','B1','B2'):
                        row[name+'_'+field]=v[offset+FIELDS.index(field)]
                for i in range(2):row['geometric_Gamma'+str(i+1)]=v[4*N+i];row['matter_Gamma'+str(i+1)]=v[4*N+2+i]
                rows.append(row)
    identity=compare(ROOT/'equilibrium/experimental/chk/EMS_000000.2d.hdf5',ROOT/'equilibrium/moving_puncture/chk/EMS_000000.2d.hdf5',skip=(16,17))
    save('t22-equilibrium-native.csv',rows);save('t22-equilibrium-summary.csv',summaries);save('t22-initial-state-census.csv',states);save('t22-equilibrium-resources.csv',receipts)
    result=dict(status='E_HIGH_NATIVE_EQUILIBRIUM_AUDIT_PASS_WITH_MEASURED_K_ROUNDOFF',levels=15,boxes_equal_collected_E_high=True,physical_initial_identity=identity,preKO_MPG_shift_exact_zero=True,preKO_MPG_lapse_exact_zero=False,preKO_MPG_lapse_measured_peak=trace['MPG_preKO_lapse_valid_peak'],K_measured_peak=trace['K_valid_peak'],K_bound_eps=trace['declared_valid_bound_eps'],K_scaled_eps_peak=trace['K_valid_scaled_eps_peak'],Theta_exact_zero=True,preKO_MPG_B_equals_complete_Gamma=True,KO_B_initial_zero=True,peak_RSS_bytes=max(a['peak_RSS_bytes'] for a in receipts),scope='all valid cells K/Theta/B and native lapse rates; native sampled puncture/ray Gamma and geometric/matter/KO split on all levels 0-14; no advances or positive-time static reads',exact_zero_premise_correction='consult premise that setter supplies exact K and Theta zeros was false for K: setter computes a Float64 trace; measured roundoff bound replaces former exact-zero assertion explicitly; scientific screen unchanged')
    (HERE/'t22-equilibrium-qualification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
def launch():
    # Stream only completed native states. This is a gauge-labelled T13 serial
    # audit, not the production T21 recorder. No native replay is needed to
    # obtain the signed evolved Gamma disturbance.
    dec=module('frames',HERE/'t7-check.py');reader=module('sampling',HERE/'t13-analyze.py')
    rays={'axis':np.array([1.,0.]),'diagonal':np.ones(2)/np.sqrt(2)}
    radius=np.unique(np.r_[np.linspace(.00075,.0025,257),.0015]);receipts=[];history=[];cache={};floorrows=[]
    def norm(a,kind):return float(np.max(abs(a))) if kind=='peak' else float(np.sqrt(np.trapezoid(a*a,radius)/(radius[-1]-radius[0])))
    def interp(lo,a,h,vec,n):
        xy=radius[:,None]*vec;z=np.c_[(336+xy[:,0])/h-.5,xy[:,1]/h-.5];start=np.floor(z).astype(int)-n//2+1
        wx=reader.weights(z[:,0],start[:,0],n);wy=reader.weights(z[:,1],start[:,1],n);ans=np.zeros((len(radius),N));anchor=None
        odd=[2,7,12,15,17,22,23,25,26]
        for j in range(n):
            yy=start[:,1]+j
            for i in range(n):
                ix=start[:,0]+i-lo[0];iy=np.where(yy<0,-yy-1,yy)-lo[1]
                assert np.all(ix>=0) and np.all(ix<a.shape[1]) and np.all(iy>=0) and np.all(iy<a.shape[0]),'native frame lacks I10 stencil'
                v=a[iy,ix].copy();assert np.isfinite(v).all()
                v[np.ix_(np.flatnonzero(yy<0),odd)]*=-1
                if anchor is None:anchor=v.copy();ans=anchor.copy()
                ans+=(v-anchor)*(wx[i]*wy[j])[:,None]
        return ans
    for gauge in ('experimental','moving_puncture'):
        for half in (False,True):
            name=gauge+('-half' if half else '');d=ROOT/'launch'/name;receipts.append(own(d,name))
            assert 'T13 clean native stop' in (d/'run.log').read_text() and 'GRChombo finished.' in (d/'run.log').read_text()
            actual=float((d/'t13-stop.csv').read_text().splitlines()[1]);target=56*(7/12/4096/4)
            assert abs(actual-target)<64*EPS,'wrong physical launch stop'
            assert not list(d.rglob('*.hdf5'))
            for p in d.glob('t13-floors-L*.csv'):
                rr=list(csv.DictReader(p.open()));assert rr and all(int(r['nonfinite'])==0 for r in rr)
                row=dict(run=name,file=p.name,rows=len(rr),chi_activations=sum(int(r['chi_activations']) for r in rr),lapse_activations=sum(int(r['lapse_activations']) for r in rr),nonfinite=0,chi_min=min(float(r['chi_min']) for r in rr),lapse_min=min(float(r['lapse_min']) for r in rr));floorrows.append(row)
            group=[];previous=None;times=[];profiles={ray:{n:[] for n in (8,10)} for ray in rays}
            def consume(group):
                lo,a=reader.dense(group);h=group[0]['meta'][1];times.append(group[0]['meta'][0])
                for ray,vec in rays.items():
                    for n in (8,10):
                        sampled=interp(lo,a,h,vec,n);profiles[ray][n].append(sampled[:,11:13]@vec)
            for fr in dec.frames(d/'t13-t7-stage-L12.xz'):
                if fr['phase']!=50:continue
                t=fr['meta'][0]
                if previous is not None and t!=previous:consume(group);group=[]
                group.append(fr);previous=t
            if group:consume(group)
            assert len(times)==(113 if half else 57),len(times)
            for ray in rays:
                a=np.asarray(profiles[ray][8]);b=np.asarray(profiles[ray][10]);delta=a-a[0];err=abs(a-b)+abs(a[0]-b[0]);floor=128*EPS*np.maximum(1,abs(a)+abs(a[0]))
                cache[name,ray]=dict(times=np.array(times),delta=delta,err=err,floor=floor)
                for i,t in enumerate(times):
                    for kind in ('peak','RMS'):
                        history.append(dict(run=name,gauge=gauge,ray=ray,step=i,time_M=t,norm=kind,Gamma_disturbance=norm(delta[i],kind),interpolation_spread=norm(err[i],kind),Float64_floor=norm(floor[i],kind),sampling_rule='I8 central and abs(I8-I10); fivefold',B_compared=False))
    ratios=[]
    for ray in rays:
        old=cache['experimental',ray];new=cache['moving_puncture',ray]
        old_half=cache['experimental-half',ray];new_half=cache['moving_puncture-half',ray]
        assert np.allclose(old['times'],new['times'],atol=64*EPS,rtol=0)
        assert np.allclose(old['times'],old_half['times'][::2],atol=64*EPS,rtol=0)
        for kind in ('peak','RMS'):
            for i,t in enumerate(old['times']):
                aa=norm(old['delta'][i],kind);bb=norm(new['delta'][i],kind)
                eo=norm(old['err'][i]+old['floor'][i],kind);en=norm(new['err'][i]+new['floor'][i],kind)
                to=norm(old['delta'][i]-old_half['delta'][2*i],kind);tn=norm(new['delta'][i]-new_half['delta'][2*i],kind)
                e=5*(eo+en)+to+tn;qualified=t>0 and aa>5*eo+to and bb>5*en+tn
                ratios.append(dict(ray=ray,norm=kind,step=i,time_M=t,old_Gamma=aa,new_Gamma=bb,new_over_old=bb/aa if aa else '',old_fivefold_sampling_margin=5*eo,new_fivefold_sampling_margin=5*en,old_temporal_difference=to,new_temporal_difference=tn,fivefold_joint_sampling_margin=5*(eo+en),temporal_joint_margin=to+tn,reduction=aa-bb,reduction_lower=aa-bb-e,reduction_upper=aa-bb+e,ratio_lower=max(0,bb-5*en-tn)/(aa+5*eo+to) if aa else '',ratio_upper=(bb+5*en+tn)/(aa-5*eo-to) if aa>5*eo+to else '',old_significance_margin=aa-5*eo-to,new_significance_margin=bb-5*en-tn,qualified=qualified,status='qualified' if qualified else 'uncertainty_or_initial_floor'))
    save('t22-launch-history.csv',history);save('t22-launch-ratios.csv',ratios);save('t22-launch-floors.csv',floorrows);save('t22-launch-resources.csv',receipts)
    endpoint=[r for r in ratios if r['step']==56];save('t22-launch-endpoint.csv',endpoint)
    result=dict(status='LAUNCH_PAIR_COMPLETE',endpoint=endpoint,clock_count=57,native_stop_M=target,peak_RSS_bytes=max(a['peak_RSS_bytes'] for a in receipts),any_floor_activations_including_ghosts=any(r['chi_activations'] or r['lapse_activations'] for r in floorrows),nonfinite=0,scope='E-mid levels 0-12, 0.002 M native window, both gauges plus temporal controls; not a 10.5 M verdict')
    (HERE/'t22-launch-qualification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
    if sys.argv[1]=='restart-probe':restart_probe(*sys.argv[2:])
    else:globals()[sys.argv[1].replace('-','_')]()
