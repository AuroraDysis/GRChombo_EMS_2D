#!/usr/bin/env python3
"""Stored numerical initialization census, scoped traceless witness and native rates."""
import csv, hashlib, importlib.util, json, math, resource, subprocess, sys
from pathlib import Path
import h5py, numpy as np, sympy as s, mpmath as mp
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t22');OUT=ROOT/'trace-audit';OUT.mkdir(exist_ok=True)
EPS=np.finfo(float).eps;BUDGET_EPS=64
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(name,rows):
    with (HERE/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def own(d,name):
    q=json.loads((d/(name+'.resources.json')).read_text());assert q['returncode']==q['child_measurement']['returncode']==0 and not q['gate_reason'] and q['peak_rss_bytes']<6e9
    assert (d/'done.exit').read_text().strip()=='0';return q
def witness():
    # Q[P,k,nx,ny,nz], P !=0, exact equality modulo n.n=1.
    # At zero rapidity the Lorentz matrix is identity, gamma=P^2 I and
    # Kij=P^2 k(3 ni nj-deltaij). The smallest witness is the trace remainder.
    P,k,nx,ny,nz=s.symbols('P k nx ny nz',real=True)
    matrix=P**2*k*(3*s.Matrix([nx,ny,nz])*s.Matrix([[nx,ny,nz]])-s.eye(3))
    trace=s.cancel(s.trace(matrix/P**2));residual=s.Poly(trace-3*k*(nx**2+ny**2+nz**2-1),P,k,nx,ny,nz).as_expr();assert residual==0
    mp.mp.dps=50;worst=mp.mpf(0)
    # Current-coordinate directions and arbitrary amplitudes, not static data.
    for j in range(1,33):
        x,y=mp.mpf(j)/17,mp.mpf(j+3)/29;r=mp.sqrt(x*x+y*y);n=mp.matrix([x/r,y/r,0]);p=mp.mpf(j+2)/11;kk=mp.mpf(j)/13
        gamma=p*p*mp.eye(3);K=p*p*kk*(3*n*n.T-mp.eye(3));value=mp.fsum((gamma**-1)[i,z]*K[i,z] for i in range(3) for z in range(3))
        worst=max(worst,abs(value)/max(1,abs(kk)))
    assert worst<mp.mpf('1e-47')
    return dict(exact_status='PROVED',domain='P != 0, real k and unit real direction n; zero rapidity; no companion',witness=str(residual),identity='gamma^ij Kij = 3 k (n.n-1) = 0',numerical_status='CORROBORATED',samples=32,precision_digits=50,worst_scaled_residual=str(worst),static_reference_used=False)
def compile_rates():
    commands=json.loads((ROOT/'build-commands.json').read_text());cmd=next(c for c in commands if any(x.endswith('/src/Tests/EMSNative/t22-audit.cpp') for x in c))
    cmd=[str(HERE/'t22-initial-gauge-rates.cpp') if x.endswith('/src/Tests/EMSNative/t22-audit.cpp') else x for x in cmd];cmd[-1]=str(OUT/'rates.ex')
    (OUT/'build-command.json').write_text(json.dumps(cmd,indent=2)+'\n');subprocess.run(cmd,check=True)
def stats(a,weights):
    return dict(peak=float(np.max(abs(a))),cell_RMS=float(np.sqrt(np.mean(a*a))),volume_RMS=float(np.sqrt(np.sum(weights*a*a)/np.sum(weights))),nonzero=int(np.count_nonzero(a)))
def analyze():
    witness_out=witness();rows=[];gauge_rows=[];native={};receipts=[];aliases=[];storage=[]
    for gauge in ('experimental','moving_puncture'):
        d=ROOT/'equilibrium'/gauge;receipts.append(own(d,'equilibrium-'+gauge));datafile=d/'chk/EMS_000000.2d.hdf5'
        with h5py.File(datafile) as f:
            assert float(f.attrs['time'])==0.;assert int(f.attrs['num_levels'])==15
            input_file=OUT/(gauge+'-blocks.bin');output_file=OUT/(gauge+'-rates.bin');blocks=[]
            with input_file.open('wb') as stream:
                for l in range(15):
                    group=f[f'level_{l}'];h=float(group.attrs['dx']);offset=group['data:offsets=0'][:];values=group['data:datatype=0'];centre=336.;c=round(centre/h)
                    points={(c+i,j) for j in range(4) for i in range(-4,4)}
                    for r in (.0004,.001,.003,.0055,.0060286129392158804,.0075,.012,.02):
                        for angle in (0.,np.pi/2,np.pi/4):points.add((round((centre+r*np.cos(angle))/h-.5),max(0,round(r*np.sin(angle)/h-.5))))
                    accum={scope:dict(K=[],Theta=[],scale=[],x=[],y=[],weights=[]) for scope in ('valid','stored_ghost_copies','puncture_4_cells')}
                    for box,b in enumerate(group['boxes'][:]):
                        x0,y0,x1,y1=[int(b[key]) for key in ('lo_i','lo_j','hi_i','hi_j')];nx=x1-x0+1;ny=y1-y0+1
                        a=values[int(offset[box]):int(offset[box+1])].reshape(28,ny+6,nx+6)
                        if gauge=='moving_puncture':assert np.all(a[16:18]==0),'MPG driver must be zero on valid cells and stored ghosts'
                        p=np.array(sorted([(i,j) for i,j in points if x0<=i<=x1 and y0<=j<=y1],key=lambda z:(z[1],z[0])),dtype='i4').reshape(-1,2)
                        info=np.array([l,box,h,x0,y0,x1,y1,centre,len(p)],dtype='f8');info.tofile(stream);a.tofile(stream);p.tofile(stream)
                        xx,yy=np.meshgrid(np.arange(x0-3,x1+4),np.arange(y0-3,y1+4));x=(xx+.5)*h-centre;y=(yy+.5)*h
                        valid=(xx>=x0)&(xx<=x1)&(yy>=y0)&(yy<=y1);puncture=valid&(abs(xx-c+.5)<=4)&(yy<4)&(yy>=0)
                        det=a[1]*a[3]-a[2]*a[2]
                        scale=abs(a[3]*a[6]/det)+2*abs(a[2]*a[7]/det)+abs(a[1]*a[8]/det)+abs(a[9]/a[4])
                        for scope,take in [('valid',valid),('stored_ghost_copies',~valid),('puncture_4_cells',puncture)]:
                            for key,value in [('K',a[5]),('Theta',a[10]),('scale',scale),('x',x),('y',y),('weights',abs(y))]:accum[scope][key].append(value[take])
                        blocks.append(dict(info=info,points=p,lapse=a[13,3:-3,3:-3].copy(),K=a[5,3:-3,3:-3].copy(),Theta=a[10,3:-3,3:-3].copy(),weights=abs(y[3:-3,3:-3]).copy()))
                    for scope,raw in accum.items():
                        arrays={k:np.concatenate(v) for k,v in raw.items()};kk=arrays['K'];tt=arrays['Theta'];scale=arrays['scale'];w=arrays['weights'];assert len(kk)>0
                        ratio=abs(kk)/(EPS*np.maximum(scale,np.finfo(float).tiny));loc=int(np.argmax(abs(kk)));ks=stats(kk,w);ts=stats(tt,w)
                        row=dict(gauge=gauge,level=l,scope=scope,cells=len(kk),K_peak=ks['peak'],K_cell_RMS=ks['cell_RMS'],K_volume_RMS=ks['volume_RMS'],K_nonzero=ks['nonzero'],Theta_peak=ts['peak'],Theta_cell_RMS=ts['cell_RMS'],Theta_volume_RMS=ts['volume_RMS'],Theta_nonzero=ts['nonzero'],K_max_x_M=arrays['x'][loc],K_max_y_M=arrays['y'][loc],K_max_radius_M=np.hypot(arrays['x'][loc],arrays['y'][loc]),K_over_eps_contraction_scale_peak=float(ratio.max()),declared_bound_eps=BUDGET_EPS,roundoff_bound_pass=bool(np.all(ratio<=BUDGET_EPS)),Theta_exact_zero=bool(np.all(tt==0)))
                        rows.append(row)
                        # Ghost K may be a point-interpolated rounding residual.
                        # This check is measured and reported separately.
                        assert np.all(tt==0),row
                        if scope=='valid':assert np.all(ratio<=BUDGET_EPS),row
            subprocess.run([str(OUT/'rates.ex'),gauge,str(input_file),str(output_file)],check=True)
            local=[];rates_by_level={l:[] for l in range(15)};weights_by_level={l:[] for l in range(15)};formula_bound=0.
            with output_file.open('rb') as stream:
                for block in blocks:
                    info=np.frombuffer(stream.read(9*8),dtype='f8');assert np.array_equal(info,block['info']);l=int(info[0]);nx=int(info[5]-info[3]+1);ny=int(info[6]-info[4]+1)
                    rr=np.frombuffer(stream.read(nx*ny*2*8),dtype='f8').reshape(ny,nx,2).copy();rates_by_level[l].append(rr[...,1].ravel());weights_by_level[l].append(block['weights'].ravel())
                    if gauge=='moving_puncture':
                        expected=-2*block['lapse']*(block['K']-2*block['Theta']);scale=abs(expected)+abs(rr[...,1]);frac=abs(expected-rr[...,1])/(128*EPS*np.maximum(scale,np.finfo(float).tiny));formula_bound=max(formula_bound,float(frac.max()));assert np.all(frac<=1)
                    else:
                        expected=rr[...,0]-1.8*block['lapse']*(block['K']-2*block['Theta']);assert np.allclose(rr[...,1],expected,rtol=128*EPS,atol=np.finfo(float).tiny)
                    for i,j in block['points']:
                        q=np.frombuffer(stream.read(124*8),dtype='f8').copy();assert len(q)==124 and np.isfinite(q).all();local.append(np.r_[l,i,j,info[2],0.,q])
                assert not stream.read(1)
            native[gauge]=np.array(local)
            alias=list(csv.DictReader(Path(str(output_file)+'.tensor-aliases.csv').open()))
            assert all(int(a['component'])==7 for a in alias)
            aliases.extend(dict(gauge=gauge,field='A12',**a) for a in alias)
            # Native enum_mapping stores the lower off-diagonal component last.
            # Its overwritten upper alias is not an additional evolved field.
            assert np.all(native[gauge][:,-1]==0)
            assert not list(csv.DictReader(Path(str(output_file)+'.mismatches.csv').open()))
            storage.append(dict(gauge=gauge,native_points=len(local),stored_Float64_components=len(local)*28,stored_bit_mismatches=0,overwritten_tensor_alias_discrepancies=len(alias),alias_field='A12',comparison='native last-write storage against reconstructed last-write storage; bitwise, no tolerance'))
            for l in range(15):
                rate=np.concatenate(rates_by_level[l]);w=np.concatenate(weights_by_level[l]);st=stats(rate,w)
                kr=next(a for a in rows if a['gauge']==gauge and a['level']==l and a['scope']=='valid')
                gauge_rows.append(dict(gauge=gauge,level=l,valid_cells=len(rate),preKO_lapse_peak=st['peak'],preKO_lapse_cell_RMS=st['cell_RMS'],preKO_lapse_volume_RMS=st['volume_RMS'],preKO_lapse_nonzero_cells=st['nonzero'],K_peak=kr['K_peak'],Theta_peak=kr['Theta_peak'],lapse_formula='-2 alpha (K-2Theta)' if gauge=='moving_puncture' else 'native upwind beta.grad(alpha)-1.8 alpha (K-2Theta)',native_formula_budget_fraction=formula_bound if gauge=='moving_puncture' else '',source='unchanged gauge class + native FourthOrderDerivatives on stored t0 checkpoint ghosts'))
            # Large regenerable input is removed only after native success.
            digest=sha(input_file);input_file.unlink();(OUT/(gauge+'-blocks.sha256')).write_text(digest+'  '+str(input_file)+' (regenerable input; removed after success)\n')
    save('t22-KTheta-census.csv',rows);save('t22-initial-lapse-rates.csv',gauge_rows)
    save('t22-tensor-aliases.csv',aliases);save('t22-native-storage-check.csv',storage)
    for gauge,data in native.items():
        np.savez_compressed(OUT/(gauge+'-native-all-levels.npz'),data=data)
    valid=[r for r in rows if r['scope']=='valid'];moving=[r for r in gauge_rows if r['gauge']=='moving_puncture'];ghost=[r for r in rows if r['scope']=='stored_ghost_copies']
    result=dict(status='K_ROUNDOFF_THETA_EXACT_ZERO',K_valid_peak=max(r['K_peak'] for r in valid),K_valid_scaled_eps_peak=max(r['K_over_eps_contraction_scale_peak'] for r in valid),K_ghost_peak=max(r['K_peak'] for r in ghost),K_ghost_scaled_eps_peak=max(r['K_over_eps_contraction_scale_peak'] for r in ghost),declared_valid_bound_eps=BUDGET_EPS,Theta_peak=0.,MPG_preKO_lapse_valid_peak=max(r['preKO_lapse_peak'] for r in moving),setter_K_literal_zero=False,setter_K_source='EMSBH_trumpet_read.impl.hpp:419-426; full metric trace',setter_Theta_source='VarsTools::assign(vars,0) at line424; no overwrite',exact_zero_premise_false=True,scientific_screen_changed=False,witness=witness_out,process_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,completed_native_audit_peak_RSS_bytes=max(r['peak_rss_bytes'] for r in receipts))
    (HERE/'t22-trace-qualification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':globals()[sys.argv[1].replace('-','_')]()
