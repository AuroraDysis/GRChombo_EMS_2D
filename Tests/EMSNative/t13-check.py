#!/usr/bin/env python3
"""Streaming native replay and RK/projection closure; no static-file access."""
import sys
sys.dont_write_bytecode=True
import csv,hashlib,importlib.util,subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;TMP=Path('/private/tmp/ems-t13')
spec=importlib.util.spec_from_file_location('frames',HERE/'t7-check.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
N=28;COL=180;EPS=np.finfo(float).eps

def closure(directory,level):
    directory=Path(directory);states={};parts={};before={};rhs={}
    count=0;bit_mismatch=0;sum_error=0.;update_error=0.;projection_mismatch=0;term_error=0.;samples=0
    raw_sum_fraction=0.;raw_gamma_fraction=0.;sum_abs_error=0.;gamma_abs_error=0.;frame_index=0
    cache=TMP/'analysis'/(directory.name+'-stages');cache.mkdir(parents=True,exist_ok=True)
    payload=TMP/'closure-input.bin';output=TMP/'closure-output.bin'
    for fr in mod.frames(directory/f't13-t7-stage-L{level}.xz'):
        p,src=fr['phase'],fr['source'];v=fr['values'];iv=fr['cells']
        if p in (2,3,14,15):states[src,p]=fr
        if p in (20,21):parts[src,p]=fr
        if p==22:
            pre=states[src,2];lookup={tuple(k):a for k,a in zip(pre['cells'],pre['values'])}
            records=np.empty((len(iv),3+49*N));records[:,0]=fr['meta'][1]
            records[:,1]=(iv[:,1]+.5)*fr['meta'][1];records[:,2]=0.
            for j in range(-3,4):
                for i in range(-3,4):
                    records[:,3+(j+3)*7+i+3::49]=np.array([lookup[(k[0]+i,k[1]+j)] for k in iv])
            records.tofile(payload)
            (cache/f'frame-{frame_index:03}-input.sha256').write_text(hashlib.sha256(records.tobytes()).hexdigest()+'  native replay payload\n')
            label=f"replay-{directory.name}-{samples}-{src}"
            subprocess.run([sys.executable,str(HERE/'t13-run.py'),'--measure',label,
                '--directory',str(TMP/'replay-resources'),'--',str(TMP/'replay.ex'),
                '--stage',str(payload),str(output)],check=True)
            q=np.fromfile(output).reshape(-1,COL);assert q.shape[0]==len(iv)
            actual=q[:,152:180];bit_mismatch+=np.count_nonzero(actual.view('u8')!=v.view('u8'));count+=v.size
            post=states[src,3];after={tuple(k):a for k,a in zip(post['cells'],post['values'])}
            expected=np.array([after[tuple(k)] for k in iv])
            projection_mismatch+=np.count_nonzero(expected.view('u8')!=q[:,:N].copy().view('u8'))
            a,b=parts[src,20]['values'],parts[src,21]['values']
            raw=records[:,3:].reshape(-1,N,49);amp=np.max(abs(raw),axis=2);h=fr['meta'][1];y=records[:,1]
            # Before-KO + solo KO versus two sequential directional additions.
            # Use uncollapsed stencil magnitudes, not a cancelled final RHS.
            direct=abs(v-(a+b));budget=128*EPS*(abs(a)+abs(b)+2*amp/h+np.finfo(float).tiny)
            raw_sum_fraction=max(raw_sum_fraction,float(np.max(direct/budget)))
            sum_abs_error=max(sum_abs_error,float(direct.max()))
            terms=q[:,112:138].reshape(-1,13,2);gamma_defect=abs(terms[:,:12].sum(1)-a[:,11:13])
            scale=amp[:,14:16].max(1)*(1/h**2+1/y**2)
            scale+=amp[:,13]*amp[:,6:10].max(1)*(1/h+1/y)*(amp[:,1:5].max(1)+amp[:,0]/raw[:,0].min(1))
            scale+=amp[:,13]*(amp[:,5]+32*np.pi*amp[:,19]*amp[:,18])/h
            scale+=amp[:,10]*amp[:,13]/h+amp[:,14:16].max(1)*amp[:,11:13].max(1)/h
            scale+=np.abs(terms[:,:12]).sum((1,2))+1
            gamma_budget=128*EPS*scale
            raw_gamma_fraction=max(raw_gamma_fraction,float(np.max(gamma_defect/gamma_budget[:,None])))
            gamma_abs_error=max(gamma_abs_error,float(gamma_defect.max()))
            sum_error=max(sum_error,float(np.max(abs(v-(a+b))/(EPS*np.maximum(1.,abs(a)+abs(b))))))
            sum_error=max(sum_error,float(np.max(abs(q[:,28:56]-a)/(EPS*np.maximum(1.,abs(a))))))
            sum_error=max(sum_error,float(np.max(abs(q[:,56:84]-b)/(EPS*np.maximum(1.,abs(b))))))
            term_error=max(term_error,float(q[:,143].max()));samples+=len(iv)
            rhs[src]={tuple(k):a for k,a in zip(iv,v)}
            np.savez_compressed(cache/f'frame-{frame_index:03}.npz',cells=iv,q=q,
                physical=a,KO=b,total=v,meta=fr['meta'],stage=fr['stage'],source=src,
                pre_cells=pre['cells'],pre_values=pre['values'])
            frame_index+=1
        if p==12:before[src]=fr
        if p==13:
            old=before[src];assert np.array_equal(old['cells'],iv)
            rv=np.array([rhs[src][tuple(k)] for k in iv]);dt=fr['meta'][7]
            predicted=old['values']+dt*rv
            update_error=max(update_error,float(np.max(abs(v-predicted)/(EPS*np.maximum(1.,abs(v)+abs(old['values'])+abs(dt*rv))))))
        if p==3:
            assert np.array_equal(states[src,15]['values'].view('u8'),v.view('u8'))
    payload.unlink(missing_ok=True);output.unlink(missing_ok=True)
    row=dict(directory=str(directory),level=level,native_samples=samples,
        RHS_Float64_values=count,RHS_bit_mismatches=int(bit_mismatch),
        projected_state_bit_mismatches=int(projection_mismatch),
        physical_plus_KO_max_eps=sum_error,RK_update_max_eps=update_error,
        Gamma_term_sum_max_eps=term_error,recorded_physical_plus_KO_max_absolute=sum_abs_error,
        recorded_sum_raw_stencil_budget_fraction=raw_sum_fraction,
        Gamma_term_sum_max_absolute=gamma_abs_error,Gamma_raw_stencil_budget_fraction=raw_gamma_fraction)
    assert count>0 and bit_mismatch==0 and projection_mismatch==0,row
    assert raw_sum_fraction<1 and raw_gamma_fraction<1 and update_error<8,row
    print(row,flush=True);return row

if __name__=='__main__':
    rows=[]
    # This replay pins the E coupling f2=-0.8; reference f2=-20 is covered by
    # the binary identity controls, not by this E-specific term audit.
    for name in ('E',):
        rows.append(closure(TMP/'controls'/f'{name}-point-capture',1))
    with (HERE/'t13-replay-controls.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
