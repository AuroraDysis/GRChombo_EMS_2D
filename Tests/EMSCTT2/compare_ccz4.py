from pathlib import Path
import subprocess,os,json
import numpy as np
R=Path(__file__).resolve().parent;out=R/'results';data=out/'input-snapshots';exe=str(next(R.glob('EMSCTT2Readback*.ex')))
summary={}
for name in ('zero','n24'):
    p=subprocess.run([exe,str(data/'a0.9-e7.trumpet'),str(data/(name+'-final.ctt')),'--ccz4'],input=(out/f'ccz4-{name}-points.txt').read_text(),capture_output=True,text=True,check=True,timeout=240,env={**os.environ,'OMP_NUM_THREADS':'1'})
    a=np.array([[float(v) for v in line.split()] for line in p.stdout.splitlines() if len(line.split())==21])
    b=np.loadtxt(out/f'ccz4-{name}-julia.txt');sc=np.loadtxt(out/f'ccz4-{name}-scales.txt')
    assert a.shape==b.shape==(5*len(sc),21)
    a=a.reshape(-1,5,21)*sc[:,None,:];b=b.reshape(-1,5,21)*sc[:,None,:]
    err=abs(a-b)/(1+abs(b));assert err.max()<5e-13,(name,err.max())
    result={'values':int(a.shape[0]*5),'max_value_error':float(err.max())}
    for order,w,ids in [(1,[1,-8,0,8,-1],range(21)),(2,[-1,16,-30,16,-1],[0,1,2,3,4,10,11,18,19,20])]:
        w=np.array(w)[None,:,None];ids=list(ids)
        f=(a*w).sum(axis=1)/(12*.001**order);g=(b*w).sum(axis=1)/(12*.001**order)
        err=abs(f[:,ids]-g[:,ids])/(1+abs(g[:,ids]))
        budget=5e-13*((1+abs(b))*abs(w)).sum(axis=1)/(12*.001**order)
        budget=budget[:,ids]/(1+abs(g[:,ids]))
        assert np.all(np.isfinite(err)) and np.all(err<=budget)
        result['D'+str(order)]={'max_error':float(err.max()),'max_allowance':float(budget.max()),'max_fraction_of_allowance':float((err/budget).max())}
    summary[name]=result
(out/'ccz4-readback.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
