"""C++/Julia transfer through the consumed fourth-order Cartesian stencils.
Pre-data allowance: 5e-13 value-transfer scale, amplified by each stencil's
absolute coefficient sum / h^order. This tests transfer, not FD truncation.
"""
from pathlib import Path
from evidence import RAW, OUT, reference
import json,os,subprocess,time
import numpy as np
R=Path(__file__).resolve().parent;out=R/'results'
exe=str(next(R.glob('EMSCTT2Readback*.ex')))
args=[exe,str(out/'input-snapshots/a0.9-e7.trumpet'),str(out/'input-snapshots/n24-final.ctt')]
start=time.monotonic()
p=subprocess.run(args,input=(reference('stencil-points.txt')).read_text(),text=True,capture_output=True,check=True,timeout=240,env={**os.environ,'OMP_NUM_THREADS':'1'})
rows=[]
for line in p.stdout.splitlines():
    try:r=[float(x) for x in line.split()]
    except ValueError:continue
    if len(r)==29:rows.append(r)
cpp=np.array(rows);ref=np.loadtxt(reference('stencil-julia.txt'));sc=np.loadtxt(reference('jet-scales.txt'))
assert cpp.shape==ref.shape and cpp.shape==(len(sc)*10,29)
cpp=cpp.reshape(len(sc),2,5,29)*sc[:,None,None,:]
ref=ref.reshape(len(sc),2,5,29)*sc[:,None,None,:]
value_error=abs(cpp-ref)/(1+abs(ref));assert value_error.max()<=5e-13
h=np.array([.002,.001])[None,:,None]
summary={'jets':len(sc),'stencil_values':len(rows),'max_value_error':float(value_error.max()),'directions':'x,y,z,(x+y)/sqrt(2),(x+z)/sqrt(2),(y+z)/sqrt(2)','selection':'automatic crossing plus both explicit panel sides'}
for order,weights,divisor,ids in [(1,[1,-8,0,8,-1],12,range(29)),(2,[-1,16,-30,16,-1],12,[*range(9),18,19,26,27,28])]:
    ids=list(ids);w=np.array(weights)[None,None,:,None]
    a=(cpp*w).sum(axis=2)/(divisor*h**order)
    b=(ref*w).sum(axis=2)/(divisor*h**order)
    error=abs(a[:,:,ids]-b[:,:,ids])/(1+abs(b[:,:,ids]))
    # Each value has its own scale; propagate the retained 5e-13 allowance.
    budget=5e-13*((1+abs(ref))*abs(w)).sum(axis=2)/(divisor*h**order)
    budget=budget[:,:,ids]/(1+abs(b[:,:,ids]))
    assert np.all(np.isfinite(error)) and np.all(error<=budget),(order,error.max(),(error/budget).max())
    summary['D'+str(order)]={'max_error':float(error.max()),'max_allowance':float(budget.max()),'max_fraction_of_allowance':float((error/budget).max())}
summary['seconds']=time.monotonic()-start
(OUT/'stencil-readback.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
