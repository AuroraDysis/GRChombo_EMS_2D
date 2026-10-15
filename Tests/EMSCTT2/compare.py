"""Independent readback, zero identity, and consumed-order transfer controls."""
from pathlib import Path
import subprocess,json,math,time
import numpy as np
R=Path(__file__).resolve().parent
out=R/'results';data=out/'input-snapshots'
profile=str(out/'input-snapshots/a0.9-e7.trumpet')
exe=str(next(R.glob('EMSCTT2Readback*.ex')))
def run(name,points):
    t=time.monotonic()
    p=subprocess.run([exe,profile,str(data/(name+'-final.ctt')) if name!='plain' else 'plain'],input=points,text=True,capture_output=True,check=False,timeout=240,env={**__import__('os').environ,'OMP_NUM_THREADS':'1'})
    assert p.returncode==0,(p.returncode,p.stdout[-300:],p.stderr)
    rows=[]
    for line in p.stdout.splitlines():
        try: row=[float(x) for x in line.split()]
        except ValueError: continue
        if len(row)==29: rows.append(row)
    return np.array(rows),time.monotonic()-t
points=(out/'points.txt').read_text();summary={}
# P^2 and P are fixed at the sample, obtained from the uncorrected setter.
plain,_=run('plain',points)
X=1/np.sqrt(plain[:,0]);scale=np.ones_like(plain);scale[:,:18]=X[:,None]**2;scale[:,20:26]=X[:,None]
for name in ['zero','half','n12','n18','n24']:
    cpp,seconds=run(name,points);ref=np.loadtxt(out/(name+'-julia.txt'))
    err=np.abs((cpp-ref)*scale)/(1+np.abs(ref*scale))
    assert cpp.shape==ref.shape and np.all(np.isfinite(err))
    # 5e-13: retained trumpet/CTT Float64 transfer allowance.
    assert err.max()<=5e-13,(name,err.max(),np.unravel_index(err.argmax(),err.shape))
    summary[name]={'points':len(cpp),'max_scaled_value_error':float(err.max()),'seconds':seconds}
    if name=='zero':
        identity=np.abs((cpp-plain)*scale)/(1+np.abs(plain*scale))
        assert identity.max()<=5e-13
        summary[name]['plain_identity']=float(identity.max())
    np.savetxt(out/(name+'-cpp.txt'),cpp)
# 5-point Cartesian differences of independently executed C++ values, compared
# with Julia AD jets, and h/2 differences for an independent truncation estimate.
jets=np.loadtxt(out/'jet-points.txt');jref=np.loadtxt(out/'jets-julia.txt')
queries=[]
for row in jets:
    x=row[:3];r=np.linalg.norm(x);panel=int(row[3]);axis=int(row[4])-1
    direction=np.zeros(3)
    if axis<3: direction[axis]=1
    else:
        i,j=[(0,1),(0,2),(1,2)][axis-3]
        direction[i]=direction[j]=1/math.sqrt(2)
    for h in (.002,.001):
        for k in (-2,-1,0,1,2):
            y=x+k*h*r*direction
            queries.append(' '.join(map(str,[*y,panel])))
a,seconds=run('n24','\n'.join(queries)+'\n');a=a.reshape(len(jets),2,5,29)
sc=jref[:,:29];v=jref[:,29:58];d=jref[:,58:87];dd=jref[:,87:116]
a*=sc[:,None,None,:]
d1=(a[:,:,0]-8*a[:,:,1]+8*a[:,:,3]-a[:,:,4])/(12*np.array([.002,.001])[None,:,None])
d2=(-a[:,:,4]+16*a[:,:,3]-30*a[:,:,2]+16*a[:,:,1]-a[:,:,0])/(12*np.array([.002,.001])[None,:,None]**2)
# D1 all physical fields. D2 gamma and scalars (Pi is consumed by tagging).
# Tolerance fixed by measured independent Richardson estimate plus Float64
# differentiation conditioning: 128 eps/h^j, normalized by 1+|AD jet|.
report={}
for order,fd,ref,ids in [(1,d1,d,range(29)),(2,d2,dd,[*range(9),18,19,26,27,28])]:
    ids=list(ids);err=np.abs(fd[:,1,ids]-ref[:,ids])/(1+np.abs(ref[:,ids]))
    trunc=np.abs(fd[:,0,ids]-fd[:,1,ids])/(15*(1+np.abs(ref[:,ids])))
    floor=128*np.finfo(float).eps/.001**order
    tolerance=4*trunc+floor
    assert np.all(np.isfinite(err)) and np.all(err<=tolerance),(order,err.max(),(err/tolerance).max(),np.unravel_index((err/tolerance).argmax(),err.shape))
    report['D'+str(order)]={'samples':len(jets),'max_error':float(err.max()),'roundoff_allowance':float(floor),'max_error_over_allowance':float((err/tolerance).max()),'max_allowance':float(tolerance.max())}
summary['jets']=dict(report,seconds=seconds)
(out/'readback.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
