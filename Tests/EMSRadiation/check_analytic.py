"""Production surface extraction checks and postprocessor regression."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('rad',ROOT/'Examples/EMS/ems_radiation_postprocess.py')
rad=importlib.util.module_from_spec(spec); spec.loader.exec_module(rad)
OUT=Path(sys.argv[1]).resolve()
EXE=next((ROOT/'Tests/EMSRadiation').glob('EMSRadiationTest2d.*.ex'))
records=[]


def run(name,n,kind,mode=0):
    path=OUT/name; path.mkdir(parents=True,exist_ok=True)
    with (path/'run.log').open('w') as log:
        subprocess.run([str(EXE),str(path),str(n),str(kind),str(mode)],
                       stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120,
                       env=dict(os.environ,OMP_NUM_THREADS='1'))
    return np.genfromtxt(path/'ems_radiation.csv',delimiter=',',names=True)


for n in (33,65,129):
    a=run(f'pulse-{n}',n,0)
    d=rad.load(OUT/f'pulse-{n}/ems_radiation.csv')[20]
    v=rad.prepare(d,0)
    lum=rad.luminosities(v['news'],v['ds'],v['P'],1)[:,1]
    analytic=2*(.03*(v['t']-4)/.5**2*np.exp(-(v['t']-4)**2/(2*.5**2)))**2
    energy=rad.cumulative(lum,v['u'])[-1]
    expected=.03**2*np.sqrt(np.pi)/.5
    stress_energy=rad.cumulative(v['L_scalar_stress'],v['t'])[-1]
    error=abs(energy/expected-1)
    assert min(lum)>=0 and error<3e-5
    assert abs(stress_energy/expected-1)<1e-6
    records.append(dict(test='pulse',theta=n,expected=expected,measured=energy,
                        relative_error=error,peak_absolute_error=max(abs(lum-analytic)),
                        stress_energy=stress_energy))
for l in range(7):
    errors=[]
    for n in (65,129,257):
        a=run(f'mode-{l}-{n}',n,1,l)
        target=np.zeros(7); target[l]=.03
        err=max(abs(a['S_re']-target)); errors.append(err)
        assert err<1e-4
        records.append(dict(test='mode',l=l,theta=n,expected=.03,measured=a['S_re'][l],
                            absolute_error=err))
    assert errors[-1]<1e-7 and errors[-2]/errors[-1]>8
    assert max(abs(a['D_re']-target*.07/.03))<3e-7
for kind,name in ((2,'outgoing-EM'),(3,'incoming-EM')):
    a=run(name,129,kind)
    L=2*np.exp(-.4)*sum(a['P_re']**2+a['P_im']**2)
    expected=16*np.pi/3*np.exp(-.4)*.02**2
    assert abs(a['L_EM_stress'][0]/expected-(1 if kind==2 else -1))<1e-6
    assert abs(L/expected-(1 if kind==2 else 0))<1e-6
    records.append(dict(test=name,expected=expected if kind==2 else 0,measured=L,
                        stress=a['L_EM_stress'][0]))

# Fixed-u polynomial data: distinguish fixed-time fitting and preserve a
# nonzero scalar monopole. Constant news and scalar derivative are not filtered.
series=[]
for r in (50.,75.,100.,150.):
    u=np.linspace(-20,20,801); t=u+r
    base=np.exp(-u*u/16)
    d={'u':u,'R':np.full_like(u,r),'F_inf':np.ones_like(u)}
    for name in ('news','ds','P'):
        a=np.zeros((len(u),7),complex)
        l={'news':2,'ds':0,'P':1}[name]
        a[:,l]=base*(1+3/r+7/r**2)*(1+.2j)
        d[name]=a
    series.append(d)
u,quad=rad.extrapolate(series,2,.05)
expected=np.exp(-2*u*u/16)*1.04
assert max(abs(quad[:,1]-2*expected))<2e-12
u,linear=rad.extrapolate(series,1,.05)
assert max(abs(linear[:,1]-quad[:,1]))>1e-5
linear_only=[]
for r,d in zip((50.,75.,100.,150.),series):
    new=dict(d)
    for key in ('news','ds','P'):
        new[key]=d[key]*(1+3/r)/(1+3/r+7/r**2)
    linear_only.append(new)
_,linear_exact=rad.extrapolate(linear_only,1,.05)
assert max(abs(linear_exact[:,1]-2*expected))<2e-12
# Nonuniform retarded time and growing areal radius: D alone is insufficient.
t=np.linspace(0,3,101)**1.2; R=20+.01*t
d={'t':t,'R':R,'alpha2':np.ones_like(t),'F_inf':np.ones_like(t),
   'D':np.zeros((len(t),7),complex),'S':np.zeros((len(t),7),complex),
   'W':np.zeros((len(t),7),complex)}
d['S'][:,0]=R*.1
v=rad.prepare(d,0)
assert max(abs(v['ds'][:,0]-.001/.99))<1e-12
records.append(dict(test='postprocessor',quadratic_max_absolute_error=float(max(abs(quad[:,1]-2*expected))),
                    linear_max_absolute_error=float(max(abs(linear[:,1]-2*expected)))))
# Restart rewind discards the superseded future; incomplete blocks fail.
lines=(OUT/'pulse-129/ems_radiation.csv').read_text().splitlines()
restart=OUT/'restart.csv'
restart.write_text('\n'.join(lines+lines[1+400*7:1+411*7])+'\n')
assert abs(rad.load(restart)[20]['t'][-1]-4.1)<1e-12
restart.write_text('\n'.join(lines[:-1])+'\n')
try:
    rad.load(restart)
except ValueError:
    pass
else:
    raise AssertionError('truncated block accepted')
restart.unlink()
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'analytic.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2))
