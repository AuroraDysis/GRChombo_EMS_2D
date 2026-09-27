"""Assert local radiation gates and retain compact, reproducible evidence."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import re
import shutil
import sys
import subprocess

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
BASE=Path(sys.argv[1]).resolve()
OUT=ROOT/'Tests/EMSRadiation/results'
spec=importlib.util.spec_from_file_location('rad',ROOT/'Examples/EMS/ems_radiation_postprocess.py')
rad=importlib.util.module_from_spec(spec); spec.loader.exec_module(rad)
OUT.mkdir(exist_ok=True)
report={'pulse_grid':[],'stationary':[],'production_gauge_control':[],
        'threshold':[],'provenance':{},'passed':False}


def authenticate(path):
    report['provenance'][str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()


exact=.03**2*np.sqrt(np.pi)/.5
for n in (64,128,256):
    path=BASE/f'pulse-grid2/pulse-grid-{n}/ems_radiation.csv'; authenticate(path)
    for r,d in rad.load(path).items():
        v=rad.prepare(d,0); L=rad.luminosities(v['news'],v['ds'],v['P'],1)[:,1]
        E=rad.cumulative(L,v['u'])[-1]; stress=rad.cumulative(d['L_scalar_stress'],d['t'])[-1]
        u=d['t']-r+1.5; expectedL=2*(.03*u/.25*np.exp(-u*u/.5))**2
        error=abs(E/exact-1); peak_error=max(abs(L-expectedL))/max(expectedL)
        assert min(L)>=0 and E>0 and error<.03 and peak_error<.06
        if n==256: assert error<2e-4 and peak_error<5e-4 and abs(stress/exact-1)<2e-4
        report['pulse_grid'].append(dict(n=n,dx=24/n,r=r,expected_energy=exact,energy=E,
            relative_error=error,peak_relative_error=peak_error,stress_energy=stress))
for r in (4,6,8):
    e=[v['relative_error'] for v in report['pulse_grid'] if v['r']==r]
    order=float(np.log2(e[-2]/e[-1])); assert order>3.5
    next(v for v in report['pulse_grid'] if v['n']==256 and v['r']==r)['order']=order

for group,label in (('killing','stationary'),('evolution2','production_gauge_control')):
    for n in (128,256,512):
        path=BASE/f'{group}/stationary-{n}-{n*3//8}/ems_radiation.csv'; authenticate(path)
        raw=np.genfromtxt(path,delimiter=',',names=True)
        for r,d in rad.load(path).items():
            assert d['t'][0]==0 and d['t'][-1]==3
            v=rad.prepare(d,1); L=rad.luminosities(v['news'],v['ds'],v['P'],1)
            energy=rad.cumulative(L,v['u'])[-1]
            pi=raw['Pi_rms'][raw['r']==r][0]; assert pi>1e-6
            if label=='stationary':
                assert max(L[:,1])<1e-10 and max(L[:,2])<1e-10
                assert max(energy[1:])<1e-10
            report[label].append(dict(n=n,dx=32/n,r=r,expected=0 if label=='stationary' else None,
                max_L_phi=float(max(L[:,1])),max_L_EM=float(max(L[:,2])),
                E_phi=float(energy[1]),E_EM=float(energy[2]),Pi_rms_initial=float(pi)))
for r in (4,6,8):
    e=[v['max_L_phi'] for v in report['stationary'] if v['r']==r]
    assert np.log2(e[0]/e[1])>6

for threshold in ('1e-7','1e-10','1e-12'):
    path=BASE/f'threshold/{threshold}/static-n48-dx32'
    authenticate(path/'params.txt'); authenticate(path/'rh_surf_0.dat')
    rows=[line.split() for line in (path/'rh_surf_0.dat').read_text().splitlines()
          if line.strip() and not line.lstrip().startswith('#')]
    last=rows[-1]; assert last[-1]=='found'
    residual=float(last[9]); assert residual<=float(threshold)
    fresh=float(re.findall(r'EMSRH_FINAL_FRESH 0 (\S+)',(path/'run.log').read_text())[-1])
    assert fresh<=float(threshold)*1.001
    report['threshold'].append(dict(threshold=float(threshold),residual=residual,fresh=fresh,
        area=float(last[4]),updates=(path/'run.log').read_text().count('RHFinder::update'),
        seconds=json.loads((path/'status.json').read_text())['seconds']))
assert report['threshold'][-1]['area']<report['threshold'][0]['area']

amr=rad.load(BASE/'amr/amr/ems_radiation.csv')
assert all(len(d['t'])==9 and d['t'][0]==0 and d['t'][-1]==.5 for d in amr.values())
report['amr']={'radii':list(amr),'samples_per_radius':9,'window':[0,.5],'dt':.0625,'passed':True}
report['off_identity']=json.loads((BASE/'off-final/off-identity.json').read_text())
report['baseline']=json.loads((BASE/'off-final/baseline.json').read_text())
report['parameter_checks']=json.loads((BASE/'parameters/parameter-checks.json').read_text())
report['regression_status']=json.loads((BASE/'regressions/grid/status.json').read_text())
assert len(report['regression_status'])==26 and all(v['exit']==0 for v in report['regression_status'])
fixture=(BASE/'regressions/trumpet/run.log').read_text()
for fingerprint in ('e1107d452abc4e4b','e68444ed5f33ba7b','024d0092dba552bc'):
    assert fingerprint.lstrip('0') in fixture
assert 'b3b32ee0d0a729dc' in (BASE/'regressions/ctt/run.log').read_text()
report['fixture_fingerprints']=['e1107d452abc4e4b','e68444ed5f33ba7b','024d0092dba552bc','b3b32ee0d0a729dc']
for src,dst in [('analytic-verified/analytic.json','analytic.json'),
                ('cas-final.log','cas.log'),('pilot-dryrun.log','pilot-dryrun.log'),
                ('parameters/parameter-checks.json','parameter-checks.json'),
                ('regressions/trumpet/run.log','trumpet.log'),('regressions/ctt/run.log','ctt.log'),
                ('postprocess/summary.json','postprocess.json'),
                ('analytic-verified/pulse-33/run.log','tensor-check.log')]:
    shutil.copyfile(BASE/src,OUT/dst)
native=[json.loads(p.read_text()) for p in sorted((BASE/'native-rh-pilot').glob('*/status.json'))]
assert len(native)==4 and all(p['exit']==0 for p in native)
(OUT/'native-rh-pilot.json').write_text(json.dumps(native,indent=2)+'\n')
for path in [Path('/Users/auroradysis/Workspace/EMS/artifacts/reference-alpha20-qfile.trumpet'),
             *ROOT.glob('Examples/EMS/*.ex'),*ROOT.glob('Tests/EMSRadiation/*.ex'),
             *ROOT.glob('Tests/EMSRHFinder/*.ex')]: authenticate(path)
# The single production callback retains every legacy expression verbatim.
original=subprocess.check_output(['git','show',report['baseline']['revision']+
    ':Examples/EMS/EMSBH2DLevel.cpp'],cwd=ROOT,text=True).split('void EMSBH2DLevel::specificPostTimeStep()',1)[1]
current=(ROOT/'Examples/EMS/EMSBH2DLevel.cpp').read_text().split('void EMSBH2DLevel::specificPostTimeStep()',1)[1]
assert current.split('    if (m_radiation.active &&',1)[0]+'}\n'==original
assert not (ROOT/'Examples/EMSRadiation/EMSRadiationLevel.cpp').exists()
report['legacy_callback_preserved']=True
report['passed']=True
(OUT/'gates.json').write_text(json.dumps(report,indent=2)+'\n')
with (OUT/'pulse-grid.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(report['pulse_grid'][-1]))
    writer.writeheader(); writer.writerows(report['pulse_grid'])
print('PASS: pulse, stationary, threshold, AMR, default-off identity, fixture fingerprints, 26 regressions')
