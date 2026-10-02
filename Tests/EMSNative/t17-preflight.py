#!/usr/bin/env python3
"""Small native tracking/restart check; no horizon or merger qualification."""
import sys
sys.dont_write_bytecode=True
import importlib.util, json, shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('work',HERE/'t17-work.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
text=w.p.replace(w.base(1,56,64),dict(max_steps=2,stop_time=1,
    checkpoint_interval=1,plot_interval=1,verbosity=0,
    ems_binary_refinement='true',ems_track_punctures='true',
    ems_puncture_tracking_level=1,regrid_interval='1'))
d=w.ROOT/'preflight/tracker-initial'
assert w.run_case(d.name,d,w.ROOT/'production-t17.ex',text)==0
# Match the regrid cadence too: regrid_interval=1 changes point-transfer calls,
# even on a uniformly covered hierarchy. The default identity controls use 0.
baseline=w.ROOT/'preflight/matched-regrid-default'
baseline_text=w.p.replace(text,dict(ems_binary_refinement='false',ems_track_punctures='false'))
assert w.run_case(baseline.name,baseline,w.ROOT/'production-t17.ex',baseline_text)==0
count,bits=w.p.compare(baseline/'chk/EMS_000002.2d.hdf5',
    d/'chk/EMS_000002.2d.hdf5');assert count>0 and bits==0
records=[x.split() for x in (d/'punctures.dat').read_text().splitlines()
    if x.strip() and not x.lstrip().startswith('#')]
assert [float(x[0]) for x in records]==[0.,.21875,.4375]
r=w.ROOT/'preflight/tracker-restart';r.mkdir(parents=True,exist_ok=True)
shutil.copyfile(d/'punctures.dat',r/'punctures.dat')
restart=w.p.replace(text,dict(restart_file=d/'chk/EMS_000001.2d.hdf5'))
assert w.run_case(r.name,r,w.ROOT/'production-t17.ex',restart)==0
count2,bits2=w.p.compare(d/'chk/EMS_000002.2d.hdf5',r/'chk/EMS_000002.2d.hdf5')
assert count2>0 and bits2==0
assert 'Read EMSTRUMPET' not in (r/'run.log').read_text()
w.save('t17-tracker-check.csv',[
    dict(check='dual-centre tracking versus original uniform hierarchy',Float64_values=count,bit_mismatches=bits,result='PASS'),
    dict(check='native checkpoint restart versus uninterrupted evolution',Float64_values=count2,bit_mismatches=bits2,result='PASS')])
print('PASS: dual-centre tracking fields and native restart identity; small uniform smoke only.')
