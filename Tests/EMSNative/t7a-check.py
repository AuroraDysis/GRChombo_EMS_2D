#!/usr/bin/env python3
"""One bounded check for T7 A sampling and retained numerical evidence."""
import importlib.util,csv
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('a',HERE/'t7a-analyze.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
def rows(name):return list(csv.DictReader((HERE/name).open()))
a.t6.selfcheck()
h=1/48;xx,yy=np.meshgrid(np.arange(int(256/h)-144,int(256/h)+144),np.arange(144));iv=np.column_stack((xx.ravel(),yy.ravel()));x=(iv[:,0]+.5)*h-256;y=(iv[:,1]+.5)*h
for n in (6,8):
 for degree in range(n):
  v=np.zeros((len(iv),33));v[:,0]=(np.sqrt(2)*x/4)**degree
  got=a.line_sample(iv,v,h,6,n)[2,0];use=(a.R>.75)&(a.R<4)
  assert np.max(abs(got[use]-(a.R[use]/4)**degree))<2e-12
  v[:,0]=(x/8)**degree*(y/16)**2;v[:,2]=(x/8)**degree*y/16
  got=a.tensor_sample(iv,v,h,6,n);use=np.isfinite(got[0,0]);assert use.any()
  assert np.max(abs(got[0,0,use]))<2e-12 and np.max(abs(got[0,2,use]))<2e-12
  use=np.isfinite(got[2,0]);q=a.R[use]/np.sqrt(2);expect=(q/8)**degree*(q/16)**2
  assert np.max(abs(got[2,0,use]-expect))<2e-12
for r in rows('t7a-run-audit.csv'):
 assert r['exit']=='0' and r['sigma']=='1' and r['amr_transfer']=='point'
 assert r['floor_cells_captured']=='0' and r['GaussB_nonzero']=='0'
 assert float(r['stage_time_error_max'])<1e-12 and float(r['ghost_replay_error_max'])<32*np.finfo(float).eps
 assert float(r['native_rhs_replay_error_max'])<1e-11
for r in rows('t7a-ratios.csv'):
 assert abs(float(r['spatial_apparent_order'])-np.log2(float(r['clock'])/float(r['space'])))<1e-12
for r in rows('t7a-projection-summary.csv'):
 assert float(r['stage_projection_Theta_absmax'])==float(r['end_projection_Theta_absmax'])==0
 assert float(r['stage_projection_Ham_absmax'])<1e-15
for r in rows('t7a-restriction-summary.csv'):
 assert float(r['replay_max_error'])<32*np.finfo(float).eps and r['condition'].startswith('COVERED_COARSE_ONLY')
for case in ('clock','space'):
 for lev in (4,5,6):
  p=HERE/f't7a-rk-budget-{case}-L{lev}.csv'
  if not p.exists():p=Path('/private/tmp/ems-t7a/tables')/p.name
  rs=list(csv.DictReader(p.open()));assert rs and max(float(r['RK4_replay_error']) for r in rs)<1e-18
  for r in csv.DictReader(a.artifact(f't7a-global-operation-deltas-{case}-L{lev}.csv').open()):
   if r['operation'] in ('stage_projection','update_projection','end_projection') and r['field']=='Theta':assert float(r['abs_max'])==0
print('PASS: six/eight-node polynomial sampling and parity; complete-run finite audits; stage times; ghost/RHS/RK replay; zero direct Theta projection; ratios')
