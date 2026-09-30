#!/usr/bin/env python3
"""Small end-to-end numerical checks of the T7-E evidence tables."""
import csv,math,json
from pathlib import Path
H=Path(__file__).resolve().parent
R=lambda n:list(csv.DictReader((H/('t7e-'+n+'.csv')).open()))
q=R('qualified-horizons');assert len(q)==12
assert all(x['status']=='QUALIFIED_FROZEN_GRID' and max(float(x['theta_plus_rms48']),float(x['theta_plus_rms96']))<=1e-6 for x in q)
assert len(R('finder-runs'))==8 and all(int(x['exit'])==0 and float(x['wall_seconds'])<120 for x in R('finder-runs'))
assert all(x['pass_control']=='True' and int(x['differences'])==0 for x in R('controls'))
orders=R('constraint-orders');assert len(orders)==720
assert len(set((x['time_M'],x['mask'],x['constraint']) for x in orders))==720
for x in orders:
 for a,b,field in [('low','mid','order_low_mid'),('mid','high','order_mid_high')]:
  assert math.isclose(float(x[field]),math.log(float(x[a+'_rms'])/float(x[b+'_rms']))/math.log(1.5),abs_tol=1e-12)
 if float(x['time_M']) in (0,4.375,8.75):assert x['time_match']=='EXACT' and float(x['interpolation_order_sensitivity'])<1e-12
for x in R('coordinate-errors'):
 raw=next(r for r in R('constraint-rms') if r['scale']==x['scale'] and r['time_M']=='0.0' and r['mask']==x['mask'] and r['constraint']=='Ham')
 assert int(raw['cells'])==int(x['cells'])
 if x['scale']=='low':assert float(x['coordinate_error_max_M'])==0
for x in R('coordinate-controls'):
 if x['scale']=='high' and x['mask'] in ('inside_inner_ring','cavity','between_rings','cavity_core'):assert float(x['native_Ham'])>10*float(x['chi_fused_Ham'])
local=R('checkpoint-localization')
for x in local:
 if x['mask']=='far':assert x['sampled_cells']==x['estimated_cells']
 if float(x['time_M'])>0:assert x['continuum_Ham_rms']==''  # no static evaluator output after t=0
assert max(int(x['parent_peak_rss_bytes']) for x in R('analysis-performance'))<6*1024**3
assert all(int(x['threads'])==1 for x in R('analysis-performance'))
for x in R('finder-input-audit'):
 if x['input_path'].endswith('RHSurf.hpp'):assert x['sha256']=='3c9e3e946f95375a27477466f4e51b1cafbea2f14392acb2d19128582c98cfd6'
 if x['input_path'].endswith('RHUnion.hpp'):assert x['sha256']=='1174747f4ee8f2e18a61cd2fb0618d1a4f3c6974c9929bcae85fd8f64519f885'
print('PASS: qualification, T6 numerical control, exact time matching, orders, counts, coordinate control, t=0 guard, and load evidence.')
