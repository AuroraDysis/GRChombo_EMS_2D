#!/usr/bin/env python3
"""T15 current-field reductions; no evolution, static data, or remote access."""
import csv,importlib.util,json,math,os,resource,sys
from pathlib import Path
sys.dont_write_bytecode=True
import h5py
import numpy as np
HERE=Path(__file__).resolve().parent
DATA=Path('/Users/auroradysis/Workspace/EMS/.data')
SUB=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0021/submissions/exp-0021')
STAGE=Path('/private/tmp/ems-t15/E-4th')
LEGS=('E-low','E-mid','E-high','E-4th')
TIMES=np.arange(7)*1.75
sys.path.insert(0,str(SUB))
import audit
spec=importlib.util.spec_from_file_location('sealed_budget',SUB/'budget.py')
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f))
def save(name,rows):audit.save(HERE/name,rows)
def slope(t,v):
 t=np.asarray(t);w=(t-t.mean())/np.sum((t-t.mean())**2)
 return float(w@np.asarray(v)),w
def directory(leg):return DATA/('exp-0021' if leg=='E-4th' else 'exp-0020')/leg

def shell():
 from scipy.interpolate import CubicSpline
 b.selfcheck();rows=[]
 for leg in LEGS[:3]:
  plots=sorted((directory(leg)/'plt').glob('*.hdf5'))
  for rung,t in enumerate(TIMES):
   path=plots[rung*2]
   with h5py.File(path) as f:
    g0=f['level_12'];boxes=g0['boxes'][:];xy=np.array([[int(q[k]) for k in audit.BOX_KEYS] for q in boxes])
    lo=xy[:,:2].min(0);hi=xy[:,2:].max(0)
    a=np.full((34,hi[1]-lo[1]+1,hi[0]-lo[0]+1),np.nan)
    for (x0,y0,x1,y1),raw in audit.blocks(f,12):a[:,y0-lo[1]:y1-lo[1]+1,x0-lo[0]:x1-lo[0]+1]=raw
    assert abs(float(g0.attrs['time'])-t)<1e-9 and audit.names(f)[32]=='GaussE'
    g=dict(a=a,lo=lo,hi=hi,h=float(g0.attrs['dx']),centre=336.,time=t)
   assert np.isfinite(a).all()
   s=np.fromstring((directory(leg)/'qualified'/f't{rung*2:03d}'/'n96'/'shape-0-2.dat').read_text(),sep=' ')
   assert len(s)==98 and abs(s[0]-t)<1e-9 and abs(s[1]-336)<1e-10
   th=(np.arange(96)+.5)*np.pi/96
   shape=CubicSpline(np.r_[-th[:4][::-1],th,2*np.pi-th[-4:][::-1]],np.r_[s[2:6][::-1],s[2:],s[-4:][::-1]])
   inv,ff,vol=b.geometry(a);D=b.displacement(a);density=vol*a[32];sphere=CubicSpline([0,np.pi],[.02,.02])
   for nt,nr in ((96,64),(192,128)):
    theta,rh,wt,r,w=b.quadrature(g,shape,nt,nr);x,y=r*np.cos(theta[:,None]),r*np.sin(theta[:,None])
    for n in (6,8):
     v=b.NORM*np.sum(w*b.sample(density[None],g,x,y,n,odd=[])[0]);qo=b.flux(D,g,sphere,nt,n);qh=b.flux(D,g,shape,nt,n)
     rows.append(dict(leg=leg,checkpoint=path.name,time_M=t,N_theta=nt,N_r=nr,point_nodes=n,Q_H=qh,Q_0p02=qo,flux_gap=qo-qh,signed_Gauss_volume=v,volume_minus_gap=v-(qo-qh),source='CURRENT_PLOT;QUALIFIED_NUMERICAL_SHAPE;SEALED_BUDGET_KERNEL'))
   print('SHELL',leg,t,flush=True)
 save('t15-shell-reference-raw.csv',rows)


def tables():
 q=read(DATA/'exp-0021/qualified-horizons-four-rung.csv');rates=read(DATA/'exp-0021/charge-rates-four-rung.csv')
 hist=[];out=[];projections=[]
 for leg in LEGS:
  by={n:sorted([r for r in q if r['leg']==leg and int(r['N_theta'])==n],key=lambda r:float(r['time_M'])) for n in (48,96)}
  assert all(len(v)==7 for v in by.values())
  raw=read(directory(leg)/'qualified-horizons.csv')
  assert all(r['status']=='FOUND' and int(r['stage'])==2 and float(r['expansion_squared'])<=1e-12 and float(r['theta_minus'])<0 for r in raw)
  base=by[96][0];rows=[]
  for rung,(r,s) in enumerate(zip(by[96],by[48])):
   row=dict(leg=leg,time_M=float(r['time_M']),A=float(r['A']),Q=float(r['Q']),expansion_squared=float(r['expansion_squared']))
   for f in ('A','Q'):
    z=float(r[f]);z0=float(base[f]);delta=z/z0-1
    row.update({f'delta_{f}':z-z0,f'relative_{f}':delta,f'angular_drift_{f}':abs(delta-float(s[f])/float(by[48][0][f])+1),f'stopping_drift_{f}':(float(r[f'stopping_delta_{f}'])+abs(z/z0)*float(base[f'stopping_delta_{f}']))/z0,f'angular_absolute_{f}':float(r[f'angular_delta_{f}']),f'stopping_absolute_{f}':float(r[f'stopping_delta_{f}'])})
   rows.append(row);hist.append(row)
  late=[r for r in rows if r['time_M']>=5.25-1e-9];tt=np.array([r['time_M'] for r in late]);assert len(late)==4
  record=dict(leg=leg,h0_M=(7/8)/(1.5**LEGS.index(leg)),A0=float(base['A']),Q0=float(base['Q']),late_samples=4)
  for f in ('A','Q'):
   vv=np.array([r[f] for r in late]);rate,w=slope(tt,vv)
   n48=[r for r in by[48] if float(r['time_M'])>=5.25-1e-9]
   rate48,_=slope(tt,[float(r[f]) for r in n48]);angular=abs(rate-rate48)
   stopping=float(sum(abs(w[i])*float(by[96][i+3][f'stopping_delta_{f}']) for i in range(4)))
   fitted=vv.mean()+rate*(tt-tt.mean());rms=float(np.sqrt(np.mean((vv-fitted)**2)));se=float(np.sqrt(np.sum((vv-fitted)**2)/2/np.sum((tt-tt.mean())**2)))
   adjacent=np.diff(vv)/np.diff(tt)
   native=next(r for r in rates if r['leg']==leg and r['source']=='native_every_level0_step')
   native_abs=float(native['dQ_per_M']) if f=='Q' else float(native['relative_A_rate_per_M'])*float(base['A'])
   record.update({f'd{f}_per_M':rate,f'relative_{f}_rate_per_M':rate/float(base[f]),f'angular_rate_{f}':angular,f'stopping_rate_{f}':stopping,f'angular_stopping_rate_{f}':angular+stopping,f'fit_rms_{f}':rms,f'fit_slope_SE_{f}':se,f'interval_rate_min_{f}':float(adjacent.min()),f'interval_rate_max_{f}':float(adjacent.max()),f'native_d{f}_per_M':native_abs,f'native_relative_{f}_rate_per_M':native_abs/float(base[f]),f'native_samples':int(native['samples']),f'max_abs_relative_{f}':max(abs(r[f'relative_{f}']) for r in rows),f'endpoint_relative_{f}':rows[-1][f'relative_{f}'],f'max_angular_drift_{f}':max(r[f'angular_drift_{f}'] for r in rows),f'max_stopping_drift_{f}':max(r[f'stopping_drift_{f}'] for r in rows)})
   sealed=next(r for r in rates if r['leg']==leg and r['source']=='qualified_n96')
   assert abs((rate if f=='Q' else rate/float(base[f]))-float(sealed['dQ_per_M' if f=='Q' else 'relative_A_rate_per_M']))<1e-14
  out.append(record)
  if leg in LEGS[2:]:
   delta=record['endpoint_relative_Q'];rr=record['relative_Q_rate_per_M'];u=(record['angular_stopping_rate_Q']+record['fit_slope_SE_Q'])/record['Q0']
   projections.append(dict(leg=leg,endpoint_M=100,measured_relative_Q_t10p5=delta,late_relative_rate_per_M=rr,angular_stopping_and_fit_sensitivity_per_M=u,anchored_late_projection=delta+89.5*rr,projection_sensitivity=rows[-1]['angular_drift_Q']+rows[-1]['stopping_drift_Q']+89.5*u,endpoint_average_projection=delta*100/10.5,charge_budget=2e-4,conditional_crossing_M=10.5+(2e-4-delta)/rr if rr>0 else '',assumption='LATE_5p25_TO_10p5_SLOPE_PERSISTS_TO_100;NOT_AN_EVOLUTION_OR_PREDICTION'))
 save('t15-horizon-history.csv',hist);save('t15-late-rates.csv',out);save('t15-charge-projections.csv',projections)
 ratios=[]
 for x,y in zip(out,out[1:]):
  ratios.append(dict(coarse=x['leg'],fine=y['leg'],signed_fine_over_coarse=y['dQ_per_M']/x['dQ_per_M'],absolute_coarse_over_fine=abs(x['dQ_per_M']/y['dQ_per_M']),ratio_relative_angular_stopping_sensitivity=x['angular_stopping_rate_Q']/abs(x['dQ_per_M'])+y['angular_stopping_rate_Q']/abs(y['dQ_per_M'])))
 save('t15-rate-ratios.csv',ratios)


def orders():
 fields=read(DATA/'exp-0021/self-differences.csv');constraints=read(DATA/'exp-0021/constraint-orders.csv')
 assert len(fields)==37*11*7 and len(constraints)==12*11*7
 groups=[]
 for mask in audit.MASKS:
  for t in TIMES:
   for kind,rows,columns in (('evolved28',[r for r in fields if r['field'] in audit.EVOLVED],('measured_order_low_triplet','measured_order_high_triplet')),('all37',fields,('measured_order_low_triplet','measured_order_high_triplet')),('constraints12',constraints,('order_low_mid','order_mid_high','order_high_fourth'))):
    selected=[r for r in rows if r['mask']==mask and float(r['time_M'])==t]
    rec=dict(mask=mask,time_M=t,kind=kind,rows=len(selected))
    for col in columns:
     values=[float(r[col]) for r in selected if r[col]!='']
     rec.update({col+'_min':min(values) if values else '',col+'_max':max(values) if values else '',col+'_negative':sum(v<0 for v in values),col+'_low_positive':sum(0<=v<3 for v in values),col+'_undefined':len(selected)-len(values)})
    if kind!='constraints12':
     spreads=[float(r['interpolation_spread_4_vs_6'])/min(float(r[k]) for k in ('difference_low_mid','difference_mid_high','difference_high_fourth')) for r in selected if min(float(r[k]) for k in ('difference_low_mid','difference_mid_high','difference_high_fourth'))>float(r['roundoff_floor'])]
     rec.update(interpolation_ratio_max=max(spreads) if spreads else '',interpolation_above_smallest_difference=sum(v>=1 for v in spreads),cells=selected[0]['cells'],coordinate_volume=selected[0]['coordinate_volume'])
    else:rec.update(interpolation_ratio_max='',interpolation_above_smallest_difference='',cells='',coordinate_volume='')
    groups.append(rec)
 save('t15-field-order-summary.csv',[r for r in groups if r['kind']!='constraints12'])
 save('t15-constraint-order-summary.csv',[r for r in groups if r['kind']=='constraints12'])
 # Preserve every negative and undefined value verbatim, with all raw differences/norms.
 save('t15-field-orders.csv',fields);save('t15-constraint-orders.csv',constraints)


def check():
 b.selfcheck()
 rows=read(HERE/'t15-horizon-history.csv');assert len(rows)==28
 r=read(HERE/'t15-late-rates.csv');assert [x['leg'] for x in r]==list(LEGS)
 assert all(len([x for x in rows if x['leg']==leg])==7 for leg in LEGS)
 for x in read(HERE/'t15-constraint-orders.csv'):
  for k,a,c in (('order_low_mid','low','mid'),('order_mid_high','mid','high'),('order_high_fourth','high','fourth')):
   if x[k]:assert abs(float(x[k])-math.log(float(x[a])/float(x[c]))/math.log(1.5))<1e-12
 budget=read(HERE/'t15-checkpoint-budget.csv')
 assert len(budget)==56 and len({r['checkpoint'] for r in budget})==14
 for r in budget:
  assert all(math.isfinite(float(r[k])) for k in ('time_M','Q_H','Q_0p02','flux_gap','signed_Gauss_volume','volume_minus_gap'))
  assert abs(float(r['flux_gap'])-(float(r['Q_0p02'])-float(r['Q_H'])))<1e-16
  assert abs(float(r['volume_minus_gap'])-(float(r['signed_Gauss_volume'])-float(r['flux_gap'])))<1e-16
 assert {r['checkpoint'] for r in budget}=={r['checkpoint'] for r in read(HERE/'t15-checkpoint-census.csv')}
 for cp in {r['checkpoint'] for r in budget}:
  assert {(int(r['N_theta']),int(r['N_r']),int(r['point_nodes'])) for r in budget if r['checkpoint']==cp}=={(96,64,6),(96,64,8),(192,128,6),(192,128,8)}
 print('T15_CHECK_PASS numerical baselines, four late samples, all constraint orders, fourteen complete checkpoint budgets, frozen budget selfcheck')
if __name__=='__main__':
 {'shell':shell,'tables':tables,'orders':orders,'check':check}[sys.argv[1]]()
 print('SELF_PEAK_RSS_BYTES',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)
