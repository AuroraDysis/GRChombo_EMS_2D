#!/usr/bin/env python3
"""Qualified numerical horizons and published cylindrical E norms, no static targets."""
import csv,json,math,hashlib,importlib.util
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scienceplots
HERE=Path(__file__).resolve().parent;INPUT=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0019');FIG=HERE/'figures'
plt.style.use(['science','ieee','no-latex']);FIG.mkdir(exist_ok=True)
SCALES=('low','mid','high');C=('Ham','Mom','GaussE');MASKS=('horizon','inside_inner_ring','cavity','cavity_core','between_rings','outer_ring','outer_ring_core','far','far_core','outer_boundary_shell')
def read(p):return list(csv.DictReader(p.open()))
def save(name,rows):
 with (HERE/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def order(a,b):return math.log(a/b)/math.log(1.5) if a>0 and b>0 else math.nan
def bracket(t,times,values):
 exact=np.where(abs(times-t)<2e-10)[0]
 if len(exact):i=int(exact[0]);return float(values[i]),float(times[i]),float(times[i]),0.
 assert times[0]<=t<=times[-1],(t,times[[0,-1]])
 i=np.searchsorted(times,t)-1;w=(t-times[i])/(times[i+1]-times[i]);return float((1-w)*values[i]+w*values[i+1]),float(times[i]),float(times[i+1]),float(w)
def horizons():
 old=read(HERE/'t6-finder-final.csv');finals=[];runs=[];stages=[];audits=[]
 for r in old:finals.append({k:(float(v) if k not in ('scale','checkpoint','status') else v) for k,v in r.items()})
 for cp in sorted((INPUT/'E-high/chk').glob('*.hdf5')):
  pair=[]
  for n in (48,96):
   d=Path('/private/tmp/ems-t7e/horizons/E-high')/cp.stem/f'n{n}';state=json.loads((d/'status.json').read_text());assert state['exit']==0
   runs.append({k:v for k,v in state.items() if k!='sha256'})
   for p,h in state['sha256'].items():audits.append(dict(checkpoint=cp.name,N_theta=n,input_path=p,sha256=h))
   rr=read(d/'surfaces.csv');last=rr[-1];assert last['status']=='FOUND' and last['stage']=='2' and float(last['expansion_squared'])<=1e-12
   pair.append(rr);stages.extend(dict(scale='E-high',checkpoint=cp.name,**x) for x in rr)
  a,b=pair[0][-1],pair[1][-1];prev=next(x for x in pair[1] if x['stage']=='1');seed=pair[1][0]
  finals.append(dict(scale='E-high',checkpoint=cp.name,time_M=float(b['time']),status='QUALIFIED_FROZEN_GRID',A48=float(a['A']),A96=float(b['A']),Q48=float(a['Q']),Q96=float(b['Q']),angular_delta_A=abs(float(a['A'])-float(b['A'])),angular_delta_Q=abs(float(a['Q'])-float(b['Q'])),stopping_delta_A=abs(float(prev['A'])-float(b['A'])),stopping_delta_Q=abs(float(prev['Q'])-float(b['Q'])),replay_A=float(seed['A']),replay_Q=float(seed['Q']),replay_theta_plus_rms=math.sqrt(float(seed['expansion_squared'])),theta_plus_rms48=math.sqrt(float(a['expansion_squared'])),theta_plus_rms96=math.sqrt(float(b['expansion_squared'])),theta_minus96=float(b['theta_minus'])))
 save('t7e-finder-runs.csv',runs);save('t7e-finder-stages.csv',stages);save('t7e-finder-input-audit.csv',audits)
 histories={};qualified=[];drifts=[];budgets=[]
 for scale in SCALES:
  native=read(INPUT/('E-'+scale)/'finder-values.csv');times=np.array([float(x['t']) for x in native]);histories[scale]=native
  q=sorted((x for x in finals if x['scale']=='E-'+scale),key=lambda x:x['time_M']);base=q[0]
  for r in q:
   A=np.interp(r['time_M'],times,[float(x['Area']) for x in native]) if r['time_M']>0 else math.nan
   Q=np.interp(r['time_M'],times,[float(x['Q_charge']) for x in native]) if r['time_M']>0 else math.nan
   err=np.interp(r['time_M'],times,[float(x['err']) for x in native]) if r['time_M']>0 else math.nan
   native_exact=[x for x in native if float(x['t'])==float(format(r['time_M'],'.8e'))]
   if native_exact:
    v=native_exact[-1];A=float(v['Area']);Q=float(v['Q_charge']);err=float(v['err']);nlo=nhi=float(v['t'])
   elif r['time_M']>0:
    _,nlo,nhi,_=bracket(r['time_M'],times,np.array([float(x['Area']) for x in native]))
   else:nlo=nhi=math.nan
   qualified.append(dict(**r,inrun_time_lo=nlo,inrun_time_hi=nhi,inrun_A=A,inrun_Q=Q,inrun_expansion_rms=math.sqrt(err),inrun_row='NO_T0_ROW' if r['time_M']==0 else ('TIME_INTERPOLATED' if scale=='high' and '000097' in r['checkpoint'] else 'EXACT_PRINTED_TIME')))
   drifts.append(dict(scale=scale,time_M=r['time_M'],A96=r['A96'],Q96=r['Q96'],relative_A=r['A96']/base['A96']-1,relative_Q=r['Q96']/base['Q96']-1,relative_A48=r['A48']/base['A48']-1,relative_Q48=r['Q48']/base['Q48']-1))
  d=[x for x in drifts if x['scale']==scale];last=d[-1]
  budgets.append(dict(scale=scale,terminal_time_M=last['time_M'],terminal_relative_A=last['relative_A'],terminal_relative_Q=last['relative_Q'],max_checkpoint_abs_relative_A=max(abs(x['relative_A']) for x in d),max_checkpoint_abs_relative_Q=max(abs(x['relative_Q']) for x in d),area_budget=1e-3,charge_budget=2e-4,checkpoint_area_pass=max(abs(x['relative_A']) for x in d)<=1e-3,checkpoint_charge_pass=max(abs(x['relative_Q']) for x in d)<=2e-4,interval_certification='ONLY_CHECKPOINTS_QUALIFIED',max_inrun_abs_relative_A=max(abs(float(x['Area'])/base['A96']-1) for x in native),max_inrun_abs_relative_Q=max(abs(float(x['Q_charge'])/base['Q96']-1) for x in native)))
 save('t7e-qualified-horizons.csv',qualified);save('t7e-horizon-drifts.csv',drifts);save('t7e-horizon-budgets.csv',budgets)
 rich=[]
 # Qualified checkpoint interpolation at common times; linear and native-trend sensitivity are both recorded.
 for t in (0.,4.5,9.,10.):
  for field in ('A','Q'):
   vals=[];trend=[]
   for scale in SCALES:
    q=sorted([x for x in finals if x['scale']=='E-'+scale],key=lambda x:x['time_M']);ts=np.array([x['time_M'] for x in q]);vs=np.array([x[field+'96'] for x in q]);base=vs[0]
    value,lo,hi,w=bracket(t,ts,vs);vals.append(value)
    native=histories[scale];nt=np.array([float(x['t']) for x in native]);nv=np.array([float(x['Area' if field=='A' else 'Q_charge']) for x in native])
    anchor=int(np.argmin(abs(ts-t)));shift=0 if t==0 else np.interp(t,nt,nv)-np.interp(ts[anchor],nt,nv)
    trend.append(vs[anchor]+shift)
   da,db=vals[0]-vals[1],vals[1]-vals[2];p=order(abs(da),abs(db));ext=vals[2]+(vals[2]-vals[1])/(1.5**p-1) if da*db>0 and p>0 else math.nan
   rich.append(dict(time_M=t,observable=field,low=vals[0],mid=vals[1],high=vals[2],p=p,extrapolated=ext,signed_difference_ratio=da/db if db else math.nan,time_interpolation='LINEAR_QUALIFIED_CHECKPOINTS',native_trend_low=trend[0],native_trend_mid=trend[1],native_trend_high=trend[2],native_trend_p=order(abs(trend[0]-trend[1]),abs(trend[1]-trend[2])),time_matching_sensitivity_max=max(abs(np.array(vals)-trend)),interpretation='CONDITIONAL_NOT_ORDER_QUALIFICATION'))
 for row in rich:
  field=row['observable'];baselines=[next(x[field+'96'] for x in finals if x['scale']=='E-'+scale and x['time_M']==0) for scale in SCALES]
  delta=[row[scale]/base-1 for scale,base in zip(SCALES,baselines)]
  d1,d2=delta[0]-delta[1],delta[1]-delta[2];p=order(abs(d1),abs(d2));ext=delta[2]+(delta[2]-delta[1])/(1.5**p-1) if d1*d2>0 and p>0 else math.nan
  row.update(relative_low=delta[0],relative_mid=delta[1],relative_high=delta[2],relative_p=p,relative_extrapolated=ext,monotone=(row['low']-row['mid'])*(row['mid']-row['high'])>0)
 save('t7e-horizon-richardson.csv',rich)
 fig,axes=plt.subplots(2,1,figsize=(7,4.8),layout='constrained',sharex=True)
 for scale,color in zip(SCALES,('C0','C1','C2')):
  q=[x for x in finals if x['scale']=='E-'+scale];base=next(x for x in q if x['time_M']==0);native=histories[scale]
  for ax,key,nk,budget in zip(axes,('A','Q'),('Area','Q_charge'),(1e-3,2e-4)):
   ts=[float(x['t']) for x in native];ax.plot(ts,[float(x[nk])/base[key+'96']-1 for x in native],color=color,alpha=.35,lw=.8)
   d=[x for x in drifts if x['scale']==scale];ax.plot([x['time_M'] for x in d],[x['relative_'+key] for x in d],'o-',color=color,label=scale+' qualified')
   ax.axhline(budget,color='.5',ls=':',lw=.7);ax.axhline(-budget,color='.5',ls=':',lw=.7);ax.set_ylabel(r'$\Delta '+key+'/'+key+'(0)$');ax.grid(alpha=.2)
 axes[0].legend(ncol=3);axes[1].set_xlabel(r'$t/M$');axes[0].set_title('48/96-point qualification; pale curves: unqualified in-run finder')
 for suffix in ('png','pdf'):fig.savefig(FIG/f't7e-horizons.{suffix}',dpi=300)
 plt.close(fig);return budgets

def constraints():
 raw=[];groups={}
 for scale in SCALES:
  for x in read(INPUT/('E-'+scale)/'constraint-rms.csv'):
   if x['cell_class']!='all':continue
   r=dict(scale=scale,time_M=float(x['time_M']),mask=x['mask'],constraint=x['constraint'],cells=int(x['cells']),cylindrical_rms=float(x['cylindrical_rms']),maximum=float(x['maximum']),norm='SUM_Y_C2_OVER_SUM_Y_UNCOVERED_NO_H2');raw.append(r)
 for scale in SCALES:
  for mask in MASKS:
   for field in C:
    rr=sorted([x for x in raw if (x['scale'],x['mask'],x['constraint'])==(scale,mask,field)],key=lambda x:x['time_M']);groups[scale,mask,field]=(np.array([x['time_M'] for x in rr]),np.array([x['cylindrical_rms'] for x in rr]),rr)
 common=[float(x) for x in groups['low','far','Ham'][0] if x<=10]+[10.]
 orders=[]
 for t in common:
  for mask in MASKS:
   for field in C:
    vals=[];row=dict(time_M=t,mask=mask,constraint=field)
    for scale in SCALES:
     ts,vs,rr=groups[scale,mask,field];v,lo,hi,w=bracket(t,ts,vs);vals.append(v);row.update({scale+'_rms':v,scale+'_time_lo':lo,scale+'_time_hi':hi,scale+'_time_weight':w,scale+'_cells_lo':next(x['cells'] for x in rr if x['time_M']==lo),scale+'_cells_hi':next(x['cells'] for x in rr if x['time_M']==hi)})
    p,q=order(vals[0],vals[1]),order(vals[1],vals[2]);hi_ts,hi_vs,_=groups['high',mask,field]
    near=float(hi_vs[np.argmin(abs(hi_ts-t))]);p_nearest=order(vals[1],near)
    hi_lo,hi_hi=row['high_time_lo'],row['high_time_hi'];w=row['high_time_weight']
    a=float(hi_vs[np.argmin(abs(hi_ts-hi_lo))]);b=float(hi_vs[np.argmin(abs(hi_ts-hi_hi))]);squared=math.sqrt((1-w)*a*a+w*b*b);p_squared=order(vals[1],squared)
    row.update(order_mid_high_nearest=p_nearest,order_mid_high_interpolate_squared_rms=p_squared,interpolation_order_sensitivity=max(abs(q-p_nearest),abs(q-p_squared)))
    row.update(order_low_mid=p,order_mid_high=q,flag='NEGATIVE' if min(p,q)<0 else ('BELOW_3' if min(p,q)<3 else 'GE_3'),time_match='EXACT' if all(row[s+'_time_lo']==row[s+'_time_hi'] for s in SCALES) else 'LINEAR_RMS_INTERPOLATION',norm='LEGACY_CYLINDRICAL_NO_H2')
    orders.append(row)
 save('t7e-constraint-rms.csv',raw);save('t7e-constraint-orders.csv',orders)
 fig,axes=plt.subplots(len(MASKS),3,figsize=(9,17),layout='constrained',sharex=True)
 for i,mask in enumerate(MASKS):
  for j,field in enumerate(C):
   ax=axes[i,j];values=[]
   for scale in SCALES:
    ts,vs,_=groups[scale,mask,field];ax.semilogy(ts,vs,label=scale,lw=.8);values.extend(vs)
   positive=[v for v in values if v>0];ax.set_ylim(min(positive)/2,max(positive)*2);ax.grid(alpha=.2)
   if i==0:ax.set_title(field);ax.legend(ncol=3,fontsize=5)
   if j==0:ax.set_ylabel(mask.replace('_','\n'),fontsize=6)
   if i==len(MASKS)-1:ax.set_xlabel(r'$t/M$')
 fig.suptitle('E: published cylindrical RMS over uncovered cells (without level-volume factor)')
 for suffix in ('png','pdf'):fig.savefig(FIG/f't7e-constraint-rms.{suffix}',dpi=240)
 plt.close(fig)
 # One simple invariant checks exact matching and the interpolation machinery.
 exact=[r for r in orders if r['time_match']=='EXACT'];assert {round(r['time_M'],8) for r in exact}>={0.,4.375,8.75}
 assert bracket(1.,np.array([0.,2.]),np.array([2.,4.]))[0]==3.
 return orders
def merge_native():
 folder=Path('/private/tmp/ems-t7e/tables')
 for prefix,name in [('local-','checkpoint-localization'),('profile-','far-profiles'),('coordinate-control-','coordinate-controls'),('gamma-control-','gamma-controls'),('performance-','analysis-performance')]:
  rows=[r for path in sorted(folder.glob('t7e-'+prefix+'*.csv')) for r in read(path)]
  assert rows,(prefix,folder)
  columns=list(rows[0]);columns+=sorted(set().union(*(r.keys() for r in rows))-set(columns))
  save('t7e-'+name+'.csv',[{k:r.get(k,'') for k in columns} for r in rows])

def far_figure():
 rows=read(HERE/'t7e-far-profiles.csv');fig,axes=plt.subplots(2,3,figsize=(9,5),layout='constrained',sharex=True,sharey=True)
 for j,(m,h) in enumerate((('000032','000051'),('000064','000097'),('000069','000103'))):
  for scale,step,color in (('mid',m,'C1'),('high',h,'C2')):
   selected=[x for x in rows if x['scale']==scale and step in x['checkpoint']]
   for lev in (3,4):
    data=sorted([x for x in selected if int(x['level'])==lev],key=lambda x:float(x['r_mid_M']))
    for i,c in enumerate(('Ham','Mom')):
     axes[i,j].semilogy([float(x['r_mid_M']) for x in data],[float(x[c+'_rms']) for x in data],color=color,ls='-' if lev==4 else '--',label=scale+f' L{lev}')
   axes[0,j].set_title(f"mid {float(next(x['time_M'] for x in selected)):.3f}" if scale=='mid' else axes[0,j].get_title()+f" / high {float(selected[0]['time_M']):.3f}")
  for ax in axes[:,j]:
   ax.axvline(5.4444444444,color='.5',lw=.7,ls=':');ax.axvline(5.5416666667,color='.5',lw=.7,ls=':');ax.set_ylim(1e-12,1e-8);ax.set_xlim(4,8);ax.grid(alpha=.2)
  axes[1,j].set_xlabel(r'$r/M$')
 axes[0,0].set_ylabel('Hamiltonian RMS');axes[1,0].set_ylabel('Momentum RMS');axes[0,0].legend(fontsize=6)
 fig.suptitle('Current checkpoint fields: all far cells; solid L4, dashed L3; dots mark axial face radii')
 for suffix in ('png','pdf'):fig.savefig(FIG/f't7e-far-profiles.{suffix}',dpi=300)
 plt.close(fig)

if __name__=='__main__':
 b=horizons();o=constraints();merge_native();far_figure();print(json.dumps(b,indent=2))
 for r in o:
  if r['mask']=='far' and r['time_match']=='EXACT':print(r['time_M'],r['constraint'],r['low_rms'],r['mid_rms'],r['high_rms'],r['order_low_mid'],r['order_mid_high'],r['flag'])
