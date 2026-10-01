#!/usr/bin/env python3
"""Independent analytic derivatives of Float64 Chebyshev coefficients, 60 digits."""
import csv,json,hashlib
from pathlib import Path
import numpy as np
import mpmath as mp
mp.mp.dps=60
root=Path(__file__).resolve().parents[2];here=root/'Tests/EMSNative';tmp=Path('/private/tmp/ems-t12')
profile=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet');lines=profile.read_text().splitlines()
meta={a[2:].split('=',1)[0]:a.split('=',1)[1] for a in lines if a.startswith('# ')}
cs=[]
for label in 'YDS':
 i=lines.index(label+' 97');cs.append([mp.mpf(float(a)) for a in lines[i+1:i+98]])
def cheb(a,s):
 x=2*s-1;t0=mp.mpf(1);t1=x;v=a[0]+a[1]*x
 for c in a[2:]:t0,t1=t1,2*x*t1-t0;v+=c*t1
 return v
ell,nu,C,phinf=[mp.mpf(float(meta[a])) for a in ('ell','nu','C','phi_inf')]
def Y(s):return cheb(cs[0],s)
def delta(s):return (1-s)**2*cheb(cs[1],s)
def radius(s):return ell*s**(1/nu)/(1-s)
def X(s):return ell*s**(1/nu)/Y(s)
def alphaK(s):return mp.exp(-delta(s))*nu*s*(Y(s)+(1-s)*mp.diff(Y,s))/(Y(s)*(1+(nu-1)*s))
def beta(s):return -C*ell*s**(1/nu)*(1-s)**2/Y(s)**3
def curvature(s):return C*mp.exp(delta(s))*(1-s)**3/Y(s)**3
def phi(s):return phinf+(1-s)*cheb(cs[2],s)
def Pi(s):return -C*mp.exp(delta(s))*(1-s)**4*mp.diff(phi,s)/(Y(s)**2*(Y(s)+(1-s)*mp.diff(Y,s)))
def independent(R):
 guess=(R/ell)**nu if R<ell else 1-ell/R
 s=mp.findroot(lambda z:mp.log(radius(z)/R),(guess/2,min(mp.mpf('.99999'),guess)),solver='secant') if R<ell else mp.findroot(lambda z:mp.log(radius(z)/R),(guess,(guess+1)/2),solver='secant')
 sr=1/mp.diff(radius,s);srr=-mp.diff(radius,s,2)*sr**3
 def d(f):return mp.diff(f,s)*sr
 def dd(f):return mp.diff(f,s,2)*sr**2+mp.diff(f,s)*srr
 b=beta(s);bp=d(beta);bpp=dd(beta);k=curvature(s);kp=d(curvature);ak=alphaK(s);akp=d(alphaK);xx=X(s);xp=d(X);c=2*xp/xx;p=Pi(s);phip=d(phi)
 mb=bp-b/R-3*ak*k;mbp=bpp-bp/R+b/R**2-3*(akp*k+ak*kp);cm=2*kp+6*k/R-3*k*c-16*mp.pi*p*phip
 values={}
 for mode,al,ap in (('original',xx,xp),('maximal',ak,akp)):
  shift=mp.mpf(4)/3*(bpp+2*bp/R-2*b/R**2);lapse=-4*k*ap;chi=-6*al*k*c;matter=-32*mp.pi*al*p*phip
  source=shift+lapse+chi+matter;mismatch=4*(kp*(ak-al)+k*(akp-ap))+12*k*(ak-al)/R
  mr=mp.mpf(4)/3*(mbp+3*mb/R);cr=2*al*cm
  # Budget the operations before their radial cancellations, not the small collapsed shift term.
  operation_scale=mp.mpf(4)/3*(abs(bpp)+2*abs(bp)/R+2*abs(b)/R**2)+abs(lapse)+abs(chi)+abs(matter)
  values[mode]=dict(source=source,mismatch=mismatch,M_beta=mb,M_beta_prime=mbp,C_M=cm,M_beta_source=mr,C_M_source=cr,alpha_K=ak,alpha=al,alpha_dot=b*ap,beta_dot=b*bp,B_dot=mp.mpf('.1')*b,metric_coefficient=2*k*(ak-al),identity3_defect=source-mismatch-mr-cr,term_scale=abs(shift)+abs(lapse)+abs(chi)+abs(matter),operation_scale=operation_scale)
 return values
zo=np.load(tmp/'radial-original-continuum.npz');zk=np.load(tmp/'radial-maximal-continuum.npz');r=np.hypot(*zo['xy'].T)
rows=[];allrows=[]
for mode,z in (('original',zo),('maximal',zk)):
 j=z['radial_jets'];x,ak,b,k,ch,p,ph=[j[:,i] for i in range(7)];a=x if mode=='original' else ak
 mb=b[:,1]-b[:,0]/r-3*ak[:,0]*k[:,0];mbp=b[:,2]-b[:,1]/r+b[:,0]/r**2-3*(ak[:,1]*k[:,0]+ak[:,0]*k[:,1])
 cm=2*k[:,1]+6*k[:,0]/r-3*k[:,0]*ch[:,1]/ch[:,0]-16*np.pi*p[:,0]*ph[:,1]
 mismatch=4*(k[:,1]*(ak[:,0]-a[:,0])+k[:,0]*(ak[:,1]-a[:,1]))+12*k[:,0]*(ak[:,0]-a[:,0])/r
 residual=4/3*(mbp+3*mb/r)+2*a[:,0]*cm;source=z['q'][:,39:41].sum(1)/np.sqrt(2)
 scale=4/3*(abs(b[:,2])+2*abs(b[:,1])/r+2*abs(b[:,0])/r**2)+4*abs(k[:,0]*a[:,1])+6*abs(a[:,0]*k[:,0]*ch[:,1]/ch[:,0])+32*np.pi*abs(a[:,0]*p[:,0]*ph[:,1])
 for i,R in enumerate(r):
  defect=source[i]-mismatch[i]-residual[i];budget=128*np.finfo(float).eps*max(1,scale[i]);assert abs(defect)<=budget
  allrows.append(dict(mode=mode,r_M=R,source=source[i],mismatch=mismatch[i],M_beta=mb[i],M_beta_prime=mbp[i],C_M=cm[i],M_beta_source=4/3*(mbp[i]+3*mb[i]/R),C_M_source=2*a[i,0]*cm[i],residual_budget=residual[i],identity3_defect=defect,Float64_term_budget=budget))
 for target in (1e-6,1e-5,.0001,.0002,.001,.005,.01,.02,.05,.1,.2,.5,1.,10.):
  i=int(abs(r-target).argmin());R=mp.mpf(float(r[i]));v=independent(R)[mode]
  row=dict(mode=mode,r_M=float(R),source_Float64=source[i],source_60digit=float(v['source']),mismatch_Float64=mismatch[i],residual_Float64=residual[i],identity3_Float64_defect=source[i]-mismatch[i]-residual[i],evaluation_discrepancy=abs(source[i]-float(v['source'])),collapsed_term_budget=128*np.finfo(float).eps*max(1,float(v['term_scale'])),term_roundoff_budget=128*np.finfo(float).eps*max(1,float(v['operation_scale'])))
  row.update({name+'_60digit':float(value) for name,value in v.items() if name!='source'})
  assert abs(v['identity3_defect'])<mp.mpf('1e-40')*max(1,v['term_scale'])
  assert abs(row['identity3_Float64_defect'])<row['term_roundoff_budget'],row
  if mode=='maximal':assert abs(row['source_Float64'])<=abs(float(v['M_beta_source']+v['C_M_source']))+row['evaluation_discrepancy']+row['term_roundoff_budget']
  rows.append(row)
for filename,data in (('t12-identity-budget.csv',rows),('t12-identity-profiles.csv',allrows)):
 with (here/filename).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
Y0=Y(mp.mpf(0));G0=Y0+mp.diff(Y,mp.mpf(0));aK0=mp.exp(-delta(mp.mpf(0)))*nu*G0/(Y0*ell**nu)
result=dict(status='CORROBORATED',scope='transfer layer — production physics not certified',definition='Forward Chebyshev recurrence of exactly lifted Float64 coefficients; independent implicit-map differentiation at 60 digits; no constraint imposed',points=len(rows),r_min=min(x['r_M'] for x in rows),r_max=max(x['r_M'] for x in rows),precision_digits=60,profile_sha256=hashlib.sha256(profile.read_bytes()).hexdigest(),nu=str(nu),Y0=str(Y0),G0=str(G0),a_K=str(aK0),maximal_continuum_source_max=max(abs(x['source_Float64']) for x in rows if x['mode']=='maximal'),worst_evaluation_discrepancy=max(x['evaluation_discrepancy'] for x in rows),worst_identity3_Float64_defect=max(abs(x['identity3_Float64_defect']) for x in rows),reader_budget='independently measured residual contribution + Float64/60-digit evaluation discrepancy + 128 epsilon sum of pre-cancellation radial operation magnitudes; multiplier unchanged, collapsed shift-term estimate corrected to include bprime/R and b/R^2')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
