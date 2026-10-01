#!/usr/bin/env python3
"""Small tables and scientific figures for the completed T12 t=0 audit."""
import sys
sys.dont_write_bytecode=True
import csv,json,math,resource
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;TMP=Path('/private/tmp/ems-t12');LEGS=('E-low','E-mid','E-high')
def read(name):return list(csv.DictReader((HERE/('t12-'+name+'.csv')).open()))
def save(name,rows):
 with (HERE/('t12-'+name+'.csv')).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def summaries():
 p=read('native-profiles');rows=[]
 for leg in LEGS:
  for ray in ('axis','diagonal'):
   for mode in ('original','maximal'):
    for l in range(7,13):
     a=[x for x in p if x['leg']==leg and x['ray']==ray and x['mode']==mode and x['level']==str(l)];g=np.array([float(x['Gamma_rhs']) for x in a]);k=int(abs(g).argmax())
     row=dict(leg=leg,mode=mode,ray=ray,level=l,points=len(a),peak=g[k],peak_r_M=a[k]['r_M'],RMS=float(np.sqrt(np.mean(g*g))))
     for term in ('shift_laplacian_2D','shift_graddiv_2D','cartoon_shift','lapse_A','chi_A','matter','KO'):row[term+'_peak_abs']=max(abs(float(x[term])) for x in a)
     rows.append(row)
 save('source-summary',rows)
 rows=[]
 for leg in LEGS:
  for mode in ('original','maximal'):
   for ray in ('axis','diagonal'):
    z=np.load(TMP/f'{leg}-{mode}-native.npz');m=z['meta'];q=z['q'];ids=np.where((m[:,0]==12)&(m[:,2]>0)&((m[:,3]==0) if ray=='axis' else (abs(m[:,4]-m[:,5])<m[:,1]/100)))[0]
    i=ids[np.argmin(np.hypot(m[ids,4],m[ids,5]))];cc=np.load(TMP/f'{leg}-{mode}-continuum.npz')['q'];proj=lambda a,c:(a[c] if ray=='axis' else (a[c]+a[c+1])/np.sqrt(2))
    rows.append(dict(leg=leg,mode=mode,ray=ray,x_M=m[i,4],y_M=m[i,5],R_M=np.hypot(m[i,4],m[i,5]),native_alpha_advection=q[i,84+13],continuum_alpha_dot=cc[i,28+13],native_beta_advection=proj(q[i,84:112],14),continuum_beta_dot=proj(cc[i,28:56],14),native_B_decay=proj(q[i,141:143],0),continuum_B_dot=proj(cc[i,28:56],16),native_h11_dot=q[i,29],native_h12_dot=q[i,30],continuum_h11_dot=cc[i,29],continuum_h12_dot=cc[i,30],KO_A12=q[i,56+7],KO_alpha=q[i,56+13],KO_phi=q[i,56+18],KO_Pi=q[i,56+19],KO_Ex=q[i,56+24]))
 save('stage-b-inputs',rows)
 rows=[]
 for leg in LEGS:
  a=read('census-'+leg+'-maximal');rows.append(dict(leg=leg,boxes=len(a),points=sum(int(x['points']) for x in a),**{key:min(float(x[key]) for x in a) for key in ('min_R','min_Y','min_G','min_alpha_K','min_chi','alpha_floor_margin','chi_floor_margin')}))
 save('reader-summary',rows)
 with (TMP/'continuum-memory.csv').open() as f:a=list(csv.DictReader(f))
 save('resources',[dict(stage='continuum',chunk_points=1024,peak_RSS_bytes=max(int(x['peak_RSS_bytes']) for x in a if not x['comparison'].startswith(('bounded','fixed'))),wall_s=15.0277945,completed=True),dict(stage='measure',chunk_points=1024,peak_RSS_bytes=486227968,wall_s=5.7125087,completed=True),dict(stage='native_high',chunk_points=0,peak_RSS_bytes=548700160,wall_s=8.0072943,completed=True)])
def figures():
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 try:
  import scienceplots;plt.style.use(['science','no-latex'])
 except ImportError:pass
 plt.rcParams.update({'font.size':9,'figure.dpi':150,'axes.grid':True,'grid.alpha':.2})
 colors=['#245685','#cc7134','#4e855e'];out=HERE/'figures';out.mkdir(exist_ok=True)
 def finish(fig,name):
  fig.savefig(out/(name+'.png'),dpi=240);fig.savefig(out/(name+'.pdf'));plt.close(fig)
 p=read('native-profiles');fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained')
 for ir,ray in enumerate(('axis','diagonal')):
  for im,mode in enumerate(('original','maximal')):
   ax=axs[ir,im]
   for leg,color in zip(LEGS,colors):
    a=[x for x in p if x['ray']==ray and x['mode']==mode and x['leg']==leg and x['covered']=='False'];a.sort(key=lambda x:float(x['r_M']))
    ax.plot([float(x['r_M']) for x in a],[float(x['Gamma_rhs']) for x in a],color=color,label=leg)
   ax.set_xscale('log');ax.set_yscale('symlog',linthresh=1e-5);ax.set_xlim(4e-5,.5);ax.set_ylim((-1,400) if mode=='original' else (-1,5));ax.set_title(ray+(' (first native row)' if ray=='axis' else '')+', '+mode);ax.set_xlabel('Ray coordinate [M]');ax.set_ylabel('Native longitudinal Gamma RHS');ax.legend()
 finish(fig,'t12-gamma-source')
 p=read('reader-line');r=np.array([float(x['r_M']) for x in p]);a=np.array([float(x['alpha_K']) for x in p]);chi=np.array([float(x['chi']) for x in p]);Y=np.array([float(x['Y']) for x in p]);G=np.array([float(x['G']) for x in p]);ak=float(json.loads((HERE.parents[1]/'scripts/cas/t12-reader-budget.json').read_text())['a_K'])
 fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained')
 axs[0,0].loglog(r,a,label='alpha K');axs[0,0].loglog(r,np.sqrt(chi),label='sqrt chi');axs[0,0].set_ylim(1e-7,2);axs[0,0].legend();axs[0,0].set_ylabel('Initial lapse')
 use=r<=1e-4;axs[0,1].semilogx(r[use],a[use]/r[use]**1.3372155112113842);axs[0,1].axhline(ak,color=colors[1],ls='--',label=f'a K = {ak:.6f}');axs[0,1].set_ylabel('alpha K / R^nu');axs[0,1].legend()
 axs[1,0].loglog(r,Y,label='Y');axs[1,0].loglog(r,G,label='G');axs[1,0].set_ylim(.04,2);axs[1,0].set_ylabel('Reader positivity');axs[1,0].legend()
 axs[1,1].loglog(r,a/1e-12,label='alpha K / floor');axs[1,1].loglog(r,chi/1e-12,label='chi / floor');axs[1,1].axhline(1,color='black',ls='--');axs[1,1].set_ylim(1,2e12);axs[1,1].set_ylabel('Configured floor ratio');axs[1,1].legend()
 for ax in axs.flat:ax.set_xlabel('Physical R [M]')
 finish(fig,'t12-reader-checks')
 p=read('identity-profiles');fig,axs=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
 for mode,color in zip(('original','maximal'),colors):
  a=[x for x in p if x['mode']==mode and float(x['r_M'])<=.5];r=[float(x['r_M']) for x in a]
  axs[0].loglog(r,np.maximum([abs(float(x['source'])) for x in a],1e-25),color=color,label=mode)
  axs[1].loglog(r,np.maximum([abs(float(x['identity3_defect'])) for x in a],1e-25),color=color,label=mode+' defect')
  if mode=='maximal':axs[1].loglog(r,[float(x['Float64_term_budget']) for x in a],color='black',ls='--',label='128 eps operation scale')
 axs[0].set_ylim(1e-20,1e3);axs[1].set_ylim(1e-22,1e-6)
 axs[0].set_ylabel('|Continuum Gamma RHS| (Float64)');axs[1].set_ylabel('Absolute identity defect and budget')
 for ax in axs:ax.set_xlabel('Physical R [M]');ax.legend()
 finish(fig,'t12-identity-budget')
def report():
 parts=['''
## T12 — matched maximal initial lapse, Stage A

**READY.** The original initialized continuum longitudinal Gamma source is reproduced. Assigning the reader's own Killing lapse at t=0 removes that source to the independently measured representation/evaluation budget. The finite-difference RHS remains nonzero; this is the registered Stage A reading, not evidence that the emitted pulse has been removed. No evolution was launched, no equations/gauge/KO/transfers were changed, and no commit was made.

### Change and controls

`ems_use_maximal_initial_lapse = true` selects the already computed `adm_vars.lapse` only at the final evolved-lapse assignment. Omitted/false retains the original `sqrt(chi)` assignment. The internal maximal lapse used for curvature, scalar momentum, Maxwell data and boost construction is untouched. The parameter parser and setter reject boosted flags, nonzero boost rapidity (even with `boosted=false`), binary flags, nonzero separation and a CTT companion; direct boosted/binary helper calls are also rejected. The parser additionally requires EMS `emstrumpet1` data. This card supports an unboosted single hole only.

The frozen `abb0a93` setter headers were compiled with the same flags and harness as the current setter. All 13-level, per-box raw initialized-state byte fingerprints match with the option off; every non-lapse field fingerprint matches between the two options. There are 6,461,700 initialization/three-ghost-cell samples including repeated box overlaps. Gamma calculation and ExperimentalGauge initialization have no lapse dependency; their actual saved fields, including gauge B, are retained in the native replay. Every component of the T11 native input/physical RHS/KO/advection/term output is bit-identical with the option off. Between options, direct array bit comparisons find zero changes to all 27 non-lapse native state components and to KO(A). The exact native driver combination is zero in both modes.

''']
 def table(headers,rows):
  parts.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n')
  for row in rows:parts.append('| '+' | '.join(str(x) for x in row)+' |\n')
  parts.append('\n')
 c=read('controls');table(['Rung','Native cells','Default output bit mismatches','Non-lapse state bit mismatches','Driver cancellation','Maximal Gamma split defect / epsilon'],[(x['leg'],x['cells'],x['default_rhs_bit_mismatches'],x['nonlapse_bits_changed'],x['driver_cancel_maximal'],f"{float(x['maximal_term_defect_eps']):.4f}") for x in c])
 parts.append('''The native comparison covers levels 7–12 and the complete saved seven-by-seven stencil support on both rays; 82,807 cells times 28 evolved variables gives 2,318,596 components. The census fingerprints cover all levels, not only these ray windows. Full-state fingerprints are byte witnesses, not an interval proof. The non-lapse invariance also follows directly from the unchanged ADM construction and sole changed final assignment. CSVs: [controls](t12-controls.csv), [parser controls](t12-parser-controls.csv), [setter rejections](t12-rejections.csv), and the nine per-box `t12-census-*` tables.

The rebuilt EMSTrumpet fixture test passes the reference and echo B/E tables, asserted fingerprints and malformed-file controls. Its reference CCZ4 fingerprint remains `e68444ed5f33ba7b`; echo E remains `bffc958d04850e8f`. The existing reference N=64 unboosted grid check passes (full-mask H RMS 7.07603635e-7, Mx/My RMS 1.68953132e-5, GaussE RMS 1.12292660e-8, GaussB=0, zero floors). The rebuilt EMSNative operator checks pass with default and maximal parameter settings, and the unsupported parser settings fail with the intended MayDay message. Build EMSNative with `_base_srcs=PointTransferTest.cpp`: Chombo's generic one-main rule otherwise attempts to link/compile the independently built offline harnesses too. No build-system or library change was needed.

### Actual Float64 reader census and floor margins

Every valid initialization point and three-cell box halo was checked on all levels 0–12, using the actual separate multiply/subtract coordinates and the reader's Float64 inverse map. Every sampled Y, G and alpha_K is positive and finite. The minimum alpha/chi floor margins are the actual values minus the configured 1e-12 thresholds passed to `PositiveChiAndAlpha`; the class's default 1e-4 thresholds are not used by these E parameters. No clipping or replacement profile was introduced. Axis parity and same-level exchange are used for native ghosts; new coarse-fine lapse ghosts use the existing centred degree-five point operator on the current initialized coarse lapse. Every native maximal lapse after those fills also exceeds its floor. Ancestor L6 values supply valid support for L7; none of the evaluated supports reaches L6's outer face.

''')
 r=read('reader-summary');table(['Rung','Points including box halos','Minimum Y','Minimum G','Minimum alpha_K','Minimum chi','alpha_K minus floor','chi minus floor'],[(x['leg'],x['points'],*[f'{float(x[k]):.9g}' for k in ('min_Y','min_G','min_alpha_K','min_chi','alpha_floor_margin','chi_floor_margin')]) for x in r])
 parts.append('''The independent reader line has 8,009 radii: logarithmic 1e-6–1e3 M, oversampled 1e-4–0.5 M, and fixed comparison radii. Its smallest radius is a diagnostic approach to the puncture, not a new grid placement; R=0 is excluded. The floor margins there remain positive (alpha_K≈4.01195e-7, chi≈5.43021e-11 at 1e-6 M). The reader-to-analytic-jet lapse value discrepancy is at most 2.22045e-15. The positive endpoint constants are Y0=0.135703294942, G0=0.911204279678 and a_K=42.330259894. Fitting alpha_K over 1e-6–1e-5 M gives exponent 1.33717675 against nu=1.3372155112113842; alpha_K/R^nu ranges 42.3257207–42.3300504 and approaches a_K. It is C1, not C2, at the puncture; the option does not regularize all fractional-power data.

![Reader positivity, fractional asymptotics and configured floor ratios](figures/t12-reader-checks.png)

The upper-right panel approaches the positive asymptotic coefficient; the lower panels retain the entire radial sample range. These plotted checks and the complete census establish the finite sampled claims, not positivity at every real radius. Data: [reader line](t12-reader-line.csv), [census summary](t12-reader-summary.csv).

### Native source: puncture cells, all levels and fixed physical radii

The following values are the closest positive diagonal cell, physical R=sqrt(2)h12/2, with the radial longitudinal projection. Gamma RHS excludes KO, which is separately retained in every term profile. R changes with resolution here; these rows are not fixed-radius convergence tests.

''')
 a=read('puncture-terms');rows=[]
 for leg in LEGS:
  o=next(x for x in a if x['leg']==leg and x['ray']=='diagonal' and x['mode']=='original');k=next(x for x in a if x['leg']==leg and x['ray']=='diagonal' and x['mode']=='maximal')
  rows.append([leg,f"{float(k['physical_R_M']):.9g}",f"{float(o['lapse']):.9g}",f"{float(k['lapse']):.9g}",f"{float(o['Gamma_rhs']):.9f}",f"{float(k['Gamma_rhs']):.9f}",f"{float(k['continuum_Gamma']):.3e}",f"{float(k['A12_KO']):.6f}"])
 table(['Rung','Physical R/M','sqrt(chi)','alpha_K','Original Gamma RHS','Maximal Gamma RHS','Maximal continuum RHS','KO(A12), both modes'],rows)
 o=next(x for x in a if x['leg']=='E-mid' and x['ray']=='diagonal' and x['mode']=='original');k=next(x for x in a if x['leg']=='E-mid' and x['ray']=='diagonal' and x['mode']=='maximal')
 table(['E-mid same-cell term','Original','Maximal'],[(label,f'{float(o[key]):.9g}',f'{float(k[key]):.9g}') for label,key in [('2D shift Laplacian','shift_laplacian_2D'),('2D grad-div shift','shift_graddiv_2D'),('Cartoon shift terms','cartoon_shift'),('Lapse-gradient × A','lapse_A'),('chi-gradient × A','chi_A'),('Matter','matter'),('Gamma KO','KO')]])
 parts.append('''The cancellation involves all shift second derivatives, lapse/chi terms and scalar momentum. The remaining native puncture source is a discretization residual. The first native row at y=h/2 is an axis proxy, not a true-axis sample; its coordinate column is x and physical R is stored separately. Its maximal-lapse Gamma component peak is 3.71577814 / 3.24023015 / 2.82535081, giving shrinking-cell peak orders 0.337745 / 0.337913, close to nu−1. The diagonal maximum magnitudes are 0.38984269 / 0.33024429 / 0.29777945; their peak cells differ. Do not interpret either moving-cell sequence as a fixed-domain pulse convergence result. The original first-row peaks are 225.427199 / 249.704926 / 270.341072. Covered coarse cells are labeled and excluded from the plotted active hierarchy; all native level 7–12 term rows, including covered cells, are retained for inspection.

![Native original and maximal-lapse Gamma sources, all rungs](figures/t12-gamma-source.png)

The original physical source remains large on the same physical profile. With the maximal lapse, the residual becomes concentrated in shrinking puncture-end cells and decreases at fixed radii. The entire signed term split is retained in [native profiles](t12-native-profiles.csv), its compact [source summary](t12-source-summary.csv), and the [puncture table](t12-puncture-terms.csv).

The following are true fixed-physical-radius P6 samples on the diagonal, with P8 sensitivity retained in [fixed-radii](t12-fixed-radii.csv). At R=0.001 M, the low-rung P6/P8 spread is 0.033808 for the original source and 0.003806 for the maximal one, so that point is not a qualified fine-difference convergence witness. At 0.01 M the new spreads are 4.67e-13 / 1.72e-10 / 1.34e-10; at 0.1 M the native residual is a small nonmonotone floor.

''')
 a=read('fixed-radii');rows=[]
 for R in (.001,.01,.1):
  for leg in LEGS:
   o=next(x for x in a if x['leg']==leg and x['ray']=='diagonal' and x['mode']=='original' and float(x['requested_r_M'])==R);k=next(x for x in a if x['leg']==leg and x['ray']=='diagonal' and x['mode']=='maximal' and float(x['requested_r_M'])==R)
   rows.append([R,leg,*[f'{float(x):.9g}' for x in (o['lapse'],k['lapse'],o['Gamma_rhs'],k['Gamma_rhs'])]])
 table(['R/M','Rung','Original lapse','alpha_K','Original Gamma RHS','Maximal Gamma RHS'],rows)
 a=read('rhs-orders');table(['Fixed-radius window/M','Mode','Axis order','Diagonal order'],[(lo+'–'+hi,mode,*[f"{float(next(x['order'] for x in a if x['mode']==mode and x['ray']==ray and x['r_min_M']==lo)):.6f}" for ray in ('axis','diagonal')]) for lo,hi in [('0.0002','0.001'),('0.001','0.005'),('0.005','0.02'),('0.02','0.05'),('0.05','0.1'),('0.1','0.49')] for mode in ('original','maximal')])
 parts.append('''Orders are log(D_low_mid/D_mid_high)/log(1.5), on identical physical points, excluding every face's five-low-cell interpolation collar. Added exact fixed radii account for the slight difference from T11's window orders; the underlying original native output is bit-identical. The qualified 0.005–0.02 M maximal audit gives low/mid differences 1.138876e-4 axis and 2.528499e-5 diagonal, mid/high differences 2.258638e-5 and 4.987690e-6; interpolation fractions are 0.000404 and 0.002493. Inner-window orders have substantial interpolation uncertainty (up to 0.783 of the diagonal fine difference); the negative and low outer diagonal orders remain in the table. Those outer differences are 4.4e-8 down to 9.5e-11, consistent with the previously measured coordinate/metric derivative floor. This floor interpretation is inferred; no rounding-method change was tested. None of these orders is silently replaced by four.

### Continuum identity (3), with independently measured residuals

Use the actual reader's b=-CR/r³, k=C exp(delta)/r³, chi=(R/r)² and dR/R=exp(-delta)dr/(r alpha_K). Define M_beta=b'−b/R−3alpha_K k and C_M=2k'+6k/R−3k chi'/chi−16pi Pi phi'. The verified off-constraint identity is

```
S_Gamma[alpha] = 4[k(alpha_K-alpha)]' + 12k(alpha_K-alpha)/R
               + (4/3)(M_beta' + 3 M_beta/R) + 2 alpha C_M.
```

The three requested identities close exactly in a rational jet algebra, together with the map/shift relation and metric-strain residual. The independent arbitrary-function check has 24 seeded points, 60 digits, maximum normalized defect 2.878e-61. Separately, 14 radii per lapse use a forward Chebyshev recurrence and independent high-precision differentiation of the exact Float64 coefficients. Neither M_beta nor C_M is set to zero. The high-precision C_M near 1e-4/1e-3 M is 5.84007e-14 / 1.05674e-13; M_beta is numerically approximately 1e-61 there. The remaining alpha_K continuum source is exactly the independently measured residual contribution to that numerical precision.

''')
 a=read('identity-budget');rows=[]
 for R in (.0001,.001,.01,.1):
  o=min([x for x in a if x['mode']=='original'],key=lambda x:abs(float(x['r_M'])-R));k=min([x for x in a if x['mode']=='maximal'],key=lambda x:abs(float(x['r_M'])-R))
  rows.append([R,*[f'{float(x):.6e}' for x in (o['source_Float64'],k['source_Float64'],k['source_60digit'],k['identity3_Float64_defect'],k['term_roundoff_budget'])]])
 table(['R/M','Original source, Float64','Maximal source, Float64','Maximal residual source, 60 digits','Identity (3) defect, maximal','128 epsilon operation budget'],rows)
 parts.append('''At the 28 independent comparisons, the worst Float64/60-digit source discrepancy is 5.25517e-10 (at R=1e-6 M), and the worst identity defect is 1.51907e-10. Across all 8,009 radii per mode, identity defect/budget is below 0.021. The dense maximal source maximum is 7.28552e-10, dominated by cancellation arithmetic at very small radii, versus original 394.889720 at the line's minimum radius. The original puncture limit remains 419.671110806. At 1e-4 M, the alpha_K source is 1.94e-12 in Float64 and 2.21e-17 with independently evaluated derivatives.

Identity (3) was also evaluated at every native comparison point, with radial projection using its actual x,y direction. For low/mid/high, the maximal-lapse continuum source maxima are 1.79144e-12 / 6.36564e-12 / 4.80133e-12; independently evaluated Float64 residual-source maxima are 1.61244e-12 / 4.72847e-12 / 3.32123e-12. The identity defects are at most 1.74211e-12 / 3.29286e-12 / 1.96240e-12, with defect/budget ratios below 0.028. The measured Float64 C_M maxima are 4.09273e-11 / 3.99751e-11 / 4.45368e-11; unlike the 60-digit represented residuals, these include cancellation error in the reader's arithmetic. [Same-point identity summary](t12-native-identity.csv) keeps both lapse modes and the measured M_beta and its derivative.

The roundoff estimate retains the magnitudes of b'', b'/R and b/R² before their cancellations, with multiplier 128 unchanged. The first collapsed-shift estimate underestimated the R=1e-6 arithmetic error (1.519e-10 defect versus 1.263e-11 collapsed budget); both estimates are retained. This is a measured/estimated evaluation budget, not an interval certificate. The exact CAS status proves the algebraic identity; the runtime status corroborates only its transfer. [CAS evidence and scope](../../scripts/cas/t12-evidence.md), [independent budget](t12-identity-budget.csv), [dense identity profiles](t12-identity-profiles.csv).

![Continuum source cancellation and identity-defect budget](figures/t12-identity-budget.png)

The matched continuum source reaches the residual/evaluation scale throughout the puncture region; its native counterpart need not be bitwise zero. All four stop conditions pass: no positivity/floor failure; no changed non-lapse data; no loss of default identity; and no failure of identity (3) on independently validated inputs. This closes Stage A's source test only.

### Nonzero initial tendencies and Stage B capture

The continuum tendencies are alpha_dot=b alpha', beta_dot=b b', B_dot=−0.1 B (B=−beta in the ideal initial state), and h_dot_ij=2k(alpha_K−alpha)(3n_i n_j−delta_ij), with the measured M_beta strain residual retained. At the closest E-mid diagonal cell, continuum alpha_dot is 2.687573e-4 original / 9.254105e-5 maximal, beta_dot=1.322076e-5 and B_dot=3.659024e-6 in both modes. Continuum h11_dot is 1.948965e-3 / 1.11e-16 and h12_dot 5.846895e-3 / 1.48e-16. Native maximal h11/h12 tendencies remain −7.81320e-5 / 4.10976e-4. Neither the initial lapse nor the driver is stationary under the frozen gauge. Data for every rung/ray and both modes: [Stage B inputs](t12-stage-b-inputs.csv).

KO(A12) remains 19311.828 / 29010.314 / 43552.439, bit-identical between options. At the mid cell, lapse KO changes 0.342139 → 0.0940075, while scalar KO −0.0837899, Pi KO 0.943445 and Ex KO −119608.378 are unchanged. Large local KO is not yet an attribution of emitted characteristic amplitude.

The next experiment is consult 7's **matched E-mid original/alpha_K launch pair**, not a new source ladder or transport enlargement. It must capture the actual first native RK stages' lapse, A, metric, Gamma, shift/driver B, K/Theta, EMS fields and carried constraints, with geometric RHS, advection, KO and projection contributions separated. Use synchronized fixed-time axis/diagonal profiles without phase alignment or background fitting; include a bounded dt/2 control. Stop at 0.002 M; the fixed window is 0.00075–0.0025 M and probe radius 0.0015 M. Require incoming-characteristic margins from every face, complete ghost-support collars and a measured interface-error bound. The unchanged finest face is at 0.02734375 M, but a time-step snapshot cannot replace first-stage history.

The registered promotion screen is a factor-two reduction in **both peak and RMS Gamma disturbance on both rays**, significant against interpolation/temporal errors, with no new floors or failure of the existing absolute constraint screens. A reduction could be delayed emission; later exterior tests remain necessary. Stage B was not launched here. Consult 7 supersedes T11's previous Stage-C amplitude-stability requirement: a pulse converging to zero must not fail for lacking a stable nonzero amplitude.

### Reproduction, resource incident and provenance

Use `t12-analyze.py build`, `census E-low|E-mid|E-high`, `guards`, `native E-low|E-mid|E-high`, `continuum`, `measure`, `identity-points`, `check`; then the two `scripts/cas/t12-*` scripts and `t12-report.py summaries`, `figures`, `report`. Python is `/Users/auroradysis/miniconda3/bin/python`, with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` and `MPLCONFIGDIR=/private/tmp/ems-t12/mpl` for figures. The build reuses T7-E's standalone serial Chombo helper; no library source was modified. Frozen headers from `git show abb0a93:...` are retained under `/private/tmp/ems-t12/baseline`. Regression logs, native/continuum caches, exact box lists and per-chunk memory samples are retained under `/private/tmp/ems-t12`; static inputs and all EMS repository data remain read-only. The manifest pins those retained inputs and artifacts. Three dense regenerable CSVs total 13.3 MB and are registered in `LARGE-OUTPUTS-NOT-COMMITTED.txt`; the remaining tables are small.

**Resource incident:** the first continuum attempt was killed by the controller after reaching 22 GB (20 GB compressed), violating this card's 8 GB limit. NumPy array-left multiplication with the lightweight Jet class caused quadratic object allocations. The corrected operation order and 1024-point disk-backed chunks remove that path. The complete rerun took 15.028 s with measured peak RSS 304,218,112 bytes (0.304 GB); measure peaked at 486,227,968 bytes and native high at 548,700,160 bytes. Per-chunk measurements are retained and hashed. No other completed census/native result was restarted. Every command remained below the detached threshold, no process is pending, and no files belonging to another tranche were deleted. [Resource table](t12-resources.csv), [COMMIT-MANIFEST-T12.txt](COMMIT-MANIFEST-T12.txt).
''')
 text=''.join(parts);path=HERE/'README.md';old=path.read_text();assert '\n## T12 —' not in old;path.write_text(old.rstrip()+'\n\n'+text.lstrip())
if __name__=='__main__':
 if sys.argv[1]=='summaries':summaries()
 elif sys.argv[1]=='figures':figures()
 elif sys.argv[1]=='report':report()
 print(sys.argv[1],'peak RSS bytes',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
