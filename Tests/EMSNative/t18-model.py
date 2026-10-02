#!/usr/bin/env python3
"""Exact census checks; two-anchor critical-path model, explicitly predictive."""
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import mpmath as mp
import sympy as sp

HERE=Path(__file__).resolve().parent
T17=HERE.parents[1].parent/'wt-native-t4/Tests/EMSNative'
EVIDENCE=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0023/submission/phase1-evidence/smoke/census')

def csvrows(path):
    return list(csv.DictReader(path.open()))

def save(name,rows):
    with (HERE/name).open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def strips(boxes):
    # Exact coordinate set witness without allocating the gaps between holes.
    rows=collections.defaultdict(list)
    for x0,y0,x1,y1 in boxes:
        for y in range(y0,y1+1):rows[y].append((x0,x1))
    result=[]
    for y,intervals in sorted(rows.items()):
        merged=[]
        for lo,hi in sorted(intervals):
            if merged and lo<=merged[-1][1]+1:
                assert lo>merged[-1][1], 'boxes overlap'
                merged[-1]=(merged[-1][0],hi)
            else:merged.append((lo,hi))
        result.append((y,tuple(merged)))
    return result

def calibration_loads(name,gamma):
    loads=collections.defaultdict(lambda:collections.defaultdict(lambda:[0,0]))
    for path in (EVIDENCE/name).glob('floor-census.csv.rank*'):
        rank=int(path.suffix.split('rank')[-1])
        for r in csvrows(path):
            a=loads[int(r['level'])][rank]
            a[0]+=int(r['valid_cells']);a[1]+=int(r['ghost_cells'])
    assert set(loads)==set(range(15))
    result=0
    for l in range(15):
        v,g=max(loads[l].values(),key=lambda x:sum(x))
        # Affine differences are nonnegative throughout [0,1] if they are
        # nonnegative at both endpoints. Retain that exact dominance witness.
        assert all(v>=w and v+g>=w+h for w,h in loads[l].values())
        result+=(v+gamma*g)*2**l
    return result

def model():
    boxes=csvrows(HERE/'t18-census-boxes.csv')
    census=csvrows(HERE/'t18-census-loads.csv')
    baseline=csvrows(T17/'t17-census.csv')
    check=[]
    for size in (16,24,32,48,128):
        for l in range(14):
            pred=[tuple(int(r[k]) for k in ('x0','y0','x1','y1')) for r in boxes
                  if int(r['max_box_size'])==size and int(r['level'])==l]
            actual=[tuple(int(r[k]) for k in ('x0','y0','x1','y1')) for r in baseline
                    if r['rung']=='fine' and int(r['level'])==l]
            assert strips(pred)==strips(actual),(size,l,'coverage')
            if size==128:assert sorted(pred)==sorted(actual),(l,'original boxes')
            check.append(dict(max_box_size=size,level=l,coordinate_coverage_identical=True,
                              baseline_boxes_identical=size==128,valid_cells=sum((x1-x0+1)*(y1-y0+1) for x0,y0,x1,y1 in pred)))
    save('t18-census-check.csv',check)
    # Claim contract: T_i = a C_i(gamma) + tau S, gamma in [0,1].
    # The witness is the two exact rational residuals and nonzero determinant.
    # This proves algebraic calibration only, never predictive accuracy.
    gamma=sp.symbols('gamma',nonnegative=True)
    S=sp.Integer(2**15-1)
    low=calibration_loads('E-low-dense',gamma)
    high=calibration_loads('E-high-dense',gamma)
    a=sp.cancel((sp.Integer(1030)-418)/(low-high))
    tau=sp.cancel((418-a*high)/S)
    assert sp.cancel(a*low+tau*S-1030)==0
    assert sp.cancel(a*high+tau*S-418)==0
    determinant=sp.factor(S*(low-high))
    assert sp.Poly(determinant,gamma).degree()==1
    assert all(determinant.subs(gamma,g)>0 and a.subs(gamma,g)>0 and tau.subs(gamma,g)>0 for g in (0,1))
    mp.mp.dps=80;worst=mp.mpf(0)
    for g in (mp.mpf(0),mp.mpf('.25'),mp.mpf(1)):
        rational=sp.Rational(str(g))
        lo=mp.mpf(str(low.subs(gamma,rational)))
        hi=mp.mpf(str(high.subs(gamma,rational)))
        solved=mp.lu_solve(mp.matrix([[lo,int(S)],[hi,int(S)]]),mp.matrix([1030,418]))
        exact=[mp.mpf(str(v.subs(gamma,rational).evalf(80))) for v in (a,tau)]
        worst=max(worst,*[abs(v-w) for v,w in zip(solved,exact)])
    assert worst<mp.mpf('1e-70')
    record=dict(claim='The declared two-anchor affine model reproduces 1030 and 418 s/step',
        algebra='Q(gamma), 0<=gamma<=1; effective halo weight, identical 15-level substep count',
        witness='two cancelled rational residuals are zero; determinant is positive at both endpoints of its affine polynomial',
        determinant=str(determinant),a_seconds_per_weighted_cell=str(a),
        tau_seconds_per_level_substep=str(tau),status='PROVED (algebra only)',
        independent_check='mpmath lu_solve, 80 digits, gamma=0,1/4,1',
        worst_absolute_residual=str(worst),gamma_central=.25,
        caveat='Two rates do not identify halo cost separately or isolate synchronization. '
            'tau is an effective non-cell remainder; optimistic ideal OpenMP 4/t scaling, '
            'zero unmeasured multi-node penalty. Calibration layouts use the saved actual rank census.')
    (HERE/'t18-model-cas.json').write_text(json.dumps(record,indent=2)+'\n')
    coeff={g:(float(a.subs(gamma,g)),float(tau.subs(gamma,g))) for g in (0,sp.Rational(1,4),1)}
    def predict(size,ranks,threads,finest,g=sp.Rational(1,4)):
        rs=[r for r in census if int(r['max_box_size'])==size and int(r['ranks'])==ranks and int(r['level'])<=finest]
        key={0:'max_rank_valid',sp.Rational(1,4):'max_rank_weight_quarter',1:'max_rank_halo'}[g]
        critical=sum(float(r[key])*2**int(r['level']) for r in rs)
        aa,tt=coeff[g]
        return aa*critical*4/threads+tt*(2**(finest+1)-1),critical
    decision_path=HERE/'t18-dt-recommendation.json'
    decision=json.loads(decision_path.read_text()) if decision_path.exists() else {}
    recommended_dt=decision.get('dt_multiplier',.25)
    calibration_rung=decision.get('calibration_rung','fine')
    results=[];matrix=[]
    for finest,rung in ((12,'coarse'),(13,'fine')):
        trial_dt=recommended_dt if rung==calibration_rung else .25
        current=predict(128,32,4,finest)[0]
        for size in (16,24,32,48,128):
            for ranks in (32,64,128,256,512):
                t,c=predict(size,ranks,4,finest)
                sensitivity=[predict(size,ranks,4,finest,g)[0] for g in coeff]
                results.append(dict(rung=rung,max_box_size=size,ranks=ranks,threads=4,
                    weighted_critical_cell_substeps=c,sync_remainder_s=coeff[sp.Rational(1,4)][1]*(2**(finest+1)-1),
                    predicted_s_per_step=t,halo_weight_low_s=min(sensitivity),halo_weight_high_s=max(sensitivity),
                    speedup_vs_same_model_draft=current/t,
                    speedup_vs_T17_1474_659=t and (1474.6591101911486/t if rung=='fine' else 'not supplied')))
            for ranks_per_node,threads in ((32,4),(64,2),(128,1)):
                for nodes in (1,2,4):
                    total=ranks_per_node*nodes
                    t,c=predict(size,total,threads,finest)
                    sensitivity=[predict(size,total,threads,finest,g)[0] for g in coeff]
                    matrix.append(dict(rung=rung,max_box_size=size,ranks_per_node=ranks_per_node,
                        threads=threads,nodes=nodes,total_ranks=total,predicted_s_per_step=t,
                        halo_weight_low_s=min(sensitivity),halo_weight_high_s=max(sensitivity),
                        dt_multiplier=trial_dt,calibration_target=rung==calibration_rung,
                        sample_physical_duration_M=5*1.75*trial_dt,
                        sample_coarse_steps=5,warmup_steps=1,retain_steps='2,3,4,5; step5 includes interval-4 coarse regrid',
                        expected_timed_wall_minutes=5*t/60,
                        additional_inter_node_sync_penalty='unmeasured; each +1 ms/substep adds 8.19/16.38 s per coarse/fine step'))
    save('t18-model.csv',results);save('t18-cluster-matrix.csv',matrix)
    summary=[]
    for r in census:
        if int(r['ranks'])==32:
            v=int(r['valid_cells']);h=int(r['halo_cells'])
            summary.append(dict(max_box_size=r['max_box_size'],level=r['level'],boxes=r['boxes'],
                valid_cells=v,halo_cells=h,ghost_over_valid=(h-v)/v,
                ghost_exchange_doubles_upper=(h-v)*28*4*2**int(r['level'])))
    save('t18-grid-summary.csv',summary)
    print('a,tau at gamma=.25',coeff[sp.Rational(1,4)])
    for size in (16,24,32,48,128):
        print('fine size',size,[round(predict(size,r,4,13)[0],1) for r in (32,64,128,256,512)])
    for rung in ('coarse','fine'):
        best=min((r for r in matrix if r['rung']==rung),key=lambda r:r['predicted_s_per_step'])
        print('best',best)
    return results,matrix,summary

def figure(results,summary):
    os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ems-t18/mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    try:
        import scienceplots
        plt.style.use(['science','no-latex'])
    except ImportError:
        plt.rcParams.update({'font.family':'serif','font.size':9})
    fig,axes=plt.subplots(1,2,figsize=(7.2,2.7),layout='constrained')
    for size in (16,24,32,48,128):
        rs=[r for r in results if r['rung']=='fine' and r['max_box_size']==size]
        axes[0].plot([r['ranks'] for r in rs],[r['predicted_s_per_step'] for r in rs],'-o',label=str(size),ms=3)
    axes[0].set_xscale('log',base=2);axes[0].set_xticks([32,64,128,256,512],['32','64','128','256','512'])
    axes[0].set_ylim(50,250);axes[0].set(xlabel='MPI ranks (4 threads/rank)',ylabel='Predicted fine-grid seconds/step')
    axes[0].legend(title='Max box',ncol=2,fontsize=7)
    rs=[r for r in summary if int(r['level'])==13]
    axes[1].bar([str(r['max_box_size']) for r in rs],[r['ghost_over_valid']*100 for r in rs])
    axes[1].set(xlabel='Maximum box side',ylabel='Ghost cells / valid cells (%)')
    fig.savefig(HERE/'t18-performance.png',dpi=300)
    plt.close(fig)

if __name__=='__main__':
    r,m,s=model()
    figure(r,s)
