#!/usr/bin/env python3
"""T14 scientific figures, evidence checks and registered decision."""
import sys
sys.dont_write_bytecode=True
import csv, json, math, os, re
from pathlib import Path
import numpy as np
import importlib.util
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t14');OUT=ROOT/'analysis'
spec=importlib.util.spec_from_file_location('t14',HERE/'t14-analyze.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)

def read(name):return list(csv.DictReader((HERE/name).open()))
def f(row,key):return float(row[key])

def baseline():
    for leg in ('max13','max14','max14-half'):a.closure_hashes(leg)
    for suffix in ('stage-attribution','stage-budgets','projections','profile-update-closure','predictor-feedback'):
        rows=[dict(x,run='max12') for x in read('t13-'+suffix+'.csv') if x['run']=='maximal']
        a.save('t14-'+suffix+'.csv',rows,merge=True)
    row=next(x for x in read('t13-stage-closure.csv') if x['directory'].endswith('/maximal'))
    a.save('t14-stage-closure.csv',[dict(run='max12',**row)],merge=True)
    a.puncture('max12')

def figures():
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'mpl'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    try:
        import scienceplots
        plt.style.use(['science','no-latex'])
    except ImportError:pass
    plt.rcParams.update({'font.size':8,'legend.fontsize':7,'axes.titlesize':9})
    target=HERE/'figures';target.mkdir(exist_ok=True)
    def finish(fig,name):
        fig.savefig(target/(name+'.pdf'));fig.savefig(target/(name+'.png'),dpi=600);plt.close(fig)
    hist=read('t14-history.csv');screen=read('t14-screen-history.csv')
    colors={'max12':'C0','max13':'C1','max14':'C2'}
    fig,axes=plt.subplots(2,2,figsize=(7,4.8),layout='constrained')
    for row,(ray,vec) in enumerate(a.a.NV.items()):
        for leg,color in colors.items():
            z=np.load(a.cache(leg,ray));times=z['times'];index=int(np.argmin(abs(times-.001993815104166667)))
            profiles={n:(z[f'p{n}'][index,:,11:13]-z[f'p{n}'][0,:,11:13])@vec for n in (6,8)}
            axes[row,0].plot(z['r'],profiles[6],label=leg,color=color)
            axes[row,0].fill_between(z['r'],np.minimum(profiles[6],profiles[8]),np.maximum(profiles[6],profiles[8]),color=color,alpha=.2)
            h=[x for x in hist if x['run']==leg and x['ray']==ray and f(x,'time_M')<=.001993815104166668]
            axes[row,1].plot([f(x,'time_M') for x in h],[f(x,'Gamma_RMS') for x in h],label=leg,color=color)
        axes[row,0].set(title=ray+', common endpoint',xlabel='R / M',ylabel=r'$\Delta\widetilde\Gamma_n\ [M^{-1}]$')
        axes[row,1].set(title=ray+', amplitude history',xlabel='t / M',ylabel=r'RMS $\Delta\widetilde\Gamma_n\ [M^{-1}]$')
        for ax in axes[row]:ax.legend();ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    finish(fig,'t14-gamma-disturbance')
    fig,axes=plt.subplots(2,2,figsize=(7,4.8),layout='constrained')
    for ax,field in zip(axes.flat,a.FIELDS):
        for ray,color in [('axis','C0'),('diagonal','C1')]:
            for kind,style in [('RMS','-'),('peak','--')]:
                rows=[x for x in screen if x['field']==field and x['ray']==ray and x['norm']==kind and f(x,'time_M')>0]
                ax.plot([f(x,'time_M') for x in rows],[f(x,'rho_D') for x in rows],color=color,ls=style,label=ray+' '+kind)
                bad=[x for x in rows if x['status']=='uncertainty-dominated']
                ax.scatter([f(x,'time_M') for x in bad],[f(x,'rho_D') for x in bad],s=13,facecolors='none',edgecolors=color)
        ax.axhline(.8,color='black',ls=':',lw=1);ax.set(title=field,xlabel='t / M',ylabel=r'$\rho_D$')
        ax.set_ylim(0,.85);ax.legend();ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    finish(fig,'t14-contraction')
    fig,axes=plt.subplots(2,4,figsize=(9,4.8),layout='constrained')
    for row,ray in enumerate(('axis','diagonal')):
        for col,field in enumerate(a.FIELDS):
            rows=[x for x in screen if x['field']==field and x['ray']==ray and x['norm']=='RMS' and f(x,'time_M')>0]
            ax=axes[row,col];t=[f(x,'time_M') for x in rows]
            curves={'D12,13':[f(x,'D12_13') for x in rows],'D13,14':[f(x,'D13_14') for x in rows],
                'Δt/2':[f(x,'temporal_difference') for x in rows],
                '5× spread13,14':[5*f(x,'interpolation13_14') for x in rows]}
            for label,v in curves.items():ax.plot(t,v,label=label,ls=':' if 'spread' in label else '-')
            positive=[v for curve in curves.values() for v in curve if v>0]
            ax.set_yscale('log');ax.set_ylim(max(min(positive)*.5,1e-20),max(positive)*2)
            ax.set(title=ray+', '+field,xlabel='t / M',ylabel='RMS profile norm')
            if row==0 and col==0:ax.legend()
            ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    finish(fig,'t14-temporal-control')
    puncture=read('t14-puncture-budgets.csv');budgets=read('t14-stage-budgets.csv')
    fig,axes=plt.subplots(2,2,figsize=(7,4.8),layout='constrained')
    for col,group in enumerate(('total','KO')):
        for leg,color in colors.items():
            rows=[x for x in puncture if x['run']==leg and x['x_cell']==x['y_cell']=='0']
            axes[0,col].plot([int(x['step']) for x in rows],[f(x,group) for x in rows],marker='o',color=color,label=leg)
        axes[0,col].set(title='First puncture cell: signed '+group,xlabel='own native step',ylabel=r'Weighted $\Delta\widetilde\Gamma_n$')
        axes[0,col].legend()
    for col,ray in enumerate(('axis','diagonal')):
        for leg,color in colors.items():
            rows=[x for x in budgets if x['run']==leg and x['ray']==ray and x['group']=='KO']
            axes[1,col].plot([int(x['step']) for x in rows],[f(x,'KO_over_total_RMS') for x in rows],marker='o',color=color,label=leg)
        axes[1,col].set(title=ray+': direct KO norm ratio on W',xlabel='own native step',ylabel='KO norm / total norm')
        axes[1,col].legend()
    finish(fig,'t14-stage-attribution')
    print('Four T14 scientific figure pairs written.',flush=True)

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
        ['| '+' | '.join(map(str,row))+' |' for row in rows])

def report():
    rows=read('t14-screen.csv');history=read('t14-screen-history.csv');audits=read('t14-evolved-audit.csv')
    endpoint_pass=all(x['endpoint_status']=='contracting' for x in rows)
    unresolved=sum(int(x['uncertainty_dominated_times'])+int(x['numerical_floor_times']) for x in rows)
    noncontracting=sum(int(x['noncontracting_times']) for x in rows)
    marginal=sum(int(x['inconclusive_times']) for x in rows)
    temporal_fail=sum(int(x['temporal_unqualified_times']) for x in rows)
    verdict='insufficient' if noncontracting else 'inconclusive' if unresolved or marginal or temporal_fail else 'pass'
    result=dict(verdict=verdict,endpoint_pass=endpoint_pass,unresolved_field_ray_norm_time_rows=unresolved,
        significant_noncontracting_rows=noncontracting,marginal_rows=marginal,temporal_unqualified_rows=temporal_fail,
        stage_D_drafted=verdict=='pass',production_admitted=False)
    (OUT/'verdict.json').write_text(json.dumps(result,indent=2)+'\n')
    assert verdict!='pass','A passing screen requires a separately completed Stage D design.'
    text=['\n### Completed Stage C results\n',
        '**Endpoint contraction PASS; full registered history INCONCLUSIVE.** All 16 endpoint field/ray/norm entries have significant rho_D<=0.8, and every interpolation-qualified nonzero history entry contracts. However, '+str(unresolved)+' field/ray/norm/time entries fail the registered >5x sampling-spread requirement. These are retained as unresolved rather than reclassified by their small point ratios. There is no significant noncontraction, no 0.8–1 marginal ratio and no temporal failure. Under the registered unresolved-uncertainty decision, this mixed screen does not advance Stage D or admit production. No additional rung or changed KO/gauge method is automatically authorized.\n',
        'The restarted control completed with both new markers equal to zero in 1002.798 s, peak RSS **1,854,324,736 bytes**, and sampled RSS/footprint active-tree maximum **1.997461 GB**. All 449 control steps are retained; step 448 is compared to the 56/112/224-step source ladder at **0.001993815104166667 M**. Its actual final endpoint is 0.001998265584309896 M. The old gate-stopped output remains retained and is excluded from this completed control. The restart output is approximately 3.9 GiB, below its isolated 6 GB gate.\n',
        'All 57 common clocks, including t=0, match within **4.33681e-19 M**. No temporal interpolation, pulse alignment or fitted background was applied. Signed raw Gamma, native metric Gamma including ww/hww, beta_n and lapse use the unchanged anchor-subtracted tensor P6/P8 sampling at the same 258 physical points on W. Peak is the fixed-window maximum; RMS is the unweighted radial integral norm defined in T13. The Gamma disturbance is the change from each rung\'s captured numerical t=0 state, followed by fixed interpolation. Neither this numerical baseline nor any native replay opens the static solution.\n',
        'The following endpoint pairs are **peak / RMS**. A sampling margin is D/(5 times the sum of the two P6/P8 profile-difference norms), so >1 qualifies. The listed margin is the smaller of the two spatial-difference margins; temporal fraction is the larger of the dt/2 differences divided by either spatial norm. These are measured sensitivity spreads, not rigorous continuum-error intervals.\n']
    compact=[]
    for ray in ('axis','diagonal'):
        for field in a.FIELDS:
            pair=[next(x for x in rows if x['ray']==ray and x['field']==field and x['norm']==kind) for kind in ('peak','RMS')]
            compact.append([ray,field,' / '.join(f"{f(x,'endpoint_rho'):.6f}" for x in pair),
                ' / '.join(f"{min(f(x,'endpoint_margin12_13'),f(x,'endpoint_margin13_14')):.2f}" for x in pair),
                ' / '.join(f"{max(f(x,'endpoint_temporal_over_D12_13'),f(x,'endpoint_temporal_over_D13_14')):.3g}" for x in pair),
                ' / '.join(str(56-int(x['uncertainty_dominated_times']))+'/56' for x in pair)])
    text += [table(['Ray','Field','Endpoint rho_D peak / RMS','Minimum endpoint sampling margin','Finest temporal fraction','Qualified history times'],compact)+'\n',
        'Endpoint rho_D spans **0.017546–0.166404**. Every endpoint sampling margin is at least **5.297** times the already multiplied 5x threshold; every endpoint temporal fraction is <=**1.79176e-4**, far below 0.2. Across all common nonzero times, the largest finest temporal fraction is **0.00390641**. Across all qualified history entries, the largest point ratio is **0.343967**. History-supremum ratios are **0.017546–0.202392**. Thus there is substantial resolved contraction at the endpoint and over the qualified histories; this is neither a fourth-order claim nor qualification of the unresolved rows. The max12 dt/2 profile error is carried separately in [temporal CSV](t14-temporal.csv); its maximum fraction of D12,13 is **0.0265074**, also subdominant.\n']
    uncertain=[]
    for row in rows:
        bad=[x for x in history if x['ray']==row['ray'] and x['field']==row['field'] and x['norm']==row['norm'] and x['status']=='uncertainty-dominated' and f(x,'time_M')>0]
        if bad:uncertain.append([row['ray'],row['field'],row['norm'],len(bad),f"{min(f(x,'time_M') for x in bad):.9g}–{max(f(x,'time_M') for x in bad):.9g}"])
    text += [table(['Ray','Field','Norm','Unresolved clocks','Range of unresolved times / M'],uncertain)+'\n',
        'Initial Gamma/metric-Gamma contain small numerical derivative residuals, not bitwise zeros: measured axis Gamma profile peaks are approximately **0.91e-12 / 1.48e-12 / 3.14e-12** on the three rungs. Initial raw Gamma D12,13 / D13,14 peaks are **1.80e-12 / 3.68e-12 axis** and **1.42e-12 / 3.93e-12 diagonal**. The sum of each pair\'s numerical t=0 zero-field norms supplies a retained empirical floor; differences must also exceed five times that floor. Their initial ratios are explicitly undefined rather than interpreted as growth or an order. No nonzero common-time entry is dominated by this initialization floor. A first verification assumed exact-zero raw Gamma and failed; that assumption was corrected, with its resource metadata retained. Disturbances remain differences from the actual numerical initial state.\n',
        'Nonzero sampling-unresolved entries have interpolation margins below one, down to **0.033884** for the finer axis shift peak difference; temporal error is subdominant even there. They are spatial sampling limitations, not activated chi/lapse floors. Axis lapse peak differences remain unqualified at 23 clocks through 0.000818888 M; the two metric-Gamma peak exceptions occur at 0.000320435 and 0.000356038 M. There is no registered exception permitting these rows to be dropped from an unconditional history pass. [Endpoint/history summaries](t14-screen.csv), [every common-time screen](t14-screen-history.csv).\n',
        '![Signed Gamma disturbance by rung and RMS history](figures/t14-gamma-disturbance.png)\n',
        '![Contraction histories; hollow markers fail the sampling threshold](figures/t14-contraction.png)\n',
        '![Finest dt/2 error, both spatial differences and sampling spreads](figures/t14-temporal-control.png)\n']
    amplitude=read('t14-history.csv');amp=[]
    for ray in ('axis','diagonal'):
        for leg in ('max12','max13','max14'):
            candidates=[x for x in amplitude if x['run']==leg and x['ray']==ray and f(x,'time_M')<=.001993815104166668]
            row=max(candidates,key=lambda x:f(x,'time_M'))
            amp.append([ray,leg,f"{f(row,'Gamma_peak'):.9g}",f"{f(row,'Gamma_RMS'):.9g}"])
    text += [table(['Ray','Rung','Endpoint Gamma disturbance peak / M^-1','RMS / M^-1'],amp)+'\n',
        'The residual disturbance **stabilizes on these three rungs**, consistently with a nonzero resolved profile rather than decrease toward zero. max13→max14 peak changes are -0.36% axis/-0.40% diagonal; RMS changes are +1.88%/+1.25%. Both rays approach peak about 9.82e-5 M^-1 and RMS about 8.99e-5 M^-1. Coarse→middle changes are larger; no amplitude-stability tolerance was imposed as an admission condition. The common-time signed profiles and probe history are retained in [profiles](t14-profiles.csv) and [amplitude history](t14-history.csv).\n']
    lobes=[x for x in read('t14-native-lobes.csv') if x['run']!='max14-half']
    text += [table(['Ray','Rung','Half-height span / native cells','Span / M','Clipped by W'],
        [[x['ray'],x['run'],f"{f(x,'width_cells'):.3f}",f"{f(x,'width_M'):.9g}",x['clipped']] for x in lobes])+'\n',
        'Every dominant half-height lobe is **clipped by W**. These spans are lower bounds on full lobe width, not complete FWHMs. The recorded span increases from roughly 11 to 48 axial cells and 7 to 34 diagonal cells while its physical span approaches the width of W. Axis native values use the first cartoon row y=h/2, explicitly a proxy; diagonal values use actual centres, with spacing sqrt(2)h. [Native amplitudes and lobe widths](t14-native-lobes.csv).\n',
        table(['Run','Steps','Actual endpoint / M','Minimum chi','Minimum lapse','Peak RSS / GB'],
            [[x['run'],x['finest_steps'],f"{f(x,'actual_endpoint_M'):.15g}",f"{f(x,'chi_min'):.7g}",f"{f(x,'lapse_min'):.7g}",f"{f(x,'peak_rss_bytes')/1e9:.6f}"] for x in audits])+'\n',
        'Every all-level valid/ghost stage/full-step census has **zero chi activations, zero lapse activations and zero nonfinite values**. Configured floors stay 1e-12; the smallest fine chi/lapse remain 3.44027e-8 / 2.99648e-5. All completed logs have zero reader messages after the first advance. The evolved shift-speed proxy was checked on **every finest recorder cell at every native snapshot**, not only the common-time support subset: max13/max14/control maxima are 1.0010670495 / 1.0010670509 / 1.0010670509. This preserves the limiting axis/diagonal incoming margins **0.001715191465 / 0.002380857872 M** beyond 0.002 M, with the declared 7h support buffer and 3h ghosts. W reads no coarse–fine ghosts directly; characteristic separation does not imply compact FD/KO support. [Evolved audit](t14-evolved-audit.csv).\n']
    constraints=read('t14-constraints.csv');compact=[]
    for ray in ('axis','diagonal'):
        for field in ('C_Gamma','Ham','Mom'):
            matches=[next(x for x in constraints if x['run']==leg and x['ray']==ray and x['field']==field and x['norm']=='RMS') for leg in ('max12','max13','max14')]
            compact.append([ray,field,*[f"{f(x,'endpoint'):.8g}" for x in matches]])
    text += [table(['Ray','Constraint RMS','max12','max13','max14'],compact)+'\n',
        'C_Gamma is signed longitudinal evolved Gamma minus native metric Gamma; momentum is the Euclidean norm of its Cartesian components. The tables retain initial, common-endpoint and history-supremum peak/RMS values, P6/P8 spreads, and GaussE for both rays and the control. The [constraint history](t14-constraint-history.csv) retains every selected clock. GaussB is exactly zero at every replayed snapshot point. Refinement reduces the puncture-near residuals shown here, but these launch-window measurements do not replace the full exterior constraint battery. [Constraints](t14-constraints.csv).\n']
    closure=read('t14-stage-closure.csv')
    text += ['### First-stage attribution and closure\n',
        'The unchanged T13Replay executable checks every recorded finest ROI RHS/input for the first four native steps. The max12 results are reused from its pinned T13 run; all new rungs and the control were replayed. **Every actual total RHS and projected input is bit-identical**. The following split-term budgets use raw stencil magnitudes, rather than a cancelled final RHS.\n',
        table(['Rung','Native stage points','RHS floats','RHS / projection bit mismatches','RK scaled eps','Gamma raw-budget fraction'],
            [[x['run'],x['native_samples'],x['RHS_Float64_values'],x['RHS_bit_mismatches']+' / '+x['projected_state_bit_mismatches'],f"{f(x,'RK_update_max_eps'):.6g}",f"{f(x,'Gamma_raw_stencil_budget_fraction'):.6g}"] for x in closure])+'\n']
    predictor=[x for x in read('t14-predictor-feedback.csv') if x['region']=='first puncture cell' and x['run']!='max14-half']
    compact=[]
    for leg in ('max12','max13','max14'):
        matches=[next(x for x in predictor if x['run']==leg and x['variant']==variant) for variant in ('recorded_predictor','without_all_KO','without_A_KO','without_EMS_KO')]
        compact.append([leg,*[f"{float(re.search(r'signed RHS ([^;]+);',x['condition']).group(1)):.8g}" for x in matches]])
    text += [table(['First predictor Gamma physical RHS','Recorded','Remove all KO increments','Remove A KO only','Remove EMS KO only'],compact)+'\n',
        'These offline dependency tests remove only indicated KO increments from the captured first half-step predictor, reapply native projections and evaluate the frozen physical RHS. They never evolve an altered state or call a static reader. The signed changes identify sensitivity of the renewed **local puncture-cell** geometric source to KO feedback, principally A; they do not assign an additive share of the later outgoing pulse to KO. The continuing fixed-window disturbance stabilizes under refinement, so local regeneration does not establish significant noncontraction or warrant changing the numerical method under this card. [Dependency tests](t14-predictor-feedback.csv).\n']
    projections=read('t14-projections.csv');pupdate=read('t14-profile-update-closure.csv')
    max_defect=max(f(x,'absolute_defect') for x in pupdate)
    text += [f"Fixed P6 weighted step-profile closure is within **{max_defect:.6g} M^-1**. Trace removal changes only A components; floors produce zero updates. Gamma, shift, driver, lapse, chi, Theta and EMS fields have exactly zero direct projection update. Every per-operation/component maximum is retained in [projection table](t14-projections.csv); every stage's geometric/advection/direct-KO contributions on W and at the first puncture cell are in [stage attribution](t14-stage-attribution.csv). The additional [puncture stage table](t14-puncture-stages.csv) records all 28 components plus longitudinal Gamma at four cells (0,0), (1,0), (0,1), (1,1); [weighted puncture budgets](t14-puncture-budgets.csv) retain signed increments. KO/total norm ratios can exceed one because signed contributions cancel; they are not energy fractions.\n",
        '![Signed puncture-cell totals/direct KO and W direct-KO norm ratios; each rung uses its own first four native steps](figures/t14-stage-attribution.png)\n',
        '**Verdict:** the strong endpoint and qualified-history contraction plus stabilizing resolved profile are promising, but the unconditional registered history screen remains **inconclusive** because sampling obscures 53 required entries. No Stage D draft is promoted under this result; no automatic ladder extension, production admission, gauge/KO/finder change or 100 M claim follows. The short window cannot exclude delayed emission. A future controller decision would need to qualify those sampling-limited differences or explicitly predeclare a narrower history claim before another qualifying screen.\n',
        'Final analysis peak RSS is **1.198047 GB**, sampled active-tree peak **1.312571 GB**; every recorded process stays below 3 GB. The stage-attribution helper now streams one ordered stage group instead of retaining all sixteen; numerical groups and accumulation order are unchanged, with sixteen-group and native/RK closure checks. The highest evolution peak is **1.854325 GB**. All six per-rung analysis job markers are zero and no job remains pending. Exact hashes of overwritten closure outputs were reconstructed from the retained lossless q arrays and retained per frame; no native RHS rerun was needed for those hashes. Original failed census, output-gated run and first verification failure remain recorded. [Completed resource table](t14-resources.csv), [restart completion](t14-restart-completion.csv), [verification](t14-verification.csv).\n',
        'Reproduce with the measured `t14-run.py` wrapper, `T14_OUTPUT_ROOT=/private/tmp/ems-t14/analysis T14_DISK_CAP_BYTES=8000000000`: `t14-analyze.py extract max13|max14|max14-half`, `screen`, and `stages max13|max14|max14-half`; then `t14-report.py baseline`, `figures`, `check`, `report`, and `t14-prepare.py manifest`. Per-rung jobs and their atomic markers are under `analysis/jobs`; all native caches and exact removed-payload hashes are retained under `analysis`. They run serially, <=4 threads, under measured 3 GB per-process and active-tree gates. Production sources, binaries, parameters, equations, gauge, KO, transfers, precision, puncture formulation and the author\'s finder remain unchanged. No SSH or commit.\n']
    path=HERE/'README.md';old=path.read_text().split('\n### Completed Stage C results\n')[0]
    path.write_text(old+'\n'.join(text))
    print(json.dumps(result,indent=2),flush=True)

def check():
    rows=[]
    def yes(name,value,condition):
        assert condition,(name,value);rows.append(dict(check=name,value=value,result='PASS'))
    audit=read('t14-evolved-audit.csv');screen=read('t14-screen.csv');hist=read('t14-screen-history.csv')
    yes('screen summary rows',len(screen),len(screen)==16)
    yes('common-time screen rows',len(hist),len(hist)==16*57)
    yes('all endpoints contract',max(f(x,'endpoint_rho') for x in screen),all(x['endpoint_status']=='contracting' for x in screen))
    yes('initial numerical-floor Gamma ratios remain undefined',0,all(math.isnan(f(x,'rho_D')) and x['status']=='numerical-floor' for x in hist if x['field'] in ('Gamma','metric_Gamma') and f(x,'time_M')==0))
    yes('clock matching defect',max(f(x,'maximum_clock_defect_M') for x in screen),all(f(x,'maximum_clock_defect_M')<2e-18 for x in screen))
    yes('sampling-unqualified rows retained',sum(int(x['uncertainty_dominated_times']) for x in screen),sum(int(x['uncertainty_dominated_times']) for x in screen)==53)
    yes('zero floor/nonfinite/late reader counts',0,all(int(x[k])==0 for x in audit for k in ('chi_activations','lapse_activations','nonfinite','reader_messages_after_first_advance')))
    closure=read('t14-stage-closure.csv')
    yes('four-rung stage replay identity',sum(int(x['RHS_bit_mismatches']) for x in closure),len(closure)==4 and all(int(x['RHS_bit_mismatches'])==int(x['projected_state_bit_mismatches'])==0 for x in closure))
    yes('stage raw budgets',max(f(x,'Gamma_raw_stencil_budget_fraction') for x in closure),all(f(x,'Gamma_raw_stencil_budget_fraction')<1 and f(x,'recorded_sum_raw_stencil_budget_fraction')<1 and f(x,'RK_update_max_eps')<8 for x in closure))
    yes('removed native closure output hashes',sum(len(list((OUT/(leg+'-stages')).glob('*-output.sha256'))) for leg in ('max13','max14','max14-half')),
        all(len(list((OUT/(leg+'-stages')).glob('*-output.sha256')))==len(list((OUT/(leg+'-stages')).glob('frame-*.npz'))) for leg in ('max13','max14','max14-half')))
    projected=read('t14-projections.csv');protected=set('chi h11 h12 h22 hww K Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split())
    yes('protected projection updates',max(f(x,'max_abs_update') for x in projected if x['component'] in protected),all(f(x,'max_abs_update')==0 for x in projected if x['component'] in protected))
    yes('four scientific figure pairs',4,len(list((HERE/'figures').glob('t14-*.png')))==len(list((HERE/'figures').glob('t14-*.pdf')))==4)
    resources=[json.loads(p.read_text()) for p in OUT.rglob('*.resources.json')]
    yes('analysis memory / bytes',max(x['peak_rss_bytes'] for x in resources),all(x['peak_rss_bytes']<3e9 and x['tree_peak_bytes']<3e9 and not x['gate_reason'] for x in resources))
    a.save('t14-verification.csv',rows);print(json.dumps(rows,indent=2),flush=True)

if __name__=='__main__':
    {'baseline':baseline,'figures':figures,'report':report,'check':check}[sys.argv[1]]()
