#!/usr/bin/env python3
"""T14d evidence, bounded checks and an unexecuted Stage D admission draft."""
import sys
sys.dont_write_bytecode=True
import ast,csv,importlib.util,json,math,os
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('t14d',HERE/'t14d-analyze.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
ROOT=d.ROOT
def f(x,k):return float(x[k])
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
        ['| '+' | '.join(map(str,row))+' |' for row in rows])
def pairs(rows,ray,field,key):
    return ' / '.join(f"{f(next(x for x in rows if x['ray']==ray and x['field']==field and x['norm']==kind),key):.6f}" for kind in ('peak','RMS'))
def resources():
    rows=[]
    for path in sorted(ROOT.rglob('*.resources.json')):
        x=json.loads(path.read_text());sampled=max((v[0] for v in x.get('per_pid_peaks',{}).values()),default=0)
        rows.append(dict(process=x['process'],peak_rss_bytes=max(x['peak_rss_bytes'],sampled),
            sampled_tree_bytes=x['tree_peak_bytes'],wall_s=x['wall_seconds'],returncode=x['returncode'],
            gate_reason=x['gate_reason'],record=str(path)))
    return rows

def design():
    v=json.loads((HERE/'t14d-verdict.json').read_text());rh=.0060286129392158804
    masks=dict(puncture_excluded=(.005,20.),horizon=(rh,2*rh),inside_inner_ring=(2*rh,.026758927819264387),
        cavity=(.026758927819264387,.11726522703740201),between_rings=(.026758927819264387,.7385298961059745),
        outer_ring=(.5908239168847796,.8862358753271694),far=(4.,8.),cavity_core=(.055,.075),
        outer_ring_core=(.85,1.),far_core=(6.5,8.),exterior_wake=(1.,20.),receiving_side_3p5=(3.5,4.375))
    d.save('t14d-stageD-masks.csv',[dict(mask=k,R_min_M=lo,R_max_M=hi,geometry='fixed isotropic radial shell') for k,(lo,hi) in masks.items()]+
        [dict(mask='outer_boundary_shell',R_min_M=0.,R_max_M=1.,geometry='distance from x=+-336 or y=336; y>1 M')])
    census=[x for x in d.read('t14-census.csv') if x['run']=='max14']
    old=sum(int(x['valid_cells'])*2**int(x['level']) for x in census if int(x['level'])<=12)
    new=sum(int(x['valid_cells'])*2**int(x['level']) for x in census)
    assert old==604495872 and new>old
    factor=new/old;rates={x['leg']:f(x,'steady_median_s') for x in d.read('t8-run-audit.csv')}
    cost=[]
    for leg,h0,steps in (('E-low',7/8,48),('E-mid',7/12,72),('E-high',7/18,108)):
        cost.append(dict(leg=leg,h0_M=h0,h14_M=h0/2**14,dt0_M=h0/4,steps_2p625=steps//4,
            steps_4p375=steps*5//12,steps_10p5=steps,max_level=14,work_factor=factor,
            measured_L12_seconds_per_coarse_step=rates[leg],model_L14_seconds_per_coarse_step=factor*rates[leg],
            model_node_hours=steps*factor*rates[leg]/3600,range_node_hours_low=steps*factor*100/3600,
            range_node_hours_high=steps*factor*260/3600,rank_threads='32 MPI x 4 OpenMP'))
    nominal=sum(x['model_node_hours'] for x in cost);low=sum(x['range_node_hours_low'] for x in cost);high=sum(x['range_node_hours_high'] for x in cost)
    control=2*cost[-1]['model_node_hours'];clow=2*cost[-1]['range_node_hours_low'];chigh=2*cost[-1]['range_node_hours_high']
    d.save('t14d-stageD-cost.csv',cost)
    masks_table=table(['Mask','Fixed radial interval / M'],[[k,f'{lo:.15g}–{hi:.15g}'] for k,(lo,hi) in masks.items()])
    chain_table=table(['Rung','h0 / M','h14 / M','Steps at 2.625 / 4.375 / 10.5 M','Model node hours'],
        [[x['leg'],f"{x['h0_M']:.9g}",f"{x['h14_M']:.9g}",f"{x['steps_2p625']} / {x['steps_4p375']} / {x['steps_10p5']}",f"{x['model_node_hours']:.2f}"] for x in cost])
    text=f'''# Stage D draft — for the operator's admission decision

This is a proposed fresh alpha_K E global three-grid chain, not an authorization or a launched run. T14d's complete registered history is {v['verdict'].upper()}. Its unresolved launch entries remain part of that verdict. This draft is required independently of admission; it proposes no additional interpolation variant and no change to the gauge, KO, equations, transfers, Float64 precision, puncture treatment, boundary prescription or author's RHFinder.

Use exp-0020's domain x relative to the hole in [-336,336] M, y in [0,336] M, centre (336,0), point transfers, sigma=1 and dt_multiplier=0.25. Use fresh E-low/mid/high initializations with ems_use_maximal_initial_lapse=true. Keep every exp-0020 level-0–12 physical face R_l=112/2^l M and add the same two puncture levels to **all three** global rungs: level 13 at 0.013671875 M and level 14 at 0.0068359375 M. This is a uniform 3/2 global-spacing chain, not the factor-two source-only T14 ladder. Extra levels put all three finest spacings within or below T14's tested source-refinement range; they do not establish exterior convergence. Using max_level=12 instead would cost less but would place E-low's source spacing outside T14's three-rung range. Do not splice exp-0020/0021 original-lapse outputs into this new chain.

{chain_table}

Base N1/N2 are 768/384, 1152/576 and 1728/864. No regridding is allowed. Before admission to execution, run the real initialization-only census on every proposed grid, prove inherited unions/faces unchanged, check block alignment and all valid/ghost floor margins, and qualify the Tests-only dense-tag initializer workaround on all three grids. Reuse T14's workaround only after these actual censuses; it addresses initialization set storage without changing evolution. Retain absent-static-file complete-hierarchy restart and recorder/default-path identity controls, and measured <=3 GB RSS per process with <=4 OpenMP threads. Initialization uses the static EMS file only at t=0; no positive-time initialization or static ghost substitution is permitted.

## Fixed domain, masks and time claims

Predeclare the smooth puncture-excluded interior domain **0.005<=R<=20 M**. R_excl=0.005 M was fixed in t14d-registration.json before this re-screen. T14's longest measured half-height span is approximately 0.001712 M and every lobe is clipped by W=[0.00075,0.0025] M. Thus these spans are lower bounds, not complete widths. The chosen exclusion is twice W's outer radius and roughly three times the observed span, allowing a conservative fixed source neighbourhood while remaining inside the fixed horizon shell's 0.0060286 M inner radius. It is not a claim that the whole disturbance fits inside that radius. Once an outgoing or reflected disturbance enters R>=0.005 M, it participates in every declared test; the mask never moves with the pulse.

Retain every exp-0020 mask below, including full masks containing faces/corners and the outer boundary shell. The additional puncture-excluded mask and receiving-side mask do not replace failed masks. Use composite coordinate-volume RMS with weights 2pi*y*h_l^2 over uncovered cells, plus peak norms. Retain independent axis/equator/diagonal rays and fixed patch-edge, convex-corner, axis-junction, same-level seam and smooth-interior classes. Diagnostic support must be declared from the coarsest rung, and it must be identical across rungs. Report absent common classes explicitly.

{masks_table}

The outer_boundary_shell remains distance 0–1 M from x=+-336 or y=336 with y>1 M. The old Sommerfeld shell failure remains an admission issue. A restricted interior result additionally requires an all-relevant-mode arrival bound and the previously required boundary-position control. Without that evidence this draft admits neither a domain-wide claim nor a 100 M extension.

The predeclared history start is the **first positive T14d common clock at which all four fields, both rays and both norms exceed 5x interpolation**, computed as step {v['first_all_significant_step']}, t_start={v['first_all_significant_time_M']:.17g} M. This numerical start is frozen for any fresh Stage D qualifying run. Later T14d exceptions at steps {v['later_sampling_exceptions']} are retained: a first qualifying clock does not establish a fully qualified suffix and cannot change the T14d verdict. In Stage D report all t=0 and early histories, but any field/constraint convergence claim beginning at this start is explicitly narrower than convergence over the whole evolution. Do not move the start if a later clock fails. Preserve the 0.875 M synchronized exterior cadence and exact 2.625, 4.375 and 10.5 M diagnostics. Add bounded early current-state captures and diagnostic-only samples at the fixed start through the first coarse interval; if sampling these subcycled times requires time interpolation, predeclare it and measure its error with the control rather than silently comparing asynchronous data.

## Expected order and complete evidence

Design order remains **four** on the fixed smooth domain. No independent regularity argument licenses a lower expected order here. For the 3/2 chain, report log(||D_LM||/||D_MH||)/log(1.5), the signed fourth-order residual D_LM-(3/2)^4 D_MH, and both direct-to-zero constraint orders. Retain all 28 evolved fields, not only Gamma or lapse. Use the fixed I8 central diagnostic reconstruction with I8/I10 sensitivity, >5x pairwise significance, and a fresh E-high dt/2 temporal control through 10.5 M requiring temporal error <=0.2 of both spatial differences. At faces, use declared valid one-sided support; verify coverage and parity before qualifying results. A sampling/floor-limited row is unresolved, never a fourth-order pass. These samplers modify diagnostics only; evolution transfers and native derivative formulas stay frozen.

The constraint battery includes Hamiltonian, Mom1/Mom2 and their joint norm, C_Gamma components and norm, spatial Z1/Z2/Z from the code's native current metric/Gamma, Theta, det(h)-1, tr(A), electric and magnetic Gauss constraints, and the carried cleaning variables Lambda and Xi. Preserve the implementation's prescribed dimensional scalings, and report raw absolute norms, maxima, cell counts and quadrature weights alongside scaled norms. There are no invented first-order-reduction constraints. Exact zero constraints have undefined logarithmic order. Require both pair orders consistent with design four after measured interpolation, temporal and Float64 uncertainty, and decreasing residuals throughout the declared histories. No absolute constraint admission tolerances are supplied by this task or the exp-0020 contract: the operator must seal those physical tolerances for each scaled constraint before execution; measured orders alone cannot admit the chain. Do not choose tolerances after seeing its outputs.

Horizon retention always covers numerical **t=0 through 10.5 M**, independently of t_start. Re-find each rung's numerical t=0 surface using only saved current fields and a numerical seed; use it for its own A_H(0), Q_H(0) baseline. Keep the author's finder unchanged. Require negative inward expansion, fresh final RMS outgoing expansion <=1e-6, angular-resolution and stopping sensitivities below the drift budgets and inter-rung differences, and qualified coverage of between-cadence excursions. Retain native candidate histories every coarse step and numerical checkpoints at excursions; a FAR/CLOSE label or finite candidate area is not a horizon certificate. The budgets remain max|Delta A/A(0)|<=1e-3 and max|Delta Q/Q(0)|<=2e-4, with decreasing drifts under global refinement. Independent current-field sphere charges at 20,50,100 M and their Gauss/flux/cleaner budgets are reported separately. A later 100 M decision must preserve the same numerical t=0 baselines and full horizon interval; this pilot supplies no 100 M claim.

## Cluster work estimate

Use one node per leg, **32 MPI x 4 OpenMP** (128 allocated cores), as exp-0020. Its measured level-12 median seconds/coarse-step were 256/162/104, within the controller's approximately 100–260 s range. The valid-cell subcycle work model W=sum(N_l*2^l), including evolved covered cells, gives {old:,} -> {new:,} cell-steps per coarse step for E-mid, a factor **{factor:.6f}** with the extra levels. The same factor applies to low/high since every refined level has 32768/73728/165888 valid cells and each base has nine times that count. This is a pricing model, not a measured max14 cluster rate; communication, load balance, setup, finder and I/O overhead require fresh bounded cluster smokes before a sealed submit budget.

The nominal three-rung 10.5 M chain costs approximately **{nominal:.2f} node-hours / {nominal*128:.0f} core-hours**, or **{low:.2f}–{high:.2f} node-hours** using the 100–260 s baseline range scaled by the work factor. A full finest dt/2 control adds approximately **{control:.2f} node-hours**, with range **{clow:.2f}–{chigh:.2f}**. Total nominal plus control is **{nominal+control:.2f} node-hours / {(nominal+control)*128:.0f} core-hours**, with planning range **{low+clow:.2f}–{high+chigh:.2f} node-hours** before uncalibrated overhead. Reuse one physical chain per rung through the earlier diagnostic times rather than paying for three separate endpoints.

Propose four nominal segments of 2.625 M each (12/18/27 coarse steps), with complete-hierarchy checkpoints at every segment boundary and all inherited 0.875 M diagnostics. Split the finest dt/2 control into eight 1.3125 M segments (27 control coarse steps each); retain its exact 4.375 M diagnostic inside the corresponding segment. The range model prices the longest proposed segment below eight node-hours before overhead, but seal allocation limits only after new smokes. Resume only from the complete numerical checkpoint with the positive-time static guard enabled. Keep physical completion distinct from diagnostic qualification; a consumed segment budget, failed restart or missing qualification marker cannot imply admission. No SSH, submission, evolution, finder call or commit is performed by T14d.
'''
    (HERE/'t14d-stageD-design.md').write_text(text)
    print('Stage D draft and fixed masks/cost tables written; no run launched.',flush=True)

def figures():
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'t14d-mpl'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scienceplots
    plt.style.use(['science','no-latex']);plt.rcParams.update({'font.size':8,'legend.fontsize':7})
    target=HERE/'figures'
    rows=d.read('t14d-screen-history.csv')
    fig,axes=plt.subplots(2,4,figsize=(9,4.7),layout='constrained')
    for i,ray in enumerate(('axis','diagonal')):
        for j,field in enumerate(d.FIELDS):
            ax=axes[i,j]
            for kind,style in (('peak','--'),('RMS','-')):
                z=[x for x in rows if x['ray']==ray and x['field']==field and x['norm']==kind and int(x['step'])>0]
                for pair,color in (('12_13','C0'),('13_14','C1')):
                    ax.plot([f(x,'time_M') for x in z],[f(x,'margin'+pair) for x in z],style,color=color,label=pair.replace('_',',')+' '+kind)
            ax.axhline(1,color='k',ls=':');ax.set_yscale('log');ax.set_ylim(.1,1e5)
            ax.set(title=ray+', '+field,xlabel='t / M',ylabel='D / (5 x spread)')
            ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            if i==0 and j==0:ax.legend()
    for ext in ('png','pdf'):fig.savefig(target/f't14d-significance.{ext}',dpi=600)
    plt.close(fig)
    rows=d.read('t14d-constraint-orders-history.csv')
    fig,axes=plt.subplots(2,4,figsize=(9,4.7),layout='constrained')
    for i,ray in enumerate(('axis','diagonal')):
        for j,field in enumerate(d.CONSTRAINTS):
            ax=axes[i,j]
            for kind,style in (('peak','--'),('RMS','-')):
                for pair,color in (('12,13','C0'),('13,14','C1')):
                    z=[x for x in rows if x['ray']==ray and x['field']==field and x['norm']==kind and x['pair']==pair and int(x['step'])>0]
                    ax.plot([f(x,'time_M') for x in z],[f(x,'p') for x in z],style,color=color,label=pair+' '+kind)
                    bad=[x for x in z if x['qualified']=='False']
                    ax.scatter([f(x,'time_M') for x in bad],[f(x,'p') for x in bad],s=15,facecolors='none',edgecolors=color)
            ax.axhline(4,color='k',ls=':');ax.set(ylim=(2,13),title=ray+', '+field,xlabel='t / M',ylabel='Observed p')
            ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            if i==0 and j==0:ax.legend()
    for ext in ('png','pdf'):fig.savefig(target/f't14d-constraint-orders.{ext}',dpi=600)
    plt.close(fig)
    print('Two scientific figure pairs written.',flush=True)

def report():
    v=json.loads((HERE/'t14d-verdict.json').read_text());s=d.read('t14d-screen.csv');history=d.read('t14d-screen-history.csv')
    u=d.read('t14d-unqualified.csv');orders=d.read('t14d-constraint-orders.csv');rr=resources()
    compact=[]
    for ray in ('axis','diagonal'):
        for field in d.FIELDS:
            z=[next(x for x in s if x['ray']==ray and x['field']==field and x['norm']==kind) for kind in ('peak','RMS')]
            compact.append([ray,field,pairs(s,ray,field,'endpoint_rho'),
                ' / '.join(f"{min(f(x,'endpoint_margin12_13'),f(x,'endpoint_margin13_14')):.2f}" for x in z),
                ' / '.join(x['qualified_times']+'/'+x['contracting_times']+'/'+x['uncertainty_dominated_times']+'/'+x['numerical_floor_times'] for x in z)])
    unqualified=table(['Step','Time / M','Ray','Field','Norm','Status','Interpolation-failed pair','Floor-failed pair','Margin12,13','Margin13,14'],
        [[x['step'],f"{f(x,'time_M'):.12g}",x['ray'],x['field'],x['norm'],x['status'],x['significance_failed_pairs'] or '—',x['floor_failed_pairs'] or '—',
            f"{f(x,'margin12_13'):.6g}",f"{f(x,'margin13_14'):.6g}"] for x in sorted(u,key=lambda x:(int(x['step']),x['ray'],x['field'],x['norm']))])
    ot=[]
    for ray in ('axis','diagonal'):
        for field in d.CONSTRAINTS:
            vals=[]
            for column in ('endpoint_p','history_sup_p'):
                for kind in ('peak','RMS'):
                    vals.append(' / '.join(f"{f(next(x for x in orders if x['ray']==ray and x['field']==field and x['norm']==kind and x['pair']==pair),column):.3f}" for pair in ('12,13','13,14')))
            ot.append([ray,field,*vals])
    text=f'''\n## T14d — predeclared I8/I10 analysis-only follow-up

**Endpoint PASS; full registered history {v['verdict'].upper()}.** The operator predeclared one change before this result: central I8 and spread |I8-I10| replace I6 and |I6-I8|. Rungs, all 57 clocks, W=[0.00075,0.0025] M, 258 samples, peak/radial RMS norms, >5x significance, temporal <=0.2 of both spatial differences, initial numerical-floor treatment and the decision table remain frozen. No further operator variants are tried. [Registration](t14d-registration.json) pins all required stored inputs, exact common clocks, the history-start rule and Stage D exclusion radius before re-screening. Baseline HEAD is ba95217 on t14-stagec.

All four completed native recorder streams were re-extracted, including T13 alpha_K max12; the stopped control is excluded. Frozen T13Replay --state evaluates the same native current-state metric connection and constraints on 7x7 stencils. I10 uses tensor ten-point nodes, current-state reflection parity, and the same anchor-subtracted sampling routine. Stored support is complete at all 57 clocks. Every re-extracted I8 value, including constraint components, is bit-identical to the retained T13/T14 I8 cache. No static EMS file is opened, no state is evolved, and no T13/T14 artifact is overwritten.

Common endpoint remains 0.001993815104166667 M; the original control's step 448 supplies that physical time. Maximum clock defect is 4.33681e-19 M. I10 plus native derivatives requires a conservative 8h analysis support collar instead of I8's 7h. This changes only the diagnostic support accounting; all cells remain in the stored finest recorder interior. With the same 1.1 speed bound, the limiting L14 axis/diagonal arrival margins beyond 0.002 M remain approximately 0.001682824 / 0.002348491 M. Characteristic separation does not assert compact FD/KO numerical support.

### Screen and remaining entries

All endpoint pairs below are peak / RMS. Margin is min(D12,13/(5E12,13), D13,14/(5E13,14)), where E is the sum of each pair's |I8-I10| norm spreads; >1 qualifies. History counts are **qualified/contracting/uncertainty/floor**, each out of 56 positive clocks. Initial status is reported separately below.

{table(['Ray','Field','Endpoint rho_D peak / RMS','Minimum margin peak / RMS','History counts peak / RMS'],compact)}

The full screen retains **{v['positive_unqualified']} positive sampling-unqualified entries**, compared with T14's 53. All are on the axis. There are no positive numerical-floor, temporal-unqualified, marginal or significant noncontracting entries. Every qualified history entry contracts, with maximum rho_D **{v['max_qualified_rho']:.9g}**. The largest finest temporal fraction is **{v['max_temporal_fraction']:.9g}**, far below 0.2. The endpoint minimum sampling margin is **{min(min(f(x,'endpoint_margin12_13'),f(x,'endpoint_margin13_14')) for x in s):.6g}**. Recomputed disturbance amplitudes retain T14's nonzero stabilizing trend; [amplitude history](t14d-amplitude-history.csv) gives every rung at every common clock. No amplitude gate was introduced.

The previously long axis-lapse peak band shrinks to steps 1–9 and 13–15; steps 10–12 and 16–56 qualify. I8/I10 need not improve every early entry: the Gamma peak exceptions move to steps 2–4, and axis metric-Gamma peak remains unresolved at step 4. The complete table below lists every unqualified entry, including all 16 initial entries: eight Gamma/metric-Gamma numerical-floor rows and eight shift/lapse sampling-unqualified rows. The numerical t=0 baseline is not treated as bitwise zero. Pair labels refer to max12,max13 and max13,max14; native steps are 2/4/8 times the listed max12 step for max13/max14/control.

{unqualified}

[Every common-time screen](t14d-screen-history.csv) retains both 1x and **5x** rho sensitivity intervals. The 5x upper ratio is infinite whenever D12,13<=5E12,13; finite 1x bounds do not qualify these rows. Numerical-floor ratios remain undefined. Bounds describe empirical diagnostic sensitivity, not rigorous continuum-error intervals. The original significance/status gate is unchanged; no extra ratio-bound admission threshold is added.

**Verdict: {v['verdict'].upper()}.** The residual positive-time uncertainty prevents the full registered-history pass. Stop tightening interpolation. No automatic ladder extension, method change or production admission follows. The independently requested Stage D draft is supplied below for the operator's admission decision, irrespective of this verdict.

![Fixed 5x sampling margins across all fields and clocks](figures/t14d-significance.png)

### Constraint-to-zero orders on W

For c=C_Gamma,Ham,Mom,GaussE, observed p is log2(||c_h||/||c_h/2||) and log2(||c_h/2||/||c_h/4||). These are direct constraint residual orders, distinct from evolved-field self-difference orders. Mom is the norm of its sampled Cartesian components; C_Gamma is signed longitudinal evolved-minus-native metric Gamma. The same I8/I10 diagnostic operator and radial peak/RMS norms are used. **Design order is four**, with no lower regularity order assumed. Each cell below is max12→max13 / max13→max14. History sup compares each rung's supremum norm over the same 57 clocks; it is not the maximum instantaneous order.

{table(['Ray','Constraint','Endpoint peak p','Endpoint RMS p','History-sup peak p','History-sup RMS p'],ot)}

All endpoint constraint norm pairs exceed 5x their own interpolation spreads. The histories do not exhibit a uniform asymptotic fourth-order regime: axis coarse-pair endpoint orders are about 2.53–3.55, while fine-pair orders often exceed six. Large orders are measured cancellation/transient ratios, not proof of a higher design order. [All 1824 instantaneous pair orders](t14d-constraint-orders-history.csv) retain initial/positive clocks, raw norms, sampling qualification and 5x order sensitivity bounds; [summary](t14d-constraint-orders.csv) retains initial, endpoint, history-sup and positive-time min/max orders and qualified counts. Initial zero C_Gamma orders are undefined. Some positive constraint rows are sampling-limited and remain visibly marked, without a roundoff waiver.

![Constraint order histories; hollow markers are sampling-limited](figures/t14d-constraint-orders.png)

Blue curves show max12→max13 and green curves max13→max14; dashed curves are peak norms and solid curves RMS. Hollow markers retain sampling-unqualified orders.

### Stage D — for the operator's admission decision

The complete [unexecuted design](t14d-stageD-design.md) specifies the fresh E 3/2 global chain with common exp-0020 faces plus levels 13–14 on all rungs, the fixed R>=0.005 M exclusion, all inherited masks and the boundary issue, time-resolved field/constraint evidence, unchanged numerical-t=0 horizon budgets, independent sphere charges and a finest dt/2 control. The history-start rule evaluates to **step {v['first_all_significant_step']}, t={v['first_all_significant_time_M']:.17g} M**. Steps {v['later_sampling_exceptions']} later fail: that first significant clock does not retroactively pass T14d or justify dropping any row. The exclusion radius was fixed before this result and the history-start rule is frozen for a prospective fresh qualifying run.

The extra source levels multiply the valid-cell subcycle work by approximately four; the full nominal chain plus finest dt/2 control is priced in the draft using measured exp-0020 256/162/104 s coarse-step medians on one node per leg, 32 MPI x 4 OpenMP. The estimate is a cell-step model, not a measured new cluster rate. [Fixed masks](t14d-stageD-masks.csv), [rungs and cost model](t14d-stageD-cost.csv). No run, cluster contact or finder invocation is authorized or performed.

### Verification and resources

Final measured analysis peak RSS is **{max(int(x['peak_rss_bytes']) for x in rr):,} bytes**, sampled active-tree peak **{max(int(x['sampled_tree_bytes']) for x in rr):,} bytes**. Numerical processes run serially, OMP=2 and BLAS=1, under the existing 3 GB per-process/tree gates. The sandbox timer's kern.clockrate failure is handled by the existing wait4/libproc measurement; actual child return codes are retained. An initial register invocation used a relative script path under the wrapper's changed directory and exited 2; the corrected absolute-path registration completed before any extraction. Both records remain in [resources](t14d-resources.csv). Native payloads removed were only T14d-owned intermediates, with exact hashes retained; native q and I8/I10 caches remain under /private/tmp/ems-t14d. [Checks](t14d-verification.csv), [manifest](t14d-manifest.txt).

Reproduce serially with the measured t14-run.py wrapper, T14_OUTPUT_ROOT=/private/tmp/ems-t14d and T14_DISK_CAP_BYTES=8000000000, using absolute script paths: t14d-analyze.py register; extract max12|max13|max14|max14-half; screen; t14d-report.py design; figures; check; report; manifest. Never reuse the stopped control or overwrite T13/T14 outputs. No production source, evolved state, gauge, KO, equation, transfer, precision, finder or initial-data reader was changed. No SSH or commit.
'''
    path=HERE/'README.md';old=path.read_text().split('\n## T14d —')[0];path.write_text(old+text)
    print('README T14d written.',flush=True)

def check():
    rows=[]
    def yes(name,value,condition):
        assert condition,(name,value);rows.append(dict(check=name,value=value,result='PASS'))
    for n in (8,10):
        z=np.array([-.41,.13,.5,1.87]);start=np.floor(z).astype(int)-n//2+1;w=d.a.weights(z,start,n)
        for power in range(n):
            actual=np.sum(w*((start[None,:]+np.arange(n)[:,None]-z[None,:])/n)**power,axis=0)
            target=np.ones(len(z)) if power==0 else np.zeros(len(z))
            assert np.max(abs(actual-target))<2e-15,(n,power,actual)
    yes('I8/I10 centred polynomial moments',20,True)
    identities=[x for leg in d.LEGS for x in d.read(f't14d-identity-{leg}.csv')]
    yes('retained I8 bit identity',sum(int(x['Float64_values']) for x in identities),all(int(x['I8_bit_mismatches'])==0 for x in identities))
    hist=d.read('t14d-screen-history.csv');s=d.read('t14d-screen.csv');u=d.read('t14d-unqualified.csv')
    yes('57 x 16 screen rows',len(hist),len(hist)==912 and len(s)==16)
    yes('unqualified table retains each row',len(u),len(u)==sum(x['status'] in ('uncertainty-dominated','numerical-floor','temporal-unqualified') for x in hist))
    yes('initial floor versus sampling counts',16,sum(x['step']=='0' and x['status']=='numerical-floor' for x in u)==8 and sum(x['step']=='0' and x['status']=='uncertainty-dominated' for x in u)==8)
    yes('frozen endpoint contraction',max(f(x,'endpoint_rho') for x in s),all(x['endpoint_status']=='contracting' for x in s))
    yes('5x denominator bounds unbounded when required',0,all(math.isinf(f(x,'rho_upper_5x')) for x in hist if x['status']!='numerical-floor' and f(x,'D12_13')<=5*f(x,'interpolation12_13')))
    orders=d.read('t14d-constraint-orders-history.csv')
    yes('32 x 57 direct constraint orders',len(orders),len(orders)==1824)
    for x in orders:
        c=f(x,'coarse_norm');v=f(x,'fine_norm')
        if c and v:assert abs(f(x,'p')-math.log2(c/v))<1e-14
        else:assert math.isnan(f(x,'p'))
    yes('direct orders reproduce norm ratios',len(orders),True)
    rr=resources();yes('measured 3 GB process/tree caps',max(int(x['peak_rss_bytes']) for x in rr),
        all(int(x['peak_rss_bytes'])<3e9 and int(x['sampled_tree_bytes'])<3e9 and not x['gate_reason'] for x in rr))
    yes('extraction markers',4,all((ROOT/'t14d-extract'/leg/'done.exit').read_text().strip()=='0' for leg in d.LEGS))
    yes('unexecuted Stage D draft',1,(HERE/'t14d-stageD-design.md').is_file())
    d.save('t14d-verification.csv',rows);print(json.dumps(rows,indent=2),flush=True)

def manifest():
    rr=resources();d.save('t14d-resources.csv',rr)
    paths=list(HERE.glob('t14d-*'))+list((HERE/'figures').glob('t14d-*'))+[HERE/'README.md']
    paths+=list(ROOT.rglob('*'))
    reg=json.loads((HERE/'t14d-registration.json').read_text())
    lines=['T14d at ba95217; I8 and |I8-I10| only; no further variants, evolution, SSH or commit.',
        'All numerical processes measured; process/tree cap 3 GB; serial OMP=2 BLAS=1.',
        'Input hashes below were pinned before extraction and rechecked by the manifest.']
    for x in reg['inputs']:
        p=Path(x['path']);assert p.stat().st_size==x['bytes'] and d.digest(p)==x['sha256'],p
        lines.append(f"{x['sha256']}  {x['bytes']}  {p}")
    paths += [HERE/'t14-run.py',HERE/'t13-run.py',HERE/'t14d-registration.json',HERE/'t8-run-audit.csv',HERE/'t14-census.csv',HERE/'t14-native-lobes.csv',HERE/'COMMIT-MANIFEST-T14.txt']
    submission=Path('/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0020/submissions/exp-0020')
    paths += [submission/'submit-contract.md',submission/'audit.py']+[submission/f'params-E-{leg}.txt' for leg in ('low','mid','high')]
    for p in sorted(set(paths)):
        if not p.is_file() or p.name=='t14d-manifest.txt':continue
        if p.name.startswith(os.environ.get('T13_CURRENT_MEASURE','NONE')+'.'):continue
        lines.append(f'{d.digest(p)}  {p.stat().st_size}  {p}')
    lines += [f"RESOURCE {x['process']} peak_RSS_bytes={x['peak_rss_bytes']} tree_bytes={x['sampled_tree_bytes']} returncode={x['returncode']} record={x['record']}" for x in rr]
    (HERE/'t14d-manifest.txt').write_text('\n'.join(lines)+'\n');print('T14d manifest written.',flush=True)

if __name__=='__main__':
    {'design':design,'figures':figures,'report':report,'check':check,'manifest':manifest}[sys.argv[1]]()
