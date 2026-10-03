#!/usr/bin/env python3
"""Summarize completed native launch receipts and fivefold-qualified histories."""
import csv,json,resource
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
def rows(name):return list(csv.DictReader((HERE/name).open()))
def save(name,data):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
def table(headers,data):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(x) for x in r)+' |' for r in data])
a=rows('t22-launch-ratios.csv');h=rows('t22-launch-history.csv');q=json.loads((HERE/'t22-launch-qualification.json').read_text());assert q['status']=='LAUNCH_PAIR_COMPLETE'
endpoint=rows('t22-launch-endpoint.csv');history=[]
for ray in ('axis','diagonal'):
    for norm in ('peak','RMS'):
        b=[r for r in a if r['ray']==ray and r['norm']==norm and int(r['step'])>0]
        history.append(dict(ray=ray,norm=norm,positive_clocks=len(b),qualified=sum(r['qualified']=='True' for r in b),unqualified_steps=' '.join(r['step'] for r in b if r['qualified']!='True'),new_over_old_min=min(float(r['new_over_old']) for r in b),new_over_old_max=max(float(r['new_over_old']) for r in b),resolved_new_larger_steps=' '.join(r['step'] for r in b if float(r['reduction_upper'])<0),joint_temporal_over_old_peak=max(float(r['temporal_joint_margin'])/float(r['old_Gamma']) for r in b),old_temporal_over_old_peak=max(float(r['old_temporal_difference'])/float(r['old_Gamma']) for r in b),new_temporal_over_new_peak=max(float(r['new_temporal_difference'])/float(r['new_Gamma']) for r in b)))
save('t22-launch-history-summary.csv',history)
unqualified=[dict(r,failed_signal='initial zero disturbance' if int(r['step'])==0 else ' '.join(g for g in ('old','new') if float(r[g+'_significance_margin'])<=0)) for r in a if r['qualified']!='True']
save('t22-launch-unqualified.csv',unqualified)
et=table(['ray','norm','old disturbance','new disturbance','new / old [bound]','5x joint sampling','joint dt/2 difference'],[[r['ray'],r['norm'],f"{float(r['old_Gamma']):.8e}",f"{float(r['new_Gamma']):.8e}",f"{float(r['new_over_old']):.6f} [{float(r['ratio_lower']):.6f}, {float(r['ratio_upper']):.6f}]",f"{float(r['fivefold_joint_sampling_margin']):.8e}",f"{float(r['temporal_joint_margin']):.8e}"] for r in endpoint])
ht=table(['ray','norm','qualified / positive clocks','unqualified positive steps','ratio range','max joint dt/2 / old'],[[r['ray'],r['norm'],f"{r['qualified']} / {r['positive_clocks']}",r['unqualified_steps'] or 'none',f"{r['new_over_old_min']:.6f}–{r['new_over_old_max']:.6f}",f"{r['joint_temporal_over_old_peak']:.6f}"] for r in history])
ut=table(['ray','norm','step','time M','failed signal','old significance margin','new significance margin'],[[r['ray'],r['norm'],r['step'],f"{float(r['time_M']):.10e}",r['failed_signal'],f"{float(r['old_significance_margin']):.8e}",f"{float(r['new_significance_margin']):.8e}"] for r in unqualified if int(r['step'])>0])
rt=table(['case','own child exit','wall seconds','peak RSS bytes'],[[r['case'],r['returncode'],f"{float(r['wall_seconds']):.3f}",r['peak_RSS_bytes']] for r in rows('t22-launch-resources.csv')])
text=f'''# T22 short E-mid launch: qualified endpoint reductions, mixed early history

All four previously unrun trajectories completed their own native clean stops and passed receipt-first analysis. This is a local early-window audit, **not** the registered 10.5 M decision or a convergence statement. No completed initialization/evolution was repeated. The old and new packages keep identical delivered physical initial data, hierarchy, KO, transfers and t2 static guard; only gauge package and driver initialization differ. B variables are not compared across gauges.

The hierarchy has levels 0–12, h0=7/12 M, finest h=7/(12*4096)=0.00014241536458333333 M and dt/h=0.25 (0.125 for each same-gauge control). The nominal native stop is t=56*dt=0.001993815104166667 M, within the registered 0.002 M early window. The nominal trajectories have 57 saved level-12 clocks including t=0; dt/2 has 113, compared at every second clock. No unsynchronized plot/checkpoint is written. The native recorder is the qualified serial T13 stopping audit, not production recording.

Signed Gamma_n = (Gamma1,Gamma2).n is sampled on the axis and diagonal over the frozen radius window 0.00075<=R<=0.0025 M. The disturbance is the current numerical value minus its own stored numerical t=0 value. No static solution is read, subtracted or used as a reference. Peak is max absolute disturbance; RMS is its radial-integral RMS over the same window. I8 is central, with |I8-I10| at both current and initial clocks and the unchanged 128-eps Float64 floor. Each gauge must exceed five times its sampling/floor bound plus its own nominal-minus-dt/2 difference. Ratio bounds propagate these separately; the joint reduction margin is their sum.

## Endpoint at the common native stop

{et}

Every endpoint is qualified; even the upper ratio bounds are below one half. These are reductions of the measured Gamma disturbance in this window, not accuracy certificates. At the endpoint the joint temporal difference is smaller than the fivefold joint sampling margin in all four reads. [t22-launch-endpoint.csv](t22-launch-endpoint.csv) records the individual gauge temporal and sampling bounds, significance margins and reduction bounds.

## Entire common-time history

{ht}

The new gauge is **not uniformly quieter in the early history**. Both peak reads and RMS reads have resolved intervals where new exceeds old; all such steps are listed in [t22-launch-history-summary.csv](t22-launch-history-summary.csv). Axis peak reaches ratio 1.846309 and diagonal peak 2.508489 before their later decrease. In the first four steps, axis RMS ratios are 1.0723, 1.1520, 1.2594, 1.3872; diagonal RMS ratios are 1.0356, 1.1430, 1.1462, 1.1367. This short record has finite, bounded transients and no observed first-step runaway, but cannot diagnose a slow instability or decide the long-time activity screen.

{ut}

The three remaining positive-time unqualified reads are axis peak steps 2–4; their failed signal is reported explicitly. All four t=0 reads are initial zero-disturbance floor entries and are also retained in [t22-launch-unqualified.csv](t22-launch-unqualified.csv). None is relabelled as a pass. [t22-launch-ratios.csv](t22-launch-ratios.csv) has every one of the 228 nominal ray/norm/clock reads, with old/new sampling and temporal margins separately. [t22-launch-history.csv](t22-launch-history.csv) preserves both full nominal and half-step histories. Joint temporal/old disturbance is at most 3.11% over positive clocks; per-gauge maxima are in the summary, without suppressing early sampling-dominated entries.

## Floors, finite values and receipts

No chi/lapse floor activation and no nonfinite value was recorded at any inspected stage/level, including ghost copies. [t22-launch-floors.csv](t22-launch-floors.csv) lists every level and run. The smallest recorded chi is 5.49151383e-7 and lapse 1.90911861e-4, both above their 1e-12 floors. Every sampled state and I10 stencil is finite. Each run's own done marker, actual child returncode zero, empty gate reason, native stop clock, state count, matching temporal clocks, floor rows and `GRChombo finished.` log are checked. The sandboxed time-wrapper sysctl warning is not a failure.

{rt}

The maximum native launch RSS is {q['peak_RSS_bytes']:,} bytes ({q['peak_RSS_bytes']/1e9:.3f} GB), below the 6 GB per-process cap. Two OpenMP threads are used. The original evidence pipeline exit 1 and false exact-zero assertion remain retained; only the unfinished cases ran under `/private/tmp/ems-t22/continuation/plan.json`. The continuation's own receipt and marker confirm completion, but its exit alone is not used as a scientific pass.

[t22-launch-history.pdf](t22-launch-history.pdf) plots peak and RMS disturbance norms against physical time, with both same-gauge dt/2 histories. The registered activity reading remains frozen in [t22-gauge-control-design.md](t22-gauge-control-design.md); H is an engineering comparator. New-gauge 10.5 M and temporal-control legs, target compiler/MPI identity preflight and whole-node calibration remain for the submission worker/controller. No local 10.5 M evolution, cluster operation or commit is made.
'''
(HERE/'t22-launch-report.md').write_text(text)
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,2,figsize=(7.0,5.0),sharex=True,layout='constrained')
for iy,ray in enumerate(('axis','diagonal')):
    for ix,norm in enumerate(('peak','RMS')):
        ax=axs[iy,ix]
        for run,label,color,style in [('experimental','Experimental','#AA4465','-'),('experimental-half','Experimental dt/2','#AA4465',':'),('moving_puncture','Moving puncture','#2166AC','-'),('moving_puncture-half','Moving puncture dt/2','#2166AC',':')]:
            b=[r for r in h if r['run']==run and r['ray']==ray and r['norm']==norm]
            ax.plot([float(r['time_M']) for r in b],[float(r['Gamma_disturbance']) for r in b],label=label,color=color,linestyle=style,lw=1.15)
        ax.set_title(ray.capitalize()+' · '+norm);ax.set_ylabel('Gamma disturbance');ax.ticklabel_format(style='sci',axis='both',scilimits=(0,0));ax.grid(alpha=.2)
        if iy==1:ax.set_xlabel('t / M')
axs[0,0].legend(frameon=False)
fig.savefig(HERE/'t22-launch-history.pdf')
fig.savefig('/private/tmp/ems-t22/launch-history-preview.png',dpi=140);plt.close(fig)
(HERE/'t22-launch-report-resources.json').write_text(json.dumps(dict(peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,process_cap_bytes=6000000000),indent=2)+'\n')
print('T22_LAUNCH_REPORT_WRITTEN')
