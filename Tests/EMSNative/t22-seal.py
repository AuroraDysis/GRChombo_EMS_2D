#!/usr/bin/env python3
"""Seal completed T22 evidence or mark pending evidence without claiming a pass."""
import csv, hashlib, json, resource, subprocess, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];ROOT=Path('/private/tmp/ems-t22')
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
resources=[]
for p in sorted(ROOT.rglob('*.resources.json')):
    a=json.loads(p.read_text());resources.append(dict(case=a['process'],receipt=str(p),returncode=a['returncode'],gate_reason=a['gate_reason'],wall_seconds=a['wall_seconds'],peak_RSS_bytes=a['peak_rss_bytes'],own_child_returncode=a.get('child_measurement',{}).get('returncode','')))
if resources:
    with (HERE/'t22-resources.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=resources[0]);w.writeheader();w.writerows(resources)
quick=HERE/'t22-quick-qualification.json';equilibrium=HERE/'t22-equilibrium-qualification.json';launch=HERE/'t22-launch-qualification.json'
status=[]
for file,title in [(quick,'Default/generic/restart audits'),(equilibrium,'E-high native equilibrium audit'),(launch,'E-mid launch pair')]:
    if file.exists():
        a=json.loads(file.read_text());status.append(f'{title}: **{a["status"]}**, measured native/build RSS ≤{a["peak_RSS_bytes"]/1e9:.3f} GB. See [{file.name}]({file.name}).')
    else:status.append(f'{title}: **PENDING**; no scientific pass is inferred from a launcher exit.')
if equilibrium.exists():
    a=json.loads(equilibrium.read_text())
    status.append(f'**Explicit factual correction:** the consult premise that the setter supplies exact K and Theta zeros is false for K. Theta is exactly zero; K is a Float64 metric contraction, with valid-cell peak {a["K_measured_peak"]:.8e} M^-1 and maximum {a["K_scaled_eps_peak"]:.6f} eps times the local absolute contraction scale, within the explicitly declared 64-eps bound. The former exact-zero gate is replaced by this measured roundoff gate. The all-valid-cell MPG pre-KO lapse rate peaks at {a["preKO_MPG_lapse_measured_peak"]:.8e} M^-1; the shift rate remains exactly zero. [The equilibrium report](t22-equilibrium-report.md) gives every level, ghosts, locations and the complete geometric/matter/KO split. The frozen activity screen, masks, clocks and uncertainty rules remain unchanged.')
if launch.exists():
    a=json.loads(launch.read_text())
    status.append('**Short launch facts:** at t=0.0019938151 M, new/old Gamma disturbance ratios (peak/RMS) are axis 0.274335/0.138185 and diagonal 0.376254/0.175699. All endpoint reads qualify with fivefold sampling and same-gauge dt/2 margins. Early history is mixed: the new gauge has resolved larger intervals; axis peak steps 2–4 remain uncertainty-dominated. No recorded nonfinite values or chi/lapse floor activations occur. [The launch report](t22-launch-report.md), [full ratios](t22-launch-ratios.csv), [unqualified entries](t22-launch-unqualified.csv) and [norm histories](t22-launch-history.pdf) retain the limits. This is no long-time decision or convergence claim.')
status.append('The original evidence queue stopped at its exact-zero assertion (exit 1); its failed receipt is retained. The remaining-job queue marker is `/private/tmp/ems-t22/continuation/done.exit`; each native job has its own receipt and marker. No completed initialization is repeated. The 10.5 M control is designed only, not run or admitted.')
block='\n\n'.join(status)
path=HERE/'t22-gauge-control-design.md';text=path.read_text();begin='<!-- T22 qualification -->';end='<!-- /T22 qualification -->'
left,tail=text.split(begin,1);_,right=tail.split(end,1)
original=(ROOT/'design-frozen.md').read_text();ol,ot=original.split(begin,1);_,orr=ot.split(end,1)
assert left==ol and right==orr,'registered design changed outside qualification note'
path.write_text(left+begin+'\n'+block+'\n'+end+right)
readme=HERE/'README.md';text=readme.read_text();title='## T22 — opt-in nonadvective differential Gamma-driver control'
suffix=''
if title in text:
    start=text.index(title);stop=text.find('\n## ',start+len(title))
    if stop>=0:suffix=text[stop:]
    text=text[:start].rstrip()
note=f'''{title}

Consult 9b and the operator's ruling prioritize a gauge control; exp-0023 remains **FAIL**, and old-gauge level 15 remains superseded. [The design](t22-gauge-control-design.md) was frozen before new numerical output. [G parameters](t22-params-E-high.txt), [same-gauge dt/2 parameters](t22-params-E-high-dt2.txt), [fixed masks](t22-masks.csv) and [gauge-labelled T21 recording](t22-recording.json) form the submission design, with four-step whole-node calibration required before production. H is an engineering comparator; no spatial convergence or binary-gauge recommendation is made.

`ems_gauge` selects `experimental` (default) or `moving_puncture` at initialization and native RHS construction. The new differential driver receives the complete pre-KO Gamma RHS including EMS momentum; ordinary KO(B) follows, without extra KO(Gamma). It initializes B=0 while retaining delivered lapse and shift. The historical call, initialization and fixed legacy coefficients are preserved. Run records and checkpoint headers identify the gauge/driver semantics; a cross-gauge restart is refused, while an untagged old checkpoint belongs to experimental.

{block}

[The generic native harness](t22-audit.cpp), [E-high initialization harness](t22-initial-audit.cpp), [bounded build/plan preparer](t22-prepare.py) and [receipt-first verifier](t22-verify.py) expose all tests. [Exact scoped witnesses](t22-gauge-witness.json) and [asymptotic CFL checks](t22-asymptotic-CFL.csv) give directional Courant numbers 0.35355 (lapse), 0.28868 (driver) and 0.25 (light), at dt/h=0.25 on every level. Full curved CCZ4 hyperbolicity is not certified by these wave proxies. Unit-speed Sommerfeld remains unchanged; the asymptotic boundary return to R=20 is about 447 M, beyond the 10.5 M control. No run parameter is silently altered.

[The K/Theta census](t22-KTheta-census.csv) and [all-valid-cell lapse rates](t22-initial-lapse-rates.csv) come only from completed numerical t=0 checkpoints. [The native storage check](t22-native-storage-check.csv) remains bitwise: its former one-component flag compared an overwritten upper A12 tensor alias with the last-written lower alias, rather than the 28 evolved doubles. The corrected test agrees exactly on all stored components. [The alias ledger](t22-tensor-aliases.csv) retains every overwritten-value discrepancy; no production equation is changed by this test-harness repair.

The predicted G evolution cost is 12.54 node-hours at H's measured 418 s/step ×108; G dt/2 predicts 25.08. Use 72+36 step segments for G and three 72-step segments for dt/2 under 12 h limits, recalibrating actual rates and reserving measured initialization, finder and recording overhead. [Resources](t22-resources.csv) and [manifest](t22-manifest.txt) retain local process measurements. No commit, SSH, cluster work, positive-time static read or new production evolution is performed.
'''
readme.write_text(text+'\n\n'+note+suffix)
# The T21 packet has archived pre-T22 source inputs. Refresh its README hash
# after adding this section; no T21 numerical proof or frozen input is changed.
subprocess.run([sys.executable,str(HERE/'t21-manifest.py')],check=True)
registration=json.loads((HERE/'t22-registration.json').read_text())
assert sha(ROOT/'design-frozen.md')==registration['readout_files']['t22-gauge-control-design.md']
for name,digest in registration['readout_files'].items():
    if name!='t22-gauge-control-design.md':assert sha(HERE/name)==digest,name+' changed after registration'
own=HERE/'t22-seal-resources.json';own.write_text(json.dumps(dict(peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,process_cap_bytes=6000000000,threads=1),indent=2)+'\n')
paths={p for p in HERE.glob('t22-*') if p.is_file() and p.name!='t22-manifest.txt'}
paths.update([readme,HERE/'t21-manifest.txt',ROOT/'design-frozen.md'])
for f in ['Examples/EMS/EMSGaugeSelection.hpp','Examples/EMS/EMSBH2DLevel.hpp','Examples/EMS/EMSBH2DLevel.cpp','Examples/EMS/EMSBH2DMovingGauge.cpp','Examples/EMS/EMSBH2DRHS.impl.hpp','Examples/EMS/Main_EMSBH2DBH.cpp','Source/Cartoon/CCZ4Cartoon.hpp','Source/Cartoon/CCZ4Cartoon.impl.hpp','Source/CCZ4/MovingPunctureGauge.hpp','Source/CCZ4/ExperimentalGauge.hpp']:
    paths.add(REPO/f)
for f in ['/Users/auroradysis/Workspace/EMS-deps/worktrees/consult-9b-reply.md','/Users/auroradysis/Workspace/EMS/artifacts/echo/E.trumpet','/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0023/submissions/exp-0023/production/params-E-high.txt','/Users/auroradysis/Workspace/EMS/state/threads/ems-spectral-solver/runs/exp-0023/submissions/exp-0023/submit-contract.md']:
    paths.add(Path(f))
for p in ROOT.rglob('*'):
    if p.is_file() and (p.name in ('done.exit','plan.json','probe.receipt.json','params.txt','ems-gauge.txt','run.log','probe.log') or p.suffix in ('.bin','.hdf5','.xz','.csv','.ex','.npz','.sha256') or p.name in ('source-snapshots.json','reused-objects.json','build-commands.json')):
        if p.name=='run.log' and not (p.parent/'done.exit').exists():continue
        if p.suffix in ('.bin','.hdf5','.xz','.csv') and not (p.parent/'done.exit').exists() and p.parent.name!='trace-audit':continue
        paths.add(p)
paths.update(ROOT.rglob('*.resources.json'));paths.update(ROOT.rglob('*.child.json'))
with (HERE/'t22-manifest.txt').open('w') as f:
    f.write('# T22 uncommitted gauge control; frozen design, own receipts; no production admission\n# SHA256 bytes absolute_path\n')
    for p in sorted(paths):f.write(f'{sha(p)}  {p.stat().st_size}  {p}\n')
print(json.dumps(dict(status=status,sealed_files=len(paths),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),indent=2))
