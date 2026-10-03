#!/usr/bin/env python3
"""Seal the local T21 packet; update only its qualification paragraph and README section."""
import csv, hashlib, json, resource
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t21')
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
q=HERE/'t21-qualification.json';qualified=False
if q.exists():
    a=json.loads(q.read_text())
    receipt=json.loads((ROOT/'qualification/qualification.resources.json').read_text())
    qualified=a['status']=='RECORDING_BIT_IDENTITY_PASS' and receipt['returncode']==0 and not receipt['gate_reason']
    if not qualified:raise RuntimeError('qualification file is not supported by its own receipt')
    a['qualifier_measured_peak_RSS_bytes']=receipt['peak_rss_bytes']
    q.write_text(json.dumps(a,indent=2)+'\n')
    rows=[r for r in csv.DictReader((HERE/'t21-local-resources.csv').open()) if r['case'] not in ('qualification','seal','queue')]
    rows.append(dict(case='qualification',returncode=receipt['returncode'],gate_reason=receipt['gate_reason'],peak_RSS_bytes=receipt['peak_rss_bytes'],wall_seconds=receipt['wall_seconds']))
    for name,directory in [('seal',ROOT/'seal'),('queue',ROOT/'pipeline')]:
        p=directory/(name+'.resources.json')
        if p.exists():
            own=json.loads(p.read_text())
            if own['returncode']!=0 or own['gate_reason']:raise RuntimeError(f'{p}: unsuccessful own receipt')
            rows.append(dict(case=name,returncode=own['returncode'],gate_reason=own['gate_reason'],peak_RSS_bytes=own['peak_rss_bytes'],wall_seconds=own['wall_seconds']))
    with (HERE/'t21-local-resources.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    status=(f'**Local recording qualification PASS:** {a["files"]} plot/checkpoint comparisons, '
            f'{a["checked_values"]:,} Float64 values, zero bit mismatches; omitted/off/on recorder paths agree. '
            f'Numerical t=0 and positive-time checkpoint replays succeed with the static path absent. '
            f'Peak measured RSS across native/build cases {a["peak_RSS_bytes"]/1e9:.3f} GB; '
            f'qualifier {receipt["peak_rss_bytes"]/1e9:.3f} GB. '
            'See [t21-qualification.json](t21-qualification.json), [identity comparisons](t21-recording-identity.csv) '
            'and [own resource receipts](t21-local-resources.csv). MPI/target-node preflight remains required.')
else:
    status=('**Qualification pending:** the bounded local recording-only build/identity pipeline has not yet supplied its final receipts. '
            'No new-hook evolved-bit identity or production admission is claimed. '
            'Marker: `/private/tmp/ems-t21/pipeline/done.exit`; verify each child receipt before accepting it.')
path=HERE/'t21-source-level-design.md';text=path.read_text();begin='<!-- T21 qualification -->';end='<!-- /T21 qualification -->'
left,tail=text.split(begin,1);_,right=tail.split(end,1);path.write_text(left+begin+'\n'+status+'\n'+end+right)
readme=HERE/'README.md';text=readme.read_text();title='## T21 — gauge-independent recording; old-gauge level 15 superseded'
suffix=''
if title in text:
    start=text.index(title);stop=text.find('\n## ',start+len(title))
    if stop>=0:suffix=text[stop:]
    text=text[:start].rstrip()
note=f'''{title}

The operator's ruling and consult 9b supersede the old-gauge level-15 submission plan; it was not executed or failed, and its P parameters/census were not completed. Exp-0023 remains **FAIL**. [The design](t21-source-level-design.md) now delivers recording for the next fresh E-high gauge-control card. No gauge selector, gauge-equation repair or level-15 cost model is implemented here.

[The opt-in recorder overlay](t21-recorder.patch), [native checkpoint replay](t21-native-replay.cpp), [extractor](t21-record.py) and [configuration](t21-recording.json) capture thin axis/equator/diagonal profiles, complete native pre-KO/KO/total RHS near the puncture, current-numerical relational profiles and fixed-radius histories at native coarse cadence. B1/B2 carry explicit gauge/driver-semantic labels and are not cross-gauge observables. First RK states require fresh recording; the retained M/H checkpoints recover all other requested snapshot/cadence quantities on the compute node. No T13 early-stop recorder is used.

{status}

[Codec/sampling checks](t21-record-test.json) and [exact geometry witnesses](t21-geometry-witness.json) pass. Angular area-density radius and mass have the spherical areal/Misner-Sharp interpretation only when angular disagreement is qualified; ray disagreement, non-monotonic relations, invalid geometry, KO dependence and interpolation margins remain visible. [Signed histories](t21-signed-history-check.csv) reproduce A11 correlations −0.106/−0.220/−0.293 and the K/lapse anti-correlation. [Recording storage](t21-recording-storage.json) predicts about 6.8 GB of raw streams plus uncompressed profiles for an E-high history; reserve 32 GB separately from retained checkpoints/plots. [The manifest](t21-manifest.txt) seals source inputs and this packet. No SSH, cluster work, commit, static subtraction or positive-time static read is performed.
'''
readme.write_text(text+'\n\n'+note+suffix)
paths=set(p for p in HERE.glob('t21-*') if p.is_file() and p.name!='t21-manifest.txt');paths.add(readme)
own=HERE/'t21-manifest-resources.json'
own.write_text(json.dumps(dict(peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,process_cap_bytes=6000000000,threads=1),indent=2)+'\n')
paths.add(own)
for r in json.loads((HERE/'t21-registration.json').read_text())['inputs']:paths.add(Path(r['path']))
paths.update(p for p in ROOT.rglob('*.resources.json') if p.is_file())
paths.update(p for p in ROOT.rglob('*.child.json') if p.is_file())
paths.update(p for p in ROOT.rglob('done.exit') if p.is_file())
for s in ('build-commands.json','reused-objects.json','original.ex','recording.ex','replay.ex'):
    if (ROOT/s).exists():paths.add(ROOT/s)
for s in ('qualification/t21-stage-summary.csv','qualification/t21-stage-global-summary.csv'):
    if (ROOT/s).exists():paths.add(ROOT/s)
with (HERE/'t21-manifest.txt').open('w') as f:
    f.write('# T21 recording packet; old-gauge level15 superseded; no commit\n# SHA256 bytes absolute_path\n')
    for p in sorted(paths):f.write(f'{sha(p)}  {p.stat().st_size}  {p}\n')
print(json.dumps(dict(qualified=qualified,sealed_files=len(paths),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),indent=2))
