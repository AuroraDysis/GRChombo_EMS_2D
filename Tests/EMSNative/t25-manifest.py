#!/usr/bin/env python3
"""Seal the small T25 packet and index external numerical evidence by hash.

No commits; no bulk copies. Lines contain role, SHA-256, size and absolute path.
Exclude this manifest's own self-referential hash and generated build folders.
"""
import hashlib,json,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t25')
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    audit=json.loads((HERE/'t25-audit.json').read_text());reg=json.loads((HERE/'t25-regression.json').read_text())
    assert audit['all_probes_verified'] and reg['status']=='PASS'
    early_reg=HERE/'t25-early-regression.json'
    if early_reg.exists():assert json.loads(early_reg.read_text())['status']=='PASS'
    files=[p for p in sorted(HERE.glob('t25-*')) if p.is_file() and p.name!='t25-manifest.txt']
    files.append(HERE/'README.md')
    lines=['T25 signed expansion map; no commit; no evolution; no remote work.',
        'Author production source archive d50743867e48bcb466a60ca16869b07fe1546c27.',
        'Pending NaN-abort source edits are excluded from every T25 build.',
        'Workspace HEAD (read, not changed): '+subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        'Original T25 peak map RSS: 1999863808 bytes; process cap 3000000000 bytes; serial, 1 thread.',
        'Original d32 binary initial-data peak RSS: 1267761152 bytes; process cap 6000000000 bytes.',
        'The regenerated single-hole input has header A_H=0.22451077798968849; card benchmark .2218 differs.',
        'Small files only below; raw angular arrays and checkpoints are external and are NOT commit candidates.',
        'FORMAT: role sha256 bytes absolute_path']
    early=ROOT/'early'
    if early.exists():
        receipts=[json.loads(p.read_text()) for p in early.rglob('*.resources.json')]
        if receipts:lines.append('T25 early continuation observed peak RSS: '+str(max(r['peak_rss_bytes'] for r in receipts))+' bytes; see per-job caps/receipts.')
        marker=early/'pipeline/done.exit'
        lines.append('T25 early pipeline marker: '+str(marker)+'; '+(marker.read_text().strip() if marker.exists() else 'PENDING (not a numerical verdict)'))
        marker=early/'postmap-job/done.exit'
        if marker.exists():lines.append('T25 early postmap marker: '+str(marker)+'; '+marker.read_text().strip())
        final=HERE/'t25-early-final-audit.json'
        if final.exists():
            f=json.loads(final.read_text())
            assert f['status']=='PASS' and f['postmap_native_receipts_verified']==20 and f['endpoint_N96_sign_agreements']==10
            lines.append('T25 early final: 18 hole/time rows; 10 native bounded searches; 0 strict horizon passes. Maps cap 3 GB, finders cap 6 GB, serial 1 thread.')
    for p in files:lines.append(f'PACKET {sha(p)} {p.stat().st_size} {p}')
    # All raw numerical and resource evidence, including failed attempts.
    suffixes={'.json','.csv','.txt','.log','.dat','.time'}
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or p.suffix not in suffixes and p.name!='done.exit':continue
        if 'build' in p.relative_to(ROOT).parts or 'portable-build' in p.relative_to(ROOT).parts:continue
        lines.append(f'EXTERNAL {sha(p)} {p.stat().st_size} {p}')
    for item in audit['inputs']:
        lines.append(f'INPUT {item["sha256"]} {item["bytes"]} {item["path"]}')
    if (early/'inputs.json').exists():
        for item in json.loads((early/'inputs.json').read_text()):
            lines.append(f'INPUT_EARLY {item["sha256"]} {item["bytes"]} {item["path"]}')
    if (early/'initialization-d16/inputs.json').exists():
        for item in json.loads((early/'initialization-d16/inputs.json').read_text())['inputs']:
            lines.append(f'INPUT_T0_ONLY {item["sha256"]} {item["bytes"]} {item["path"]}')
    for p in [ROOT/'build/t25-expansion-map.ex',Path(reg['portable_executable'])]:
        lines.append(f'EXTERNAL_EXE {sha(p)} {p.stat().st_size} {p}')
    for p in [ROOT/'initialization'/k/'chk/EMS_000000.2d.hdf5' for k in ['single','binary']]:
        assert sha(p)==next(i['sha256'] for i in audit['inputs'] if i['path']==str(p))
    (HERE/'t25-manifest.txt').write_text('\n'.join(lines)+'\n')
    for p in files:
        if p.suffix=='.csv':assert p.stat().st_size<1000000,p
    print('T25_PACKET_SEALED',len(files),'packet files;',len(lines),'manifest entries')
if __name__=='__main__':main()
