"""Authenticate supplied bytes and make identical local readback snapshots."""
from pathlib import Path
import hashlib,json,time
root=Path('/Users/auroradysis/Workspace/EMS')
src=root/'.ariadne/worktrees/wt-pulse'
manifest={line.split(maxsplit=1)[1]:line.split()[0] for line in (src/'artifacts/echo-pulse/data.sha256').read_text().splitlines()}
snapshot=Path('Tests/EMSCTT2/results/input-snapshots');snapshot.mkdir(exist_ok=True)
paths=[root/'artifacts/echo/a0.9-e7.trumpet',src/'docs/binary-ctt-format.md',src/'src/ctt_pulse_io.jl',src/'src/ctt_pulse.jl',src/'src/binary_ctt.jl',src/'src/binary_ctt_end.jl']
paths += [Path('/Users/auroradysis/SciData/EMS/.data/echo-pulse')/f'{n}-final.ctt' for n in ['zero','half','n12','n18','n24']]
def read(path):
    # The SciData mount intermittently reports ENOENT/EMFILE; retry transport,
    # and still require the supplied digest before using any returned bytes.
    for attempt in range(3):
        try:return path.read_bytes()
        except OSError:
            if attempt==2:raise
            time.sleep(attempt+1)
rows=[]
for p in paths:
    b=read(p);h=hashlib.sha256(b).hexdigest()
    expected=manifest.get(str(p),manifest.get(str(root/'.data/echo-pulse'/p.name)))
    row=dict(path=str(p),bytes=len(b),sha256=h,manifest_sha256=expected)
    if p.suffix=='.ctt':
        payload,trailer=b.rsplit(b'SHA256 ',1)
        assert hashlib.sha256(payload).hexdigest()==trailer.decode().strip()
        assert expected==h
        meta=dict(line[2:].split('=',1) for line in payload.decode().splitlines() if line.startswith('# '))
        provenance={}
        for key,source in [('profile_sha256',root/'artifacts/echo/a0.9-e7.trumpet'),('solver_sha256',src/'src/binary_ctt.jl'),('end_solver_sha256',src/'src/binary_ctt_end.jl'),('pulse_sha256',src/'src/ctt_pulse.jl'),('checkpoint_sha256',p.with_name(p.name.replace('-final.ctt','.bin')))]:
            try:actual=hashlib.sha256(read(source)).hexdigest()
            except OSError as e:
                provenance[key]=dict(path=str(source),expected=meta[key],read_error=str(e));continue
            assert actual==meta[key],(p,key,source)
            provenance[key]=dict(path=str(source),expected=meta[key],actual=actual)
        row['provenance']=provenance
    if p.suffix=='.trumpet':assert h=='74f85f39aa652b48d179a96b8dbe57cd4468d7905e39732aaf8bd1b75931b8fa'
    if p.suffix in ('.ctt','.trumpet'):
        dest=snapshot/p.name;dest.write_bytes(b)
        assert hashlib.sha256(dest.read_bytes()).hexdigest()==h
        row['local_snapshot']=str(dest.resolve())
    rows.append(row)
Path('Tests/EMSCTT2/results/authentication.json').write_text(json.dumps(rows,indent=2)+'\n')
print('authenticated profile, five payloads and manifest hashes; current spec/source hashes and provenance checks recorded')
