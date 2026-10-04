#!/usr/bin/env python3
"""Local MPI abort adapter regression; distinguishes Init failure from tested abort."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t24/mpi-abort')
MPI=Path('/Users/auroradysis/.julia/artifacts/2e0bbdf5bae18755b0ebba03bc02d6fb568c05bc')

def main():
    if len(sys.argv)>1 and sys.argv[1]=='--launch':
        # macOS removes DYLD_* at /usr/bin/time; restore the dependency search
        # path inside its child, before launching Hydra (not in the time parent).
        os.environ.update(json.loads(Path('/private/tmp/ems-t24/mpi-probe.json').read_text())['env'])
        raise SystemExit(subprocess.run(sys.argv[2:],env=os.environ).returncode)
    s=importlib.util.spec_from_file_location('controls',HERE/'t24-controls.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    ROOT.mkdir(parents=True,exist_ok=True)
    spec=json.loads(Path('/private/tmp/ems-t24/controls/build-spec.json').read_text())
    flags=spec['flags']+['-DCH_MPI','-I'+str(MPI/'include')]
    libs=spec['libraries']+['-L'+str(MPI/'lib'),'-Wl,-rpath,'+str(MPI/'lib'),'-lmpi']
    # Use the dependency's own MPI SPMD unit, ahead of its serial archive; no
    # redefinition of procID and no mixed serial/MPI communicator metadata.
    m.measured('build-mpi-spmd',flags+['-c','/private/tmp/ems-t24/Chombo/lib/src/BaseTools/SPMD.cpp','-o',str(ROOT/'SPMD-MPI.o')],ROOT/'spmd')
    m.measured('build-mpi-probe',flags+[str(HERE/'t24-safe-nan-probe.cpp'),str(ROOT/'SPMD-MPI.o'),*libs,'-o',str(ROOT/'probe.ex')],ROOT/'build')
    env=json.loads(Path('/private/tmp/ems-t24/mpi-probe.json').read_text())['env'];os.environ.update(env)
    rows=[]
    for label,n,mode,bad,expect in [('one-finite',1,'finite',0,0),('one-abort',1,'nan',0,86),('two-finite',2,'finite',1,0),('two-abort',2,'nan',1,86)]:
        d=ROOT/label;d.mkdir(exist_ok=True)
        cmd=[MPI/'bin/mpiexec.hydra','-launcher','fork','-n',str(n),ROOT/'probe.ex',mode,d/'abort',str(bad)]
        with (d/'run.log').open('wb') as log:
            p=subprocess.run([sys.executable,str(HERE/'t13-run.py'),'--measure',label,'--directory',str(d),'--',sys.executable,str(Path(__file__).resolve()),'--launch',*map(str,cmd)],stdout=log,stderr=subprocess.STDOUT,timeout=30)
        (d/'done.exit').write_text(str(p.returncode)+'\n')
        q=json.loads((d/(label+'.resources.json')).read_text());text=(d/'run.log').read_text()
        record=d/f'abort.rank{bad}.json';bits=list(d.glob('abort.bits.rank*.json'))
        row=dict(case=label,ranks=n,native_returncode=q['returncode'],expected_returncode=expect,
            MPI_Init_failed='internal_Init' in text,abort_record_closed=record.exists(),
            field_bits_identical=all(json.loads(f.read_text())['bit_identical'] for f in bits) if bits else None,
            peak_RSS_bytes=q['peak_rss_bytes'],wall_s=q['wall_seconds'],receipt=str(d/(label+'.resources.json')))
        row['passed']=(q['returncode']!=0 if mode=='nan' else q['returncode']==0) and not q['gate_reason'] and len(bits)==n and row['field_bits_identical'] and (record.exists() if mode=='nan' else not record.exists())
        rows.append(row)
    result=dict(status='MPI_ADAPTER_PASS' if all(r['passed'] for r in rows) else 'MPI_QUALIFICATION_BLOCKED',cases=rows,
        limitation='Tests-only serial-Chombo FAB adapter with actual MPI calls; no multi-rank AMR qualification',
        MPI= str(MPI))
    (HERE/'t24-mpi-abort-qualification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
