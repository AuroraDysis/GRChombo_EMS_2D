#!/usr/bin/env python3
"""Strict local build against the clean frozen d507438 dependency archive.

No dirty Source/Examples files (including SafeNanAbort) enter this build.
The portable cluster build is t25.make. Local dependency/build provenance is
written to t25-build.json; executable and objects remain under /private/tmp.
"""
import hashlib,json,os,resource,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/ems-t25/build')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    spec=json.loads(Path('/private/tmp/ems-t24/build/build-spec.json').read_text())
    objects=[str(Path('/private/tmp/ems-t24/build')/(s+'.o')) for s in
        ('GRAMR','GRAMRLevel','GRLevelData','BoundaryConditions','SmallDataIO','PETScCommunicator','PunctureTracker')]
    obj=ROOT/'map.o';exe=ROOT/'t25-expansion-map.ex';temporary=ROOT/'t25-expansion-map.next.ex'
    commands=[spec['flags']+['-fno-access-control','-c',str(HERE/'t25-expansion-map.cpp'),'-o',str(obj)],
        spec['flags']+[str(obj),*objects,*spec['libraries'],'-o',str(temporary)]]
    for command in commands:
        print('BUILD',command,flush=True);subprocess.run(command,check=True)
    temporary.replace(exe)
    files=[HERE/'t25-expansion-map.cpp',*map(Path,objects),exe]
    source=Path('/private/tmp/ems-t24/source')
    files += [source/'Source/RHFinder'/s for s in ('RHSurf.hpp','RHUnion.hpp')]
    files += list(Path('/private/tmp/ems-t24/Chombo/lib').glob('*.a'))
    result=dict(source_commit='d50743867e48bcb466a60ca16869b07fe1546c27',
        clean_source=str(source),commands=commands,compiler=spec['compiler'],
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        serial=True,max_threads=1,files={str(p):sha(p) for p in files},
        no_dirty_source_or_nan_abort=True)
    (HERE/'t25-build.json').write_text(json.dumps(result,indent=2)+'\n')
    print('T25_BUILD_COMPLETE',result['peak_rss_bytes'],flush=True)
if __name__=='__main__': main()
