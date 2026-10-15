"""Run bounded local controls sequentially; record wall time and child peak RSS."""
from pathlib import Path
import os,subprocess,time,json,resource
R=Path(__file__).resolve().parent
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',JULIA_NUM_THREADS='1',CHOMBO_HOME='/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
records=[]
def run(name,args,cap=240):
    start=time.monotonic()
    with (R/'results'/f'{name}.log').open('w') as log:
        p=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,timeout=cap)
    rss=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    record=dict(name=name,command=args,seconds=time.monotonic()-start,exit=p.returncode,child_peak_rss_upper_bound_bytes=rss)
    records.append(record)
    (R/'results'/'cost.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    assert p.returncode==0,name
    assert rss<=6_000_000_000,(name,rss)
if __name__=='__main__':
    import sys
    phase=sys.argv[1]
    if phase in ('hierarchy','near'):
        exe=str(next(R.glob('EMSCTT2Hierarchy*.ex')))
        for n in (80,160,320):
            stem='near-' if phase=='near' else ''
            run('hierarchy-'+stem+str(n),[exe,str(R/f'params-{stem}{n}.txt')])
    elif phase=='regression':
        for path in ('Tests/EMSCTT','Tests/EMSTrumpet','Examples/EMS'):
            run('build-'+Path(path).name,['make','-C',path,'all','DIM=2','-j1'])
        root=Path('/Users/auroradysis/Workspace/EMS')
        run('emsctt1',[str(next(Path('Tests/EMSCTT').glob('EMSCTTFixtureTest*.ex'))),str(root/'artifacts/reference-alpha20-qfile.trumpet'),str(root/'.data/binary-ctt-t6/n28-r6.ctt'),str(root/'test/fixtures/emsctt1'),str(R/'results')])
        run('trumpet',[str(next(Path('Tests/EMSTrumpet').glob('EMSTrumpetFixtureTest*.ex'))),'Tests/EMSTrumpet/fixtures','Tests/EMSTrumpet/fixtures-echo/B','Tests/EMSTrumpet/fixtures-echo/E'])
    elif phase=='reference-ccz4':
        run('reference-ccz4',['julia','--heap-size-hint=2G','--project=/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-pulse',str(R/'reference_ccz4.jl'),str(R/'results/input-snapshots/a0.9-e7.trumpet'),str(R/'results/input-snapshots'),str(R/'results')],cap=290)
    elif phase=='reference-ad':
        run('reference-ad',['julia','--heap-size-hint=2G','--project=/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-pulse',str(R/'reference_ad.jl'),str(R/'results/input-snapshots/a0.9-e7.trumpet'),str(R/'results/input-snapshots/n24-final.ctt'),str(R/'results/jet-points.txt'),str(R/'results/jets-replay.txt')],cap=290)
    elif phase=='ghost-control':
        exe=str(next(R.glob('EMSCTT2Hierarchy*.ex')))
        run('analytic-ghosts',[exe,str(R/'params-320.txt'),'ctt2_test_analytic_ghosts=1'])
    else:raise ValueError(phase)
    (R/'results'/f'cost-{phase}.json').write_text(json.dumps(records,indent=2)+'\n')
