#!/usr/bin/env python3
"""Receipt-first source/build bisect; classify every saved valid/ghost difference."""
import csv,hashlib,json,os,re,resource,subprocess,sys
from pathlib import Path
import h5py,numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t23');REPO=HERE.parents[1];PY=sys.executable
FINDING=Path('/var/folders/t9/yv5hdkxs0tdf_n6b11fwrwgr0000gn/T/ariadne-explore-t8T3iZ/stage-exp-0025/state/threads/ems-spectral-solver/runs/exp-0025/submissions/exp-0025/blocked-finding.json')
CLUSTER=Path('/Users/auroradysis/Workspace/EMS/.data/exp-0025/submission/attempt-4-failed-evidence')
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(name,data):
    if not data:return
    with (HERE/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
def own(d,name):
    d=Path(d);q=json.loads((d/(name+'.resources.json')).read_text())
    assert (d/'done.exit').read_text().strip()=='0' and q['returncode']==q['child_measurement']['returncode']==0 and not q['gate_reason'] and q['peak_rss_bytes']<6e9,(d,q)
    return dict(case=name,receipt=str(d/(name+'.resources.json')),returncode=0,peak_RSS_bytes=q['peak_rss_bytes'],wall_seconds=q['wall_seconds'])
def names(f):
    def string(a):return a.decode() if isinstance(a,bytes) else str(a)
    return [string(f.attrs['component_'+str(c)]) for c in range(int(f.attrs['num_components']))]
def comparisons(a,b,case,ledger=None):
    result=[]
    params=(Path(a).parents[1]/'params.txt').read_text()
    centre=float(re.search(r'^center\s*=\s*(\S+)',params,re.M)[1])
    with h5py.File(a) as f,h5py.File(b) as g:
        assert int(f.attrs['num_levels'])==int(g.attrs['num_levels']);cs=names(f);assert cs==names(g)
        for l in range(int(f.attrs['num_levels'])):
            x,y=f[f'level_{l}'],g[f'level_{l}'];assert x.attrs['dx']==y.attrs['dx'] and x.attrs['time']==y.attrs['time']
            boxes=x['boxes'][:];assert np.array_equal(boxes,y['boxes'][:]);offset=x['data:offsets=0'][:];assert np.array_equal(offset,y['data:offsets=0'][:])
            gx,gy=map(int,x['data_attributes'].attrs['outputGhost']);assert np.array_equal(x['data_attributes'].attrs['outputGhost'],y['data_attributes'].attrs['outputGhost'])
            h=float(x.attrs['dx']);bounds=[tuple(int(z[k]) for k in ('lo_i','lo_j','hi_i','hi_j')) for z in boxes]
            for bi,box in enumerate(bounds):
                x0,y0,x1,y1=box;nx=x1-x0+1;ny=y1-y0+1;shape=(len(cs),ny+2*gy,nx+2*gx)
                aa=x['data:datatype=0'][int(offset[bi]):int(offset[bi+1])].reshape(shape);bb=y['data:datatype=0'][int(offset[bi]):int(offset[bi+1])].reshape(shape)
                xx,yy=np.meshgrid(np.arange(x0-gx,x1+gx+1),np.arange(y0-gy,y1+gy+1));valid=(xx>=x0)&(xx<=x1)&(yy>=y0)&(yy<=y1)
                covered=np.zeros(xx.shape,dtype=bool)
                for q0,r0,q1,r1 in bounds:covered|=(xx>=q0)&(xx<=q1)&(yy>=r0)&(yy<=r1)
                domain=x.attrs['prob_domain'];d0,d1,d2,d3=[int(domain[k]) for k in ('lo_i','lo_j','hi_i','hi_j')]
                physical=(xx<d0)|(xx>d2)|(yy<d1)|(yy>d3)
                regions={'valid':valid,'same_level_ghost':~valid&covered,'physical_boundary_ghost':~valid&~covered&physical,'uncovered_refinement_ghost':~valid&~covered&~physical}
                for c,field in enumerate(cs):
                    bits=aa[c].copy().view('u8')!=bb[c].copy().view('u8')
                    for region,take in regions.items():
                        changed=bits&take;count=int(changed.sum());values=int(take.sum())
                        if not values:continue
                        finite=np.isfinite(aa[c])&np.isfinite(bb[c]);diff=np.abs(aa[c]-bb[c]);maximum=float(diff[changed&finite].max()) if np.any(changed&finite) else 0.
                        result.append(dict(case=case,file=Path(a).name,level=l,box=bi,field=field,region=region,values=values,bit_mismatches=count,max_abs_difference=maximum,nonfinite_reference=int((~np.isfinite(aa[c])&take).sum()),nonfinite_candidate=int((~np.isfinite(bb[c])&take).sum())))
                        if ledger is not None and count:
                            iy,ix=np.where(changed)
                            for v,u in zip(iy,ix):ledger.writerow(dict(case=case,file=Path(a).name,level=l,box=bi,field=field,region=region,i=int(xx[v,u]),j=int(yy[v,u]),x_M=(xx[v,u]+.5)*h-centre,y_M=(yy[v,u]+.5)*h,axis_ghost=bool(yy[v,u]<0),reference=aa[c,v,u],candidate=bb[c,v,u]))
    return result
def fallback():
    a=ROOT/'runs/old-production/plt/EMS_Plot_000000.2d.hdf5';b=ROOT/'runs/candidate-production/plt/EMS_Plot_000000.2d.hdf5'
    own(a.parents[1],'old-production');own(b.parents[1],'candidate-production')
    rows=comparisons(a,b,'small-prescreen');n=sum(r['bit_mismatches'] for r in rows if r['field'] in ('Ham','Mom1','Mom2','Mom','GaussE','GaussB'))
    if n:
        (ROOT/'high-fallback/decision.json').write_text(json.dumps(dict(status='SMALL_REPRODUCTION_FOUND_NO_HIGH_RUN_NEEDED',diagnostic_differences=n),indent=2)+'\n');print('T23_SMALL_REPRODUCTION_FOUND',n,flush=True);return
    # Only initialization and plot0 on the real hierarchy; no fine evolution.
    from importlib.util import spec_from_file_location,module_from_spec
    s=spec_from_file_location('prepare',HERE/'t23-prepare.py');m=module_from_spec(s);s.loader.exec_module(m)
    seed=Path('/private/tmp/ems-t22/equilibrium/experimental/params.txt').read_text()
    seed=re.sub(r'^ems_gauge\s*=.*\n?','',seed,flags=re.M)
    for name in ('old','candidate'):
        d=ROOT/'runs'/('high-'+name)
        for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
        p=d/'params.txt';p.write_text(m.replace(seed,dict(max_steps=0,checkpoint_interval=1,plot_interval=1,RH_activate='false',t21_rhs_capture='false')))
        assert 'ems_gauge' not in p.read_text()
        with (d/'run.log').open('wb') as log:rc=subprocess.run([PY,str(HERE/'t13-run.py'),'--measure','high-'+name,'--directory',str(d),'--',str(ROOT/'builds'/name/'dense.ex'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
        (d/'done.exit').write_text(str(rc)+'\n');assert rc==0,(name,rc)
    (ROOT/'high-fallback/decision.json').write_text(json.dumps(dict(status='HIGH_INITIAL_ONLY_RUNS_COMPLETED',small_diagnostic_differences=0),indent=2)+'\n')
def cluster():
    finding=json.loads(FINDING.read_text());receipt=json.loads((CLUSTER/'smoke/identity/identity-mismatch.receipt.json').read_text());assert finding['total_bit_mismatches']==sum(r['values'] for r in receipt['fields'].values())==57806
    rows=[]
    for case in ('baseline','omitted','explicit','capture'):
        d=CLUSTER/'smoke/identity'/case;q=json.loads((d/'native.receipt.json').read_text());assert q['status']=='NATIVE_SUCCESS' and q['OMP_NUM_THREADS']==4 and q['native_exits']==[0]*32
        exits=[int(p.read_text()) for p in d.glob('native-exit-rank-*.txt')];assert sorted(exits)==[0]*32
        rows.append(dict(case=case,native_rank_count=32,native_rank_exit_max=0,OpenMP_threads=4,scope='completed cluster native cases; identity assertion failed before one-step comparison',receipt_sha256=sha(d/'native.receipt.json')))
    save('t23-cluster-native-receipts.csv',rows);save('t23-cluster-differences.csv',[dict(field=k,**v) for k,v in receipt['fields'].items()])
    return finding
def analyze():
    finding=cluster();plan=json.loads((ROOT/'pipeline/plan.json').read_text());receipts=[]
    for job in plan['jobs']:
        if job['name']=='qualification':continue
        receipts.append(own(job['directory'],job['name']))
    runs=ROOT/'runs';ledger=ROOT/'ledgers';ledger.mkdir(exist_ok=True);allrows=[]
    with (ledger/'cell-mismatches.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames='case file level box field region i j x_M y_M axis_ghost reference candidate'.split());w.writeheader()
        cases=[]
        for name in ('setter','tracking','recording','candidate','candidate-native-flags','old-all-access','old-extra'):
            for factory in ('production','dense'):cases.append(('old-'+factory,name+'-'+factory))
        cases += [('old-production','old-dense'),('candidate-production','candidate-dense'),('candidate-production','candidate-native-flags-production'),('candidate-dense','candidate-native-flags-dense'),('poison-a','poison-b'),('old-production','poison-a')]
        for left,right in cases:
            for prefix in ('chk/EMS_','plt/EMS_Plot_'):
                for step in (0,1):
                    file=f'{prefix}{step:06}.2d.hdf5';a=runs/left/file;b=runs/right/file
                    assert a.exists() and b.exists(),(a,b)
                    allrows.extend(comparisons(a,b,left+' vs '+right,w))
        if (runs/'high-old/done.exit').exists():
            for name in ('old','candidate'):receipts.append(own(runs/('high-'+name),'high-'+name))
            allrows.extend(comparisons(runs/'high-old/plt/EMS_Plot_000000.2d.hdf5',runs/'high-candidate/plt/EMS_Plot_000000.2d.hdf5','high-old vs high-candidate',w))
    from collections import defaultdict
    summary=defaultdict(lambda:dict(values=0,bit_mismatches=0,max_abs_difference=0.,nonfinite_reference=0,nonfinite_candidate=0))
    for r in allrows:
        key=tuple(r[k] for k in ('case','file','field','region'));a=summary[key]
        for k in ('values','bit_mismatches','nonfinite_reference','nonfinite_candidate'):a[k]+=r[k]
        a['max_abs_difference']=max(a['max_abs_difference'],r['max_abs_difference'])
    compact=[dict(zip(('case','file','field','region'),key),**value) for key,value in summary.items()]
    save('t23-comparisons.csv',compact);save('t23-location-summary.csv',[r for r in allrows if r['bit_mismatches']]);save('t23-local-resources.csv',receipts)
    # No undocumented tolerance replaces a bitwise comparison. These are
    # measured counts, including any compiler-induced valid-cell differences.
    evolved={'chi','h11','h12','h22','hww','K','A11','A12','A22','Aww','Theta','Gamma1','Gamma2','lapse','shift1','shift2','B1','B2','phi','Pi','Lambda','Bx','By','Bz','Ex','Ey','Ez','Xi'}
    physics=sum(r['bit_mismatches'] for r in compact if r['field'] in evolved)
    diagnostic_valid=sum(r['bit_mismatches'] for r in compact if r['field'] not in evolved and r['region']=='valid')
    poison=[r for r in compact if r['case']=='poison-a vs poison-b' and r['bit_mismatches']]
    poison_outside=all(r['field'] not in evolved and r['region']!='valid' for r in poison)
    # The witness must be the two deliberately supplied values themselves,
    # not merely some difference caused by running an instrumented binary.
    witness={}
    with (ledger/'cell-mismatches.csv').open() as f:
        for r in csv.DictReader(f):
            if r['case']!='poison-a vs poison-b':continue
            key=tuple(r[k] for k in ('file','level','field','region'))
            w=witness.setdefault(key,dict(changed_values=0,exact_seed_pair=0))
            w['changed_values']+=1
            w['exact_seed_pair']+=float(r['reference'])==1.25 and float(r['candidate'])==2.5
    witness_rows=[dict(zip(('file','level','field','region'),k),**v) for k,v in witness.items()]
    save('t23-poison-witness.csv',witness_rows)
    poison_exact=bool(witness_rows) and all(w['changed_values']==w['exact_seed_pair'] for w in witness_rows)
    untouched_halos=[r for r in compact if r['case']=='poison-a vs poison-b' and r['field'] not in evolved and r['region'] in ('physical_boundary_ghost','uncovered_refinement_ghost')]
    poison_complete=bool(untouched_halos) and all(r['bit_mismatches']==r['values'] for r in untouched_halos)
    source=json.loads((HERE/'t23-source-inputs.json').read_text());same_writer=len({source[n]['source_sha256']['Source/GRChomboCore/GRAMRLevel.cpp'] for n in ('old','setter','tracking','recording','candidate')})==1
    build_rows=[]
    for left,right in [('old','old-all-access'),('candidate','candidate-native-flags')]:
        for p in sorted((ROOT/'builds'/left).iterdir()):
            if p.suffix not in ('.o','.ex'):continue
            other=ROOT/'builds'/right/p.name;assert other.exists()
            build_rows.append(dict(reference=left,candidate=right,artifact=p.name,reference_sha256=sha(p),candidate_sha256=sha(other),bit_identical=sha(p)==sha(other)))
    save('t23-build-identity.csv',build_rows)
    repeat_rows=[r for r in compact if r['case']=='old-production vs old-all-access-production' and r['file'].startswith('EMS_Plot_')]
    repeat_counts={file:sum(r['bit_mismatches'] for r in repeat_rows if r['file']==file) for file in {r['file'] for r in repeat_rows}}
    result=dict(status='QUALIFIED_UNWRITTEN_PLOT_DIAGNOSTIC_GHOSTS' if physics==diagnostic_valid==0 and poison and poison_outside and poison_exact and poison_complete and same_writer else 'AS_MEASURED_FURTHER_DIAGNOSIS_REQUIRED',evolved_bit_mismatches=physics,valid_diagnostic_bit_mismatches=diagnostic_valid,evolved_Float64_comparisons=sum(r['values'] for r in compact if r['field'] in evolved),valid_diagnostic_Float64_comparisons=sum(r['values'] for r in compact if r['field'] not in evolved and r['region']=='valid'),poison_differences=sum(r['bit_mismatches'] for r in poison),poison_only_diagnostic_ghosts=poison_outside,poison_exact_seed_pair=poison_exact,poison_fills_all_uncovered_diagnostic_halos=poison_complete,writer_identical_across_source_commits=same_writer,access_flag_artifacts_identical=all(r['bit_identical'] for r in build_rows),identical_old_executable_repeat_ghost_mismatches=repeat_counts,peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in receipts),full_cell_ledger=str(ledger/'cell-mismatches.csv'),source_commits=source,cluster_location_unavailable='cluster receipt is per box only; raw HDF5 files were not collected',production_source_changed=False,identity_check_changed=False)
    (HERE/'t23-qualification.json').write_text(json.dumps(result,indent=2)+'\n');report(result,compact)
    seal();print(result['status'],flush=True)
def split_prepare():
    from importlib.util import spec_from_file_location,module_from_spec
    s=spec_from_file_location('prepare',HERE/'t23-prepare.py');m=module_from_spec(s);s.loader.exec_module(m)
    root=ROOT/'split-probe';root.mkdir(exist_ok=True);jobs=[]
    seed=(ROOT/'runs/poison-a/params.txt').read_text()
    for label,value in [('split-poison-a','1.25'),('split-poison-b','2.5')]:
        d=ROOT/'runs'/label
        for sub in ('chk','plt'):(d/sub).mkdir(parents=True,exist_ok=True)
        p=d/'params.txt';p.write_text(m.replace(seed,dict(max_box_size=32)))
        jobs.append(dict(name=label,directory=str(d),command=['/usr/bin/env','T23_PLOT_BUFFER_SEED='+value,str(ROOT/'builds/poison/production.ex'),str(p)]))
    (root/'plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    print('T23_SPLIT_PROBE_READY',flush=True)
def split_analyze():
    a=ROOT/'runs/split-poison-a';b=ROOT/'runs/split-poison-b';rows=[]
    with (ROOT/'ledgers/split-cell-mismatches.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames='case file level box field region i j x_M y_M axis_ghost reference candidate'.split());w.writeheader()
        for prefix in ('chk/EMS_','plt/EMS_Plot_'):
            for step in (0,1):
                file=f'{prefix}{step:06}.2d.hdf5';rows.extend(comparisons(a/file,b/file,'split-poison-a vs split-poison-b',w))
    save('t23-split-comparisons.csv',rows)
    with h5py.File(a/'chk/EMS_000000.2d.hdf5') as f:evolved=set(names(f))
    defined=[r for r in rows if r['field'] in evolved or r['region'] in ('valid','same_level_ghost')]
    assert all(r['bit_mismatches']==0 for r in defined)
    internal=[r for r in rows if r['field'] in ('Ham','Mom1','Mom2','Mom','GaussE','GaussB') and r['region']=='same_level_ghost']
    assert internal and sum(r['values'] for r in internal)>0
    poisoned=[r for r in rows if r['bit_mismatches']]
    assert poisoned and all(r['bit_mismatches']==r['values'] and r['max_abs_difference']==1.25 for r in poisoned)
    with (ROOT/'ledgers/split-cell-mismatches.csv').open() as f:
        for r in csv.DictReader(f):assert float(r['reference'])==1.25 and float(r['candidate'])==2.5
    return dict(status='EXACT_POISON_ONLY_UNCOVERED_HALOS',poisoned_values=sum(r['bit_mismatches'] for r in poisoned),same_level_diagnostic_ghost_values=sum(r['values'] for r in internal),same_level_diagnostic_ghost_mismatches=0,defined_value_comparisons=sum(r['values'] for r in defined),defined_value_mismatches=0)
def finalize():
    plan=json.loads((ROOT/'pipeline/plan.json').read_text())
    receipts=[own(j['directory'],j['name']) for j in plan['jobs']]
    queue=own(ROOT/'continuation','queue');queue['case']='continuation-queue';receipts.append(queue)
    split_plan=json.loads((ROOT/'split-probe/plan.json').read_text())
    receipts.extend(own(j['directory'],j['name']) for j in split_plan['jobs'])
    queue=own(ROOT/'split-probe','queue');queue['case']='split-probe-queue';receipts.append(queue)
    split=split_analyze();(HERE/'t23-split-qualification.json').write_text(json.dumps(split,indent=2)+'\n')
    failed=ROOT/'builds/candidate-native-flags.failed-private-base';f=json.loads((failed/'build-candidate-native-flags.resources.json').read_text())
    assert (failed/'done.exit').read_text().strip()=='1' and f['returncode']==f['child_measurement']['returncode']==1 and not f['gate_reason']
    save('t23-failed-build.csv',[dict(case='build-candidate-native-flags-original',returncode=1,peak_RSS_bytes=f['peak_rss_bytes'],wall_seconds=f['wall_seconds'],receipt=str(failed/'build-candidate-native-flags.resources.json'),cause='EMSBH2DLevel.cpp:712,729 calls private GRAMRLevel header methods',repair='access flag on EMSBH2DLevel and dense main only; Tests build only')])
    validations=[]
    for j in plan['jobs']+split_plan['jobs']:
        if Path(j['directory']).parent!=ROOT/'runs':continue
        d=Path(j['directory']);outputs=[]
        for prefix in ('chk/EMS_','plt/EMS_Plot_'):
            for step in (0,1):
                p=d/f'{prefix}{step:06}.2d.hdf5'
                with h5py.File(p) as out:
                    assert int(out.attrs['num_levels'])==2
                    for level in range(2):assert float(out[f'level_{level}'].attrs['time'])==step*.03125
                    assert len(names(out))==(28 if prefix.startswith('chk') else 34)
                outputs.append(sha(p))
        validations.append(dict(case=j['name'],own_returncode=0,levels=2,coarse_steps=1,t0_M=0.,t1_M=.03125,output_files=4,guard='t2_guard_initial_data_after_t0=true',output_sha256=';'.join(outputs)))
    save('t23-run-validation.csv',validations);save('t23-local-resources.csv',receipts)
    q=json.loads((HERE/'t23-qualification.json').read_text());q.update(completed_jobs=len(plan['jobs'])+len(split_plan['jobs']),completed_native_runs=len(validations),continuation_marker=str(ROOT/'continuation/done.exit'),split_probe=split,all_own_receipts_verified=True,peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in receipts))
    (HERE/'t23-qualification.json').write_text(json.dumps(q,indent=2)+'\n')
    with (HERE/'t23-comparisons.csv').open() as f:rows=list(csv.DictReader(f))
    for r in rows:r['bit_mismatches']=int(r['bit_mismatches'])
    report(q,rows);seal();print('T23_PACKET_SEALED',q['status'],flush=True)
def report(result=None,rows=()):
    q=result or dict(status='RUNS_PENDING',evolved_bit_mismatches='pending',valid_diagnostic_bit_mismatches='pending',poison_differences='pending',peak_RSS_bytes='pending')
    table=['| comparison | file | evolved differences | valid diagnostic differences | diagnostic ghost differences |','| --- | --- | ---: | ---: | ---: |']
    keys=sorted({(r['case'],r['file']) for r in rows});evolved=set('chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split())
    for case,file in keys:
        a=[r for r in rows if r['case']==case and r['file']==file];e=sum(r['bit_mismatches'] for r in a if r['field'] in evolved);v=sum(r['bit_mismatches'] for r in a if r['field'] not in evolved and r['region']=='valid');g=sum(r['bit_mismatches'] for r in a if r['field'] not in evolved and r['region']!='valid');table.append(f'| {case} | {file} | {e} | {v} | {g} |')
    mechanism='The controlled local receipts qualify this as an existing unwritten-output defect; they do not identify a deterministic first bad evolution commit.' if q['status']=='QUALIFIED_UNWRITTEN_PLOT_DIAGNOSTIC_GHOSTS' else 'This is a candidate explanation until the local receipts and poison witnesses qualify it.'
    text=f'''# T23 default-path diagnostic difference

Status: **{q['status']}**. Exp-0025 admission remains blocked; no cluster retry, production fix or identity-check edit is performed. The cluster finding is verified against its own four native receipts (32 zero rank exits each) and failed identity receipt: 57,806 differences, six diagnostic fields, no evolved-field differences at plot0. Cluster HDF5 files were not collected, so the original per-cell locations cannot be reconstructed from the per-box receipt alone. Local controlled evidence supplies the cell classification, not invented cluster coordinates.

The source inspection finds a pre-existing hole in plot output. `Source/GRChomboCore/GRAMRLevel.cpp:907` allocates `plot_data` with three ghost cells and no initialization. With `write_plot_ghosts=true`, evolved fields use the boundary copier at lines 927–928, while diagnostic copyTo at lines 936–937 omits it. The final same-level exchange at line 967 cannot populate uncovered refinement ghosts or physical-boundary diagnostic ghosts. `Examples/EMS/EMSBH2DLevel.cpp:302–325` computes native constraints on valid cells only. This plot writer is unchanged from 5da576b through 758f0d2. A difference in heap reuse can therefore expose saved diagnostic values with no defined numerical meaning, even when evolved fields agree exactly. {mechanism}

The earlier T22 plot regression requested only 28 evolved components; it did not establish six diagnostic-ghost identities. The T23 matrix explicitly plots all 34 components, with plot ghosts enabled, in both the production main and dense-tag factory. It isolates 31176be, d507438, a pre-selector source with only the qualified T21 hook, and 758f0d2, plus every-unit access flags and the actual unused MovingGauge object linked into the old executable. Every native unit is rebuilt against its own source/parameter layout; cross-commit GRAMRLevel.o reuse would be invalid. Local GCC is 16 on macOS, whereas the cluster compiler is GCC 11.5 on Linux: local results test the source mechanism and do not certify that compiler/MPI combination.

The original matrix stopped at `build-candidate-native-flags`: own exit 1, no resource gate, because the candidate's checkpoint overrides call private base methods at `Examples/EMS/EMSBH2DLevel.cpp:712,729` (`GRAMRLevel.hpp:98,102`). The failed build is preserved under `/private/tmp/ems-t23/builds/candidate-native-flags.failed-private-base/`, with its own [failure record](t23-failed-build.csv). The Tests-only repair adds `-fno-access-control` to EMSBH2DLevel as well as the dense main; every other native unit retains normal access checks. Thus this control is explicitly **main plus the necessary level unit**, not a claim that the pinned candidate compiles with the old main-only scope. Production sources are not repaired. The five successful original builds were reused; all remaining cases completed under `/private/tmp/ems-t23/continuation/done.exit = 0`, each verified independently.

There is no deterministic first bad numerical commit in this matrix. The earliest source-stage contrast, 31176be versus 5da576b, already shows 3,144 undefined diagnostic-ghost differences at t=0; nevertheless the old/all-access executables and all their native objects have identical SHA-256 hashes. Running that **same executable content** twice produces 5,232 plot0 and 5,841 plot1 diagnostic-ghost differences in the production-main cases, with zero defined-state differences. The candidate all-access and narrow-access executables are also byte-identical. [Build identities](t23-build-identity.csv) therefore rule out an access-flag arithmetic change locally; source/build layout is not needed to generate the mismatch. The unused extra object can change undefined ghost output too, without changing any checked evolved or valid diagnostic value. The demonstrated cause is the unwritten buffer, while the precise allocator history of the original cluster process is not established.

All ordinary small-grid cases use levels 0–1, N=64x32, L=8 M, center=(4,0), h0=0.125 M, dt0=0.03125 M, one coarse step, the same E data, alpha_K, unchanged ExperimentalGauge/KO/point transfers and the positive-time static guard. If the ordinary small grid has no raw diagnostic difference, the matrix additionally initializes the real E-high max-level-14 hierarchy in both old/new dense variants with zero advances. It never launches an expensive max-level-14 coarse step. A separate private, environment-parametric buffer probe seeds only the newly allocated plot buffer with two distinct finite values. If these appear in saved diagnostic ghosts but nowhere in valid constraints/evolved data/checkpoints, that demonstrates the unwritten-output mechanism without subtracting or evaluating static data after initialization.

{chr(10).join(table)}

Across completed comparisons: evolved bit mismatches **{q['evolved_bit_mismatches']}**, valid diagnostic bit mismatches **{q['valid_diagnostic_bit_mismatches']}**, poison-dependent saved differences **{q['poison_differences']}**. The 20 matrix comparison cases cover 9,587,200 evolved Float64 value comparisons and 768,000 valid-diagnostic comparisons, with zero mismatches. These are array identities, not whole HDF5 file identities: the candidate intentionally adds gauge metadata. All 18 matrix native runs and two supplementary split-box poison runs produced both checkpoints and plots at t=0 and t=0.03125 M; [run validation](t23-run-validation.csv) seals the four outputs of each case. [t23-comparisons.csv](t23-comparisons.csv) records every field and region including zeros; [t23-location-summary.csv](t23-location-summary.csv) records differing boxes, levels, categories and maxima. The full per-cell ledger remains under `/private/tmp/ems-t23/ledgers/`; it is not committed. Physical-boundary/axis ghosts, uncovered refinement ghosts and covered same-level ghosts are distinguished using each saved hierarchy, not a chosen radial cut. Checkpoint0, plot0, checkpoint1 and plot1 are all compared; no early failure stops the later comparisons.

The [poison witness](t23-poison-witness.csv) is exact: at **each** plot clock all 6,480 changed values are 1.25 in one run and 2.5 in the other (difference exactly 1.25), in all six diagnostic fields. There are 612 coarse physical-boundary halo cells per field and, on level 1, 162 physical/axis halo cells plus 306 uncovered-refinement halo cells per field. The coarse valid box is i=0…63, j=0…31; the fine valid box is i=40…87, j=0…23 with physical faces x=±1.5 M, y=0…1.5 M. The poison occupies only the three-cell outside halos: coarse i=−3…66, j=−3…34 outside the valid box; fine axis j=−3…−1 and the lateral/top halos outside its valid box. No valid puncture cell changes. Ordinary old/new production plots have 3,144 ghost differences at t=0 and 6,177 at t=0.03125 M; the corresponding old/new dense plots happen to have zero differences at both clocks. Accidental equality of unwritten ghosts in one case is not a defined-output guarantee.

The smallest matrix grid has one box on each level and therefore has no internal same-level box interface. A supplementary, independently receipted poison pair uses the same private executable and grid with `max_box_size=32` to create internal interfaces. [Split-box comparisons](t23-split-comparisons.csv) and [qualification](t23-split-qualification.json) check **{q.get('split_probe',{}).get('same_level_diagnostic_ghost_values','pending')}** covered same-level diagnostic ghost values with **zero** poison differences. Its **{q.get('split_probe',{}).get('poisoned_values','pending')}** changed values over both plots again consist exclusively of exact 1.25/2.5 pairs in uncovered physical/refinement halos; all **{q.get('split_probe',{}).get('defined_value_comparisons','pending')}** defined-value comparisons are identical. Thus an internal box edge with neighbouring valid support is distinguished experimentally from an uncovered level face or physical/axis boundary. Its per-cell ledger is `/private/tmp/ems-t23/ledgers/split-cell-mismatches.csv`.

The production evolution RHS consumes evolved state, not the private plot buffer. An unwritten plot diagnostic ghost cannot feed it. This structural statement does not waive bitwise checks: any actual evolved or valid-diagnostic mismatch in the local matrix is reported as requiring further diagnosis. Undefined ghost output can occur at positive-time plots as well as t=0; whether it happens in the ordinary cases is a measured result in the tables. Defined valid constraints, their native formulas and term scales must be distinguished from these output halos.

For exp-0025, do not compare undefined diagnostic ghost values or interpolate them into the normalized magnitude gate. The registered masks, clocks and thresholds do not change. The exp-0023 reducer already strips output ghosts in `phase-3/audit.py:57–66`. Its existing phase-4 `EMSConstraintPieces.cpp` / `NativeConstraintPieces.impl.hpp` / `NativeGaussPieces.hpp`, with `production/pieces.py`, recomputes native valid-cell totals and individual term scales from saved numerical states without advances or static reads. `pieces.py:50–55` explicitly crops saved halos and checks evolved valid bits; lines 60–68 validate native totals against saved valid diagnostics and term sums with the existing 2048-epsilon scaled arithmetic allowance. Lines 70–92 use AMR-uncovered valid cells, the frozen masks, coordinate-volume weights, and the largest individual term's RMS on that same mask and clock. The T19 magnitudes consume these native tables; T21 supplies gauge-labelled checkpoint replay and current-state RHS/geometry, not a replacement normalized-constraint reducer by itself. One common, source/compiler-pinned native constraint tool can replay H and G at all 13 clocks, checking untouched evolved valid bits and the inherited normalization and Gauss shell integral. Existing H replay validation is evidence for H, not a waived validation for G. The reducer must fill private evolved stencils consistently and use valid/common-supported native constraints; no plot-ghost constraint is an input. Gauge metadata must be honoured when using a gauge-aware restart wrapper, although the native constraint formulas do not depend on the driver package. This preserves comparability without changing the screen or declaring the new gauge accurate.

The old writer is defective for saved diagnostic halos; the old and new valid constraints are the defined numerical quantities. Proposed controller action: retain strict bitwise evolved checkpoint/plot comparisons, including their defined ghosts, and strict bitwise comparison of all six diagnostic fields on every **valid** box cell and covered same-level ghost at both t=0 and one step. Check layouts, clocks, component names and gauge metadata separately, recognising the intentional new metadata. Report diagnostic ghosts outside same-level valid support as undefined rather than passing or failing them numerically. The cluster should classify its original differing indices before changing admission: the collected per-box receipt cannot establish those indices, so this local result does not itself lift that blocker. A future diagnostic-output repair would need explicit ghost definitions and separate qualification; copying uncomputed diagnostic source ghosts or silently zeroing them would not establish meaningful constraints. No production patch is applied, especially none to the in-flight exp-0024 fork d507438. [t23-plot-buffer-probe.patch](t23-plot-buffer-probe.patch) is private instrumentation only, not a proposed evolution or admission-check change.

Measured per-process RSS maximum: **{q['peak_RSS_bytes']} bytes** (0.611 GB, below 6 GB). All builds/runs are serial with two OpenMP threads and single-threaded numerical libraries, a 6 GB per-process cap, 7 GB process-tree cap and 16 GB output gate. Every completed child needs its own exit/resource receipt and native output checks; a launcher exit is insufficient. The original `/private/tmp/ems-t23/pipeline/done.exit = 1` remains as the failed-build history; the successful continuation is `/private/tmp/ems-t23/continuation/done.exit = 0`. No SSH, cluster work or commit is performed. Regenerate numerical comparisons with `python Tests/EMSNative/t23-compare.py analyze`; once the continuation and qualification receipts are complete, `python Tests/EMSNative/t23-compare.py finalize` checks all own receipts, refreshes the report/README and seals the packet. Tests instrumentation and builds are confined to `/private/tmp/ems-t23/`.
'''
    (HERE/'t23-default-path-diagnosis.md').write_text(text)
    readme=HERE/'README.md';s=readme.read_text();title='## T23 — default-path constraint diagnostic identity'
    if title in s:s=s[:s.index(title)].rstrip()
    readme.write_text(s+'\n\n'+title+'\n\n'+f'Local source/build and valid/ghost diagnosis: **{q["status"]}**. [Report](t23-default-path-diagnosis.md), [isolated preparer](t23-prepare.py), [receipt-first comparator](t23-compare.py), [source hashes](t23-source-inputs.json), [cluster receipts](t23-cluster-native-receipts.csv) and [manifest](t23-manifest.txt). The unchanged old plot writer leaves uncovered diagnostic ghosts unwritten. The same executable content gives differing ghost output, and two parametric poison values survive in exactly 6,480 diagnostic ghosts per plot; 9,587,200 evolved and 768,000 valid diagnostic comparisons have zero bit mismatches. This occurs at t=0 and one step. The failed narrow-access build was repaired in Tests compilation only; the completed five builds were reused. Every job is qualified from its own receipt, peak RSS 0.611 GB. The registered exp-0025 screen and its admission blocker are unchanged pending controller classification of the original cluster indices; the proposed check compares all defined evolved values and valid diagnostics strictly. No production source or check is modified; exp-0024 stays untouched. Successful marker: `/private/tmp/ems-t23/continuation/done.exit`.\n')
def seal():
    paths=[p for p in HERE.glob('t23-*') if p.is_file() and p.name!='t23-manifest.txt']+[HERE/'README.md',HERE/'t13-run.py',HERE/'t21-recorder.patch',HERE/'t21-recorder.hpp',FINDING]
    for p in ROOT.rglob('*'):
        if p.is_file() and (p.name in ('plan.json','params.txt','commands.json','build-spec.json','decision.json','repair.json','done.exit') or p.name.endswith('.resources.json') or p.name.endswith('.child.json') or p==ROOT/'builds/candidate-native-flags.failed-private-base/run.log'):paths.append(p)
    (HERE/'t23-manifest.txt').write_text('# T23 local-only diagnostic identity investigation; no commit\n# SHA256 bytes absolute_path\n'+''.join(f'{sha(p)}  {p.stat().st_size}  {p}\n' for p in sorted(set(paths))))
if __name__=='__main__':
    if sys.argv[1]=='pending':cluster();report();seal()
    else:globals()[sys.argv[1]]()
