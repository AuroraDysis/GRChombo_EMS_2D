#!/usr/bin/env python3
"""T17 registration and serial local evidence; never starts the merger."""
import sys
sys.dont_write_bytecode = True
import csv, hashlib, importlib.util, json, shutil, subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent; REPO=HERE.parents[1]
ROOT=Path('/private/tmp/ems-t17'); PY='/Users/auroradysis/miniconda3/bin/python'
PROFILE=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo-binary/a0.9-e8.trumpet')
CTT=Path('/Users/auroradysis/Workspace/EMS/.data/echo-binary/production-outer16-a20-r32-a10/n32-r6.ctt')
PRODUCTION=REPO/'Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex'
RH=.0060404520035922523
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
p=module('prepare',HERE/'t13-prepare.py')
def save(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def base(level,L=672,n=768):
    text=(HERE/'t16-L13-geometric.txt').read_text()
    changes=dict(N1=n,N2=n//2,L=L,center=f'{L/2:g} 0',star_centre=f'{L/2:g} 0',
        ems_data_path=PROFILE,ems_ctt_data_path=CTT,binary='true',boosted='true',
        bh_mass=1,bh_charge=1.0693100593480418,separation=32,boost_rapidity=.05238600881798246,
        max_level=level,regrid_interval=' '.join(['0']*level),max_steps=0,stop_time=0,
        ems_use_maximal_initial_lapse='false',ems_use_geometric_initial_lapse='true',
        t14_dense_initial_tags='false',t13_launch_stop_time=0,
        checkpoint_interval=-1,plot_interval=0,num_plot_vars=28,
        plot_vars='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi',
        mass_extraction_center=f'{L/2:g} 0',num_mass_extraction_radii=level,
        mass_extraction_levels=' '.join(str(i) for i in range(1,level+1)),
        mass_extraction_radii=' '.join(format(112/2**i/1.2,'.17g') for i in range(1,level+1)),verbosity=1)
    return p.replace(text,changes)
def register_first():
    ROOT.mkdir(exist_ok=True)
    assert digest(PROFILE)=='4d7cc0b975d4dde0c4da8400ef0bbad8c3e57716c9708c160d19e56288a76f90'
    assert digest(CTT)=='46253cb8d3bf6255b149b44bed7d4899cd97c697fe677a9846649bfbb2bdc4ee'
    frozen=ROOT/'production-77c5b6f.ex'; frozen.write_bytes(PRODUCTION.read_bytes()); frozen.chmod(0o755)
    (HERE/'t17-inputs.json').write_text(json.dumps([dict(path=str(f),sha256=digest(f),bytes=f.stat().st_size)
        for f in (PROFILE,CTT,frozen)],indent=2)+'\n')
    jobs=[]
    for level in (11,12):
        name=f'production-original-L{level}'; d=ROOT/'census'/name
        for sub in ('plt','chk'):(d/sub).mkdir(parents=True,exist_ok=True)
        param=HERE/f't17-{name}.txt'; param.write_text(base(level))
        jobs.append(dict(name=name,directory=str(d),command=[str(frozen),str(param)]))
    (ROOT/'census/first-plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    note='''## T17 — registered low-resolution loss-measurement design (merger not run)

Baseline: fork main 77c5b6f. The supplied published EMSTRUMPET/1 and EMSCTT/1 files are pinned in [t17-inputs.json](t17-inputs.json). The companion remains DOCUMENTED DIAGNOSTIC DATA with C1 panel interfaces, not data admitted to a convergence campaign. No SSH, commit or merger is authorized here.

**Registration before numerical tests.** The first numerical jobs use the saved production executable, its original DefaultLevelFactory and native tagging, with no Tests factory. Candidate finest levels are 11 and 12 at h0=0.875 M_i: h=0.00042724609375 and 0.000213623046875 M_i, R_h/h=14.138 and 28.276; dt0=0.21875 M_i, dt_finest=h/4. This first census also checks whether the existing single extraction centre actually refines both holes. The native EMS criterion has one extraction centre; fixed concentric tags at the binary midpoint cannot resolve the holes at high levels. Any necessary opt-in dual-centre tagging must therefore be disclosed and tested on production main. Gradient thresholds remain disabled for the controlled initialization/launch hierarchy. The intended per-hole tag radii are 112/2^level M_i before native buffer/block rounding. Actual box faces, cell counts and puncture coverage must come from the census, never an assumed factory geometry. If the production TreeIntVectSet fails at the tested hierarchy, only the existing T14 dense representation workaround is permitted, with actual tag/box/state identity evidence on a working rung. A failed initialization or RSS gate is a failure receipt.

The resolution decision is frozen before seeing the new results. For each rung and each hole, the T16b launch screen uses W=[0.00075,0.0025] M_i, inward axis/diagonal rays, numerical-t0 subtraction, central I8, spread |I8-I10|, peak/trapezoidal RMS and five-spread qualification. The common endpoint is 18 coarse finest steps, t=0.001922607421875 M_i, with all 19 common clocks including t=0; the fine rung uses even steps through 36. Native stops are 18 and 37 steps. Both sqrt(chi) and geometric variants use unchanged native RK4, gauge, KO=1, point transfers and floors. A launch PASS requires both Gamma norms on both rays to exceed five interpolation spreads in both variants, their conservative five-spread ratio upper bounds below one, and smaller total native puncture Gamma RHS. The other three fields, every history entry, sampling failures, early enhancements and floors are reported, not erased. Metric Gamma qualification is reported separately. No same-grid dt/2 control is claimed: spatial and temporal refinement change together.

The unchanged RHFinder is tested on each geometric-lapse hierarchy at numerical t=0 and after one synchronized coarse step (t=0.21875 M_i). N_theta=48 and 96 use the existing three-stage thresholds expansion_squared <=1e-7,1e-10,1e-12, with fresh final interpolation. Both individual surfaces must be FOUND at the last threshold at both times; N48/N96 area and charge differences must be <=1e-3 fractionally. The squared residual and its square root, stopping differences, angular differences and caps are retained. Failed/capped probes cannot qualify a rung. Minimum-cost recommendation requires BOTH launch and finder PASS, and measured RSS below 6 GB. Local timing uses three synchronized coarse steps with extraction/finder/plot/checkpoint work excluded from the step interval, on both binary and the same-machine single-hole control at the eventual recommended rung; rate claims require actual step clocks, not whole-run initialization wall time.

Run 1 measures losses: approximate infall 306 M_i followed by 100 M_f after the first common horizon, restartable to 200–300 M_f. The notebook's 2026-10-01T17:55Z and 2026-10-02T01:10Z readings apply: zero-loss tuning near e_f=7 is not a promised endpoint; losses shift e upward, and the retuned second run carries the +/-0.5 admission test. Continue if the late 50 M_f window has not reached delta e_H<0.1 with settled exterior diagnostics. The merger parameter file is a DRAFT until the resolution evidence is complete.

All local processes are serially watched and measured: 6e9 bytes RSS/footprint per process, 6.5e9 tree gate, OMP=2 and BLAS=1, <=4 threads, 16 GB output ceiling. Long evidence work is detached with atomic done.exit markers; exit zero alone is not acceptance. Static input is initialization-only (t2 guard); capture only writes native state. The author's RHFinder, equations, gauge, KO, transfers, Float64 and puncture treatment remain unchanged.
'''
    readme=HERE/'README.md'; text=readme.read_text(); assert '## T17 —' not in text
    readme.write_text(text.rstrip()+'\n\n'+note)
    print(ROOT/'census/first-plan.json',flush=True)
def measured(label,directory,command):
    directory.mkdir(parents=True,exist_ok=True)
    for sub in ('plt','chk'):(directory/sub).mkdir(exist_ok=True)
    with (directory/'run.log').open('wb') as log:
        rc=subprocess.run([PY,str(HERE/'t17-run.py'),'--measure',label,'--directory',str(directory),
            '--',*map(str,command)],stdout=log,stderr=subprocess.STDOUT).returncode
    (directory/'done.exit').write_text(str(rc)+'\n')
    return rc
def study_base(level):
    faces={1:224.,2:144.}
    return p.replace(base(level,4480,2560),dict(ems_binary_refinement='true',
        ems_track_punctures='false',ems_puncture_tracking_level=6,
        checkpoint_interval=0,plot_interval=-1,
        mass_extraction_radii=' '.join(format(faces.get(i,224/2**i)/1.2,'.17g') for i in range(1,level+1)),verbosity=0))
def register_study():
    assert all((ROOT/'census'/f'production-original-L{l}'/'done.exit').read_text().strip()=='0' for l in (11,12))
    initial=[]
    for l in (11,12):
        d=ROOT/'census'/f'production-original-L{l}';r=json.loads((d/f'production-original-L{l}.resources.json').read_text())
        assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason']
        assert 'GRChombo finished.' in (d/'run.log').read_text()
        initial.append(dict(run=d.name,production_main=True,original_factory=True,initialized=True,
            maximum_level=l,requested_h_M=.875/2**l,peak_RSS_bytes=r['peak_rss_bytes'],wall_seconds=r['wall_seconds'],
            puncture_coverage='not measured by this first exit-only census; fixed midpoint tags do not ensure requested puncture resolution'))
    save('t17-first-census.csv',initial)
    (ROOT/'production-t17.ex').write_bytes(PRODUCTION.read_bytes());(ROOT/'production-t17.ex').chmod(0o755)
    # One representation change in the existing dense workaround, chosen before
    # candidate results: the larger global indices need dense tags from level 11.
    dense=(HERE/'T14DenseTags.hpp').read_text();assert dense.count('m_level<13')==1
    (HERE/'t17-dense.hpp').write_text('// T17 large-domain variant: identical T14 tagging, dense from level 11.\n'+dense.replace('m_level<13','m_level<11'))
    (HERE/'t17-dense-main.cpp').write_text('#include "t17-dense.hpp"\n#define DefaultLevelFactory T14LevelFactory\n#include "../../Examples/EMS/Main_EMSBH2DBH.cpp"\n#undef DefaultLevelFactory\n')
    main=(REPO/'Examples/EMS/Main_EMSBH2DBH.cpp').read_text().replace('DefaultLevelFactory<EMSBH2DLevel>','T17Factory<EMSBH2DLevel>')
    needle='    using Clock = std::chrono::steady_clock;';assert main.count(needle)==1
    main=main.replace(needle,'    t17_initial_audit(bh_amr,sim_params);\n'+needle)
    (HERE/'t17-main.hpp').write_text(main)
    rh=(REPO/'Tests/EMSRHFinder/harness/EMSRHCheckpoint.cpp').read_text()
    needle='time == 0. ? std::nextafter(0., 1.) : time);';assert rh.count(needle)==1
    rh=rh.replace(needle,'0.); // T17 fresh spheres on the current frozen checkpoint, not history replay')
    needle='        s = resample(s, n);';assert rh.count(needle)==1
    rh=rh.replace(needle,needle+'\n        s.m_chase_speed = p.m_RH_chase_speeds.at(s.m_index); // native configured parameter\n        s.m_time_step_freq = p.m_RH_time_step_freq.at(s.m_index); // native per-call quota')
    (HERE/'t17-rh.cpp').write_text('// T17 diagnostic wrapper; RHUnion/RHSurf/interpolator are unchanged.\n'+rh)
    rows=[]
    for rung,level in [('coarse',12),('fine',13)]:
        h=1.75/2**level
        for mode in ('sqrt','geometric'):
            name=rung+'-'+mode;param=HERE/f't17-{name}.txt'
            param.write_text(p.replace(study_base(level),dict(ems_use_geometric_initial_lapse=str(mode=='geometric').lower())))
            rows.append(dict(run=name,base_N1=2560,base_N2=1280,L_M=4480,h0_M=1.75,dt0_M=.4375,
                finest_level=level,h_M=h,Rh_over_h=RH/h,dt_finest_M=h/4,native_steps=int(.002/(h/4)),
                common_clocks=19,common_stop_M=18*.875/2048/4,
                face_1_tag_M=224,face_2_tag_M=144,inner_face_tag_M=224/2**level,
                actual_faces='pending native census; buffer/block rounding applies',parameters=str(param)))
    save('t17-registration.csv',rows)
    path=HERE/'README.md';path.write_text(path.read_text()+'''

### Domain/grid amendment before two-centre results

The first original-production jobs have actual child exits 0 and no gate: max11 85.09 s / 0.332562 GB RSS, max12 94.46 s / 0.389054 GB RSS. Their main/factory are unchanged. Their plot_interval=0 did not retain plots because setup leaves Chombo's disabled default in that case; these receipts establish initialization only, not box coverage. The subsequent original/default-off identity control retains full initial plots and final checkpoints.

Before any dual-centre numerical result, the production study domain is frozen to L=4480 M_i, N1=2560, N2=1280, center=(2240,0), holes (2224,0)/(2256,0), h0=1.75 M_i, levels 0–12 and 0–13. Thus the originally declared finest h, R_h/h and launch clocks are identical, while dt0=0.4375 M_i and the synchronized short-evolution horizon clock is now 0.4375 M_i. The two timing samples are two full coarse steps, at 0.4375 and 0.875 M_i (replacing the original three-step proposal before timing results). Per-hole tag faces are 224 and 144 M_i on levels 1/2, then 224/2^l for l>=3; actual native buffer/block-rounded faces are recorded. These broad first two levels resolve all 50/75/100 M_i extraction spheres at h2=0.4375 M_i without a separate wave tagging variant. The opt-in production ems_binary_refinement copies the existing criterion at both centres; ems_track_punctures reuses the existing shift-integrating PunctureTracker. Both default off. The native factory is tested first at each actual candidate level. Only an observed native TreeIntVectSet failure enables the existing T14 dense representation on the same cells; the large-domain T17 copy changes its threshold from level 13 to level 11 and requires whole-hierarchy initialization identity on a working candidate. There is no new fixed-grid Tests factory.

The moving-patch evolution uses native regrid intervals 4 on levels 0–5, 1 on level 6 and 0 above, tracking on level 6. The short frozen-face launch keeps regridding/tracking off. The timing hook begins at native level-0 advance and ends after its complete post-step: it excludes initialization and plot/checkpoint writes. It includes native lower-level regridding and puncture tracking; extraction and live finder are disabled. Coarsest-level tag work before advance is not included and must be budgeted separately for the long run. RHFinder probes use the author's unchanged implementation and fresh spherical seeds at the recorded native puncture positions, chase_speed=1 (a native parameter), 48/96 points, the original three-stage residual schedule, and 1795 s per case. TIME_CAP/FLOOR/failed probes cannot qualify a rung. This is a predeclared resource bound, not evidence of horizon convergence.

Boundary design uses an assumed characteristic envelope V<=2 throughout evolution, checked against the actual native metric/lapse/shift census and monitored during the future merger. The outer radius is 2240 M_i. It protects the largest extraction sphere even against a boundary disturbance present at t=0 until (2240-100)/2=1070 M_i. This is stricter than a pulse round trip. The initial-M_ADM upper clock 2.002174604 M_i gives planned run-1/300-M_f horizons 506.2174604/906.6523812 M_i using the 306-M_i estimate. Choosing the larger domain initially permits continuation on an unchanged restart hierarchy. If measured speeds violate V=2, or infall/continuation lasts beyond this budget, the draft fails its boundary criterion; no new assumption is substituted after results.
''')
    jobs=[dict(name='local-evidence',directory=str(ROOT/'pipeline/job'),command=[PY,str(HERE/'t17-work.py'),'pipeline'])]
    (ROOT/'pipeline').mkdir(exist_ok=True)
    (ROOT/'pipeline/plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    print('Actual-domain candidate study registered before native runs.',flush=True)
def build_harnesses():
    builder=module('builder',HERE/'t14-prepare.py');builder.ROOT=ROOT
    for source,exe in [('t17-launch.cpp','launch.ex'),('t17-dense-main.cpp','dense-production.ex'),('t17-rh.cpp','rh.ex')]:builder.build(source,exe)
    shutil.copyfile('/private/tmp/ems-t16b/replay.ex',ROOT/'replay.ex');(ROOT/'replay.ex').chmod(0o755)
def run_case(name,directory,exe,text):
    param=directory/'params.txt';directory.mkdir(parents=True,exist_ok=True)
    if (directory/'done.exit').exists():
        assert param.read_text()==text,f'Completed case parameter drift: {directory}'
        return int((directory/'done.exit').read_text())
    param.write_text(text)
    return measured(name,directory,[exe,param])
def controls():
    rows=[]
    for mode in ('sqrt','geometric'):
        text=p.replace(base(1,56,64),dict(max_steps=2,stop_time=1,checkpoint_interval=1,plot_interval=1,
            mass_extraction_center='28 0',ems_use_geometric_initial_lapse=str(mode=='geometric').lower(),verbosity=0))
        dirs=[]
        for label,exe in [('baseline',ROOT/'production-77c5b6f.ex'),('default-off',ROOT/'production-t17.ex'),('capture',ROOT/'launch.ex')]:
            d=ROOT/'controls'/f'{mode}-{label}';assert run_case(d.name,d,exe,text)==0;dirs.append(d)
        for current in dirs[1:]:
            for relative in ('plt/EMS_Plot_000000.2d.hdf5','chk/EMS_000002.2d.hdf5'):
                count,bits=p.compare(dirs[0]/relative,current/relative)
                rows.append(dict(mode=mode,variant=current.name,file=relative,Float64_values=count,bit_mismatches=bits));assert bits==0
    save('t17-controls.csv',rows)
def census():
    choices={}
    for rung,level in [('coarse',12),('fine',13)]:
        text=(HERE/f't17-{rung}-geometric.txt').read_text();d=ROOT/'census'/(rung+'-native')
        rc=run_case(d.name,d,ROOT/'production-t17.ex',text)
        if rc==0:choices[rung]=dict(dense=False,geometric_checkpoint=str(d/'chk/EMS_000000.2d.hdf5'))
        else:
            resource=json.loads((d/(d.name+'.resources.json')).read_text())
            if resource['gate_reason'] or resource['returncode'] not in (-11,139):
                raise RuntimeError(f'Native census failed without the registered tree-segfault signature: {d}; {resource}')
            d=ROOT/'census'/(rung+'-dense')
            assert run_case(d.name,d,ROOT/'dense-production.ex',p.replace(text,dict(t14_dense_initial_tags='true')))==0
            choices[rung]=dict(dense=True,geometric_checkpoint=str(d/'chk/EMS_000000.2d.hdf5'))
    # Qualify the dense representation against the same actual coarse hierarchy.
    if any(x['dense'] for x in choices.values()):
        assert not choices['coarse']['dense'],'No working actual rung for dense identity qualification'
        d=ROOT/'census/coarse-dense-control';text=(HERE/'t17-coarse-geometric.txt').read_text()
        assert run_case(d.name,d,ROOT/'dense-production.ex',p.replace(text,dict(t14_dense_initial_tags='true')))==0
        count,bits=p.compare(Path(choices['coarse']['geometric_checkpoint']),d/'chk/EMS_000000.2d.hdf5');assert bits==0
        save('t17-dense-identity.csv',[dict(rung='coarse',Float64_values=count,bit_mismatches=bits,scope='all hierarchy boxes, native t0 components and ghosts')])
    for rung in choices:
        d=ROOT/'census'/(rung+'-sqrt');text=(HERE/f't17-{rung}-sqrt.txt').read_text()
        text=p.replace(text,dict(t14_dense_initial_tags=str(choices[rung]['dense']).lower()))
        exe=ROOT/('dense-production.ex' if choices[rung]['dense'] else 'production-t17.ex')
        assert run_case(d.name,d,exe,text)==0
        choices[rung]['sqrt_checkpoint']=str(d/'chk/EMS_000000.2d.hdf5')
    (ROOT/'choices.json').write_text(json.dumps(choices,indent=2)+'\n');return choices
def launch(choices):
    for rung,choice in choices.items():
        for mode in ('sqrt','geometric'):
            name=rung+'-'+mode;d=ROOT/'evolution'/name
            text=p.replace((HERE/f't17-{name}.txt').read_text(),dict(restart_file=choice[mode+'_checkpoint'],
                t14_dense_initial_tags=str(choice['dense']).lower(),checkpoint_interval=-1,
                t13_launch_stop_time=.002,max_steps=2,stop_time=1))
            assert run_case(name,d,ROOT/'launch.ex',text)==0
def timed_binary(choices,rung):
    level=12 if rung=='coarse' else 13;d=ROOT/'timing'/(rung+'-binary')
    text=p.replace((HERE/f't17-{rung}-geometric.txt').read_text(),dict(restart_file=choices[rung]['geometric_checkpoint'],
        t14_dense_initial_tags=str(choices[rung]['dense']).lower(),checkpoint_interval=1,
        max_steps=2,stop_time=1,ems_track_punctures='true',
        regrid_interval=' '.join('4' if l<6 else '1' if l==6 else '0' for l in range(level))))
    return run_case(d.name,d,ROOT/'launch.ex',text)
def horizons(rung):
    d=ROOT/'timing'/(rung+'-binary');text=(d/'params.txt').read_text();report=[]
    for step in (0,1):
        cp=Path(json.loads((ROOT/'choices.json').read_text())[rung]['geometric_checkpoint']) if step==0 else d/'chk/EMS_000001.2d.hdf5'
        centres=[2224.,2256.]
        if step:
            records=[line.split() for line in (d/'punctures.dat').read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
            row=next(x for x in records if abs(float(x[0])-.4375)<1e-8);centres=[float(row[1]),float(row[3])]
        for n in (48,96):
            target=ROOT/'finder'/f'{rung}-step{step}-n{n}'
            params=p.replace(text,dict(restart_file=cp,RH_activate='false',RH_num_horizons=2,
                RH_initial_centre=' '.join(map(str,centres)),RH_initial_radii=f'{RH:.17g} {RH:.17g}',
                RH_num_points=f'{n} {n}',RH_level='0 0',RH_start_times='0 0',
                RH_time_step_freq='1 1',RH_chase_speeds='1 1',RH_newton_crit='0 0',
                offline_points=n,offline_seconds=1795,offline_max_updates=2000,offline_floor_window=64,
                checkpoint_interval=-1,ems_track_punctures='false',t13_launch_stop_time=0))
            rc=run_case(target.name,target,ROOT/'rh.ex',params)
            report.append(dict(rung=rung,step=step,time_M=step*.4375,N_theta=n,returncode=rc,directory=str(target)))
            if rc not in (0,1):raise RuntimeError(f'Finder process/infrastructure failure: {target}; exit {rc}')
    save(f't17-{rung}-finder-receipts.csv',report)
def timed_single(choices,rung):
    level=12 if rung=='coarse' else 13;d=ROOT/'timing'/(rung+'-single')
    text=p.replace(study_base(level),dict(binary='false',ems_ctt_data_path='""',separation=0,
        ems_binary_refinement='false',ems_track_punctures='false',max_steps=2,stop_time=1,
        t14_dense_initial_tags=str(choices[rung]['dense']).lower(),checkpoint_interval=-1,plot_interval=-1,
        regrid_interval=' '.join('4' if l<6 else '1' if l==6 else '0' for l in range(level))))
    assert run_case(d.name,d,ROOT/'launch-rate.ex',text)==0
def build_rate():
    builder=module('builder',HERE/'t14-prepare.py');builder.ROOT=ROOT
    builder.build('t17-launch.cpp','launch-rate.ex')
def pipeline():
    for item in json.loads((HERE/'t17-inputs.json').read_text()):
        assert digest(Path(item['path']))==item['sha256'],f'Input changed: {item["path"]}'
    controls();choices=census();launch(choices)
    analysis=module('analysis',HERE/'t17-analyze.py');analysis.census();analysis.launch()
    for rung in ('coarse','fine'):
        if timed_binary(choices,rung)!=0:raise RuntimeError(f'{rung} short evolution failed')
        horizons(rung)
    recommendation=analysis.resolution()
    # Price a qualified rung; if neither qualifies, the fine rung is a cost
    # reference only and remains explicitly unadmitted.
    timed_single(choices,recommendation or 'fine')
    analysis.finish();print('T17 local evidence finished; merger remains a draft.',flush=True)
def draft():
    changes=dict(ems_track_punctures='true',max_steps=32,stop_time=506.625,checkpoint_interval=8,plot_interval=-1,
        regrid_interval='4 4 4 4 4 4 1 0 0 0 0 0 0',activate_extraction=1,
        extraction_center='2240 0',num_extraction_radii=3,extraction_radii='50 75 100',extraction_levels='2 2 2',
        num_points_theta=97,num_points_phi=64,num_modes=5,modes='2 0 3 0 4 0 5 0 6 0',
        ems_radiation_activate='true',ems_radiation_phi_inf=0,ems_radiation_wave_level=0,ems_radiation_wave_radius=0,
        activate_mq_extraction=1,mq_extraction_center='2240 0 0',mq_num_extraction_radii=3,
        mq_extraction_radii='50 75 100',mq_extraction_levels='2 2 2',mq_num_points_theta=97,mq_num_points_phi=64,
        mq_num_modes=1,mq_modes='0 0',activate_rs_extraction=0,activate_em_extraction=0,
        RH_activate='false',RH_num_horizons=2,RH_initial_radii=f'{RH:.17g} {RH:.17g}',
        RH_initial_centre='2224 2256',RH_num_points='48 48',RH_level='0 0',RH_start_times='0 0',
        RH_time_step_freq='1 1',RH_chase_speeds='1 1',RH_newton_crit='0 0',verbosity=0)
    header='''# T17 DRAFT ONLY: fork 77c5b6f plus t17-production.patch; merger not run.
# Fine level is a provisional placeholder until t17-resolution.json qualifies a rung.
# Existing radiation outputs: GW, coupling-weighted EM and scalar including l=0.
# Needed write-only hooks: full-metric ADM surface mass and t0 MQ; fixed-mask
# constraints/speeds; bounded checkpoint horizon driver and segment completion checks.
# Native RH calls are off: RH_time_step_freq is an iteration quota, not cadence.
# max_steps=32 is the FIRST SEGMENT limit. Raise the absolute limit on restart.
# Set stop_time from ACTUAL first qualified common horizon +100*Mref, Mref=2.002174604.
# Map pinned published inputs to cluster-local copies before any future build/run.
# t14_dense_initial_tags requires the qualified DENSE MAIN if native census fails.
'''
    text=p.replace(study_base(13),changes)
    text=text.replace('# Charge extraction at 20, 50 and 100 M; radiative extractions off.',
        '# Radiation and charge extraction at the registered 50, 75 and 100 M_i spheres.')
    (HERE/'t17-merger-draft.txt').write_text(header+text)
    path=HERE/'t17-rh.cpp';text=path.read_text();needle='        s.m_chase_speed = p.m_RH_chase_speeds.at(s.m_index); // native configured parameter'
    if 's.m_time_step_freq = p.m_RH_time_step_freq' not in text:
        assert text.count(needle)==1
        path.write_text(text.replace(needle,needle+'\n        s.m_time_step_freq = p.m_RH_time_step_freq.at(s.m_index); // native per-call quota'))
if __name__=='__main__':
    {'register-first':register_first,'register-study':register_study,'build':build_harnesses,
        'controls':controls,'pipeline':pipeline,'draft':draft,'build-rate':build_rate}[sys.argv[1]]()
