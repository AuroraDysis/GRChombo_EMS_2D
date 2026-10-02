#!/usr/bin/env python3
"""Register and run the d16 T16b audit using native kernels and existing helpers."""
import sys
sys.dont_write_bytecode=True
import csv, hashlib, importlib.util, json, subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent; REPO=HERE.parents[1]
ROOT=Path('/private/tmp/ems-t18c'); PY='/Users/auroradysis/miniconda3/bin/python'
DATA=Path('/Users/auroradysis/Workspace/EMS/.data/echo-binary/tranche3/d16/r32-a10/n32-r6.ctt')
PROFILE=Path('/Users/auroradysis/Workspace/EMS/artifacts/echo-binary/a0.9-e8.trumpet')
DOC=Path('/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-echobin/artifacts/echo-binary/tranche3/d16/r32-a10')
SHA='107370ae00d8b69dc3122dc023e34bbff23b90afc8401bb233743c1182e6083c'
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
params=module('t18cparams',HERE/'t13-prepare.py')
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def measured(name,directory,command):
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/(name+'.log')).open('wb') as log:
        subprocess.run([PY,str(HERE/'t18c-run.py'),'--measure',name,'--directory',str(directory),'--',*map(str,command)],stdout=log,stderr=subprocess.STDOUT,check=True)
def generate():
    ROOT.mkdir(exist_ok=True)
    assert digest(DATA)==SHA
    inputs=[DATA,PROFILE,DOC/'params-initial-data.txt',DOC/'audit.md',HERE/'t16b-native.hpp',HERE/'t16b-analyze.py']
    (HERE/'t18c-inputs.json').write_text(json.dumps([dict(path=str(p),sha256=digest(p),bytes=p.stat().st_size) for p in inputs],indent=2)+'\n')
    census=(HERE/'t18-census.cpp').read_text().replace('finest=13','finest=12').replace('{16,24,32,48,128}','{24}')
    census=census.replace('2224','2232').replace('2256','2248')
    (HERE/'t18c-census.cpp').write_text(census)
    # Keep the T16b numerical audit verbatim except configuration/geometry.
    native=(HERE/'t16b-native.hpp').read_text().replace('T16B','T18C').replace('t16b','t18c')
    native=native.replace('-16.','-8.').replace('16.','8.')
    native=native.replace('std::abs(x)>3.', 'std::abs(x)>3*p.coarsest_dx')
    start=native.index('inline double t18c_coarse_h('); end=native.index('\ninline void t18c_initial_audit',start)
    geometry='''inline std::vector<std::unique_ptr<DenseIntVectSet>> t18c_cover;
inline double t18c_mid=0.;
inline double t18c_coarse_h(double x,double y)
{
    for(int l=int(t18c_cover.size())-1;l>0;--l)
    {
        double h=std::ldexp(1.75,-l);
        IntVect iv(int(std::floor((x+t18c_mid)/h)),int(std::floor(std::abs(y)/h)));
        if(t18c_cover[l]->contains(Box(iv,iv)))return h;
    }
    return 1.75;
}
'''
    native=native[:start]+geometry+native[end:]
    needle='    auto iso=p.emsbh_params;'
    native=native.replace(needle,'''    t18c_mid=p.center[0];t18c_cover.clear();
    const auto all_levels=amr.getAMRLevels();
    for(int l=0;l<all_levels.size();++l)
    {
        auto *level=dynamic_cast<T18CLevel *>(all_levels[l]);
        const auto &layout=level->getLevelData().disjointBoxLayout();Box bounds;
        for(LayoutIterator it=layout.layoutIterator();it.ok();++it)bounds.minBox(layout[it()]);
        auto cover=std::make_unique<DenseIntVectSet>(bounds);cover->makeEmptyBits();
        for(LayoutIterator it=layout.layoutIterator();it.ok();++it)(*cover)|=layout[it()];
        t18c_cover.push_back(std::move(cover));
    }
'''+needle)
    needle='    void tagCells(IntVectSet &tags) override'
    native=native.replace(needle,'''    void initialGrid(const Vector<Box> &boxes) override
    {
        GRAMRLevel::initialGrid(boxes);
        GRParmParse pp;bool audit=false,compact=true;
        pp.load("t18c_audit_only",audit,false);pp.load("t18c_compact_initial_audit",compact,true);
        // InitialData/fillAllEvolutionGhosts read only new-state arrays.
        // Ordinary-path smoke equality checks this memory-only optimization.
        if(audit&&compact){m_state_old.clear();m_state_diagnostics.clear();}
    }
'''+needle)
    (HERE/'t18c-native.hpp').write_text(native)
    panels=[]
    for line in DATA.open():
        if line.startswith('PANEL '):
            q=line.split();panels.append(dict(panel=int(q[1]),kind=q[2],lo=float(q[3]),hi=float(q[4]),center=float(q[5]),stretch=float(q[6]),mu_lo=float(q[7]),mu_hi=float(q[8]),boundary=q[9],nr=int(q[10]),na=int(q[11])))
    assert len(panels)==34
    with (HERE/'t18c-panels.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=panels[0]);w.writeheader();w.writerows(panels)
    registration=dict(baseline='dd84403',companion_sha256=SHA,companion_status='documented diagnostic only',separation_M=16,
        rapidity=0.05778205303580913,production_shape=[2560,1280],production_L_M=4480,h0_M=1.75,
        levels=list(range(13)),max_box_size=24,block_factor=8,tag_buffer=3,grid_buffer=8,fill_ratio=.7,
        h_finest_M=1.75/4096,Rh_over_h=.0060404520035922523/(1.75/4096),dt_multiplier=.5,
        requested_window_M=.002,native_clocks_M=[k*1.75/4096*.5 for k in range(10)],
        T16b_endpoint_M=.001983642578125,clock_caveat='The requested low rung/dt has 9 positive native clocks ending .001922607421875, not the original T16b 65 clocks. No time interpolation or alignment.',
        masks=dict(collars_M=[.01,.02],interhole='abs(x-midpoint)<=8, 0<=y<=8',far_M=[80,128],
            sheets='native panel selectors on +/-3 h_c stencils; h_c from actual LOW Chombo box union; AMR-face-excluded subset retained'),
        launch_W_M=[.00075,.0025],rays=['inward-axis','inward-diagonal'],holes=['left','right'],sampling='I8 with |I8-I10|, five-spread margins; native t0 subtracted',
        recommendation_rule='resolved endpoint evolved-Gamma reduction on both rays at both holes and smaller native total puncture Gamma source at both holes',
        convergence_admission=False,process_RSS_cap_bytes=4000000000,threads=2,
        launch_domain='Integer coarse-cell translation by 1792 M to L896; exact full-production refined boxes and cropped production base boxes. All captured t0 fields must match full-production audit before launch results qualify.',
        unchanged='equations, gauge, sigma=1, Float64, point transfers, puncture, RHFinder; static reader t0 only')
    (HERE/'t18c-registration.json').write_text(json.dumps(registration,indent=2)+'\n')
def build_census():
    old=module('t18cbuildcensus',HERE/'t18-work.py');old.ROOT=ROOT
    # Existing helper emits the native Chombo-only command, then pin its source.
    old.build('census')
    cmd=json.loads((ROOT/'census-build-command.json').read_text())
    cmd=[str(HERE/'t18c-census.cpp') if a.endswith('/t18-census.cpp') else a for a in cmd]
    (ROOT/'census-build-command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    subprocess.run(cmd,check=True)
    subprocess.run([str(ROOT/'census.ex'),str(HERE/'t18c-census')],check=True)
def prepare():
    boxes=list(csv.DictReader((HERE/'t18c-census-boxes.csv').open()))
    base=(HERE/'t18-launch-coarse-dt0.5.txt').read_text()
    data=params.parameters((DOC/'params-initial-data.txt').read_text());jobs=[]
    for scope in ('full','launch'):
        file=HERE/f't18c-{scope}-boxes.txt'; lines=[]
        for b in boxes:
            l=int(b['level']);x0,y0,x1,y1=(int(b[k]) for k in ('x0','y0','x1','y1'))
            if scope=='launch':
                shift=1024*2**l;x0-=shift;x1-=shift
                if l==0:
                    x0=max(x0,0);x1=min(x1,511);y1=min(y1,255)
                    if x0>x1 or y0>y1:continue
                assert 0<=x0<=x1<512*2**l and 0<=y0<=y1<256*2**l
            lines.append(' '.join(map(str,[l,x0,y0,x1,y1])))
        file.write_text('\n'.join(lines)+'\n')
        for mode in ('sqrt','geometric'):
            name=scope+'-'+mode;d=ROOT/'evolution'/name;d.mkdir(parents=True,exist_ok=True)
            changes=dict(data);changes.update(ems_data_path=PROFILE,ems_ctt_data_path=DATA,
                N1=2560 if scope=='full' else 512,N2=1280 if scope=='full' else 256,L=4480 if scope=='full' else 896,
                center='2240 0' if scope=='full' else '448 0',star_centre='2240 0' if scope=='full' else '448 0',
                mass_extraction_center='2240 0' if scope=='full' else '448 0',max_box_size=24,
                t18c_boxes=file,t18c_audit_only=str(scope=='full').lower(),t13_launch_stop_time=.002,
                ems_use_geometric_initial_lapse=str(mode=='geometric').lower(),
                t18c_compare_nonlapse=ROOT/'evolution'/(scope+'-sqrt')/'t18c-nonlapse.bin' if mode=='geometric' else '""')
            param=HERE/f't18c-{name}.txt';param.write_text(params.replace(base,changes))
            jobs.append(dict(name=name,directory=str(d),command=[str(ROOT/'launch.ex'),str(param)]))
    pipeline=ROOT/'pipeline';pipeline.mkdir(exist_ok=True)
    # Analysis runs only after all individual native receipts have completed.
    jobs.append(dict(name='analysis',directory=str(ROOT/'analysis-job'),command=[PY,str(HERE/'t18c-analyze.py')]))
    (pipeline/'plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    return pipeline/'plan.json'
def build_native():
    record=ROOT/'build-production/job/production-build.resources.json'
    if not list((REPO/'Examples/EMS/o').glob('**/EMSBH2DLevel.o')):
        measured('production-rebuild',ROOT/'production-rebuild',['make','-C',REPO/'Examples/EMS','all','DIM=2','-j1'])
        record=ROOT/'production-rebuild/production-rebuild.resources.json'
    r=json.loads(record.read_text())
    assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason'] and r['peak_rss_bytes']<4e9
    b=module('t18cnativebuilder',HERE/'t14-prepare.py');b.ROOT=ROOT;b.build('t18c-launch.cpp','launch.ex')
    cmd=json.loads((ROOT/'launch.ex-build-command.json').read_text())
    cmd=[a for a in cmd if not a.endswith('.o')]
    cmd=[str(HERE/'t16-replay.cpp') if a.endswith('/t18c-launch.cpp') else str(ROOT/'replay.ex') if a==str(ROOT/'launch.ex') else a for a in cmd]
    (ROOT/'replay.ex-build-command.json').write_text(json.dumps(cmd,indent=2)+'\n');subprocess.run(cmd,check=True)
def smoke():
    original=(HERE/'t18c-full-sqrt.txt').read_text()
    boxfile=ROOT/'smoke-boxes.txt'
    boxfile.write_text('1 32 0 95 31\n2 96 0 159 31\n')
    for variant in ('ordinary','compact','geometric'):
        d=ROOT/'smoke'/variant;d.mkdir(parents=True,exist_ok=True)
        p=d/'params.txt';p.write_text(params.replace(original,dict(N1=64,N2=32,L=112,center='56 0',star_centre='56 0',mass_extraction_center='56 0',
            max_level=2,t18c_boxes=boxfile,regrid_interval='0 0',t18c_compact_initial_audit=str(variant!='ordinary').lower(),
            ems_use_geometric_initial_lapse=str(variant=='geometric').lower(),
            t18c_compare_nonlapse=ROOT/'smoke/compact/t18c-nonlapse.bin' if variant=='geometric' else '""')))
        measured('smoke-'+variant,d,[ROOT/'launch.ex',p])
    for file in ('t18c-nonlapse.bin','t18c-t0-native.csv','t18c-ranges.csv','t18c-source.csv'):
        assert (ROOT/'smoke/ordinary'/file).read_bytes()==(ROOT/'smoke/compact'/file).read_bytes(),file
    assert (ROOT/'smoke/compact/t18c-t0-native.csv').read_bytes()==(ROOT/'smoke/geometric/t18c-t0-native.csv').read_bytes()
    (HERE/'t18c-software-check.json').write_text(json.dumps(dict(compact_initialization='BIT_IDENTICAL',nonlapse='BIT_IDENTICAL',native_constraints='BIT_IDENTICAL',scope='L0-2 software fixture including point transfers; physical results remain pending'),indent=2)+'\n')
if __name__=='__main__':
    {'generate':generate,'census':build_census,'prepare':prepare,'build':build_native,'smoke':smoke}[sys.argv[1]]()
