#!/usr/bin/env python3
"""Register/build the documented-diagnostic binary and its serial pipeline."""
import sys
sys.dont_write_bytecode=True
import csv,hashlib,importlib.util,json,math,shutil,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];ROOT=Path('/private/tmp/ems-t16b')
PY='/Users/auroradysis/miniconda3/bin/python'
DATA=Path('/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-echobin/.data/echo-binary/production-outer16-a20-r32-a10/n32-r6.ctt')
DOC=Path('/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-echobin/artifacts/echo-binary/production/outer16-a20/r32-a10')
PROFILE=Path('/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-echobin/artifacts/echo-binary/a0.9-e8.trumpet')
SHA='46253cb8d3bf6255b149b44bed7d4899cd97c697fe677a9846649bfbb2bdc4ee'
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
p=module('prepare',HERE/'t13-prepare.py')
def save(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def measured(label,directory,cmd):
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/(label+'.log')).open('wb') as log:
        subprocess.run([PY,str(HERE/'t16b-run.py'),'--measure',label,'--directory',str(directory),'--',*map(str,cmd)],stdout=log,stderr=subprocess.STDOUT,check=True)
def register():
    ROOT.mkdir(exist_ok=True);assert hashlib.sha256(DATA.read_bytes()).hexdigest()==SHA
    inputs=[DATA,PROFILE,DOC/'params-initial-data.txt',DOC/'audit.md']
    (HERE/'t16b-inputs.json').write_text(json.dumps([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),bytes=f.stat().st_size) for f in inputs],indent=2)+'\n')
    panels=[]
    with DATA.open() as f:
        for line in f:
            if line.startswith('PANEL '):
                q=line.split();panels.append(dict(panel=int(q[1]),kind=q[2],lo=float(q[3]),hi=float(q[4]),center=float(q[5]),stretch=float(q[6]),mu_lo=float(q[7]),mu_hi=float(q[8]),boundary=q[9],nr=int(q[10]),na=int(q[11])))
    assert len(panels)==34;save('t16b-panels.csv',panels)
    base=(HERE/'t16-L13-geometric.txt').read_text();audit_params=p.parameters((DOC/'params-initial-data.txt').read_text())
    jobs=[];registration=[]
    for rung,n in [('coarse',672),('fine',1344)]:
        h0=672/n;h=h0/8192;dt=h/4;steps=int(.002/dt)
        for mode in ('sqrt','geometric'):
            name=rung+'-'+mode;d=ROOT/'evolution'/name
            for sub in ('plt','chk'):(d/sub).mkdir(parents=True,exist_ok=True)
            changes=dict(audit_params);changes.update(N1=n,N2=n//2,L=672,center='336 0',star_centre='336 0',
                ems_use_geometric_initial_lapse=str(mode=='geometric').lower(),ems_use_maximal_initial_lapse='false',
                t16b_audit_only='false',t16b_compare_nonlapse=ROOT/'evolution'/(rung+'-sqrt')/'t16b-nonlapse.bin' if mode=='geometric' else '""',
                t14_dense_initial_tags='false',activate_mass_extraction=0,t13_launch_stop_time=.002,max_level=13,
                RH_activate='false',RH_num_horizons=0,verbosity=0,regrid_interval=' '.join(['0']*13))
            param=HERE/f't16b-{name}.txt';param.write_text(p.replace(base,changes))
            jobs.append(dict(name=name,directory=str(d),command=[PY,str(HERE/'t16b-work.py'),'leg',name]))
            registration.append(dict(run=name,base_N1=n,base_N2=n//2,h0_M=h0,finest_level=13,h_M=h,Rh_over_h=.0060404520035922523/h,
                face_each_hole_M=128/8192,dt_M=dt,native_steps=steps,native_stop_M=steps*dt,common_clocks=66,
                common_stop_M=65/8192/4,cap_RSS_bytes=6000000000,parameters=str(param)))
    save('t16b-registration.csv',registration)
    (ROOT/'evolution/plan.json').write_text(json.dumps(dict(jobs=jobs),indent=2)+'\n')
    d=ROOT/'pipeline/job';d.mkdir(parents=True,exist_ok=True)
    (ROOT/'pipeline/plan.json').write_text(json.dumps(dict(jobs=[dict(name='binary-pipeline',directory=str(d),command=[PY,str(HERE/'t16b-work.py'),'pipeline'])]),indent=2)+'\n')
    note='''## T16b — documented-diagnostic binary companion, preregistered native audit and launch

Baseline is fork main 31176be. The operator authorizes n32-r6.ctt (SHA-256 46253cb8d3bf6255b149b44bed7d4899cd97c697fe677a9846649bfbb2bdc4ee) as DOCUMENTED DIAGNOSTIC DATA for this binary test, not certified data. Its audit reports only value/first-derivative interface matching, a D2(logpsi) jump about 8.89e-6, and crossing/within ratios up to 306. No companion is constructed or altered. The exact supplied mass, separation, rapidity and coupling lines are consumed; the centres are translated together to (320,0) and (352,0) in the native [0,672]x[0,336] domain.

Before compilation or evolution, the two GLOBAL rungs are fixed: base 672x336 with h0=1 M_i and 1344x672 with h0=0.5 M_i. Both have levels 0–13 and per-hole square faces 128/2^l M_i at l>=1; overlapping outer boxes merge. Finest h=1/8192 and 1/16384, R_h/h=49.4834 and 98.9668 at each hole. Dense Tests-only fixed tags select those boxes; production equations, gauge, KO=1, Float64, puncture treatment, transfers and the unmodified RHFinder remain unchanged. Regridding and RHFinder are disabled and the t2 reader guard is on. Global refinement is necessary to change the native interface-sheet spacing; adding only a puncture level would leave the exterior test unchanged.

At numerical t=0, the final native CTT state is audited before any advance. Every actual cell including ghosts is checked for each isolated geometric a_A and combined alpha_bin: finite, 0<a<=1, without clipping. The stable existing setter formula is used. The geometric state is compared directly, bit by bit, to the default's stored non-lapse stream on the same hierarchy before it advances. Native physical and KO Gamma RHS, its existing 13-term attribution, K/lapse RHS and initialized-driver cancellation are recorded in the four innermost cells at each hole, with all 28 RHS components checked against Base. There is no binary translation or Killing-source identity. Puncture-limit samples use R=10^-2,...,10^-10 and report alpha_bin/a_own for each hole, excluding the puncture itself.

Composite native constraint regions use uncovered valid cells and cylindrical-coordinate-volume weights y*h^2. Fixed regions are each Euclidean hole collar 0.01<=r_A<=0.02, the inter-hole box |x-336|<=8, 0<=y<=8, far field 80<=r_global<=128, and elsewhere with r_near>=0.01 and r_global<=128 excluding interface sheets. Ham, vector Mom, GaussE, GaussB and vector C_Gamma are reported in peak/RMS. No continuum solver residual is substituted for native derivatives.

Interface sheets are defined geometrically from all 34 companion PANEL records and the native selector, not from measured constraints. At each physical point, the declared COARSE hierarchy fixes h_c(x,y). Every distinct panel selected at offsets (i h_c,j h_c), i,j=-3,...,3, labels the associated sheet with the ordered panel-ID pair. This defines the same physical sheet mask on both rungs and covers the native FD/KO stencil. The selector is cross-checked against native forced-panel evaluation on sampled actual cells. Actual native-stencil crossings are recorded separately, as are coarse-AMR-face proximity flags. Each pair is reported both in the full sheet and with coarse AMR faces excluded. Interfaces that no actual native stencil crosses are listed as unobserved. Fine/coarse constraint norm ratios report growth, decrease or sampling absence; absence of growth on two coarse rungs does not certify C2 matching or admit a convergence campaign despite the known interface failure.

The launch uses dt=h/4 and the isolated T16 W=[0.00075,0.0025] M_i at both holes. Axis/diagonal rays point toward the other hole: n=(+1,0),(+1,+1)/sqrt(2) on the left and (-1,0),(-1,+1)/sqrt(2) on the right. The numerical t=0 profile is subtracted without alignment. Central I8 and |I8-I10| spreads, peak/trapezoidal RMS, and the 5x significance rule are frozen. All 66 common clocks (t=0 plus 65 positive clocks) through 0.001983642578125 M_i are retained; native fine stop is 0.0019989013671875 M_i and its unmatched last half-clock is excluded. Ratios are geometric/sqrt(chi), with five-spread sensitivity intervals and early enhancement entries retained. Native puncture cells and floors corroborate the profiles. Both h and dt change, so no same-grid temporal bound is claimed. Copy-free extra recorder calls only write current native arrays; the first four steps also retain pre-RHS state, enforced state and resulting RHS at both holes. Snapshot compression is lossless and crops only stored output to the same I10/native-stencil support.

Numerical work is serial with OMP=2 and BLAS=1, <=4 active threads. The per-process RSS/footprint watchdog is 6,000,000,000 bytes, tree gate 6,500,000,000 bytes, output ceiling 8 GB. The potentially long four-leg queue and analysis run detached, with own done.exit and actual-child resource receipts. No SSH or commit. A low-resolution merger recommendation will distinguish an initial-lapse choice from the companion's unresolved convergence admission.
'''
    path=HERE/'README.md';text=path.read_text();assert '## T16b —' not in text
    path.write_text(text.rstrip()+'\n\n'+note)
    print('Registered diagnostic binary, global two-rung chain and fixed masks before runs.',flush=True)
def build():
    subprocess.run(['make','-j1','all'],cwd=REPO/'Examples/EMS',check=True)
    text=(REPO/'Examples/EMS/Main_EMSBH2DBH.cpp').read_text().replace('DefaultLevelFactory<EMSBH2DLevel>','T16BFactory<EMSBH2DLevel>')
    needle='    using Clock = std::chrono::steady_clock;'
    assert text.count(needle)==1
    text=text.replace(needle,'    t16b_initial_audit(bh_amr,sim_params);\n    bool audit_only=false;pp.load("t16b_audit_only",audit_only,false);\n    if(audit_only)return 0;\n    pout()<<"T16b evolution starts after completed t=0 audit"<<std::endl;\n'+needle)
    (HERE/'t16b-main.hpp').write_text(text)
    b=module('builder',HERE/'t14-prepare.py');b.ROOT=ROOT;b.build('t16b-launch.cpp','launch.ex')
    shutil.copyfile('/private/tmp/ems-t16/replay.ex',ROOT/'replay.ex');(ROOT/'replay.ex').chmod(0o755)
    # Existing lossless filter, with the binary physical centres and common clock.
    text=(HERE/'t16-xz.py').read_text().replace('clock=.875/8192/4','clock=1./8192/4')
    text=text.replace('(abs((iv[:,0]+.5)*h-336)<=.0035)',
        '((abs((iv[:,0]+.5)*h-320)<=.0035)|(abs((iv[:,0]+.5)*h-352)<=.0035))')
    (HERE/'t16b-xz.py').write_text(text)
    d=ROOT/'bin';d.mkdir(exist_ok=True);(d/'xz').write_text(text);(d/'xz').chmod(0o755)
def leg(name):
    import os
    os.environ['PATH']=str(ROOT/'bin')+':'+os.environ['PATH']
    os.execv(str(ROOT/'launch.ex'),[str(ROOT/'launch.ex'),str(HERE/f't16b-{name}.txt')])
def pipeline():
    measured('evolution-worker',ROOT/'evolution',[PY,HERE/'t16b-run.py','--worker',ROOT/'evolution/plan.json'])
    measured('binary-analysis',ROOT/'analysis',[PY,HERE/'t16b-analyze.py'])
    measured('binary-manifest',ROOT/'manifest',[PY,HERE/'t16b-analyze.py','manifest'])
    print('T16b diagnostic binary pipeline finished.',flush=True)
def smoke():
    base=(HERE/'t16b-coarse-sqrt.txt').read_text()
    for mode in ('sqrt','geometric'):
        d=ROOT/'smoke'/mode
        for sub in ('plt','chk'):(d/sub).mkdir(parents=True,exist_ok=True)
        changes=dict(N1=336,N2=168,max_level=4,regrid_interval='0 0 0 0',t16b_audit_only='true',
            ems_use_geometric_initial_lapse=str(mode=='geometric').lower())
        if mode=='geometric':changes['t16b_compare_nonlapse']=ROOT/'smoke/sqrt/t16b-nonlapse.bin'
        param=d/'params.txt';param.write_text(p.replace(base,changes))
        measured('smoke-'+mode,d,[ROOT/'launch.ex',param])
    check_smoke()
def check_smoke():
    for mode in ('sqrt','geometric'):
        d=ROOT/'smoke'/mode
        r=json.loads((d/('smoke-'+mode+'.resources.json')).read_text())
        assert r['returncode']==r['child_measurement']['returncode']==0 and not r['gate_reason'] and r['peak_rss_bytes']<6e9
        ranges=list(csv.DictReader((d/'t16b-ranges.csv').open()))
        for key in ('a_left','a_right','alpha_bin','native_lapse'):
            assert 0<min(float(x[key+'_min']) for x in ranges)<=max(float(x[key+'_max']) for x in ranges)<=1
        source=list(csv.DictReader((d/'t16b-source.csv').open()))
        assert len(source)==8 and all(int(x['Base_bit_mismatches'])==0 and float(x['driver_cancel'])==0 for x in source)
    assert (ROOT/'smoke/sqrt/t16b-t0-native.csv').read_bytes()==(ROOT/'smoke/geometric/t16b-t0-native.csv').read_bytes()
    save('t16b-software-check.csv',[dict(check='coarse binary t=0 software fixture; not a physical result',result='PASS',
        default_actual_child=0,geometric_actual_child=0,nonlapse_bit_mismatches=0)])
def filter_check():
    import struct
    import numpy as np
    decoder=module('decoder',HERE/'t7-check.py');h=1/16384;clock=1/8192/4
    cells=np.array([[round(c/h)+i,0] for c in (320,352) for i in (0,1,64)]+[[round(336/h),0]],dtype='<i4')
    frames=[];previous=None;raw=bytearray(b'T7OP0002')
    for index,(phase,time) in enumerate([(50,0.),(3,clock/2),(50,clock/2),(50,clock),(50,2*clock)]):
        values=(np.arange(len(cells)*28).reshape(len(cells),28)+index/8).astype('<f8')
        payload=values.view('u1').reshape(-1,8).T.reshape(-1).copy();encoded=payload.copy()
        if previous is not None:encoded^=previous
        raw.extend(struct.pack('<11i',phase,13,0,index,28,len(cells),int(previous is not None),0,0,1,1))
        raw.extend(struct.pack('<10d',time,h,clock/2,0,0,0,0,0,0,0));raw.extend(cells.tobytes());raw.extend(encoded.tobytes())
        previous=payload;frames.append((phase,time,values))
    d=ROOT/'software';d.mkdir(exist_ok=True);path=d/'filter.xz'
    result=subprocess.run([PY,HERE/'t16b-xz.py'],input=raw,capture_output=True,check=True);path.write_bytes(result.stdout)
    got=list(decoder.frames(path));expected=[x for x in frames if x[0]!=50 or x[1]!=clock/2]
    assert len(got)==len(expected)
    for actual,(phase,time,values) in zip(got,expected):
        keep=np.array([0,1,3,4]) if phase==50 else np.arange(len(cells))
        assert actual['phase']==phase and actual['meta'][0]==time
        assert actual['cells'].tobytes()==cells[keep].tobytes() and actual['values'].tobytes()==values[keep].tobytes()
    save('t16b-filter-check.csv',[dict(check='both-hole crop, half-clock exclusion, first-stage retention, XOR round trip',result='BIT_IDENTICAL',input_frames=5,output_frames=4)])
if __name__=='__main__':
    if sys.argv[1]=='leg':leg(sys.argv[2])
    else:{'register':register,'build':build,'pipeline':pipeline,'smoke':smoke,'check-smoke':check_smoke,'filter-check':filter_check}[sys.argv[1]]()
