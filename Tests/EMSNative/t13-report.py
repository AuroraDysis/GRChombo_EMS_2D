#!/usr/bin/env python3
"""T13 measured-stage budgets, figures and conditional (unrun) source ladder."""
import sys
sys.dont_write_bytecode=True
import csv, importlib.util, json, math, os
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent; TMP=Path('/private/tmp/ems-t13'); OUT=TMP/'analysis'
spec=importlib.util.spec_from_file_location('analysis',HERE/'t13-analyze.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
GROUPS=('advection','shift_laplacian_2D','shift_grad_div_2D','lapse_gradient_A',
        'chi_gradient_A','K_gradient','Theta_gradient','geometry_2D',
        'cartoon_geometry','cartoon_shift','reduction_constraint','matter','KO')
VARS='chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi'.split()

def read(name):
    with (HERE/name).open() as f:return list(csv.DictReader(f))

def attribution():
    rows=[];budgets=[];projections=[];closure=[]
    # The same fixed P6/P8 spatial operator is applied to each signed RHS group.
    # Groups are summed as profiles before taking norms (norms are not additive).
    a.SELECT=list(range(a.COL));a.ODD_C=a.ODD
    for leg in a.LEGS:
        dt=a.DT/(2 if leg.endswith('half') else 1);groups={};integrals={}
        for f in sorted((OUT/(leg+'-stages')).glob('frame-*.npz')):
            with np.load(f) as z:
                key=(float(z['meta'][0]),int(z['stage']))
                groups.setdefault(key,[]).append((z['cells'].copy(),z['q'].copy(),z['meta'].copy()))
        assert len(groups)==16,(leg,len(groups))
        for (start,stage),frames in sorted(groups.items()):
            cells=np.concatenate([z[0] for z in frames]);q=np.concatenate([z[1] for z in frames]);meta=frames[0][2]
            assert len(set(map(tuple,cells)))==len(cells)
            step=round(start/dt)+1
            for ray,vec in a.NV.items():
                p6=a.sample(cells,q,vec,6);p8=a.sample(cells,q,vec,8)
                terms6=p6[:,112:138].reshape(-1,13,2)@vec
                terms8=p8[:,112:138].reshape(-1,13,2)@vec
                total=p6[:,163:165]@vec;adv=p6[:,95:97]@vec;ko=p6[:,67:69]@vec
                physical=p6[:,39:41]@vec;geometric=physical-adv
                profiles=dict(geometric=geometric,advection=adv,KO=ko,total=total)
                weight=(1/6,1/3,1/3,1/6)[stage]*dt
                for group,v in profiles.items():
                    key=(leg,step,ray,group)
                    integrals[key]=integrals.get(key,np.zeros(len(a.R)))+weight*v
                for i,group in enumerate(GROUPS):
                    v=terms6[:,i];spread=terms6[:,i]-terms8[:,i]
                    rows.append(dict(run=leg,step=step,stage=stage,stage_time_M=float(meta[3]),
                        ray=ray,region='fixed W P6',group=group,peak_RHS=float(abs(v).max()),
                        RMS_RHS=a.rms(v),probe_RHS=float(np.interp(.0015,a.R,v)),
                        interpolation_peak=float(abs(spread).max()),interpolation_RMS=a.rms(spread)))
                for group,v in profiles.items():
                    rows.append(dict(run=leg,step=step,stage=stage,stage_time_M=float(meta[3]),
                        ray=ray,region='fixed W P6',group=group,peak_RHS=float(abs(v).max()),
                        RMS_RHS=a.rms(v),probe_RHS=float(np.interp(.0015,a.R,v)),
                        interpolation_peak='',interpolation_RMS=''))
            # Native innermost diagonal cell: no ray interpolation or background.
            i=np.flatnonzero((cells[:,0]==a.CENTRE)&(cells[:,1]==0));assert len(i)==1
            v=q[i[0]];vec=a.NV['diagonal'];terms=v[112:138].reshape(13,2)@vec
            for group,value in zip(GROUPS,terms):
                rows.append(dict(run=leg,step=step,stage=stage,stage_time_M=float(meta[3]),
                    ray='diagonal',region='first puncture cell',group=group,peak_RHS=abs(value),
                    RMS_RHS=abs(value),probe_RHS=value,interpolation_peak=0.,interpolation_RMS=0.))
            for var,index in [('lapse',13),('shift_n',None),('B_n',None),('A11',6),('A12',7),
                              ('A22',8),('Aww',9),('phi',18),('Pi',19),('Ex',24),('Ey',25),('Xi',27)]:
                take=lambda offset: float(v[offset+index]) if index is not None else float(v[offset+(14 if var=='shift_n' else 16):offset+(16 if var=='shift_n' else 18)]@vec)
                for group,offset in [('total',152),('physical',28),('native_advection',84),('KO',56)]:
                    value=take(offset)
                    rows.append(dict(run=leg,step=step,stage=stage,stage_time_M=float(meta[3]),
                        ray='diagonal',region='first puncture cell '+var,group=group,peak_RHS=abs(value),
                        RMS_RHS=abs(value),probe_RHS=value,interpolation_peak=0.,interpolation_RMS=0.))
        for ray in a.NV:
            z=np.load(OUT/(leg+'-'+ray+'-profiles.npz'));v=z['p6'][:,:,11:13]@a.NV[ray]
            for step in range(1,5):
                actual=v[step]-v[step-1]
                summed=integrals[leg,step,ray,'total']
                defect=float(abs(actual-summed).max())
                closure.append(dict(run=leg,ray=ray,step=step,
                    peak_update=float(abs(actual).max()),weighted_RHS_peak=float(abs(summed).max()),
                    absolute_defect=defect,condition='fixed P6 on native states; interpolation evaluation roundoff'))
                assert defect<1e-13,(leg,ray,step,defect)
                for group in ('geometric','advection','KO','total'):
                    value=integrals[leg,step,ray,group]
                    budgets.append(dict(run=leg,step=step,ray=ray,group=group,
                        peak_step_delta=float(abs(value).max()),RMS_step_delta=a.rms(value),
                        signed_probe_step_delta=float(np.interp(.0015,a.R,value)),
                        KO_over_total_peak=float(abs(integrals[leg,step,ray,'KO']).max()/abs(summed).max()),
                        KO_over_total_RMS=a.rms(integrals[leg,step,ray,'KO'])/a.rms(summed)))
        pairs={14:(2,'input_trace'),15:(14,'input_floors'),5:(4,'RK_trace'),
               16:(6,'full_step_trace'),17:(16,'full_step_floors')};before={};summary={}
        for frame in a.decoder.frames(TMP/'evolution'/leg/'t13-t7-stage-L12.xz'):
            p=frame['phase'];src=frame['source']
            if p in (2,4,6,14,16):before[src,p]=frame
            if p not in pairs:continue
            oldp,op=pairs[p];old=before[src,oldp]
            assert np.array_equal(old['cells'],frame['cells'])
            d=frame['values']-old['values'];mx=abs(d).max(0)
            current=summary.get(op,np.zeros(28));summary[op]=np.maximum(current,mx)
            assert not np.count_nonzero(d[:,[0,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]]), (leg,op)
        for op,mx in summary.items():
            for var,value in zip(VARS,mx):
                projections.append(dict(run=leg,operation=op,component=var,
                    max_abs_update=float(value),condition='all recorded native cells incl. ghost collar; first four steps'))
    a.save('t13-stage-attribution.csv',rows);a.save('t13-stage-budgets.csv',budgets)
    a.save('t13-projections.csv',projections);a.save('t13-profile-update-closure.csv',closure)
    print('Attribution:',len(rows),'rows; max profile-update defect',max(z['absolute_defect'] for z in closure),flush=True)

def stage_c():
    # Reuse exact parameter editing/census, no initialization or evolution run.
    s=importlib.util.spec_from_file_location('prepare',HERE/'t13-prepare.py');p=importlib.util.module_from_spec(s);s.loader.exec_module(p)
    screen=read('t13-screen.csv');cs=read('t13-constraint-screen.csv');audit=read('t13-run-audit.csv')
    assert all(float(z['worst_nonzero_time_ratio'])<=.5 and float(z['reduction_over_interpolation'])>5 and float(z['reduction_over_temporal'])>5
               and float(z['minimum_history_reduction_over_interpolation'])>5 and float(z['minimum_history_reduction_over_temporal'])>5 for z in screen)
    assert all(int(z['times_significantly_worse'])==0 for z in cs)
    assert all(int(z['chi_activations'])==int(z['lapse_activations'])==int(z['nonfinite_count'])==0 for z in audit)
    source=(HERE/'params/t13-E-mid-maximal.txt').read_text();base=p.parameters(source);radii=base['mass_extraction_radii'].split()
    levels=base['mass_extraction_levels'].split();census=read('t13-box-census.csv');oldcounts=[int(z['valid_cells']) for z in census]
    cellsteps0=sum(c*max(1,math.ceil(56/2**(12-l))) for l,c in enumerate(oldcounts))
    measured=float(next(z for z in audit if z['run']=='maximal')['wall_s'])
    rss=float(next(z for z in audit if z['run']=='maximal')['peak_rss_GB']);extra=int(census[-1]['estimated_array_bytes'])
    steps=[];faces=[];storage=[]
    for finest in (13,14):
        rr=radii+[format(float(radii[-1])/2**j,'.17g') for j in range(1,finest-11)]
        ll=levels+[str(j) for j in range(13,finest+1)]
        name=f't13-stageC-E-mid-max{finest}.txt'
        text=p.replace(source,dict(max_level=finest,regrid_interval=' '.join(['0']*finest),
             num_mass_extraction_radii=finest,mass_extraction_radii=' '.join(rr),mass_extraction_levels=' '.join(ll)))
        (HERE/'params'/name).write_text(text)
        parsed=p.parameters(text)
        assert parsed['mass_extraction_radii'].split()[:12]==radii
        assert parsed['mass_extraction_levels'].split()[:12]==levels
        for key in base:
            if key not in ('max_level','regrid_interval','num_mass_extraction_radii','mass_extraction_radii','mass_extraction_levels'):
                assert parsed[key]==base[key],key
        counts=oldcounts+[oldcounts[-1]]*(finest-12);nt=56*2**(finest-12)
        cost=sum(c*max(1,math.ceil(nt/2**(finest-l))) for l,c in enumerate(counts))
        # Preserve measured wall including startup; the ratio is conservative
        # because startup does not scale with the number of native cell steps.
        planned=measured*cost/cellsteps0
        peak=rss+(finest-12)*extra*1.35/1e9+.35 # recorder/copy and allocator headroom
        assert peak<6
        steps.append(dict(finest_level=finest,h_M=(7/12)/2**finest,finest_steps=nt,
            dt0_M=7/48,dt_finest_M=(7/48)/2**finest,endpoint_M=56*a.DT,
            cell_steps=cost,cell_step_ratio=cost/cellsteps0,wall_estimate_s=planned,
            planned_wall_range_s=f'{planned*.7:.0f}–{planned*1.5:.0f}',estimated_peak_RSS_GB=peak,
            condition='NOT RUN; baseline startup retained; eight aligned boxes per added level predicted, actual census mandatory'))
        for lev in range(finest+1):
            h=(7/12)/2**lev;face=336 if lev==0 else 112/2**lev
            faces.append(dict(ladder=finest,level=lev,h_M=h,face_M=face,
                face_cells=round(face/h),old_face_unchanged=lev<=12,
                condition='old actual census retained; added unions predicted x±R,y[0,R]; 192 fine cells to face'))
        # Output estimate from the actual nominal recorder lengths, linear in
        # recorded ROI area and number of snapshots; compression ratio retained.
        d=TMP/'evolution/maximal';old_bytes=sum(f.stat().st_size for f in d.glob('*.xz'))
        roi0=math.ceil(.006/a.H);roi=math.ceil(.006/((7/12)/2**finest))
        estimate=old_bytes*(roi+3)**2/(roi0+3)**2*(nt+1)/57
        storage.append(dict(finest_level=finest,existing_compressed_bytes=old_bytes,
            ROI_cells_per_side=roi,estimated_compressed_bytes=round(estimate),
            gate_bytes=5_700_000_000,condition='upper scaling estimate includes fixed first-stage capture; actual xz ratio may differ'))
    a.save('t13-stageC-plan.csv',steps);a.save('t13-stageC-faces.csv',faces);a.save('t13-stageC-storage.csv',storage)
    (HERE/'params/t13-stageC-plan.json').write_text(json.dumps({'jobs':[
        dict(name=f'max{level}',directory=f'/private/tmp/ems-t13-stageC/evolution/max{level}',
             command=[str(TMP/'launch.ex'),str(HERE/'params'/f't13-stageC-E-mid-max{level}.txt')])
        for level in (13,14)]},indent=2)+'\n')
    (HERE/'params/t13-stageC-reading.md').write_text('''Prepared only; no Stage C evolution has been run.

Run independent max13 and max14 alpha_K initializations, alongside the completed max12 branch. Keep dt0=7/48 M, all 12 existing tags/faces, sigma=1, cleaners, transfers, guard, recorder and W unchanged. The native endpoint is 0.001993815104166667 M in all three runs. Actual t=0 box unions, positive lapse/floor margins and incoming speed/support bounds must be checked before advancing; stop if any old face changes. Added face coordinates are 0.013671875 and 0.0068359375 M. Added square unions and their eight-box decomposition are block-aligned predictions, not a measured new census.

Compare raw signed Gamma, native metric Gamma and longitudinal shift profiles on the same fixed physical W and common max12 sample times, on axis and diagonal, without phase alignment or background fitting. Use P6 and P8 on current evolved fields. Primary contraction norms: fixed-window peak and unweighted radial RMS, both endpoint and history. rho_D=norm(u_h/2-u_h/4)/norm(u_h-u_h/2) <=0.8; each difference must exceed 5 times its joint P6/P8 spread. Temporal error must be subdominant; carry the Stage B temporal spread conservatively and require a bounded dt/2 check at the refined rung if it could decide the screen. Do not impose amplitude stability. Report differences obscured by uncertainty rather than dividing by numerical zero.

Significant contraction promotes to the existing exterior/near-hole pilot; 0.8<rho_D<1 or uncertainty-dominated differences are inconclusive; significant growth/noncontraction with rho_D>=1 kills lapse alone as sufficient. Fresh positivity/floor screens and the absolute constraint battery remain required. The short window does not rule out delayed emission or admit 100 M production. A separate dissipation card is needed if KO regeneration prevents contraction; this card does not alter KO, projections or gauge.

Local estimate: 80/156 s for max13/max14 (ranges 56-120/109-233 s), 2.15/2.33 GB peak RSS; estimates include startup and headroom, not measurements. Estimated new compressed outputs 0.57/4.35 GB; use a fresh, separately gated output root capped at 5.7 GB, retain the completed Stage B evidence, and stop if actual storage exceeds the gate. Set T13_OUTPUT_ROOT=/private/tmp/ems-t13-stageC and copy t13-stageC-plan.json into its evolution directory before using the existing detached runner. Do not reuse the completed Stage B directories or markers. This plan is prepared, not launched.
''')
    print('Prepared, NOT RUN:',steps,storage,flush=True)

def predictor():
    """One recorded predictor dependency test, NOT an evolution or initialization."""
    rows=[];cells=a.selected();a.SELECT=list(range(a.COL));a.ODD_C=a.ODD
    for leg in a.LEGS:
        dt=a.DT/(2 if leg.endswith('half') else 1);cached={0:[],1:[]}
        for f in sorted((OUT/(leg+'-stages')).glob('frame-*.npz')):
            with np.load(f) as z:
                if float(z['meta'][0])==0 and int(z['stage']) in cached:
                    cached[int(z['stage'])].append((z['cells'].copy(),z['q'].copy()))
        iv0=np.concatenate([x[0] for x in cached[0]]);q0=np.concatenate([x[1] for x in cached[0]])
        lookup={tuple(k):v for k,v in zip(iv0,q0)}
        iv1=np.concatenate([x[0] for x in cached[1]]);q1=np.concatenate([x[1] for x in cached[1]])
        for variant in ('recorded_predictor','without_all_KO','without_A_KO','without_EMS_KO'):
            record=np.empty((len(cells),3+49*28));record[:,0]=a.H
            record[:,1]=(cells[:,1]+.5)*a.H;record[:,2]=(cells[:,0]+.5)*a.H-336
            for j in range(-3,4):
                for i in range(-3,4):
                    values=[]
                    for x,y in cells+np.array([i,j]):
                        yy=-y-1 if y<0 else y;source=lookup[x,yy]
                        state=source[:28]+dt/2*source[152:180]
                        if variant!='recorded_predictor':
                            keep=list(range(28)) if variant=='without_all_KO' else list(range(6,10)) if variant=='without_A_KO' else list(range(18,28))
                            state[keep]-=dt/2*source[56:84][keep]
                        if y<0:state[a.ODD_FIELDS]*=-1
                        values.append(state)
                    record[:,3+(j+3)*7+i+3::49]=np.array(values)
            payload=OUT/(leg+'-'+variant+'-input.bin');output=OUT/(leg+'-'+variant+'-q.bin');record.tofile(payload)
            a.timed('predictor-'+leg+'-'+variant,[str(TMP/'replay.ex'),'--stage',str(payload),str(output)])
            qq=np.fromfile(output).reshape(-1,a.COL)
            for ray,vec in a.NV.items():
                rhs=a.sample(cells,qq,vec,6)[:,39:41]@vec
                actual=a.sample(iv1,q1,vec,6)[:,39:41]@vec
                mismatch=float(abs(rhs-actual).max())
                if variant=='recorded_predictor':assert mismatch<1e-8,(leg,ray,mismatch)
                rows.append(dict(run=leg,variant=variant,ray=ray,region='W P6',
                    physical_Gamma_peak=float(abs(rhs).max()),physical_Gamma_RMS=a.rms(rhs),
                    change_from_recorded_peak=mismatch,change_from_recorded_RMS=a.rms(rhs-actual),
                    condition='first dt/2 predictor; only indicated recorded KO increment subtracted; native projections and RHS reapplied'))
            i=np.flatnonzero((cells[:,0]==a.CENTRE)&(cells[:,1]==0));j=np.flatnonzero((iv1[:,0]==a.CENTRE)&(iv1[:,1]==0));assert len(i)==len(j)==1
            rhs=float(qq[i[0],39:41]@a.NV['diagonal']);actual=float(q1[j[0],39:41]@a.NV['diagonal'])
            rows.append(dict(run=leg,variant=variant,ray='diagonal',region='first puncture cell',
                physical_Gamma_peak=abs(rhs),physical_Gamma_RMS=abs(rhs),
                change_from_recorded_peak=abs(rhs-actual),change_from_recorded_RMS=abs(rhs-actual),
                condition=f'signed RHS {rhs:.17g}; first predictor only; not outgoing amplitude attribution'))
            (OUT/(leg+'-'+variant+'-input.sha256')).write_text(a.hashlib.sha256(record.tobytes()).hexdigest()+' '+str(payload)+'\n')
            payload.unlink() # Regenerable dependency-test payload; hash and native output retained.
    a.save('t13-predictor-feedback.csv',rows)
    print('First-predictor dependency test:',json.dumps(rows,indent=2),flush=True)

def input_hashes():
    # Recover the exact earlier regenerable payload hashes from retained frames;
    # no RHS rerun, no static input, and one bounded frame at a time.
    count=0
    for leg in a.LEGS:
        for f in sorted((OUT/(leg+'-stages')).glob('frame-*.npz')):
            with np.load(f) as z:
                cells=z['cells'];pre={tuple(k):v for k,v in zip(z['pre_cells'],z['pre_values'])}
                h=z['meta'][1];record=np.empty((len(cells),3+49*28))
                record[:,0]=h;record[:,1]=(cells[:,1]+.5)*h;record[:,2]=0.
                for j in range(-3,4):
                    for i in range(-3,4):
                        record[:,3+(j+3)*7+i+3::49]=np.array([pre[k[0]+i,k[1]+j] for k in cells])
                f.with_name(f.stem+'-input.sha256').write_text(a.hashlib.sha256(record.tobytes()).hexdigest()+'  reconstructed native replay payload\n')
                count+=1
    assert count==128,count
    print('Exact replay payload hashes recovered from retained frames:',count,flush=True)

def verify():
    rows=[]
    def check(name,value,condition):
        assert condition,(name,value);rows.append(dict(check=name,value=value,result='PASS'))
    controls=read('t13-controls.csv');identity=read('t13-initial-identity.csv');stages=read('t13-stage-closure.csv')
    check('frozen default/control files',len(controls),len(controls)==12 and all(int(x['bit_mismatches'])==0 for x in controls))
    check('non-lapse initialized mismatches',sum(int(x['bit_mismatches']) for x in identity),all(int(x['bit_mismatches'])==0 for x in identity))
    check('actual production RHS mismatches',sum(int(x['RHS_bit_mismatches']) for x in stages),len(stages)==4 and all(int(x['RHS_bit_mismatches'])==int(x['projected_state_bit_mismatches'])==0 for x in stages))
    for leg in a.LEGS:
        d=TMP/'evolution'/leg
        check(leg+' done.exit',d.joinpath('done.exit').read_text().strip(),d.joinpath('done.exit').read_text().strip()=='0')
    audit=read('t13-run-audit.csv')
    check('floor/nonfinite/reader counts',0,all(int(x[k])==0 for x in audit for k in ('chi_activations','lapse_activations','nonfinite_count','reader_messages_after_first_advance')))
    check('common native endpoints',float(audit[0]['end_time_M']),max(float(x['end_time_M']) for x in audit)-min(float(x['end_time_M']) for x in audit)<1e-18)
    screen=read('t13-screen.csv')
    check('largest history disturbance ratio',max(float(x['worst_nonzero_time_ratio']) for x in screen),all(float(x['worst_nonzero_time_ratio'])<.5 for x in screen))
    check('resolved absolute constraint worsening',0,all(int(x['times_significantly_worse'])==0 for x in read('t13-constraint-screen.csv')))
    check('all exact replay payload hashes',128,len(list(OUT.glob('*-stages/*-input.sha256')))==128)
    resource=[json.loads(p.read_text()) for p in TMP.rglob('*.resources.json')]
    check('peak process RSS / GB',max(x['peak_rss_bytes'] for x in resource)/1e9,all(x['peak_rss_bytes']<6e9 for x in resource))
    check('peak sampled total / GB',max(x['tree_peak_bytes'] for x in resource)/1e9,all(x['tree_peak_bytes']<8e9 and not x['gate_reason'] for x in resource))
    check('figure artifact pairs',4,len(list((HERE/'figures').glob('t13-*.png')))==len(list((HERE/'figures').glob('t13-*.pdf')))==4)
    check('innermost cell R / M',a.H/math.sqrt(2),a.H/math.sqrt(2)>.0001)
    a.save('t13-verification.csv',rows);print(json.dumps(rows,indent=2),flush=True)

def figures():
    os.environ.setdefault('MPLCONFIGDIR',str(TMP/'mpl'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    try:
        import scienceplots
        plt.style.use(['science','no-latex'])
    except ImportError:pass
    plt.rcParams.update({'font.size':8,'axes.labelsize':8,'legend.fontsize':7,'axes.titlesize':9})
    figdir=HERE/'figures';figdir.mkdir(exist_ok=True)
    def finish(fig,name):
        fig.savefig(figdir/(name+'.pdf'));fig.savefig(figdir/(name+'.png'),dpi=600);plt.close(fig)
    # Separate columns preserve the sign and show the much smaller maximal signal.
    fig,ax=plt.subplots(2,2,figsize=(7,5),layout='constrained')
    for row,(ray,vec) in enumerate(a.NV.items()):
        for col,leg in enumerate(('original','maximal')):
            z=np.load(OUT/(leg+'-'+ray+'-profiles.npz'));d=(z['p6'][:,:,11:13]-z['p6'][0,:,11:13])@vec
            e=(z['p8'][:,:,11:13]-z['p8'][0,:,11:13])@vec
            for ti in (1,14,28,56):
                line=ax[row,col].plot(a.R,d[ti],label=f"{z['times'][ti]*1000:.3f}")[0]
                ax[row,col].fill_between(a.R,np.minimum(d[ti],e[ti]),np.maximum(d[ti],e[ti]),alpha=.2,color=line.get_color())
            ax[row,col].set(title=f'{ray}, '+('√χ' if leg=='original' else 'αK'),xlabel='R / M',ylabel=r'$\Delta\widetilde\Gamma_n / M^{-1}$')
            ax[row,col].legend(title=r'$t / (10^{-3} M)$');ax[row,col].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    finish(fig,'t13-gamma-disturbance')
    hist=read('t13-history.csv');fig,ax=plt.subplots(2,2,figsize=(7,4.6),layout='constrained')
    for row,ray in enumerate(a.NV):
        for leg in a.LEGS:
            z=[x for x in hist if x['run']==leg and x['ray']==ray];t=[float(x['time_M']) for x in z]
            for col,key in enumerate(('probe_Gamma','probe_C_Gamma')):
                ax[row,col].plot(t,[float(x[key]) for x in z],label=leg,ls='--' if leg.endswith('half') else '-')
                ax[row,col].set(xlabel='t / M',ylabel=key+r' / M$^{-1}$',title=ray)
                ax[row,col].set_yscale('symlog',linthresh=1e-6);ax[row,col].set_ylim(*((-.003,.15) if col==0 else (-.0005,.001)))
                ax[row,col].legend();ax[row,col].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    finish(fig,'t13-probe-history')
    fig,ax=plt.subplots(2,2,figsize=(7,4.5),layout='constrained')
    for row,(ray,vec) in enumerate(a.NV.items()):
        for leg in ('original','maximal'):
            b=np.load(OUT/(leg+'-'+ray+'-profiles.npz'));c=np.load(OUT/(leg+'-half-'+ray+'-profiles.npz'))
            d=(b['p6'][:,:,11:13]-b['p6'][0,:,11:13])@vec;e=(c['p6'][::2,:,11:13]-c['p6'][0,:,11:13])@vec
            error=e-d;spread=((b['p6'][:,:,11:13]-b['p6'][0,:,11:13])-(b['p8'][:,:,11:13]-b['p8'][0,:,11:13]))@vec
            for col,calc in enumerate((lambda x:np.max(abs(x),1),lambda x:np.sqrt(np.trapezoid(x*x,a.R,axis=1)/(a.R[-1]-a.R[0])))):
                ax[row,col].plot(b['times'][1:],calc(error)[1:],label=leg+' Δt/2')
                ax[row,col].plot(b['times'][1:],calc(spread)[1:],ls=':',label=leg+' P6/P8')
                ax[row,col].set(title=ray,xlabel='t / M',ylabel=('peak' if col==0 else 'RMS')+r' uncertainty / M$^{-1}$')
                ax[row,col].set_yscale('log');ax[row,col].set_ylim(1e-10,1e-4 if col==0 else 1e-5);ax[row,col].legend()
    finish(fig,'t13-temporal-control')
    budgets=read('t13-stage-budgets.csv');fig,ax=plt.subplots(2,2,figsize=(7,4.6),layout='constrained')
    for row,ray in enumerate(a.NV):
        for col,leg in enumerate(('original','maximal')):
            for group in ('geometric','advection','KO','total'):
                z=[x for x in budgets if x['run']==leg and x['ray']==ray and x['group']==group]
                ax[row,col].plot([int(x['step']) for x in z],[float(x['RMS_step_delta']) for x in z],marker='o',label=group)
            ax[row,col].set(title=ray+', '+leg,xlabel='native step',ylabel=r'RMS signed step profile / M$^{-1}$')
            ax[row,col].set_yscale('log');ax[row,col].set_ylim(1e-12,1e-2 if col==0 else 1e-5);ax[row,col].legend()
    finish(fig,'t13-stage-attribution')
    print('Four scientific figure pairs written.',flush=True)

if __name__=='__main__':
    {'attribution':attribution,'stage-c':stage_c,'predictor':predictor,'input-hashes':input_hashes,'verify':verify,'figures':figures}[sys.argv[1]]()
