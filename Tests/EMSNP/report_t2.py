"""Assess the requested controls without converting a failed gate into a pass.

Roundoff budgets are fixed by the contraction arithmetic. Static spatial
error uses fourth-order Richardson scaling and the 1/dx^2 conditioning of
second differences. The Weyl4 comparison uses a separate sixth-order cell
calculation as an error estimate, with a factor-eight margin for the two
fourth-order evaluation orders (interpolate/reconstruct and its reverse).
"""
import collections
import csv
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results-t2'
EPS=sys.float_info.epsilon
NAMES=['Psi0','Psi1','Psi2','Psi3','Psi4','Phi0','Phi1','Phi2','X_out','X_in']


def rows(path):
    with Path(path).open() as stream:
        result=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]
    assert result and all(math.isfinite(x) for row in result for x in row.values()), path
    return result


def value(row,name):
    return complex(row[name+'_re'],row[name+'_im'])


def grouped(data):
    result=collections.defaultdict(list)
    for row in data:result[row['r']].append(row)
    return result


def rms(data,f):
    return math.sqrt(sum(d['weight']*abs(f(d))**2 for d in data)/sum(d['weight'] for d in data))


def ratio(a,b):
    return a/b if b else None


def rotations():
    report=[]
    for kind in ('checkpoint','cavity'):
        # The vacuum form is retained as a comparison, not called matter-complete.
        for form in ('vacuum','matter'):
            for r,data in grouped(rows(OUT/kind/form/'np_points.csv')).items():
                entry={'case':kind,'form':form,'r':r,'points':len(data),
                       'R':math.sqrt(sum(d['weight']*r*r*d['area'] for d in data)/(4*math.pi)),
                       'max_coordinate_trace':max(d['trace'] for d in data),'rotations':{}}
                for num,den in [('Phi2','Phi1'),('Psi3','Psi2'),('Phi0','Phi1'),('Psi1','Psi2')]:
                    entry[num+'/'+den+'_ratio_of_rms']=ratio(rms(data,lambda d:value(d,num)),rms(data,lambda d:value(d,den)))
                    entry[num+'/'+den+'_rms_of_pointwise_ratio']=rms(data,lambda d:value(d,num)/value(d,den)) if all(value(d,den)!=0 for d in data) else None
                for pre in ('l_out','l_in','k_out','k_in'):
                    errors=[abs(value(d,pre+'_fd')-value(d,pre+'_cas')) for d in data]
                    undefined=sum(d[pre+'_den']==0 for d in data)
                    entry['rotations'][pre]={
                        'max_fd_cas_error':max(errors),
                        'max_error_over_roundoff_budget':max(e/d[pre+'_budget'] for e,d in zip(errors,data)),
                        'failed_points':sum(e>d[pre+'_budget'] for e,d in zip(errors,data)),
                        'max_roundoff_budget':max(d[pre+'_budget'] for d in data),
                        'suppression_ratio_of_rms':ratio(rms(data,lambda d:value(d,pre+'_fd')),rms(data,lambda d:d[pre+'_den'])),
                        'suppression_cas_ratio_of_rms':ratio(rms(data,lambda d:value(d,pre+'_cas')),rms(data,lambda d:d[pre+'_den'])),
                        'suppression_rms_of_pointwise_ratio':rms(data,lambda d:value(d,pre+'_fd')/d[pre+'_den']) if not undefined else None,
                        'undefined_pointwise_ratios':undefined,
                        'status':'PASS' if all(e<=d[pre+'_budget'] for e,d in zip(errors,data)) else 'FAIL'}
                report.append(entry)
    return report


def static():
    report=[]
    for case,N,L in [('rn',256,16),('ems',1024,4)]:
        oracle={d['r']:d for d in rows(OUT/(case+'-oracle.csv'))}
        coarse=grouped(rows(OUT/f'{case}-n{N}'/'np_points.csv'))
        fine=grouped(rows(OUT/f'{case}-n{2*N}'/'np_points.csv'))
        dx=L/(2*N)
        for r,ds in fine.items():
            cs=coarse[r];ref=oracle[r]
            assert [d['theta'] for d in cs]==[d['theta'] for d in ds]
            entry={'case':case,'r':r,'dx':dx,'oracle':ref,'scalars':{}}
            # A Float64 second derivative has absolute error O(eps/dx^2),
            # irrespective of convergence once its truncation error falls below it.
            psi_roundoff=128*EPS/dx**2
            phi_roundoff=128*EPS*abs(ref['Phi1'])
            x_roundoff=2*abs(ref['Phi1'])*psi_roundoff+3*abs(ref['Psi2'])*phi_roundoff
            for name in NAMES:
                expected=ref[name] if name in ('Psi2','Phi1') else 0.
                ce=max(abs(value(d,name)-expected) for d in cs)
                fe=max(abs(value(d,name)-expected) for d in ds)
                difference=max(abs(value(c,name)-value(f,name)) for c,f in zip(cs,ds))
                arithmetic=psi_roundoff if name.startswith('Psi') else phi_roundoff if name.startswith('Phi') else x_roundoff
                discretization=8*difference/15
                tolerance=discretization+arithmetic
                entry['scalars'][name]={'max_error_coarse':ce,'max_error_fine':fe,
                    'max_normalized_fine_error':fe/max(abs(ref['Psi2']) if name.startswith('Psi') else abs(ref['Phi1']) if name.startswith('Phi') else abs(ref['Psi2']*ref['Phi1']),1e-300),
                    'observed_order':math.log2(ce/fe) if ce and fe else None,
                    'richardson_tolerance':discretization,'roundoff_tolerance':arithmetic,
                    'total_tolerance':tolerance,'status':'PASS' if fe<=tolerance else 'FAIL'}
            report.append(entry)
    return report


def weyl4():
    report=[]
    for form in ('vacuum','matter'):
        p=OUT/'checkpoint'
        np=rows(p/form/'ems_np.csv');w=rows(p/form/'ems_radiation.csv')
        high=rows(p/('sixth-'+form)/'ems_radiation.csv')
        points=grouped(rows(p/form/'np_points.csv'))
        assert len(np)==len(w)==len(high)==21
        for a,b,c in zip(np,w,high):
            assert (a['t'],a['r'],a['l'],a['m'])==(b['t'],b['r'],b['l'],b['m'])==(c['t'],c['r'],c['l'],c['m'])
            if a['l']<2:continue
            error=abs(value(a,'Psi4')-value(b,'W'))
            discretization=8*abs(value(b,'W')-value(c,'W'))
            arithmetic=4096*EPS*a['R']*rms(points[a['r']],lambda d:value(d,'Psi2'))
            tolerance=discretization+arithmetic
            report.append({'form':form,'r':a['r'],'l':a['l'],
                'NP':value(a,'Psi4').real,'W':value(b,'W').real,'W_sixth':value(c,'W').real,
                'absolute_error':error,'relative_error':ratio(error,abs(value(b,'W'))),
                'tolerance':tolerance,'order_estimate':discretization,'roundoff':arithmetic,
                'same_jet_legacy_max_error':max(abs(value(d,'Psi4')-d['legacy']) for d in points[a['r']]),
                'status':'PASS' if error<=tolerance else 'FAIL'})
    return report


def default():
    worst=0;count=0
    for form in ('vacuum','matter'):
        ds=rows(OUT/'default'/form/'ems_radiation.csv')
        original=rows(ROOT/'results'/form/'ems_radiation.csv')
        enabled=rows(OUT/'checkpoint'/form/'ems_radiation.csv')
        assert len(ds)==len(original)==len(enabled)==21
        for a,b,c in zip(ds,original,enabled):
            assert a.keys()==b.keys()==c.keys()
            for k in a:
                for ref in (b[k],c[k]):
                    error=abs(a[k]-ref);scale=max(abs(a[k]),abs(ref))
                    assert error<=1e-12*scale+1e-100,(form,k,a[k],ref)
                    worst=max(worst,error/scale if scale else 0)
            count+=1
        assert not (OUT/'default'/form/'ems_np.csv').exists()
        assert not (OUT/'default'/form/'np_points.csv').exists()
    return {'status':'PASS','rows':count,'relative_tolerance':1e-12,'absolute_tolerance':1e-100,'max_relative_difference':worst,'new_files_with_switch_off':False}


def harmonic(l,spin,theta):
    # Associated Legendre polynomial from its factorial coefficients, separate
    # from the production half-angle/binomial implementation. Its m=0 phase
    # is _1Y10=+sqrt(3/(8pi))*sin(theta), _-1Y10=-_1Y10.
    s=abs(spin)
    if l<s:return 0.
    x=math.cos(theta)
    derivative=sum((-1)**k*math.factorial(2*l-2*k)*x**(l-2*k-s)/
                   (2**l*math.factorial(k)*math.factorial(l-k)*math.factorial(l-2*k-s))
                   for k in range((l-s)//2+1))
    positive=math.sqrt((2*l+1)/(4*math.pi)*math.factorial(l-s)/math.factorial(l+s))*math.sin(theta)**s*derivative
    return positive*((-1)**s if spin<0 else 1)


def projections():
    spins=dict(zip(NAMES,[2,1,0,-1,-2,1,0,-1,-1,1]))
    assert abs(harmonic(1,1,.7)-math.sqrt(3/(8*math.pi))*math.sin(.7))<16*EPS
    count=0;worst=0
    for kind in ('checkpoint','cavity'):
        for form in ('vacuum','matter'):
            points=grouped(rows(OUT/kind/form/'np_points.csv'))
            for m in rows(OUT/kind/form/'ems_np.csv'):
                for name in NAMES:
                    terms=[d['weight']*harmonic(int(m['l']),spins[name],d['theta'])*value(d,name) for d in points[m['r']]]
                    expected=m['R']*sum(terms)
                    budget=2048*EPS*m['R']*sum(abs(x) for x in terms)
                    error=abs(value(m,name)-expected)
                    assert error<=budget,(kind,form,m['r'],m['l'],name,error,budget)
                    worst=max(worst,error/budget if budget else 0);count+=1
            for data in points.values():
                for d in data:
                    for name,a,b in [('X_out','Psi3','Phi2'),('X_in','Psi1','Phi0')]:
                        p=2*value(d,'Phi1')*value(d,a);q=3*value(d,'Psi2')*value(d,b)
                        assert abs(value(d,name)-(p-q))<=64*EPS*(abs(p)+abs(q))
    return {'status':'PASS','mode_components':count,'max_error_over_roundoff_budget':worst,
            'pointwise_X_definitions':'PASS','independent_harmonics':'associated Legendre factorial polynomial'}


def summary(report):
    lines=['# NPTET T2 measured controls','',
           '**Not ready for a simulation:** (b)(i) fails and two (c) modes exceed their declared error estimates.','',
           '## Static metric/potential oracles','',
           '| case | coordinate r | areal R | exact Psi2 | exact native Phi1 |',
           '|---|---:|---:|---:|---:|']
    for s in report['static']:
        o=s['oracle'];lines.append(f"| {s['case']} | {s['r']:g} | {o['R']:.9g} | {o['Psi2']:.9g} | {o['Phi1']:.9g} |")
    lines += ['', '## Fine-grid static residuals (70 scalar/sphere checks pass)','',
              'Psi2 and Phi1 cells show maximum absolute error / total tolerance. The other columns normalize maximum absolute residuals by their Coulomb field scale.', '',
              '| case | r | dx | Psi2 error / tol | Phi1 error / tol | max radiative Psi / abs(Psi2) | max Phi0,2 / abs(Phi1) | max X / abs(Phi1 Psi2) |',
              '|---|---:|---:|---:|---:|---:|---:|---:|']
    for s in report['static']:
        v=s['scalars'];a=v['Psi2'];b=v['Phi1']
        wp=max(v[x]['max_normalized_fine_error'] for x in ('Psi0','Psi1','Psi3','Psi4'))
        wf=max(v[x]['max_normalized_fine_error'] for x in ('Phi0','Phi2'))
        wx=max(v[x]['max_normalized_fine_error'] for x in ('X_out','X_in'))
        lines.append(f"| {s['case']} | {s['r']:g} | {s['dx']:g} | {a['max_error_fine']:.3e} / {a['total_tolerance']:.3e} | {b['max_error_fine']:.3e} / {b['total_tolerance']:.3e} | {wp:.3e} | {wf:.3e} | {wx:.3e} |")
    lines += ['', '## Matter checkpoint suppression: RMS of pointwise ratios','',
              'Both classes and both combinations are retained. k is the second null vector (called n in the task). Opposite-class denominators have no Coulomb leading term, so their O(1) ratios do not test the protected cancellation.', '',
              '| r | areal R | Phi2/Phi1 | Psi3/Psi2 | l-fixed Xout | l-fixed Xin | k-fixed Xout | k-fixed Xin |',
              '|---|---:|---:|---:|---:|---:|---:|---:|']
    selected=[s for s in report['rotation'] if s['form']=='matter']
    for s in selected:
        vals=[s['Phi2/Phi1_rms_of_pointwise_ratio'],s['Psi3/Psi2_rms_of_pointwise_ratio']]+[s['rotations'][p]['suppression_rms_of_pointwise_ratio'] for p in ('l_out','l_in','k_out','k_in')]
        lines.append(f"| {s['r']:g} | {s['R']:.9g} | "+' | '.join(f'{x:.3e}' for x in vals)+' |')
    lines += ['', '## Pointwise derivative roundoff gate (b)(i): FAIL','',
              '| r | maximum FD/CAS absolute error | maximum error / budget | failed comparisons / 380 | coordinate Weyl trace max |',
              '|---|---:|---:|---:|---:|']
    for s in selected:
        v=list(s['rotations'].values())
        lines.append(f"| {s['r']:g} | {max(x['max_fd_cas_error'] for x in v):.3e} | {max(x['max_error_over_roundoff_budget'] for x in v):.3e} | {sum(x['failed_points'] for x in v)} | {s['max_coordinate_trace']:.3e} |")
    lines += ['', '## Weyl4: l=2 matter modes pass; 28/30 total modes pass','',
              '| r | R Psi4 | old W | absolute error | tolerance | relative error |',
              '|---|---:|---:|---:|---:|---:|']
    for s in report['weyl4']:
        if s['form']=='matter' and s['l']==2:
            lines.append(f"| {s['r']:g} | {s['NP']:.12g} | {s['W']:.12g} | {s['absolute_error']:.3e} | {s['tolerance']:.3e} | {s['relative_error']:.3e} |")
    lines += ['', 'Failed weak modes:','', '| form | r | l | absolute error | tolerance |', '|---|---:|---:|---:|---:|']
    for s in report['weyl4']:
        if s['status']=='FAIL':lines.append(f"| {s['form']} | {s['r']:g} | {s['l']:g} | {s['absolute_error']:.3e} | {s['tolerance']:.3e} |")
    lines += ['', '## Output checks','',
              f"- Independent harmonic re-projection: {report['projections']['mode_components']} components pass; worst residual/budget {report['projections']['max_error_over_roundoff_budget']:.4g}. Both pointwise X expressions pass.",
              '- Default off and NP enabled: 42 old radiation rows unchanged within 1e-12 relative; maximum difference 0. No new file with NP off.',
              '- Exact CAS witnesses and production syntax check pass. No production link/evolution or continuum certification is claimed.',
              '- Complete budgets, observed convergence orders, incoming ratios, both Weyl forms and all modes are in validation.json.']
    (OUT/'summary.md').write_text('\n'.join(lines)+'\n')


def main():
    # Self-check the RMS and undefined-ratio rules independently of runtime data.
    assert rms([{'weight':1,'x':3},{'weight':3,'x':1}],lambda d:d['x'])==math.sqrt(3)
    assert ratio(1,0) is None
    report={'rotation':rotations(),'static':static(),'weyl4':weyl4(),'default':default(),'projections':projections()}
    report['ready_for_simulation']=all(x['status']=='PASS' for s in report['static'] for x in s['scalars'].values()) and all(x['status']=='PASS' for x in report['weyl4']) and all(x['status']=='PASS' for s in report['rotation'] if s['form']=='matter' for x in s['rotations'].values())
    (OUT/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    summary(report)
    print('Static controls:',collections.Counter(x['status'] for s in report['static'] for x in s['scalars'].values()))
    print('Weyl4 modes:',collections.Counter(x['status'] for x in report['weyl4']))
    print('Matter rotation sphere/classes:',collections.Counter(x['status'] for s in report['rotation'] if s['form']=='matter' for x in s['rotations'].values()))
    print('Default:',report['default'])
    print('Projections:',report['projections'])
    print('Ready for simulation:',report['ready_for_simulation'])
    return 0 if report['ready_for_simulation'] else 1


if __name__=='__main__':raise SystemExit(main())
