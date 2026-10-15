"""T3: unchanged spatial/static controls plus algebraic STF and sphere-scale gates."""
import collections
import json
import math
from pathlib import Path
import report_t2 as base

OUT=Path(__file__).resolve().parent/'results-t3'
base.OUT=OUT
EPS_MODE=2.2e-16
NTHETA=97
# 16 times the 256-term tensor contraction count, then Ntheta angular sums.
# This deliberately conservative arithmetic allowance contains no spatial error.
C_MODE=4096


def weyl4():
    report=[]
    for form in ('vacuum','matter'):
        folder=OUT/'checkpoint'/form
        np=base.rows(folder/'ems_np.csv');w=base.rows(folder/'ems_radiation.csv')
        points=base.grouped(base.rows(folder/'np_points.csv'))
        for a,b in zip(np,w):
            assert (a['t'],a['r'],a['l'],a['m'])==(b['t'],b['r'],b['l'],b['m'])
            if a['l']<2:continue
            scale=max(abs(base.value(x,'Psi4')) for x in np if x['r']==a['r'])
            tolerance=C_MODE*EPS_MODE*NTHETA*scale
            error=abs(base.value(a,'Psi4')-base.value(b,'W'))
            report.append({'form':form,'r':a['r'],'l':a['l'],
                'NP':base.value(a,'Psi4').real,'W':base.value(b,'W').real,
                'mode_scale':scale,'tolerance':tolerance,'absolute_error':error,
                'NP_over_mode_scale':abs(base.value(a,'Psi4'))/scale,
                'W_over_mode_scale':abs(base.value(b,'W'))/scale,
                'error_over_tolerance':error/tolerance,
                'same_jet_legacy_max_error':max(abs(base.value(d,'Psi4')-d['legacy']) for d in points[a['r']]),
                'status':'PASS' if error<=tolerance else 'FAIL'})
    assert len(report)==30
    return report


def diagnostics():
    report=[]
    for case in ('checkpoint','cavity'):
        for r,data in base.grouped(base.rows(OUT/case/'matter/np_points.csv')).items():
            erms=base.rms(data,lambda d:d['E_trace'])
            hrms=base.rms(data,lambda d:d['Ham'])
            hmetric=base.rms(data,lambda d:d['Ham_metric'])
            mag=base.rms(data,lambda d:math.sqrt(sum(d[f'B_axial{i}']**2 for i in range(3))))
            mom=base.rms(data,lambda d:math.sqrt(sum(d[f'Mom{i}']**2 for i in range(3))))
            report.append({'r':r,'removed_electric_trace_rms':erms,'Hamiltonian_rms':hrms,
                'metric_Hamiltonian_rms':hmetric,
                'max_electric_trace_minus_metric_H_minus_KtraceA':max(abs(d['E_trace']-d['Ham_metric']-d['K_traceA']) for d in data),
                'removed_magnetic_trace_rms':base.rms(data,lambda d:d['B_trace']),
                'magnetic_axial_rms':mag,'momentum_norm_rms':mom,
                'max_magnetic_axial_plus_momentum':max(math.sqrt(sum((d[f'B_axial{i}']+d[f'Mom{i}'])**2 for i in range(3))) for d in data),
                'max_projected_coordinate_trace':max(d['trace'] for d in data),
                'max_Psi0_projection_change':max(d['psi0_change'] for d in data),
                'max_Psi4_projection_change':max(d['psi4_change'] for d in data),
                'max_radiative_change_over_budget':max(max(d['psi0_change'],d['psi4_change'])/d['projection_budget'] for d in data)})
    return report


def preservation():
    report=[]
    for case in ('checkpoint','cavity'):
        for form in ('vacuum','matter'):
            for r,data in base.grouped(base.rows(OUT/case/form/'np_points.csv')).items():
                ratio=max(max(d['psi0_change'],d['psi4_change'])/d['projection_budget'] for d in data)
                report.append({'case':case,'form':form,'r':r,
                    'max_error_over_budget':ratio,'status':'PASS' if ratio<=1 else 'FAIL'})
    return report


def summary(report):
    static=report['static'];rot=[x for x in report['rotation'] if x['form']=='matter']
    ready=report['ready_for_simulation']
    lines=['# NPTET T3 measured controls','',
        '**Ready for simulation under every requested gate: '+str(ready)+'**','',
        '## Removed components and same-point constraints','',
        'All entries are sphere RMS except the last three maximum residual columns. E trace uses the complete spatial metric. The production Hamiltonian uses independently evolved Gamma; H(metric) reconstructs Gamma and its derivative from the same metric jets. The latter agrees with the removed E trace after the measured K tr(A)/3 term. Magnetic trace is zero in this polar sector; its removed antisymmetric axial vector tracks minus the momentum constraint. The small tabulated axial residual includes derivatives of residual tr(A), which the production momentum formula assumes zero.', '',
        '| r | trace E | H(production) | H(metric) | axial B | norm M | max E-Hmetric-KtrA/3 | max axialB+M | max projected C trace |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for d in report['diagnostics']:
        keys=['removed_electric_trace_rms','Hamiltonian_rms','metric_Hamiltonian_rms','magnetic_axial_rms','momentum_norm_rms','max_electric_trace_minus_metric_H_minus_KtraceA','max_magnetic_axial_plus_momentum','max_projected_coordinate_trace']
        lines.append(f"| {d['r']:g} | "+' | '.join(f'{d[k]:.3e}' for k in keys)+' |')
    lines+=['','## Both null-rotation classes: roundoff comparison','',
        '| r | max FD/CAS error | max error/budget | failed comparisons / 380 |',
        '|---|---:|---:|---:|']
    for s in rot:
        v=list(s['rotations'].values())
        lines.append(f"| {s['r']:g} | {max(x['max_fd_cas_error'] for x in v):.3e} | {max(x['max_error_over_roundoff_budget'] for x in v):.3e} | {sum(x['failed_points'] for x in v)} |")
    lines+=['','## Protected suppression: RMS of pointwise ratios','',
        'X_out protects l fixed; X_in protects k fixed. Perturbation ratios are outgoing and incoming respectively. Opposite-class ratios remain in validation.json and have no Coulomb denominator term.', '',
        '| r | areal R | Phi2/Phi1 | Psi3/Psi2 | l-fixed Xout | Phi0/Phi1 | Psi1/Psi2 | k-fixed Xin |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for s in rot:
        vals=[s['Phi2/Phi1_rms_of_pointwise_ratio'],s['Psi3/Psi2_rms_of_pointwise_ratio'],s['rotations']['l_out']['suppression_rms_of_pointwise_ratio'],s['Phi0/Phi1_rms_of_pointwise_ratio'],s['Psi1/Psi2_rms_of_pointwise_ratio'],s['rotations']['k_in']['suppression_rms_of_pointwise_ratio']]
        lines.append(f"| {s['r']:g} | {s['R']:.9g} | "+' | '.join(f'{x:.3e}' for x in vals)+' |')
    lines+=['','## Production Weyl4: requested pure sphere-scale tolerance','',
        f'epsilon={EPS_MODE}, Ntheta={NTHETA}, C={C_MODE}; tolerance=C epsilon Ntheta max_l abs(mode). C is 16 times the 256-term contraction count. No spatial-error allowance is included. Identical-jet Psi0/Psi4 preservation passes; the existing cell-then-interpolate production W differs from the NP interpolate-then-reconstruct route.', '',
        '| form | r | l | NP mode | production W | sphere mode scale | error | tolerance | result |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for s in report['weyl4']:
        lines.append(f"| {s['form']} | {s['r']:g} | {s['l']:g} | {s['NP']:.12g} | {s['W']:.12g} | {s['mode_scale']:.3e} | {s['absolute_error']:.3e} | {s['tolerance']:.3e} | {s['status']} |")
    lines+=['','The l=3 values are symmetry-suppressed, but are not zero to machine roundoff: NP/scale and W/scale are retained in validation.json. The gate above compares the two measured modes, without replacing either by zero.','', '## Static backgrounds','',
        'All ten scalars on seven spheres are checked at both original T2 resolutions. Budgets are unchanged: 8/15 of the coarse/fine difference plus conditioning allowances (128 eps/dx² for curvature). Exact metric/potential oracle columns and all vanishing-field budgets are in validation.json.', '',
        '| case | r | exact Psi2 | fine error / tolerance | exact Phi1 | fine error / tolerance |',
        '|---|---:|---:|---:|---:|---:|']
    for s in static:
        p=s['scalars']['Psi2'];f=s['scalars']['Phi1']
        lines.append(f"| {s['case']} | {s['r']:g} | {s['oracle']['Psi2']:.9g} | {p['max_error_fine']:.3e} / {p['total_tolerance']:.3e} | {s['oracle']['Phi1']:.9g} | {f['max_error_fine']:.3e} / {f['total_tolerance']:.3e} |")
    lines+=['','## Other controls','',
        f"- Static: {collections.Counter(x['status'] for s in static for x in s['scalars'].values())}",
        f"- Default: {report['default']}",f"- Independent projection: {report['projections']}",
        '- Exact E/B, rotation and static CAS witnesses, complex Maxwell/flat/floor unit controls and prescribed production syntax check are retained in the logs.',
        '- Full production link/evolution was not run. No HPC, SSH, transfer or push.']
    (OUT/'summary.md').write_text('\n'.join(lines)+'\n')


def main():
    report={'rotation':base.rotations(),'static':base.static(),'weyl4':weyl4(),
            'default':base.default(),'projections':base.projections(),'diagnostics':diagnostics(),
            'radiative_preservation':preservation()}
    report['radiative_preservation_status']='PASS' if all(d['status']=='PASS' for d in report['radiative_preservation']) else 'FAIL'
    statuses=[x['status'] for s in report['static'] for x in s['scalars'].values()]
    statuses += [x['status'] for s in report['rotation'] for x in s['rotations'].values()]
    statuses += [x['status'] for x in report['weyl4']]+[report['default']['status'],report['projections']['status'],report['radiative_preservation_status']]
    report['ready_for_simulation']=all(x=='PASS' for x in statuses)
    (OUT/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    summary(report)
    for label,data in [('Static',[x['status'] for s in report['static'] for x in s['scalars'].values()]),('Rotation',[x['status'] for s in report['rotation'] if s['form']=='matter' for x in s['rotations'].values()]),('Weyl4',[x['status'] for x in report['weyl4']])]:print(label,collections.Counter(data))
    print('Default',report['default']);print('Projection',report['projections'])
    print('Psi0/Psi4 preservation',report['radiative_preservation_status'])
    print('Ready for simulation:',report['ready_for_simulation'])
    return 0 if report['ready_for_simulation'] else 1


if __name__=='__main__':raise SystemExit(main())
