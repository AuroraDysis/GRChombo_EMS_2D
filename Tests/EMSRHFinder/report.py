"""Tabulate native finds against separate known-Killing-surface controls."""
import csv
import json
import math
from pathlib import Path
import sys


def load(path):
    return {(r['case'],int(r['N_theta']),int(r['resolution'])):r
            for r in csv.DictReader(path.open())}


def report(native_dir, control_dir, output):
    native,control=load(native_dir/'gates.csv'),load(control_dir/'gates.csv')
    expected={(case,n,dx) for case in ('static','plus','minus','rn') for n in (48,96,192) for dx in (32,64)}
    assert set(native)==expected
    assert set(control)=={(case,n,dx) for case in ('static','plus','minus','rn') for n in (48,96,192) for dx in (32,64,128)}
    rows=[]
    for key in sorted(expected):
        r,t=native[key],control[key]
        c32,c64,c128=(control[key[0],key[1],dx] for dx in (32,64,128))
        stop=abs(float(r['M_eq'])-float(t['M_eq']))
        grid=abs(float(c32['M_eq'])-float(c64['M_eq']))
        fine=abs(float(c64['M_eq'])-float(c128['M_eq']))
        order=math.log2(grid/fine) if grid and fine else float('nan')
        # No mass correction. 4/3 times the full coarse-fine change is the
        # second-order coarse-grid estimate, used conservatively for both grids.
        # The third grid tests whether that assumed order is supported.
        total=float(r['delta_M'])+stop+4*grid/3
        row=dict(r,search_offset_vs_known_surface=stop,grid_refinement_change=grid,grid_order=order,
                 delta_M_empirical=total,empirical_covers=float(r['mass_error'])<=total,
                 control_M_eq=t['M_eq'],control_mass_error=t['mass_error'],
                 control_quad_budget=t['delta_M'],control_A=t['A'],control_Q=t['Q'],
                 control_expansion_residual=t['expansion_residual'])
        rows.append(row)
    bykey={(r['case'],int(r['N_theta']),int(r['resolution'])):r for r in rows}
    for r in rows:
        control=bykey['static',int(r['N_theta']),int(r['resolution'])]
        r['boost_difference']=abs(float(r['M_eq'])-float(control['M_eq'])) if r['case'] in ('plus','minus') else 0.
        r['boost_combined_error']=r['delta_M_empirical']+control['delta_M_empirical']
        r['boost_agrees']=r['boost_difference']<=r['boost_combined_error']
    with output.with_suffix('.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    text=['# Isolated RHFinder: measured gates', '',
          'All rows use unchanged RHUnion/RHSurf, production EMSTRUMPET grid data and AMR Lagrange<4> interpolation at t=0. '
          'Expansion residuals are fresh test-only replays of the saved final shapes, without chase or re-centring. '
          'Eight Nθ=192 serialized native shapes marginally exceed the native 1e-7 mean-square threshold; the CSV retains their explicit status. '
          'Masses below come only from Julia `horizon_mass`. The displayed ± allowance is the requested round-sphere '
          'quadrature plus table allowance; it excludes the search stopping error. `delta_M_empirical` in the CSV adds '
          'the measured offset to a separate known-Killing-surface control and 4/3 of its full M/32–M/64 mass difference. M/128 tests the assumed at-least-second-order interpolation behavior. These controls are seeded with known geometry; they are not independent blind horizon finds.', '',
          '| Case | Nθ | dx/M | A | Q | e | M_eq ± δ_quad+table | mass error | M_RN | φ_mean | φ_rms (fluctuation) | model φ_H | R_A√〈Θ²〉 |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        f=lambda k:float(r[k])
        text.append(f"| {r['case']} | {r['N_theta']} | 1/{r['resolution']} | {f('A'):.8f} | {f('Q'):.9f} | {f('e'):.9f} | {f('M_eq'):.9f} ± {f('delta_M'):.2e} | {f('mass_error'):.2e} | {f('M_RN_legacy'):.9f} | {f('phi_mean'):.9f} | {f('phi_rms'):.2e} | {f('model_phi_H'):.9f} | {f('expansion_residual'):.2e} |")
    text += ['', '| Case | Nθ | dx/M | Control A | Control Q | Control M_eq ± δ_quad+table | Control mass error | ΔM_search | ΔM_grid | δM_empirical | Covers native error |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in rows:
        f=lambda k:float(r[k])
        text.append(f"| {r['case']} | {r['N_theta']} | 1/{r['resolution']} | {f('control_A'):.8f} | {f('control_Q'):.9f} | {f('control_M_eq'):.9f} ± {f('control_quad_budget'):.2e} | {f('control_mass_error'):.2e} | {f('search_offset_vs_known_surface'):.2e} | {f('grid_refinement_change'):.2e} | {f('delta_M_empirical'):.2e} | {r['empirical_covers']} |")
    text += ['', 'The empirical allowance is an engineering estimate from a known stationary surface and three spacings, '
             'not a rigorous bound or a branch/equilibrium certificate. All settled-remnant masses remain unset.']
    output.with_suffix('.md').write_text('\n'.join(text)+'\n')
    summary=dict(records=len(rows),mass_gates=all(r['mass_gate']=='true' for r in rows),
                 charge_gates=all(r['charge_gate']=='true' for r in rows),
                 quad_budget_covers=sum(r['budget_covers']=='true' for r in rows),
                 fresh_expansion_passes=sum(r['expansion_status']=='FRESH_REPLAY' for r in rows),
                 empirical_covers=sum(r['empirical_covers'] for r in rows),
                 boosted_agreement=all(r['boost_agrees'] for r in rows),
                 max_mass_error=max(float(r['mass_error']) for r in rows),
                 max_charge_error=max(float(r['charge_relative_error']) for r in rows),
                 control_max_mass_error=max(float(r['control_mass_error']) for r in rows),
                 minimum_grid_order=min(r['grid_order'] for r in rows),
                 max_area_readback_relative=max(abs(float(r['area_readback_relative'])) for r in rows))
    output.with_suffix('.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    report(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]))
