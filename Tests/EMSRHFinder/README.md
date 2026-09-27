# Isolated RHFinder and Julia mass post-processing

**The 24 native finds recover unit rest mass within the requested absolute mass and charge tolerances. The quadrature-plus-table error allowance alone fails to cover their mass errors.** The cause is the unchanged finder's stopping tolerance, with a smaller grid interpolation contribution. A separately measured engineering allowance covers all 24 discrepancies. This is not an unconditional pass of the original error-budget requirement.

This records the operator's original offline validation. No mass calculation is wired into C++. During that study, the author's finder, geometry, initial-data setter, production example and parameters were unchanged. The subsequent interpolation performance, restart and pilot-parameter changes are documented in [PERFORMANCE.md](PERFORMANCE.md). All masses below use the unchanged Julia `horizon_mass` from the read-only EMS worktree. That original study made no commits, SSH/HPC calls, evolution runs, radiation extraction, table adjustments or EMS.jl edits.

The complete requested case × angular resolution × spacing table is in [results/gates.md](results/gates.md), with full precision and all diagnostic columns in [results/gates.csv](results/gates.csv). [results/gates.json](results/gates.json) records the gate counts. Known-geometry controls are kept separate in [results/killing/gates.csv](results/killing/gates.csv).

## Route and inputs

The new test harness [harness/EMSRHSurfaceTest.cpp](harness/EMSRHSurfaceTest.cpp) builds a real, single-level Chombo hierarchy, following AMRInterpolatorTest, and calls the author's public `RHUnion` directly. The production `EMSBH_trumpet_read` setter fills the grid; the production `AMRInterpolator<Lagrange<4>>` supplies the finder. The evolution RHS is disabled. Repeated finder updates operate on the same t=0 snapshot until `FOUND` or the cap.

The reduced domain is x ∈ [0,4], y ∈ [0,2], with the object at (2,0), reflection at the cartoon axis, and uniform spacing M/32 or M/64. There is one active AMR level; the requested finest spacing is therefore the uniform spacing. The native matrix starts from spheres of coordinate radius 0.64 for EMS and 0.575 for RN, uses the author's mean-square expansion threshold 1e-7 and chase speed 1, and takes at most 1000 updates of 400 chase steps. These are searches from approximate radii, not exact-horizon seeds.

EMS uses the existing certified `artifacts/reference-alpha20-qfile.trumpet`: M=1, Q=0.700784945779, A=38.510457444941309, and native φ_H=0.052357655979104342. The three velocities are 0 and ±tanh(0.20273255), approximately ±0.2. RN was exported in seconds using Julia `reissner_nordstrom` → `trumpet_profile` → `write_trumpet`, then passed through exactly the same C++ setter. Its scalar is zero, with the same quadratic coupling and Q, and its exact horizon area is 36.890411235758208. The exporter and input digests are in [fixtures.jl](fixtures.jl) and [fixtures/manifest.sha256](fixtures/manifest.sha256). Extended precision was used only to construct the RN fixture; the finder and post-processing calculations use Float64.

Both production tables have 81 nodes on e ∈ [0.30,0.50]. The Julia reader binds each run's `ems_alpha=4π`, `f0=f1=0`, `f2=-20`, `phi_inf=0` and branch, and pins the approved payload hashes:

| Branch | Payload SHA-256 |
|---|---|
| scalarized-positive-n0 | `dd4013ca19a746b5bab847fcaaf7e9ece4d5ddf7f4ba5b29d84cb86386b23a07` |
| rn | `b98d6dc915b60656e848a6ca45d38bfb2f7d5ba2a5b6273395ab3ee3fe412c05` |

The Julia reader supplies its full format, convention, Smarr, ordering and hash validation. Its out-of-range errors propagate; the scripts do not clamp, extrapolate, substitute RN or change a table. Complete-file hashes, which differ from payload hashes, are in [results/provenance.json](results/provenance.json).

## Offline measurement and output record

[postprocess.py](postprocess.py) selects the last native row per search and time and requires `found`. It matches the corresponding `rh_f` row and t=0 plot file. Native `rh_surf` A and Q are used for the mass; Q is never reconstructed with a different normalization. Four-point tensor Lagrange interpolation of the plotted metric and φ, including the axis parity, supplies the scalar integral. Shape derivatives and midpoint area weights reproduce the author's `dA` expression. Recomputed A agrees with native A to at most 1.18e-9 relative, consistent with the native nine-significant-digit serialization.

`phi_mean` is ∑φ dA/∑dA. Following ruling §2.5, `phi_rms` denotes the RMS fluctuation about that mean, sqrt(∑(φ−phi_mean)² dA/∑dA), not sqrt(∑φ² dA/∑dA). The latter can be obtained as hypot(phi_mean,phi_rms). This reader deliberately supports only the uniform t=0 test hierarchy; it rejects multiple AMR levels, unsupported boundaries, incomplete plots and ambiguous time matches.

The standalone Julia script [postprocess.jl](postprocess.jl) computes R_A=sqrt(A/(4π)), e=Q/R_A and M_eq=R_A μ(e). For b_N=(π/(2N))/sin(π/(2N))−1, it passes δA=|A|b_N and δQ=|Q|b_N to `horizon_mass`. The reported primary allowance is exactly the requested first-law propagation plus the table allowance. `model_phi_H` and `model_Qs_over_R_A` are the table's linearly interpolated branch hints; their errors are not part of the certified mass allowance.

Each CSV has a self-describing first row. Its columns are:

| Columns | Meaning |
|---|---|
| case, N_theta, resolution, time | Case; angular points; inverse spacing in M units; snapshot time |
| horizon_id, search_index, duplicate_count | Geometric identity at this time, representative search and number merged |
| A, Q, R_A, M_irr, e | Native area/charge, area radius, R_A/2, charge ratio |
| M_eq, delta_M, delta_M_measurement, delta_M_table, delta_A, delta_Q | Julia mass and the unmodified first-law uncertainty components |
| phi_mean, phi_rms, model_phi_H, phi_difference, model_Qs_over_R_A | Scalar branch measurements and model hints |
| M_RN_legacy | Author's original RN-style mass, retained for comparison |
| expansion_squared_native, expansion_squared, expansion_residual, theta_minus, expansion_status, threshold | Native and freshly replayed mean-square outgoing expansion; dimensionless R_A sqrt(〈Θ+²〉); inward expansion; verification status; test threshold |
| area_from_plot, area_readback_relative | Independent shape/plot quadrature check |
| spread_A, spread_Q, duplicate_A_Q, spread_M_eq | Max–min spreads and constituent area/charge pairs |
| ems_alpha, f0, f1, f2, phi_inf, branch, table_hash | Run binding and table identity |
| mass_error, charge_relative_error, area_relative_error, RN_relation_error | Errors against the isolated reference; RN formula check |
| mass_gate, charge_gate, budget_covers, area_bias_covers | Explicit primary gate booleans, without hidden slack |
| equilibrium_status, M_settled, seconds | Static projection status; settled mass always unset; run duration |
| search_offset_vs_known_surface, grid_refinement_change, grid_order, delta_M_empirical, empirical_covers | Separate control-based error estimate, added only by report.py |
| control_M_eq, control_mass_error, control_quad_budget, control_A, control_Q, control_expansion_residual | Corresponding known-geometry control |
| boost_difference, boost_combined_error, boost_agrees | Difference from the unboosted mass at the same N and spacing |

Coincident searches are grouped by mutual geometric enclosure of their saved surfaces within 0.2% radial tolerance, not by search index. Distinct nested or separated surfaces are retained. IDs are local to a snapshot; this is not a worldline tracker. The dedicated two-search test produces one physical horizon with ΔA=0.0025204, ΔQ=3.0e-9 and ΔM_eq=2.65566e-5; see [results/duplicates/gates.csv](results/duplicates/gates.csv). No duplicate spread is silently absorbed into the quadrature allowance.

## Gates and cause of the failed allowance

| Check | Result |
|---|---|
| Native finds | 24/24 `found` at t=0 |
| Mass: ≤1e-3 at N=48, ≤3e-4 at N=96/192 | 24/24 pass; maximum error 2.54272e-4 |
| Relative charge error ≤5e-4 | 24/24 pass; maximum 1.78990e-4 |
| Boosted/unboosted rest mass | Pass even with the primary combined allowances; maximum difference 4.44349e-6 |
| RN table gives R_A(1+e²)/2 | Exact at the stored Float64 native inputs; control discrepancy at most Float64 roundoff |
| Scalarized case distinguished from RN | Native legacy mass 1.01570–1.01581, far outside the unit-mass gates |
| Requested quadrature-plus-table allowance covers mass error | **0/24: fail** |
| Separately measured engineering allowance covers native error | 24/24, with the limitations below |
| Fresh expansion on serialized geometry | 16/24 below 1e-7; eight N=192 rows marginally above, explicitly flagged |

At M/64, the native static mass errors for N=48,96,192 are 2.53751e-4, 1.74202e-4, 1.54324e-4. They decrease, then approach a stopping-error floor. Changing M/32 to M/64 changes native mass by only 1.83e-8 to 5.33e-7. The dominant floor therefore does not come from spatial resolution. The default stopping rule is 〈Θ+²〉≤1e-7, not a mass tolerance; the final dimensionless residuals are approximately 5.4–5.5e-4.

To isolate quadrature and interpolation, 36 additional controls use the known stationary Killing-horizon surface, at all three N values and spacings M/32, M/64 and M/128. Their positive radial seed is f(θ)=R_h/sqrt(cosh²(η)cos²θ+sin²θ), with R_h=0.63593977642346233 for EMS and 0.5701934158668871 for RN. Each is submitted to the unchanged finder and accepted in one update. **These are known-geometry seeded controls, not independent blind finds.**

At M/64 the static control mass errors are 1.06064e-4, 2.65224e-5 and 6.63605e-6; the boosted errors are 1.01814e-4, 2.54639e-5 and 6.37177e-6. Each angular doubling reduces these errors by approximately four. The RN control behaves similarly. Spatial differences of these controls are 1.01e-7 to 1.75e-7 in mass on the first pair; the third spacing supports at least second order (minimum measured order 2.12, limited by the native decimal output for tiny differences). Seven of the 36 known-surface controls also slightly exceed the bare allowance, so interpolation and the limitations of first-order propagation cannot be omitted even after stopping error is removed.

The reported engineering estimate is

    δM_empirical = δM_quad+table
                  + |M_native − M_known_surface|
                  + (4/3) |M_known_surface(M/32) − M_known_surface(M/64)|.

It leaves the native mass uncorrected. The last term uses the conservative second-order coarse-grid estimate for both spacings; M/128 checks the assumed order. It covers all native errors, with allowances from 1.40546e-4 to 2.54353e-4. It requires this known stationary reference and is not a rigorous bound or a general uncertainty estimator for an unknown binary horizon. It does not turn the failed bare-budget gate into a pass.

A stricter static N=48, M/32 pilot at the original chase speed reached an oscillatory expansion-error floor and failed its 1e-12 target after the bounded update limit. Reducing the existing chase-speed parameter to 0.125 reached 〈Θ+²〉=9.998e-17 in 15.56 s: A=38.5173783, Q=0.700910361, M_eq=1.000106630. Thus the dominant native bias can be reduced using the existing finder settings. The residual excess above the 1.06081e-4 primary allowance is about 5.5e-7, involving interpolation and finite angular representation. This pilot is separate from the native matrix. A full sweep at this strict setting is **BLOCKED by the declared runtime tier**: the measured chase scaling projects beyond 30 minutes (roughly 17 minutes for a single N=192 case); it was not launched. No upstream algorithm or table was adjusted.

## Consequences for evolution runs

The default stopping rule 〈Θ+²〉 ≤ 1e-7 is controlled by RHUnion's public member `m_thresh_super_low`, not an evolution run parameter, and leaves a positive relative mass bias of about 1.5e-4 on an isolated hole. Chase speed 0.125 reached 〈Θ+²〉 ≈ 1e-16 in 15.6 s at N=48; the saved successful pilot used a threshold of **1e-16**, not 1e-12 (the earlier speed-1 pilot failed its 1e-12 target). Quadrature bias falls approximately fourfold for each doubling of `RH_num_points`.

## Fresh expansion and decimal serialization

The test's verification mode reloads each saved centre and shape, interpolates fields afresh at those final points, and calls the author's public `expansion_error()` without a subsequent chase or re-centring. All 60 native/control surfaces were replayed; every inward expansion is negative. The gate table uses these fresh residuals, preserving the original native residual in a separate column. The raw checks are [results/fresh-native.json](results/fresh-native.json) and [results/fresh-killing.json](results/fresh-killing.json).

Eight serialized N=192 native surfaces exceed 1e-7 by less than 1.0e-11 absolute; the largest relative change from the native residual is 1.17e-4. A new static N=192, M/64 search confirmed the cause: the fresh in-memory residual and full-precision shape replay both give 9.99927888e-8, whereas replay of the author's nine-digit shape gives 1.0000390015e-7. Angular differentiation amplifies shape rounding. The full-precision check is [results/serialization-check.json](results/serialization-check.json). This directly verifies the mechanism for that case; the other seven small crossings remain explicitly flagged rather than waived. The native output precision was not changed.

## Area discriminator and scalar branch

The native default finds give approximately A=38.525 at N=192: a **positive** stopping bias. They do not satisfy an area error estimate containing only midpoint bias. They also do not produce A≈38.435.

The known-horizon controls give, at N=192 and M/64:

| Case | A | Q | φ_mean | model φ_H at measured e | φ_rms |
|---|---:|---:|---:|---:|---:|
| static | 38.5108874 | 0.700792780 | 0.05235765610 | 0.05235768797 | 2.28e-10 |
| +0.2 | 38.5108702 | 0.700792471 | 0.05235765609 | 0.05235767495 | 2.50e-10 |
| −0.2 | 38.5108702 | 0.700792471 | 0.05235765609 | 0.05235767495 | 2.50e-10 |

The midpoint relative biases at N=48,96,192 are +1.78509e-4, +4.46231e-5 and +1.11555e-5. The static M/64 area differs from this bias by only about 1e-8 relative; boosted angular quadrature has a slightly different coefficient and the same second-order reduction. All three isolated controls converge to the rest-frame A=38.510457444941309 and Q=0.700784945779. The exp-0004 deficit A≈38.435 (about −0.2%, persisting with spatial refinement) is absent here and has the opposite sign to the native finder bias. This supports attributing that deficit to the binary CTT data rather than an isolated setter/finder normalization error. The binary exp-0004 data were not rerun in this task, so this is a discriminator result, not a new binary validation.

Native scalar means are approximately 0.0523496, reflecting the small outward surface displacement. Native φ_mean−model φ_H ranges from −2.19e-6 to +2.85e-6, and the largest RMS fluctuation is 8.80e-8. On the known surface, φ_mean agrees with the true stationary value and is unchanged by the boost. The table's linear scalar hint and the quadrature-induced shift in measured e explain why model φ_H approaches the measured mean with N. These hints are not a branch proof. All records retain `BRANCH_UNVERIFIED_STATIC_PROJECTION`; `M_settled` remains unset.

## Validation

Cleanup validation rebuilt the harness and reran a native static case and an RN known-surface control at N=48, M/32 through fresh expansion replay and Julia post-processing. Both reproduce the previous A, Q, M_eq, uncertainty and scalar diagnostics; see [results/cleanup-check.json](results/cleanup-check.json).

- The offline interpolator's runnable checks cover even/odd polynomial interpolation through the axis, exact midpoint sphere quadrature at all three N values, and coincident/separated/nested surface grouping.
- The existing EMSTrumpet fixture suite passes the reference and B/E echo tables, rejection controls and asserted fingerprints. Reference objects/single CCZ4/binary CCZ4 are `e1107d452abc4e4b`, `e68444ed5f33ba7b`, `024d0092dba552bc`; see [results/trumpet-regression.log](results/trumpet-regression.log).
- EMSCTT passes all 42 physical and 42 CCZ4 rows, nine SHA controls, 30 reader/binding rejections and its CCZ4 fingerprint `b3b32ee0d0a729dc`; see [results/ctt-regression.log](results/ctt-regression.log).
- All 26 existing t=0 grid regression jobs exit zero. Their 116 norm rows are finite with zero chi/lapse floor hits. Production grid/kernel discrepancy is at most 1.69e-16; the intentionally old field-superposition control retains its expected finite difference. CTT fine collar H/M/GaussE orders are 3.99578/3.99340/3.99621. See [results/regression-summary.json](results/regression-summary.json) and [results/regressions](results/regressions). No evolution regressions were run.
- The static rebuild checks preserve the native shape output and accepted diagnostic row; a few intermediate near-zero quadratic diagnostics differ at roundoff in the N=192 control. See [results/final-binary-identity.json](results/final-binary-identity.json) and [results/serialization-check.json](results/serialization-check.json). No production build path changed.

An in-run C++ front end was dropped because post-processing of the unchanged finder's output suffices (operator, 2026-09-27).

The files added by the original study were this README, `GNUmakefile`, `harness/EMSRHSurfaceTest.cpp`, `run_cases.py`, `postprocess.py`, `postprocess.jl`, `report.py`, `regressions.py`, `fixtures.jl`, their fixture/result artifacts, and `ems_rh_seed.py`. That study did not modify the then-existing tracked files.

## Reproduction, provenance and runtime

Build with the local serial Chombo toolchain documented in `EMS-deps/BUILD.md`, `CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib`, and `OMP_NUM_THREADS=1`. From the fork root, the commands are:

```sh
make -C Tests/EMSRHFinder all DIM=2 -j 2
python3 Tests/EMSRHFinder/run_cases.py /private/tmp/emsrh-repeat/native matrix
python3 Tests/EMSRHFinder/run_cases.py /private/tmp/emsrh-repeat/controls killing
python3 Tests/EMSRHFinder/run_cases.py /private/tmp/emsrh-repeat/fresh-native verify /private/tmp/emsrh-repeat/native
python3 Tests/EMSRHFinder/run_cases.py /private/tmp/emsrh-repeat/fresh-controls verify /private/tmp/emsrh-repeat/controls
uv run --with h5py --with numpy python Tests/EMSRHFinder/postprocess.py /private/tmp/emsrh-repeat/native /private/tmp/emsrh-repeat/native-results /Users/auroradysis/Workspace/EMS /private/tmp/emsrh-repeat/fresh-native/fresh.json
uv run --with h5py --with numpy python Tests/EMSRHFinder/postprocess.py /private/tmp/emsrh-repeat/controls /private/tmp/emsrh-repeat/control-results /Users/auroradysis/Workspace/EMS /private/tmp/emsrh-repeat/fresh-controls/fresh.json
python3 Tests/EMSRHFinder/report.py /private/tmp/emsrh-repeat/native-results /private/tmp/emsrh-repeat/control-results /private/tmp/emsrh-repeat/gates
```

Run builds and matrices in the background and poll with backoff, as required by the operator. The driver caps each executable at 240 seconds. A nonzero verification result for the eight flagged rounded shapes is expected and is retained in `fresh.json`; the mass post-processor still reports their native `found` rows with that failure status. A run is reused only if both its parameter text and executable SHA match. Known-geometry controls additionally require an explicit seed acknowledgement in the log, preventing silent use of a stale executable.

The delivered native matrix took 613.49 seconds total (0.81–78.94 seconds per run); the four pilot cases provided the dry run for this 5–30 minute batch. The 36 controls took 4.35 seconds total. Existing grid regressions took 57.06 seconds total, at most 23.86 seconds per job. The final serialization diagnostic took 53.02 seconds. All actual individual runs fit the under-five-minute tier. The projected full strict sweep was not run.

Compact accepted native rows, shapes and parameters are retained under `results/native/accepted` and `results/killing/accepted`. Their `inputs.json` files record SHA-256 hashes of the complete original files, including plots. Full plots/logs/statuses remain at `/private/tmp/emsrh-t2/postprocess-native` and `/private/tmp/emsrh-t2/killing-controls`; these temporary locations must be retained or the runs regenerated to repeat the scalar integration. Strict/duplicate/serialization controls have separate sibling directories. There are no outstanding background runs.

The algebra-only verification card is: fixed-branch first-law chain rule and Hermite endpoint identities use cancelled rational residuals/zero polynomial coefficients; the transformed-sphere seed uses the exact residual γ²f²cos²θ+f²sin²θ−R_h²=0 and D=1+(γ²−1)cos²θ≥1 on R_h>0, γ≥1. The positive radial root is selected. Independent numerical checks include 2016 deterministic Float64 seed points and direct midpoint sphere quadrature. Status: proved algebraic identities and passed numerical transfer checks. See [results/cas.log](results/cas.log), [results/seed-cas.log](results/seed-cas.log), and [ems_rh_seed.py](ems_rh_seed.py). **CAS scope: algebraic identity — production physics not certified.**
