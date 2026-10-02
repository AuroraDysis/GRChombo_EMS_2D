# Final T17 finder diagnosis and the exp-0024 policy

**Verdict: no tested boosted-single or binary individual-horizon setting qualifies strictly. The exp-0024 interim policy stands unchanged.** The final confirmation ladder admits neither T17 rung. The completed diagnosis is evidence about bounded cold searches of native snapshots; it does not establish nonexistence of the horizons or failure of every possible author-finder parameter choice. No new probe is started for this final report.

This report supersedes the pending interpretation in [the interim note](t17-finder-interim.md), while preserving that dated snapshot. The complete evidence is in [all stages and both individual surfaces](t17-finder-final-stages.csv), [one receipt-backed row per probe](t17-finder-final-probes.csv), [final-stage minima and residual histories](t17-finder-final-residuals.csv), [N48/N96 pair checks](t17-finder-final-pairs.csv), [initialization controls](t17-finder-final-controls.csv), [cost extrapolation](t17-finder-final-cost.csv), and [the audit result](t17-finder-final.json). The table at the end lists every probe. Pipeline exit zero is never used as a horizon pass.

The final audit checks all **58 diagnosis probes and eight original binary probes** against their own marker, resources and wait4 child receipts. It also verifies the four isolated production-native initialization receipts and both completed pipeline receipts. Actual return codes are one scientific success and 65 bounded scientific failures; every probe has a clear gate and peak RSS below 6e9 bytes. The former interrupted fine N48 accounting recovery is retained with explicit provenance and an unavailable tree peak, not a claimed zero measurement. The pinned data, finder/executable and author-source hashes match. The initialization controls have the same finest spacing, faces and rational subcell phase as the binary at R_h/h=14.138/28.276. Static input only initialized those controls at t=0; all horizon probes read frozen native checkpoints. No evolution, field setter, author-finder code, KO, gauge, transfer, convention or floor treatment is changed by this analysis.

All residuals are the author's area-weighted **mean squared outgoing expansion**, in native M_i units ([RHSurf.hpp:228](../../Source/RHFinder/RHSurf.hpp:228)). The stages are squared residual <=1e-7, 1e-10 and 1e-12, measured after fresh final interpolation. Passing all stages at each N is necessary; a qualified N48/N96 pair additionally needs area and charge agreement <=1e-3. Coarse/fine comparisons of failed surfaces cannot establish the spatial horizon test. The original T17 resolution verdict and launch results remain intact. The all-stage CSV contains each seed and each attempted stage per hole; unattempted stages are left unattempted rather than inferred from the final residual.

The **long confirmation** uses the predeclared selected accelerated setting, chase=32/quota=1, 100000 updates per stage and a 1795-second numerical case budget. It was selected by the smallest worst final residual among the accelerated ladder, not by outperforming chase=1. Boosted single t=0 finishes at squared residual 0.00239103/0.00353514 on coarse and 0.00239014/0.00353371 on fine (N48/N96), all UPDATE_CAP. Their finder times are 912.55/1069.33 seconds and 1178.66/1317.54 seconds, respectively. These results remain near the 595-second stalled values.

Binary long-confirmation results below are the worst of the two surfaces. Each reaches TIME_CAP at about 1795 finder seconds, fails the first stage and therefore never attempts the latter stages.

| Rung | Physical time M_i | N48 squared residual | N96 squared residual | Strict angular pair |
|---|---:|---:|---:|---|
| coarse | 0 | 0.0557856 | 0.0724400 | unavailable: both failed |
| fine | 0 | 0.0557409 | 0.0558764 | unavailable: both failed |
| coarse | 0.4375 | 0.131328 | 0.181089 | unavailable: both failed |
| fine | 0.4375 | 0.132639 | 0.140096 | unavailable: both failed |

Both binary surfaces are recorded separately at full precision in the CSVs. Close stopped area/charge values, where present, are not angular qualification. The physically later binary snapshots give larger residuals with this setting, but these cold searches do not diagnose a warm-track error or long-time physical growth rate.

## What limits the searches

**Cap and wrapper FLOOR.** The original chase=1 settings hit the 2000-update cap on boosted N48/N96 at both rungs; binary N48 also hits that cap. Original binary N96 instead triggers `FLOOR`, the local wrapper's comparison of extrema in successive 64-update windows ([t17-rh.cpp:238](t17-rh.cpp:238)). It can stop a residual that is still decreasing. This is neither an author convergence status nor an established physical/discretization floor. Making FLOOR inert within the 100000-update cap allows further reduction at chase=1: the 235-second coarse boosted controls reach 6.92364e-7/6.65118e-6, and binary reaches 4.62711e-4/4.89482e-3 (N48/N96, worst binary surface). They still decline when their time budget ends. Thus longer chase=1 convergence for boosted/binary is **not ruled out**, but no completed probe demonstrates it.

**Chase multiplier.** Chase=32 is demonstrably unsuitable for the strict schedule on these snapshots. Its residual stalls at a nonzero value even in unboosted controls, and extending to the full confirmation budget does not cure it. Chase=128/512 and the tested chase=128/quota=16 controls finish with very large residuals, commonly 1e4–1e6 for boosted/binary. More chasing cannot be assumed to mean faster convergence. The author's update contains the literal factor `courant*dtheta*dtheta*f*f` ([RHSurf.hpp:828](../../Source/RHFinder/RHSurf.hpp:828)); N96 takes smaller steps at a fixed multiplier. The smooth-flow stability region, the exact nonlinear mechanism of the chase=32 stalled iteration, and the optimal stable multiplier are not proved by this finite ladder. The author's Newton polish is commented out ([RHUnion.hpp:321](../../Source/RHFinder/RHUnion.hpp:321)); a Newton parameter does not activate a solver that is not executed.

**Seed.** Coarse binary N96 changes to spherical radius 0.9/1.1 R_h or centre offsets -/+0.1 R_h all end near squared residual 0.07243999 at chase=32. Those modest changes do not cure that regime. The source checks also confirm that author recentering moves the area-weighted centre ([RHSurf.hpp:617](../../Source/RHFinder/RHSurf.hpp:617)); the search subsequently changes its angular radius. The completed failed chase=1 boosted shape is nonspherical, so a spherical seed is not itself a demand that the final horizon be spherical. An improved independently supplied CTT shape or a different basin could still help. The current evidence does not exclude arbitrary seeds or isolate seed error from the full nonlinear iteration. Failed fitted shapes/centres are retained as iteration diagnostics, not qualified horizons.

**Tolerance, spatial resolution and data.** The unboosted coarse N48 chase=1 pass proves that 1e-12 is attainable by the unchanged finder on one native grid. No universal Float64 floor is established, and the finite cap/time failures do not show that the strict tolerances are physically impossible. Doubling spatial resolution barely changes original cap-limited boosted/binary residuals, while the unboosted fine controls reach much tighter residuals at that cap. The long chase=32 results also remain poor on both grids. These observations argue against simply buying an extra puncture level to cure the tested controls; they do not establish a spatial convergence order or irreducible resolution floor for a fully converged search.

The geometric lapse cannot directly cause the t=0 failure: lapse is absent from RHUnion's interpolated field list and the expansion formula ([RHUnion.hpp:160](../../Source/RHFinder/RHUnion.hpp:160), [RHSurf.hpp:430](../../Source/RHFinder/RHSurf.hpp:430)), and the stored geometric/sqrt non-lapse states are bit-identical. Boosted isolated controls fail without a CTT companion, so CTT is not the sole explanation. The CTT correction does change the spatial geometry and extrinsic curvature, hence can change the binary trapped surface and its search difficulty. These probes neither invalidate the data nor repair its documented interfaces. T17 tests d=32/v=0.0523381404947, whereas exp-0024 uses its separately audited d=16/v=0.05771783189084103 state.

## The exp-0024 contract remains unchanged

The controller's policy is retained **verbatim**:

> individual = APPROXIMATE warm tracking (N48, chase 1, quota 1, every 16 coarse steps; squared residual ≤ 1e-3 → APPROXIMATE else UNRESOLVED; CTT-side A_i = 0.2220848223, Q_i = 1.06931005935); common/remnant = STRICT (N48 + N96, stages 1e-7/1e-10/1e-12) before the 100 M_f clock and horizon-based losses; chase 32/128/512 excluded.

No completed probe argues for relaxing the common/remnant thresholds or reinstating an excluded multiplier. The approximate individual target is a reporting tolerance, **not a measured residual floor or A/Q error bound**. Warm tracking and its cadence were not tested by these cold frozen probes. A track remains approximate even if its stopped area/charge looks stable; values above the stated target remain UNRESOLVED. Preserve shapes, centres, residuals, finder state and checkpoints so a later strict search can revisit them. The unchanged author's FOUND flag must not be manufactured from the approximate label.

The independent initial values are the qualified d16 **CTT-side** reference in [horizons.toml](/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-echobin/artifacts/echo-binary/tranche3/d16/r32-a10/horizons.toml): A_1=0.22208482228786347, A_2=0.22208482228786353 M_i^2; Q_1=1.0693100593480518, Q_2=1.0693100593480516 M_i. They characterize that audited reconstructed CTT state, not a native-grid t=0 certification. They enter only offline initial-reference budgets, never an evolution RHS or post-t0 static target. The legacy RN mass column of failed searches is not the scalarized EMS loss/remnant mass.

A settled, nearly at-rest remnant remains **plausibly findable** at chase=1, supported by the unboosted strict N48 pass. This is an expectation, not a demonstrated common/remnant pair. The initial common surface may be highly distorted; a strict N48/N96 pair and the existing area/charge comparison remain mandatory. Neither an approximate individual track, an infall estimate nor a failed common search starts the 100 M_f clock. If no strictly qualified common/remnant surface is obtained, a horizon-based final loss/remnant verdict remains unavailable.

## Cost for strict common/remnant qualification

The only measured complete strict search is **unboosted coarse N48, chase=1/quota=1**: radius 0.0060404520035922523 M_i at the puncture, floor_window=100001 and cap=100000. Its stages take 1, 9797 and 9311 updates, with cumulative finder times 0.009375, 101.001626 and **201.542791 seconds**. Final squared residual is **9.9931491e-13**. Process wall time, including loading, is 201.843421 seconds. This is a local frozen isolated control, not a measured common horizon or exclusive cluster-node timing.

Unboosted coarse N96 at chase=1 reaches only stage 1, squared residual **3.8418643e-10**, after 235.003402 seconds. Thus the measured costs establish **more than 436.546 seconds** for the sequential local N48/N96 pair, not a completed strict-pair cost. No boosted/binary success supplies an alternative pair calibration.

For a provisional model, [t17-finder-final-cost.csv](t17-finder-final-cost.csv) fits log(squared residual) against update number over the last 50%, 25% and 10% of that N96 stage-1 history and converts updates using its observed mean time per update. Carrying those late slopes unchanged to 1e-12 projects N96 completion around **1256–1281 seconds**; adding the measured N48 gives **1458–1482 seconds, approximately 24–25 minutes per sequential local strict pair**. This is a conditional extrapolation, **not a confidence interval, rigorous bound or measured convergence**. A check on the completed N48 history projects 252–268 seconds from its stage-1 late slope, versus the measured 201.54 seconds, overpredicting by about 25–33%. N96's future slope and a common/remnant geometry may differ in either direction.

Use **about 1500 local seconds per pair as a provisional chase=1 cost term**, retaining retry/failed-search costs separately. It cannot be substituted directly for node seconds: native checkpoint size, interpolation refresh cost, angular dynamics, MPI layout and machine speed differ. These local controls do not justify replacing the controller's existing 5522–6164-second cluster cold-pair planning allowance. Until a chase=1 common/remnant pair is timed on the actual exclusive-node layout, retain that cluster allowance as a planning reserve, explicitly uncalibrated for this finder setting. No SSH, cluster job or multi-week merger run is started by this report.

## Complete probe table

E0/E1/E2 are the final fresh squared residuals at the three attempted stages, with the **worst of both surfaces** shown for binary cases. A dash means the stage was not attempted. The full-precision stage CSV contains each surface separately, its own status, A/Q, update count and cumulative seconds. Updates in the table refer to the last attempted stage; total updates and the actual child/gate receipts are in the probe CSV. Finder seconds are cumulative numerical-search time; process seconds include loading. O identifiers are original binary probes, D identifiers are the separate diagnosis. A capped or timed-out result remains a failure regardless of its stopped A/Q agreement.

<!-- ALL-PROBES -->

| ID | Data/rung/t/N | Chase/quota | Seed R/Rh, offset/Rh | Cap, floor, seconds | E0 | E1 | E2 | Stop, updates | Finder / process seconds | RSS GB |
|---|---|---|---|---|---:|---:|---:|---|---:|---:|
| O01 | binary/coarse/0/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.00348404 | — | — | UPDATE_CAP, 2000 | 88.31 / 89.07 | 2.0206 |
| O02 | binary/coarse/0/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.00927842 | — | — | FLOOR, 320 | 15.39 / 16.20 | 1.9397 |
| O03 | binary/coarse/0.4375/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.00547209 | — | — | UPDATE_CAP, 2000 | 94.44 / 95.16 | 2.0279 |
| O04 | binary/coarse/0.4375/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.0135505 | — | — | FLOOR, 512 | 31.56 / 32.36 | 1.9373 |
| O05 | binary/fine/0/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.00348156 | — | — | UPDATE_CAP, 2000 | 79.41 / 80.02 | 2.0729 |
| O06 | binary/fine/0/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.00927368 | — | — | FLOOR, 320 | 12.28 / 12.89 | 2.0736 |
| O07 | binary/fine/0.4375/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.0054507 | — | — | UPDATE_CAP, 2000 | 76.00 / 76.54 | 2.0732 |
| O08 | binary/fine/0.4375/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 0.0146038 | — | — | FLOOR, 192 | 7.56 / 7.96 | 2.0737 |
| D09 | binary/coarse/0/48 | 1/1 | 1, 0 | 100000, 100001, 235 | 0.000462711 | — | — | TIME_CAP, 6019 | 235.02 / 236.30 | 1.9573 |
| D10 | binary/coarse/0/48 | 128/1 | 1, 0 | 100000, 100001, 595 | 223510 | — | — | TIME_CAP, 13987 | 595.00 / 596.31 | 1.9554 |
| D11 | binary/coarse/0/48 | 128/16 | 1, 0 | 100000, 100001, 595 | 177451 | — | — | TIME_CAP, 9640 | 595.03 / 596.66 | 1.9575 |
| D12 | binary/coarse/0/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.0557856 | — | — | TIME_CAP, 48666 | 1795.03 / 1796.23 | 1.9551 |
| D13 | binary/coarse/0/48 | 32/1 | 1, 0 | 100000, 100001, 595 | 0.0557856 | — | — | TIME_CAP, 12699 | 595.01 / 596.33 | 1.9521 |
| D14 | binary/coarse/0/48 | 512/1 | 1, 0 | 100000, 100001, 595 | 87290 | — | — | TIME_CAP, 9741 | 595.05 / 596.36 | 1.9506 |
| D15 | binary/coarse/0/96 | 1/1 | 1, 0 | 100000, 100001, 235 | 0.00489482 | — | — | TIME_CAP, 5252 | 235.01 / 235.61 | 2.0303 |
| D16 | binary/coarse/0/96 | 128/1 | 1, 0 | 100000, 100001, 595 | 1.33119e+06 | — | — | TIME_CAP, 12810 | 595.04 / 595.73 | 1.9645 |
| D17 | binary/coarse/0/96 | 128/16 | 1, 0 | 100000, 100001, 595 | 1.13956e+06 | — | — | TIME_CAP, 10830 | 595.04 / 595.81 | 1.9326 |
| D18 | binary/coarse/0/96 | 32/1 | 0.9, 0 | 100000, 100001, 595 | 0.07244 | — | — | TIME_CAP, 11998 | 595.02 / 595.77 | 1.9325 |
| D19 | binary/coarse/0/96 | 32/1 | 1, -0.1 | 100000, 100001, 595 | 0.07244 | — | — | TIME_CAP, 9967 | 595.01 / 597.05 | 1.9471 |
| D20 | binary/coarse/0/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.07244 | — | — | TIME_CAP, 47956 | 1795.01 / 1795.88 | 1.9419 |
| D21 | binary/coarse/0/96 | 32/1 | 1, 0 | 100000, 100001, 595 | 0.07244 | — | — | TIME_CAP, 11530 | 595.03 / 595.94 | 1.9355 |
| D22 | binary/coarse/0/96 | 32/1 | 1, 0.1 | 100000, 100001, 595 | 0.07244 | — | — | TIME_CAP, 11891 | 595.03 / 596.01 | 1.9356 |
| D23 | binary/coarse/0/96 | 32/1 | 1.1, 0 | 100000, 100001, 595 | 0.07244 | — | — | TIME_CAP, 9972 | 595.02 / 595.79 | 2.0290 |
| D24 | binary/coarse/0/96 | 512/1 | 1, 0 | 100000, 100001, 595 | 1.24726e+06 | — | — | TIME_CAP, 9332 | 595.04 / 595.97 | 1.9412 |
| D25 | binary/coarse/0.4375/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.131328 | — | — | TIME_CAP, 38625 | 1795.02 / 1796.07 | 1.9376 |
| D26 | binary/coarse/0.4375/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.181089 | — | — | TIME_CAP, 31552 | 1795.04 / 1796.30 | 1.9404 |
| D27 | boosted/coarse/0/48 | 1/1 | 1, 0 | 100000, 100001, 235 | 6.92364e-07 | — | — | TIME_CAP, 23060 | 235.00 / 235.37 | 0.2952 |
| D28 | boosted/coarse/0/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 3.38884e-05 | — | — | UPDATE_CAP, 2000 | 19.92 / 20.25 | 0.2961 |
| D29 | boosted/coarse/0/48 | 128/1 | 1, 0 | 100000, 100001, 595 | 103651 | — | — | TIME_CAP, 33761 | 595.01 / 595.33 | 0.2888 |
| D30 | boosted/coarse/0/48 | 128/16 | 1, 0 | 100000, 100001, 595 | 538562 | — | — | TIME_CAP, 27642 | 595.02 / 595.44 | 0.2973 |
| D31 | boosted/coarse/0/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.00239103 | — | — | UPDATE_CAP, 100000 | 912.55 / 913.18 | 0.2951 |
| D32 | boosted/coarse/0/48 | 32/1 | 1, 0 | 100000, 100001, 595 | 0.00239103 | — | — | TIME_CAP, 54965 | 595.00 / 595.28 | 0.2951 |
| D33 | boosted/coarse/0/48 | 512/1 | 1, 0 | 100000, 100001, 595 | 124717 | — | — | TIME_CAP, 34358 | 595.01 / 595.58 | 0.2899 |
| D34 | boosted/coarse/0/96 | 1/1 | 1, 0 | 100000, 100001, 235 | 6.65118e-06 | — | — | TIME_CAP, 20932 | 235.00 / 235.25 | 0.2898 |
| D35 | boosted/coarse/0/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 9.08512e-05 | — | — | UPDATE_CAP, 2000 | 20.87 / 21.03 | 0.2968 |
| D36 | boosted/coarse/0/96 | 128/1 | 1, 0 | 100000, 100001, 595 | 336084 | — | — | TIME_CAP, 43762 | 595.01 / 595.45 | 0.2937 |
| D37 | boosted/coarse/0/96 | 128/16 | 1, 0 | 100000, 100001, 595 | 194190 | — | — | TIME_CAP, 26367 | 595.00 / 595.28 | 0.2893 |
| D38 | boosted/coarse/0/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.00353514 | — | — | UPDATE_CAP, 100000 | 1069.33 / 1069.49 | 0.2941 |
| D39 | boosted/coarse/0/96 | 32/1 | 1, 0 | 100000, 100001, 595 | 0.00353514 | — | — | TIME_CAP, 50941 | 595.01 / 595.34 | 0.2909 |
| D40 | boosted/coarse/0/96 | 512/1 | 1, 0 | 100000, 100001, 595 | 587517 | — | — | TIME_CAP, 28796 | 595.01 / 595.53 | 0.2861 |
| D41 | unboosted/coarse/0/48 | 1/1 | 1, 0 | 100000, 100001, 235 | 9.6357e-09 | 9.99692e-11 | 9.99315e-13 | FOUND, 9311 | 201.54 / 201.84 | 0.2953 |
| D42 | unboosted/coarse/0/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 9.6357e-09 | 8.88698e-10 | — | UPDATE_CAP, 2000 | 18.79 / 19.01 | 0.2958 |
| D43 | unboosted/coarse/0/48 | 128/1 | 1, 0 | 100000, 100001, 595 | 9.6357e-09 | 48643.5 | — | TIME_CAP, 39825 | 595.00 / 595.38 | 0.2963 |
| D44 | unboosted/coarse/0/48 | 128/16 | 1, 0 | 100000, 100001, 595 | 9.6357e-09 | 101919 | — | TIME_CAP, 32202 | 595.00 / 595.59 | 0.2915 |
| D45 | unboosted/coarse/0/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 9.6357e-09 | 0.00435793 | — | UPDATE_CAP, 100000 | 1293.01 / 1293.52 | 0.2941 |
| D46 | unboosted/coarse/0/48 | 32/1 | 1, 0 | 100000, 100001, 595 | 9.6357e-09 | 0.00435793 | — | TIME_CAP, 53384 | 595.00 / 595.34 | 0.2960 |
| D47 | unboosted/coarse/0/48 | 512/1 | 1, 0 | 100000, 100001, 595 | 9.6357e-09 | 1.03953e+06 | — | TIME_CAP, 40302 | 595.01 / 595.26 | 0.2892 |
| D48 | unboosted/coarse/0/96 | 1/1 | 1, 0 | 100000, 100001, 235 | 1.08424e-08 | 3.84186e-10 | — | TIME_CAP, 20883 | 235.00 / 235.28 | 0.2956 |
| D49 | unboosted/coarse/0/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 1.08424e-08 | 2.6232e-09 | — | UPDATE_CAP, 2000 | 20.37 / 20.59 | 0.2952 |
| D50 | unboosted/coarse/0/96 | 128/1 | 1, 0 | 100000, 100001, 595 | 1.08424e-08 | 25206.1 | — | TIME_CAP, 34802 | 595.01 / 595.76 | 0.2957 |
| D51 | unboosted/coarse/0/96 | 128/16 | 1, 0 | 100000, 100001, 595 | 1.08424e-08 | 627420 | — | TIME_CAP, 26900 | 595.02 / 595.52 | 0.2902 |
| D52 | unboosted/coarse/0/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 1.08424e-08 | 0.00496734 | — | UPDATE_CAP, 100000 | 985.93 / 986.31 | 0.2947 |
| D53 | unboosted/coarse/0/96 | 32/1 | 1, 0 | 100000, 100001, 595 | 1.08424e-08 | 0.00496734 | — | TIME_CAP, 49404 | 595.00 / 595.37 | 0.2959 |
| D54 | unboosted/coarse/0/96 | 512/1 | 1, 0 | 100000, 100001, 595 | 1.08424e-08 | 669204 | — | TIME_CAP, 32475 | 595.00 / 595.39 | 0.2953 |
| D55 | binary/fine/0/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.0557409 | — | — | TIME_CAP, 40958 | 1795.01 / 1796.67 | 2.0140 |
| D56 | binary/fine/0/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.0558764 | — | — | TIME_CAP, 40711 | 1795.01 / 1795.60 | 2.0720 |
| D57 | binary/fine/0.4375/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.132639 | — | — | TIME_CAP, 37581 | 1795.05 / 1795.62 | 2.0474 |
| D58 | binary/fine/0.4375/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.140096 | — | — | TIME_CAP, 37800 | 1795.03 / 1796.03 | 1.9744 |
| D59 | boosted/fine/0/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 3.4011e-05 | — | — | UPDATE_CAP, 2000 | 21.43 / 21.79 | 0.3203 |
| D60 | boosted/fine/0/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.00239014 | — | — | UPDATE_CAP, 100000 | 1178.66 / 1179.00 | 0.3202 |
| D61 | boosted/fine/0/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 9.09436e-05 | — | — | UPDATE_CAP, 2000 | 22.05 / 22.40 | 0.3198 |
| D62 | boosted/fine/0/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 0.00353371 | — | — | UPDATE_CAP, 100000 | 1317.54 / 1318.04 | 0.3201 |
| D63 | unboosted/fine/0/48 | 1/1 | 1, 0 | 2000, 64, 1795 | 4.59674e-11 | 4.59585e-11 | 2.16778e-12 | UPDATE_CAP, 2000 | 21.05 / 21.28 | 0.3191 |
| D64 | unboosted/fine/0/48 | 32/1 | 1, 0 | 100000, 100001, 1795 | 4.59674e-11 | 4.59585e-11 | 0.00435629 | UPDATE_CAP, 100000 | 1167.24 / 1167.74 | 0.3200 |
| D65 | unboosted/fine/0/96 | 1/1 | 1, 0 | 2000, 64, 1795 | 4.26398e-11 | 4.26397e-11 | 5.2744e-12 | UPDATE_CAP, 2000 | 21.98 / 22.18 | 0.3201 |
| D66 | unboosted/fine/0/96 | 32/1 | 1, 0 | 100000, 100001, 1795 | 4.26398e-11 | 4.26397e-11 | 0.00496537 | UPDATE_CAP, 100000 | 1166.84 / 1167.15 | 0.3196 |
