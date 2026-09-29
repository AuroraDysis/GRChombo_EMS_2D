# T1: local qualification

**Status: BLOCKED.** The E data file is absent during a successful bitwise restart with a complete hierarchy, but an incomplete-hierarchy restart still calls the EMSTRUMPET reader. All B level-15 grid builds segfault before the finest level exists. The E t = 0 constraints converge at fourth order only in the horizon shell; the full cavity, outer-ring and far-field masks do not. No gauge, evolution equation, AMR transfer, reader, or finder code was changed.

Build: `CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib make -C Examples/EMS all DIM=2 -j4` (serial Chombo, Float64). Runs used four OpenMP threads except the single-thread restart check. The two source files are read-only inputs: E SHA-256 `2a8de074ae17c0b11d323d4b0933a6bdb7430a5305473c8cc7d4ce37f39fa793`, B `405a2b0e1f2f28b7dea2e3757ee1780cafdcb91538aba397b251cbbc10d5ba86`. Every supplied parameter file pins `sigma = 1`. Final run logs have no `Variable with name ... not found` or `Parameter: sigma not found` line; other optional parameters report their documented defaults.

## A. Data-file reachability and restart

| Call site (tree-relative unless noted) | Role and reachability |
| --- | --- |
| `Examples/EMS/EMSBH2DLevel.cpp:91-102` | Only production construction, `compute_1d_solution()`, and grid setter; all inside `initialData()`. |
| `Source/InitialConditions/EMSBH/EMSBH_trumpet_read.impl.hpp:29-31,50` | Initial-data reader delegates to profile `main()` and optional binary CTT companion. |
| `Source/InitialConditions/EMSBH/1D_SOL/EMSTrumpetSolution_read.impl.hpp:97-123` | `main → check_file → read_from_file → ifstream`; actual file access. |
| `Tests/EMSTrumpet/EMSTrumpetFixtureTest.cpp:230-237,293`; `EMSTrumpetGridConvergence.cpp:107-125`; `EMSTrumpetBinaryGridConvergence.cpp:136-137` | Standalone test calls, never evolution callbacks. |
| `Tests/EMSCTT/EMSCTTFixtureTest.cpp:331-332`; `EMSCTTGridConvergence.cpp:87-88` | Standalone CTT test calls, never evolution callbacks. |
| `Chombo/lib/src/AMRTimeDependent/AMR.cpp:456,1656,1757,1764` | Calls `initialData()` during new-run initialization. |
| `Chombo/lib/src/AMRTimeDependent/AMR.cpp:693-698` | **Failure path:** after checkpoint read, calls `initialData()` on every undefined level above the checkpoint's finest level, so the data file is required after t = 0 in that case. That path was stopped, not exercised. |

The exhaustive `rg -n` symbol and method audit is in `reader-sites.csv` (170 file:line records, including declarations and tests). No reader construction or call occurs in production tagging (`Examples/EMS/EMSBH2DLevel.cpp:259-281`), RHS, diagnostic hooks, RHFinder or AMR transfer. The T1 hierarchies use geometric tagging only during initialization, with all `regrid_interval` entries zero and derivative-tag thresholds `1e100`.

| E restart experiment | Result |
| --- | --- |
| `a-first.txt`: 2 coarse steps, then `a-restart.txt` to step 4 with the copied file renamed away; `a-continuous.txt`: direct 4 steps | Both exit 0; restart log contains no `Read EMSTRUMPET` line. |
| Final raw HDF5 checkpoint comparison | Byte-identical, SHA-256 `7066f6208aff3f3795663867dcd1e05f445ca3aeaffd8c7f7717f32c4e92edd5`; `h5diff` also reports no differences. |
| Scope | `max_level = finest_level = 0`, `h0 = 0.125 M`, `dt0 = 0.03125 M`, 1 rank × 1 thread. This proves the exercised complete-hierarchy restart path, not the incomplete-hierarchy path above. |

## B. Actual hierarchy census and cost

Domain `x ∈ [-224,224] M`, cartoon `y ∈ [0,224] M`; `512 × 256` base, `h0 = 0.875 M`, `dt0 = 0.21875 M`, 458 coarse steps modeled to 100 M, fixed ratio 2 and no radial map. Patch width is the nominal x-cell extent on each refined level; real boxes include Chombo tag and nesting buffers. `census-levels.csv` gives **every real level**: composite valid cells, total box cells, number of boxes, 3-layer ghost cells, `dx`, `dt`, and level steps. `W = 458 Σ_l (box cells)_l 2^l`, matching Chombo's updated-cell convention. Planned peak is `2 × (4 × 28 + 18) × 8 × Σ(box cells + ghost cells)` bytes: four RK evolution arrays, one diagnostic array, plus 100% headroom for patchers and allocator overhead. These are estimates, not measured RSS.

| Member | Max level | Patch width | W (cell steps) | Finest cells/R_h | Planned peak GiB | Inner ring level (cells/R) | Outer ring level (cells/R) | Real build |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| E | 12 | 384 | 3.0014e11 | 28.22 | 2.368 | 12 (125.3) | 7 (108.0) | yes |
| E | 12 | 192 | 8.1202e10 | 28.22 | 0.841 | 11 (62.6) | 6 (54.0) | yes |
| E | 12 | 96 | 2.3586e10 | 28.22 | 0.448 | 10 (31.3) | 6 (54.0) | yes |
| E | 13 | 384 | 6.0030e11 | 56.44 | 2.542 | 12 (125.3) | 7 (108.0) | yes |
| E | 13 | 192 | 1.6236e11 | 56.44 | 0.888 | 11 (62.6) | 6 (54.0) | yes |
| E | 13 | 96 | 4.7119e10 | 56.44 | 0.463 | 10 (31.3) | 6 (54.0) | yes |
| B | 14 | 384 | 1.2006e12 | 26.57 | 2.716 | 14 (136.3) | 7 (134.3) | yes |
| B | 14 | 192 | 3.2469e11 | 26.57 | 0.935 | 13 (68.1) | 6 (67.1) | yes |
| B | 14 | 96 | 9.4183e10 | 26.57 | 0.477 | 12 (34.1) | 5 (33.6) | yes |
| B | 15 | 384 | ≈2.4012e12 | 53.14 | ≈2.890 | unavailable | unavailable | SIGSEGV |
| B | 15 | 192 | ≈6.4934e11 | 53.14 | ≈0.981 | unavailable | unavailable | SIGSEGV |
| B | 15 | 96 | ≈1.8831e11 | 53.14 | ≈0.491 | unavailable | unavailable | SIGSEGV |

All B level-15 attempts stop while Chombo tags level 14 (`AMR.cpp:1656-1692`); the p96 case was repeated at one thread and exited 139. The level-15 `W` and memory entries are **unverified extensions** of the corresponding measured level-14 hierarchy with one equal-size top patch, not a real box census. Full-ring coverage is tested over 721 angles against the actual boxes. `rings.csv` derives isotropic light-ring radii from the delivered file's Chebyshev `Y(s)` and the independent areal-ring table. It is a pre-run resolution census only.

## C. t = 0 constraints

E uses `c-E-low/mid/high.txt` at `h0 = 0.875, 0.5833333333, 0.3888888889 M` (factor 3/2) and the **same physical** fixed tagging radii on all levels. The reference patch is nominal 192 × 96 cells on the low-resolution hierarchy. Norms are composite-grid cylindrical RMS: `sqrt(Σ y |C|² / Σ y)`, with covered coarse cells removed. Momentum is `sqrt(Mom1² + Mom2²)`; GaussE and GaussB are separate. Masks are fixed coordinate shells: horizon `[R_h,2R_h]`; inside inner ring `[2R_h,R_inner]`; inner-side cavity `[R_inner,R_stable]`; between unstable rings `[R_inner,R_outer]`; outer ring `[0.8 R_outer,1.2 R_outer]`; far `[4,8] M`. These cover both common meanings of “inner cavity.” `constraints-t0.csv` includes counts and norms; `orders-t0.csv` includes both observed orders. GaussB is identically zero at t = 0, so no order is assigned.

| E mask | Hamiltonian order (low→mid, mid→high) | Momentum order | GaussE order |
| --- | --- | --- | --- |
| Horizon | 4.007, 4.004 | 4.005, 4.006 | 4.017, 4.008 |
| Inside inner ring | 0.348, 0.391 | 1.115, 1.239 | 1.308, 1.368 |
| Inner-side cavity | 0.464, 0.478 | 1.469, 1.532 | 1.413, 1.445 |
| Between unstable rings | 0.477, 0.488 | 1.480, 1.540 | 1.424, 1.452 |
| Outer ring | 0.275, 0.353 | 1.173, 1.284 | 1.356, 1.277 |
| Far field | 0.424, 0.409 | 1.295, 1.314 | 1.309, 1.342 |

The CSV also records narrower registered core shells; they do not repair the low orders. The full-mask result is a **failed convergence gate**. Large t = 0 Hamiltonian values occur in cells near refinement patch corners; for example, E low's level-11 cavity cells reach `|Ham| = 2.19` while level-10 cavity interior reaches only `8.03e-7`. This is an observed AMR interface residual, not a fourth-order pass.

| B mask | p96/p384 Hamiltonian | p96/p384 momentum | p96/p384 GaussE |
| --- | ---: | ---: | ---: |
| Horizon | 1.00 | 1.00 | 1.00 |
| Inside inner ring | 1.02e6 | 2.06e3 | 2.52e4 |
| Inner-side cavity | 1.98 | 20.03 | 6.99 |
| Between unstable rings | 1.91 | 19.55 | 6.95 |
| Outer ring | 2.00 | 7.61 | 9.58 |
| Far field | 1.92 | 6.35 | 6.23 |

`patch-effects-t0.csv` also gives the p192 ratios and raw B norms. The compact B patch degrades constraints outside the horizon shell, catastrophically between `2R_h` and the inner light ring: Hamiltonian RMS `1.75e-6 → 1.79` from p384 to p96. The tested E p192 hierarchy fails the full-mask convergence gate; no qualifying E patch was established.

## D. Bounded E smoke (diagnostic only)

| Hierarchy | Planned peak | Wall watchdog | Completed coarse steps | Peak RSS | Finite fields / floor counts | Constraint change |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| E level 12, nominal 192 × 96 | 0.841 GiB | 1200 s | 1 (4096 finest steps) | 0.432 GiB | 0 nonfinite evolved values; 0 saved floor-contact cells | Ham 0.766×, Mom 2.975×, GaussE 15.752× |

The run uses `d-E12-p192.txt` and is explicitly **not** C-qualified. It has fixed geometric boxes, no regridding, 28 evolved variables and 5 constraints in t = 0/end plots, checkpoint every coarse step, and no static profile use after initialization. End diagnostics use global composite norms without any profile-derived mask. `smoke-levels.csv` reports each level's exact completed step count (1, 2, …, 4096); `smoke.csv` reports t = 0 and t = 0.21875 M norms and field checks. The run exited 0 after 739.61 s and 177,297,152 updated cells, yielding 0.18057 s per finest-level step and 239,717 updated cells/s including initialization and I/O. On the 325,760 composite cells, global RMS changed: Ham `2.50574e-3 → 1.91843e-3`, Mom `0.116349 → 0.346125`, GaussE `4.63098e-4 → 7.29487e-3`; GaussB remained zero. End minima are χ `1.6385e-6` and lapse `2.8705e-3`, above their `1e-12` floors. The existing clamp has no activation counter; saved-state floor contacts are zero, while transient activation events are **unmeasured**.

## E. Throughput

| Source | Updated cells | Wall time | Rate |
| --- | ---: | ---: | ---: |
| `Chombo/lib/src/AMRTimeDependent/AMR.cpp:747-763,1186-1195` | Prints per-level and total counts only if `verbosity >= 2`; counts all box cells at every advance. | — | — |
| exp-0002 n64 `pout.0` | Absent (`verbosity = 0`); no recoverable counter. | 1.95 h for 100 M reported by operator | Cannot infer cell updates without the moving-grid history. |
| Local E diagnostic smoke, 4 threads, includes setup/plots | 177,297,152 | 739.61 s | 239,717 cells/s |
| Short exp-0002-like reference, 2 steps, 4 threads, includes setup/regrid | 15,161,344 | 74.08 s | 204,658 cells/s |

The local reference peak RSS was 0.744 GiB. These two short whole-process rates include different setup and output costs; they are measured local rates, not steady-state cluster throughput. `throughput.csv` records the conditions and exact counts.

## F. Lower-x mirror plane (read-only)

| Question | Finding |
| --- | --- |
| Does the implementation support it? | `BoundaryConditions.cpp:649-673` mirrors `i → -i-1` at the lower boundary and multiplies each component by its `vars_parity` sign. The supplied 28-entry array covers every evolved EMS variable. |
| Odd under x reflection | `h12`, `A12`, `Gamma1`, `shift1`, `B1`, `Lambda`, `By`, `Bz`, `Ex`; the remaining evolved variables are even under x. Magnetic `B` is axial, so `Bx` is even and `By/Bz` odd. |
| What a separate study would require | Use x-domain `[0,L_x]` with puncture/coordinate, tagging, extraction and finder centers on x=0; set `lo_boundary = 2 2`, retain reflective y=0, and verify output/diagnostic parities. The supplied diagnostic parity array leaves `Sx` even although an x-directed momentum density should be odd, so it needs an audit before diagnostic ghost use. No lower-x run was made here. |

## Artifacts

`params/` contains every run parameter file. `census-levels.csv` contains real built levels; `census-cases.csv` marks B15 model-only rows as failed; `rings.csv`, `reader-sites.csv`, `constraints-t0.csv`, `orders-t0.csv`, `patch-effects-t0.csv`, `smoke.csv`, `smoke-levels.csv`, and `throughput.csv` hold the call-site, geometric, t = 0, smoke, and rate measurements. `COMMIT-MANIFEST-T1.txt` hashes all delivered files except itself. Raw run logs and HDF5 outputs are under `/private/tmp/ems-t1-{a,b,c,d,e}` and are not included in the worktree.

# T2: local qualification

**Status: BLOCKED for E evolution.** The t = 0 defect is localized to coarse–fine ghost interpolation, but no tested remedy meets both the convergence and evolution gates. The opt-in after-t=0 initial-data guard and a parameter-only B15 grid remedy pass. No equation, gauge, RHFinder, Chombo library, radial coordinate, or data file was changed.

All runs are Float64 with `sigma = 1`. Input SHA-256 prefixes are E `2a8de074ae17c0b1`, B `405a2b0e1f2f28b7`, reference `6f0820a576620f1f` (full hashes in `COMMIT-MANIFEST-T2.txt`). E uses `max_level = 12`, `h0 = 0.875, 0.583333333333, 0.388888888889 M` on the registered 3/2 sequence, the same physical geometric tagging radii, and `regrid_interval = 0`. The extra `h0 = 0.259259259259 M` grid checks the Float64 limit. All norms use T1's fixed coordinate masks and composite cylindrical RMS. `t2-layout.csv` records 59 actual E/reference levels; every measured box union is rectangular.

## T2 A. Ghost source and cell localization

| Cells whose constraint stencil reads ghosts | t = 0 ghost source |
| --- | --- |
| Fine patch edge and convex corner | `GRAMRLevel::fillAllEvolutionGhosts` calls `FourthOrderFillPatch::fillInterp` first. `FourthOrderFillPatch.cpp:254-300` selects uncovered coarsened ghosts; `FourthOrderInterpStencil.cpp:145-195` builds weights with **cell-average** polynomial moments (`power1dcoarseind0avg`, `power1dfineind0avg`). EMS stores point values. |
| Internal seam between fine boxes | `m_state_new.exchange(a_comps, m_exchange_copier)` overwrites interpolated values where another fine box has valid cells. |
| Cartoon axis `y = 0` | `BoundaryConditions::fill_reflective_cell` mirrors ghosts with `vars_parity`; cells where the axis meets a patch x edge also need coarse interpolation. |
| Outer x and upper-y boundary | Sommerfeld boundaries: `fill_solution_boundaries` skips them at t = 0. The direct `initialData()` setter supplied these ghosts at t = 0; later RHS boundaries use Sommerfeld. |
| Re-entrant union corner | None in the measured E p192/p384 or reference hierarchies; no residual or order is assigned. |

The E p192 low level-11 union is a 208 × 104-cell rectangle in two boxes, with an internal same-level seam and exterior convex corners. Its block factor is 8, maximum box size 128, original grid buffer 8; the reference block factor is 16. Classes use two valid cell layers and remove covered coarse cells. This table gives **Hamiltonian RMS in the E cavity** (`R_inner ≤ ρ ≤ R_stable`) except the last boundary row. Orders are low→mid / mid→high. `t2-class-orders.csv` gives every field and mask.

| Cell class | Old low RMS | Old orders | Pointwise low RMS | Pointwise orders |
| --- | ---: | ---: | ---: | ---: |
| Straight patch edge | 0.777621 | 0.00679 / 0.01013 | 2.14927e-7 | 3.696 / 3.739 |
| Convex corner | 1.70053 | -0.0545 / -0.0376 | 2.51568e-7 | 3.650 / 3.785 |
| Axis × patch edge | 0.758064 | -0.0546 / -0.0381 | 6.72841e-6 | 3.559 / 3.715 |
| Cartoon axis away from interface | 8.55735e-7 | 3.763 / 3.843 | 8.55735e-7 | 3.763 / 3.843 |
| Internal box seam | 1.37890e-7 | 2.841 / 3.664 | 1.37890e-7 | 2.841 / 3.664 |
| Outer boundary, 1 M shell | 1.79808e-14 | roundoff; no order | 1.79808e-14 | roundoff; no order |
| Re-entrant corner | absent | — | absent | — |

The diagnostic-only `t2_exact_initial_ghost_diagnostic = true` rereads EMSTRUMPET **only for a `max_steps = 0`, `t = 0` plot**, computes the connection on scratch boxes, and substitutes only fine ghosts before constraint evaluation. It aborts before any file read at later times or in an evolving run. E low cavity patch-edge Hamiltonian falls `0.777621 → 2.35430e-8`, convex-corner `1.70053 → 5.51065e-9`, and axis-junction `0.758064 → 1.16421e-7`. This isolates interpolation from valid-cell data. The working reference member's exp-0002-like seven-level hierarchy has the same defect: patch-edge Hamiltonian in `2 ≤ ρ ≤ 4 M` is `1.76816e-3` old, `4.43832e-10` with the pointwise experiment, and `5.08968e-11` with exact t = 0 ghosts. All three reference runs have identical boxes; momentum, GaussE, and `[4,8] M` values are in `t2-localization.csv`.

## T2 B. Candidate tests

The rejected pointwise candidate uses a sixth-order tensor-product polynomial of **current coarse fields** and Chombo's RK4 coarse-time polynomial. It obeys the no-static-solution rule. The C++ code is archived only in `t2-pointwise-experiment.patch`, **not applied** to the active build. Its three-grid t = 0 composite constraint orders are measured, not rounded to fourth:

| E mask | Ham low→mid / mid→high | Mom | GaussE |
| --- | ---: | ---: | ---: |
| Horizon | 4.007 / 4.004 | 4.005 / 4.006 | 4.017 / 4.008 |
| Inside inner ring | 3.889 / 3.892 | 3.952 / 3.958 | 3.870 / 3.892 |
| Inner-side cavity | 3.843 / 3.780 | 3.898 / 3.932 | 3.791 / 3.849 |
| Between unstable rings | 3.850 / 3.788 | 3.909 / 3.940 | 3.788 / 3.849 |
| Outer ring | 3.773 / 3.829 | 3.569 / 3.701 | 3.556 / 3.696 |
| Far `[4,8] M` | 3.785 / 3.561 | 3.716 / 3.799 | 3.719 / 3.801 |

The minimum requested three-grid order is **3.556**, so a strict fourth-order gate is unproved. On the extra grid, momentum and GaussE reach about fourth order, but Hamiltonian drops to order 1.22 in the cavity and -0.047 in the far shell. This is consistent with Float64 second-derivative cancellation; the cause of that extra-grid limit was not isolated. Exact norms and the extra-grid orders are in `t2-orders.csv`.

Two coarse E steps reach `t = 0.4375 M` with p192, one rank × four threads. The unfixed first step is T1's run, restarted from its complete step-1 checkpoint for step 2; its step-1 plot data before and after restart are bitwise equal. The pointwise experiment ran continuously for two steps. Patch-edge Hamiltonian RMS is:

| Mask | Unfixed step 1 | Pointwise step 1 | Unfixed step 2 | Pointwise step 2 | Step-2 ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Inside inner ring | 0.0276881 | 1.31938 | 0.00322416 | 1.03974 | 322.48 |
| Cavity | 0.0192615 | 0.413183 | 0.00348899 | 0.348143 | 99.78 |
| Outer ring | 0.000613773 | 0.00507743 | 0.0000860827 | 0.00493605 | 57.34 |
| Far field | 0.000145793 | 0.0000606907 | 0.0000114003 | 0.00000560364 | 0.492 |

`t2-smoke-comparison.csv` gives all constraints and class/mask combinations. The geometric p384 rectangles are the second substantive remedy: they recover about fourth order inside the inner ring by moving that interface, but cavity Hamiltonian orders remain `0.454 / 0.470`, outer-ring `0.388 / 0.424`, and far `0.432 / 0.454`; cavity patch-edge Hamiltonian remains about 0.72 at all three resolutions. Raising `grid_buffer_size` from 8 to 32 produced **identical boxes and norms**, so it did not move the interface. `t2-p384-orders.csv` and `t2-buffer32-orders.csv` give all fields. **Neither remedy passes the combined t = 0 and evolution gate; no production transfer change is applied.**

## T2 C. After-t=0 initial-data guard

`t2_guard_initial_data_after_t0 = true` (default false) makes `EMSBH2DLevel::initialData()` abort before constructing a reader when either the level time or AMR current time is positive. The AMR-time check covers undefined restart levels whose own time was never loaded. Starting from a checkpoint at `t = 0.0625 M`, `finest_level = 0`, with the source file absent, a normal `max_level = 0` restart exits 0 with no reader call; a constructed `max_level = 1` restart aborts with `initialData called after t=0 (incomplete restart hierarchy)` before reading. `t2-guard.csv` records both. The default-off incomplete-restart path retains the T1 problem and must not be used after t = 0.

## T2 D. B level 15

The macOS postmortem of T1's centered `N1 = 512`, p96 SIGSEGV resolves `TreeIntVectSet::clearTree ← TreeIntVectSet::grow ← GRAMRLevel::tagCells(14)`. Chombo's `TreeIntVectSet.cpp:248` and `TreeIntVectSet.H:314` use fixed 24-entry traversal arrays. The puncture's level-14 x index is `2^22`; tags cross that dyadic span edge and require an extra tree depth during `grow(3)`. The fixed-depth tree traversal is the root-cause diagnosis from the source, crash stack and layout perturbations; an instrumented memory trace was unavailable because LLDB could not launch under the sandbox. Shifting the fixed puncture and geometric tagging center **one base cell** to `x = 223.125 M` retains `N1 = 512`, `h0 = 0.875 M`, ratio 2 and all equations; level 15 initializes and the t = 0 run exits 0. Centered `N1 = 384` and `256` at the same h0 also build level 15. `t2-b15.csv` records these tests. No Chombo modification is needed for these layouts.

## T2 default-path identity and artifacts

The final build contains only opt-in guard and t = 0 diagnostic code, both default off. Against committed T1 output, the final default-off E t = 0 plot matches **all 52 datasets and 173 attributes bit for bit**; its two-step checkpoint matches all 4 datasets and 48 attributes. Normal restart with guard on also matches all 4 checkpoint datasets and 48 attributes. Raw HDF5 file bytes differ in container metadata, so raw-file identity is **not** claimed. `t2-bit-identity.csv` records each control. No run log reports a missing `sigma` or missing plot variable. All file reads in the two-step runs precede the first advance; no static solution is used after t = 0.

`t2-*.csv` files hold the layouts, residuals, orders, controls, guard, B15 and smoke results. `params/t2-*.txt` are the run parameters, including those for the rejected patch. Applying `t2-pointwise-experiment.patch` reconstructs that experiment and is **not recommended for evolution**. The active source changes are `Examples/EMS/EMSBH2DLevel.cpp`, `Examples/EMS/SimulationParameters.hpp`, and `Source/GRChomboCore/GRAMR.hpp`. `COMMIT-MANIFEST-T2.txt` hashes every delivered T2 file except itself, plus the read-only inputs. Raw HDF5 and logs are in `/private/tmp/ems-t2-*`, outside the worktree.

# T3: initial restriction and evolved convergence

**Status: BLOCKED.** The controller's proposed restriction mechanism is real, but applying restriction once at initialization does not pass the registered t = 0 constraint gate. The fixed-layout evolved controls below test whether it affects later residuals. All T3 evolution uses the existing equations, gauge, puncture, Float64, `sigma = 1`, fixed boxes, and only current evolved fields after t = 0.

## T3 A. Mechanism and t = 0 test

`GRAMRLevel::postInitialize()` previously only set `m_restart_time = 0`, whereas `postTimeStep()` calls the finer level's `CoarseAverage::averageToCoarse` before filling boundaries. Chombo calls `postInitialize()` finest to coarsest after all levels' `initialData()`. The opt-in `t3_restrict_initial_covered = true` (default false) now applies that same existing restriction once in `EMSBH2DLevel::postInitialize()`, followed by the same boundary fill. It reads no file and uses no static reference. All three E registered t = 0 hierarchies retain exactly the T1 boxes.

`t3-restriction-check.py` compares **every** covered coarse valid cell against the arithmetic 2 × 2 fine average in the HDF5 checkpoints: 64,896 cells on 12 interfaces, all 28 evolved fields. “Near interface” means two coarse layers by either x edge or the upper y edge of the fine rectangle (4,896 cells). All values below are maximum absolute differences, not constraint norms.

| E low state | All covered components | Near-interface components |
| --- | ---: | ---: |
| Unfixed t = 0 setter state | 1.05408 (`A12`, level 11) | 1.09814e-5 (`Ey`, level 9) |
| Unfixed after first coarse step | 8.88178e-16 | 2.22045e-16 |
| Opt-in restricted t = 0 | 1.77636e-15 | 2.22045e-16 |

At the low level-11 interface, the unfixed near-interface differences are 5.21261e-8 in χ, 1.19214e-7 in lapse, and 4.04574e-7 in φ. `t3-mechanism.csv` gives all 1,008 field/level/state checks. This establishes the change in the covered coarse representation and verifies that the opt-in hook actually performs the intended averaging.

Composite cylindrical RMS constraint orders on E's same low/mid/high 3/2 spacing sequence are below (low→mid / mid→high). Masks and light-ring landmarks were fixed before the runs as in T1; none enter evolution.

| E mask | Hamiltonian | Momentum | GaussE |
| --- | ---: | ---: | ---: |
| Horizon | 4.007 / 4.004 | 4.005 / 4.006 | 4.017 / 4.008 |
| Inside inner ring | 0.346 / 0.391 | 1.107 / 1.236 | 1.303 / 1.368 |
| Inner-side cavity | 0.465 / 0.478 | 1.487 / 1.527 | 1.420 / 1.444 |
| Between unstable rings | 0.478 / 0.488 | 1.498 / 1.536 | 1.430 / 1.451 |
| Outer ring | 0.273 / 0.352 | 1.169 / 1.281 | 1.359 / 1.273 |
| Far `[4,8] M` | 0.424 / 0.408 | 1.298 / 1.314 | 1.311 / 1.343 |

The cavity (`R_inner ≤ ρ ≤ R_stable`) classes show the interface failure directly:

| Cell class | Low Hamiltonian RMS | Hamiltonian order | Momentum order | GaussE order |
| --- | ---: | ---: | ---: | ---: |
| Straight patch edge | 0.646157 | 0.00591 / 0.00953 | 1.015 / 1.065 | 0.959 / 0.979 |
| Convex corner | 1.42344 | -0.0552 / -0.0381 | 0.704 / 0.798 | 0.848 / 0.897 |
| Axis × patch edge | 0.599889 | -0.0556 / -0.0389 | 2.264 / 1.845 | 0.844 / 0.891 |
| Cartoon axis away from edge | 0.0474497 | 0.480 / 0.485 | 2.313 / 1.808 | 1.345 / 1.393 |

The old low cavity edge and corner Hamiltonian RMS values were 0.777621 and 1.70053, so the one-time restriction makes a small improvement there; it makes the cartoon-axis-away class much worse than T2's 8.55735e-7. The order ≥ 3.5 acceptance fails in five of six full masks and the key interface classes. `t3-e-t0-localization.csv` and `t3-e-t0-orders.csv` include every measured T2 class, mask, field, cell count and both orders; absent or roundoff-only classes receive no pass claim. GaussB is identically zero at t = 0 and has no order. The remaining low orders after covered cells become exact fine averages show that this initialization mismatch alone cannot explain the t = 0 interface residual.

## T3 B. Evolved controls

The reference sequence was registered before the final runs: `h0 = 2, 4/3, 8/9 M`, ratio 3/2, `max_level = 6`, fixed `regrid_interval = 0`, `dt0 = h0/16`, and one plot per 0.5 M through 5 M. Separate pre-run geometric radii produce **identical physical patch unions** on every resolution: x half width and y height 64, 32, 16, 8, 4, 2 M on levels 1–6; `t3-reference-grid-check.csv` records all 18 verified rectangular extents and box counts. The individual same-level box splits differ. Low/mid/high outputs at 1, 2.5, 5 M are steps 8/12/18, 20/30/45, 40/60/90. Fixed masks are the T2 reference horizon, `[0.1,2]`, `[2,4]`, and `[4,8] M` shells; they are offline diagnostics only. All six final parameters are `params/t3-ref-match-*.txt`; `t3-grid-probes.csv` records exploratory geometry tuning, and the unmatched first low pair is excluded from the convergence fit.

The E low and mid opt-in p192 hierarchies each completed two coarse steps under their 90-minute per-run watchdogs: **2,385.18 s** (4 threads) and **4,686.08 s** (8 threads). The pre-run mid estimate was about 65 minutes in isolation from T1's low smoke cost and the measured low→mid setup ratio; the actual 78.10 minutes includes overlap with four-thread reference runs. Low step 2 is at `t = 0.4375 M`; mid step 2 is at `t = 0.291667 M`, so their values are a class comparison and **not** an observed spatial order. Patch-edge Hamiltonian RMS at step 2 is:

| E mask | T3 low restricted | T2 low unfixed | T2 low pointwise | T3 mid restricted |
| --- | ---: | ---: | ---: | ---: |
| Inside inner ring | 0.00330973 | 0.00322416 | 1.03974 | 0.00461742 |
| Inner-side cavity | 0.00347211 | 0.00348899 | 0.348143 | 0.00340290 |
| Outer ring | 8.31418e-5 | 8.60827e-5 | 0.00493605 | 1.11350e-4 |
| Far `[4,8] M` | 1.10388e-5 | 1.14003e-5 | 5.60364e-6 | 1.04956e-5 |

The one-time restriction stays close to the unfixed low run and does not reproduce the rejected pointwise candidate's interface growth in the inner/cavity/ring masks. The pointwise candidate is smaller in the far mask, as the table shows. `t3-e-smoke-localization.csv` gives all masks, classes, fields and steps; `t3-e-class-comparison.csv` joins the low values to T2. T2 contains no evolved mid-resolution control, so no mid/T2 ratio is assigned.

The four completed reference controls have these measured wall times: low off **365.75 s**, low on **392.05 s**, mid off **902.31 s**, mid on **1,175.59 s**. Each reached 5 M. The table gives observed **low→mid only** orders of Hamiltonian/momentum/GaussE composite cylindrical RMS. “Near hole” is `[0.1,2] M`, “ring” is `[2,4] M`; the other masks are the reference horizon shell and `[4,8] M`. These values alone are not a three-grid convergence claim.

| t/M | Mask | Default off H/M/G | Initial restriction on H/M/G |
| ---: | --- | --- | --- |
| 1 | Horizon | 2.625/3.675/7.876 | 2.625/3.676/7.876 |
| 1 | Near hole | 2.694/1.899/2.029 | 2.695/1.899/2.029 |
| 1 | Ring | 0.730/0.704/1.630 | 0.733/0.707/1.629 |
| 1 | Far | 0.752/0.720/1.667 | 0.759/0.724/1.665 |
| 2.5 | Horizon | 0.728/1.017/4.044 | 0.728/1.026/4.036 |
| 2.5 | Near hole | 2.450/2.095/1.961 | 2.451/2.095/1.961 |
| 2.5 | Ring | 1.199/0.912/1.630 | 1.213/0.923/1.629 |
| 2.5 | Far | 0.686/0.758/1.823 | 0.701/0.761/1.822 |
| 5 | Horizon | 0.920/0.856/6.215 | 0.929/0.864/6.221 |
| 5 | Near hole | 2.882/1.965/1.985 | 2.882/1.965/1.985 |
| 5 | Ring | 1.023/0.921/2.013 | 1.033/0.932/2.015 |
| 5 | Far | 1.170/1.023/1.961 | 1.178/1.031/1.960 |

For the ring's Hamiltonian cell classes, low→mid orders off/on at t = 1, 2.5, 5 M are: patch edge **0.994/0.997, 1.381/1.379, 1.264/1.288**; convex corner **0.982/0.985, 2.890/2.829, 1.329/1.343**; axis × patch edge **1.797/1.801, 1.956/1.941, 1.228/1.186**. `t3-reference-timeseries.csv` has every measured mask/class/field/half-M output; `t3-reference-orders.csv` has the first-pair orders and blank high-grid columns. GaussB remains zero in all 44 completed plots (`t3-reference-gaussb.csv`). There is no material first-pair improvement from the one-time restriction at the coarse–fine interface.

**PENDING: high-resolution reference off/on and all mid→high evolved orders.** The initially attached high-off run was interrupted at about 1.52 M solely to replace it with a detached run; its partial output is kept at `/private/tmp/ems-t3-ref/t3-ref-match-high-off-interrupted` and is excluded from analysis. The detached launcher `t3-run-pending.sh` runs high off then high on sequentially with eight OpenMP threads, stdout and stderr in each run's `run.log`. It was started with Python `subprocess.Popen(..., start_new_session=True, stdin=DEVNULL, close_fds=True)`; launcher PID and log are `/private/tmp/ems-t3-ref/t3-pending-launcher.pid` and `/private/tmp/ems-t3-ref/t3-pending-launcher.log`. The running off simulation has an open log descriptor independently of this agent session.

| Pending run | Command executed by detached launcher (from run directory) | Run directory | Log | Exit-code marker |
| --- | --- | --- | --- | --- |
| High off | `OMP_NUM_THREADS=8 /Users/auroradysis/Workspace/EMS-deps/worktrees/wt-native-t1/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex /Users/auroradysis/Workspace/EMS-deps/worktrees/wt-native-t1/Tests/EMSNative/params/t3-ref-match-high-off.txt` | `/private/tmp/ems-t3-ref/t3-ref-match-high-off` | `run.log` | `done.exit` |
| High on | `OMP_NUM_THREADS=8 /Users/auroradysis/Workspace/EMS-deps/worktrees/wt-native-t1/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex /Users/auroradysis/Workspace/EMS-deps/worktrees/wt-native-t1/Tests/EMSNative/params/t3-ref-match-high-on.txt` | `/private/tmp/ems-t3-ref/t3-ref-match-high-on` | `run.log` | `done.exit` |

Each command redirects stdout and stderr to the stated log; the launcher atomically writes the numeric exit code to the stated marker on return. The queued on run starts after a successful off run. After both markers read `0`, rerun `PYTHONDONTWRITEBYTECODE=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t3-reference-analyze.py /private/tmp/ems-t3-ref`; this fills the high-grid columns of `t3-reference-orders.csv`, updates the time-series/layout/GaussB CSVs, and permits the two-pair orders at 1, 2.5, 5 M and the two high wall times to be reported. Recheck `pout` for missing variables and `sigma`, finish the read audit, and refresh the T3 manifest then. The current T3 status remains **BLOCKED** by the completed t = 0 gate even if the pending evolved orders improve.

## T3 C. Active-source clean-up and identity

The exact-ghost diagnostic was removed from `EMSBH2DLevel::prePlotLevel()` and its parameter removed from `SimulationParameters.hpp`. `t2-exact-ghost-diagnostic.patch` archives precisely that T2 experiment and passes `git apply --check` against the active source. Its recorded T2 measurements remain historical evidence. The active production EMSTRUMPET reader construction is now only inside `initialData()`; the opt-in after-t = 0 guard is retained. All T3 evolution parameter files enable that guard.

With both T3 flags absent/default off, the final build matches T2's E low t = 0 plot in all 52 datasets and 173 attributes, and its two-step complete-hierarchy checkpoint in all 4 datasets and 48 attributes, bit for bit (`t3-bit-identity.csv`). The raw HDF5 container bytes are not asserted equal. No source equation, gauge, RHFinder, Chombo library or radial coordinate was changed.

# T4: point transfer with RK4 stage ghosts

**Status: PARTIAL evolved convergence; the registered t = 0 acceptance gate remains BLOCKED.** The operator and wave checks pass, and default-off / restart controls are bitwise equal. The all-mask/all-class order ≥ 3.5 gate does not pass: 24 measured mask/class/constraint combinations fail at least one pair. All three detached reference runs to 2.5 M have completed with exit 0; the completed legacy/point comparison is in T4c below.

## T4 A. Scheme and scope

`amr_transfer = point` is the single new opt-in switch; absent or `legacy` retains the original paths. Ratio 2, DIM=2, three fine ghosts and Float64 are explicitly checked. Every delivered physics parameter pins `sigma = 1`. No equation, gauge, floor, RHFinder, coordinate, Chombo library or read-only EMS input is modified. No commit was made.

`Source/GRChomboCore/PointAMRTransfer.hpp` implements six-point tensor Lagrange prolongation, degree five in each coordinate, O(h^6). Ordinary evolution ghosts, diagnostic ghosts and regrid valid/boundary cells use the same spatial operator. At reflective boundaries the stencil remains centred through parity-filled coarse ghosts; four coarse buffer layers cover the three fine ghosts plus the six-point stencil. Same-level exchange wins at seams and periodic images. Point restriction uses the six-point midpoint weights `(3,-25,150,150,-25,3)/256` in each coordinate, with current fine ghosts filled before covered coarse cells are replaced. Accumulation anchors on a current source cell to preserve constants exactly. Initial valid point samples are retained; restriction occurs at synchronization, not as a second initial-data prescription.

The coarse RK4 initial value and four RHSs are saved into Chombo's existing dense polynomial. Each fine RHS uses the *start fraction of its fine step* and a distinct stage index 0–3; the two half-time stages are deliberately different. A local `PointRK4Interpolator` wrapper corrects Chombo's two `diff12` coefficients from r^2 to r^3, where r=dt_fine/dt_coarse. The exact generic vector-ODE jet witness gives old stage residuals `-H^3 F'^2 F/128` and `+H^3 F'^2 F/64` at ratio 2; the corrected residuals through H^3 are zero. The library file is untouched. Sources: [Chombo implementation](https://github.com/GRTLCollaboration/Chombo/blob/main/lib/src/AMRTimeDependent/TimeInterpolatorRK4.cpp) and [McCorquodale–Colella 2011, equations 43–49](https://msp.org/camcos/2011/6-1/camcos-v6-n1-p01-p.pdf); the local source/PDF text was inspected. Context7 did not index Chombo or GRChombo, so the pinned local source supplied the API contract.

`t3_restrict_initial_covered` is removed from the active parser and hook. Historical T3 parameter files and measurements retain their original provenance; the retired parameter has no active implementation. In point mode, `initialData()` aborts before constructing any reader after t=0, including an undefined restart level. The static data reader exists only in `initialData()`; transfer, RK RHS, regrid and synchronization take current arrays and grid/boundary metadata only. All masks below are the requested registered offline diagnostic shells; none enters the evolution, tagging, transfer or damping. `t4-static-read-audit.csv` and the absent-file restart controls record this scope.

Builds: `CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib make -C Examples/EMS all DIM=2 -j4`, and the same command with `-C Tests/EMSNative` for `PointTransferTest`. The pristine cd92a4b executable is `/private/tmp/ems-t4/baseline.ex`; the current executable is in `Examples/EMS`. The runtime queues are serial; all measured run RSS bounds are below 12 GB. Peak RSS is recorded with Python `getrusage`, because the sandbox denies the sysctl used by `/usr/bin/time -l`. Within a queue the recorded RSS is conservative: it retains the largest earlier child's peak.

## T4 B. Operator checks

`PointTransferTest.cpp` exercises the compiled production operators, not a second Python implementation. It tests all 36 tensor monomials x^a y^b, a,b≤5, at N=32,64,128, with even/odd axis parity and two fine boxes. Direct analytic fine samples independently test restriction. Maximum polynomial errors are 1.3877787807814457e-17 for prolongation and 3.4694469519536142e-18 for restriction. The non-polynomial field is exp(2x) cos(3y), or exp(2x) sin(3y) for odd parity. The following orders use max errors at 32→64 / 64→128.

| Operator / parity | Whole tested patch | Axis, first two fine/coarse rows |
| --- | ---: | ---: |
| Prolongation, even | 6.056 / 6.028 | 6.056 / 6.028 |
| Restriction, even | 5.976 / 5.990 | 5.976 / 5.990 |
| Prolongation, odd | 6.056 / 6.029 | 7.054 / 7.027 |
| Restriction, odd | 5.973 / 5.988 | 6.973 / 7.010 |

The compiled RK stage wrapper is compared against independently constructed fine RK stages for u'=0.7u, fine-step start theta=0,1/2, H=1/8,1/16,1/32. Its largest stage errors are 1.0732744670782779e-6, 6.6931720610341472e-8, 4.1786385462927456e-9, with residual orders 4.003/4.002. `t4-operator-checks.csv` holds every error. The exact moment/vector-jet proof is rerunnable in `scripts/cas/t4-transfer-verify.py` with its JSON transcript and evidence card. CAS scope: algebraic identity — production physics not certified.

## T4 C. Registered t = 0 constraints

E uses the unchanged T1 3/2 sequence h0=0.875,0.583333333333,0.388888888889 M, max_level=12, fixed geometric tags, no regrid. Reference uses the T3 matched exp-0002-like h0=2,4/3,8/9 M, max_level=6. Reference physical patch half widths/heights remain 64,32,16,8,4,2 M; same-level splits differ. `t4-t0-layout.csv` records all 60 actual levels. Norms and classes reuse T2's classifier, covered coarse removal and cylindrical RMS weighting. The T1 core shells are included: cavity_core [0.055,0.075], outer_ring_core [0.85,1.0], far_core [6.5,8] M. `t2-localize.py` gained only an optional extra-mask argument, preserving existing callers.

Every full physical mask passes ≥3.5 at both pairs; this does not pass the required cell-class gate. Entries are low→mid / mid→high.

| Member / mask | Class | Ham | Mom | GaussE |
| --- | --- | ---: | ---: | ---: |
| E / horizon | all | 4.007/4.004 | 4.005/4.006 | 4.017/4.008 |
| E / inside_inner_ring | all | 3.886/3.890 | 3.952/3.958 | 3.870/3.892 |
| E / cavity | all | 3.834/3.776 | 3.898/3.932 | 3.791/3.849 |
| E / between_rings | all | 3.840/3.784 | 3.909/3.940 | 3.788/3.849 |
| E / outer_ring | all | 3.761/3.825 | 3.568/3.701 | 3.556/3.696 |
| E / far | all | 3.746/3.543 | 3.716/3.799 | 3.719/3.801 |
| E / cavity_core | all | 3.948/3.843 | 3.941/3.966 | 3.962/3.979 |
| E / outer_ring_core | all | 3.908/3.940 | 3.829/3.914 | 3.866/3.932 |
| E / far_core | all | 3.875/3.600 | 3.882/3.910 | 3.896/3.914 |

| Member / mask | Class | Ham | Mom | GaussE |
| --- | --- | ---: | ---: | ---: |
| ref / horizon | all | 4.013/4.011 | 4.033/4.012 | 4.032/4.023 |
| ref / near_hole | all | 3.898/4.096 | 3.923/3.999 | 5.863/6.025 |
| ref / ring | all | 4.053/4.045 | 3.989/3.995 | 3.997/3.999 |
| ref / far | all | 4.036/4.024 | 4.013/4.006 | 4.001/4.001 |

The requested interface classes show the remaining failures explicitly. “Absent” is a cell census result, not a convergence pass.

| Member / mask | Class | Ham | Mom | GaussE |
| --- | --- | ---: | ---: | ---: |
| E / inside_inner_ring | patch_edge | 3.716/3.806 | 3.730/3.832 | 3.881/3.932 |
| E / inside_inner_ring | axis_patch_junction | 3.602/3.437 | 3.632/3.763 | 3.790/3.876 |
| E / cavity | patch_edge | 3.686/3.732 | 4.123/4.090 | 3.849/3.900 |
| E / cavity | convex_corner | 3.652/3.775 | 3.638/3.741 | 4.121/4.120 |
| E / cavity | axis_patch_junction | 3.575/3.894 | 3.963/3.952 | 3.811/3.881 |
| E / between_rings | patch_edge | 3.729/3.762 | 4.152/4.102 | 3.877/3.912 |
| E / between_rings | convex_corner | 3.666/3.796 | 3.638/3.741 | 3.999/4.010 |
| E / between_rings | axis_patch_junction | 3.601/3.865 | 3.956/3.947 | 3.809/3.878 |
| E / outer_ring | patch_edge | 3.786/3.844 | 3.757/3.832 | 3.799/3.864 |
| E / outer_ring | axis_patch_junction | 3.711/3.781 | 3.525/3.672 | 3.691/3.800 |
| E / far | patch_edge | 3.827/2.489 | 3.726/3.816 | 3.734/3.821 |
| E / far | convex_corner | 2.551/2.937 | 5.727/3.793 | 5.495/3.812 |
| E / far | axis_patch_junction | 4.114/3.323 | 3.730/3.818 | 3.737/3.823 |
| E / cavity_core | patch_edge | 3.765/3.768 | 3.800/3.874 | 4.100/4.080 |
| E / cavity_core | convex_corner | 3.695/3.716 | 3.703/3.810 | 4.070/4.071 |
| E / outer_ring_core | patch_edge | 3.743/3.785 | 4.129/4.060 | 4.029/4.008 |
| E / outer_ring_core | convex_corner | 2.644/3.608 | 3.695/3.877 | 3.938/3.929 |
| E / far_core | patch_edge | 3.759/2.493 | 3.819/3.873 | 3.670/3.814 |
| E / far_core | convex_corner | 2.320/2.937 | 3.601/3.793 | 3.634/3.812 |

| Member / mask | Class | Ham | Mom | GaussE |
| --- | --- | ---: | ---: | ---: |
| ref / near_hole | patch_edge | 3.910/3.972 | 4.175/4.135 | 4.056/4.046 |
| ref / near_hole | axis_patch_junction | 3.952/3.971 | 4.244/4.181 | 4.146/4.106 |
| ref / ring | patch_edge | 3.811/3.818 | 4.087/4.034 | 3.961/3.924 |
| ref / ring | convex_corner | 4.101/4.065 | 5.520/5.311 | 4.393/4.269 |
| ref / ring | axis_patch_junction | 3.854/3.909 | 4.160/4.082 | 4.272/4.202 |
| ref / far | patch_edge | 3.964/3.917 | 3.979/3.941 | 4.119/4.055 |
| ref / far | convex_corner | 3.907/3.942 | 4.136/4.108 | 4.065/4.058 |
| ref / far | axis_patch_junction | 4.093/4.057 | 4.232/4.153 | 4.164/4.106 |

All 252 mask/class/constraint cases are in `t4-t0-orders.csv`: 197 pass, 24 fail, 21 have a class absent on at least one grid, and 10 have max RMS <1e-13 (reported as roundoff, no pass). This roundoff label was fixed in the analyzer and does not erase the 24 failures. All measured same-level seam, interior, axis-away and outer-boundary cases are retained in the CSV, including the failures. Re-entrant corners are absent from these rectangular unions. GaussB is identically zero and all saved plot values are finite (`t4-t0-field-checks.csv`). A cause for the remaining low class orders is not established; do not reinterpret the near-zero far-field Hamiltonian values or changing class membership as a passed gate.

## T4 D. Time-dependent two-level test

The test-only WaveLevel evolves u_t=v, v_t=∂xx u+∂yy u with the code's native fourth-order derivatives and native Kreiss–Oliger dissipation, sigma=1 on both fields. Initial data are a right-moving Gaussian exp(-((x-1.5)/0.4)^2), periodic in x, independent of y with reflective y boundaries. Domain [0,8]×[0,4], refinement strip [2,6]×[0,4], two fine boxes, ratio 2. The pulse crosses x=2 before t=0.75. No analytic value is used by the RHS or transfer after initialization.

N1=64,128,256, dt0=h0/8, 48/96/192 coarse steps to t=0.75. Fine steps are half the coarse dt. Independent six-point interpolation on the composite nonuniform output nodes evaluates 481 fixed x points in [1,4]. Self-convergence is log2(||U64-U128||/||U128-U256||); it is one three-grid order, not two separate constraint-to-zero orders.

| Scheme | Field | Low–mid RMS | Mid–high RMS | RMS order | Max order |
| --- | --- | ---: | ---: | ---: | ---: |
| legacy | u | 0.00242473734 | 0.000446288819 | 2.442 | 2.595 |
| legacy | v | 0.0116612692 | 0.00192840359 | 2.596 | 2.953 |
| point | u | 0.00197328007 | 0.0000931218346 | 4.405 | 4.306 |
| point | v | 0.0158387419 | 0.000785863332 | 4.333 | 4.333 |

`t4-wave-orders.csv` holds unrounded values. The final legacy runs took 3.998,12.981,47.489 s; point runs took 0.930,2.245,8.278 s, four threads. These include output and setup and are not a benchmark claim. Preliminary no-dissipation outputs are archived under `/private/tmp/ems-t4/wave-no-ko` and are excluded from this final sigma=1 qualification.

## T4 E. Bit identity and restart

The baseline executable was built directly from detached HEAD cd92a4b before any production edit. Comparisons check the bits of all numeric HDF5 values and attributes (compound-field padding and HDF5 container metadata are outside the claim). All ten comparisons pass; `t4-bit-identity.csv` records exact paths and counts.

| Control | E datasets / attributes | Reference datasets / attributes | Result |
| --- | ---: | ---: | --- |
| Registered low t=0 plot, default vs cd92a4b | 52 / 158 | 28 / 86 | bitwise equal |
| Small two-level t=0 plot | 8 / 26 | 8 / 26 | bitwise equal |
| Small two-level checkpoint after 2 steps | 8 / 26 | 8 / 26 | bitwise equal |
| Same checkpoint, regrid_interval=1 | 8 / 26 | 8 / 26 | bitwise equal; regrid callback exercised |
| Point continuous 4 steps vs restart at step 2 | 8 / 26 | 8 / 26 | bitwise equal; source path absent |

The small controls use domain [−4,4]×[0,4] M, N1=64,N2=32, max_level=1, one thread. Coarse dx=0.125 M; E dt0=0.03125 M, reference dt0=0.0078125 M. The 2-step checkpoint equality is scoped to these complete two-level hierarchies, not an E level-12 evolution. Point regrid smoke runs also exit 0. Default identity is not claimed for the removed T3-on path. No post-t=0 static-profile read occurs in the complete point restarts.

## T4 F. Detached reference evolution to 2.5 M

The registered T3 match inputs are copied to `params/t4-ref-{low,mid,high}-2p5M.txt`, with `amr_transfer=point`, sigma=1, max_steps=20/30/45, checkpoint at the endpoint, original plot cadence 4/6/9 (0.5 M), no regrid or run-time profile-driven placement. The time sequence is unchanged: dt0=h0/16 and each level subcycles by 2. One detached serial queue executes low then mid then high. Per-job allocation bound 4 GB is below the 12 GB gate; the completed t=0 queue's conservative peak was 1.0361 GiB. The launcher uses a new session, disconnected stdin, separate stdout/stderr log and atomic done.exit per run. It stops on a nonzero run; later markers then remain absent.

Launch command: `python3 Tests/EMSNative/t4-run.py --detach /private/tmp/ems-t4/evolution/plan.json`. Launcher PID and log are `/private/tmp/ems-t4/evolution/launcher.pid` and `/private/tmp/ems-t4/evolution/launcher.log`. Commands below run with OMP_NUM_THREADS=4 from their stated directories; the executable is the frozen snapshot `/private/tmp/ems-t4/point.ex` of this worktree's final EMS build, hashed in `t4-builds.csv`.

| Run (at launch) | Parameter path (tree relative) | Directory | stdout/stderr | Atomic marker |
| --- | --- | --- | --- | --- |
| ref-low | Tests/EMSNative/params/t4-ref-low-2p5M.txt | `/private/tmp/ems-t4/evolution/ref-low` | `/private/tmp/ems-t4/evolution/ref-low/run.log` | `/private/tmp/ems-t4/evolution/ref-low/done.exit` |
| ref-mid | Tests/EMSNative/params/t4-ref-mid-2p5M.txt | `/private/tmp/ems-t4/evolution/ref-mid` | `/private/tmp/ems-t4/evolution/ref-mid/run.log` | `/private/tmp/ems-t4/evolution/ref-mid/done.exit` |
| ref-high | Tests/EMSNative/params/t4-ref-high-2p5M.txt | `/private/tmp/ems-t4/evolution/ref-high` | `/private/tmp/ems-t4/evolution/ref-high/run.log` | `/private/tmp/ems-t4/evolution/ref-high/done.exit` |

The full absolute commands, memory bounds and remaining statuses are in `t4-runs.csv` and the detached plan JSON. This turn ends after launch; the controller watches the markers.

The launch-time remaining analysis is now completed in T4c below: all markers are 0, both times and both resolution pairs are measured, logs/floors/RHFinder availability are audited, and the manifest is refreshed. The T4c common-support analyzer replaces the launch-time proposed independent-class analysis; the original t=0 acceptance gate remains BLOCKED.

`COMMIT-MANIFEST-T4.txt` hashes the changed production source, reused analyzer change, every new T4 test/parameter/CSV, CAS witness/transcript/evidence and this README. Raw outputs, the baseline binary, logs, plan JSONs and markers remain under `/private/tmp/ems-t4`; local make products are excluded from the manifest. The read-only E/reference input hashes are in `t4-inputs.csv`.

# T4b: E failing-class census and Float64 sensitivity

**Result: MIXED.** The finest-grid Hamiltonian residuals in these 19 failed cases are compatible with the measured-field initialization sensitivity envelope; the momentum and electric-Gauss seam residuals remain above it. A floor-compatible envelope is not proof of floor saturation. The low/mid edge/corner signal has a real truncation component localized to coarse–fine **chi ghost prolongation**, while the seam orders compare substantially different physical cell samples. The registered T4 gate remains BLOCKED; no pass is inferred from this diagnosis.

The quoted minima are Ham −2.629 (box_seam/far_core), Mom 0.391 (box_seam/outer_ring_core), and GaussE 2.695 (box_seam/outer_ring). The three minima are different constraints, not three Hamiltonian orders. This section covers every FAIL row for E in box_seam, patch_edge and convex_corner, including the added T1 core masks.

All RMS values below come from the existing `/private/tmp/ems-t4/t0/E-{low,mid,high}/plt/EMS_Plot_000000.2d.hdf5` outputs, with T2's cylindrical weights, registered offline shells and covered-coarse removal. They have only ten saved variables. To estimate sensitivity without reconstructing a solution, three short **max_steps=0, stop_time=0** snapshots of the same frozen `point.ex` saved all 28 current evolution fields and the same diagnostics, with three ghost layers. Every original saved valid value is bitwise identical. The reader only initialized t=0. No evolution, static diagnostic ghost substitution, production edit, or reference-run inspection was performed.

## Raw failures and spacings

Triples are low / mid / high. Units are M for spacing. `h_min` is the finest spacing actually represented in that class/mask, not an unused innermost level. The whole hierarchy's level-12 spacings are 0.000213623046875 / 0.00014241536458333334 / 0.00009494357638888889 M. `t4b-roundoff.csv` also gives h_max, exact full-precision RMS, radial range, feature intersections and operation probes.

| Class | Mask | Constraint | Raw RMS, low / mid / high | Cells, low / mid / high | h_min, low / mid / high |
| --- | --- | --- | --- | --- | --- |
| box_seam | between_rings | GaussE | 3.701133e-08 / 4.193214e-09 / 1.159257e-09 | 964 / 8502 / 14784 | 4.272461e-04 / 2.848307e-04 / 1.898872e-04 |
| box_seam | between_rings | Ham | 5.025867e-08 / 1.572869e-08 / 3.673741e-09 | 964 / 8502 / 14784 | 4.272461e-04 / 2.848307e-04 / 1.898872e-04 |
| box_seam | cavity | GaussE | 9.120828e-08 / 1.064261e-08 / 2.875782e-09 | 424 / 3618 / 6332 | 4.272461e-04 / 2.848307e-04 / 1.898872e-04 |
| box_seam | cavity | Ham | 1.378898e-07 / 4.357875e-08 / 9.862898e-09 | 424 / 3618 / 6332 | 4.272461e-04 / 2.848307e-04 / 1.898872e-04 |
| box_seam | cavity_core | Ham | 1.571321e-08 / 1.747339e-08 / 3.491047e-09 | 96 / 737 / 1168 | 8.544922e-04 / 5.696615e-04 / 3.797743e-04 |
| box_seam | far | GaussE | 3.329830e-10 / 4.164935e-11 / 1.202889e-11 | 200 / 1760 / 3056 | 5.468750e-02 / 3.645833e-02 / 2.430556e-02 |
| box_seam | far | Ham | 2.071924e-10 / 1.512836e-10 / 3.856196e-11 | 200 / 1760 / 3056 | 5.468750e-02 / 3.645833e-02 / 2.430556e-02 |
| box_seam | far | Mom | 8.884828e-13 / 1.155746e-13 / 3.398316e-14 | 200 / 1760 / 3056 | 5.468750e-02 / 3.645833e-02 / 2.430556e-02 |
| box_seam | far_core | Ham | 2.676280e-11 / 7.770931e-11 / 9.230340e-12 | 56 / 292 / 560 | 1.093750e-01 / 7.291667e-02 / 4.861111e-02 |
| box_seam | inside_inner_ring | GaussE | 5.784140e-07 / 6.123720e-08 / 1.524239e-08 | 228 / 2245 / 3716 | 2.136230e-04 / 1.424154e-04 / 9.494358e-05 |
| box_seam | outer_ring | GaussE | 5.170542e-09 / 6.748822e-10 / 2.262721e-10 | 116 / 1076 / 1980 | 6.835938e-03 / 4.557292e-03 / 3.038194e-03 |
| box_seam | outer_ring | Ham | 2.756318e-08 / 7.084270e-09 / 1.550870e-09 | 116 / 1076 / 1980 | 6.835938e-03 / 4.557292e-03 / 3.038194e-03 |
| box_seam | outer_ring | Mom | 5.398567e-11 / 1.207649e-11 / 5.011810e-12 | 116 / 1076 / 1980 | 6.835938e-03 / 4.557292e-03 / 3.038194e-03 |
| box_seam | outer_ring_core | Mom | 8.217183e-12 / 7.012644e-12 / 1.404427e-12 | 44 / 241 / 448 | 1.367188e-02 / 9.114583e-03 / 6.076389e-03 |
| convex_corner | far | Ham | 1.766813e-10 / 6.281243e-11 / 1.909186e-11 | 8 / 8 / 8 | 2.734375e-02 / 3.645833e-02 / 2.430556e-02 |
| convex_corner | far_core | Ham | 1.608889e-10 / 6.281243e-11 / 1.909186e-11 | 6 / 8 / 8 | 5.468750e-02 / 3.645833e-02 / 2.430556e-02 |
| convex_corner | outer_ring_core | Ham | 1.618437e-08 / 5.540547e-09 / 1.283139e-09 | 6 / 8 / 8 | 6.835938e-03 / 4.557292e-03 / 3.038194e-03 |
| patch_edge | far | Ham | 2.255775e-10 / 4.780205e-11 / 1.742588e-11 | 808 / 1192 / 1768 | 5.468750e-02 / 3.645833e-02 / 2.430556e-02 |
| patch_edge | far_core | Ham | 2.343003e-10 / 5.102712e-11 / 1.856720e-11 | 344 / 440 / 596 | 5.468750e-02 / 3.645833e-02 / 2.430556e-02 |

## Constraint-specific Float64 estimates

epsilon = 2.220446049250313e-16. `B_cell` propagates epsilon times the magnitude of every saved stencil input and arithmetic result through the native Hamiltonian/momentum equations, the native fourth-order derivative weights, metric inverse and cartoon terms. It includes the native GammaCartoonCalculator sensitivity before differentiating the stored Gamma, so it includes both differentiation operations. The electric-Gauss contractions are the unchanged source expression, evaluated on those same saved fields. For a second derivative this contains epsilon times the absolute weighted field sum divided by h^2; `eps_chi_over_h2_rms` records the simpler epsilon*abs(chi)/h^2 proxy separately. No static field, reference or profile derivative supplies a scale.

`B_initial` adds a Float64 initialization sensitivity that a one-ulp field estimate omits. `EMSTrumpetSolution_read.impl.hpp:285` accepts its log-coordinate root residual below **100 epsilon**. A residual in that equation is a perturbation in log(R/ell). The extra current-field input scale is therefore `100*epsilon*abs(x*D_x U + y*D_y U)`, with slopes measured by native differences on the saved fields; one epsilon of local coordinate sensitivity is added. This is a first-order local estimate and uses the central slope as the neighbouring-stencil envelope. Native GCC 16 -O3 emits `fnmsub` for `(i+.5)*h-center` on this ARM build (`/private/tmp/ems-t4/t4b/coordinate.s`), so the estimate uses the local coordinate result rather than an unfused 224 M subtraction bound. No compactification or coordinate change is applied to evolution.

These are conservative L1 operation/sensitivity envelopes, **not certified intervals or stochastic RMS noise measurements**. RMS below B_initial is marked compatible with the Float64/100-epsilon initialization floor, not proven saturated. Above it is a resolved component under this estimate. The eight-point probe and current-metric term separation below distinguish resolved interpolation effects from mere envelope compatibility. `scripts/cas/t4b-verify.py` verifies the log-residual identity and P8 moments exactly, with independent Float64 checks; its scope does not certify the physical floor.

| Class | Mask | Constraint | B_cell, low / mid / high | B_initial, low / mid / high | High RMS / B_initial | High assessment |
| --- | --- | --- | --- | --- | ---: | --- |
| box_seam | between_rings | GaussE | 8.648e-13 / 1.678e-12 / 2.379e-12 | 9.306e-12 / 1.836e-11 / 2.601e-11 | 44.564 | above floor |
| box_seam | between_rings | Ham | 4.283e-10 / 1.085e-09 / 2.217e-09 | 5.644e-09 / 1.431e-08 / 2.926e-08 | 0.126 | floor-compatible |
| box_seam | cavity | GaussE | 2.437e-12 / 4.938e-12 / 6.897e-12 | 2.583e-11 / 5.314e-11 / 7.424e-11 | 38.737 | above floor |
| box_seam | cavity | Ham | 1.084e-09 / 2.861e-09 / 5.763e-09 | 1.384e-08 / 3.651e-08 / 7.361e-08 | 0.134 | floor-compatible |
| box_seam | cavity_core | Ham | 8.187e-10 / 1.866e-09 / 4.160e-09 | 1.044e-08 / 2.383e-08 / 5.312e-08 | 0.066 | floor-compatible |
| box_seam | far | GaussE | 5.088e-16 / 1.066e-15 / 1.460e-15 | 7.010e-15 / 1.490e-14 / 2.043e-14 | 588.655 | above floor |
| box_seam | far | Ham | 6.309e-12 / 1.661e-11 / 3.350e-11 | 3.295e-11 / 8.521e-11 / 1.721e-10 | 0.224 | floor-compatible |
| box_seam | far | Mom | 1.260e-18 / 2.259e-18 / 3.076e-18 | 1.645e-17 / 3.214e-17 / 4.387e-17 | 774.567 | above floor |
| box_seam | far_core | Ham | 2.436e-12 / 5.503e-12 / 1.233e-11 | 9.874e-12 / 2.202e-11 / 4.971e-11 | 0.186 | floor-compatible |
| box_seam | inside_inner_ring | GaussE | 1.441e-11 / 2.944e-11 / 4.136e-11 | 1.392e-10 / 2.858e-10 / 4.020e-10 | 37.919 | above floor |
| box_seam | outer_ring | GaussE | 3.135e-14 / 6.736e-14 / 8.969e-14 | 6.114e-13 / 1.332e-12 / 1.777e-12 | 127.348 | above floor |
| box_seam | outer_ring | Ham | 1.035e-10 / 2.906e-10 / 5.588e-10 | 1.462e-09 / 4.072e-09 / 7.838e-09 | 0.198 | floor-compatible |
| box_seam | outer_ring | Mom | 6.392e-16 / 1.201e-15 / 1.592e-15 | 9.337e-15 / 1.883e-14 / 2.500e-14 | 200.467 | above floor |
| box_seam | outer_ring_core | Mom | 2.272e-16 / 3.486e-16 / 5.297e-16 | 3.162e-15 / 5.194e-15 / 8.000e-15 | 175.564 | above floor |
| convex_corner | far | Ham | 1.534e-11 / 2.223e-11 / 4.981e-11 | 8.125e-11 / 8.595e-11 / 1.946e-10 | 0.098 | floor-compatible |
| convex_corner | far_core | Ham | 9.939e-12 / 2.223e-11 / 4.981e-11 | 3.789e-11 / 8.595e-11 / 1.946e-10 | 0.098 | floor-compatible |
| convex_corner | outer_ring_core | Ham | 2.132e-10 / 4.697e-10 / 1.041e-09 | 2.670e-09 / 5.929e-09 / 1.321e-08 | 0.097 | floor-compatible |
| patch_edge | far | Ham | 9.490e-12 / 2.120e-11 / 4.746e-11 | 4.111e-11 / 9.320e-11 / 2.108e-10 | 0.083 | floor-compatible |
| patch_edge | far_core | Ham | 9.712e-12 / 2.177e-11 / 4.886e-11 | 3.956e-11 / 8.927e-11 / 2.014e-10 | 0.092 | floor-compatible |

For the finer grids specifically: all 11 targeted Ham rows are compatible on **high** (RMS/B_initial = 0.0657–0.2240). On **mid**, cavity_core seam and all targeted patch/corner rows are compatible; between_rings, cavity, far, far_core and outer_ring seam Ham remain 1.10, 1.19, 1.78, 3.53 and 1.74 times their estimates. Every targeted Mom/GaussE row remains above the estimate on both finer grids; even on high the smallest ratio is 37.9. The earlier universal max-RMS<1e-13 label cannot diagnose this mixture: far seam Mom at 3.3983e-14 is still 775 times its initial sensitivity scale.

## Seam membership and feature coincidences

The seam class excludes the first two axis rows, the current level's outer patch edge/corner and covered coarse cells. Every targeted seam row has **zero covered cells and zero axis-class cells**. A separate nested Gamma stencil may touch axis parity ghosts in rows 2–3; those counts are retained below. “Near finer” means an uncovered coarse cell within two coarse cells of the next finer valid union, so its stencil can reach covered coarse cells. It is not a fine ghost fill on the cell being measured. Thus a seam can coincide with a *next* level's boundary even though it excludes its *own* level's boundary.

| Seam mask | All cells, low / mid / high | Near next-finer boundary, low / mid / high | Nested stencil reaches axis halo, low / mid / high |
| --- | --- | --- | --- |
| between_rings | 964 / 8502 / 14784 | 40 / 160 / 4160 | 0 / 38 / 40 |
| cavity | 424 / 3618 / 6332 | 16 / 64 / 1960 | 0 / 16 / 16 |
| cavity_core | 96 / 737 / 1168 | 0 / 24 / 180 | 0 / 0 / 0 |
| far | 200 / 1760 / 3056 | 8 / 32 / 896 | 0 / 8 / 8 |
| far_core | 56 / 292 / 560 | 0 / 24 / 332 | 0 / 0 / 0 |
| inside_inner_ring | 228 / 2245 / 3716 | 8 / 8 / 668 | 0 / 8 / 8 |
| outer_ring | 116 / 1076 / 1980 | 8 / 8 / 732 | 0 / 8 / 8 |
| outer_ring_core | 44 / 241 / 448 | 0 / 24 / 244 | 0 / 0 / 0 |

Every refined E union has 2 boxes on low and 8 on mid/high. At level 4 the physical patch half widths are 5.6875 / 5.541666666667 / 5.444444444444 M. Low has only the central x=0 seam. Mid has x seams at −2.916666666667, 0, +2.625 and y=2.625 M. High has x seams at −2.722222222222, 0, +2.722222222222 and y=2.722222222222 M. High's off-centre seams coincide exactly with the next finer union's x boundaries and upper-y boundary. Mid's asymmetric seams do not. The same pattern scales on the other refined levels; `t4b-layout.csv` records all 39 levels. Both the box splitting and snapped physical patch extent change, so the aggregate class does compare different locations, levels, angles and weights. Far_core corner membership is 6 / 8 / 8 cells rather than a fixed stencil sample.

As a diagnostic only, restrict the existing seam census to the central strip shared by all three grids, abs(x)<2h. This does not redefine the registered gate. Orders use the original log(3/2) denominator.

| Mask | Constraint | Original seam orders | Common-central seam orders | Common-central counts |
| --- | --- | --- | --- | --- |
| between_rings | GaussE | 5.371 / 3.171 | 3.759 / 3.847 | 964 / 1428 / 2124 |
| between_rings | Ham | 2.865 / 3.587 | 3.805 / 3.861 | 964 / 1428 / 2124 |
| cavity | GaussE | 5.298 / 3.227 | 3.766 / 3.849 | 424 / 628 / 932 |
| cavity | Ham | 2.841 / 3.664 | 3.827 / 3.881 | 424 / 628 / 932 |
| cavity_core | Ham | -0.262 / 3.972 | 4.206 / 1.761 | 96 / 140 / 208 |
| far | GaussE | 5.127 / 3.063 | 3.726 / 3.797 | 200 / 296 / 440 |
| far | Ham | 0.776 / 3.371 | 3.863 / 3.258 | 200 / 296 / 440 |
| far | Mom | 5.030 / 3.019 | 3.719 / 3.791 | 200 / 296 / 440 |
| far_core | Ham | -2.629 / 5.254 | 3.939 / 0.830 | 56 / 84 / 124 |
| inside_inner_ring | GaussE | 5.538 / 3.430 | 3.911 / 3.868 | 228 / 332 / 496 |
| outer_ring | GaussE | 5.022 / 2.695 | 3.510 / 3.718 | 116 / 164 / 248 |
| outer_ring | Ham | 3.351 / 3.746 | 3.716 / 3.849 | 116 / 164 / 248 |
| outer_ring | Mom | 3.693 / 2.169 | 3.429 / 3.593 | 116 / 164 / 248 |
| outer_ring_core | Mom | 0.391 / 3.966 | 4.128 / 4.061 | 44 / 68 / 100 |

All failing electric-Gauss seam cases regain ≥3.5 at both pairs on the common central seam. The outer_ring_core momentum minimum 0.391 becomes 4.128 / 4.061. The original failed outer_ring momentum mid→high pair becomes 3.593; the central low→mid value is 3.429, so this subset is not a new all-pair qualification. Negative cavity_core/far_core Ham low→mid orders become 4.206 / 3.939, while their finer central Ham orders deteriorate as the initialization sensitivity becomes comparable. The common-strip result and extra-seam census localize the anomalous aggregate orders to class membership and spatially varying native truncation/sensitivity, not to a bad same-level copy.

## Operation localization from current fields

1. **Same-level exchange:** every saved ghost overlap with another valid owner, including covered coarse owners, is identical in all 28 evolution fields. Numeric comparisons are 427728 / 2159136 / 3874752 values with **0 / 0 / 0** bit differences. The original ten valid plot fields compare 3906560 / 8494080 / 18677760 values with **0 / 0 / 0** bit differences (`t4b-identity.csv`). Replaying the native constraints differs from the plotted value by at most 0.02371 B_cell in RMS across these cases (scalar-expression/FMA ordering); this is inside the arithmetic envelope.
2. **Seam source:** its native derivative stencil uses same-level current valid values. The offline ghost probe changes *only* coarse–fine chi ghosts; every seam Hamiltonian change is exactly zero. The standard-mask seam Mom/GaussE values already match T2's legacy measurements to the CSV's final rounding (relative differences ≤2.23e-16). This is direct evidence that the seam failures predate the new transfer and are not caused by its prolongation or restriction.
3. **Patch/corner source:** independently interpolate current coarse chi with eight-point degree-seven Lagrange weights and replace *only* coarse–fine chi stencil entries in the offline constraint replay. No valid sample, other variable, axis ghost or same-level ghost is changed. This isolates chi prolongation feeding native second derivatives. Exact P8 moments through degree seven and 64 tensor-monomial Float64 checks pass; maximum error is 2.168404344971009e-19. This probe is not installed in production.

| Class / mask, Ham | Existing RMS, low / mid / high | Offline P8-chi-ghost replay RMS, low / mid / high |
| --- | --- | --- |
| convex_corner / far | 1.766813e-10 / 6.281243e-11 / 1.909186e-11 | 1.353243e-11 / 6.605970e-12 / 9.265502e-12 |
| convex_corner / far_core | 1.608889e-10 / 6.281243e-11 / 1.909186e-11 | 1.760231e-12 / 6.605970e-12 / 9.265502e-12 |
| convex_corner / outer_ring_core | 1.618437e-08 / 5.540547e-09 / 1.283139e-09 | 4.433796e-10 / 9.242201e-11 / 1.893993e-10 |
| patch_edge / far | 2.255775e-10 / 4.780205e-11 / 1.742588e-11 | 9.206305e-11 / 2.192481e-11 / 1.509934e-11 |
| patch_edge / far_core | 2.343003e-10 / 5.102712e-11 / 1.856720e-11 | 2.227455e-11 / 7.326503e-12 / 1.563158e-11 |

The low/mid reductions isolate a resolved chi-prolongation truncation contribution, especially the outer_ring_core corner (1.6184e-8 → 4.4339e-10 on low, 5.5405e-9 → 9.242e-11 on mid). They do not demonstrate an algebraic bug in the degree-five operator, whose polynomial and smooth-field checks pass. On high the far corner P8 replay is 9.266e-12, still compatible with the initialization sensitivity; neither its small RMS nor a ratio through that transition establishes a fourth-order interface gate. A separate native Ricci term split uses the same current metric/Gamma with chi-derivative terms omitted solely to report the existing metric contribution: high far corner is 2.793e-13 versus total 1.909e-11, and high far patch edge is 5.993e-13 versus total 1.743e-11. Current diagonal-metric stencil spans are 0 to 4.441e-16. Thus the observed Hamiltonian is not explained just by metric/Gamma arithmetic; chi derivatives and the t=0 sampling precision matter. Both term values are in `t4b-roundoff.csv`; no altered diagnostic quantity feeds evolution.

No new production non-roundoff copy/restriction defect is established. The remaining resolved interface component is localized to chi prolongation followed by second differentiation; the seam failures above the floor are dominated by changing cell-class membership and native truncation. The finest Hamiltonian floor statement remains sensitivity-compatible, not a saturation proof or an acceptance waiver.

## Reproduction and completed short runs

`PYTHONDONTWRITEBYTECODE=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t4b-diagnose.py /private/tmp/ems-t4` builds the standalone native audit, checks original-value identity, verifies every original RMS/census, compares same-level ghost bits, and regenerates `t4b-roundoff.csv`, `t4b-layout.csv`, `t4b-identity.csv`, and `t4b-localization.csv`. It reads only the completed E t=0 plots and their completed diagnostic snapshots. `t4b-case-summary.csv` joins all 19 failures, spacings and sensitivity ratios. The CAS witness is `scripts/cas/t4b-verify.py`, with transcript `t4b-verify.json`.

The three snapshots were queued with `python Tests/EMSNative/t4-run.py --detach /private/tmp/ems-t4/t4b/state-plan.json`; all markers are 0. Each command was `/private/tmp/ems-t4/point.ex <absolute parameter path below>`, OMP_NUM_THREADS=1, from the stated directory. Each run has stdout/stderr in run.log, an atomic done.exit and no evolution steps. The maximum conservative child RSS is 545652736 bytes, below 12 GB. The detached reference directories and markers were not read, polled, changed or waited for in T4b.

| Scale | Parameter (tree relative) | Run directory | Log / atomic marker | Seconds | Conservative peak RSS bytes |
| --- | --- | --- | --- | ---: | ---: |
| low | Tests/EMSNative/params/t4b-E-low-state.txt | /private/tmp/ems-t4/t4b/E-low | run.log / done.exit=0 | 11.835895 | 208977920 |
| mid | Tests/EMSNative/params/t4b-E-mid-state.txt | /private/tmp/ems-t4/t4b/E-mid | run.log / done.exit=0 | 27.136814 | 386760704 |
| high | Tests/EMSNative/params/t4b-E-high-state.txt | /private/tmp/ems-t4/t4b/E-high | run.log / done.exit=0 | 62.459082 | 545652736 |

`COMMIT-MANIFEST-T4.txt` is refreshed for this analysis and its artifacts. No commit or production change was made in T4b.

# T4c: completed reference evolution and matched class samples

**Result: PARTIAL.** All three point-mode runs reached 2.5 M and exited 0. Ring and far pass order ≥ 3 for Hamiltonian, momentum and GaussE at both pairs at 1 M; far also passes at 2.5 M. The horizon shell passes at 2.5 M but fails Hamiltonian and the second momentum pair at 1 M. At 2.5 M the ring remains below 3 for Hamiltonian and momentum. The near-hole mask also has failures. Thus the requested exterior gate at both times does not pass. Of the 48 full-mask individual pair orders, 31 reach 3; 14 of 24 full-mask constraint cases pass both pairs. The earlier registered T4 t=0 class gate remains BLOCKED.

## Completed runs and log/field audit

No simulation was started in this analysis. The directories, commands and detached markers are those recorded in T4 F; `t4-runs.csv` now records completion. There is no separate pout file: `run.log` contains the serial pout stream and redirected stdout/stderr. Every log ends with `GRChombo finished.`. All runs use the frozen `point.ex`, `amr_transfer=point`, sigma=1, four threads, nan_check=1, and chi/lapse floors 1e-12.

| Scale | Atomic marker | Exit | Wall seconds | Conservative peak RSS bytes |
| --- | --- | ---: | ---: | ---: |
| low | /private/tmp/ems-t4/evolution/ref-low/done.exit | 0 | 108.320360 | 121913344 |
| mid | /private/tmp/ems-t4/evolution/ref-mid/done.exit | 0 | 202.674329 | 251084800 |
| high | /private/tmp/ems-t4/evolution/ref-high/done.exit | 0 | 704.000926 | 609697792 |

The peak child RSS is 609697792 bytes, below the 12 GB gate. Missing required variables: 0. The logs have 16 distinct optional parameters using their printed defaults (listed in `t4c-run-audit.csv`), including num_ghosts=3 and write_plot_ghosts=0; sigma and both floors are explicitly set in all three input files. Error/warning/nonfinite log lines: 0; nonfinite tokens in the three numeric .dat outputs per run: 0. Each run has 28 initial reader messages, all before the first `GRAMRLevel::advance`; reader messages after that point: 0. The analysis itself never opens the initial-data profile.

All 18 plots (t=0,0.5,1,1.5,2,2.5 M) contain the expected ten components, and the three final checkpoints contain all 28 evolution components. All stored numeric field arrays are Float64 and finite, including checkpoint ghosts and covered coarse cells. GaussB is exactly zero in every valid cell of all 18 plots; checkpoints do not store GaussB as a diagnostic. Saved valid and ghost chi/lapse values at or below their 1e-12 floors: 0. Saved valid chi below the existing momentum diagnostic regularizer 1e-6: 0. There are no runtime floor-activation counters, so this audit does not establish that every intermediate RK stage avoided clamping.

| Time M | Scale | Minimum saved valid chi | Minimum saved valid lapse |
| ---: | --- | ---: | ---: |
| 1 | low | 2.779341924279e-04 | 3.775059102610e-02 |
| 1 | mid | 1.222116801195e-04 | 2.680310769147e-02 |
| 1 | high | 5.402409064070e-05 | 1.893966465868e-02 |
| 2.5 | low | 2.876798629178e-04 | 7.154022880252e-02 |
| 2.5 | mid | 1.210610607483e-04 | 5.835141702248e-02 |
| 2.5 | high | 5.101475625633e-05 | 4.689607281897e-02 |

The pre-existing `min_chi.dat` is mislabeled: `EMSBH2DLevel.cpp` constructs `AMRReductions<VariableType::diagnostic>` and calls `min(c_chi)`. Evolution c_chi and diagnostic c_mod_F both have index 0, so that file reduces mod_F rather than chi. Its negative entries are not evidence of negative evolved chi. The table and floor counts above use the actual HDF5 chi/lapse arrays. No production diagnostic was changed.

## Masks, baseline and common support

The hierarchy has max_level=6, h0=2,4/3,8/9 M and finest spacings 1/32,1/48,1/72 M. The offline masks are exactly the registered T3 shells: horizon [0.63593977642346233,1.2718795528469247], near_hole [0.1,2], ring [2,4], far [4,8] M. Covered coarse cells are removed. These shell constants only select analysis samples; they are not used by evolution, transfer, tagging or a horizon finder. Constraints are measured from the saved current fields, with no static field, subtraction or constraint target supplied to the code.

The legacy baseline is restriction=off in the completed, read-only `/Users/auroradysis/Workspace/EMS-deps/worktrees/wt-native-t1/Tests/EMSNative/t3-reference-orders.csv`, SHA-256 `f017043b463094e52ad0dfcf98282f064bcfaa3d4f6eb85dfffbb831dc38e360`. Both legacy pairs are available; nothing is pending. The full-mask legacy orders below are the supplied CSV values. Recomputed native legacy RMS values and cell counts agree with that CSV on every available row at both times and all three resolutions (relative RMS tolerance 2e-13).

The cell-centred 3/2 grids have no coincident centres. Here “common cells” means the intersection of eligible physical cell footprints at the same composite AMR level, requiring each grid’s native radial-mask and T2-class membership. A shared lattice of spacing h_low/9 subdivides low/mid/high footprints into 9/6/4 cells per direction. Only intersections eligible on all three grids contribute. A native field value is retained for its footprint, with a fractional cylindrical weight; no evolved value is interpolated or resampled. The original T2 per-cell y weighting, without an h² factor between levels, is preserved. Assembly of every unfiltered mask/class reproduces T2’s native counts and RMS to relative 2e-13; full-footprint weight sums are also checked. The class columns labeled legacy-common and point-common use precisely the same support and weights. The raw legacy and point box arrays and physical bounds match at each scale/time.

The requested class rows therefore use shared footprints rather than their original independent cell censuses. The CSV retains the supplied legacy class orders in legacy_published_* columns alongside recalculated legacy-common and point-common orders, so a change in sampling can be distinguished from a scheme change. `t4c-common-support.csv` records each contributing level, overlay census, cylindrical physical measure and native cell counts; `t4c-common-norms.csv` records the RMS values and total weights. Participating native cells may be only partially weighted, and are not coincident point samples.

## Full-mask constraint orders

Each entry is low→mid / mid→high. These are constraint-to-zero RMS orders, log(RMS_coarse/RMS_fine)/log(1.5), not evolved-field self-difference orders. H=Hamiltonian, M=the joint norm of Mom1 and Mom2, G=GaussE. A large order through a small residual does not certify asymptotic convergence.

| Time M | Mask | Legacy H | Point H | Legacy M | Point M | Legacy G | Point G |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | horizon | 2.625/2.388 | 2.641/2.435 | 3.675/2.919 | 3.676/2.919 | 7.876/4.070 | 7.876/4.070 |
| 1 | near_hole | 2.694/4.612 | 2.702/4.867 | 1.899/2.758 | 1.899/2.758 | 2.029/2.436 | 2.029/2.436 |
| 1 | ring | 0.730/0.747 | 3.700/3.664 | 0.704/0.711 | 3.485/3.471 | 1.630/1.612 | 4.037/4.023 |
| 1 | far | 0.752/0.774 | 3.789/3.719 | 0.720/0.725 | 3.452/3.477 | 1.667/1.628 | 4.055/4.027 |
| 2.5 | horizon | 0.728/0.717 | 4.616/4.244 | 1.017/0.820 | 4.689/3.968 | 4.044/1.617 | 9.307/4.079 |
| 2.5 | near_hole | 2.450/4.547 | 2.461/4.833 | 2.095/2.729 | 2.095/2.729 | 1.961/2.707 | 1.961/2.707 |
| 2.5 | ring | 1.199/1.168 | 1.615/1.759 | 0.912/0.903 | 1.364/1.510 | 1.630/1.654 | 3.097/3.315 |
| 2.5 | far | 0.686/0.541 | 3.763/3.728 | 0.758/0.822 | 3.388/3.397 | 1.823/1.758 | 4.087/4.053 |

## Common-support cell-class orders

Within each scheme column the three entries are H, M, G, each low→mid / mid→high. Every present class is shown; the absent classes are specified below. Of 96 present class/constraint cases, 66 pass ≥3 at both pairs and 30 fail. No floor-based acceptance waiver is applied.

At t=1 M:

| Mask | Class | Legacy-common H, M, G | Point-common H, M, G |
| --- | --- | --- | --- |
| horizon | cartoon_axis | 3.632/3.119, 3.674/3.566, 10.153/3.615 | 3.640/3.167, 3.674/3.566, 10.153/3.615 |
| horizon | interior | 2.578/2.351, 3.558/2.812, 5.715/4.055 | 2.583/2.382, 3.558/2.812, 5.715/4.055 |
| near_hole | patch_edge | 1.576/1.634, 2.026/1.732, 2.314/2.161 | 3.051/3.399, 3.553/3.641, 2.921/3.456 |
| near_hole | axis_patch_junction | 2.463/2.766, 1.774/2.015, 3.379/3.914 | 2.794/3.169, 3.270/3.407, 2.655/3.335 |
| near_hole | cartoon_axis | 1.630/4.559, 1.538/3.022, 1.823/2.464 | 1.630/4.559, 1.538/3.022, 1.823/2.464 |
| near_hole | interior | 3.310/4.393, 3.246/3.941, 2.544/2.762 | 3.331/4.953, 3.246/3.942, 2.544/2.762 |
| ring | patch_edge | 1.020/1.080, 1.199/1.182, 2.212/1.940 | 3.004/3.287, 3.283/3.439, 3.155/3.336 |
| ring | convex_corner | 1.024/1.419, 1.284/1.535, 2.355/1.634 | 3.662/3.783, 3.560/3.707, 4.254/4.043 |
| ring | axis_patch_junction | 2.686/2.484, 2.249/2.342, 3.771/2.565 | -5.894/2.811, 3.273/3.399, 3.580/3.577 |
| ring | cartoon_axis | 0.743/0.753, 0.716/0.714, 1.580/1.585 | 3.937/3.799, 3.389/3.582, 4.121/4.068 |
| ring | interior | 0.732/0.748, 0.706/0.712, 1.609/1.604 | 3.802/3.743, 3.982/3.647, 4.008/3.996 |
| far | patch_edge | 0.342/1.236, 0.919/0.962, 2.966/1.762 | 3.167/3.111, 3.247/3.458, 3.294/3.401 |
| far | convex_corner | 3.456/1.712, 1.648/1.411, 2.165/1.771 | 2.796/3.558, 3.352/3.517, 0.932/3.129 |
| far | axis_patch_junction | -9.062/3.165, 1.449/2.890, -3.050/3.157 | 2.049/4.889, 3.009/3.244, 3.189/3.241 |
| far | cartoon_axis | 0.863/0.812, 0.707/0.718, 1.705/1.626 | 4.493/4.182, 3.760/3.816, 4.075/4.037 |
| far | interior | 0.754/0.774, 0.717/0.725, 1.637/1.628 | 3.914/3.866, 3.838/3.685, 3.946/4.006 |

At t=2.5 M:

| Mask | Class | Legacy-common H, M, G | Point-common H, M, G |
| --- | --- | --- | --- |
| horizon | cartoon_axis | 2.933/1.287, 3.326/2.977, 6.172/1.588 | 4.317/4.533, 3.508/3.453, 9.847/4.069 |
| horizon | interior | 0.619/0.666, 0.928/0.779, 3.007/1.588 | 4.597/4.327, 4.623/3.977, 7.934/4.083 |
| near_hole | patch_edge | 4.368/3.675, 1.802/1.671, 2.413/2.242 | 3.140/5.753, 2.728/5.993, 4.868/5.328 |
| near_hole | axis_patch_junction | 2.691/2.900, 3.225/1.127, 3.339/3.746 | 3.819/5.223, 5.490/4.812, 5.075/5.214 |
| near_hole | cartoon_axis | 1.382/4.333, 1.723/2.973, 1.731/2.711 | 1.382/4.334, 1.723/2.973, 1.731/2.711 |
| near_hole | interior | 3.001/4.276, 3.361/4.051, 2.528/2.900 | 3.025/4.772, 3.361/4.051, 2.528/2.900 |
| ring | patch_edge | 1.349/1.739, 1.173/1.456, 2.088/2.046 | 1.543/1.916, 1.154/1.382, 3.689/3.691 |
| ring | convex_corner | 3.119/4.406, 1.660/1.447, 2.281/0.905 | 3.671/3.821, 2.982/3.205, 7.787/1.457 |
| ring | axis_patch_junction | 2.976/0.883, 1.372/1.519, 3.055/2.845 | 3.625/3.828, 3.079/3.283, 3.668/3.568 |
| ring | cartoon_axis | 1.147/1.038, 0.726/0.927, 1.392/1.583 | 2.079/2.857, 2.985/3.167, 3.090/3.453 |
| ring | interior | 0.885/0.906, 0.831/0.795, 1.579/1.621 | 1.590/1.602, 1.767/1.644, 3.050/3.275 |
| far | patch_edge | 0.608/0.663, 0.874/0.930, 1.880/1.917 | 3.161/3.390, 3.175/3.357, 3.400/3.445 |
| far | convex_corner | 2.043/2.182, 1.614/1.724, 1.869/1.772 | 2.670/3.166, 3.418/3.546, 2.814/3.197 |
| far | axis_patch_junction | 2.281/2.417, 1.820/2.177, 2.363/1.706 | 4.236/4.604, 3.086/3.290, -0.696/3.184 |
| far | cartoon_axis | -0.076/0.276, 1.034/1.007, 1.781/1.687 | 2.813/3.118, 3.351/3.555, 4.107/4.063 |
| far | interior | 0.683/0.528, 0.743/0.824, 1.814/1.744 | 3.980/3.858, 4.184/3.695, 4.018/4.018 |

Shared native-cell census (the hierarchy and masks are unchanged between the two times):

| Mask | Class | Low | Mid | High |
| --- | --- | ---: | ---: | ---: |
| horizon | cartoon_axis | 42 | 122 | 184 |
| horizon | interior | 1836 | 3792 | 8344 |
| near_hole | patch_edge | 28 | 76 | 104 |
| near_hole | axis_patch_junction | 2 | 8 | 8 |
| near_hole | cartoon_axis | 118 | 352 | 532 |
| near_hole | interior | 5988 | 12728 | 28346 |
| ring | patch_edge | 244 | 704 | 1068 |
| ring | convex_corner | 2 | 8 | 8 |
| ring | axis_patch_junction | 2 | 8 | 8 |
| ring | cartoon_axis | 60 | 180 | 272 |
| ring | interior | 5460 | 11692 | 26196 |
| far | patch_edge | 244 | 704 | 1068 |
| far | convex_corner | 2 | 8 | 8 |
| far | axis_patch_junction | 2 | 8 | 8 |
| far | cartoon_axis | 60 | 180 | 272 |
| far | interior | 5460 | 11692 | 26196 |

box_seam and outer_boundary have no common support in any requested mask at either time. Low has one box on each refined level, whereas mid has two and high eight; its only base-level seam is covered within these masks. Independent mid/high seam samples consequently cannot form a three-grid seam class. The outer domain boundary is outside all four shells. Patch edges, convex corners and axis × patch junctions are absent from the horizon shell; convex corners are absent from near_hole. These are 24 absent mask/class/time combinations (72 constraint rows), reported as ABSENT_COMMON with zero participating counts, not as pending or convergence passes.

The remaining low ring orders at 2.5 M occur on common patch-edge samples (H 1.543/1.916, M 1.154/1.382) and common interior samples (H 1.590/1.602, M 1.767/1.644). Thus they survive matching the sampled physical cells and extend into the interior; these outputs alone do not identify a faulty transfer operation. Tiny classes have only 2/8/8 participating cells: their negative or large orders are reported as measured and are not taken as proof of a roundoff floor. Raw RMS values and counts for every constraint are in the CSVs.

## Author’s RHFinder output and reproduction

RH_activate=false and AH_activate=0 in all three point runs. No rh_surf_* or rh_f* output exists; horizon area and horizon charge are unavailable on low, mid and high, and their resolution trend cannot be measured. Qscalar and the fixed offline horizon shell are not substituted for RHFinder measurements. The author’s RHFinder source is unchanged.

Reproduce the analysis (short postprocessing only):

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/auroradysis/miniconda3/bin/python Tests/EMSNative/t4c-analyze.py /private/tmp/ems-t4 /private/tmp/ems-t3-ref /Users/auroradysis/Workspace/EMS-deps/worktrees/wt-native-t1/Tests/EMSNative/t3-reference-orders.csv
```

Authoritative T4c artifacts: `t4c-run-audit.csv`, `t4c-field-checks.csv`, `t4c-native-norms.csv`, `t4c-common-norms.csv`, `t4c-common-support.csv`, and `t4c-orders.csv`. These supersede the empty launch-time t4-evolved-* stubs for this completed analysis; use t4c-analyze.py for the common-support class comparison. The analysis verifies native T2 norms/censuses, all supplied legacy rows, common-support eligibility, identical legacy/point layouts and full-footprint weighting. The manifest is refreshed. No production source, run output, Chombo file, input repository or wt-native-t1 file was modified; no commit was made.
