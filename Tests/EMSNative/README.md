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
