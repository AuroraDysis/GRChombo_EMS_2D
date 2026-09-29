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
