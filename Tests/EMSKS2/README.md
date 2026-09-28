# EMSKS 2 C++ readback and t = 0 grid result

**Current T4C outcome: PRE-ASYMPTOTIC.** At the authorized fourth resolution, the joining Hamiltonian L2 order is 3.771 (B), 3.891 (E) and 3.818 (reference), with both momentum components above 3.77. The h⁻⁴-scaled radial profiles tighten on each successive pair. The original three-level NOT-CONVERGENT result and the T4B floor-margin stop remain documented below as earlier measurements. Only the new finest diagnostic runs use `min_lapse = min_chi = 1e-9`; every floor remains inactive and every active and ghost cell passes the unchanged factor-100 check. This is a t = 0 grid result, not an evolution claim.

## Inputs and fingerprints

The read-only Julia worktree is commit `552ab65acdcf4bddc359f672c2ca465a9e7b5322`; this C++ worktree started at `34f2ef0b307428209c92b8e3b33d6f8df02ed45c`. All 17 fixture files match the supplied `manifest.txt` byte counts and SHA-256 values. The three top-level profiles are byte-identical to the fixture copies. The contract `docs/ks-puncture-format.md` has SHA-256 `82947f2ceaadc9330040898c89c17c7e1c9c86f98c6ea8d93d0c675dfacb5914`. The source `ks_puncture.jl` and `ks.jl` hashes are `9edf76b2cc72c16d274e7d97530451c00807c7124088314192e664b5c4cc4e03` and `23ebcc9888d60c0dc72727b256d4ddda136d086e6f9de9a766c16956eb93fb52`.

| Member | complete `.ks2` SHA-256 | comparisons | maximum normalized fixture error |
|---|---|---:|---:|
| B | `6fe2e42892c3842a78d15ed8b1cdaccacf75d545aad8917085b8049a55a17d01` | 730 | 2.337e-15 |
| E | `b2faa29e1c81cf0c871e6a4a6c8eab39d0ab01f26f0a6c0675534b8895768912` | 730 | 7.251e-16 |
| reference | `fb46ceecd7565c5afb5b04b4ace4dd0790a832844d3f57a28333700230870533` | 730 | 1.343e-15 |

`EMSKS2FixtureTest.cpp` reads all source element, radial, Cartesian CCZ4 and reference gauge rows. It checks `|C++ − Julia|/(1 + |Julia|)` against `1e-12` for values and `1e-10` for first derivatives. The supplied tables have no second-derivative columns; the registered `1e-8` second-derivative class therefore has no rows to test. All seven supplied malformed files are rejected; exact reasons are in [`fixture-run.log`](fixture-run.log). The reader checks the SHA-256 trailer, sorted mandatory metadata, domains, coupling, mass/charge/time normalization, critical cylinder, five source-derived endpoint Taylor coefficients, chart antiderivative, chart derivative and source consistency.

## Grid method and results

The cartoon boxes span approximately `[-2.5 r_h, 2.5 r_h] × [0, 2.5 r_h]`, with the puncture on a shared vertex. The four fixed physical masks are join `[ρ_a,r_m]`, KS collar `[r_m,r_h]`, exterior `[r_h,2r_h]`, and far `[2r_h,2.3r_h]`; each additionally requires `ρ−4h ≥ ρ_a/2`. The far upper limit leaves at least `0.2 r_h` before the domain boundary. This cut is fixed across each resolution ladder. The test samples the setter on all ghost cells, calculates discrete Γ̃ on the interior plus two ghost layers, and invokes the production `Constraints<CouplingFunction>` and `EMSCartoonGaussConstraints` operators. Electric Gauss is the signed diagnostic `D_i(F E^i)`; no reference RHS is subtracted. Magnetic Gauss is identically zero. Analytic fixture Γ̃ and discrete grid Γ̃ are both retained; their difference is tabulated.

The [grid norms](grid-results.csv) contain L2 and L∞ of H, Mx, My, electric Gauss, magnetic Gauss, and Γ̃ error for every member, mask and resolution. [Observed orders](grid-orders.csv) contain both adjacent pairs. The table below gives the finer-pair L2/L∞ orders; full norms and mask counts are in those CSVs.

| Member | Mask | H | Mx | My | electric Gauss | Γ̃ error |
|---|---|---:|---:|---:|---:|---:|
| B | join | 3.18 / 2.80 | 3.28 / 3.12 | 3.28 / 3.12 | 3.58 / 3.37 | 3.57 / 3.39 |
| B | KS collar | 7.22 / 6.59 | 6.80 / 6.34 | 6.80 / 6.34 | 7.51 / 7.93 | 7.71 / 7.21 |
| B | exterior | 3.96 / 3.92 | 4.00 / 3.99 | 4.00 / 3.99 | 4.00 / 4.00 | 4.00 / 3.99 |
| B | far | 0.74 / 0.32 | 4.00 / 3.99 | 4.00 / 3.99 | 4.00 / 3.98 | 4.00 / 4.00 |
| E | join | 3.56 / 3.31 | 3.60 / 3.49 | 3.60 / 3.49 | 3.76 / 3.62 | 3.76 / 3.66 |
| E | KS collar | 7.50 / 6.89 | 7.03 / 6.49 | 7.03 / 6.49 | 5.41 / 7.58 | 8.12 / 7.39 |
| E | exterior | 3.16 / 3.21 | 4.00 / 3.99 | 4.00 / 3.99 | 4.00 / 3.98 | 4.00 / 4.00 |
| E | far | −1.44 / −1.48 | 4.00 / 3.95 | 4.00 / 3.95 | 4.00 / 3.92 | 4.00 / 3.97 |
| reference | join | 3.31 / 2.96 | 3.39 / 3.23 | 3.39 / 3.23 | 3.63 / 3.40 | 3.63 / 3.49 |
| reference | KS collar | 7.02 / 6.27 | 6.89 / 6.29 | 6.89 / 6.29 | 5.39 / 7.03 | 7.55 / 6.75 |
| reference | exterior | 4.00 / 3.96 | 4.00 / 3.99 | 4.00 / 3.99 | 4.00 / 3.97 | 4.00 / 4.00 |
| reference | far | 2.57 / 1.86 | 4.00 / 3.99 | 4.00 / 3.99 | 4.00 / 3.99 | 4.00 / 4.00 |

Far-mask Hamiltonian reaches an arithmetic/readback floor in the fine pair (for E, L2 grows from `3.651e-11` to `9.917e-11`); its poor fine-pair order is reported, not used to adjust the mask. The B and reference join failures occur at L2 `1e-1` and `9e-2`, far above that floor. The joining mask's low L∞ orders are an additional failure. Electric and magnetic Gauss meet their expected behavior where above the floor.

The [first puncture cell audit](puncture-audit.csv) reports α, χ, α/ρ^ν, χr₀²/ρ², E_iE^i and its limiting ratio. It makes no order claim at the puncture. At the finest B/E grids the smallest lapse/floor factors are 184/356 and smallest chi/floor factors are 1208/647; the corresponding factors for the finest reference grid are 5558/1751. All active and sampled ghost cells pass the required factor 100. Every [grid row](grid-results.csv) has zero chi and lapse floor activation, as do the three production single-level initialization runs. The new path leaves gauge B, Θ, Ξ, Λ and the magnetic field zero and uses the file's α_* and β_*.

## Baseline regression

[`regression-results.csv`](regression-results.csv) compares this worktree against a separately compiled archive of commit `34f2ef0`. The identical `EMSKS2Regression.cpp` is compiled against each source tree. It hashes all 28 evolved variables in every cell of the test FArrayBox, including five ghost layers, after the unchanged setter, discrete Γ̃ pass and old gauge initialization. SHA-256 matches bit for bit for `legacy_dat`, trumpet reference, trumpet echo-B and EMSCTT. The harness tests the direct initialized ghost values; it does not execute the full production boundary-condition fill. The old production branches and their `fillAllGhosts()` calls are unchanged by this patch.

| Case | cells including ghosts | SHA-256 on both builds |
|---|---:|---|
| legacy_dat | 10,212 | `d34ebb047cf788af8a5a137562349d2b5cdd5dc010284eecbfb81221bd89eeb1` |
| trumpet reference | 10,212 | `27983762ef33929ffa6081e1f166f01fc1d46d0e05cfbaacbb1cad89fdcdd212` |
| trumpet echo-B | 10,212 | `43e6f7c54b60b061093f2f004dbe428700d8b9f611edb8e2e6c59820c8c30c64` |
| EMSCTT | 2,484 | `e5a8040030b6ee6093ee4f66764e03939f0c42e91ab36a2ca3ce5fdb1930cc5b` |

## Reproduce and resource record

Run from the C++ worktree root. `JL` is the supplied read-only Julia worktree.

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
JL=/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-LSb44u
make -C Examples/EMS all DIM=2 -j 8
make -C Tests/EMSKS2 all DIM=2 -j 8
/opt/homebrew/bin/g++-16 -std=c++17 -O2 -I Source/InitialConditions/EMSBH/1D_SOL -I Source/utils Tests/EMSKS2/EMSKS2FixtureTest.cpp -o /private/tmp/EMSKS2FixtureTest
/private/tmp/EMSKS2FixtureTest "$JL/artifacts/echo-evolution/t3/fixtures/emsks2"
/opt/homebrew/bin/g++-16 -std=c++17 -O2 -g -fsanitize=address,undefined -fno-omit-frame-pointer -I Source/InitialConditions/EMSBH/1D_SOL -I Source/utils Tests/EMSKS2/EMSKS2FixtureTest.cpp -o /private/tmp/EMSKS2FixtureTest.san
/private/tmp/EMSKS2FixtureTest.san "$JL/artifacts/echo-evolution/t3/fixtures/emsks2"
python3 Tests/EMSKS2/run_grid.py "$JL/artifacts/echo-evolution/t3"
```

The B M/512 real-path pilot took 0.825 s. Extrapolating by measured cell count gave less than 90 s for all nine; the completed nine-run ladder took about 43 s of stage time. The [stage log](grid-stages.log) records elapsed seconds, measured peak RSS and items at every stage boundary. The slowest unit was E M/2048: 17.755 s, 501,579,776 byte peak RSS. `memory_pressure -Q` reported 32% free before the ladder; this maximum process RSS is well below half of 24 GiB and leaves more than one fifth free. The final production and test builds both returned 0 in 14.1 s combined; the O2 fixture run took about 2.2 s including compilation, ASan/UBSan about 5.3 s including compilation, baseline regression build 4.2 s, and the four paired regression cases each took under 0.6 s per process. The last production `max_steps=0` initializations from `Examples/EMS/params-ks2-{B,E,reference}.txt` took 0.926, 0.655 and 0.250 s with exit 0 and `chi=0 lapse=0` printed on level 0. The baseline build used `git archive 34f2ef0` for `Source` and required EMS headers under `/private/tmp/emsks2-baseline`; the exact regression source and compiler flags were identical across builds.

The claim boundary is t = 0 readback and initialization on the evolution code's own grid. No evolution, gauge stability or global puncture norm is established.

## T4B finer-grid feasibility stop (2026-09-28; resolved by T4C authorization below)

**Outcome: BLOCKED by the specified lapse-floor margin, before the requested fourth-grid diagnostic.** None of the registered numerical outcomes PRE-ASYMPTOTIC, STRUCTURAL or FLOOR applies: the B grid fails an initialization precondition, not a Float64 arithmetic-floor test. The reader and setter, floors, masks, constraints and evolution code were left unchanged. The earlier three-level NOT-CONVERGENT finding remains the last measured join order.

The [finer-grid preflight](finer-feasibility.csv) samples the exact `EMSKS2Profile::object(h/2,h/2)` path at an active cell nearest the puncture in the existing vertex-centred box. It applies the setter's `alpha >= 100 * min_lapse && chi >= 100 * min_chi` check with both floors fixed at `1e-8`. A single failed active cell proves the full B grid cannot initialize under the unchanged contract; passing this one-cell check does not certify all cells of E or the reference.

| Member | h | nearest ρ | nearest α | α / min_lapse | nearest χ | χ / min_chi | one-cell check |
|---|---:|---:|---:|---:|---:|---:|---|
| B, M/4096 | 0.0019185238278559367 | 0.0013566012085449053 | 7.351893830993478e-7 | 73.51893831 | 3.019663379799114e-6 | 301.96633798 | reject |
| E, M/4096 | 0.001397210020898587 | 0.0009879766805191885 | 1.408617602306571e-6 | 140.86176023 | 1.618305079288644e-6 | 161.83050793 | pass |
| reference, M/256 | 0.002230811200079717 | 0.001577421727123268 | 2.102715486282423e-5 | 2102.71548628 | 4.378440188942360e-6 | 437.84401889 | pass |

The B lapse exceeds the floor itself but falls below the required `1e-6` margin. The setter's [check](../../Source/InitialConditions/EMSBH/EMSBH_ks2_read.hpp) would reject that active cell even with zero floor activations. The previous B M/2048 nearest-cell lapse was `1.8414241101576945e-6` (factor 184.14); the new exact readback follows the file's `ν=1.3246940667627503` puncture scaling. The supplied Julia report places the B peak of `M|K|` at `ρ/M ≈ 0.1096` and its M/512 10–90% joining width at about 12.57 cells, but no new localization is inferred from those source quantities.

The existing [stage log](grid-stages.log) measured B M/2048 at 9.054 s and 271,450,112 byte peak RSS, E M/2048 at 17.755 s and 501,579,776 bytes, and reference M/128 at 5.741 s and 203,653,120 bytes. Fourfold cell-count extrapolations to the requested grids are about 36, 71 and 23 s, with conservative fourfold peak RSS bounds 1.09, 2.01 and 0.81 GB. `memory_pressure -Q` reported 72% free on the 24 GiB host before the preflight. Resource limits therefore did not trigger a domain reduction. The preflight compile plus three reads took 1.75 s wall; the printed per-read times were 0.314, 0.011 and 0.012 s. No fourth-grid production run was launched.

To make the exact full-box B M/4096 diagnostic possible, the contract would need to authorize a lower `min_lapse`: at this cell it must be at most `7.3518938309934775e-9` for the factor-100 check, with a new all-cell and ghost audit. That floor change is outside this tranche. Accordingly there are no new fourth-resolution results, third-pair orders, radial profiles or relative-size ratios; the earlier [grid norms](grid-results.csv) and [orders](grid-orders.csv) remain unchanged. Reproduce the preflight from the worktree root with:

```sh
/opt/homebrew/bin/g++-16 -std=c++17 -O2 -I Source/InitialConditions/EMSBH/1D_SOL -I Source/utils Tests/EMSKS2/EMSKS2FinerFeasibility.cpp -o Tests/EMSKS2/EMSKS2FinerFeasibility
JL=/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-LSb44u
for spec in B:4096 E:4096 reference:256; do member=${spec%%:*}; n=${spec#*:}; Tests/EMSKS2/EMSKS2FinerFeasibility "$member" "$n" "$JL/artifacts/echo-evolution/t3/$member.ks2"; done
```

## T4C: fourth-grid join diagnosis (2026-09-28)

**Pre-registered outcome: PRE-ASYMPTOTIC.** The finest-pair join L2 orders exceed 3.7 for H, both momentum components, signed electric Gauss and discrete-minus-analytic Γ̃ for B, E and reference. The h⁻⁴-scaled H and |M| radial profiles move closer on the later pair. The B Hamiltonian join L∞ order is 3.641, below 3.7 but above the structural threshold of approximately 3.5; the registered pivot concerns the L2 order. Magnetic Gauss is exactly zero at all four resolutions. The earlier three-level NOT-CONVERGENT classification was a finite-resolution finding, now resolved by the authorized finer diagnostic.

The [four-level norms](grid-results-four.csv) give L2 and L∞ for H, Mx, My, signed G_E, G_B and Γ̃ error on the **unchanged join and KS collar masks** at every one of four resolutions. The [four-level orders](grid-orders-four.csv) give each of the three consecutive pair orders in both norms. Blank G_B orders mean zero divided by zero, not a failed residual. The table gives the newest pair's L2/L∞ orders; Mx and My agree to printed precision.

| Member | Mask | H | Mx = My | G_E | Γ̃ error |
|---|---|---:|---:|---:|---:|
| B | join | 3.771 / 3.641 | 3.779 / 3.714 | 3.873 / 3.765 | 3.870 / 3.825 |
| B | KS collar | 7.649 / 7.176 | 7.341 / 6.630 | 4.253 / 6.612 | 8.337 / 7.690 |
| E | join | 3.891 / 3.835 | 3.888 / 3.853 | 3.935 / 3.888 | 3.934 / 3.894 |
| E | KS collar | 7.212 / 7.220 | 7.115 / 6.659 | 4.012 / 5.067 | 7.493 / 7.706 |
| reference | join | 3.818 / 3.707 | 3.819 / 3.764 | 3.892 / 3.850 | 3.894 / 3.833 |
| reference | KS collar | 7.524 / 7.338 | 7.266 / 6.801 | 4.028 / 5.112 | 7.856 / 7.919 |

For the decision residual, the four join H L2 norms and three consecutive orders are:

| Member | H L2 at four resolutions, coarse to fine | Pair orders, coarse to fine |
|---|---|---|
| B | 5.394829, 1.301674, 0.144063, 0.010556 | 2.051, 3.176, 3.771 |
| E | 2.840976, 0.442413, 0.037535, 0.002530 | 2.683, 3.559, 3.891 |
| reference | 4.038296, 0.878050, 0.088343, 0.006265 | 2.201, 3.313, 3.818 |

The [join profiles](join-profiles.csv) contain 40 fixed ρ bins from ρ_a to r_m for each of three diagnostic resolutions per member. Each row gives the angular maximum and RMS of |H| and |M| multiplied by h⁻⁴, its sampled areal r and blend z, and relative-size ratios. The required two finest profiles are B/E M/2048 and M/4096 and reference M/128 and M/256. The additional earlier profile measures whether overlay disagreement decreases across successive pairs. The relative L2 distance between two h⁻⁴-scaled 40-bin RMS curves is `sqrt(sum((coarse−fine)^2)/sum(fine^2))`:

| Member | H earlier → later discrepancy | |M| earlier → later discrepancy | Finest peak ρ, areal r, z |
|---|---:|---:|---|
| B | 0.4510 → 0.1516 | 0.4119 → 0.1464 | 0.87154, 0.87208, 0.75374 |
| E | 0.2685 → 0.0746 | 0.2464 → 0.0762 | 0.87162, 0.87205, 0.76154 |
| reference | 0.3933 → 0.1217 | 0.3501 → 0.1184 | 0.85339, 0.85557, 0.68001 |

Both H and |M| RMS and angular maximum peak in bin 38 for B/E and bin 37 for reference; these lie on the outer shoulder of the blend near the supplied peak-|K| radii, rather than at the puncture or r_m endpoint. Bins 37–39 hold respectively 99.974% / 99.989% (B), 99.961% / 99.969% (E), and 99.210% / 99.509% (reference) of the join's H / |M| L2 power at the finest grid. The exact [concentration fractions](join-concentration.csv), [peak values and locations](join-localization.csv), and [pair discrepancies](join-overlay.csv) are saved separately.

The relative-size denominator follows the production `Constraints::constraint_equations` decomposition. For H it is `|R| + |(2/3)K²| + |Ã_ijÃ^ij| + 16π|ρ_EMS|` (cosmological constant zero). For each M_i it sums the absolute values of its `−(2/3)∂_iK` term, the cartoon `reg_07` term, each `Γ_ww Ã` term, each `h^jk(∇_kÃ_ji − 3Ã_ij ∂_kχ/(2χ_reg))` term, and `−8πS_i`; `χ_reg=max(1e-6,χ)` as in the production operator. The vector M denominator is the Euclidean norm of the two component denominators. Each bin's reported ratio is RMS residual divided by RMS denominator over its cells. The test-only decomposition uses the production Ricci and EMS tensor routines and checks that it reconstructs H, Mx and My at every join cell.

| Finest member | H ratio at residual peak | |M| ratio at residual peak | Largest H ratio in any bin | Largest |M| ratio in any bin |
|---|---:|---:|---:|---:|
| B | 1.6013e-4 | 1.0453e-4 | 1.5083e-3 | 1.2944e-4 |
| E | 4.9617e-5 | 2.4741e-5 | 3.1717e-4 | 2.9943e-5 |
| reference | 1.1055e-4 | 8.6040e-5 | 8.7676e-4 | 1.8819e-4 |

The largest relative ratios occur in bin 39, closest to r_m, where the absolute residual RMS is lower than at the bin-37/38 peak. All 40 bin ratios for H, M, Mx and My are in [join-profiles.csv](join-profiles.csv).

The fine-grid [parameter files](params-fine-B.txt) set both floors to `1e-9` only for B/E M/4096 and reference M/256. The earlier ladder still uses `1e-8`. The harness scans the complete state box including five ghost layers after the unchanged setter; minimum α/`min_lapse` and χ/`min_chi` were respectively **735.189 / 3019.66** (B), **1408.62 / 1618.31** (E), and **21027.2 / 4378.44** (reference). All are above 100; all floor counts are zero. The old-grid diagnostic reruns reproduce the saved join and collar norms to at worst `1.293e-11` relative. The stored-field term reconstruction check passed at every join cell.

The full boxes and masks were unchanged. Before the largest run, `memory_pressure -Q` reported 60% free on the 24 GiB host. The B M/2048 real-path diagnostic pilot took 9.198 s and 276,856,832 bytes peak RSS; earlier member-specific grid timings and fourfold cell-count extrapolation put the full bundle below five minutes. Nine diagnostic runs took 173.536 s in total stage time. The three new finest units took 37.856, 72.206 and 22.541 s, with peak RSS 1,035,010,048, 1,851,293,696 and 796,000,256 bytes. Per-run margins, wall times and RSS are in [grid-finer-resources.csv](grid-finer-resources.csv); every stage boundary with elapsed time, RSS and items done/total is in [grid-finer-stages.log](grid-finer-stages.log). The test build returned 0 in 4.275 s; `analyze_finer.py` returned 0 with 144 result rows, 36 order rows, 360 profile rows, maximum old-grid rerun difference `1.293e-11` relative, and `join_L2_gate=True overlay_tightens=True outcome=PRE-ASYMPTOTIC`.

The three `.ks2` and contract SHA-256 values above were rechecked before the T4C runs; the Julia report SHA-256 is `7b9052a7b24b9b2b1b71f3745410ce9ad71b992ba3433d6aa5765a8d3e87da8b`. Output SHA-256: `grid-results-four.csv` `a2364781dc6a74d6df03c785c792490dcc3d6de7309d81ce44e0626fe47fe164`; `grid-orders-four.csv` `ddf6bcbc0e609658b5467af18cc4a580b6bf879a26da5f0a6e42ed729c535b0a`; `join-profiles.csv` `da56fcb9cb021dfab532127fcbdece33b3f10245e3c1b12e0bb63eeba1346369`; `join-localization.csv` `be142a3a4a7fc22e44d3c504180ba43fcb054d9e8458307f66832a30f9b55999`; `join-overlay.csv` `e6ea6e3dd22c29545d90187b47c18b4ad0008aa08a69332057871298159d8058`; `join-concentration.csv` `8a63a539d22bdf79ce21068a56fd5e94708665e7c79b3f7ec55fa012b3a923f6`.

Reproduce from the worktree root after setting `CHOMBO_HOME` and `JL` as above. The nine runs use the same executable; for B/E substitute N=`1024 2048 4096`, for reference N=`64 128 256`. Pass the fine parameter file as the fourth argument only on the largest N, followed by the profile CSV path. Redirect each run's stdout to `grid-MEMBER-N.out` and stderr to `grid-MEMBER-N.log`, then run `python3 Tests/EMSKS2/analyze_finer.py`. For example:

```sh
env CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib make -C Tests/EMSKS2 all DIM=2 -j 4
EXE=Tests/EMSKS2/EMSKS2GridConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex
env OMP_NUM_THREADS=1 "$EXE" B 2048 "$JL/artifacts/echo-evolution/t3/B.ks2" Tests/EMSKS2/profile-B-2048.csv > Tests/EMSKS2/grid-B-2048.out 2> Tests/EMSKS2/grid-B-2048.log
env OMP_NUM_THREADS=1 "$EXE" B 4096 "$JL/artifacts/echo-evolution/t3/B.ks2" Tests/EMSKS2/params-fine-B.txt Tests/EMSKS2/profile-B-4096.csv > Tests/EMSKS2/grid-B-4096.out 2> Tests/EMSKS2/grid-B-4096.log
python3 Tests/EMSKS2/analyze_finer.py
```

The claim boundary is t = 0 initial data on the evolution code's own grid and these fixed masks. No evolution or gauge stability follows.

## T5: opt-in reference-stationary gauge and short evolution (2026-09-28)

**G2 outcome: NOT-STATIONARY under the registered all-mask fourth-order test.** The complete production CCZ4 + EMS RHS, with reference-stationary gauge and no reference-RHS subtraction, reaches approximately fourth order in the join and KS collar. Exterior and far-mask components reach a `1e-12`–`1e-10` arithmetic/readback floor at the finest grids and do not show fourth-order orders there. For example, KO-on far K has finest-pair L2 order −2.020 (B), −2.007 (E), −1.841 (reference), at respective fine norms `3.303e-11`, `1.259e-10`, `1.088e-10`. The roughly `h^-2` behavior for several second-derivative components is consistent with Float64 subtraction noise, an inference from the ladder rather than a corrected residual. The fixed masks, files, floors, CCZ4 equations and EMS sector were not changed to make these numbers pass.

The supplied contract, Julia report and B/E/reference profiles were rehashed before T5: SHA-256 values remain `82947f2c…acb5914`, `7b9052a7…e87da8b`, `6fe2e428…17d01`, `b2faa29e…768912` and `fb46ceec…70533`, respectively, matching the complete T4C fingerprints above. Output SHA-256: `rhs-results.csv` `0939ff2b…fb1a284`, `rhs-orders.csv` `d62a2007…37412ad6`, `rhs-cost.csv` `001188ea…399c5fc`, `smoke-B-harmonic.csv` `15828d40…4b0bf3e`, `smoke-cost.csv` `b77ba697…9aa0263`, `gauge-regression-results.csv` `35d0deb5…20f86a`.

**G3 outcome: SMOKE-FINITE for B with harmonic slicing through `t/M=0.046875`.** This is 12 coarse steps and 96 finest steps on a fixed four-level hierarchy, with no nonfinite values or floor activations and an inward-pointing sampled outgoing characteristic at `r_s=0.95r_h` throughout. It does not establish long-term gauge stability. The optional onepluslog smoke was not launched because the measured local compute total, including pilots, reached approximately 9.2 minutes; its projected additional 57 s would cross the approximately 10-minute local limit.

### Gauge path and direct-evaluation control

`gauge_type=reference_stationary` is legal only with `ems_data_format=emsks2`; the default is still `experimental`. The new branch prescribes the file's Cartesian shift, sets its shift and gauge-B RHS exactly to zero after KO, and uses equation (13) with `reference_f=harmonic` by default or `reference_f=onepluslog` and `reference_onepluslog_n=2`. The reference profiles are sampled from the `.ks2` file at every patch cell and ghost centre during initialization, new-grid creation, and restart; no coarse-grid interpolation is used for the reference cache. The full RHS uses `sigma_KO=1` times the file's taper; the old gauge branch and its `sigma` remain unchanged. The new path checks lapse and chi at every stage and post-step active and ghost cell before positivity enforcement, aborting on a floor activation or nonfinite lapse/chi. The existing active-cell NanCheck covers every evolved variable after a full step.

The RHS harness samples about 1024 cells from every full reference cache plus all 81 cells of a newly created box, using a freshly loaded profile. All six cached components (α_*, q_*, K_*, β_*^x, β_*^y, d) agree **bit for bit** with direct `profile.object` evaluation at those samples. Its algebraic fixed-point check gives lapse RHS `≤1e-15 max(1,α_*)` for both harmonic and onepluslog; the production RHS asserts shift and gauge-B RHS equal exactly zero after KO on every masked cell. Production startup created all four levels and repeatedly rebuilt the cache on new grids, with level-0/1/2/3 final cache fills of `627200/677448/182408/270848` cells. A checkpoint restart was not executed, so the `postRestart()` hook is compiled but lacks an end-to-end restart control.

At the largest B/E/reference diagnostic levels, cache fill times were `18.236/37.069/11.208 s` and one KO-on RHS pass took `2.323/5.152/1.585 s`, respectively. The cache is a one-time new-grid/restart cost, about `7.1–8.0` RHS passes at these levels. The [per-resolution cost and peak-RSS table](rhs-cost.csv) includes the measured floor margins. The largest measured process peak was `3.093 GB` for E M/4096; `memory_pressure -Q` showed 47% free before that dispatch, so the projected peak stayed below half the 24 GiB RAM and left over one fifth free. Fine B/E/reference runs alone used the authorized `1e-9` floors; earlier levels and the smoke used `1e-8`.

### Full t = 0 RHS ladder

The [RHS norms](rhs-results.csv) contain L2 and L∞ for all 15 groups, four masks, both KO settings, and four resolutions per member; [pair orders](rhs-orders.csv) contain all three consecutive pairs. Vector/tensor group magnitudes are Euclidean norms of their listed components. `B_gauge` is a synthesized zero row backed by the harness's exact-zero assertion on every masked cell; all other rows come directly from the production `CCZ4Cartoon<ReferenceStationaryGauge>` operator. The boxes, centres, masks and floor check match the T4C ladder. The no-KO and KO-on finest-pair **ranges** below exclude exactly zero groups; full per-group and L∞ orders are in the CSV.

| Member | Mask | KO off L2 order range | KO on L2 order range |
|---|---|---:|---:|
| B | join | 3.588–3.860 | 3.767–4.628 |
| B | KS collar | 3.995–8.586 | 4.404–9.257 |
| E | join | 3.774–3.928 | 3.889–4.736 |
| E | KS collar | 3.985–8.102 | 4.015–9.011 |
| reference | join | 3.664–3.890 | 3.771–4.167 |
| reference | KS collar | 4.000–8.307 | 4.044–8.834 |

The KO-on join K and Θ L2 orders show the same rise through the steep joining layer as the T4C constraints. The far K row shows the acceptance failure; the exterior also has low or negative orders in K, Θ and Π at the finest pair.

| Member | Join K orders, three pairs | Join Θ orders, three pairs | KS-collar K orders, three pairs | Far K orders, three pairs |
|---|---|---|---|---|
| B | 1.758, 3.384, 4.410 | 2.000, 3.157, 3.767 | 3.837, 6.210, 7.069 | 3.649, 0.111, −2.020 |
| E | 2.629, 3.994, 4.645 | 2.652, 3.551, 3.889 | 5.184, 6.718, 7.011 | 3.395, −1.537, −2.007 |
| reference | 2.052, 3.427, 4.074 | 2.140, 3.299, 3.815 | 4.317, 6.224, 6.913 | 3.963, 2.990, −1.841 |

The no-KO far K orders likewise end at `−2.01/−2.01/−1.84`. All floor counts are zero and all setter margins are above 100. The [stage logs](rhs-B-4096.log) include elapsed time, peak RSS and cells at every boundary; the 12 full diagnostic runs used about 200 s in total. The first-puncture-cell audit in those logs records χ, K and α RHS with no order claim; it does not yet cover every evolved group. No h⁻⁴ binned join-RHS profile was saved, so the location of the join **RHS** peak is unresolved; the T4C outer-shoulder localization applies only to constraints.

### Fixed-hierarchy B smoke and cost

The [smoke parameter file](../../Examples/EMS/params-ks2-reference-stationary-B-smoke.txt) uses a base `M/64` grid on `[-8M,8M]×[0,8M]`, with forced fixed radii through the existing mass-extraction tagging parameters: level 1 `M/128` to `4M`, level 2 `M/256` to `1M`, and level 3 `M/512` to `0.6M`. All four masks lie inside the finest region; their measured cell counts are `5274/1274/19984/8608`, identical to the B M/512 single-grid mask counts. Regridding during evolution is disabled. The outer boundary is more than `7M` from the horizon, far outside the causal reach of this `0.046875M` smoke.

The [96-step time series](smoke-B-harmonic.csv) reports H, vector momentum and signed electric-Gauss L2, maximum `|α/α_*−1|`, maximum absolute finite-difference time derivative over all evolved components, and the maximum outgoing radial characteristic over 382 cells in an `h`-wide ring at `0.95r_h`. The characteristic is computed from the evolved `−β^r+α/√γ_rr`. The table shows the first and final finest steps; the file contains all intermediate steps.

| `t/M` | Mask | H L2 | M L2 | G_E L2 | max lapse drift | max component rate |
|---:|---|---:|---:|---:|---:|---:|
| 0.00048828125 | join | 4.20448 | 20.3854 | 0.00110014 | 0.0100824 | 195.360 |
| 0.00048828125 | KS collar | 0.843510 | 1.57905 | 0.000130189 | 0.000631562 | 54.1889 |
| 0.00048828125 | exterior | 5.54468e-7 | 7.92308e-7 | 4.74002e-8 | 1.33894e-9 | 7.68488e-5 |
| 0.00048828125 | far | 1.68994e-9 | 2.83050e-9 | 7.29295e-10 | 5.70655e-14 | 1.35688e-8 |
| 0.046875 | join | 1.10645 | 3.45147 | 0.0135972 | 0.125244 | 3.74419 |
| 0.046875 | KS collar | 1.34567 | 1.12918 | 0.00160684 | 0.0290288 | 2.99186 |
| 0.046875 | exterior | 0.0243871 | 0.00953301 | 3.39657e-5 | 0.000342083 | 0.240568 |
| 0.046875 | far | 6.83722e-9 | 1.79978e-8 | 7.39330e-10 | 5.70033e-12 | 4.57298e-9 |

The sampled outgoing characteristic maximum moves from `−4.93515e-4` to `−4.73377e-4`, always negative. The minimum finest-level lapse/chi floor factors over the smoke were `1155.0/19316.6`; each level's stage scans found zero activations, and the diagnostic nonfinite count stayed zero. The smoke uses `min_lapse=min_chi=1e-8`. It started recording after the first finest step; the production AMR t=0 mask values were not separately written.

The exact initialization-only pilot took `39.71 s`; the one-coarse-step real-path pilot took `57.17 s`; the 12-coarse-step run took `227.059 s`, of which the timed evolution was `182.371 s`. Dividing by 96 finest steps gives `1.900 s/finest step`, or `3891 s/M` at M/512. [Cost extrapolation](smoke-cost.csv) gives **108.1 h per 100M** at M/512 on this serial host. Holding the physical hierarchy and using twice the resolution in 2D would multiply cell count by about four and steps per M by two, giving **864.6 h per 100M** at M/1024; this is an extrapolation, not a measured M/1024 evolution. No cluster run was submitted.

`memory_pressure -Q` showed 39% free before the 12-step dispatch; the four-level state and reference-cache cell counts gave a projected process peak below half of 24 GiB and left more than one fifth free. Peak RSS of that production smoke was not captured. The preceding finest-grid RHS run measured 3.093 GB peak RSS.

### Regression controls and reproduction

The [T5 initialization hash comparison](gauge-regression-results.csv) reran the same 34f2ef0 and current `EMSKS2Regression` binaries for legacy_dat, EMSTRUMPET reference and echo-B, and EMSCTT. All 28 evolved variables including five ghost layers are bit-identical in each paired run. The CTT input in this rerun is `.../.data/binary-ctt-t6/n16-r6.ctt`, so its hash differs from the earlier T4C CTT row, while baseline and current agree with each other. The old `ExperimentalGauge` source branch remains the default and its arithmetic/KO call is unchanged. A few-step baseline-versus-current ExperimentalGauge **evolution** comparison and an EMSKS2 default-initialization full-state hash were not completed, so G4 as a whole remains unverified; there is no observed regression.

Reproduce after setting `CHOMBO_HOME` and `JL` as in the earlier section:

```sh
env CHOMBO_HOME=$CHOMBO_HOME OMP_NUM_THREADS=1 make -C Examples/EMS all DIM=2 -j4
env CHOMBO_HOME=$CHOMBO_HOME OMP_NUM_THREADS=1 make -C Tests/EMSKS2 all DIM=2 -j4
EXE=Tests/EMSKS2/EMSKS2RHSConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex
env OMP_NUM_THREADS=1 "$EXE" B 512 "$JL/artifacts/echo-evolution/t3/B.ks2" > Tests/EMSKS2/rhs-B-512.csv 2> Tests/EMSKS2/rhs-B-512.log
env OMP_NUM_THREADS=1 "$EXE" B 4096 "$JL/artifacts/echo-evolution/t3/B.ks2" Tests/EMSKS2/params-fine-B.txt > Tests/EMSKS2/rhs-B-4096.csv 2> Tests/EMSKS2/rhs-B-4096.log
python3 Tests/EMSKS2/analyze_rhs.py
env OMP_NUM_THREADS=1 Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex Examples/EMS/params-ks2-reference-stationary-B-smoke.txt > Tests/EMSKS2/smoke-B-harmonic-12.log 2>&1
```

For the full RHS ladder, run B/E `N=512,1024,2048,4096` and reference `N=32,64,128,256`; pass `params-fine-MEMBER.txt` only at the largest N. `analyze_rhs.py` returned `results=1440 orders=1080 finest_L2_below_3.5=100 outcome=NOT-STATIONARY`. The smoke command returned exit 0. Builds returned exit 0 in 3.5 s (RHS harness) and 9.4/12.2 s (production stages). The local compute ledger, including RHS, initialization and smoke pilots, and the final smoke, is approximately 9.2 minutes; the optional slicing smoke and end-to-end restart/evolution regression were not launched within the hard local limit. This section is limited to t = 0 stationary-RHS behavior and a `0.046875M` harmonic smoke on this grid, with no claim about long evolutions or gauge stability.

## T5b: evolution controls and fixed-time reference smoke (2026-09-28)

**Registered outcomes:** G4 is **REGRESSION** for the pre-existing ExperimentalGauge evolution path (for example, `emstrumpet1` reference `chi`, level 0, step 4). Its t = 0 state is bit-identical and its production checkpoint restart control is bit-identical. S1 is **NOT-CLOSED / BLOCKED** at M/1024 by the stated memory reserve, so neither GRID-SCALE nor RESOLUTION-INDEPENDENT nor MIXED is assigned. S2 is **GROWS**: on the M/32 reference member, exterior H rises through 1.5 M, although the join lapse drift peaks near 1 M and then declines. The optional M/32 onepluslog run also reaches 1.5 M and grows in exterior H.

### Inputs, edits, and provenance

The supplied Julia contract and B/E/reference `.ks2` files were rehashed and match the complete SHA-256 values in the opening table. [T5b fingerprints](t5b-fingerprints.csv) include those inputs, the earlier B smoke and cost CSVs, and the main new output CSVs. The current production executable is SHA-256 `b5ffe95174e4945a426b16e2041ab4e0ddc3d523646473ac163821d33978e848`; the separately built 34f2ef0 executable is `64a01543c25fee87563120822695cd19d34200fdfd9dd6696ca4f4b3d4a2ca8e`; the one-line T4-initializer comparison build is `22254d8969d3f8ba3738b4c6ce66f99c1a01be50149bfcb1a5945d5de6f5053a`.

The reference branch now passes `m_p.sigma` to its CCZ4 RHS constructor. All four existing reference-stationary example parameter files already explicitly set `sigma = 1`; the new S1/S2/S3 files inherit it. The smoke mask tests are byte-for-byte the same inequalities as the t = 0 RHS harness: `ρ−4h ≥ ρ_a/2`, join `[ρ_a,r_m)`, KS collar `[r_m,r_h)`, exterior `[r_h,2r_h)`, and far `[2r_h,2.3r_h]`. No mask change was needed. The finest-level smoke now also writes `diagnostics.csv.radial.csv` at each finest step: 60 fixed ρ bins from `ρ_a/2` to `3r_h`, with midpoint ρ, sampled areal r and join z, cell count, and RMS/max of `|H|`, `|M_i|`, `|G_E|`, `|Θ|`, `|Ξ|`, `|Λ|`, `|K−K_*|`, and `|α/α_*−1|`. Empty bins are reported with count zero. Existing mask CSV columns remain intact. The [analysis script](analyze_t5b.py) asserts zero nonfinite/floor counts and verifies all four B M/512 smoke mask counts against the t = 0 RHS CSV before writing comparison tables.

[All defaults printed by pout](parameter-defaults.csv) are retained per run, including repetitions across production cases. In reference runs the consequential defaults were `dt_multiplier=0.25`, `num_ghosts=3`, `kappa1=0.1`, `kappa2=0`, `kappa3=1`, `covariantZ4=1`, `formulation=0`, and `nan_check=1`; `reference_onepluslog_n=2` defaulted in harmonic runs and is explicit in S3. No reference run defaulted `sigma`. The pre-existing ExperimentalGauge controls defaulted `sigma=0.1` where their parameter files did not set it. These are the actual pout values, not substituted assumptions.

### G4: full-state evolution, restart, initialization, and join RHS

The separate 34f2ef0 production build ran four steps on exactly matched physics parameter files under [g4-controls](g4-controls/): reference `emstrumpet1` with levels 0–2, echo-B `emstrumpet1`, one boosted EMSCTT binary, and `legacy_dat`. Each paired run used the same parameter-file bytes in separate output directories. [Per-variable SHA-256 comparisons](g4-evolution-hashes.csv) cover all 28 evolved components and all checkpointed active and ghost cells on every level. Every t = 0 component matched. At step 4:

| ExperimentalGauge path | levels | matched component-level rows | largest absolute float difference |
|---|---:|---:|---:|
| trumpet reference | 0–2 | 23 / 84 | 2.90e-14 |
| trumpet echo-B | 0 | 14 / 28 | 1.11e-16 |
| EMSCTT | 0 | 7 / 28 | 1.95e-15 |
| legacy_dat | 0 | 11 / 28 | 5.06e-16 |

These differences occur in active cells as well as ghosts: for reference level 0, 46,180 of 229,376 active component values differ; for echo-B, 43,856 of 229,376. A second 34f2ef0 reference run matched the first at step 4 for all 84 component-level rows ([self-repeat hashes](g4-baseline-repeat-hashes.csv)). Thus the requested **bit identity fails** despite the small numerical differences. Every affected variable and level is recorded in the per-variable CSV; the cause is unresolved. The old gauge arithmetic and matter equations were not adjusted to hide this registered result.

For reference-stationary B, the cheap fixed three-level hierarchy (finest M/64) ran four coarse steps continuously and as two steps plus checkpoint and two restarted steps. [Restart hashes](g4-restart-hashes.csv) match **84/84** evolved component-level rows bit for bit, including checkpointed ghosts. The restart log shows fresh reference cache fills on levels 0, 1, and 2 (`39,200`, `67,712`, and `23,328` cells); the prescribed shift was reimposed. A separate M/64 reference-member two-step versus one-step-plus-restart control also matched **84/84** ([schedule hashes](s2-reference-64-schedule-hashes.csv)), establishing bitwise equivalence for the split S2 schedule.

For EMSKS2 B with `gauge_type=experimental`, a production t = 0 checkpoint was compared against a copied current production build with only the old unconditional `ExperimentalGauge` initializer restored by [this patch](g4-t4-init.patch). [All-variable hashes](g4-ks2-init-hashes.csv) match **26/28**; **B1 and B2 are the only differences**. This is the exact T4 initializer behavior comparison on the same production grid and checkpoint path. The current initializer keeps B1/B2 zero as set by EMSKS2.

The [KO-on B join RHS profiles](join-rhs-B-1024.csv) and [fine profile](join-rhs-B-2048.csv) provide 40 fixed bins for K, Ã, Γ̃, Θ and α at M/1024 and M/2048, with RMS and max multiplied by `h⁻⁴`. All five M/2048 RMS peaks occur in bin 38 at `ρ=0.8715406 r_h`, `z=0.751339`: K `2.2993e10`, Ã `2.5878e10`, Γ̃ `1.0811e9`, Θ `3.3322e7`, α `3.0213e6`. The matching M/1024 K peak is `1.4258e10`; no overlay or stationarity claim is inferred from this pair. The raw [norm CSVs](join-rhs-B-2048-norms.csv) and stage logs retain the unchanged full t = 0 RHS controls. These two RHS programs returned 0 in 3.46 and 11.56 s.

### S1: B fixed-time refinement, stopped before M/1024

The B pilots used the T5 fixed physical hierarchy `[-8M,8M]×[0,8M]`, refinement radii `4M`, `1M`, `0.6M`, and finest spacings M/256 and M/512. The puncture lies on the shared vertex at both resolutions. Their two-coarse-step programs returned 0 in 19.038 s (0.603 GB peak RSS) and 76.290 s (2.733 GB peak RSS). [Time series](s1-time-series.csv), [all common-time first-pair orders](s1-first-pair-orders.csv), and [radial snapshots at `t=M/128`](s1-radial-common.csv) preserve the controls. The table shows H/M/`G_E` L2 as **M/256 / M/512**:

| t/M | mask | H L2 | M L2 | G_E L2 |
|---:|---|---:|---:|---:|
| 0.00390625 | join | 3.611e+00 / 1.677e+00 | 2.241e+01 / 7.624e+00 | 5.171e-03 / 1.770e-03 |
| 0.00390625 | KS collar | 2.565e+00 / 9.353e-01 | 6.177e+00 / 1.959e+00 | 2.780e-03 / 5.156e-04 |
| 0.00390625 | exterior | 1.491e-02 / 1.616e-04 | 2.925e-02 / 1.785e-04 | 1.276e-05 / 7.460e-08 |
| 0.00390625 | far | 2.922e-08 / 2.299e-09 | 4.803e-08 / 4.434e-09 | 1.164e-08 / 7.301e-10 |
| 0.00781250 | join | 2.220e+00 / 1.134e+00 | 1.510e+01 / 5.019e+00 | 9.470e-03 / 2.487e-03 |
| 0.00781250 | KS collar | 1.895e+00 / 1.274e+00 | 7.424e+00 / 1.875e+00 | 3.459e-03 / 8.125e-04 |
| 0.00781250 | exterior | 1.733e-02 / 6.960e-04 | 5.732e-02 / 5.735e-04 | 4.667e-05 / 4.974e-07 |
| 0.00781250 | far | 3.295e-08 / 3.026e-09 | 5.738e-08 / 6.534e-09 | 1.165e-08 / 7.309e-10 |

The first-pair exterior H orders are 6.528 and 4.638 at those times; corresponding M orders are 7.356 and 6.643, and `G_E` orders 7.418 and 6.552. At `t=M/128`, both radial H profiles peak immediately outside the horizon in bin 19, `ρ=1.0226r_h`: RMS `9.694e-2` at M/256 and `4.567e-3` at M/512. At about `1.95r_h` the values are `5.04e-8` and `4.26e-9`, localizing the measured exterior contamination near `r_h`.

The requested M/1024 program was **not launched**; its [exact parameter file](s1-B-1024-blocked/params.txt) is saved for review. Immediately after the M/512 pilot, `memory_pressure -Q` reported 52% free of 24 GiB. Fourfold cell scaling projects at least 10.93 GB from its measured 2.73 GB peak, while preserving 20% free allows only 8.25 GB for one dispatch ([resource decision](s1-resource-decision.csv)); the earlier M/256-to-M/512 measured RSS ratio was 4.53. The M/512 pilot's 76.290 s wall included 32.775 s evolution; scaling setup by four and each coarse step by eight projects 436 s for two M/1024 coarse steps or 567 s for three, both before any unmeasured cache penalty. The exact M/1024 hierarchy therefore violates the hard memory reserve even though this time estimate is below ten minutes. Its measured runtime is absent because memory preflight refused the run. No finest-pair order, longest common three-resolution time, or S1 registered numerical classification is claimed. Re-cutting onto an authorized larger-memory host is the concrete next route; it is not an authorized local or cluster submission in this tranche.

### S2 and S3: reference member through and past lapse-lock time

The reference member used `M=0.5710876672204076`, a base domain `[-16M,16M]×[0,16M]`, forced refinement radii `12M` and `9M`, and finest coverage past `5r_h=8.75M`. The final checkpoint boxes cover all `123,296`/`493,168` finest cell centres inside `5r_h` at M/32/M/64, with zero missing ([coverage audit](s2-coverage.csv)); the outer boundary lies `7.25M` beyond this radius. M/32 ran to `1.5M` (48 coarse steps) and M/64 to `0.5M` (32 coarse steps). Both S2 resolutions used the harmonic equation, `sigma=1`, and existing CCZ4/matter parameters; S3 differs only by the authorized `reference_f=onepluslog`, `reference_onepluslog_n=2`. All recorded floor and nonfinite counts are zero, and sampled outgoing characteristics at `0.95r_h` remain negative. [Full mask time series](s2-s3-time-series.csv), [common-time ratios](s2-common-ratios.csv), and [60-bin radial snapshots](s2-s3-radial-snapshots.csv) preserve the measured comparisons; complete every-step profiles remain in each final stage directory.

| t/M | exterior H M/32 harmonic | exterior H M/64 harmonic | H ratio M/32 : M/64 | join lapse drift M/32 | join lapse drift M/64 |
|---:|---:|---:|---:|---:|---:|
| 0.125 | 1.9306e-2 | 6.3904e-5 | 302.1 | 5.807% | 1.681% |
| 0.25 | 7.4538e-2 | 7.0973e-4 | 105.0 | 7.034% | 2.369% |
| 0.5 | 1.5060e-1 | 2.8971e-3 | 51.98 | 8.806% | 3.027% |
| 1.0 | 2.6204e-1 | — | — | 10.855% | — |
| 1.5 | 3.9492e-1 | — | — | 9.802% | — |

The M/32 exterior H secant slope over `[1M,1.5M]` is `+0.2658` per M; exterior M and `G_E` also increase there, from `0.1556` to `0.2081` and `1.11e-5` to `2.13e-5`. Join lapse drift peaks near `1M` and falls by 1.5M, so its lock does not settle exterior H on this grid. The M/64 exterior H also rises through its available 0.5M, while its H is 52 times smaller at that common time. This is **GROWS**, with a strong resolution dependence and no claim about the finer run after 0.5M. The final M/32 exterior radial peak is bin 19 at `ρ=1.0271r_h`, H RMS `1.923` (M/64 at 0.5M: `0.01344` in the same bin); at `ρ=1.953r_h`, M/32 final H RMS is only `2.92e-8`. The growth is attached to the horizon-side exterior, not the far mask. No floor, barrier crossing, or nonfinite event occurred by either endpoint.

The optional onepluslog M/32 run reached 1.5M with zero floor/nonfinite counts. [Comparison rows](s3-comparison.csv) show exterior H `0.05910`, `0.11362`, `0.16091` at 0.5, 1.0, 1.5M, versus harmonic `0.15060`, `0.26204`, `0.39492`. Its join lapse drift is `13.04%`, `16.04%`, `20.71%` at those times, versus harmonic `8.81%`, `10.86%`, `9.80%`. Onepluslog reduces H on this coarse grid but still grows; no gauge was selected or tuned from this comparison. Its exterior radial H peak at 1.5M is also bin 19 (`0.8152` RMS).

### Commands and resource ledger

With `CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib` and `OMP_NUM_THREADS=1`, the current EMS production and RHS harness builds returned 0 in 10.9 and 4.6 s. The separate 34f2ef0 production archive build returned 0 in 12.3 s. The T4 comparison copy build returned 0 in 16.8 s after applying [the saved one-line patch](g4-t4-init.patch). Reproduce the builds and a representative real-path program from this worktree root with:

```sh
env CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib OMP_NUM_THREADS=1 make -C Examples/EMS all DIM=2 -j4
env CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib OMP_NUM_THREADS=1 make -C Tests/EMSKS2 all DIM=2 -j4
mkdir -p /private/tmp/emsks2-t5b-baseline
git archive 34f2ef0 | tar -x -C /private/tmp/emsks2-t5b-baseline
env CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib OMP_NUM_THREADS=1 make -C /private/tmp/emsks2-t5b-baseline/Examples/EMS all DIM=2 -j4
python3 Tests/EMSKS2/run_g4_controls.py
env OMP_NUM_THREADS=1 Tests/EMSKS2/EMSKS2RHSConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex B 2048 /Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-LSb44u/artifacts/echo-evolution/t3/B.ks2 - Tests/EMSKS2/join-rhs-B-2048.csv > Tests/EMSKS2/join-rhs-B-2048-norms.csv 2> Tests/EMSKS2/join-rhs-B-2048.log
python3 Tests/EMSKS2/run_stage.py Tests/EMSKS2/s1-B-512-pilot B-M512-pilot 2 3
python3 Tests/EMSKS2/run_stage.py Tests/EMSKS2/s2-reference-32-stage48 reference-M32-to-1.5M 1 3
/Users/auroradysis/.cache/ariadne/venvs/T601/bin/python Tests/EMSKS2/compare_state.py Tests/EMSKS2/g4-restart/continuous/chk/EMS_000004.2d.hdf5 Tests/EMSKS2/g4-restart/split-restart/chk/EMS_000004.2d.hdf5 reference-B-restart
/Users/auroradysis/.cache/ariadne/venvs/T601/bin/python Tests/EMSKS2/check_coverage.py
python3 Tests/EMSKS2/analyze_t5b.py
```

Each stage directory contains its exact `params.txt`, `run.log`, checkpoint, time-series CSV, and every-step radial CSV. [The cost ledger](t5b-cost.csv) records each measured wall time, peak process RSS and exit status; its 11 simulation stages total `1,690.081 s` (28.17 minutes), excluding short builds, G4 controls, and the two RHS runs, so the local tranche remains below approximately 45 minutes. The longest stage was M/64 step 2→17 at `461.599 s`, `3.599 GB` peak RSS; the next was step 17→32 at `421.166 s`, `3.809 GB`. The onepluslog continuation took `261.122 s`, `1.070 GB`; harmonic M/32 continuation took `259.639 s`, `1.079 GB`. Every program stayed below ten minutes and every dispatch below half of 24 GiB. The M/1024 B memory refusal is separate from the completed cost ledger.

The T4 initialization comparison build was made by copying `Source/` and `Examples/EMS/` into `/private/tmp/emsks2-t5b-t4/` (excluding `initial_data/`, `o/`, `d/`, and executables), applying `g4-t4-init.patch` there, and running the same `make -C .../Examples/EMS all DIM=2 -j4` command. The four G4 checkpoint pairs and the separate B restart checkpoints can each be re-compared with `compare_state.py` using the corresponding files under their stage directories. The staged S2 runs take the prior stage's named checkpoint and both diagnostic CSVs as inputs; those copied inputs remain alongside the final outputs. The measured M/64 restart-equivalence control is in `s2-reference-64-one/` and `s2-reference-64-restart1/`.

This tranche establishes a bitwise failure in pre-existing evolution despite matching t = 0 states; a bitwise-safe reference restart on three levels; the exact EMSKS2 B initializer delta; two-resolution early B refinement evidence without the required finest pair; and reference-member growth through 1.5M at M/32. It makes no continuum-mode claim for B, no finer-grid stability claim past 0.5M for the reference member, and no change to the gauge, taper, floors, CCZ4 or matter equations. Binary/boost EMSKS2 evolution and cluster submission remain archive-only.

## T5c: ExperimentalGauge evolution contraction test (2026-09-28)

**Registered outcome: CONTRACTION.** The two builds of `34f2ef0b307428209c92b8e3b33d6f8df02ed45c` and this worktree, compiled with identical production settings plus `cxxoptflags='-O3 -ffp-contract=off'`, match at step 4 for all **168/168** evolved component-level checkpoint hashes. The original production-flag comparison was 55/168. Thus the observed difference depends on floating-point contraction under the production compilation, without evidence of a changed arithmetic expression on the old path. D2 and the T4 comparison build were not run because D1 met its registered bit-identity predicate.

The exact input parameter files under `g4-controls/CASE/{baseline,current}/params.txt` were byte-identical within every pair and were copied unchanged into the scratch runs. The four data inputs were SHA-256 `6f0820a576620f1f7230c701131af56a2312b130e31d5de2a752c32cefe6d24f` (reference trumpet), `d54cc1141dc90536dbb6fda07d408b8b143294791269f65dd5b270b9dac91251` (echo-B trumpet), `a875dd55cf628efe7db59fc2e2eb867c8794434fd4e31bd00324f544fec844c1` (CTT), and `77f43e34138acf9ed20f9ea03f8bbea12bda59e81bbe40babed0885f29a32845` (legacy data). The earlier `g4-evolution-hashes.csv` rehashed to its recorded `0363024c…92002e38`. The D1 binaries rehashed to `255a0d64…3b595cc216b50de43cb95defb43c` and `446b9054…9a9e6a6e33a5bdb1aaf9e51`; clean repeat builds produced the same executable hashes.

For every level and evolved component, the unchanged `compare_state.py` computes `SHA256` over the ordered raw bytes of each checkpoint box, including active and ghost cells. It also asserts equal iteration, time, boxes, offsets, component names and cell counts. The acceptance predicate is equality of **every** component hash; no float tolerance was introduced. [D1 hashes](t5c-d1-hashes.csv) include both step 0 and step 4:

| G4a case | levels | step 0, contraction off | step 4, contraction off | step 4, original production flags |
|---|---:|---:|---:|---:|
| trumpet reference | 0–2 | 84/84 | 84/84 | 23/84 |
| trumpet echo-B | 0 | 28/28 | 28/28 | 14/28 |
| EMSCTT | 0 | 28/28 | 28/28 | 7/28 |
| legacy_dat | 0 | 28/28 | 28/28 | 11/28 |
| **all** | | **168/168** | **168/168** | **55/168** |

D3 tested an out-of-line split of the reference-stationary member functions and the reference branch of `specificEvalRHS` in a separate scratch build at the **unchanged production flags** (`-O3 -std=c++17 -fopenmp`). The exact trial is archived as [a patch](t5c-d3-trial.patch), which passes `git apply --check` but is **not applied**. Against the 34f2ef0 production executable, its [G4a hashes](t5c-d3-hashes.csv) returned to 168/168 at steps 0 and 4. The required reference-stationary preservation control failed: relative to the pre-refactor production build, only 26/84 hashes matched after two steps and 25/84 after four. The refactor was therefore rejected. Its internal continuous-versus-restart comparison remained 84/84; that narrower result does not satisfy preservation against the pre-refactor build.

| D3 control | matched component-level hashes |
|---|---:|
| G4a old paths, step 0 / step 4 | 168/168 / 168/168 |
| EMSKS2 B experimental initialization, pre-refactor vs trial | 28/28 |
| B reference-stationary initialization, pre-refactor vs trial | 84/84 |
| B reference-stationary continuous step 2, pre-refactor vs trial | **26/84** |
| B reference-stationary continuous step 4, pre-refactor vs trial | **25/84** |
| B reference-stationary trial continuous vs trial restart, step 4 | 84/84 |
| B reference-stationary restart, pre-refactor vs trial, step 4 | **25/84** |

The full per-component controls are in [the initialization hashes](t5c-d3-init-hashes.csv) and [the reference hashes](t5c-d3-reference-hashes.csv). The trial fixture executable returned 0: 730 comparisons for each of B, E and reference, all seven malformed fixtures rejected, with the same three value-error maxima as the saved fixture results ([fixture log](t5c-d3-fixture.log)). The trial continuous, split-first and split-restart programs returned 0 in 6.75, 4.07 and 2.15 seconds. Its EMSKS2 B t = 0 program returned 0 in 0.52 seconds. D3 build time was 19.93 seconds.

| Output | SHA-256 |
|---|---|
| `t5c-d1-hashes.csv` | `27286987f517b2a1246e588324ff689afbcc9cdf1f397e800b64b22a4c013026` |
| `t5c-d3-hashes.csv` | `b4d8637723a414f82981097bc08c03f30ef99c83c19267fbeff2b038f2e6f237` |
| `t5c-d3-init-hashes.csv` | `f983c9fdfd281b3b1ed1c9565cbba0e4cb526b4a24900a287d0be8cd97c5561a` |
| `t5c-d3-reference-hashes.csv` | `7719499a319d1b522e938ccd22d1dacf4f4afb616580137511b2dace42d7da62` |
| `t5c-d3-trial.patch` | `9b12e81b08d06deba17053bfbf71aa83ac0d0e86e26d7f30c52541a1e589c4fc` |

### Reproduction and cost

The D1 build trees are `/private/tmp/emsks2-t5c-d1-baseline-measured` (fresh `git archive 34f2ef0`) and `/private/tmp/emsks2-t5c-d1-current-measured` (copy of this worktree's `Source/` and `Examples/EMS/`, excluding build outputs and `initial_data/`). The excluded data are read from the unchanged absolute paths in the parameter files. Both use the same installed Chombo/HDF5 libraries and `CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib`. The compile **and link** command lines in each build log contain `-O3 -ffp-contract=off -std=c++17 -fopenmp`. The builds returned 0 in **14.60** and **17.97** seconds, respectively. The D1 G4a bundle returned 0 in **10.29** seconds (individual stages: 4.123, 4.246, 0.230, 0.219, 0.408, 0.414, 0.330, 0.335 seconds in [stage order](t5c-d1-stages.log)). Every program was estimated from the prior G4a controls at under five minutes and completed far below the ten-minute per-program cap. `memory_pressure -Q` reported 32% free before builds and 67% free before the reference control; no peak RSS was recorded for these short programs.

```sh
CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib OMP_NUM_THREADS=1 make -C /private/tmp/emsks2-t5c-d1-baseline-measured/Examples/EMS all DIM=2 -j4 cxxoptflags='-O3 -ffp-contract=off'
CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib OMP_NUM_THREADS=1 make -C /private/tmp/emsks2-t5c-d1-current-measured/Examples/EMS all DIM=2 -j4 cxxoptflags='-O3 -ffp-contract=off'
python3 Tests/EMSKS2/run_g4_controls.py --output /private/tmp/emsks2-t5c-d1-controls --baseline /private/tmp/emsks2-t5c-d1-baseline-measured/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex --current /private/tmp/emsks2-t5c-d1-current-measured/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex
/Users/auroradysis/.cache/ariadne/venvs/T601/bin/python Tests/EMSKS2/collect_g4_hashes.py /private/tmp/emsks2-t5c-d1-controls > Tests/EMSKS2/t5c-d1-hashes.csv
```

The optional `--output`, `--baseline` and `--current` arguments added to `run_g4_controls.py` only redirect the existing eight runs and copy the exact parameter-file bytes into fresh output directories; its default behavior remains the original G4a bundle. `collect_g4_hashes.py` invokes the unchanged `compare_state.py` for the same cases at steps 0 and 4. Regenerating both D1 and D3 G4a CSVs reproduced their saved bytes exactly. The D3 bundle took **11.21** seconds ([stage log](t5c-d3-stages.log)); its old-path table is all 84/84, 28/28, 28/28 and 28/28 at each step. No production flags, gauge, CCZ4, matter, data or acceptance criterion was changed in this worktree by T5c. The conclusion is limited to these G4a paths, four coarse steps and compiler/build environment; the failed D3 reference control forbids keeping the trial split.

## T5d: adopted split, MPI diagnostics, and cluster preflight (2026-09-28)

**Registered outcome: READY-EXCEPT.** P1 and the serial and source-level parts of P2 pass. P3 has four completed local pilots and two memory refusals. The MPI smoke is required on the cluster because no local `mpicxx` exists. P4's requested archive move is blocked by the filesystem sandbox, so the worktree is not yet output-free or ready to commit as-is.

The controller adopted the [T5c D3 patch](t5c-d3-trial.patch): the reference-stationary implementation is now in `Examples/EMS/EMSBH2DLevelReference.cpp`. The production `-O3 -std=c++17 -fopenmp` example build returned 0 in **13.52 s**. [G4a hashes](t5d-g4-hashes.csv) match 34f2ef0 at step 0 and step 4 on all four old paths:

| Path | Levels | step 0 | step 4 |
|---|---:|---:|---:|
| Trumpet reference | 0–2 | 84/84 | 84/84 |
| Trumpet echo-B | 0 | 28/28 | 28/28 |
| EMSCTT | 0 | 28/28 | 28/28 |
| Legacy data | 0 | 28/28 | 28/28 |
| **Total** | | **168/168** | **168/168** |

The eight G4a stages returned 0 in **8.74 s** combined (3.135, 3.722, 0.205, 0.203, 0.376, 0.477, 0.297, 0.329 s; [stage log](t5d-g4-stages.log)). The fixture compile/run returned 0 in **2.04/0.34 s**, with 730 comparisons per member and all seven malformed fixtures rejected ([fixture log](t5d-fixture.log)). The EMSKS2 B experimental initializer returned 0 in **0.62 s**; [T4 comparison hashes](t5d-init-hashes.csv) match **26/28**, differing only in `B1` and `B2`. A four-step B reference-stationary continuous run returned 0 in **6.25 s**; split-first and split-restart returned 0 in **5.19/2.64 s**, and [restart hashes](t5d-restart-hashes.csv) match **84/84**.

The pre-split and post-split sources were also built separately with `cxxoptflags='-O3 -ffp-contract=off'` (exit 0 in **15.58/17.82 s**). Both ran the same four-coarse-step B reference-stationary parameters (exit 0 in **7.16/6.98 s**). [All 84 evolved component-level hashes](t5d-reference-contract-off-hashes.csv) match bit for bit, including ghosts. This is the controller's contraction-disabled implementation-equivalence control. It does not claim production-flag bit identity for the new reference path.

### Diagnostic reduction and cadence

`write_reference_diagnostics` accumulates locally over each rank's boxes, then under `CH_MPI` reduces mask and radial squared sums and counts with `MPI_SUM`, maxima (including barrier speed) with `MPI_MAX`, and lapse/χ margins with `MPI_MIN` on `Chombo_MPI::comm`. Only rank 0 opens either CSV. Empty masks are rejected after the global count reduction. The serial branch performs no reduction and, with both new intervals at their zero defaults, retains the old every-finest-step schedule. The serial time-series and radial CSVs from the same four-step B run are **byte-identical** to the saved pre-MPI D3 output: SHA-256 `64f7e4637e85120165c3e5b17c350b545f1a3c7dbaa88c42a7bacbf96d179b9c` and `39cc01e34ed05eab1d2133a2069eda31bd891074239982c5408f257d7fcc1e43`.

`reference_diagnostics_interval` gates both outputs; `reference_radial_interval` further selects radial rows. Zero means every finest step for either interval. Output is written at the first finest step crossing each requested time interval. An M/32 reference run with `max_steps=4` and intervals `0.01M`/`0.1M` returned 0 in **33.48 s**: 48 mask rows at 12 distinct times (`t/M=0.015625` through `0.125`), and 60 radial rows at `t/M=0.1015625`. No MPI compiler is installed locally (`mpicxx` was absent), so the `CH_MPI` branch has only source review and the serial control, not an MPI runtime pass.

**Required cluster MPI smoke before production:** stage the same `reference.ks2` (SHA-256 `fb46ceec…70533`) beside `params-cluster-ref-M32.txt`; shorten only `max_steps` to 4; run once with one MPI rank and once with 32 ranks in separate output directories. Require exit 0, exactly one header per file, the same 12 time-series timestamps and 60 radial rows, identical time/level/mask/bin/count and barrier-sample keys, and finite fields. Compare sums-derived RMS columns within reduction-order round-off (for example `64 ε × (|a|+|b|+1e-30)`); maxima and minima should agree at the same precision, with any larger difference investigated. Repeat continuous four steps versus two plus restart at 32 ranks, carrying the checkpoint and both existing diagnostic CSVs into the restart directory, and require bitwise checkpoint identity at fixed rank count. This is a required cluster gate; no cluster access was attempted here.

### Six cluster inputs and local cost

The six [parameter files](../../Examples/EMS/) use the unchanged T5b S2 reference hierarchy or T5 B smoke hierarchy. Every file explicitly sets `sigma=1`, `kappa1/2/3`, `covariantZ4`, `dt_multiplier=0.25`, `formulation=0`, all six MovingPunctureGauge parameters, `reference_onepluslog_n=2`, `gauge_type=reference_stationary`, `reference_f=harmonic`, and `min_chi=min_lapse=1e-8`. Only `ref-M64-k1` changes κ₁, to `1.3760187883491917 = 0.1 × 7.8582735988979167 / 0.5710876672204076`. Each requests `stop_time=10M`, diagnostic intervals `0.01M`/`0.1M`, and `plot_interval=-1`; the radial interval is an additional selector on the diagnostic interval. Stage the supplied byte-identical `reference.ks2` or `B.ks2` beside the parameter file at run time; the relative `ems_data_path` is intentional.

The puncture is the common vertex at the centre of the x grid and at y=0 on every level (all base N values are even, refinement ratio 2). The nearest non-cartoon outer boundary is **16M** from the reference hole and **8M** from B, giving nominal light round-trip times **32M** and **16M**; their margins beyond the 10M target are 22M and 6M. The hierarchy is fixed (`regrid_interval=0`). Checkpoint intervals below use a **16× effective speedup estimate on 32 ranks**, pending the MPI smoke and first cluster timing; their projected spacing is about 0.04, 0.33, 1.31, 0.33, 0.87 and 0.87 wall-hours in table order. Cases whose full projected wall time is shorter than an hour checkpoint at the end.

| Case | finest h | init RSS GB | one-step peak GB | measured s/finest step | projected serial core-h to 10M | checkpoint coarse steps | local pilot |
|---|---:|---:|---:|---:|---:|---:|---|
| ref-M32 | M/32 | 0.212 | 0.778 | 1.085 | 0.65 | 320 | exit 0, 4.43 + 8.77 s |
| ref-M64 | M/64 | 0.736 | 2.992 | 4.490 | 5.24 | 640 | exit 0, 18.33 + 36.29 s |
| ref-M128 | M/128 | 2.945 projected | **14.096 projected** | 29.426 projected | 41.87 | 640 | memory refused |
| ref-M64-k1 | M/64 | 0.774 | 3.057 | 4.556 | 5.24 | 640 | exit 0, 18.05 + 36.28 s |
| B-M512 | M/512 | 1.100 | 1.878 | 2.455 | 13.98 | 2560 | exit 0, 43.00 + 62.64 s |
| B-M1024 | M/1024 | 4.401 projected | **10.931 projected** | 9.819 projected | 111.77 | 640 | memory refused |

The first time in each completed pilot row is a separate `max_steps=0` initialization, the second is initialization plus one coarse step; pilot inputs changed only `max_steps` from the then-current cluster files. A redundant second `sigma=1` line was subsequently removed from all six final files. Repeating the M/32 one-step pilot with the final file returned 0 in **12.01 s** and produced a byte-identical diagnostic CSV, confirming that the cleanup did not change the parsed physics. Measured evolution seconds per finest step equal `(one-step wall − init wall)/2^max_level`. [The cost CSV](t5d-cluster-cost.csv) has the exact values. For 10M, the projection uses the larger of recent long-stage and one-step seconds per coarse step (M/32 additionally uses its new four-step interval test), scaling the unrun doubled-resolution cases by four for 2D cell count; `serial core-hours = (init + steps × projected seconds per coarse step)/3600`. This is a capacity estimate, not a measured 32-rank speedup. Before either high-resolution pilot, `memory_pressure -Q` reported **50% free** of 24 GiB. Leaving 20% free allowed only 7.73 GB of new process peak; prior two-step peaks scale to 14.10 GB (ref-M128) and 10.93 GB (B-M1024). Neither was launched. The projected init plus one-step walls, before a memory cache penalty, are 145 and 251 seconds; the refusal is memory, not time.

Every pout `not found` parameter and its printed value for each of the four run cases is retained in [parameter defaults](t5d-parameter-defaults.csv): 38 distinct keys per case, including repeated parser occurrences. None is among the explicitly required CCZ4, gauge, timestep or floor keys. There is no pout for the two memory-refused cases, so defaults for those cases are not reported as measurements.

### Commit and archive boundary

[The archive inventory](t5d-archive-plan.csv) lists **206 run-output files, 6,736,539,218 bytes**, with their worktree-relative paths and the requested destination under `/Users/auroradysis/Workspace/EMS/.data/echo-evolution/fork-tests/`. The destination could not be created from this sandbox: `mkdir` returned `Operation not permitted`. No file was moved and the existing README links still point to the worktree copies. Production and test executables plus `o/` and `d/` build directories are left alone as requested, and are excluded from the commit manifest. [COMMIT-MANIFEST.txt](COMMIT-MANIFEST.txt) names only modified or new sources, headers, tests, scripts, READMEs, parameter files and CSVs below 1 MB, with sizes; it does not claim that the worktree is output-free. The archive move and MPI smoke are the remaining release gates.

The manifest contains **167 paths**, including all five modified tracked files and the new separate translation unit. A fresh `git archive 34f2ef0` plus exactly the final source and parameter paths in the manifest built the production example and all four `Tests/EMSKS2/GNUmakefile` executables (exit 0 in **16.85/6.82 s**). The standalone `EMSKS2FinerFeasibility.cpp` compiled in **1.66 s**, and all manifest Python scripts passed `py_compile`. The worktree build artifacts deliberately left untouched and excluded from the manifest are `Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex`, `Tests/EMSKS2/{EMSKS2FixtureTest,EMSKS2GridConvergence,EMSKS2RHSConvergence,EMSKS2Regression}2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex`, `Examples/EMS/{o,d}`, `Tests/EMSKS2/{o,d}`, and `Tests/EMSKS2/__pycache__`. The clean-build scratch tree is `/private/tmp/emsks2-t5d-manifest-tree-final`.

Reproduce the controls from the repository root with `CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib` and `OMP_NUM_THREADS=1`:

```sh
make -C Examples/EMS all DIM=2 -j4
python3 Tests/EMSKS2/run_g4_controls.py --output /private/tmp/emsks2-t5d-g4 --baseline /private/tmp/emsks2-t5b-baseline/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex --current "$PWD/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex"
/Users/auroradysis/.cache/ariadne/venvs/T601/bin/python Tests/EMSKS2/collect_g4_hashes.py /private/tmp/emsks2-t5d-g4 > Tests/EMSKS2/t5d-g4-hashes.csv
python3 Tests/EMSKS2/run_stage.py /private/tmp/emsks2-t5d-pilots/ref-M32/one ref-M32-one 1 1
```

The cluster parameter validation checked explicit required keys, common puncture vertices, 10M stop times, exact coarse-step counts, diagnostic intervals and boundary distances on all six files. The manifest overlay also checked every listed file size before copying. The MPI build/runtime test and requested archive move remain outstanding; neither is inferred from the serial controls or clean build.
