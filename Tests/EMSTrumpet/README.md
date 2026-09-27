# EMSTRUMPET 1 initial data checks

## Direct CCZ4 initialization

The setter and fixture API now share one ADM-to-CCZ4 conversion. Density
assembly returns `gamma_left + gamma_right - I`; CTT assembly returns
`psi^4 * gamma_bar`. The latter no longer rounds the diagonal through
`(gamma + 1) - 1`. The shared converter retains production's
`pow(det, -1./3.)` and initial lapse `sqrt(chi)`; the former fixture-only
converter used `1/cbrt(det)`.

The reference CCZ4 assertion covers one table containing all three rapidities,
not separate boosted/unboosted assertions. Its fingerprint changes from
`c30cf22126f14d61` to `e68444ed5f33ba7b`. The density CCZ4 table changes from
`3d5a69d4cda91686` to `024d0092dba552bc` and now has an explicit assertion.
The 3D-object assertion remains `e1107d452abc4e4b`. All nine radial/jet/object
table fingerprints, every Julia fixture, and all input SHA hashes are unchanged.
Echo B/E CCZ4 helper fingerprints were printed, never asserted: they change
from `8cfc83d44aa754d9` / `59fb7d5a88e073fd` to
`0d83df086db272d1` / `bffc958d04850e8f` solely with the shared converter's
floating-point evaluation. These helper tables must not be confused with
the production-grid comparison: the unboosted reference and unboosted echo-B
grids remain bit-identical in every evolved variable and every ghost cell.

The old/new production comparison uses b397e42 as its baseline. Maximum
absolute primitive differences are `1.1102e-16` for the boosted single,
`4.3368e-19` for density, and `4.4409e-16` for CTT. Their largest differences
relative to each field's maximum magnitude are `3.1685e-15`, `2.5419e-17`,
and `4.8805e-15`, respectively. Shift, scalar, electromagnetic, Theta, Lambda
and Xi fields are bit-identical in all cases.

The controller accepted the initial STOP: applying a field-relative gate to
finite-difference-derived Gamma omitted amplification of metric round-off
by `1/dx`. Gauge B inherits Gamma through `B = 0.75 Gamma - eta shift`.
The corrected gate is `1e-14` relative to the field scale for primitive fields,
and `100 * epsilon * max(abs(h)) / dx` absolute for Gamma and gauge B, with
`epsilon = 2.220446049250313e-16`. All cells, including ghosts, pass:

| Case | dx | Max primitive scaled difference | FD allowance | Max Gamma difference | Max gauge B difference |
|---|---:|---:|---:|---:|---:|
| Reference, unboosted | 1/32 | 0 | 7.1054e-13 | 0 | 0 |
| Reference, boosted | 1/32 | 3.1685e-15 | 7.5153e-13 | 0 | 0 |
| Density binary | 0.3125 | 2.5419e-17 | 7.5082e-14 | 0 | 0 |
| CTT binary / T7 checkpoint | 0.3125 | 4.8805e-15 | 7.5082e-14 | 1.7090e-15 | 1.2820e-15 |
| Echo-B, unboosted | 1/128 | 0 | 2.8422e-12 | 0 | 0 |

All 28 fixture/grid executions pass, including both echo fixture sets, all
parser/SHA controls and all 26 grid runs. Fine-pair H/M L2 orders are
`3.992251/3.979856` for the unboosted reference,
`3.982369/3.977222` for either boost sign,
`3.209941/4.638946` for echo-B, and `3.995777/3.993403` on the CTT exterior
collar. Reference and echo orders match the baseline to the shown precision;
CTT collar H/M orders change by `8.4e-7` / `3.3e-10`. At the
inter-hole Hamiltonian floor the finest order changes from `1.801698` to
`1.833311`; its L2 changes from `5.67227e-11` to `5.54869e-11`. Density H/M
retain their continuum residuals. No numerical fixture tolerance is changed.
The resumed report, all 116 mask comparisons, raw logs and full-state dumps
are retained under `/private/tmp/ems-fixsup-92881/`.

## Echo B and E with nonunit compactification scale

The supplied B and E files were authenticated with `EMS.read_trumpet` and
copied byte for byte into `fixtures-echo/{B,E}/reference.trumpet`. Their SHA-256
hashes are `d54cc1141dc90536dbb6fda07d408b8b143294791269f65dd5b270b9dac91251`
and `7408f854a8d9e390c9db4e0e77400bdfa98a6d622fd20233c8be3f0ba4ace809`.
Each manifest records the source and table hashes, the Julia reader and original
fixture generator hashes, the radii, and the rapidities. The fixture script uses
the same `trumpet_jet`, `trumpet_sample`, `trumpet_object`, and `trumpet_ccz4`
calls and table columns as `scripts/emstrumpet_fixtures.jl`. It samples the
cylinder approach, half and twice the isotropic horizon radius, and `R=0.16M`,
where B has areal radius `r=0.99396M`. Objects and single CCZ4 fields cover
masses 1 and 2 and rapidities `0, ±0.20273255`.

```sh
julia --project=/Users/auroradysis/Workspace/EMS Tests/EMSTrumpet/echo_fixtures.jl
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
make -C Tests/EMSTrumpet all DIM=2 -j 12
OMP_NUM_THREADS=1 Tests/EMSTrumpet/EMSTrumpetFixtureTest2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
  Tests/EMSTrumpet/fixtures Tests/EMSTrumpet/fixtures-echo/B Tests/EMSTrumpet/fixtures-echo/E
```

The C++ test first validates each header and takes the coupling parameters from
that file. It runs the original reference binary table, bit fingerprints, and
nine malformed-file rejections. For the echo rows, every value is checked with
`|C++ - Julia|/(1 + |Julia|)`. The maximum in each tolerance class is:

| Member | Table | Values, limit `1e-12` | First, limit `1e-10` | Second/higher jets, limit `1e-8` |
| --- | --- | ---: | ---: | ---: |
| B | jets | 6.83e-17 | 4.72e-15 | 6.38e-14 |
| B | radial | 4.10e-15 | 4.55e-15 | — |
| B | objects | 1.50e-14 | 2.63e-14 | — |
| B | single CCZ4 | 4.48e-15 | 2.84e-15 | — |
| E | jets | 2.71e-16 | 2.36e-15 | 8.22e-15 |
| E | radial | 1.25e-15 | 1.26e-15 | — |
| E | objects | 1.78e-15 | 2.69e-15 | — |
| E | single CCZ4 | 2.78e-15 | 2.34e-15 | — |

The reference still passes: maximum normalized errors are `2.51e-16` (jets),
`7.84e-16` (radial), `8.89e-16` (objects), `1.15e-15` (single CCZ4), and
`1.11e-16` (binary CCZ4). The object and single CCZ4 fingerprints are
`e1107d452abc4e4b` and `e68444ed5f33ba7b`; all nine malformed files are
rejected for their expected reasons. The echo reader accepts `ell=0.02149952219`
for B and `ell=0.05709817337` for E with the unchanged endpoint contract.

### B single-object grid at t=0

The existing C2 test gains an `echo` mode. At `N=128,256,512`, `dx=M/N`, its
cell-centred box is `x ∈ [-0.5,0.5]M`, cartoon `y ∈ [0,0.5]M`, with
`N × (N/2)` cells. It uses the same CCZ4 setter, Gamma, constraint, and Gauss
operators as the reference mode. The fixed physical mask is
`0.01 ≤ R=√(x²+y²) ≤ 0.5M`, outside B's `R_h=0.001418901919M`. The table gives
unweighted L2 residuals on that mask and `log2(E_coarse/E_fine)`; full L∞ values
and orders are in [`echo-convergence-results.csv`](echo-convergence-results.csv).

| dx/M | Mask cells | H L2 (order) | Mx/My L2 (order) | GaussE L2 (order) | GaussB L2 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1/128 | 6,444 | 1.253e-1 | 6.671e-1 / 6.671e-1 | 4.389e-2 | 0 |
| 1/256 | 25,722 | 2.821e-3 (5.47) | 3.935e-2 / 3.935e-2 (4.08) | 6.803e-4 (6.01) | 0 |
| 1/512 | 102,906 | 3.049e-4 (3.21) | 1.579e-3 / 1.579e-3 (4.64) | 6.782e-5 (3.33) | 0 |

```sh
for n in 128 256 512; do
  OMP_NUM_THREADS=1 Tests/EMSTrumpet/EMSTrumpetGridConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
    "$n" 0 Tests/EMSTrumpet/fixtures-echo/B/reference.trumpet echo > "b-$n.log"
done
python3 Tests/EMSTrumpet/echo_grid_results.py b-{128,256,512}.log > Tests/EMSTrumpet/echo-convergence-results.csv
```

The three boxes have 8,192, 32,768, and 131,072 interior cells. C++ reported
0.076, 0.285, and 1.109 seconds respectively with one OpenMP thread. All mask
points stayed above the `1e-12` chi and lapse floors, and the largest normalized
field comparison error was `2.78e-16`. All four nonzero L2 constraints decrease,
but their fine-pair orders range from 3.21 to 4.64; this is not a uniform
fourth-order regime. At `M/512`, the inner mask radius is 5.12 grid steps from
the center and the `r≈M` feature at `R≈0.16M` is about 82 steps out. Thus the
finest box samples the stated mask, but it does not resolve the horizon:
`2R_h/dx=1.45` cells across its diameter.

### Head-on refinement estimate

B's areal horizon radius is `r_h/M=0.1272544`, while its isotropic coordinate
radius is `R_h/M=0.001418901919`: the areal-to-isotropic ratio is 89.69.
The exp-0001 layout has base `dx_0=2M` and doubles resolution on each level,
so the diameter cell count on level `L` is
`2R_h/dx_L = (R_h/M) 2^L`. The calculation is reproducible with
`python3 -c 'import math; h=0.0014189019194714417; print([(l,h*2**l) for l in range(16)]); print(math.ceil(math.log2(25/h)))'`.

| Level | dx/M | Cells across horizon diameter |
| ---: | ---: | ---: |
| 0 | 2 | 0.00142 |
| 1 | 1 | 0.00284 |
| 2 | 1/2 | 0.00568 |
| 3 | 1/4 | 0.01135 |
| 4 | 1/8 | 0.02270 |
| 5 | 1/16 | 0.04540 |
| 6 (exp-0001 finest) | 1/32 | 0.09081 |
| 7 | 1/64 | 0.18162 |
| 8 | 1/128 | 0.36324 |
| 9 | 1/256 | 0.72648 |
| 10 | 1/512 | 1.45296 |
| 11 | 1/1024 | 2.90591 |
| 12 | 1/2048 | 5.81182 |
| 13 | 1/4096 | 11.62364 |
| 14 | 1/8192 | 23.24729 |
| 15 | 1/16384 | 46.49458 |

With diameter as the meaning of "across," level 15 is the first level with
at least 25 cells: nine more levels than exp-0001. Its spacing and CFL timestep
are `1/512` of level 6 (`dx_15/M=6.1035e-5`). This ratio is a minimum local
step-count cost for a fixed physical duration; the total work also depends on
the refined box volume and subcycling. These grid results are t=0 diagnostics,
not a head-on evolution.

Fixtures come from the supplied EMS.jl worktree
`/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-binsup-QSsl2l`.
Its density-seed source is uncommitted on `f477b22`; the fixture generator is
`scripts/emstrumpet_fixtures.jl`. The copied `manifest.txt` records the fixture
hashes. `binary_ccz4_fieldsum.tsv` retains the previous binary as a control.

Build and run from the worktree root with the local serial Chombo toolchain:

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
make -C Tests/EMSTrumpet all DIM=2 -j 12
Tests/EMSTrumpet/EMSTrumpetFixtureTest2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
  Tests/EMSTrumpet/fixtures
```

`EMSTrumpetSolution_read` owns file validation and radial reconstruction;
`EMSBH_trumpet_read` uses GRChombo tensors and cartoon variables for initial
data. `check_file(path, reason)` checks malformed files without aborting. The
production `main(path)` uses the same validation and calls `MayDay::Error` on
failure.

The runner compares every numeric column in `jets.tsv`, `radial.tsv`,
`objects.tsv`, `single_ccz4.tsv`, and `binary_ccz4.tsv`, checks each malformed
file's rejection reason.
It prints the largest normalized error per table, using
`abs(actual - expected) / (1 + abs(expected))`; this gives an absolute check
near zero. The density binary limit is `5e-13`; other limits are `1e-12` for
field values and coordinates, `1e-10` for
first derivative quantities (`jets` d1, radial PR/alphaR/betaR_deriv,
scalar Pi, and curvature K/A), and `1e-8` for higher jets d2/d3.

## Density-seed binary at t = 0

The setter sums the closed electric densities and the magnetic densities
obtained from each unchanged single-object tensor. It lowers both with the
physical metric `gamma_left + gamma_right - I`, then divides by its volume
factor and, for E, the coupling at the superposed scalar. The ADM assembly
returns that physical metric. One conversion helper supplies the final CCZ4
fields to both the grid setter and the fixture API, with initial lapse
`sqrt(chi)` and zero constraint/gauge auxiliaries. The level then fills ghosts
and calculates Gamma; the shared second ghost fill and gauge setup follow.

At separation 32 and rapidities ±0.20273255, the binary fixture's largest
normalized error is `1.107e-16`, below `5e-13`. The 90 single-object CCZ4
rows have largest error `1.149e-15`; their table hash remains
`18f1187e6cafc1f9b4504603210dbb649f96503ad705efac56acdb665ebada44`.
The fixture runner also checks `e1107d452abc4e4b` for full 3D object rows,
`e68444ed5f33ba7b` for single CCZ4 rows, and `024d0092dba552bc` for density
binary CCZ4 rows. The accepted fingerprint changes are documented above.
All nine malformed files are rejected. The Julia density binary test passes
110/110 checks.

Run the binary grid test at `N=32,64,128` (spacing `1/N`) with:

```sh
for mode in density old; do
  for n in 32 64 128; do
    OMP_NUM_THREADS=1 Tests/EMSTrumpet/EMSTrumpetBinaryGridConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
      "$n" 0.20273255 Tests/EMSTrumpet/fixtures/reference.trumpet "$mode" > "$mode-$n.log"
  done
done
cd /Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-binsup-QSsl2l
julia --project=test scripts/binary_density_fd_ad.jl /path/to/{density,old}-{32,64,128}.log
```

The physical masks are each hole's `0.5 <= R_rest <= 1.5` collar and the
inter-hole rectangle `|x| <= 14, 0 < y <= 1.5`. Their axis strips have fixed
physical width `0 < y <= 0.125`. All masks exclude the punctures. The table
shows GaussE L2 at the three spacings and orders on the coarse and fine pairs;
left and right agree to printed precision. The old field sum is a control.

| Mask | Density GaussE L2 at N=32/64/128 | Density orders | Old GaussE L2 at N=128 |
| --- | --- | --- | ---: |
| Left/right collar | `3.690e-8 / 2.296e-9 / 1.434e-10` | `4.007 / 4.001` | `5.936e-4` |
| Inter-hole | `5.831e-11 / 3.646e-12 / 2.279e-13` | `4.000 / 4.000` | `1.762e-5` |
| Left/right axis strip | `3.938e-8 / 2.649e-9 / 1.675e-10` | `3.894 / 3.983` | `6.320e-4` |
| Inter-hole axis strip | `9.702e-11 / 6.071e-12 / 3.795e-13` | `3.998 / 4.000` | `1.952e-5` |

The FD−AD H, Mx and My L2 sample differences converge at fourth order. For
left collar, the fine-pair orders are `3.949, 3.962, 3.941`; for right collar
`3.911, 4.087, 3.988`; for inter-hole `3.510, 3.961, 4.001`. The inter-hole
H difference at N=128 is `8.33e-11`, approaching numerical cancellation.
The old control has the same FD−AD differences to printed precision because
only its algebraic electromagnetic source changes; its GaussE remains finite.
Raw H and M do not tend to zero for this density-only seed. All values,
including old/control norms and both resolution-pair orders, are in
[`binary-convergence-results.csv`](binary-convergence-results.csv); the exact
FD/AD sample coordinates and residuals are in [`binary-fd-ad.csv`](binary-fd-ad.csv).

The 10-step legacy smoke `constraint_norms.dat` is byte identical to the
archived file after adding the two existing Gauss diagnostic parity entries
(`0 0`) to a run copy of the old parameter file. The single-object grid
output is byte identical with the pre-change binary setter. Against the
committed `convergence-results.csv`, however, the worst relative norm
difference is `8.405e-10` (N=128, eta=-0.20273255, full My_Linf), exceeding
the requested `1e-12` gate; the pre-change setter has the same difference.

## GRChombo grid check (C2)

From the worktree root, with the local Chombo build:

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
make -C Tests/EMSTrumpet all DIM=2 -j 12
for eta in 0 0.20273255 -0.20273255; do
  for n in 32 64 128; do
    Tests/EMSTrumpet/EMSTrumpetGridConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
      "$n" "$eta" /Users/auroradysis/Workspace/EMS/artifacts/reference-alpha20-qfile.trumpet
  done
done
```

This single-box test has no AMR or file output. It initializes the mass-one
object on cell centres in `x ∈ [-2,2]`, cartoon `y ∈ [0,2]` with spacing
`1/N`, then runs the same CCZ4 setter, Gamma calculator, Hamiltonian and momentum
constraint operators used by EMS. It also checks the signed electric and
magnetic Gauss diagnostics and compares the final CCZ4 fields directly with
the accepted kernel. The fixed mask is `0.25 ≤ R_rest ≤ 1.5`, with
`R_rest = hypot(cosh(eta)*x,y)`; the separate axis subset is the first two
positive-y rows within that mask. The test fills analytic ghost cells and
calculates Gamma in two ghost rows so fourth-order stencils work at the axis.

`convergence-results.csv` records each run's L2 and L∞ norms and observed
orders `log2(E_(N/2) / E_N)` for H, Mx, My, GaussE, GaussB. The columns
`chi_floor` and `lapse_floor` count mask points below the configured `1e-12`
floors before any clipping; all are zero. `field_error` is the maximum
`abs(grid - kernel)/(1 + abs(kernel))` across the 21 sampled CCZ4 fields.
Zero GaussB is reported without an observed order. The suggested fourth-order
gate is 3.5–4.5 on the finer pair for nonzero errors; the few exceptions in
the axis subset and the unboosted full-mask GaussE L∞ are retained as measured.

The EMS example can parse and initialize the reference file using
`Examples/EMS/params-trumpet-single.txt`; the corresponding binary parameter
file is `params-trumpet-binary.txt`. Both examples use `ems_data_format =
emstrumpet1` and the reference coupling. Their `max_steps = 0` requests an
initialization only. The binary data are also covered by the fixture and
evolution checks below.

```sh
make -C Examples/EMS all DIM=2 -j 12
cd Examples/EMS
mkdir -p chk plt
OMP_NUM_THREADS=1 ./Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
  params-trumpet-single.txt
```

## Short evolution (C3)

The C2 `axis` rows cover the first two y rows and therefore move toward the
axis as N increases. C3 uses the fixed physical strip `0 < y ≤ 1/8` inside `0.25 ≤ R_rest ≤ 1.5`.
Its nine-run norms and orders are in
[`evolution/fixed-axis-strip.csv`](evolution/fixed-axis-strip.csv).

Build EMS as above. The run files in `Examples/EMS/params-trumpet-evolve-*.txt`
use the delivered reference profile, `OMP_NUM_THREADS=1`, CFL 0.25 and
`min_chi = min_lapse = 1e-12`. Run each from a distinct working directory
containing `chk/` and `plt/`, passing the absolute path to its parameter file.
The unboosted `n16`, `n32`, `n64` files cover `[-4,4] × [0,4]` to `T=1/4`;
`boosted` uses the same box at N=32; `binary` covers `[-20,20] × [0,4]` at
N=16 with separation 32. `binary-n32` and `binary-n64` add the two resolutions
needed for binary field self-convergence. `amr` has three levels and regrids
after each coarse step. The binary files use the existing opposite-rapidity
placement; none changes the evolution equations or gauge.

From the worktree root, after building EMS:

```sh
worktree=$PWD
run_root=/private/tmp/emstrumpet-evolution
for case in n16 n32 n64 boosted binary binary-n32 binary-n64 amr rh; do
  mkdir -p "$run_root/$case/chk" "$run_root/$case/plt"
  (cd "$run_root/$case" && OMP_NUM_THREADS=1 \
    "$worktree/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex" \
    "$worktree/Examples/EMS/params-trumpet-evolve-$case.txt" > run.log 2>&1)
done
```

Analyze the resulting directories with:

```sh
uv run --with h5py --with numpy python Tests/EMSTrumpet/evolution/analyze.py RUN_ROOT GRID_LOG
```

`RUN_ROOT` contains directories named `n16`, `n32`, `n64`, `boosted`, `binary`,
`binary-n32`, `binary-n64`, and `amr`, each with its `plt/` directory.
`GRID_LOG` is the captured stdout of the nine-run grid loop above; omit it
when reanalyzing evolution outputs only.
The script reconstructs Chombo boxes, checks every plotted component for
finite values, and writes the CSVs beside itself. It uses the fixed exterior
mask `0.75 ≤ R ≤ 1.5`; binary masks are the union of such shells around the
two objects. Field Richardson ratios use the N=16 cell centres common to the
analysis, with degree-five tensor-product interpolation of N=32 and N=64
values. `common-constraints.csv` uses those same physical points. The
`single-constraints.csv` table retains each grid's own mask sampling as a
separate check. Floor counts in the CSVs count plotted interior cells below
`1e-12`, not individual RK-substage clamp calls.
The measured times, convergence summary, smoke outcomes and limitations are
in [`evolution/results.md`](evolution/results.md).

## Boundary parity

`GaussE` and `GaussB` are scalar diagnostic variables with parity `0`. Parameter
files written for the previous diagnostic list must append two `0` entries
to `vars_parity_diagnostic`; the EMS example files do so. The core
`BoundaryConditions.cpp` uses the standard GRChombo array loading path.
