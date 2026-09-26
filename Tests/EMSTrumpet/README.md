# EMSTRUMPET 1 initial data checks

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
factor and, for E, the coupling at the superposed scalar. The metric returned
to `FixSuperposition_metric` remains `gamma_left + gamma_right`.

At separation 32 and rapidities ±0.20273255, the binary fixture's largest
normalized error is `1.119e-16`, below `5e-13`. The 90 single-object CCZ4
rows have largest error `1.149e-15`; their table hash remains
`18f1187e6cafc1f9b4504603210dbb649f96503ad705efac56acdb665ebada44`.
The fixture runner also checks bit fingerprints frozen from the pre-change C++
setter: `e1107d452abc4e4b` for full 3D object rows and `c30cf22126f14d61`
for single CCZ4 rows.
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
`1/N`, then runs the same finalizer, Gamma calculator, Hamiltonian and momentum
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
