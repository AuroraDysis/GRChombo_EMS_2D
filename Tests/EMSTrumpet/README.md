# EMSTRUMPET 1 initial data checks

Fixtures come from `/Users/auroradysis/Workspace/EMS/test/fixtures/emstrumpet1/`.
The EMS.jl reference revision for this port is `14fa579`; its `src/trumpet.jl`
matches the supplied read-only file byte-for-byte. The fixture generator is
`scripts/emstrumpet_fixtures.jl`. The copied `manifest.txt` records the fixture
hashes and its generating worktree revision (`5e894187`).

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
near zero. Limits: `1e-12` for field values and coordinates, `1e-10` for
first derivative quantities (`jets` d1, radial PR/alphaR/betaR_deriv,
scalar Pi, and curvature K/A), and `1e-8` for higher jets d2/d3.

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
