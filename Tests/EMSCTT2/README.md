# EMSCTT/2 single-hole pulse initialization

**The opt-in reader/setter reproduces the Julia reconstruction.** Zero-pulse,
field, consumed-derivative and old-reader controls pass. Two-level t=0 AMR
hierarchies initialize successfully. T2's production point transfer removes
the large outer error seen with T1's default legacy transfer. This is an
initialization/readback result, not an echo evolution or a production-resolution
certificate.

## Selection and defining construction

Use the existing parameter; its default remains empty:

```text
ems_not_rn = true
ems_data_format = emstrumpet1
ems_data_path = /Users/auroradysis/Workspace/EMS/artifacts/echo/a0.9-e7.trumpet
ems_ctt_data_path = /Users/auroradysis/Workspace/EMS/.data/echo-pulse/n24-final.ctt
binary = false
boosted = false
boost_rapidity = 0
separation = 0
bh_mass = 1
ems_alpha = 12.566370614359172
ems_f0 = 0
ems_f1 = 0
ems_f2 = -0.9
ems_use_maximal_initial_lapse = false
ems_use_geometric_initial_lapse = false
amr_transfer = point
num_ghosts = 3
max_steps = 0
stop_time = 0
```

The version line selects the separate `EMSCTTPulse_read` parser. EMSCTT/1
still selects its original binary reader. No-companion and legacy defaults
are unchanged. On every rank, `SimulationParameters::check_params` reads the
version before hierarchy setup and refuses EMSCTT/2 unless `amr_transfer = point`
and `num_ghosts = 3`. This is a mandatory startup guard.

The single object is positioned at `star_centre`, on the cartoon axis. The reader requires `bh_mass = 1` because the companion's coordinates,
panels and pulse parameters are already expressed in units of the hole's mass M.
In geometric units G=c=1, this massless EMS theory with dimensionless couplings
has no independent dimensional scale besides M: code coordinates are x/M and
t/M, while code K, Pi and Cartesian Maxwell fields represent M times their
physical values. Mass one selects those units; it is not a hard-coded physical
mass or a restriction to a particular mass in external units. The reader does
not perform a second runtime rescaling of the supplied representation.
Zero boost/separation and the default initial lapse are also checked. Physical
pulse amplitude, width, launch radius, panels and phase come from the file.

The format's defining expressions are retained:

```
rho = |x|, mu = x/rho, R = rho/X
phi_s = phi0 + A/R exp(-z^2/2) w(z) sqrt(5/(16pi)) (3mu^2-1)
z = xi(rho)/sigma
Pi_s = Pi0 - A Y20 X/R exp(-z^2/2) [w'(z)-z w(z)] xi'(rho)/sigma
              + k0 (phi_s-phi0)
W = (rho*a, y*b, z_cart*b)/(1+rho)
gamma = psi^4 gamma0
K = psi^-2 [K0 - gamma0 tau0/3 + gamma0 (L_gamma0 W) gamma0]
      + gamma tau0/3
phi = phi_s, Pi = psi^-6 Pi_s
E = psi^-2 F(phi0)/F(phi_s) E0, B = psi^-2 B0 = 0
```

Here `tau0` is the original trumpet's trace, zero to arithmetic accuracy.
`F(phi0) E0` fixes D0 before adding the pulse. The stored full-a0 optical
Chebyshev map and its analytic derivative define xi; no replacement integration
is performed. The specified taper has w=1 at |z|<=6, w=0 at |z|>=8, with the
specified exponential ratio in between; the declared radial support is also
respected exactly. The focusing momentum term above is the supplied Julia
expression, including its k0 term.

Radial Lobatto/angular Gauss nodal data are transformed once into tensor
Chebyshev coefficients by the corresponding cosine transforms. Analytic first
polynomial derivatives reconstruct W's Cartesian gradient without an axis
sine division. For the unboosted `gamma0=P^2 I`,
`gamma0 (L_gamma0 W) gamma0 = P^2 [dW+dW^T-2 div(W) I/3]`.
The exact executable witness is [ctt2_longitudinal.py](../../scripts/cas/ctt2_longitudinal.py);
its scope and independent runtime check are in [findings.md](findings.md).
The final physical metric reaches the existing shared ADM-to-CCZ4 packing,
with lapse sqrt(chi), original shift, and zero gauge/constraint auxiliaries.

## Authentication and parser controls

[authentication.json](results/authentication.json) records the supplied profile,
current specification/source hashes, five complete-file hashes, each payload
checksum, and all five provenance bindings per companion. The profile SHA-256 is
`74f85f39aa652b48d179a96b8dbe57cd4468d7905e39732aaf8bd1b75931b8fa`.
All five companion hashes match the supplied `data.sha256` manifest; source and
checkpoint provenance digests also match. The consumer checks the payload and
actual parsed profile bytes; other digests identify provenance, which it does
not load or execute.

The parser rejects unknown/missing/duplicate keys, unsupported reconstruction,
nonfinite data, invalid integers/dimensions/maps/topology, gaps/overlaps,
nonconverged solve metadata, incorrect phase/launch normalization, wrong profile
bytes, checksum damage, truncation, trailing data and a missing final newline.
[T3 parser controls](results/review/parser.json) record 31 rejection controls,
including checksum-valid CRLF and line-boundary truncation mutants.
[T3 binding controls](results/review/bindings.json) record eight constructor
rejections and four fatal domain controls executed inside a two-thread OpenMP
region. Evaluation failures use `MayDay::Error`; parser validation still returns
its structured failure outside BoxLoops. A CR in the version line is refused
with `EMSCTT 2 version line contains CR; use LF line endings`, in the parser,
dispatcher and startup check; checksum bytes are never silently normalized.
Below `1.2664165549094182e-20 M`, evaluation fails; there is no padding or
constant continuation. A positive explicit panel index is a fixture extension,
while production uses the first containing panel, including the inner shared face.

The SciData mount intermittently returned ENOENT/EMFILE during resumed readback.
Final comparisons use identical local input snapshots, checked against the
supplied hashes by `authenticate.py`. These regenerable snapshots are ignored
by git; no generating input or coefficient has been changed.

## Readback controls and tolerance basis

The fixed sample set has 146 off-grid points per companion, covering every
panel interior/face, both explicit face sides, exact/near axes, off-plane
locations, the inner continuation through 1e-18 M, the horizon vicinity, pulse
support/taper and compact exterior through 1e4 M. All 29 physical components
are compared for zero, half, n12, n18 and n24 companions.

Errors are `|C++-Julia|*s/(1+|Julia|*s)`, with fixed local weights X0^2 for
physical gamma/K and X0 for covariant E/B; scalars and shift have unit weights.
This removes the trumpet's coordinate blowup without changing the fields.
The retained trumpet/CTT transfer tolerance is **5e-13**, allowing Float64
profile inversion, matrix operations and spectral interpolation conditioning.
The worst field error is **1.291e-14** for each supplied companion; the zero
companion versus plain setter identity is **1.928e-15**.
[readback.json](results/readback.json) retains the measurements and timings.

An independent 512-bit Julia AD reference supplies 207 local Cartesian
D1/D2 samples, including both sides of every face and near-axis/deep-end points.
The C++ five-point fourth-order differences use h/rho=.002 and .001. D1 is
checked on every physical field; D2 on gamma, phi, Pi (tagging), and shift.
The allowance, fixed before comparison, is four times the Richardson
truncation estimate plus `128 epsilon64/h^j`, divided by `1+|AD jet|`.
The roundoff terms are 2.842e-11 / 2.842e-8 for D1/D2. Worst normalized errors
are **9.169e-11 / 1.560e-8**, at most **0.406 / 0.486** of the respective
pointwise allowances. This checks agreement with the delivered formulas;
it is not a physical waveform-error allocation.

Float64 nested AD at a |z|=6 taper face produced a spurious Pi derivative
about 26, while the actual derivative is about 6e-10. The original diagnostic
is retained in `jets-julia-float64.txt`; evaluating the unchanged formulas at
512 bits resolves the cancellation. BigFloat is needed here for its exponent
range in `exp(-1/t)` at the taper face. The consumer remains Float64 and its
values/finite differences pass; no Julia source was patched.

A second control directly compares C++ and Julia through 540 fourth-order
stencils (5,400 values), in x/y/z and three diagonal directions, exercising
mixed Hessians, automatic interface crossings and explicit panel extensions.
Value / D1 / D2 maximum errors are **1.217e-14 / 4.041e-11 / 3.122e-8**.
The retained 5e-13 pointwise tolerance is propagated by the absolute stencil
weights and h^-j; the derivative discrepancies use at most 1.044% / 0.703%
of those allowances. See [stencil-readback.json](results/stencil-readback.json).

The actual shared CCZ4 packing is also compared at 465 stencil values each
for the zero and n24 companions, including chi, h, A, K, matter, lapse and
shift. D2 explicitly covers chi/h/lapse/shift and the scalar tagging fields;
D1 covers all 21 returned fields. Maximum value / D1 / D2 errors are
**7.439e-15 / 1.163e-11 / 9.178e-9**, within the same propagated transfer
allowance. See [ccz4-readback.json](results/ccz4-readback.json).

## T1 initialized AMR constraints (legacy transfer)

T1 omitted `amr_transfer`; its parameter log explicitly records the default
`legacy`. Its historical parameter files allocated five ghosts. These measurements
are retained as the T1 baseline; T2 below supplies the production setting.

`EMSCTT2Hierarchy` builds a fresh real GRAMR hierarchy with two levels, ratio 2,
cell-centred cartoon parity, production ghost filling, Gamma initialization,
`Constraints<CouplingFunction>` and both `EMSCartoonGaussConstraints` diagnostics.
Its RHS throws. No time step, restart, horizon search, extraction or gauge
package is run. Coarse cells covered by fine boxes and all ghost cells are
excluded from norms. Each norm is area-weighted RMS in the Cartesian (x,y)
half-plane with measure dx dy, not physical 3-volume; maxima are also saved.
Momentum is `hypot(Mx,My)` in the diagnostic's covariant coordinate components,
not an orthonormal norm. H, momentum and both Gauss diagnostics have units M^-2.
The sampled mask cells never invoke the diagnostic chi/lapse floors. Grid
packing versus the direct helper stays within 2.12e-16, below 5e-13.

The global domain is [-80,80] x [0,80] M, with base N1=80/160/320 and fine
dx=1/.5/.25 M. Geometric tagging refines r<64 M; block layout follows Chombo.
The masks are r=4..8 (inner), 25..55 (pulse), 68..76 (outer), and the pulse
mask restricted to y<=2 (axis). Only resolution changes. Fine-pair orders:

| Region | H | Momentum | Electric Gauss | Magnetic Gauss |
|---|---:|---:|---:|---:|
| inner | 3.860 | 4.105 | 3.981 | exactly zero |
| pulse | 3.943 | 3.960 | 3.997 | exactly zero |
| pulse axis | 3.976 | 3.955 | 4.003 | exactly zero |
| outer | 0.723 | 1.791 | 1.774 | exactly zero |

At finest global resolution, pulse RMS H/M/GE are
**7.453e-13 / 1.543e-13 / 2.180e-13**; GB=0.
Outer H RMS is **1.701e-7** and does not show fourth-order convergence.
The paired analytic-ghost t=0 control reduces outer H/M/GE to
**1.384e-13 / 1.369e-17 / 1.127e-14**, identifying ghost transfer as the
source of that discrepancy. These are the legacy-transfer measurements;
the analytic control does not replace them. T2 below measures point transfer.

A separate small domain [-.16,.16] x [0,.16] M measures the horizon collar
r=.008...016 M, outside R_h=.007036891 M. Fine dx=.002/.001/.0005 M gives
H RMS **.12535, .0083387, .00054525** and momentum RMS
**3.9164, .23678, .015676**. Fine-pair H/M/GE orders are
**3.935 / 3.917 / 4.112**, with GB=0. These finite-spacing coordinate residuals
are substantial; no high-resolution evolution accuracy is inferred. Other
small-domain masks are the global masks scaled by .002; the inherited label
`pulse` in those logs denotes r=.05...11 M and contains no scalar packet.

For t=0 controls, the acceptance is finite diagnostics, valid-only nonempty
masks, no mask floors and correct field staging within 5e-13, followed by an
explicit report of all trends. Fourth order is the FD expectation away from
transfer/representation/roundoff errors, not a silently relaxed universal
pass gate. [constraints.csv](results/constraints.csv) and
[constraint-tables.md](results/constraint-tables.md) preserve every RMS,
maximum, cell count and both adjacent-pair orders.

## T2: production point transfer at t = 0

**Yes on the measured three-resolution trend in every mask.** Use
`amr_transfer = point` and `num_ghosts = 3` for this EMSCTT/2 initialization.
No default or transfer implementation has been changed. All nonzero point-mode
RMS and maxima decrease on both pairs: 96 positive pair orders across the
global/near domains, four masks and H/momentum/electric Gauss. Magnetic Gauss
is zero everywhere, so its order is undefined. The smallest pair order is
2.055, the global outer H RMS fine pair, already near the analytic reference
residual. Uniform fourth order at that floor is not claimed.

The original domains, two levels, refinement ratio, tagging,
box sizes, boundary conditions and masks are unchanged. The point operator
requires exactly three ghosts (`PointAMRTransfer::define`); it cannot accept
T1's five-ghost allocation. The harness therefore computes Gamma on valid
cells and fills Gamma ghosts from current arrays in the three-ghost cases,
without reading metric stencils outside the allocation. The harness retains
its five-ghost calculation for the historical reference; EMSCTT/2 startup now
refuses five ghosts. A matched `legacy3` control separates this harness adjustment
from the transfer setting: its finest global outer H RMS is 1.701e-7, the same
legacy defect as five ghosts. The legacy replay reproduces all archived T1
norms with measured zero difference, without making byte identity an
acceptance criterion. [layout-check.json](results/point-transfer/layout-check.json)
verifies identical actual boxes, uncovered valid-cell totals and mask counts
for every transfer/ghost case at all six domain/resolution combinations.

The new results keep five controls beside one another: `legacy` (five ghosts,
T1 prescription), `legacy3` (three ghosts, matched Gamma handling), `point`
(production point transfer, three ghosts), `analytic` (analytic metric/matter
ghosts, three ghosts with Gamma ghosts filled through the point operator), and
`analytic5` (the original fully analytic five-ghost reference, rerun). Every
case is a fresh hierarchy at t=0; no evolution or restriction-after-stepping
is invoked. Norm definitions and the 5e-13 staging allowance are unchanged.

Global outer RMS values and adjacent-pair orders:

| Mode | N1 | H | order | Momentum | order | Electric Gauss | order |
|---|---:|---:|---:|---:|---:|---:|---:|
| legacy | 80 | 3.41805e-7 | — | 1.55550e-13 | — | 3.51737e-10 | — |
| legacy | 160 | 2.80679e-7 | 0.284 | 6.86002e-14 | 1.181 | 1.46591e-10 | 1.263 |
| legacy | 320 | 1.70085e-7 | 0.723 | 1.98241e-14 | 1.791 | 4.28686e-11 | 1.774 |
| point | 80 | 9.05084e-12 | — | 1.60377e-15 | — | 1.44569e-12 | — |
| point | 160 | 5.76596e-13 | 3.972 | 1.98956e-16 | 3.011 | 1.65850e-13 | 3.124 |
| point | 320 | 1.38781e-13 | 2.055 | 1.36879e-17 | 3.861 | 1.12731e-14 | 3.879 |

Fine-pair point-mode RMS orders in all masks:

| Domain | Mask | H | Momentum | Electric Gauss |
|---|---|---:|---:|---:|
| global | inner | 3.860 | 4.105 | 3.981 |
| global | pulse | 3.943 | 3.960 | 3.997 |
| global | pulse axis | 3.976 | 3.955 | 4.003 |
| global | outer | 2.055 | 3.861 | 3.879 |
| near | inner | 3.935 | 3.917 | 4.112 |
| near | pulse | 3.995 | 3.999 | 3.997 |
| near | pulse axis | 3.997 | 4.022 | 3.994 |
| near | outer | 3.818 | 3.872 | 3.899 |

The finest global outer analytic reference is H=1.38374e-13,
M=1.36878e-17, GE=1.12731e-14 (three ghosts); the independently rerun original
five-ghost reference gives H=1.38373e-13 with the same displayed M/GE.
Point H is within 0.30% of these reference residuals. This is compatible
with reaching the analytic initialization/FD arithmetic residual and explains
why a fourth-order claim on the finest global H pair would overstate the data;
it does not isolate a unique source of the remaining arithmetic error.
The near-domain outer point H RMS sequence is
5.96537e-7, 6.69029e-8, 4.74331e-9, with orders 3.156/3.818;
its analytic reference is 5.75839e-7, 6.65427e-8, 4.74115e-9
(see the full-precision tables for the controls). No nonconverging outer
residual remains on these ladders, so the conditional field/boundary
localization is not required. The transfer code is untouched.

[All RMS, maxima and both pair orders](results/point-transfer/tables.md) are
reported per mask and diagnostic, including both analytic controls and the
matched legacy control. [norms.csv](results/point-transfer/norms.csv) preserves
960 full-precision rows. [authentication.json](results/point-transfer/authentication.json)
binds the consumed profile/companion snapshots to T1's supplied hashes and
records the executable, source and parameter hashes. Runtime profile/payload
checks remain active. [cost.json](results/point-transfer/cost.json) records each
command, exit, wall time and peak child RSS. The 30 hierarchies took
54.70 s total; individual costs are 0.34–5.45 s, with a conservative child
RSS peak of 122.4 MB.
All runs are serial and below five minutes individually, within the
two-process/6 GB-per-process budget.

T2 is frozen historical evidence. Its legacy runs are deliberately refused by
the T3 startup guard. The original T2 driver is archived as
`T1-T2/source-transfer_controls.py` with its hash; the current driver below
runs only supported point-transfer configurations.

## T3: review fixes and panel-interface controls

[Startup controls](results/review/startup.json) verify rejection of explicit
legacy transfer (three and five ghosts), omitted transfer (the unchanged legacy
default), point with five ghosts, and the CRLF version line. Two ghosts are
rejected earlier by the existing base derivative-stencil validation. Point with
three ghosts and no-companion legacy startup are accepted. All these checks end
before initialization; the twelve hierarchy controls below also verify actual
initialization of the accepted EMSCTT/2 configuration. The guard has no rank-zero
condition and uses the code's fatal error path on every rank. Local builds use
serial Chombo; multi-rank execution was not tested.

The physical ADM return type is retained for existing callers. In the pulse
branch, `chi` is recomputed from the final metric determinant and `lapse` is
sqrt(chi), matching shared CCZ4 packing. Boost factors J and w are explicitly
one because this branch admits only an unboosted object. The plain/binary paths
are unchanged. Consistency checks cover all 146 sample points for both zero and
n24 companions, with the retained 5e-13 normalized allowance for chi/lapse.
The field, AD-derivative, stencil and CCZ4 comparisons against the authenticated
archived Julia references pass again; [readback](results/review/readback.json),
[stencils](results/review/stencil-readback.json), and
[CCZ4](results/review/ccz4-readback.json) retain the numbers.

The global hierarchy adds four shells of half-width 2 M, centred on radii
read from the supplied n24 panel rows, rather than changing the representation:
16.72433994232569, 22.149975026273744, 56.262273380783796 and
62.06455126251464 M. They cover the omitted faces in the 8–25 M and 55–68 M
gaps. The original four masks, domains, tagging, levels and layouts are retained.
Point and analytic-ghost controls run at N1=80/160/320; the near-domain ladder
also repeats its original four masks. Valid cells only, norm definitions and
5e-13 staging allowance remain as above. Layouts match T2; original-mask norms
reproduce within 5e-13 (measured normalized difference zero).

New-shell point-transfer RMS values and pair orders (GB is zero throughout):

| Face / M | Diagnostic | N1=80 | N1=160 | N1=320 | order 80→160 | order 160→320 |
|---:|---|---:|---:|---:|---:|---:|
| 16.72434 | H | 1.62855e-08 | 1.14967e-09 | 7.00461e-11 | 3.824 | 4.037 |
| 16.72434 | M | 7.56803e-12 | 5.25861e-13 | 3.14756e-14 | 3.847 | 4.062 |
| 16.72434 | GE | 4.39380e-09 | 3.07430e-10 | 1.83412e-11 | 3.837 | 4.067 |
| 22.14998 | H | 2.23367e-09 | 1.42997e-10 | 8.95358e-12 | 3.965 | 3.997 |
| 22.14998 | M | 9.00346e-13 | 5.67821e-14 | 3.51708e-15 | 3.987 | 4.013 |
| 22.14998 | GE | 6.67083e-10 | 4.21284e-11 | 2.60501e-12 | 3.985 | 4.015 |
| 56.26227 | H | 2.93896e-12 | 1.90239e-13 | 1.84791e-13 | 3.949 | 0.042 |
| 56.26227 | M | 3.57761e-15 | 1.99754e-16 | 1.16671e-17 | 4.163 | 4.098 |
| 56.26227 | GE | 1.04988e-12 | 6.56904e-14 | 4.08736e-15 | 3.998 | 4.006 |
| 62.06455 | H | 1.46753e-12 | 1.03278e-13 | 1.90281e-13 | 3.829 | -0.882 |
| 62.06455 | M | 9.18204e-16 | 5.75537e-17 | 3.59953e-18 | 3.996 | 3.999 |
| 62.06455 | GE | 5.28466e-13 | 3.31319e-14 | 2.07044e-15 | 3.996 | 4.000 |

**The expanded mask set does not establish decreasing H in every mask.**
At 56.26 M the finest H RMS is 1.848e-13 (order .042); its maximum rises to
6.171e-13 (order −.493). At 62.06 M H RMS rises from 1.033e-13 to 1.903e-13
(order −.882), and its maximum rises to 6.476e-13 (order −.996).
The analytic-ghost control gives the same measured values and orders in all
four new shells, within the 5e-13 replay allowance (measured difference zero).
Thus these new H limitations survive replacement of transferred metric/matter
ghosts by analytic ones; they do not reproduce the legacy coarse-fine defect.
They are compatible with the analytic initialization/FD arithmetic residual,
without isolating a unique source or asserting asymptotic convergence there.
Momentum and electric Gauss decrease at both pairs in every new shell; the
first two shells show approximately fourth-order H/M/GE convergence as well.
The original outer-mask point-transfer H trend remains 3.972/2.055.

[Every RMS, maximum, cell count and pair order](results/review/tables.md),
including analytic ghosts, is preserved with a
[full-precision CSV](results/review/norms.csv) and
[decision receipt](results/review/decision.json). No transfer code changed.
Twelve fresh t=0 hierarchies cost 18.20 s total, at most 3.76 s each and
112.5 MB child RSS. [Cost receipts](results/review/cost-transfer.json) and
[authentication](results/review/authentication.json) bind inputs and executables.
The EMS production build and both old fixture suites pass again; the build
costs 10.82/13.78/29.38 s for EMSCTT/trumpet/EMS, and fixture runs .97/.20 s.
[Regression receipts](results/review/cost-regression.json) give exact costs and
child RSS (conservative maximum 751.5 MB). No job needed the over-five-minute tier.

Raw evidence, including numeric reference/sample tables, now lives under
`/Users/auroradysis/Workspace/EMS/.data/fork-ctt2-logs/`, in `T1-T2/` and `T3/`.
Only evidence tables and integrity/cost/decision receipts remain in git.
[The archive receipt](results/raw-archive.json) records every archived T1/T2
file hash and the manifest hash; [T3 raw hashes](results/review/raw-manifest.json)
cover the new logs. `evidence.py` resolves reference files from this authenticated
archive, preferring newly regenerated references in T3. Archived T1/T2 receipts
keep their original command paths as provenance; raw files moved, inputs did not.

Reproduce current supported controls:

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
make -C Tests/EMSCTT2 all DIM=2 -j2
python3 Tests/EMSCTT2/authenticate.py
python3 Tests/EMSCTT2/parser.py
python3 Tests/EMSCTT2/transfer_controls.py
uv run --no-project --with numpy python Tests/EMSCTT2/compare.py
uv run --no-project --with numpy python Tests/EMSCTT2/compare_stencils.py
uv run --no-project --with numpy python Tests/EMSCTT2/compare_ccz4.py
python3 Tests/EMSCTT2/run_controls.py regression
# Regenerate the current norm tables from completed logs:
python3 Tests/EMSCTT2/transfer_controls.py report
```

## T1 reproduction and cost

The current supported commands are above. To regenerate Julia reference values,
pass `/Users/auroradysis/Workspace/EMS/.data/fork-ctt2-logs/T3` as the output
directory to `reference.jl`; its profile/data arguments remain the authenticated
input snapshots. `run_controls.py reference-ccz4` and `reference-ad` also write
raw outputs there. Historical `constraints_report.py` reads archived T1 logs;
its legacy tables are not regenerated from the now-production parameter files.

`compare.py` consumes the archived AD reference. Regenerate it with
`python3 Tests/EMSCTT2/run_controls.py reference-ad`; the replay output is
`jets-replay.txt`, deliberately separate from the archived reference.
`reference_ad.jl` and `reference_ccz4.jl` retain the executable reference paths.
Long commands were backgrounded with watchdogs and backed-off log checks.
Every job stayed within the under-five-minute tier; none exceeded 30 minutes.
The interrupted model-service turn left completed receipts that were reused;
interrupted builds/regressions were finished after resumption.

The production EMS build and both unchanged fixture suites pass. EMSCTT/1
retains 42 physical + 42 CCZ4 points, SHA controls and parser tests; the trumpet
suite retains single/binary/radial/echo controls. Their numeric tolerances are
unchanged. No new bit-identity acceptance requirement was introduced.
[T1 regression log](/Users/auroradysis/Workspace/EMS/.data/fork-ctt2-logs/T1-T2/regression.log) and `cost-*.json` record commands,
wall times, exit status and child peak RSS. Final production build: 21.23 s;
old CTT/trumpet readbacks: .80/.18 s. Global and small-domain AMR ladders cost
about 1–4.2 s per hierarchy, with peak RSS <=118 MB. The field readback costs
.04–.08 s per companion; the 5,400-point C++ stencil comparison costs .59 s.

The independent 512-bit AD replay completed in **142.58 s**, peak RSS
**883 MB**; its normalized difference from the archived reference is zero
under a 1e-12 conditioning allowance, not a byte-identity requirement.
The CCZ4 Julia reference cost **9.66 s**, peak RSS **543 MB**. At most two
heavy processes ran together; recorded process peaks remain below 6 GB.
