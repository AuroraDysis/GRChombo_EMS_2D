# EMSCTT 1 C++ handoff

**Decision: yes on the requested fixture and t=0 checks, with the declared numerical-floor qualification.**
All 42 physical and 42 CCZ4 points pass. The original T7 handoff preserved all
13 no-companion bit fingerprints; the subsequent direct CCZ4 conversion is
documented below. On the exterior collars the fine-pair H/M/GaussE orders are
3.996/3.993/3.996; the density control tends to H≈3.626e-3 and M≈1.575e-3.
Inter-hole H shows fourth order over three levels before Float64 evaluation
noise dominates the finest spacing. All raw orders and the noise diagnostic
remain in the tables.

This is the file-only C++ consumer of the accepted degree-28, separation-32 M
CTT binary. It evaluates the exported representation; it does not solve (10).
No solver, refit, new physics or external hash library is used. The Julia worktree
only receives the controller's certification wording in the format document's
first paragraph and the fixture manifest's `certification` string.

Reader design (eight lines):

1. `EMSCTTSolution_read` checks all 35 records, exact convention/schema, finite numbers and decimal integer controls.
2. It rejects unknown/duplicate/missing keys, malformed blocks/maps, end gaps/overlaps, invalid `END` and trailing bytes.
3. The self-contained SHA-256 checks the payload trailer; the profile binding hashes the exact bytes actually parsed.
4. Both masses, relative centres, rapidities and all four couplings must equal the run and profile parameters.
5. Log/end, bridge (plane/sphere) and inverse maps follow the Julia reader, including panel order and 8-epsilon admission.
6. Hole-relative coordinates precede maps/Clenshaw; five radial-fast tensor series produce log ψ and Cartesian C using `(y,z)/R`.
7. The density seed receives the documented ψ/C reconstruction; its physical γ_final goes directly through the setter's shared ADM-to-CCZ4 conversion, with initial lapse √χ.
8. Below the minimum represented radius, evaluation fails: no clipping, extrapolated end or constant padding. The trumpet reader likewise rejects its limiting cylinder.

## Inputs and activation

The single new parameter is `ems_ctt_data_path`, default empty. It requires
`ems_data_format = emstrumpet1` and `binary = true`. The existing setter places
holes at `±separation/2` relative to `star_centre`; the companion records those
relative centres. Its effective rapidity is zero if `boosted = false`, so that
setting is rejected for the supplied boosted companion. Explicit binary helper
calls also check mass/separation/rapidity against the loaded binding.

The complete read-only companion (11,387,862 bytes) remains in the EMS.jl worktree:

```text
.data/binary-ctt-t6/n28-r6.ctt
ecba66b51c0bdb9a470afa50194dcf23bd0d7163a1af1e4b4b1fc5408b7e11b1
```

The EMSTRUMPET profile is `artifacts/reference-alpha20-qfile.trumpet` (5,464 bytes):

```text
6f0820a576620f1f7230c701131af56a2312b130e31d5de2a752c32cefe6d24f
```

These complete-file hashes and both fixture hashes match the supplied manifest;
[authentication.json](results/authentication.json) records them. The payload
SHA-256 in the companion's trailer is separately verified. No companion or Julia
fixture is copied here. The represented minimum radius is
`1.2664165549094182e-20 M`; `evaluate(dx,y,z,panel,origin)` also supports offsets
smaller than an ulp of a hole's global centre. Forced 1-based panels are for
one-sided fixture checks; production uses automatic selection.

## Direct CCZ4 conversion validation

The physical metric now reaches the shared conversion directly, removing
the CTT diagonal `(gamma + 1) - 1` round trip. The converter uses production's
`pow(det, -1./3.)`; the old fixture-only path used `1/cbrt(det)`. Against
b397e42, the CTT CCZ4 numeric fingerprint changes from `4caed5832f720e5a` to
`b3b32ee0d0a729dc` (42 rows, 28 returned fields per row). This fingerprint
now has an explicit assertion; the old CTT test had only numerical parity
and input-file SHA assertions. All four authenticated input hashes and the
`5e-13` Julia parity tolerance are unchanged. The physical-table fingerprint
is printed for diagnosis, not used to replace any input hash.

The [trumpet README](../EMSTrumpet/README.md#direct-ccz4-initialization)
records the controller's corrected primitive/FD gate, why the initial STOP
was accepted, the unchanged production single/echo snapshots, and the
no-companion fixture fingerprint changes. For CTT and the T7 checkpoint,
the maximum primitive error is `4.4409e-16` absolute and `4.8805e-15`
relative to its field scale. Gamma/gauge-B maxima are `1.7090e-15` /
`1.2820e-15`, below the `7.5082e-14` FD allowance. Every evolved variable in
every cell, including ghosts, passes. The worktree build's t=0 checkpoint
also matches the instrumented candidate checkpoint byte for byte.

The final fixture run passes all 42 physical and 42 CCZ4 points, nine SHA
vectors, 30 parser/binding rejections and the domain/translated-coordinate
controls. The largest Julia scaled error is `5.9062e-15` in physical Kij
and `2.2700e-15` across CCZ4 fields. All eight CTT/density grid runs pass.
The exterior collar's fine H/M orders are `3.995778/3.993403` before and
`3.995777/3.993403` after; at dx=1/128, H/M L2 are `3.99101e-9/1.35497e-8`.
The inter-hole H floor changes from `5.67227e-11` to `5.54869e-11`
(finest order `1.801698` to `1.833311`); M remains fourth order.
Raw baseline/candidate norms and orders for every mask are retained in
`/private/tmp/ems-fixsup-92881/convergence-comparison.csv`, alongside the
resumed report, corrected gate and smoke evidence.

## Original T7 build and reproduction record

The commands, timings and result links below record the original T7 handoff.
Its byte-identical stdout comparison and `report.py` apply to the archived
pre-cleanup logs, not the updated CCZ4 fingerprints above. The C++ fixture
and grid executables remain the current regression tests; write new logs
to a separate directory when rerunning them.

From this fork worktree root, using the serial toolchain in `EMS-deps/BUILD.md`:

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
export OMP_NUM_THREADS=1
ems_fixture_root=/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/wt-binsup-t6
make -C Tests/EMSCTT all DIM=2 -j 2
make -C Tests/EMSTrumpet all DIM=2 -j 3
make -C Examples/EMS all DIM=2 -j 4

Tests/EMSCTT/EMSCTTFixtureTest2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
  "$ems_fixture_root/artifacts/reference-alpha20-qfile.trumpet" \
  "$ems_fixture_root/.data/binary-ctt-t6/n28-r6.ctt" \
  "$ems_fixture_root/test/fixtures/emsctt1" Tests/EMSCTT/results \
  > Tests/EMSCTT/results/parity.log

Tests/EMSTrumpet/EMSTrumpetFixtureTest2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
  Tests/EMSTrumpet/fixtures Tests/EMSTrumpet/fixtures-echo/B Tests/EMSTrumpet/fixtures-echo/E \
  > Tests/EMSCTT/results/current-fixtures.log
cmp Tests/EMSCTT/results/baseline-fixtures.log Tests/EMSCTT/results/current-fixtures.log

for mode in ctt density; do
  companion=density
  if test "$mode" = ctt; then companion="$ems_fixture_root/.data/binary-ctt-t6/n28-r6.ctt"; fi
  for level in 0.5 1 2 4; do
    /opt/homebrew/bin/gtimeout -k 5s 240s \
      Tests/EMSCTT/EMSCTTGridConvergence2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
      "$level" "$ems_fixture_root/artifacts/reference-alpha20-qfile.trumpet" "$companion" \
      > "Tests/EMSCTT/results/$mode-$level.log" 2>&1
    echo $? > "Tests/EMSCTT/results/$mode-$level.exit"
  done
done
python3 Tests/EMSCTT/report.py
```

Long commands were backgrounded and their log checks backed off from 15 seconds;
each build/run had a 240-second cap. The first real CTT path measured 3.62 seconds
at refinement 1, projecting about 58 seconds for refinement 4. All programs fit
the under-five-minute tier. No 5–30 minute or over-30 minute program was needed.
Initial test build: 6.07 s; full test rebuild: 6.84 s; trumpet tests: 6.55 s;
production EMS build: 16.05 s. The final fixture program took 0.312 s; the initial run took 0.431 s. Grid timings
are in [results.md](results.md), and logs retain shell wall/user/system times.
The EMS initialization-only smoke used the new parameter and exited 0 in 0.82 s;
its parameter file and log are in [results/smoke](results/smoke). No time evolution
or HPC run is claimed.

## Controls and bit identity

The fixture executable authenticates all four consumed files before comparing
42 physical rows and 42 CCZ4 rows, including explicitly selected interface sides.
Every field uses `abs(C++ − Julia)/(1 + abs(Julia)) ≤ 5e-13`, the existing trumpet
binary fixture tolerance. The complete per-group maxima and worst point IDs are
in [results.md](results.md); [parity.log](results/parity.log) is the raw output.

The parser creates a small zero-correction synthetic file with the same bound
maps and 3×3 coefficient blocks. It accepts that control and rejects 30 mutants:
wrong magic/convention/schema; missing, duplicate and unknown keys; block counts
and dimensions; nonfinite metadata/coefficients; noninteger controls; invalid
finite topology; malformed END/trailing data; profile hash; file mass/coupling;
truncation; checksum and final newline; and ten independent run bindings
(two masses, centres, rapidities and four couplings). Two inner-domain controls
reject the puncture and `R=1e-22`; a translated-coordinate control distinguishes
opposite axial offsets of `1e-18` near x=-16.

Nine SHA known-answer vectors cover empty input, `abc`, the standard 56-byte
message, one million `a` bytes, high-bit binary input and zero-byte messages
at lengths 55/56/57/64. The one/two-block examples and padding-boundary binary
vectors come from NIST's [SHA-256 examples](https://csrc.nist.gov/CSRC/media/Projects/Cryptographic-Standards-and-Guidelines/documents/examples/SHA256.pdf)
and [additional SHA-2 data](https://csrc.nist.gov/CSRC/media/Projects/Cryptographic-Standards-and-Guidelines/documents/examples/SHA2_Additional.pdf).
The implementation also reproduces the manifest's 11 MB file digest.

Before the original T7 reader changes, the e1ee1b8 fixture executable was rebuilt with only
all-table fingerprint printing added. At that stage, all stdout compared byte-for-byte
identical: 13 numeric bit fingerprints, all fixture tolerances, and the nine old
parser rejections. This covers the single-object, density-binary, radial/jet and
B/E echo paths. The legacy `.dat` initializer is untouched. The baseline and
current logs are retained. The later CCZ4 fingerprint refresh is described
in the direct-conversion validation section above.

## Fixed-mask t=0 constraints

`EMSCTTGridConvergence` extends the existing binary harness: SetValue/CCZ4 setter,
analytic ghost fill, `GammaCartoonCalculator`,
`ExperimentalGauge`, `Constraints<CouplingFunction>` from `ConstraintsCartoon`,
and `EMSCartoonGaussConstraints`. The underlying operators are unchanged.
No chi/lapse clipping is applied. The masks have zero cells below the actual
constraint operator's chi floor (`1e-6`) or the lapse floor (`1e-12`).

For each refinement, the core box is `[-18,18] × [0,2]`, with spacing
`1/(32 refinement) M`; the outer box is `[-96,96] × [0,96]`, with spacing
`4/refinement M`. This fixed two-box spacing ratio was declared before the first
run. Only spacing changes. The masks are:

- Left/right collars: `0.75 ≤ hypot(cosh(η)(x∓16),y) ≤ 1.5 M`, outside the authenticated reference horizon `R_h=0.63593977642346233 M`.
- Inter-hole: `|x| ≤ 14 M`, `0 < y ≤ 1.5 M`.
- Axis subsets of those three masks: fixed `0 < y ≤ 0.125 M`.
- Outer shell: `64 ≤ hypot(x,y) ≤ 96 M`.

The reported L2 is the unweighted RMS over each mask, as in the original harness;
L∞ is the maximum absolute value. Momentum components are the code's covariant
coordinate diagnostics, with `M = hypot(Mx,My)` as a convenient additional norm.
The maximum difference between the initialized grid and direct CCZ4 packing
is also recorded. Both Gauss constraints are reported; magnetic Gauss is exactly
zero for this axisymmetric data.

The expected order is four before the declared Float64 and representation floors.
The companion's D¹C interface-crossing floor is about `1e-11` in the scaled measure
of the Julia format document. That is not a universal bound on these unscaled
mask RMS values. This ladder does not claim to resolve that particular momentum
plateau; it retains all measured values and orders. The unconsumed D²C limitation
and the controller's accepted Float64-level end/end log ψ residual remain part
of the input certification. No acceptance threshold has been changed.

See [results.md](results.md) for every mask's H, both momentum components, combined
M and Gauss norms and orders at all spacings, CTT and density side by side;
[convergence.csv](convergence.csv) preserves full precision and L∞ results.

The original binary harness collar began at `R=0.5 M`, inside that horizon.
The final tables use the `.75` exterior cutoff from the existing evolution harness;
the earlier `.5` diagnostic runs are retained in `results/inner-collar-control/`.
This correction satisfies the requested exterior mask; it changes no initial field
or acceptance tolerance.

## Files changed

Fork production files:

- `Source/utils/SHA256.hpp` (new).
- `Source/InitialConditions/EMSBH/1D_SOL/EMSCTTSolution_read.{hpp,impl.hpp}` (new).
- `Source/InitialConditions/EMSBH/1D_SOL/EMSTrumpetSolution_read.{hpp,impl.hpp}` (hash the parsed profile snapshot).
- `Source/InitialConditions/EMSBH/EMSBH_trumpet_read.{hpp,impl.hpp}` (optional companion loading, evaluation and assembly).
- `Source/InitialConditions/EMSBH/EMSBHParams.hpp` and `Examples/EMS/SimulationParameters.hpp` (optional parameter).

Fork tests and evidence:

- `Tests/EMSTrumpet/EMSTrumpetFixtureTest.cpp` (print all existing table fingerprints).
- `Tests/EMSCTT/GNUmakefile`, `EMSCTTFixtureTest.cpp`, `EMSCTTGridConvergence.cpp`, `report.py`, this README,
  `results.md`, `convergence.csv`, and `results/` logs/authentication/smoke input.

EMS.jl worktree, exactly two text edits:

- `docs/binary-ctt-format.md`: first paragraph replaced with the controller T6c ruling.
- `test/fixtures/emsctt1/manifest.json`: only `certification` replaced with that same ruling.

No commits, state/AGENTS/GOAL edits, companion edits or Julia source changes.

The builds also leave standard generated `o/`, `d/` and `*.ex` products under
`Tests/EMSTrumpet`, `Tests/EMSCTT` and `Examples/EMS`. These are uncommitted build
products, not additional source changes.
