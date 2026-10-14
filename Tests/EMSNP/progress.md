# Local command record — 2026-10-15

All paths below are relative to the worktree unless marked absolute.
No subagents, remote compute, transfers, or evolution steps were used.

1. `git status --short`, `git log -3 --oneline`, `rg`/`sed`/`cat` of the
   reconstruction, stress tensor, production scheduling and existing test
   harnesses. Initial worktree clean at d10d27e; branch np-scalars.
2. Python hashlib SHA-256 streaming over the five supplied absolute paths:
   identities in `results/inputs.json`; trumpet header and TOML inspection.
   `h5dump -a time -a num_levels -a num_components <checkpoint>`:
   `results/checkpoint-header.log`, time 602, 16 levels, 28 variables.
3. Wolfram/Context7 tool discovery: Wolfram MCP evaluator unavailable;
   official Context7 reference selected. Literature Eq. (3) checked at
   https://arxiv.org/html/2312.17160v2. This supports the binomial laws;
   the polynomial script and tensor replay carry the actual controls.
4. `CHOMBO_HOME=... make -C Tests/EMSWeylOffline all DIM=2 -j 8` under
   `/usr/bin/time -p`, backgrounded with `nohup`, parent waiting for exit:
   `results/build.log`, 6.81 s, exit 0. The first detached launch left an
   empty log; reran with its parent kept alive. Poll after 15 seconds;
   completion seen on the next inspection. An attempted `ps` was denied
   by the sandbox; completion is recorded by explicit exit markers.
5. In `Tests/EMSNP/results`, `OMP_NUM_THREADS=8
   ../../EMSWeylOffline/EMSWeylCheckpoint2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex
   params.txt`, timed and backgrounded with parent waiting. 7.65 s, exit 0;
   both forms, all existing production radiation spheres, one jet each.
6. `python3 Tests/EMSNP/check_rotation.py <vacuum jet> <matter jet>`:
   self-check passed; `results/rotation.json`. Also `python3 -m py_compile`
   passed; generated pycache removed.
7. `make -n -C Examples/EMS all DIM=2`: exact compiler flags in
   `results/production-dry-build.log`. A Python driver copied
   `Source/GRChomboCore/*.hpp` to a temporary directory outside the worktree,
   asserted line 51 of GRAMRLevel.hpp was private, changed only that line
   to protected in the copy, prepended its include directory and
   `-fsyntax-only`, removed `-c`/`-o`, and invoked g++-16 on
   EMSBH2DLevel.cpp. `results/production-syntax.json`: exit 0, 1.42 s.
   Temporary directory removed; source untouched.
8. Repeated the checkpoint harness from `results/default`, using the saved
   parameter copy without the three np_probe controls. 5.25 s, exit 0;
   `results/default.log`. Standard-library CSV comparison: 21 rows per
   form, all columns within 1e-12 relative (+1e-100 for zeros), maximum
   difference 0, no np_jet files. `results/default-comparison.json`.
9. `/usr/bin/time -p /opt/homebrew/bin/gtimeout -k 5s 60s
   /Applications/Wolfram.app/Contents/MacOS/wolfram -script
   scripts/cas/np_null_rotation.wls`: cold CLI kernel, final rerun .75 s,
   exit 0; exact coefficients and counterexample in `results/cas.log`.
10. `git diff --exit-code HEAD -- Source Examples/EMS/*.cpp
    Examples/EMS/*.hpp`: no production source changes. `git diff --check`
    passed for source, docs and structured results; raw program logs retain
    the program's trailing spaces. Generated object/dependency directories and executable removed;
    all evidence logs/CSVs retained. Reviewed staged diff before commit.

Stop reason: (b) cannot be the stated fixed-field quadratic-order claim
while retaining the exact defining products. No production NP integration
or unperformed static/sphere validation is claimed.

## NPTET T2 — 2026-10-15

Starting revision ef0e2d5, branch np-scalars, initially clean. Re-read the
T1 README, probe, Python contractions and CAS script, the production
reconstruction/stress/extraction/level paths, the static harness and
fixtures, BUILD notes and native Julia static/trumpet formulation and
reader. SHA-256 reauthentication of the five supplied inputs matched T1;
all five hashes were checked unchanged again after the controls.

Implementation scope: new Source/Extraction/EMSNPScalar.hpp and
EMSNPExtraction.hpp; default-false parser flag and production radiation
scheduling hook; T1 probe reduced to a wrapper around the production
geometry; existing frozen harness extended; new static/unit executables
in the same GNUmakefile; test-only derivative/order controls and scripts.
No diagnostic enum, plot/checkpoint schema, gauge, RH finder, archived
extraction class, or evolution law was edited. No new dependency.

The reproducible command is `zsh Tests/EMSNP/run_t2.sh`. It was run in
the background with nohup and its parent waiting. Polling started after
15 seconds; independent report work continued between observations.
Every individual command measured below five minutes; the final full
workflow measured 34.05 s. Final exit 1 is the reported physics-gate
failure, not a hidden harness crash. Full output: results-t2/run-final.log.

The script commands and final outputs are:

- `CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib make -C
  Tests/EMSWeylOffline all DIM=2 -j 8`: results-t2/build-release.log.
  Three executable targets build normally. Initial builds measured
  7.53–9.47 s; the final incremental rebuild was shorter.
- `python3 Tests/EMSNP/syntax_check.py`: make dry-run plus the exact
  prescribed throwaway-header -fsyntax-only invocation. Exit 0;
  production-syntax.json and syntax-driver.log retain flags/timing/output.
- `OMP_NUM_THREADS=1 EMSNPTest...ex <RN params>`: unit.log, PASS flat
  spacetime, complex Maxwell normalization/signs, chi-floor rejection.
- Fresh CLI `gtimeout -k 5s 60s /Applications/Wolfram.app/Contents/MacOS/wolfram
  -script scripts/cas/np_static_metric.wls` and np_first_order.wls:
  static-cas.log and first-order-cas.log, each under 1.1 s, all witnesses
  PROVED. The runner also requires a PROVED marker and rejects parse or
  NOT-CLOSED messages because a kernel syntax error can otherwise exit 0.
- `julia --project=/Users/auroradysis/Workspace/EMS
  Tests/EMSNP/static_oracle.jl <supplied trumpet> <oracle CSV> <radii>`:
  rn-oracle.log and ems-oracle.log, exit 0, roughly 1.5–3 s each.
  The Julia development skill was used for this test script. Context7
  confirmed generic Real/nested Dual derivative support in ForwardDiff
  v1.4.6. It also resolved the official Wolfram reference. The Wolfram
  context tool returned an internal error; the evaluator and fresh CLI
  worked. Evaluator results were written only by the sandboxed tool.
- `OMP_NUM_THREADS=8 EMSNPStatic...ex params.txt`: rn-n256/n512 and
  ems-n1024/n2048 folders, corresponding *-release.log files. t=0 only,
  exit 0; approximately .1–1.6 s per hierarchy. No low-resolution echo
  evolution or RHS step was launched.
- `OMP_NUM_THREADS=8 EMSWeylCheckpoint...ex params.txt`: checkpoint,
  cavity, default folders and *-release.log. Each loads the unchanged
  16-level t=602 hierarchy; no advancement/regrid. Both Riemann forms
  and requested radii are exercised. All three executables exit 0.
  Observed cold/cached restart wall times were roughly 5–14 s.
- Activation dependency: static executable with EMS radiation false,
  NP true, expected exit 2 and the precise dependency error; PASS in
  activation-dependency.json. No hierarchy is constructed on this path.
- `python3 Tests/EMSNP/report_t2.py`: report.log, validation.json and
  summary.md. 70 static checks PASS; 980 independent harmonic mode
  components and pointwise X checks PASS; 42 unchanged legacy rows PASS.
  28/30 Weyl4 comparisons PASS, with two l=3,r=100 failures retained.
  All 28 matter rotation groups FAIL their roundoff-only transfer gate.
  ready_for_simulation=false; report exit 1. No tolerance was widened to
  accept the failed Weyl4 modes or rotation comparisons.

Exploratory/failed attempts are retained: build-final.log had an auto
declaration error in the refactored probe and an include-order exp
declaration error in the standalone unit. Both were fixed and rebuilt.
The first unit invocation printed usage without a parameter file; the
actual unit result above is from the required parameter invocation.
first-order-cas-syntax-failure.log and its MCP failure transcript record
a missing closing bracket, corrected before the passing cold-kernel run.
default-lock-failure.log records an attempted concurrent restart rejected
by HDF5's exclusive lock; subsequent/final restart runs are sequential.

An exploratory Lagrange<6> **jet** reconstruction disagreed strongly with
the otherwise controlled curvature; its CSVs remain under
exploratory-sixth-jets and are not used as accuracy references. It was
not investigated or promoted. The final sixth-order cell-Weyl control
uses the same existing Lagrange<4> interpolation as production, isolating
the cell derivative order and avoiding that unvalidated route.

Static outermost-EMS curvature lacks fourth-order decrease once errors
reach the Float64 second-derivative floor; the report explicitly separates
the arithmetic and Richardson allowances. Required method was preserved:
no projection using reconstructed Ricci replaced field-equation Ricci.
No full production link/evolution was attempted past the specified syntax
workaround. No HPC, SSH, file transfer, push, rebase or branch deletion.
