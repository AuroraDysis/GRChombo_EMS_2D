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
