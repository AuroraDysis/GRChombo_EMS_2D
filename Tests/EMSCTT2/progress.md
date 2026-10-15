# Execution record

Original command paths and measured costs are retained in results/cost-*.json.
No remote-compute ledger applies: every operation was local.

- Read the EMSCTT/2 specification and Julia pulse/interpolation/reconstruction.
- Implemented a separate strict reader and optional single-hole setter branch.
- Executed `python3 Tests/EMSCTT2/authenticate.py`; payload/profile/provenance pass.
- Executed `uv run --no-project --with sympy python scripts/cas/ctt2_longitudinal.py`;
  eleven exact witnesses close (results/cas.log).
- Julia reference: `reference.jl` under the supplied wt-pulse project; all five
  files evaluated. Original Float64-AD taper failure and corrected 512-bit jets
  retained. `compare.py` and `compare_stencils.py` assert the numerical budgets.
- Built isolated t=0 hierarchies at three resolutions, global and near-hole;
  recorded valid-cell norms and production/analytic ghost comparison.
- Rebuilt EMS and old CTT/trumpet fixture suites after the source changes.
- No evolution, HPC, SSH, file transfer, push, upstream PR or input modification.
- Archived console logs have trailing whitespace trimmed; numerical tokens and
  output ordering are preserved. CSV uses LF line endings.
