# Native EMS radiation instruments

Local instrument gates pass. The fork has a radiation-enabled pilot path; an
energy-balance certificate still requires a resolved binary run, an independently
certified initial ADM energy, and a settled branch-identified remnant. No commit,
SSH, HPC run, or EMS.jl edit was made.

The implementation follows the supplied horizon ruling, sections 3.1–3.6. All
operational arithmetic is Float64. Raw local runs are in
`/private/tmp/ems-radiation-t3b`; compact results, hashes and gates are in
[results/gates.json](results/gates.json), with the analytic checks in
[results/analytic.json](results/analytic.json). These temporary raw runs can be
regenerated using the commands below.

## Files and activation

| Files | Purpose |
|---|---|
| `Source/Extraction/EMSNativeRadiation.hpp` | Native Phi2 and scalar/EM worldtube stress flux |
| `Source/Extraction/EMSRadiationExtraction.hpp` | AMR interpolation, x-polar harmonics, measured sphere geometry, CSV output |
| `Source/Extraction/EMSRadiationParameters.hpp` | Default-off controls and shared RH threshold reader |
| `Examples/EMS/{GNUmakefile,Main_EMSBH2DBH.cpp,EMSBH2DLevel.hpp,EMSBH2DLevel.cpp}` | Guarded hooks in the existing executable, initial sample, wave-zone refinement |
| `Examples/EMS/ems_radiation_postprocess.py` | Retarded time, luminosities, energies, radius extrapolation |
| `Examples/EMS/params-radiation-pilot.txt` | Complete proposed pilot parameters; passed parameter-only dry run |
| `Tests/EMSRadiation/{README.md,GNUmakefile,harness/*,check_analytic.py,check_runs.py,run_local.py,verify_algebra.py,results/*}` | Reproducible tests and compact evidence |
| `Tests/EMSRHFinder/harness/EMSRHSurfaceTest.cpp` | Three added lines to exercise the production threshold reader in our existing harness |

The existing EMS executable now owns the radiation hooks. The author's post-step
body is retained verbatim, followed by one guarded radiation call at the minimum
extraction level's cadence. The tagging hook adds spherical extraction/wave-zone
floors only when native radiation is active; the legacy tagging still runs first.
Main sets the RH threshold and prepares the initial radiation sample after the
interpolator is attached. All new parameters are read in our own parameter header.
The existing emstrumpet1/ems_ctt_data_path initialization and evolution stay intact.
There is no separate radiation executable or copied post-step implementation.

| New parameter | Default | Meaning |
|---|---:|---|
| `ems_radiation_activate` | 0 | Activate all new radiation measurements and the t=0 sample |
| `ems_radiation_phi_inf` | 0 | Native asymptotic scalar value |
| `ems_radiation_wave_level` | 0 | Additional wave-zone refinement floor; zero disables this extra floor |
| `ems_radiation_wave_radius` | 0 | Radius of that floor, with a 20% buffer |
| `ems_rh_expansion_threshold` | 1e-7 | Assign RHUnion's public `m_thresh_super_low`; finite values in (0,1e-7] |

Radiation activation requires `activate_extraction=1` and reuses its sphere
parameters. Each requested extraction level is also enforced as a refinement floor with the author's 20% buffer. The additional
wave-zone floor preserves the legacy physical/mass-sphere tagging.

The new output is `data_subpath/ems_radiation.csv`, with every allowed
axisymmetric mode through l=6: S=R(phi-phi_inf), D=R partial_t(phi), P=R Phi2,
and W=R Psi4. The scalar monopole is retained. All m != 0 vanish in this cartoon
model; the harmonic polar axis is x. R is measured sqrt(A/4pi), while `r` labels
the fixed coordinate sphere. Both scalar and electromagnetic stress integrals
are recorded; their sum is L_matter. Their measure is alpha sqrt(det gamma)
r^2 dOmega, including the full determinant rather than assuming det(h)=1.
Area-weighted alpha^2 and Pi RMS are also recorded. No evolution/diagnostic
variable slots were added.

## Normalization audit

Finding about the author’s code, **not fixed here**: legacy `Pheyl2` output is
undefined and cannot be assigned a reliable luminosity coefficient. Its
constructor accepts the coupling but does not initialize `m_coupling_params`
(`Source/Matter/Pheyl2.hpp`, lines 69–82); the producer reads that uninitialized
member to multiply by exp[-alpha(f0+f1 phi+f2 phi^2)] (`Pheyl2.impl.hpp`, lines
142–148). The intended factor is sqrt(F), but its actual value is undefined.
There is no justified constant correction to existing Pheyl files.

Its intended spatial frame begins with u=(x,y,0), v=(-y,x,0), and
w proportional to gamma^{-1} epsilon(v,u); Gram–Schmidt normalizes v, then u,
then w (lines 174–223). In flat space w=-e_phi. The pre-coupling expression is
the conjugate of the standard x-polar outgoing Phi2 with that opposite azimuthal
orientation. `PheylExtraction.hpp`, lines 54–70, already multiplies its output
by radius and projects with spin -1. These files are preserved and excluded
from the new energy calculation.

The fork-only producer instead uses the sphere normal
s=gamma^{-1}dr/|dr|, normalized e_theta=(-sin(theta),cos(theta),0), and
e_phi=(0,0,1/sqrt(gamma_zz)). Thus

    k = (n-s)/sqrt(2), m = (e_theta+i e_phi)/sqrt(2),
    Phi2 = F_ab conjugate(m)^a k^b
         = [(E_theta+B_phi) + i(B_theta-E_phi)]/2.

These are native covariant electric and magnetic fields. There is no coupling
rescaling inside Phi2. The implementation is pinned by
`EMSNativeRadiation.hpp`, lines 41–51. Therefore the exact luminosity multiplier
on the squared radius-multiplied modes is **2 F_inf**:

    F_inf = exp[-2 ems_alpha (ems_f0 + ems_f1 phi_inf + ems_f2 phi_inf^2)].
    L_EM = 2 F_inf sum_l |P_l0|^2.

The independent outgoing EM test uses F_inf=exp(-0.4), rather than relying on
F_inf=1 to hide a missing factor. An incoming control has vanishing outgoing
Phi2 and negative outward stress flux. General lapse/shift stress identities
were checked by exact rational witnesses and 24 independent numerical tensor
contractions; [cas.log](results/cas.log) and
[tensor-check.log](results/tensor-check.log) record the scope. The unchanged
algebra witness is retained from t3a; numerical instrument checks were rerun
after folding the hooks into EMS.

GW extraction uses the unchanged WeylOmScalar producer. Its two 1/4 curvature
contractions and transverse difference (`WeylOmScalar.impl.hpp`, lines 101,
124–126) yield Psi4=-h_plus'' for the supplied outgoing TT-curvature control
(measured -0.13 for h_plus''=0.13). This pins the magnitude convention used with
L_GW=sum|integral W du|^2/(16pi). It is a tetrad-factor test; it does not certify
finite-radius gauge or matter errors of the Weyl reconstruction.

## Measured tests

| Test | Expected | Measured / error | Convergence |
|---|---|---|---|
| Exact-surface outgoing scalar Gaussian, A=.03, sigma=.5 | E_phi=0.00319041693163, L_phi>=0 | E_phi=0.00319041695093 at 129 angles; relative error 6.05e-9; peak absolute L error 1.61e-11 | Angular energy errors 1.57e-6, 9.69e-8, 6.05e-9 relative: fourth order |
| Same pulse through production AMR interpolation, r=4,6,8 | Same energy; correct instantaneous flux | Finest dx=.09375: E_phi=0.00319010858, .00319010390, .00319004988; worst energy error 1.15e-4 relative, peak-flux error 2.61e-4 | Energy orders 4.05, 4.11, 3.73; stress-energy error at most 1.28e-4 |
| Known scalar Y_l0, l=0..6 | Coefficient .03, all other coefficients zero | Maximum coefficient/leakage error 6.27e-9 at 257 angles; D normalization also passes | All angular refinements reduce errors by more than 8; approximately fourth order |
| Outgoing transverse EM dipole, E_theta=.02 sin(theta)/20 | L=0.00449252806865 | NP .00449252791457, stress .00449252798708; relative NP error 3.43e-8 | 129-angle Simpson integration |
| Incoming EM control | Outgoing L_EM=0; negative stress flux | Zero within Float64 amplitude tolerance; stress=-.00449252798708 | Direction/sign control |
| Certified stationary trumpet, v=0, 0<=t<=3M, r=4,6,8 | L_phi=L_EM=0 despite Pi!=0 | Worst coarse peaks 1.74e-12 / 1.18e-12; at dx=.125 both <1e-14; initial Pi_RMS(r=4)=1.459e-4 | Approximately fourth order in amplitude, eighth in squared luminosity |
| Curved metric/lapse/shift flux | Agreement with direct -T^i_t contraction | 24 cases, absolute discrepancy <1e-14 | Float64 check |
| GW TT-curvature tetrad | Psi4=-.13 | -.13 | Magnitude-factor check |
| Fixed-u linear/quadratic synthetic waveform extrapolation | Known infinity waveform, including scalar monopole | Both matching polynomial orders recover the answer to Float64 tolerance | Linear fit on quadratic data deliberately differs by 0.00462 in peak luminosity |
| Two-level AMR callback/regrid smoke | Extraction at fine-level cadence from t=0 | 9 samples/radius, dt=.0625, through t=.5 | All modes and geometric fields finite |
| Default off / explicit off | Original bytes | Five files identical, including both HDF5 checkpoints and all three legacy diagnostic files | Exact SHA-256 record |
| Existing fork regressions | Prior fixture fingerprints and successful grid runs | Reference and echo B/E fixtures, CTT checks, all 26 grid runs pass | Four asserted fingerprints preserved |
| Original RH pilot, new threshold key absent | Native default stopping and found surfaces | All four existing cases pass: both boost signs, RN, and boosted N=192 at dx=1/64 | 0.85–74.3 s/case; [record](results/native-rh-pilot.json) |

The full pulse grid table is [pulse-grid.csv](results/pulse-grid.csv).

The zero stationary test is an actual CCZ4/EMS evolution in Killing coordinates:
the test-only level restores the certified physical lapse and freezes lapse,
shift and gauge B after evaluating the normal production RHS. Every physical
field evolves. This matters because the production setter deliberately stores
lapse=sqrt(chi) (`EMSBH_trumpet_read.impl.hpp:339`) and uses ExperimentalGauge.
That separate, unchanged production-gauge control shows a finite-radius scalar
signal which persists under grid refinement (raw D-based peak about 1.84e-6 at
r=4). A coordinate waveform at small radius does not make that signal physical
radiation at infinity. It is retained in `gates.json`, not relabelled as a failed
stationary cancellation or filtered away.

At N_theta=48, dx=1/32 and chase speed .125, the unchanged finder gives:

| New threshold | Native mean-square expansion | Area | Updates (at most 400 chase steps each) |
|---:|---:|---:|---:|
| 1e-7 | 9.99938984e-8 | 38.5313939 | 98 |
| 1e-10 | 9.99852902e-11 | 38.5178208 | 178 |
| 1e-12 | 9.99670894e-13 | 38.5174221 | 197 |

Fresh interpolation of the final in-memory surfaces also passes each threshold.
The remaining area bias is principally the known angular quadrature error;
changing the stopping threshold is not an area/mass accuracy certificate.
The old test-only key remains supported, and the new production key takes
precedence when supplied.

## Post-processing

From the fork root, with NumPy installed (or using the existing uv runner):

```sh
uv run --with numpy python Examples/EMS/ems_radiation_postprocess.py \
  /path/to/ems_radiation.csv /path/to/energy --mass 2.012499 --du 0.125
```

The example mass must be replaced by the independently certified ADM value of
the particular binary input. Neither M_ADM.dat nor MQ mass is read. Output is
one `r*.csv` per radius, `infinity_order1.csv`, `infinity_order2.csv`, and
`summary.json` with input hash, integration assumptions and actual u windows.

The clock is t_corr=integral sqrt(<alpha^2>_A/(1-2M/R)) dt, and
u=t_corr-[R+2M log(R/(2M)-1)]. M=0 uses r*=R. The scalar derivative is
(D+Rdot S/R)/(du/dt), so changing areal radius is retained. The processor also
writes a finite-difference derivative-of-S luminosity as a comparison. It
integrates W to news, then evaluates the three native luminosities. Matter
stress energy is integrated in coordinate time. Extrapolation fits complex
news, scalar derivative and Phi2 amplitudes in 1/R at fixed common u, then
squares the infinity amplitudes.

Use `--coordinate-time` for a clock/gauge comparison and `--radii 50 75` (or
another subset) for radius sensitivity. Three radii give an exactly determined
quadratic fit with no redundancy; a conservation run needs at least four.
`--initial-news FILE.json` supplies seven [real,imag] constants per coordinate
radius (keys such as `"50.0"`); otherwise news is explicitly assumed zero at the
first sample. No high-pass, drift subtraction or scalar-monopole removal occurs.
An appended restart supersedes the old future; incomplete seven-mode blocks
are rejected. The processor does not silently clamp invalid R<=2M or
nonmonotone retarded time.

## Proposed pilot and limits

Subsequent finder performance, restart validation and the measured box-size choice
are recorded in [the performance report](../EMSRHFinder/PERFORMANCE.md).

The complete [pilot parameter file](../../Examples/EMS/params-radiation-pilot.txt)
passed the executable's parameter-only dry run. Its essential block is:

```text
L = 1536                 # half-extent 768 M about center = 768 0
N1 = 768
N2 = 384                # coarse dx=2; cartoon y>=0
max_level = 6           # finest dx=1/32
stop_time = 250
ems_radiation_activate = 1
ems_radiation_phi_inf = 0
activate_extraction = 1
extraction_center = 768 0
num_extraction_radii = 3
extraction_radii = 50 75 100
extraction_levels = 2 2 2
num_points_theta = 257
num_points_phi = 2
ems_radiation_wave_level = 2
ems_radiation_wave_radius = 110   # dx=.5 through r=132 with buffer
RH_num_horizons = 3
RH_initial_radii = 0.63593977642346233 0.63593977642346233 4
RH_initial_centre = 752 784 768
RH_num_points = 96 96 96
RH_level = 3 3 3
RH_start_times = 0 0 70
RH_time_step_freq = 400 400 400  # chase iterations/callback, not output stride
RH_chase_speeds = 0.125 0.125 0.125
ems_rh_expansion_threshold = 1e-12
```

This proposes the certified separation-32, opposite-rapidity CTT binary and
the exp-0004 surface guesses: two individual surfaces from t=0 and a radius-4
common surface at the centre from t=70, all called at level 3. All three
use 96 points, chase speed .125 and expansion threshold 1e-12. MQ extraction
remains on for its charge; its mass column is not used. `mod_F` is included in
`plot_vars`. Both initial-data paths are `/HPC_STAGING/` placeholders for staging.
As in the native finder workflow, accepted coincident common surfaces must
subsequently be deduplicated. A later mass precision comparison should use 192 angular points; a wave comparison should
raise the wave-zone/extraction floor to level 3 independently of inner AMR.
The stricter threshold's tracking cost/lag in a moving binary is not established
by the isolated tests.

At the outer extraction sphere the nearest physical outer boundary is 668 M
away. A conservative assumed incoming coordinate-speed bound of 2 gives a
first-contamination time of 334 M, exceeding 250 M by 84 M. This addresses
initial boundary errors, not just a reflected merger wave. ExperimentalGauge's
asymptotic lapse characteristic is about sqrt(1.8); the bound 2 still needs
checking against the actual lapse, shift and characteristic propagation in
the pilot. The y=0 cartoon reflection surface is a symmetry axis, not an outer
boundary. By contrast, half-extent 512 leaves only 206 M at speed 2, too short
for this conservative 250 M window.

Local builds and each actual test run took under five minutes; the longest
individual evolution was 84.7 s and each existing grid regression was below
29 s. No >30-minute local simulation or strict full finder sweep was launched.
The proposed binary run was not initialized/evolved/submitted, and MPI was not
tested with the local serial Chombo. Finite-radius gauge effects, wave-zone and
radius convergence, l>6 truncation, GW integration constants, radiation already
outside the first spheres, late tails and remnant settling remain pilot/balance
questions. A complete ADM-at-zero ledger must account for those missing windows;
the common-u extrapolated energy is explicitly only the common recorded window.

## Reproduce

Run long commands in the background and poll with backoff. The local driver
caps every build/evolution/finder job at 240 s; its short preliminary evolution
precedes the three-M sequence.

```sh
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/builds build
uv run --with sympy python Tests/EMSRadiation/verify_algebra.py > /private/tmp/ems-radiation-t3b/cas-final.log
uv run --with numpy python Tests/EMSRadiation/check_analytic.py /private/tmp/ems-radiation-t3b/analytic-verified
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/pulse-grid2 pulse-grid
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/killing stationary
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/evolution2 evolve
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/off-final off
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/amr amr
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/threshold threshold
python3 Tests/EMSRHFinder/run_cases.py /private/tmp/ems-radiation-t3b/native-rh-pilot pilot
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/regressions fixtures
python3 Tests/EMSRadiation/run_local.py /private/tmp/ems-radiation-t3b/parameters parameters
Examples/EMS/Main_EMSBH2DBH2d.*.ex Examples/EMS/params-radiation-pilot.txt check_params=1 > /private/tmp/ems-radiation-t3b/pilot-dryrun.log
uv run --with numpy python Examples/EMS/ems_radiation_postprocess.py /private/tmp/ems-radiation-t3b/pulse-grid2/pulse-grid-256/ems_radiation.csv /private/tmp/ems-radiation-t3b/postprocess --mass 0 --du 0.05
uv run --with numpy python Tests/EMSRadiation/check_runs.py /private/tmp/ems-radiation-t3b
```

The report collector also copies the CAS log, pilot dry-run log and CLI example
summary from that directory. Production execution uses
`Examples/EMS/Main_EMSBH2DBH2d.*.ex`. The `off` check exports and builds
fork main `4bf79b8e829159f5edbfb2a769e1fcbe1ad888f7` separately before comparing
default-off and explicit-off checkpoint/legacy-diagnostic bytes.

HDF5 object metadata contains whole-second timestamps. The initial comparison
crossed a second boundary: four timestamp bytes differed in each checkpoint,
while `h5diff` found no dataset differences and the three diagnostics were
identical. The rerun started the tiny cases at a second boundary and passed raw
byte equality for all five files; no checkpoint bytes were rewritten or masked.

## Files to commit

The inventory above includes the inherited t3a work. Stage only these files:

```text
Examples/EMS/EMSBH2DLevel.cpp
Examples/EMS/EMSBH2DLevel.hpp
Examples/EMS/GNUmakefile
Examples/EMS/Main_EMSBH2DBH.cpp
Examples/EMS/ems_radiation_postprocess.py
Examples/EMS/params-radiation-pilot.txt
Source/Extraction/EMSNativeRadiation.hpp
Source/Extraction/EMSRadiationExtraction.hpp
Source/Extraction/EMSRadiationParameters.hpp
Tests/EMSRHFinder/harness/EMSRHSurfaceTest.cpp
Tests/EMSRadiation/GNUmakefile
Tests/EMSRadiation/README.md
Tests/EMSRadiation/check_analytic.py
Tests/EMSRadiation/check_runs.py
Tests/EMSRadiation/run_local.py
Tests/EMSRadiation/verify_algebra.py
Tests/EMSRadiation/harness/EMSRadiationGridTest.cpp
Tests/EMSRadiation/harness/EMSRadiationTest.cpp
Tests/EMSRadiation/harness/EMSStationaryRadiationTest.cpp
Tests/EMSRadiation/results/analytic.json
Tests/EMSRadiation/results/cas.log
Tests/EMSRadiation/results/ctt.log
Tests/EMSRadiation/results/gates.json
Tests/EMSRadiation/results/native-rh-pilot.json
Tests/EMSRadiation/results/parameter-checks.json
Tests/EMSRadiation/results/pilot-dryrun.log
Tests/EMSRadiation/results/postprocess.json
Tests/EMSRadiation/results/pulse-grid.csv
Tests/EMSRadiation/results/tensor-check.log
Tests/EMSRadiation/results/trumpet.log
```

Do not stage `.ex`, `o/`, `d/`, `__pycache__/`, or raw run outputs. The removed
`Examples/EMSRadiation` files were untracked, so no tracked deletion is required.
