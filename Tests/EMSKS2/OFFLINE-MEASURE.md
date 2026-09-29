# T13 offline checkpoint measurements

**Registered outcome: READY-EXCEPT the strict outgoing-expansion target on B and evolved surfaces.** The E (t=0) surface reaches (\langle\Theta_+^2\rangle<10^{-14}); B (t=0) and the two-step continuations retain `FLOOR` status. The failed offline Newton polish trial was removed. Its divergence is not a validation of a different root. These rows are useful measurements, but no strict-horizon acceptance or late-time energy-balance claim follows from a `FLOOR` row.

## Command

Build from this worktree with the toolchain in `EMS-deps/BUILD.md`:

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
make -C Tests/EMSRHFinder all DIM=2 -j4
uv run --with h5py --with numpy python Tests/EMSKS2/offline_measure.py RUN_DIR RUN-measure.csv \
  --points 48 96 192 --threads 4 --cap 240 --dry-run
uv run --with h5py --with numpy python Tests/EMSKS2/offline_measure.py RUN_DIR RUN-measure.csv \
  --points 48 96 192 --threads 4 --cap 240
```

`RUN_DIR` contains `params.txt` and `chk/*.hdf5`. `--checkpoints PATH ...` selects a synchronized subset, including checkpoints held in other directories. The output CSV path must be new. Its sibling `.work/` tree retains exact parameter copies, per-resolution `surfaces.csv`, `exterior.csv`, shapes, progress, logs, input hashes and status with wall time and child peak RSS. The command exits nonzero if a horizon fails its final (10^{-14}) target; the CSV and status remain available. Run once per run directory, not once per MPI rank. This is a serial offline instrument using the production AMR interpolator; MPI performance is unmeasured. Checkpoints are read without advancing the hierarchy. `RHSurf.hpp`, `RHUnion.hpp`, their native histories, and evolution code are unchanged.

Compact local rows are archived in [`t13-B.csv`](t13-B.csv) and [`t13-E.csv`](t13-E.csv): all three angular counts at (t=0) and (N=48) after two coarse steps. Their raw checkpoints and full iteration traces remain in `/private/tmp/t13-*`.

The puncture-centred seed uses `star_centre[0]` and `r_h` from the hashed KS2 file. It does not require a surviving native finder history. The underlying `offline.py` still accepts saved shape histories when `--sphere-seed` is absent. `--skip-mass` bypasses the alpha-20 horizon-family projection and enables these KS2 exterior measurements. The frozen finder uses chase speed 0.125, quota 400, stages (10^{-7},10^{-10},10^{-14}), and retains the best fresh interpolation when a stage stalls. The existing plateau detector and hard cap remain distinct statuses.

## Definitions and uncertainty

`Q` in `surfaces.csv` is exactly `RHSurf::Q_charge()`: (Q_S=\sqrt{8\pi}/(4\pi)\int_S F(\phi)E_i s^i\,dA=\int_S F(\phi)E_i s^i\,dA/\sqrt{2\pi}). The raw `Q_S` on exterior rows is the same integral. `A` and `R_areal=sqrt(A/(4π))` are the raw midpoint sums. `A_corrected` and `Q_corrected` multiply those raw integrals by `sin(π/(2N))/(π/(2N))`, removing the exact midpoint-sphere bias. Both raw and corrected values are retained; `R_corrected=sqrt(A_corrected/(4π))`. For nonspherical data the adjacent (N,2N) corrected difference provides the residual angular estimate, divided by three on the finer value (four thirds on the coarser value). It does not estimate grid or evolution error.

Each horizon row gives `expansion_squared` (area-weighted (\langle\Theta_+^2\rangle)), `theta_max`, and `theta_minus` (area-weighted ingoing expansion). The harness displaces the whole measured surface radially by (\pm\min(dx_0/8,r_{\min}/8)), reinterpolates, and computes (\delta r=\max_\theta |\Theta_+|/|\partial_r\Theta_+|), with the finite-difference derivative measured pointwise. `area_delta` and `charge_delta` multiply this worst displacement by the measured uniform radial derivatives of area and charge. The corrected versions apply the same quadrature factor. These are conservative first-order position sensitivities, not full grid error bars.

`exterior_coordinate` rows use the listed coordinate radii in units of (r_h). `exterior_areal` rows solve (R(\rho,\theta)=\rho\sqrt{h_{ww}/\chi}=R_\mathrm{target}) at each common angle, to relative root tolerance (10^{-12}). The list is (1.2,1.5,2,3,5,8,12,R_\mathrm{outer\ ring}/r_h,20,30). The optical-potential root computed from the KS2 source is 14.87445955 for B and 9.473392... for E, consistent with the EMS.jl light-ring table at the precision needed here. `m_mean` averages (m=R(1-\gamma^{ij}\partial_iR\partial_jR)/2) with physical area weights; the gradient includes both Cartesian derivatives of (h_{ww}/\chi). This is the **specified spatial-slice formula**. `phi_mean`, `lapse_mean`, `K_mean`, and `electric_E2_mean=\langle E_iE^i\rangle` use the same weights. Each `*_max_departure` is the maximum absolute difference from its weighted mean, so nonspherical structure stays visible. `Q_S` is a whole-surface flux. The angular charge diagnostic is (dQ/d\Omega=F E_i s^i(dA/d\Omega)/\sqrt{2\pi}), reported as its solid-angle mean and maximum departure. `m_static`, `phi_static`, `Q_static`, `lapse_static`, and `K_static` come from the unchanged KS2 source profile at the corresponding areal radius. Coordinate-sphere (R(t)/R(0)-1) is a separate drift observable; this test does not replace it by a fixed-coordinate field comparison.

## Local validation

The production executable wrote both four-level (t=0) checkpoints, with base spacing (M/64), finest spacing (M/512), and 1,434,624 active cells. The unchanged KS2 source has (r_h=1), hence exact (A_H=4\pi=12.566370614359172), with (Q_B=8) and (Q_E=5.9999999999966764). The table is the final (N=192) row; `δA_pos` is the radial-residual sensitivity. The area and charge errors are smaller than these sensitivities or the measured angular estimates. The position contribution is 1.14% (B) and 0.28% (E) of one tenth of the area target (10^{-3}A_H); the charge-position contribution is far below one tenth of (2\times10^{-4}).

| member | status | (A_\mathrm{corrected}) | (Q_\mathrm{corrected}) | (\langle\Theta_+^2\rangle) | max \(|\Theta_+|\) | (\delta r) | (\delta A_\mathrm{pos}) |
|---|---|---:|---:|---:|---:|---:|---:|
| B | FLOOR | 12.566370562918 | 7.999999965257 | 4.4542e-14 | 5.8104e-7 | 5.6921e-7 | 1.4306e-5 |
| E | FOUND | 12.566370602230 | 5.999999993982 | 3.3327e-15 | 1.5050e-7 | 1.4138e-7 | 3.5533e-6 |

At (N=192), the largest (t=0) absolute errors over all ten fixed-areal radii are:

| member | max \(|m-m_\mathrm{static}|\) | max \(|\phi-\phi_\mathrm{static}|\) | max \(|Q_\mathrm{corrected}-Q_\mathrm{static}|\) | max scalar nonspherical departure |
|---|---:|---:|---:|---:|
| B | 3.16e-9 | 2.60e-10 | 1.59e-8 | 1.26e-8 |
| E | 6.27e-9 | 4.07e-10 | 1.60e-8 | 1.60e-8 |

All (t=0) coordinate-sphere (R) values at the listed radii agree with the source to the shown interpolation accuracy. The (N=192) clean extraction sphere at (R=30r_h) has corrected charge error below (3\times10^{-10}) for both members. After two coarse steps, B at (t=0.06139276249) gives corrected (A=12.56608109576), corrected (Q=7.99999989726), and (\langle\Theta_+^2\rangle=6.34\times10^{-13}) (`FLOOR`); E at (t=0.04471072067) gives (A=12.56637060310), (Q=5.99999999402), and (1.55\times10^{-13}) (`FLOOR`). All reported invariant fields remain finite. Maximum coordinate-sphere areal-radius drifts across the ten radii are (3.16\times10^{-11}) (B) and (1.24\times10^{-12}) (E), below 2%. These are short local continuations, not a long-time drift test.

An independent EMS.jl `read_ks2`/`ks2_sample` readback at (R/r_h=1.2,2,30) agrees with the C++ static mass and scalar columns to Float64 roundoff. For example, at (R=1.2r_h), EMS.jl gives B ((m,\phi,Q)=(0.2821071452389381,0.5744398719572448,8)) and E ((0.28119838744303877,0.4731238464288072,5.999999999996676)). This checks the static reference reconstruction; the numerical checkpoint fields remain an independent interpolation.

## Cost and claim boundary

Observed wall times on one local process with four OpenMP threads: B (t=0) (N=48/96/192): 44.58/51.58/82.71 s, 179.37 s for the bundle; E (t=0): 0.93/1.32/0.42 s, 3.22 s including Python; B and E two-step (N=48): 68.16 and 44.36 s. Peak child RSS was at most 0.86 GB for these final runs. The large difference reflects whether the surface starts below target or reaches the plateau window. B/E production checkpoint generation took 45.02/41.15 s at (t=0); the two B continuation steps took 36.54 s and the E steps 25.32+21.80 s. These are serial local measurements, with another session using the machine, not cluster timings.

Scaling the measured B floor bundle by active-cell count relative to the 1,434,624-cell local hierarchy gives a **serial sizing estimate**, not a measured MPI rate:

| production case | active cells | ratio | floor-bundle estimate | driver's memory plan |
|---|---:|---:|---:|---:|
| B M/512 | 2,181,632 | 1.52 | 4.5 min | 5.3 GiB |
| B M/720 | 4,257,280 | 2.97 | 8.9 min | 8.5 GiB |
| B M/1024 | 8,560,640 | 5.97 | 17.8 min | 15.1 GiB |
| E | 3,639,296 | 2.54 | 7.6 min | 7.6 GiB |
| reference M/64 edge | 5,612,544 | 3.91 | 11.7 min | 10.6 GiB |

The memory plan is the existing `2 GiB + 6 × checkpoint field payload` heuristic, with payload scaled by active cells. It is not a measured production peak. B M/1024 exceeds the local half-RAM dispatch limit, so no local production checkpoint was run. A cluster pilot must time one real checkpoint before setting the batch allocation or asserting MPI scaling. The local (t=0) agreements validate formulas and units on these grids. They do not establish continuum error, late-time stationarity, a settled remnant, or the strict B horizon target.
