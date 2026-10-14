# Weyl4 on frozen checkpoints: vacuum and matter-complete

`EMSWeylCheckpoint` loads a production checkpoint (its parameter file with `restart_file` set), computes Weyl4 on every level in the
production order (ghosts, then the compute class on valid cells) twice — `WeylOmScalar`, the production form, which takes the
normal-normal Riemann projections from the vacuum field equations, and `EMSMatterWeylOmScalar`, which keeps the matter term
8 pi G [S_ij - gamma_ij (S - rho)/2] with the evolution's own EMS stress tensor — and after each runs the production
`EMSRadiationExtraction` on the parameter file's spheres. Output: `vacuum/ems_radiation.csv`, `matter/ems_radiation.csv`.
The hierarchy is never advanced or regridded.

```sh
export CHOMBO_HOME=/path/to/Chombo/lib
make -C Tests/EMSWeylOffline all DIM=2 -j 8
```

Check (2026-10-15, exp-0033 checkpoint at t = 602, r = 50 / 75 / 100): the vacuum form reproduces the production
`ems_radiation.csv` (Psi4 l = 2 to 2e-11 relative, Phi2 and phi to 1e-12 or better); matter - vacuum in Psi4 l = 2 at r = 50 is
-6.9e-11 against -2.9e-4, second order as expected. One checkpoint: 13 s, 1.3 GB, serial with 8 OpenMP threads.
Production switch: `ems_weyl_matter = true` (default false) selects the matter-complete form in `EMSBH2DLevel`.

Test-only NP null-rotation acceptance probe: set `np_probe_level` plus
`np_probe_x` and `np_probe_y` (coordinates relative to extraction center)
to export one valid cell's reconstructed tensors to each form's `np_jet.txt`.
Absent `np_probe_level`, this probe is off. See [the failed acceptance
control and convention note](../EMSNP/README.md); this is not a production
NP extraction.
