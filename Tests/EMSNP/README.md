# NP tetrad acceptance probe: stopped on control (b)

NPTET T1, run `ariadne-run:v1:51304_760a6373`, branch `np-scalars`,
starting revision `d10d27e`. **The requested production diagnostic is not
implemented and is not ready for a new simulation.** The fixed-field,
second-order-in-null-rotation requirement in (b) is stronger than the
first-order Coulomb cancellation of the defining expressions. Changing that
requirement was explicitly a stop in the task; this packet preserves the
failed control without changing the expressions, gauge, or production code.

## Convention and normalization

Signature is `(-,+,+,+)`. Use `WeylOmScalar::compute_null_tetrad` unchanged:
`u=m1` is the normalized coordinate radial direction; `v=m2` is the
Gram–Schmidt orthogonalization of `(-y,x)` in the cartoon plane; the physical
azimuthal unit vector is `w=mw ∂w`. The unused, degenerate `m3` is not used.
Let `n` be the future Eulerian unit normal and set

```
l=(n+u)/sqrt(2), k=(n-u)/sqrt(2), m=(v+i w)/sqrt(2).
l.k=-1, m.mbar=1.
Psi0=C(l,m,l,m),        Psi1=C(l,k,l,m),
Psi2=C(l,m,mbar,k),     Psi3=C(l,k,mbar,k),
Psi4=C(k,mbar,k,mbar).
Phi0=F(l,m), Phi1=[F(l,k)+F(mbar,m)]/2, Phi2=F(mbar,k).
X_out=2 Phi1 Psi3-3 Psi2 Phi2; X_in=2 Phi1 Psi1-3 Psi2 Phi0.
```

These positive Weyl contractions fix the signs of every Psi together. With
the evolved covectors, `F(i,n)=E_i` and `F_ij=sqrt(det gamma)[ijk] B^k`.
There is no coupling factor in any Phi. In an oriented orthonormal frame,
`Phi1=(E_u+i B_u)/2` and
`Phi2=[E_v+B_w+i(B_v-E_w)]/2`, agreeing with `EMSNativeRadiation` when its
transverse frame agrees. Its surface-normal frame need not equal this
radial-first Gram–Schmidt frame in a nonspherical metric.

In the reflection-even polar sector of the supplied hierarchy (`E_w=0`,
`B_u=B_v=0`, no metric component with one w index), all five Psi, all three
Phi, and both X are real. Their imaginary parts vanish; their real parts
are generally nonzero. Axisymmetry alone would not exclude an independent
magnetic sector. On an aligned spherical Coulomb solution only Psi2 and
Phi1 can remain. The requested spin weights are respectively `2,1,0,-1,-2`,
`1,0,-1`, and `-1,1`; no new harmonic projections were implemented.

The existing `c_Weyl4_Re` is **`+Re(Psi4)`**, factor **+1** in this convention.
To see this, expand the contractions:

```
Re Psi4 = [C(n,v,n,v)-2 C(n,v,u,v)+C(u,v,u,v)
          -C(n,w,n,w)+2 C(n,w,u,w)-C(u,w,u,w)]/4.
```

Replacing C by R leaves the transverse difference unchanged: each Ricci
subtraction in `C(k,e,k,e)` is `-Ric(k,k)/2`, the same for e=v,w.
The expression is exactly the `omlocij[1][1]-Omww` in
`compute_Weyl_Om`; its returned `out.Om22` already includes that subtraction.
This identity assumes the same physical metric in the reconstruction and
tetrad. The probe rejects chi at or below `1e-4`, so neither existing chi
floor changes that metric.

On the **same** t=602 hierarchy, at the valid level-2 cell
`(x-center_x,y)=(35.21875,35.21875)`, dx=0.4375, the absolute normalization
residual is `2.118e-22` for the vacuum form and `6.353e-22` for the
matter-complete form. The magnitude-based Float64 tolerance is `1.677e-16`
(`4096 epsilon` times the tensor norm and frame contraction scale).
This is a pointwise normalization check, **not** the requested sphere-mode
comparison at r=50,75,100.

Lengths are in the input's geometric units, G=c=1; `[Psi]=L^-2`,
`[Phi]=L^-1`, `[X]=L^-3`, and the null-rotation parameter is dimensionless.
No mass or field normalization is adjusted.

## Exact obstruction

For a real l-fixed null rotation of parameter a,

```
l'=l, m'=m+a l, k'=k+a(m+mbar)+a^2 l,
Phi0'=Phi0,
Phi1'=Phi1+a Phi0,
Phi2'=Phi2+2a Phi1+a^2 Phi0,
Psi1'=Psi1+a Psi0,
Psi2'=Psi2+2a Psi1+a^2 Psi0,
Psi3'=Psi3+3a Psi2+3a^2 Psi1+a^3 Psi0.
```

The binomial laws follow by multilinearity of the symmetric rank-2 and
rank-4 spinors. They also follow from direct tensor contraction; the
independent electric/magnetic Weyl self-check in `check_rotation.py`
agrees to `5.929e-17`. For a literature cross-check see Eq. (3) of
[Mao and Zhao, arXiv:2312.17160v2](https://arxiv.org/html/2312.17160v2#S2).
That paper uses the opposite metric signature; the binomial coefficients
are unaffected by the common sign conversion of the scalar definitions.

Substituting the defining expressions gives the exact linear coefficients:

| Fixed null vector | dX_out/da at 0 | dX_in/da at 0 |
|---|---|---|
| l | `2 Phi0 Psi3 - 6 Phi2 Psi1` | `2 Phi1 Psi0 - 4 Phi0 Psi1` |
| k | `2 Phi1 Psi4 - 4 Phi2 Psi3` | `2 Phi2 Psi1 - 6 Phi0 Psi3` |

These are not identically zero. For example, in mass-one units choose
`Psi2=-1, Phi1=1, Psi3=epsilon, Phi0=epsilon`, with all other Psi and Phi
zero. The l-fixed derivative of X_out is exactly `2 epsilon^2`, nonzero
for every nonzero epsilon, however small. This is an algebraic
counterexample, **not a replacement for the supplied checkpoint**.

The protected coefficients (l-fixed X_out, k-fixed X_in) are quadratic in
non-Coulomb perturbations. They cancel to first order about an aligned
Coulomb background. They do **not** in general vanish to first order in a
at fixed perturbed fields. A joint expansion in the data perturbation and
tetrad rotation, or a bound allowing the linear remainder, is a different
acceptance statement. Neither was substituted for the requested test.

## Frozen-checkpoint control

`CheckpointJet.hpp` reuses the existing Gauss/Codazzi/normal-normal
reconstruction and `add_matter`, exports its Riemann tensor, and obtains
`Ric_ab=8 pi G (T_ab-T g_ab/2)` from `compute_EMS_EM_tensor`. The test
reader subtracts that Ricci part and contracts the resulting tensor and
the rebuilt F with one tetrad before and after the rotation. It does not
apply NP transformation formulas to generate the checkpoint results.
The hierarchy is never advanced, regridded, or initialized from a static
solution. Only one selected valid cell is exported per form.

The matter-complete result gives:

| a (l fixed) | abs(delta X_in) | abs(delta Phi2) |
|---|---:|---:|
| 1e-3 | 3.391041889e-16 | 1.590913996e-7 |
| 1e-4 | 3.391040768e-17 | 1.590914453e-8 |
| 1e-5 | 3.391040826e-18 | 1.590914499e-9 |

The decade slopes for X_in are `1.000000144` and `0.999999993`; for
Phi2 they are `0.999999875` and `0.999999988`. A second-order response
would have slope 2. Over these steps, a slope tolerance of 0.003 follows
conservatively from the `<=1.452e-20` paired-change roundoff budget for X:
4096 Float64 eps times the absolute sums in the tensor contractions and
products, including both evaluations. The measured changes exceed this
budget by factors 234 to 23,400. This comparison diagnoses the rotation
order **on the frozen discretized fields**, not continuum convergence.
For a k-fixed rotation, abs(delta X_out) is `1.754826600e-13`,
`1.754823291e-14`, `1.754822960e-15` over those same steps, also linear.
Phi2 is invariant under the k-fixed rotation, as expected.

There is a further numerical qualification: using the field-equation
Ricci with discretized Gauss/Codazzi data leaves a residual Weyl trace.
Its orthonormal maximum is `3.036e-10` for the matter form (`8.059e-7`
for the vacuum reconstruction with matter Ricci). This is reported in
`rotation.json`; **no trace projection or change of reconstruction method
was made**. Consequently the ideal trace-free binomial derivatives and
the direct derivatives need not coincide exactly. In particular, the
very small l-fixed X_out changes at a<=1e-5 are below the conservative
roundoff budget and are not evidence of a resolved order. The exact
algebraic obstruction above does not depend on this discretization issue.

## Validation status and files

| Requested control | Result |
|---|---|
| (a) RN and alpha=.9,e=7 EMS static profiles | Inputs fingerprinted and metadata inspected; hierarchies, static formulas, and truncation tolerances **not run/derived** after the stop. |
| (b) Quadratic rotation response of both X on perturbed data | **Failed as written**; exact obstruction and direct checkpoint evidence above. Joint perturbative cancellation is a separate statement. |
| (c) Psi4 versus production W on all three spheres | Factor +1 derived and checked at one cell for both forms; full NP sphere projections **not implemented/tested**. Existing radiation extraction ran on r=50,75,100. |
| (d) Default production path | `EMSBH2DLevel.cpp` syntax check passed with the prescribed throwaway header. No Source/ or Examples/ source changes. Test probe disabled: 42 radiation CSV rows agree within 1e-12 relative, maximum difference 0; no jet files created. |

No production `ems_np_activate`, diagnostic variables, extraction class,
or `ems_np.csv` has been added. The only harness change is an optional
`np_probe_level` (absent means off) with required finite `np_probe_x,y`.
Existing radiation columns and production defaults remain unchanged.

`results/inputs.json` records absolute names, sizes, and SHA-256 fingerprints
of all five supplied inputs. The trumpet headers identify M=1,Q=.7 and
M=1,Q=1.0681492364093355 respectively, native conventions and the requested
couplings. The supplied TOML identifies e=7; the checkpoint header records
t=602, 16 levels, 28 variables. The original parameter file has no
`restart_file`; the saved test copy appends it and the three probe controls.
These are fingerprints of the supplied objects, not comparison with an
unprovided historical digest. The harness also verifies grid size, dx,
and level synchronization against the checkpoint.

Source changes: `CheckpointJet.hpp`, `check_rotation.py`, the optional
hook in `../EMSWeylOffline/harness/EMSWeylCheckpoint.cpp`, and the executable
exact witness `../../scripts/cas/np_null_rotation.wls`. This README,
`findings.md`, `progress.md`, and `results/` retain conventions and evidence.

## Reproduction and timings

From the worktree root:

```sh
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
make -C Tests/EMSWeylOffline all DIM=2 -j 8
cd Tests/EMSNP/results
OMP_NUM_THREADS=8 ../../EMSWeylOffline/EMSWeylCheckpoint2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex params.txt
cd ../../..
python3 Tests/EMSNP/check_rotation.py Tests/EMSNP/results/vacuum/np_jet.txt Tests/EMSNP/results/matter/np_jet.txt
/opt/homebrew/bin/gtimeout -k 5s 60s /Applications/Wolfram.app/Contents/MacOS/wolfram -script scripts/cas/np_null_rotation.wls
```

The output is append-only when reusing an existing radiation CSV; use an
empty run directory for another run. Both `vacuum/` and `matter/` must exist.
For a probe-off replay use `results/default/params.txt` in a fresh directory.

Full logs: `results/build.log` (exit 0, 6.81 s), `checkpoint.log` (exit 0,
7.65 s), `default.log` (exit 0, 5.25 s), `cas.log` (fresh Wolfram CLI kernel,
exit 0, .75 s), `production-dry-build.log`, `production-syntax.json`
(exit 0, 1.42 s, exact compiler command and temporary-header adjustment),
`checkpoint-header.log`, `rotation.json`, and `default-comparison.json`.
The two jets, both forms' existing radiation CSVs, and both parameter copies
are retained. No run was expected to exceed five minutes; no program
measuring over thirty minutes was launched. No evolution, HPC, SSH, file
transfer, upstream PR, or push was performed.
