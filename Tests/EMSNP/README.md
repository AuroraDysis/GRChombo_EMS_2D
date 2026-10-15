# One-tetrad EMS Newman–Penrose extraction

NPTET T5b, `ariadne-run:v1:83992_07d05728`, resumes 56f9ff6 on `np-scalars`.
The review adds the fixed fourth-order cell-derivative error omitted by
T4. Historical records remain in README-T1.md,
README-T2.md and [README-T3.md](README-T3.md).

**Radiation validation passes through l=6, with switches default off.**
All 18 sphere/mode/form comparisons pass the corrected budget, even
without its permitted safety factor 2. Linked production executables
also pass the non-echo, three-step smoke described below.
The accepted T3 static, rotation, algebraic, normalization and default
controls remain valid; L4 output reproduction is checked again in T4.

## Production interface

Enable `ems_radiation_activate = true` and `ems_np_activate = true` to write
`ems_np.csv` on the existing radii, center, Simpson quadrature and schedule,
including the existing t=0 call. NP defaults false and requires EMS radiation.
If any angular node has invalid geometry, an active tetrad chi floor,
or a non-finite NP value, every mode of that sphere is NaN for that time.
One warning per sphere and time is printed, and other spheres and the
evolution continue. Consumers must treat these rows as unavailable data.
Startup parameter validation remains fatal on all ranks. NP output IO
errors are reported without throwing from rank zero.
Physical/model parameters are live inputs. No gauge, evolution, KO, AMR
transfer, or static target after t=0 is introduced.

`EMSNPScalar` reuses the inherited Gauss, Codazzi, normal-normal Riemann,
and `EMSMatterWeylOmScalar::add_matter`, including the evolution's own
`compute_EMS_EM_tensor`. It builds symmetric trace-free E and B in the
complete oriented x,y,w orthonormal frame and reconstructs an algebraic
Weyl tensor. The old full Ricci-subtracted tensor is reconstructed only
in the test control for radiative preservation.
`EMSNPExtraction` evaluates it on interpolated fourth-order spatial jets:
the existing Lagrange<4> first/second derivative stencils have five/six points.
The legacy W is still computed on cells and then interpolated.
**No new diagnostic variables or checkpoint/plot layout changes were needed.**

The new production extraction is always matter-complete. Legacy W retains
its independent `ems_weyl_matter` flag; enable that flag when comparing the
same matter-complete form. The frozen harness runs both forms. Its vacuum
control omits both the normal-normal matter correction and the Ricci
terms of the electric/magnetic projection, giving E=STF(Q) consistently.
It retains the evolved Maxwell fields as a diagnostic comparison on the
same hierarchy; it is not called a matter-complete geometry.

CSV columns are `t,r,R,l,m`, then real/imaginary pairs for `Psi0`…`Psi4`,
`Phi0`…`Phi2`, `X_out`, `X_in`. Each value is **R times its harmonic
coefficient**, matching EMS radiation. Every allowed m=0 mode through l=6
is recorded; l below the spin magnitude is zero and axisymmetry removes
m!=0. Both X are formed pointwise before projection. R comes from the
proper sphere area, r is coordinate radius. The two zero-measure polar
endpoints are omitted from reconstruction and pointwise controls because
the inherited formulas contain 1/y and 1/y². All 95 nonzero-weight nodes
are tested; no offset, replacement pole value, or new frame is introduced.

## Fixed convention and normalization

Signature `(-+++)`; n is the future Eulerian normal. Use the inherited
radial-first Gram–Schmidt vectors u=m1 and v=m2, and positive cartoon
azimuthal w=mw∂w. The unused degenerate m3 is not consumed.

```
l=(n+u)/sqrt(2), k=(n-u)/sqrt(2), m=(v+i w)/sqrt(2),
l.k=-1, m.mbar=1.
Psi0=C(l,m,l,m),       Psi1=C(l,k,l,m),
Psi2=C(l,m,mbar,k),    Psi3=C(l,k,mbar,k),
Psi4=C(k,mbar,k,mbar).
Phi0=F(l,m), Phi1=[F(l,k)+F(mbar,m)]/2, Phi2=F(mbar,k).
X_out=2 Phi1 Psi3-3 Psi2 Phi2,
X_in =2 Phi1 Psi1-3 Psi2 Phi0.
```

E,B are covectors: `F(i,n)=E_i`, `F_ij=sqrt(det gamma)[ijk] B^k`.
There is no coupling multiplier in Phi. In an oriented orthonormal frame,
`Phi1=(E_u+i B_u)/2`, `Phi2=[E_v+B_w+i(B_v-E_w)]/2`.
All scalars are real in the supplied polar, nonrotating sector; their
imaginary parts vanish. Generic poloidal magnetic fields can give imaginary
Maxwell components, which the implementation retains. Only Psi2 and Phi1
remain in an aligned spherical static background. Spin weights are
`2,1,0,-1,-2`, `1,0,-1`, `-1,1`. The existing harmonic phase has
`_1Y10=+sqrt(3/(8pi)) sin(theta)`, `_(-1)Y10=-_1Y10`.

**c_Weyl4_Re = +Re(Psi4)** (factor +1):

```
Re Psi4 = [C(n,v,n,v)-2C(n,v,u,v)+C(u,v,u,v)
          -C(n,w,n,w)+2C(n,w,u,w)-C(u,w,u,w)]/4.
```

The Ricci subtraction in C(k,e,k,e) is -Ric(k,k)/2 for both e=v,w,
so it cancels in the difference. Replacing C by R gives exactly the old
omlocij[1][1]-Omww; returned out.Om22 already contains the subtraction.
The identity and preservation of Psi0/Psi4 by projection are checked on
the same jets. Comparing cell-then-interpolate W is a separate spatial test.
At chi<1e-4 the inherited Riemann and tetrad floors use inconsistent metrics;
NP rejects that case instead of redefining the normalization. No tested
sphere node hits this rejection. In these geometric units G=c=1:
[Psi]=L^-2, [Phi]=L^-1, [X]=L^-3; the rotation parameter is dimensionless.

## Authorized electric/magnetic construction

Let Q_ij=3R_ij+K K_ij-K_ik K^k_j, and r_ab=8piG(T_ab-T g_ab/2).
The inherited matter-complete normal-normal Riemann is Q_ij-r_ij.
Its electric **Weyl** projection before TF is

```
Eraw_ij=R_0i0j+r_ij/2-gamma_ij(r_00/2+rscalar/6),
E=STF(Eraw)=[Q_ij-4piG S_ij]^TF,
Braw_ij=epsilon_j^kl [R_0ikl+(gamma_ik r_0l-gamma_il r_0k)/2]/2,
B=STF(Braw), STF(A)=(A+A^T)/2-gamma tr_gamma(A)/3.
```

Using the raw matter-corrected R_injn **without** its electric Weyl
projection would instead give -8piG S_ij; that is not the consult's
bracket. The isotropic terms above are retained to diagnose the removed
trace. The Ricci contribution to Braw is antisymmetric and disappears
under STF. In oriented orthonormal components reconstruct

```
C_0i0j=E_ij, C_0ijk=epsilon_jkl B_il,
C_ijkl=delta_ik E_jl-delta_il E_jk-delta_jk E_il+delta_jl E_ik,
```

with Riemann pair symmetries, then transform back using the same spatial
frame. `np_tracefree.wls` proves all 256 component symmetries, all 16
traces, both null-rotation binomial laws from tensor contractions, and
Psi0/Psi4 preservation for a generic Gauss/Codazzi tensor. It does not
assume Hamiltonian or momentum constraints vanish. Exact arithmetic is
over a polynomial ring with I and sqrt(2); the only denominators are
nonzero constants. The C++ generic E/B test independently checks traces
and the Bianchi identity in floating arithmetic.

For matter-complete data the removed electric trace equals the geometric
Hamiltonian plus K tr(A)/3 (the latter accounts for residual trace of the
conformal extrinsic curvature). The axial vector epsilon_ijk Braw_jk
is minus the geometric momentum constraint. Its antisymmetric part has
Frobenius norm |M|/sqrt(2). Braw's scalar trace
vanishes by the Codazzi Bianchi identity. These are L^-2 quantities.
The production momentum formula assumes tr(A)=0; the small tabulated
axial residual also contains derivatives of its nonzero discrete residual.

The test evaluates the existing ConstraintsCartoon equations on identical
sphere jets. Its usual Hamiltonian consumes independently evolved Gamma;
therefore a second evaluation reconstructs Gamma and its derivative from
the metric jets. Both values are reported, rather than treating a Gamma
constraint residual as roundoff. See the same-point RMS table and maximum
residuals in results-t3/summary.md and validation.json.

## Independent static oracle and control (a)

For ds²=-A(r)dt²+dr²/N(r)+r²dΩ², A=N exp(-2δ), direct four-dimensional
Christoffel/Riemann construction gives orthonormal projections

```
x=R0101=N[A''/(2A)-A'^2/(4A²)]+N'A'/(4A),
y=R0202=N A'/(2Ar), z=R1212=-N'/(2r), w=R2323=(1-N)/r²,
Psi2=(x-y+z-w)/6
 =[N''/2-Nδ''-3N'δ'/2+Nδ'^2-N'/r+Nδ'/r-(1-N)/r²]/6,
Phi1=exp(δ)V'/2=Q/[2 sqrt(8pi) F r²].
```

V'=Q exp(-δ)/(sqrt(8pi) F r²) is the native potential first integral.
The exact executable `scripts/cas/np_static_metric.wls` proves all four
projections, Ricci subtraction and the Schwarzschild/RN limits
Psi2=-M/r³ and -M/r³+Q²/r⁴. A radial boost to the maximal trumpet normal
leaves Psi2 and Phi1 invariant. All oracle points are outside the static
horizon; poles r=0,A=0,N=0,F=0,sin(theta)=0 are excluded.

`static_oracle.jl` uses the supplied Julia trumpet reader, Y,D,S coefficients
and signed C, independently of Cartesian Gauss/Codazzi. It reconstructs
r=Y/(1-s), δ=(1-s)²D, the maximal lapse, A=alpha_K²-C²/r⁴ and N=exp(2δ)A.
Nested ForwardDiff and the s-to-r chain rule supply the derivatives. No
unprovided serialized solution, root, field rescaling or gauge is selected.
RN also passes a numerical closed-form check before the oracle is written.

The static harness builds genuine t=0 hierarchies with the production
trumpet initial-data reader; its RHS throws and it calls no RH finder.
RN uses dx=.0625,.03125 at r=1,2,4. EMS uses dx=.00390625,.001953125
at the four specified cavity coordinate radii. These are t=0 diagnostics,
not low-resolution echo evolution runs. Static data enter no evolution.

All **70 scalar/sphere controls pass**. Tolerance is
`8/15 max|f_h-f_(h/2)|` from fourth-order refinement, plus arithmetic:
128ε/dx_fine² for Psi, 128ε|Phi1| for Phi, and the corresponding product
bound for X. The 1/dx² term accounts for second-difference conditioning.
The outermost EMS curvature is at its Float64 arithmetic floor and is
not claimed to show converged fourth-order truncation there. Both budget
parts, observed orders and every error are in validation.json; summary.md
lists the independent Psi2/Phi1 comparisons and vanishing-field residuals.

## Corrected rotation controls (b)

The second null vector called n in T2 is k here; n denotes the Eulerian
normal. For l fixed, m'=m+a l, k'=k+a(m+mbar)+a²l. Interchange l,k for
the other class. The exact derivatives are:

| Fixed vector | derivative X_out | derivative X_in |
|---|---|---|
| l | 2 Phi0 Psi3-6 Phi2 Psi1 | 2 Phi1 Psi0-4 Phi0 Psi1 |
| k | 2 Phi1 Psi4-4 Phi2 Psi3 | 2 Phi2 Psi1-6 Phi0 Psi3 |

`np_first_order.wls` proves all four coefficients vanish on the aligned
Coulomb background. With a counted as first order, both changes in both
classes start at second joint order. The l-fixed X_out and k-fixed X_in
coefficients themselves start at second order in non-Coulomb fields;
the opposite-class coefficients generally start at first order.

The test rotates the tetrad itself at ±h,±2h, h=1/32. Its derivative is
[8(X(h)-X(-h))-(X(2h)-X(-2h))]/(12h). CAS proves this exact for the
rotation polynomial of degree <=4. No finite-step allowance is required.
The pointwise roundoff budget is
4096ε Cmax Fmax [tmax(1+2h)²]^6/h, where tmax is the largest tetrad
absolute component sum. **No spatial tolerance is added to (b)(i).**

The T3 roundoff results for both classes and all seven spheres are in
results-t3/summary.md. No spatial tolerance is added. The new tensor is
algebraically trace-free, including on constraint-violating input data.
All 2,660 matter point/class/combination comparisons pass. Both vacuum
and matter forms pass; across both forms maximum error/budget is below
6.9e-6. The matter projected coordinate trace is below 4e-20 on the
production spheres and 2.3e-13 on
the innermost cavity sphere, whose coordinate metric components are large.

For (b)(ii), summary.md reports sphere RMS of the pointwise suppression
and perturbation ratios, including areal radii. validation.json also gives
ratios of RMS values, CAS responses, incoming ratios, and undefined-ratio
counts (none here). The protected responses are compared with both
outgoing and incoming perturbation ratios in the T3 table. X_out protects
l fixed; X_in protects k fixed. Their derivative coefficients are quadratic
in non-Coulomb fields. Protected suppression ratios range from 4.99e-7 to
7.95e-6, consistent with products of the measured perturbation ratios.
The l-fixed coefficient obeys
`|dX|/(12|Phi1 Psi2|) <= |Phi0/Phi1||Psi3/Psi2|/6 +
|Phi2/Phi1||Psi1/Psi2|/2`; k-fixed Xin reverses the coefficients.
Opposite-class ratios are O(1), since their specified denominators have
no Coulomb leading term. They are reported and are not called suppressed.

## Weyl4, projection and default controls (c,d)

Same-jet normalization remains +1. T4 replaces T3's inappropriate
pure-roundoff W comparison with interpolation orders 2,4,6 in both
routes, on both Weyl forms and r=50,75,100, l=2,4,6. Production W cell
derivatives remain fourth order throughout. Before the numerical study,
the error contract was fixed to

```
E_NP=abs(NP4-NP6), E_W=abs(W4-W6),
E=E_NP+E_W, tolerance=2 E,
accept iff abs(NP4-W4)<=tolerance.
```

T5 corrects that incomplete estimate using the already measured
fourth-versus-sixth cell-derivative control at fixed L4 interpolation:

```
Eprime=E+abs(W4cell-W6cell),
accept iff abs(NP4-W4cell)<=2 Eprime.
```

All 18 comparisons pass within Eprime itself. The earlier r=75,l=6
failure now has discrepancy 2.02193e-10 against Eprime=3.23412e-10
(2Eprime=6.46824e-10). The complete table, including abs(NP-W6cell),
is in results-t5/cell-error.md and cell-error.json. NP is closer to the
sixth-order cell result in only 4/18 modes; discrepancy ratios range
from .260585 to 5.86966. The reviewer's claimed general strong improvement
is not observed and is not used as an acceptance condition.

Each estimate term compares two orders of the same route; their sum
accounts for both discretizations. Factor 2 is the permitted maximum.
This is an empirical order-variation estimate at this checkpoint, not a
rigorous bound on continuum error. No L2-based larger allowance replaces
it after seeing failures. Relative differences use abs(NP4), with the
areal R factor included in every mode exactly as written by each route.

All L2-to-L4 cross-route differences decrease; all L6 cross-route
differences exceed L4. Interpolation-order variation leaves the fixed
fourth-order cell-Weyl derivative error in place. The two l=6,r=75
comparisons failed under T4's incomplete budget: abs(NP4-W4)=2.02193e-10, E_NP=6.3051e-11,
E_W=5.478e-12, tolerance=1.37060e-10. The relative difference is
3.92304e-4, against a relative budget 2.65928e-4. A factor about 2.95
on E would be required, above the allowed 2. All 18 comparisons and
their raw modes are in results-t4/summary.md and validation.json. T5's
corrected budget above supersedes that readiness failure.

The checkpoint's point-transfer implementation requires three ghosts.
An attempted four-ghost restart was rejected before extraction. L6's
eight-point second-derivative stencils therefore use the existing
Lagrange search with its source bounded to the sampled FArrayBox. This
shifts a full-width stencil near box edges, retaining the original three
ghosts, transfer law and valid-cell hierarchy. Every read is bounds
checked. The unbounded L6 negative control rejects an out-of-array
stencil; it is not used as numerical evidence. L4 agrees with the
ordinary uninstrumented interpolator and the accepted T3 CSVs within
4096 epsilon times each column's scale.

The uninstrumented production L4 matter NP call costs a median 6.323 ms
(6.258–6.412 ms over three calls after one warm-up) on this checkpoint,
for the three spheres and 97 theta points. This includes construction,
interpolation, reconstruction, projection and CSV IO; it excludes restart
loading and shared ghost/cell-Weyl preparation. The existing radiation
call is 3.290 ms. These are local measurements, not an HPC wall-time
prediction. T5's corrected budget passes the two earlier failed modes.

An independent associated-Legendre factorial implementation rebuilds
all **980 mode components**, including signs, spin weights and R factors;
maximum residual is .00181 of its 2048ε absolute-sum budget. Both X
expressions are independently checked at every point before projection.

T3's EMSBH2DLevel.cpp syntax check used the prescribed temporary header
workaround. Default-off and NP-enabled frozen runs reproduce the 42 old radiation
rows to 1e-12 relative (+1e-100 at zeros), maximum difference zero. No new
NP file appears with the flag off. Flat space, complex Maxwell signs and
normalization, chi-floor rejection and activation-dependency controls pass.
T5b additionally builds, links and runs the real production executable.
Only the two checkpoint-header methods in GRAMRLevel are now protected,
allowing the EMS overrides to call their base methods on local g++-16.
Other private methods keep their original access.

## Use in a production run

The validated scope is the existing axisymmetric extraction, m=0 and
l<=6. Keep the existing sphere, center, angular quadrature and extraction
schedule parameters, and enable:

```
activate_extraction = 1
ems_radiation_activate = true
ems_np_activate = true
ems_weyl_matter = true
```

The two new/related switches default false. NP itself is always
matter-complete; ems_weyl_matter enables the matching legacy W control.
The columns and R normalization are specified above. A puncture crossing
produces unavailable NaN rows for that sphere/time and one warning;
consumers must skip these rows. On the supplied t=602 checkpoint, the
incremental NP extraction costs median 6.323 ms for three spheres and
97 theta points, including CSV IO (measurement boundaries given above).

The linked smoke uses the authenticated **non-echo** reference
alpha=20, Q/M=.7 trumpet, at dx=.25 with KO sigma=1, three RK4 steps to
t=.075, spheres r=2,4 and 33 theta points. This is a build control, not
an evolution-accuracy or stability result. Raw minimum chi is negative
on this coarse grid in both baseline and current builds; no new branch
or output difference is observed. The baseline is an authenticated
2516-blob checkout of 00d7171, with only the same two access fixes needed
to compile on this machine. Both builds link the real local Chombo/HDF5.
All three runs exit 0. Off/off comparisons pass at absolute tolerance
1e-12 for both checkpoints and nine legacy text outputs (maximum text
difference observed zero). This covers Float64 conditioning at this dx
over three RK4 steps; byte identity is not required. NP is absent off;
on produces 56 finite rows at four times. Evidence is in
results-t5/production-verification.json and the three run directories.
The final current build took 22.375 s (peak RSS .572 GB), with smoke runs
below .3 s and .019 GB. No run approached the five-minute dry-run limit.

Reproduce with `python3 Tests/EMSNP/production_smoke.py --prepare`, then
`python3 Tests/EMSNP/production_smoke.py`, then its `--verify` mode.
The script uses an external baseline checkout and archives checkpoints
under EMS_NP_LOG_ROOT/t5/smoke-checkpoints with SHA-256 receipts.

## Evidence and reproduction

All five inputs were rehashed against T1 identities before and after T4;
see results-t4/inputs.json and inputs-final.json.
The t=602 checkpoint has 16 synchronized levels and 28 evolved fields.
The full hierarchy is read without advancing, regridding or replacing data.
Parameter copies state every restart/sphere/test selector.

Run `zsh Tests/EMSNP/run_t4.sh` from the worktree. It builds with the local
serial Chombo, syntax-checks the production level with the prescribed
temporary header workaround, checks the negative stencil control, then
runs the three orders and the uninstrumented L4 cost measurement on one
unchanged hierarchy. It authenticates the inputs again and reports every
comparison. Exit 1 records the failed readiness gate. See progress.md and
results-t4/summary.md, validation.json and full logs. Accepted T3 controls
and their executable derivations remain in results-t3 and run_t3.sh.
No Julia process, evolution, HPC, SSH, transfer or push was launched in T4.
