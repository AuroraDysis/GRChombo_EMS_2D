# One-tetrad EMS Newman–Penrose extraction

NPTET T3, `ariadne-run:v1:88362_7ba22e11`, resumes bd8e890 on `np-scalars`.
Restart continuation: `ariadne-run:v1:11338_1612a245`.
T3 authorizes the consult's symmetric, trace-free electric/magnetic Weyl
construction. Historical records remain in [README-T1.md](README-T1.md)
and [README-T2.md](README-T2.md).

**Default-off implementation; not ready under all requested gates.**
The authorized trace-free construction fixes the rotation-law failure:
28/28 matter rotation groups pass, as do the 70 static checks and default
output checks. The requested pure roundoff comparison with legacy
production W passes all odd modes, but fails 18 even modes because of
the distinct evaluation orders described below.

## Production interface

Enable `ems_radiation_activate = true` and `ems_np_activate = true` to write
`ems_np.csv` on the existing radii, center, Simpson quadrature and schedule,
including the existing t=0 call. NP defaults false and requires EMS radiation.
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

Same-jet normalization remains +1 in both forms. T3's production CSV
comparison uses exactly `C epsilon Ntheta max_l|mode|`, epsilon=2.2e-16,
Ntheta=97 and C=4096. C is a conservative 16 times the 256-term tensor
contraction count, followed by the angular sum. It is an arithmetic
allowance, with **no spatial error term**. All modes, including the
symmetry-suppressed l=3 values, are reported against the sphere scale.
No relative error with respect to that tiny l=3 mode is used as a gate.
All six l=3 comparisons pass: at r=100 the matter NP and W values are
7.0429552e-12 and 7.0095248e-12, error 3.343e-14 against the sphere-scale
tolerance 3.765e-14. These odd modes are small but not zero to machine
roundoff in the checkpoint; no value was replaced by its symmetry limit.

The legacy W first reconstructs curvature on cells then interpolates it;
NP interpolates metric jets then reconstructs curvature. Those operations
do not commute at finite resolution. Projection preserves Psi4 on identical
jets, and cannot remove that pre-existing spatial difference. The strict
T3 gate is retained even where it consequently fails. The old sixth-order
cell control remains available as historical discretization evidence; it
is not an allowance in this tranche's acceptance.
The 18 l=2,4,6 comparisons fail the strict arithmetic gate. For the matter
l=2 mode, errors at r=50,75,100 are 4.094e-10,2.804e-10,2.943e-10,
against tolerances 2.515e-14,3.256e-14,3.765e-14. No C was enlarged or
spatial allowance added to turn those failures into passes.

An independent associated-Legendre factorial implementation rebuilds
all **980 mode components**, including signs, spin weights and R factors;
maximum residual is .00181 of its 2048ε absolute-sum budget. Both X
expressions are independently checked at every point before projection.

EMSBH2DLevel.cpp passes -fsyntax-only with the prescribed temporary header
copy and line-51 private-to-protected workaround. Source headers are not
patched. Default-off and NP-enabled runs reproduce the 42 old radiation
rows to 1e-12 relative (+1e-100 at zeros), maximum difference zero. No new
NP file appears with the flag off. Flat space, complex Maxwell signs and
normalization, chi-floor rejection and activation-dependency controls pass.
A full production executable link/evolution is not claimed: its known
private-header build issue is outside this task.

## Evidence and reproduction

All five inputs were rehashed against T1 identities; see results-t3/inputs.json.
The t=602 checkpoint has 16 synchronized levels and 28 evolved fields.
The full hierarchy is read without advancing, regridding or replacing data.
Parameter copies state every restart/sphere/test selector.

Run `zsh Tests/EMSNP/run_t3.sh`. It builds with the local serial Chombo,
runs syntax/unit/CAS/oracle controls, both static resolutions, both frozen
Weyl forms on production and cavity spheres, and the default replay.
HDF5 restart access is sequential because the reader takes an exclusive
file lock. All reports are written; final exit is **1** when
`ready_for_simulation=false`. This failure is intentional evidence, not
silently accepted as a passing test suite.

See results-t3/summary.md, validation.json and progress.md. Full logs,
point files, NP/radiation CSVs, parameters, compiler commands, CAS scripts
and transcripts are retained. No program was expected above five minutes;
the complete workflow measured about one minute. No >30-minute program,
evolution, HPC, SSH, file transfer or push was launched.
