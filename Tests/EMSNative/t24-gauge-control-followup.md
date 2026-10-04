# T24 follow-up: gauge stability, numerical horizons and safe abort

**READY-EXCEPT.** The collected exp-0026 interventions and the code-stencil
audit support the shift/Gamma discrete-stability candidate under the registered
reading. They do not uniquely prove its mechanism or admit a production restart.
The expanded RK4 rule, all local bounded numerical horizon probes and the
serial-qualified opt-in safe abort are complete. No strict individual or common
horizon is demonstrated. Original common-finder iteration logs are absent from
the collected packet; MPI abort execution is blocked at local MPI_Init.
No full evolution replay is repeated. No job remains pending. No gauge, KO,
transfer, initial-data or evolution-equation remedy is applied; no commit,
SSH or cluster operation is performed.

## 1. Code-level frozen-coefficient audit

All equation references here are to archived **d507438**, not the later gauge
control. `ExperimentalGauge.hpp:77–90` supplies beta_t=advec(beta)+c_Gamma Gamma
−eta beta−B and B_t=−0.1 B. Its lapse coefficient 1.8 and B decay 0.1 are
existing literals; this tranche changes neither. `CCZ4Cartoon.impl.hpp:369–408`
supplies the pure inverse-metric shift Laplacian and **1/3** longitudinal second
derivative, the cartoon reg03/reg04 contributions and the evolved Gamma row.
With the production CCZ4 formulation and kappa3=1, contracted Christoffel plus
2 Z/chi is **evolved Gamma** exactly. No contracted-Christoffel replacement is
silently made. Matter, curvature, A and damping contributions proportional to
lapse vanish in the specified alpha→0 limit. This does not remove the finite
alpha/chi terms from the actual evolution, especially once chi hits a floor.

Write H=inverse([[h11,h12],[h12,h22]]), Hww=1/hww, rho the cartoon radius,
d=partial_i beta^i+beta_rho/rho, and dw=partial_i beta^i−2 beta_rho/rho.
    The alpha=0 subsystem implemented by the audit is

```
beta_t = adv(beta) + c_Gamma Gamma - eta beta - B + KO(beta)
B_t = -0.1 B + KO(B)
Gamma_i,t = adv(Gamma_i) + (2/3) d Gamma_i - Gamma_j partial_j beta_i
            + H_jk partial_j partial_k beta_i
            + (1/3) H_ij partial_j partial_k beta_k
            + Hww[(partial_rho beta_i)/rho - delta_i,rho beta_rho/rho^2]
            + (1/3) H_ij[(partial_j beta_rho)/rho - delta_j,rho beta_rho/rho^2]
            + KO(Gamma_i)
h_ij,t = adv(h_ij) -(2/3) h_ij d
          + h_ki partial_j beta_k + h_kj partial_i beta_k + KO(h_ij)
hww_t = adv(hww) -(2/3) hww dw + KO(hww)
chi_t = adv(chi) -(2/3) chi d + KO(chi)
```

The script differentiates these expressions analytically around a homogeneous
Cartesian background, retains cylindrical source terms and delta(H)=−H delta(h) H,
and assembles an **11×11** complex matrix. Gamma and B perturbations are scaled
by dx, chi by its numerical background. The damping rates are physical rates:
eta enters as eta*dx, not eta=1 in grid units. Actual numerical Gamma and shift
are supplied. K/Theta/A have advective alpha-independent blocks; at alpha=0
they do not close an extra fast feedback loop into this subsystem. A lapse
perturbation can feed those rows, but its coupling back through alpha K vanishes
at alpha=0; this audit freezes the collapsed lapse rather than claiming a
complete finite-lapse CCZ4 spectrum.

`FourthOrderDerivatives.hpp:114–190,258–365` gives the pure second stencil
(-1,16,-30,16,-1)/(12 dx^2); mixed derivatives are the tensor product of
centered fourth-order first derivatives. The shift-dependent one-sided
advection stencil is (-3,-10,18,-6,1)/(12 dx) at offsets (-1,0,1,2,3) for
positive shift. Its Nyquist symbol is **−8/(3 dx)**, not a centered imaginary
symbol. KO is (1,-6,15,-20,15,-6,1)*sigma/(64 dx), **on every evolved row,
including B** (`CCZ4Cartoon.impl.hpp:162`). Each direction contributes
−sigma/dx at Nyquist. The compiled Chombo `LevelRK4.H:107–162` has the classical
RK4 weights and stage times; its stability polynomial is
R(z)=1+z+z^2/2+z^3/6+z^4/24.

The exact CAS witnesses in [t24-gauge-stability.json](t24-gauge-stability.json)
verify the stencil symbols, two-component beta/Gamma characteristic polynomial,
and |R(−1+i y)|^2−1. A separately evaluated 70-digit **four-stage recurrence**
checks the boundary and both sides. The reduced zero-shift, diagonal-metric,
two-direction Nyquist criterion is indeed
**S <= 4.9062387672032055**, where S=H22+0.75 H11. Without KO it is 6;
the single-rho-direction KO limit on H22 is **6.36217048021**.
These are reduced-model limits, not a complete AMR or variable-background proof.

For nonzero h12 the fastest Nyquist coefficient is
c_Gamma*(trace(H)+lambda_max(H)/3); H22+0.75 H11 is only its diagonal rho-mode
special case. Shift upwinding and eta move the RK4 argument left. At the saved
cell the Nyquist wave has Re(z)≈−1.03274, changing the bound. The resolved
129×129 wave-number scan brackets the decaying-branch transition, when all four
h components are scaled together, at **4.774408817 < S_limit < 4.774408913**.
This is a sampled numerical bracket, not a rigorous bound between wave numbers.
Thus **4.906 is the correct scalar reduction but is too permissive as this
cell's complete discrete criterion**. The saved S is about
97.658% of that sampled transition.

The independent blockwise L12 active-valid maximum is S=4.662576980611591 at
(5243409,135), rho=0.057891845703125, dx=0.00042724609375. Its h11/h12/h22/hww
are 7.75956300347 / 0.0422448889773 / 0.219248673040 / 0.588659453868;
beta=(0.0221622779682,0.00231354007068), chi=0.000455879319582,
lapse=1e-12. All 28 values are in [t24-gauge-cell.json](t24-gauge-cell.json).
The prescribed variants use the measured cell, scaling all h entries by
S_original/S_target; this choice preserves metric shape but not unit determinant.
They are operator sensitivity tests, not physically initialized states.

| Case | metric U | weighted S | max |R|, full | max |R|, Re(lambda) <= 0 | limiting decaying k/pi | Nyquist wave |R| |
| --- | --- | --- | --- | --- | --- | --- |
| step152 | 4.66257698 | 4.66257698 | 1.001982115 | 0.999999886 | (-0.1015625, -0.0859375) | 0.969970252 |
| S_4.9 | 4.90000000 | 4.90000000 | 1.035229566 | 1.035229566 | (-0.984375, -1) | 1.034498487 |
| S_4.96 | 4.96000000 | 4.96000000 | 1.052878456 | 1.052878456 | (-0.984375, 1) | 1.052098413 |
| S_5.03 | 5.03000000 | 5.03000000 | 1.074156157 | 1.074156157 | (-0.984375, 1) | 1.073316283 |
| S_5.5 | 5.50000000 | 5.50000000 | 1.237083089 | 1.237083089 | (-0.984375, 1) | 1.235772521 |
| S_6.4 | 6.40000000 | 6.40000000 | 1.652661490 | 1.652661490 | (0.9765625, -1) | 1.649994631 |
| half_dt | 4.66257698 | 4.66257698 | 1.000990748 | 0.999999943 | (-0.1015625, -0.0859375) | 0.623752049 |
| sigma_half | 4.66257698 | 4.66257698 | 1.001982215 | 0.999999369 | (-0.1171875, -0.0859375) | 0.474698120 |
| sigma_zero | 4.66257698 | 4.66257698 | 1.001982315 | 0.999999901 | (-0.171875, -0.1640625) | 0.590705757 |
| cGamma_half | 4.66257698 | 2.33128849 | 1.001370814 | 0.999999886 | (-0.1015625, -0.0859375) | 0.588826064 |
| half_dt_S_25.0 | 25.00000000 | 25.00000000 | 1.002363557 | 0.999999943 | (-0.1015625, -0.0859375) | 0.946863726 |
| half_dt_S_26.0 | 26.00000000 | 26.00000000 | 1.100245345 | 1.100245345 | (0.9609375, 1) | 1.097274428 |


The unstable high-frequency branches at S=4.9/5.5/6.4 are dominated by **beta2
and Gamma2**, with ky at Nyquist and kx near Nyquist. The exact modes and complex
RK arguments are in [t24-gauge-stability.csv](t24-gauge-stability.csv).
The neutral near-zero-wave branches mean a stable principal scan's overall
maximum is 1; the Nyquist column retains the useful margin. Half dt, half
c_Gamma, sigma=0.5 and sigma=0 have no decaying-branch RK4 violation at the
saved metric. Sigma=0 is **an offline operator test only**; KO remains 1 in
every production/evolution test and no KO removal is proposed.

The **full** frozen cartoon matrix has positive-real low-frequency eigenvalues:
its overall maximum is 1.001982115 at k/pi≈(−0.0078125,−0.0546875) for the saved
cell. These are already growing eigenvalues of that homogeneous cylindrical
semi-discrete operator, not an RK4 amplification error of a decaying mode.
Freezing Cartesian background derivatives at zero loses the actual numerical
gradients and radial wave envelopes. Their interpretation as radial amplitude
transport, physical/local source growth or a real instability remains open.
S alone cannot bound them, and they are not silently removed from the table.

### RK4 rule, expanded before use in a run design

For arbitrary positive c_Gamma define **U=H22+0.75 H11** and the weighted
rho-mode coefficient **S=c_Gamma*(H11+4 H22/3)=(c_Gamma/0.75) U**.
Earlier columns named S in the single-cell script mean U; the CSV now also
records explicit `metric_U` and `weighted_S` columns. With diagonal H and
H22 >= H11, zero shift, alpha=0 and negligible eta*dx, the two-direction
Nyquist branch has z=−2 sigma nu ± i nu sqrt(16 S/3), nu=dt/dx.
Let y_RK(a) be the first positive edge of the RK4 stable segment on Re(z)=a.
Then **S_lim=3 y_RK(−2 sigma nu)^2/(16 nu^2)**. The weighted limit is
independent of c_Gamma; the permitted metric U is 0.75 S_lim/c_Gamma.
For general off-diagonal H use S_eff=c_Gamma*(trace(H)+lambda_max(H)/3)
at Nyquist, and scan all wave numbers with the actual shift and metric.
The restriction on H's ordering is part of the rho-mode claim, not an implicit
assumption about every cell in a simulation.

| dt/dx | sigma | weighted S limit | metric U limit, c=.75 | U limit, c=.375 | U limit, c=1 | sampled S, measured cell/c=.75 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 0.0 | 6.00000000 | 6.00000000 | 12.00000000 | 4.50000000 | 5.99326177 |
| 0.5 | 0.3 | 6.46654932 | 6.46654932 | 12.93309864 | 4.84991199 | 6.36708118 |
| 0.5 | 0.5 | 6.36217048 | 6.36217048 | 12.72434096 | 4.77162786 | 6.30249035 |
| 0.5 | 1.0 | 4.90623877 | 4.90623877 | 9.81247753 | 3.67967908 | 4.77440884 |
| 0.375 | 0.0 | 10.66666667 | 10.66666667 | 21.33333333 | 8.00000000 | 10.67469706 |
| 0.375 | 0.3 | 11.42637072 | 11.42637072 | 22.85274144 | 8.56977804 | 11.24459926 |
| 0.375 | 0.5 | 11.48817724 | 11.48817724 | 22.97635449 | 8.61613293 | 11.42483155 |
| 0.375 | 1.0 | 10.35318118 | 10.35318118 | 20.70636236 | 7.76488588 | 10.21222849 |
| 0.25 | 0.0 | 24.00000000 | 24.00000000 | 48.00000000 | 18.00000000 | 24.05277652 |
| 0.25 | 0.3 | 25.36446906 | 25.36446906 | 50.72893812 | 19.02335180 | 25.04807384 |
| 0.25 | 0.5 | 25.78169746 | 25.78169746 | 51.56339492 | 19.33627310 | 25.46371285 |
| 0.25 | 1.0 | 25.44868192 | 25.44868192 | 50.89736384 | 19.08651144 | 25.34582197 |
| 0.125 | 0.0 | 96.00000000 | 96.00000000 | 192.00000000 | 72.00000000 | 96.24382195 |
| 0.125 | 0.3 | 99.23373627 | 99.23373627 | 198.46747254 | 74.42530220 | 98.77535281 |
| 0.125 | 0.5 | 100.81653105 | 100.81653105 | 201.63306210 | 75.61239829 | 99.95365413 |
| 0.125 | 1.0 | 103.12678984 | 103.12678984 | 206.25357968 | 77.34509238 | 102.14194570 |


[t24-stability-rule.py](t24-stability-rule.py) supplies all **48** requested
(nu,sigma,c_Gamma) combinations in [t24-stability-limits.csv](t24-stability-limits.csv).
For each real RK argument it closes the exact conjugate-polynomial witness
and independently checks the positive root with a 70-digit four-stage recurrence.
The last column is a 257×257 **decaying principal-branch** scan at the measured
L12 metric shape, shift, Gamma and dx, scaling the inverse metric to each S.
It is a reproducible sampled threshold, not a universal table depending only
on three parameters. `t24-stability-rule.json` records the narrow numerical
brackets, modes, witness statuses and resource receipt. Shift, h12, radial
source terms and background gradients prevent a complete universal criterion
in just (nu,sigma,c_Gamma). A safety factor is not invented after these data.
Sigma variants are offline calculations only; evolution retains sigma=1.

At nu=0.25 the measured-cell h/chi advection and B-decay branches have max
RK amplification <=1 across the full sampled domain. The full 11-row cartoon
scan at weighted S=25 has no decaying-mode RK violation; at S=26 its beta2/
Gamma2 wave is unstable (max |R|=1.10025). No other **decaying** alpha-independent
mode crosses the RK boundary first in this audit. The positive-real,
low-frequency homogeneous cylindrical modes reported above remain present
even at smaller dt: they are source/transport growth in this frozen model,
not a competing numerical CFL threshold. Their interpretation in the actual
inhomogeneous radial problem is **not closed**. K/Theta/A and lapse perturbations
do not add a fast feedback loop at fixed alpha=0: their advected/zero-order
blocks are triangular with respect to this principal subsystem. In the actually
compiled CCZ4 file, `COVARIANTZ4` is defined and the damping uses
kappa1*alpha/(0.005+alpha), which also vanishes as alpha approaches zero
(`CCZ4Cartoon.impl.hpp:13,346–351`). No finite-alpha, matter, floor-projection,
AMR-interface or boundary stability guarantee follows from this limit.

### Growth and the exp-0026 interventions

At nu=.5, sigma=1, **S=4.96 / 5.03** gives scalar Nyquist gains
**1.015262958 / 1.035814535 per fine step**, e-folds **66.01684 / 28.41869
steps**. On L10/L11/L12 these are e-fold times
0.056411/0.028205/0.014103 and 0.024284/0.012142/0.006071 M_i.
The measured-shift full cartoon variants instead give Nyquist gains
**1.052098413 / 1.073316283** (whole-grid maxima 1.052878456/1.074156157).
They are a sensitivity scan with the step-152 background held fixed, not
measured growth of a trajectory. Both sets are exposed rather than using
the scalar estimate as the full code's amplification.

The collected S153/S154 slope gives the scalar crossing **t=134.408499455**
and extrapolated S=5.031717543 at 135.173828125.
The resulting quasi-frozen Nyquist ramp predicts log10 amplitude gains
**6.903/13.806/27.611** on L10/L11/L12 over the **0.765329 M_i** interval.
This accounts for how a tiny grid-scale seed can become very large after
hundreds to thousands of fine steps; it is consistent with a delayed NaN.
It is not a fit to measured checkerboard amplitudes: the cell of the spatial
maximum can migrate, the coefficients change, and the scalar model omits
shift and AMR effects. See [t24-stability-growth.csv](t24-stability-growth.csv).

| Leg/step | t | max weighted S | cells >4.9 | native evidence |
| --- | --- | --- | --- | --- |
| exp-0024/152 | 133.0 | 4.677619954 | 0 | reduction only |
| P/153 | 133.875 | 4.818769342 | 0 | reduction only |
| P/154 | 134.75 | 4.962229183 | 1548 | reduction only |
| G/158 | 138.25 | 2.658945805 | 0 | 224/224 zero |
| D/158 | 136.5 | 5.244598436 | 5444 | 224/224 zero |


[t24-cluster-evidence.py](t24-cluster-evidence.py) checked **147** locally
collected files against their SHA256 manifests. It then checked each G/D/D0
**native** receipt: 224/224 rank exits zero, the independent finite-state
reduction has no unreadable checkpoints or valid non-finite values, and the
dt readback and absent-static-input proofs hold. The own reduction receipt
also reports all active metric inversions valid. The summary does not mistake
the scheduler exit for a numerical pass. Raw intervention checkpoints remain
on the cluster; no bulk file was copied into this worktree.

Both one-parameter interventions pass coarse 155 and the original event.
G lowers c_Gamma to .375, giving weighted S=2.658946 at t=138.25;
D keeps c_Gamma=.75 and halves dt, giving S=5.244598 at t=136.5 against the
scalar limit 25.448682. Restart already supplies the requested per-level dt.
The checkpoint-state gauge kick in G leaves B untouched, so G and D are
different gauge trajectories and do not certify each other's physical accuracy.
**Supported:** the collapsed-lapse beta/Gamma loop approaches/crosses its
discrete RK boundary and two targeted interventions prevent the original NaN.
**Not established:** a unique causal mechanism, the exact first unstable
stage, measured checkerboard growth, convergence or production admission.
The appropriate run-design rule evaluates the complete discrete symbol with
the current numerical metric/shift and leaves an explicit margin; S<4.906
alone is too permissive at the measured shift. The existing half-dt control
has a large measured margin through its stated bound, not through merger.

## 2. Numerical apparent-horizon searches

[t24-horizon-search.py](t24-horizon-search.py) loaded only the numerical step-152
checkpoint into the frozen T17 harness, with unchanged RHUnion/RHSurf and native
Lagrange<4> interpolation. Static and CTT paths were deliberately absent; no
advance or initialization was executed. Seeds 0.002,0.004,0.008,0.016,0.032,0.064,
0.096 M_i were chosen as a numerical radius sweep, independently at both
numerically tracked centres, not from a static trumpet. N48, chase 1, quota 1,
64 updates and the inherited stage-0 squared-expansion threshold 1e-7 were used.

| Hole | Seed radius | Squared expansion after 64 updates | Trial rho_max | Status |
| --- | --- | --- | --- | --- |
| 1 | 0.002 | 9.6364445 | 0.0019957945 | UPDATE_CAP |
| 1 | 0.004 | 2.9550364 | 0.0039942718 | UPDATE_CAP |
| 1 | 0.008 | 1.3454555 | 0.0080022756 | UPDATE_CAP |
| 1 | 0.016 | 0.91395261 | 0.016044113 | UPDATE_CAP |
| 1 | 0.032 | 0.58571499 | 0.032204374 | UPDATE_CAP |
| 1 | 0.064 | 0.37034439 | 0.064816098 | UPDATE_CAP |
| 1 | 0.096 | 0.29201689 | 0.097833616 | UPDATE_CAP |
| 2 | 0.002 | 9.6364445 | 0.0019957945 | UPDATE_CAP |
| 2 | 0.004 | 2.9550364 | 0.0039942718 | UPDATE_CAP |
| 2 | 0.008 | 1.3454555 | 0.0080022756 | UPDATE_CAP |
| 2 | 0.016 | 0.91395261 | 0.016044113 | UPDATE_CAP |
| 2 | 0.032 | 0.58571499 | 0.032204374 | UPDATE_CAP |
| 2 | 0.064 | 0.37034439 | 0.064816098 | UPDATE_CAP |
| 2 | 0.096 | 0.29201689 | 0.097833616 | UPDATE_CAP |


None qualifies even the first stage; the minimum squared residual is
**0.2920168877**, and every search reports UPDATE_CAP. The trial extents in the
table are **unconverged surfaces**, never apparent-horizon extents. N96 and
later strict stages were consequently not run. The native nonzero exit **1**
is the bounded finder failure, not a launcher success, memory failure or NaN.
Cost including checkpoint load: **34.418 s**,
peak RSS **2,914,631,680 B**; actual chase iteration
time is about 18.675 s across all 14 surfaces. See
[t24-horizon-search.csv](t24-horizon-search.csv) and its own JSON receipt.

### Expanded common search and the production flow

The same unchanged author's finder searched fresh spheres at numerical midpoint
x=2240 with coordinate radii **.3,.5,.75,1,1.5,2,3,4 M_i**. N48/chase1/quota1
was bounded by 128 updates; N96/chase1/quota1 by 90 s and 2048 updates.
An additional N48/chase.125/quota1 control tried .75 and 3 for 90 s.
All queries used the saved t=133 state and native interpolation, never a static
seed, target or subtraction. All retain the same first-stage squared threshold
1e-7. The N96 run ended after 587 updates; the slow N48 run after 574.

| Probe | Seed | Squared expansion | A (trial) | Q (trial) | rho_max (trial) | stop |
| --- | --- | --- | --- | --- | --- | --- |
| common N48 chase1 | 0.3 | 0.16018541 | 395.54174 | 2.0303459 | 0.28128913 | UPDATE_CAP |
| common N48 chase1 | 0.5 | 0.18537023 | 637.6207 | 2.0499265 | 0.61767767 | UPDATE_CAP |
| common N48 chase1 | 0.75 | 0.1483043 | 770.64635 | 2.0653494 | 1.0686733 | UPDATE_CAP |
| common N48 chase1 | 1.0 | 0.18804543 | 755.43279 | 2.0770537 | 1.5940408 | UPDATE_CAP |
| common N48 chase1 | 1.5 | 3.2189923 | 632.78339 | 2.0876768 | 4.0779226 | UPDATE_CAP |
| common N48 chase1 | 2.0 | 10.832789 | 678.52462 | 2.0904544 | 6.9299347 | UPDATE_CAP |
| common N48 chase1 | 3.0 | 3100.295 | 857.31553 | 2.0989033 | 9.397069 | UPDATE_CAP |
| common N48 chase1 | 4.0 | 2471.7531 | 1105.6438 | 2.8243028 | 9.9946459 | UPDATE_CAP |
| common N96 chase1 | 0.3 | 0.16690857 | 400.19329 | 2.0300099 | 0.27867363 | TIME_CAP |
| common N96 chase1 | 0.5 | 0.18228942 | 648.4717 | 2.0504248 | 0.64721531 | TIME_CAP |
| common N96 chase1 | 0.75 | 0.14441056 | 778.57754 | 2.0672071 | 1.1272081 | TIME_CAP |
| common N96 chase1 | 1.0 | 0.201245 | 738.08684 | 2.0792396 | 1.7960423 | TIME_CAP |
| common N96 chase1 | 1.5 | 35.593369 | 633.87897 | 2.0891721 | 4.9368597 | TIME_CAP |
| common N96 chase1 | 2.0 | 2394.2974 | 742.96874 | 2.082263 | 7.8781241 | TIME_CAP |
| common N96 chase1 | 3.0 | 3528.2085 | 1047.1198 | 2.5058677 | 9.7446837 | TIME_CAP |
| common N96 chase1 | 4.0 | 3082.2201 | 1101.6074 | 2.4769751 | 9.9986614 | TIME_CAP |
| common N48 chase.125 | 0.75 | 0.16351806 | 748.44706 | 2.0619649 | 0.90501074 | TIME_CAP |
| common N48 chase.125 | 3.0 | 0.050858645 | 675.88151 | 2.0910834 | 7.2195828 | TIME_CAP |


| Probe | Native exit | Wall s | RSS GB | Best final squared residual |
| --- | --- | --- | --- | --- |
| individual | 1 | 34.418 | 2.915 | 0.29201689 |
| common N48 chase1 | 1 | 38.105 | 3.099 | 0.1483043 |
| common N96 chase1 | 1 | 105.707 | 3.038 | 0.14441056 |
| common N48 chase.125 | 1 | 104.542 | 2.925 | 0.050858645 |


Every native probe exits **1**, with no resource gate: these are bounded
finder failures, not numerical evolution failures. None qualifies even stage
one, so later 1e-10/1e-12 qualification and N48/N96 A/Q agreement cannot be
claimed. **No qualified individual or common horizon, area or charge is
established at t=133.** Trial areas/charges/extents in the table describe
unconverged surfaces. Thus whether rho=.047–.08 above either puncture is
inside an apparent horizon is **undetermined**. Failure of these finite
searches does not prove the absence of a trapped surface or of a non-star-shaped
surface. All trials are reproduced by [t24-horizon-search.py](t24-horizon-search.py),
[t24-horizon-all-probes.csv](t24-horizon-all-probes.csv) and their own receipts
in [t24-horizon-summary.json](t24-horizon-summary.json).

The collected run summary, [t24-common-horizon-history.csv](t24-common-horizon-history.csv),
has **18** common rows from steps 96–128, all TIME_CAP/UNRESOLVED: areas
1016.6–1279.9, squared residuals .0100433–563.326, and about 234.8 s/call.
The cap is the wrapper's `offline_seconds`; the author's stopping test is
area-average Theta_plus^2 (`RHSurf.hpp:228`) against the specified threshold,
not stabilization of area. Its update step is
delta_f=−.25*chase*dtheta^2*f^2*(Theta_plus+.01*Theta_plus/sqrt(Theta_plus^2+1e-13)*sqrt(abs(Theta_plus)))
(`RHSurf.hpp:823–839`). FAR uses one fresh plus **five stale** steps before
another interpolation (`RHUnion.hpp:297–313`). Per-node radii are clamped
to [.0001,10]; a mean radius outside those bounds resets to 5
(`RHUnion.hpp:220–229`). Newton polishing is commented out, so changing its
parameter cannot rescue this build. The residual/time tests, step and bounds
are the unchanged author's machinery, not a new horizon solver.

The actual driver initially seeds common surfaces at d/2+.5 and d+1 and
midpoint. At later calls it copies the preceding `common.dat` into **both**
seeds (`ev1/.../runtime/production.py:231–236`). It saves the lowest-error
finite trial as `common.dat` even when UNRESOLVED (`:259–265`); thus subsequent
paired searches are not independent fresh-radius trials. Already at t=133,
fresh radius .3 has numerical area **370.837**, and fresh radii 1–4 have
areas about 777–828. A large area is neither an admission nor a measure of
proximity to a horizon on this numerical geometry.

The controlled flow-step comparison is informative: for seed 3 at accumulated
chase 72, N48/chase1 has squared residual **2.37627**; at accumulated chase
71.75, N48/chase.125 gives **.0508586** (about 46.7 times smaller).
Thus large angular-flow excursions are sensitive to the step, beyond merely
allowing less chase time. Yet the smaller step remains more than five orders
above the strict threshold and does not demonstrate a horizon. The numerical
pilot with larger seeds develops very distorted and near-bound trial radii.
This supports investigating bounded, independently reseeded flow controls;
it does not prove the area-1000 history was caused by a reset or a particular
seed. The collected **ev1 packet lacks `diagnostics/` and `shapes/`**, including
the original finder run/progress/rh_f logs. Its available rank logs are the
evolution logs (`RH_activate=false`). Therefore the exact original reset,
stale-step excursion or seed basin cannot be read from the supplied finder
log: that requested evidence is missing. The report separates these measured
flow and policy facts from an unproved diagnosis of each historical call.

## 3. Implemented opt-in joined NaN abort

`ems_safe_nan_abort=false` preserves the old NanCheck path.
`ems_safe_nan_abort=true` uses [SafeNanAbort.hpp](../../Source/BoxUtils/SafeNanAbort.hpp)
in `EMSBH2DLevel::specificAdvance`, after the existing projections and on valid
cells only. Parameters are `ems_nan_max_abs` (default 1e20, the old threshold),
`ems_nan_abort_exit_code` (default 86) and `ems_nan_abort_prefix`.
Detection collects a first local witness without a conditional OpenMP barrier.
After all workers join it closes/renames a per-rank JSON record with field,
cell, level, time, dx, dt, phase, threshold, exit code, finite values in
exact hexadecimal notation and non-finite labels (not NaN payload bits), then directly calls MPI_Abort on Chombo's communicator.
The non-MPI fallback uses the same nonzero code. There is no abort collective
that could wait for a rank stranded in an exchange. If the record cannot be
written it prints an explicit write failure and still terminates nonzero.
The hook diagnoses the post-complete-RK check; it does not claim to identify
the first generating RK stage.

The isolated patch [t24-safe-abort.patch](t24-safe-abort.patch) applies only this
feature to d507438. It contains no later gauge or recording changes. Every
native unit was rebuilt against the new parameter layout: GRAMRLevel stores
SimulationParameters by value. The initial partial-object build was rejected
and corrected, not accepted as a regression.

Omitted/off/on native checks at t=0, after two regridding coarse steps and after
restart have **3,106,560 defined Float64 comparisons, zero bit mismatches**.
All evolved checkpoint/plot values, including saved evolution ghosts, and
valid diagnostic cells are compared. Only the previously demonstrated undefined
plot-diagnostic ghost slots are excluded. A NaN on the non-master OpenMP worker
preserves all 112 input Float64 bit patterns, closes the registered record,
and exits **86**; a finite input is unchanged and exits 0. No RHS, floor,
transfer, gauge, tagging or tracker operation is changed. Evidence:
[t24-controls-qualification.json](t24-controls-qualification.json),
[t24-control-bits.csv](t24-control-bits.csv),
[t24-control-resources.csv](t24-control-resources.csv).

The CH_MPI probe and the dependency's actual MPI SPMD unit compile with the
strict flags. **Every local one/two-rank attempt fails in MPI_Init**, exit 15,
at POSIX shared-memory bootstrap before the checker; the direct-abort path
is **not MPI-qualified**. [t24-mpi-abort-qualification.json](t24-mpi-abort-qualification.json)
states this failure rather than treating its nonzero code as an abort pass.
Changing the runtime SHM selector did not cure it. The controller's submission
worker still needs the two-rank/two-thread finite/default identity test and a
worker/rank NaN injection while the healthy rank waits in MPI, on an exclusive
node. No such cluster test was run here. The old cluster wait remains unproved.

The MPI choice is supported by the [MPICH API source](https://github.com/pmodels/mpich/blob/main/src/binding/mpi_standard_api.txt)
and its [abort notes](https://github.com/pmodels/mpich/blob/main/doc/mansrc/funcnotes.txt):
terminate the communicator's processes and perform the request outside worker
threads. That API argument is separate from the missing execution qualification.

## Cancelled builds and retained evidence

The controller cancelled the restart-dt feature and the annular stage-recorder
build after the interventions. **Neither was started.** Restart already
re-derives every level's dt from dt_multiplier; the completed local dyadic
readback regression is retained in [t24-restart-dt.csv](t24-restart-dt.csv)
and agrees with the cluster's 13-level readback. The annular recorder remains
a historical **design only**, explicitly cancelled as a build in
[t24-annular-recorder-design.md](t24-annular-recorder-design.md); no hook code
exists. Its old size estimate remains evidence, not a submission plan.

## Mechanism ranking after both interventions

1. **Collapsed-lapse shift/Gamma discrete CFL crossing:** supported by the
   code's alpha-independent loop, beta2/Gamma2 unstable mode, active-S history
   and two successful targeted interventions. Exact stage/checkerboard growth
   and uniqueness remain unproved.
2. **Inherited point-transfer/resolution damage:** remains possible as a
   contributor to the large inner metric and S; a switch only at step 152 does
   not exonerate the preceding 152 steps. Negative covered lapse is excluded
   as the proximate requirement by legacy's failure with positive valid lapse.
3. **Regrid/interpolation amplification:** remains a possible seed or contributor;
   frequent descendant regrids coincide with every fine-time neighbourhood.
   The named final regrid does not change count, and a lower-CFL run with the
   same transfer/regrid policy passes it. A unique footprint-change trigger
   is unsupported.
4. **Tracker/tag nontermination:** disfavoured as the initiating cause by the
   reproduced numerical non-finite witnesses and successful interventions.
   The unsafe abort is a separate failure to terminate the MPI job cleanly.

The gauge control does not establish accuracy, convergence, horizon retention
or a unique cause. No evolution fix beyond the optional termination mechanism
is applied in this turn.
