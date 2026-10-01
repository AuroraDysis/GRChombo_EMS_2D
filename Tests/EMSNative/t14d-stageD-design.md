# Stage D draft — for the operator's admission decision

This is a proposed fresh alpha_K E global three-grid chain, not an authorization or a launched run. T14d's complete registered history is INCONCLUSIVE. Its unresolved launch entries remain part of that verdict. This draft is required independently of admission; it proposes no additional interpolation variant and no change to the gauge, KO, equations, transfers, Float64 precision, puncture treatment, boundary prescription or author's RHFinder.

Use exp-0020's domain x relative to the hole in [-336,336] M, y in [0,336] M, centre (336,0), point transfers, sigma=1 and dt_multiplier=0.25. Use fresh E-low/mid/high initializations with ems_use_maximal_initial_lapse=true. Keep every exp-0020 level-0–12 physical face R_l=112/2^l M and add the same two puncture levels to **all three** global rungs: level 13 at 0.013671875 M and level 14 at 0.0068359375 M. This is a uniform 3/2 global-spacing chain, not the factor-two source-only T14 ladder. Extra levels put all three finest spacings within or below T14's tested source-refinement range; they do not establish exterior convergence. Using max_level=12 instead would cost less but would place E-low's source spacing outside T14's three-rung range. Do not splice exp-0020/0021 original-lapse outputs into this new chain.

| Rung | h0 / M | h14 / M | Steps at 2.625 / 4.375 / 10.5 M | Model node hours |
| --- | --- | --- | --- | --- |
| E-low | 0.875 | 5.34057617e-05 | 12 / 20 / 48 | 13.64 |
| E-mid | 0.583333333 | 3.56038411e-05 | 18 / 30 / 72 | 12.95 |
| E-high | 0.388888889 | 2.37358941e-05 | 27 / 45 / 108 | 12.47 |

Base N1/N2 are 768/384, 1152/576 and 1728/864. No regridding is allowed. Before admission to execution, run the real initialization-only census on every proposed grid, prove inherited unions/faces unchanged, check block alignment and all valid/ghost floor margins, and qualify the Tests-only dense-tag initializer workaround on all three grids. Reuse T14's workaround only after these actual censuses; it addresses initialization set storage without changing evolution. Retain absent-static-file complete-hierarchy restart and recorder/default-path identity controls, and measured <=3 GB RSS per process with <=4 OpenMP threads. Initialization uses the static EMS file only at t=0; no positive-time initialization or static ghost substitution is permitted.

## Fixed domain, masks and time claims

Predeclare the smooth puncture-excluded interior domain **0.005<=R<=20 M**. R_excl=0.005 M was fixed in t14d-registration.json before this re-screen. T14's longest measured half-height span is approximately 0.001712 M and every lobe is clipped by W=[0.00075,0.0025] M. Thus these spans are lower bounds, not complete widths. The chosen exclusion is twice W's outer radius and roughly three times the observed span, allowing a conservative fixed source neighbourhood while remaining inside the fixed horizon shell's 0.0060286 M inner radius. It is not a claim that the whole disturbance fits inside that radius. Once an outgoing or reflected disturbance enters R>=0.005 M, it participates in every declared test; the mask never moves with the pulse.

Retain every exp-0020 mask below, including full masks containing faces/corners and the outer boundary shell. The additional puncture-excluded mask and receiving-side mask do not replace failed masks. Use composite coordinate-volume RMS with weights 2pi*y*h_l^2 over uncovered cells, plus peak norms. Retain independent axis/equator/diagonal rays and fixed patch-edge, convex-corner, axis-junction, same-level seam and smooth-interior classes. Diagnostic support must be declared from the coarsest rung, and it must be identical across rungs. Report absent common classes explicitly.

| Mask | Fixed radial interval / M |
| --- | --- |
| puncture_excluded | 0.005–20 |
| horizon | 0.00602861293921588–0.0120572258784318 |
| inside_inner_ring | 0.0120572258784318–0.0267589278192644 |
| cavity | 0.0267589278192644–0.117265227037402 |
| between_rings | 0.0267589278192644–0.738529896105975 |
| outer_ring | 0.59082391688478–0.886235875327169 |
| far | 4–8 |
| cavity_core | 0.055–0.075 |
| outer_ring_core | 0.85–1 |
| far_core | 6.5–8 |
| exterior_wake | 1–20 |
| receiving_side_3p5 | 3.5–4.375 |

The outer_boundary_shell remains distance 0–1 M from x=+-336 or y=336 with y>1 M. The old Sommerfeld shell failure remains an admission issue. A restricted interior result additionally requires an all-relevant-mode arrival bound and the previously required boundary-position control. Without that evidence this draft admits neither a domain-wide claim nor a 100 M extension.

The predeclared history start is the **first positive T14d common clock at which all four fields, both rays and both norms exceed 5x interpolation**, computed as step 10, t_start=0.00035603841145833332 M. This numerical start is frozen for any fresh Stage D qualifying run. Later T14d exceptions at steps [13, 14, 15] are retained: a first qualifying clock does not establish a fully qualified suffix and cannot change the T14d verdict. In Stage D report all t=0 and early histories, but any field/constraint convergence claim beginning at this start is explicitly narrower than convergence over the whole evolution. Do not move the start if a later clock fails. Preserve the 0.875 M synchronized exterior cadence and exact 2.625, 4.375 and 10.5 M diagnostics. Add bounded early current-state captures and diagnostic-only samples at the fixed start through the first coarse interval; if sampling these subcycled times requires time interpolation, predeclare it and measure its error with the control rather than silently comparing asynchronous data.

## Expected order and complete evidence

Design order remains **four** on the fixed smooth domain. No independent regularity argument licenses a lower expected order here. For the 3/2 chain, report log(||D_LM||/||D_MH||)/log(1.5), the signed fourth-order residual D_LM-(3/2)^4 D_MH, and both direct-to-zero constraint orders. Retain all 28 evolved fields, not only Gamma or lapse. Use the fixed I8 central diagnostic reconstruction with I8/I10 sensitivity, >5x pairwise significance, and a fresh E-high dt/2 temporal control through 10.5 M requiring temporal error <=0.2 of both spatial differences. At faces, use declared valid one-sided support; verify coverage and parity before qualifying results. A sampling/floor-limited row is unresolved, never a fourth-order pass. These samplers modify diagnostics only; evolution transfers and native derivative formulas stay frozen.

The constraint battery includes Hamiltonian, Mom1/Mom2 and their joint norm, C_Gamma components and norm, spatial Z1/Z2/Z from the code's native current metric/Gamma, Theta, det(h)-1, tr(A), electric and magnetic Gauss constraints, and the carried cleaning variables Lambda and Xi. Preserve the implementation's prescribed dimensional scalings, and report raw absolute norms, maxima, cell counts and quadrature weights alongside scaled norms. There are no invented first-order-reduction constraints. Exact zero constraints have undefined logarithmic order. Require both pair orders consistent with design four after measured interpolation, temporal and Float64 uncertainty, and decreasing residuals throughout the declared histories. No absolute constraint admission tolerances are supplied by this task or the exp-0020 contract: the operator must seal those physical tolerances for each scaled constraint before execution; measured orders alone cannot admit the chain. Do not choose tolerances after seeing its outputs.

Horizon retention always covers numerical **t=0 through 10.5 M**, independently of t_start. Re-find each rung's numerical t=0 surface using only saved current fields and a numerical seed; use it for its own A_H(0), Q_H(0) baseline. Keep the author's finder unchanged. Require negative inward expansion, fresh final RMS outgoing expansion <=1e-6, angular-resolution and stopping sensitivities below the drift budgets and inter-rung differences, and qualified coverage of between-cadence excursions. Retain native candidate histories every coarse step and numerical checkpoints at excursions; a FAR/CLOSE label or finite candidate area is not a horizon certificate. The budgets remain max|Delta A/A(0)|<=1e-3 and max|Delta Q/Q(0)|<=2e-4, with decreasing drifts under global refinement. Independent current-field sphere charges at 20,50,100 M and their Gauss/flux/cleaner budgets are reported separately. A later 100 M decision must preserve the same numerical t=0 baselines and full horizon interval; this pilot supplies no 100 M claim.

## Cluster work estimate

Use one node per leg, **32 MPI x 4 OpenMP** (128 allocated cores), as exp-0020. Its measured level-12 median seconds/coarse-step were 256/162/104, within the controller's approximately 100–260 s range. The valid-cell subcycle work model W=sum(N_l*2^l), including evolved covered cells, gives 604,495,872 -> 2,416,435,200 cell-steps per coarse step for E-mid, a factor **3.997439** with the extra levels. The same factor applies to low/high since every refined level has 32768/73728/165888 valid cells and each base has nine times that count. This is a pricing model, not a measured max14 cluster rate; communication, load balance, setup, finder and I/O overhead require fresh bounded cluster smokes before a sealed submit budget.

The nominal three-rung 10.5 M chain costs approximately **39.07 node-hours / 5001 core-hours**, or **25.32–65.82 node-hours** using the 100–260 s baseline range scaled by the work factor. A full finest dt/2 control adds approximately **24.94 node-hours**, with range **23.98–62.36**. Total nominal plus control is **64.01 node-hours / 8194 core-hours**, with planning range **49.30–128.18 node-hours** before uncalibrated overhead. Reuse one physical chain per rung through the earlier diagnostic times rather than paying for three separate endpoints.

Propose four nominal segments of 2.625 M each (12/18/27 coarse steps), with complete-hierarchy checkpoints at every segment boundary and all inherited 0.875 M diagnostics. Split the finest dt/2 control into eight 1.3125 M segments (27 control coarse steps each); retain its exact 4.375 M diagnostic inside the corresponding segment. The range model prices the longest proposed segment below eight node-hours before overhead, but seal allocation limits only after new smokes. Resume only from the complete numerical checkpoint with the positive-time static guard enabled. Keep physical completion distinct from diagnostic qualification; a consumed segment budget, failed restart or missing qualification marker cannot imply admission. No SSH, submission, evolution, finder call or commit is performed by T14d.
