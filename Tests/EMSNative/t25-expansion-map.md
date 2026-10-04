# T25 — signed outgoing expansion on frozen numerical checkpoints

**READY-EXCEPT.** The offline tool, native sign check, local initial-data
regeneration, three requested maps and bracket-seeded searches are complete.
The unboosted N48 search passes all three stages. Both boosted binary holes
have sampled pointwise trapped/untrapped barriers at initialization, but their
bounded N48 searches remain unqualified. At step 152 none of the sampled
individual or common families has a pointwise barrier; there are ten weaker
mean-only crossings. These observations neither locate a late horizon nor
establish that a horizon shrank or disappeared. Early positive-time checkpoints
are on the clusters and have not been inspected locally.

## What is evaluated

The map computes the signed outgoing expansion

\[
\theta_+=D_i s^i+K_{ij}s^is^j-K,\qquad
\gamma_{ij}=\widetilde h_{ij}/\chi,\qquad
K_{ij}=(\widetilde A_{ij}+K\widetilde h_{ij}/3)/\chi.
\]

The code uses the convention
\(K_{ij}=-(\partial_t\gamma_{ij}-{\cal L}_\beta\gamma_{ij})/(2\alpha)\):
the metric rate contains \(-2\alpha\widetilde A_{ij}\)
([CCZ4Cartoon.impl.hpp:299](../../Source/Cartoon/CCZ4Cartoon.impl.hpp:299)).
This gives the same outgoing sign as the author's
[RHSurf.hpp:430](../../Source/RHFinder/RHSurf.hpp:430). The lapse is reported
as a surface mean and does not enter the expansion or the horizon equation.

An axisymmetric surface has cylindrical semiaxis \(a\), axial semiaxis \(c\)
and centre \(x_c\). Its implicit equation is
\(F=(x-x_c)^2/c^2+(\rho^2+w^2)/a^2-1=0\).
The normal is the outward normalized gradient of this equation. Its covariant
Hessian gives
\[
D_i s^i={ (\gamma^{ij}-s^is^j)
 (\partial_i\partial_jF-\Gamma^k{}_{ij}\partial_kF)
 \over \sqrt{\gamma^{ij}\partial_iF\partial_jF}}.
\]
Spheres are the case \(a=c\). The geometry is exactly the requested
\(\rho=a\sin\vartheta,\ x-x_c=c\cos\vartheta\); sampling uses the equivalent
polar radial graph
\(f(\theta)=[\cos^2\theta/c^2+\sin^2\theta/a^2]^{-1/2}\).
Its tangent, normal and surface Hessian are analytic. Angular derivatives of
the surface are therefore not inferred from nearest cells or finite angular
differences in the primary map. Cell-centred polar angles avoid the coordinate
poles; the cartoon transverse metric derivatives are retained, including
\(\partial_w h_{xw}=h_{x\rho}/\rho\) and
\(\partial_w h_{\rho w}=(h_{\rho\rho}-h_{ww})/\rho\).

[t25-expansion-map.cpp](t25-expansion-map.cpp) uses the unchanged frozen
RHUnion field query and `AMRInterpolator<Lagrange<4>>`, including its `dx` and
`dy` derivative queries. Ownership is determined by the interpolator's own
finest **valid-box** search, not by sampling nearest-cell values. Stencil
ghosts are filled from the checkpoint's valid evolved fields with the native
production transfers and boundaries. No plot diagnostic ghosts are read.
The map consumes the actual checkpoint hierarchy, not a prescribed hierarchy
template. The old evolution buffer and unused diagnostics are released per
level; only evolution ghosts are filled for maps. No advance or regrid is
possible: the Tests-only level rejects those entry points. The map runs with
one process and one thread. Static and CTT paths are deliberately absent in
every map parameter file; they are used only by the separate t=0 initializers.

The physical metric is inverted directly. The author's finder instead uses
unit-determinant conformal cofactors and the corresponding volume derivative
([RHSurf.hpp:315](../../Source/RHFinder/RHSurf.hpp:315),
[RHSurf.hpp:500](../../Source/RHFinder/RHSurf.hpp:500)). No determinant or
positivity correction is made by the map. A non-positive interpolated chi or
non-positive metric is flagged `INVALID_GEOMETRY`. The initial unboosted
sphere at the linear bracket root agrees with the author's expansion to
\(1.88\times10^{-11}\) pointwise; across its sampled spheres the maximum
disagreement is \(2.94\times10^{-11}\). Initial binary spheres agree to
\(5.67\times10^{-11}\). This confirms the convention on the native initial
control. It is not a promise of identical geometry prescriptions at late
times: interpolation near the late punctures gives determinant deviations as
large as 0.203, and the largest resolved-sphere expansion disagreement with
the finder is 0.147. Both values are retained in the CSVs. No late numerical
surface is certified by the initial sign check.

The charge uses the native EMS displacement flux
\[
Q={1\over\sqrt{2\pi}}\oint
\exp[-2\alpha(f_0+f_1\phi+f_2\phi^2)]\,E_i s^i\,dA.
\]
The stored electric components are **covariant**; this contraction is
identically the requested \(E^i s_i\).
[RHSurf.hpp:584](../../Source/RHFinder/RHSurf.hpp:584) implements the same
coupling and normalization, while
[EMSBH_trumpet_read.impl.hpp:367](../../Source/InitialConditions/EMSBH/EMSBH_trumpet_read.impl.hpp:367)
lowers the electric density with the physical metric. Area, charge and all
means use physical area weights; RMS is \(\sqrt{\int\theta_+^2dA/A}\).
The negative fraction is an **area fraction**, not a count of angles.

Each row reports both local minimum and maximum spacing, x extents, cylindrical
extent, determinant error, and enclosure of each supplied numerical puncture.
Resolution is conservatively the smallest semiaxis divided by the **largest**
local spacing on the surface. A value below 3 is `UNRESOLVED` and never enters
a bracket or a finder seed. Three cells is the requested eligibility rule,
not an accuracy certificate for strongly varying geometry.

[t25-brackets.py](t25-brackets.py) emits the consecutive-input-surface table.
A pointwise pair requires negative expansion at every sampled inner angle
and positive expansion at every sampled outer angle (the reverse direction
is also identified). A mean sign change without these inequalities is
`AVERAGE_ONLY_WEAKER`. Dense sampling can insert a mixed-sign surface between
uniformly negative and uniformly positive surfaces. The additional
`t25-barriers-*` tables retain those tight enclosing barriers, record their
adjacency and the number of intervening surfaces, and discard no input rows.
All intervening geometry must be resolved. These are **sampled** pointwise
inequalities at N96 and N192, not a proof at every continuous angle or of a
continuum MOTS. A mean root is a seed or a trial surface, never a horizon
measurement by itself.

## Validation and initial control

Both controls were regenerated locally with the pinned T17 executable and
`max_steps=stop_time=0`, geometric initial lapse and the static-file guard on.
No evolution occurred. The unboosted control has L0–12, h0=1.75 and
h12=0.00042724609375 M_i, centre x=40 on its L=112 domain. The initial binary
is the T17 d=32 control, centres 2224/2256 in the L=4480 domain, with the same
per-hole spacing and approximately fourteen cells per horizon radius. It is
distinct from the d=16 merger state at t=133. The initial binary companion is
the documented diagnostic n32-r6 file with SHA-256 `46253cb8…`; its data
certification status is unchanged. Full input hashes are in
[t25-audit.json](t25-audit.json).

The implicit-surface formula has an exact flat-space algebra witness and an
independent 50-digit curved-metric check at sixteen sampled points, with worst
residual \(1.07\times10^{-50}\). These checks verify the geometry algebra,
not the accuracy of the evolved data or existence of a horizon
([t25-geometry-verification.json](t25-geometry-verification.json)).

For the single hole the N192 pointwise bracket is r=0.00604–0.00606 M_i.
Its inner expansion range is approximately
[-6.947e-4, -2.848e-4], and the outer range is [0.019128, 0.019537].
Linear interpolation of the area-weighted signed means gives
**r=0.006040488470 M_i**. Evaluating a new numerical surface at that radius
gives **A=0.2245141754, Q=1.0693233457**. It has 14.138 cells per radius,
signed RMS 9.826e-5, and mixed tiny angular signs; it is a bracket-derived
trial sphere, not yet the final fitted surface. N96 to N192 changes the
linear root by 4.87e-10 M_i and A at the linear root by approximately 7.5e-6.

**The card's area benchmark 0.2218 does not describe the supplied progenitor.**
The pinned file's header has A_H=0.22451077798968849, 1.222% higher. The map's
root area agrees with that header and the prior T17 native control. That
header is inspected only as t=0 provenance; it is never an input, seed,
subtraction or target in the expansion map. The discrepancy is reported,
not removed by changing the input or the area normalization.

The unchanged author's N48 finder, seeded from the measured bracket root,
passes squared-expansion stages 1e-7/1e-10/1e-12. Stage residuals are
8.4655e-9, 9.9780e-11 and **9.9961e-13**, using 1/1411/2019 updates and
61.162 s total finder time. Its fitted radial range is
0.006040479697–0.006040494153 M_i; final A=0.2245517229 and Q=1.0695022343
at N48. The N48/N192 area and charge differences include angular quadrature
and surface fitting; they are not spatial error bounds. This is a strict
N48 all-stage pass as requested, **not a newly qualified N48/N96 pair**.
The prior 9.99e-13 result is reproduced without using its fitted shape as a
seed. See [t25-root-single.csv](t25-root-single.csv) and
[t25-finder-single-N48.csv](t25-finder-single-N48.csv).

## Binary at t=0

The initial map contains 1,200 surfaces at each N, spanning two holes,
five centre offsets (-0.001, -0.0005, 0, +0.0005, +0.001 M_i), five axial
aspect ratios c/a (0.85, 0.95, 1, 1.05, 1.15), and twenty-four radii
0.0013–0.1 M_i. Twenty surfaces are unresolved; all remaining surfaces have
valid metric geometry. The fifty consecutive-family sign changes are
mean-only. Twenty-two families nevertheless have tight enclosing pointwise
barriers once intervening mixed-sign surfaces are retained explicitly. The
counts are unchanged at N96 and N192. A shifted centre or spheroid is
therefore **not required for a sampled barrier at t=0**.

The tight centred-sphere pair is identical by reflection at the two holes:

| surface | r (M_i) | signed minimum | signed maximum | area-weighted mean |
|---|---:|---:|---:|---:|
| inner barrier | 0.00608 | -0.0801841 | -0.00735877 | -0.0567711 |
| intervening mixed sphere | 0.00610 | -0.0601846 | +0.0125824 | -0.0367915 |
| outer barrier | 0.00620 | +0.0385480 | +0.1110235 | +0.0618324 |

Linear interpolation across the uniform-sign pair gives
r=0.006137439565 M_i. A new N192 surface at that radius has
A=0.2232824138, Q=1.0693235191, expansion range
[-0.0229732, +0.0496843], mean +3.806e-4 and RMS 0.0208971.
It has 14.365 cells per radius. This signed mean zero estimate is not a
pointwise root. The centred sphere was selected for one bounded search at
each hole because it already supplies the tight resolved barrier and the
least sampled RMS among the near-root surface families; no finder sweep was
launched. Search seeds and their numerical bracket proofs are retained.

| hole | stopping status | updates | finder seconds | final squared expansion | trial A | trial Q |
|---|---|---:|---:|---:|---:|---:|
| left | TIME_CAP | 2336 | 235.029 | 2.27087e-5 | 0.223318623 | 1.069502217 |
| right | TIME_CAP | 2349 | 235.075 | 2.24868e-5 | 0.223318619 | 1.069502217 |

Neither search reaches even the first squared-residual threshold 1e-7;
later stages were not attempted. The differences in update counts reflect
wall caps under concurrent serial searches, not different physics. Both
residuals continue decreasing at the cap. They are about 155 times below
the historical binary value 3.5e-3, and below the historical boosted-single
3.4e-5, but these are **different seeds and update budgets**. The improvement
is an observed bounded-search result, not a spatial or angular convergence
claim. No new boosted-single search was performed.

The initial failure is not absence of a resolved numerical trapped/untrapped
pair in these families. A spherical radius inherited from the unboosted hole
is an inferior binary seed here. The result supports a seed/cap limitation
of those earlier searches, but does not isolate it from the author's flow,
recentring or discretization of the fitted surface. Lapse does not enter
this equation, so changing only initial lapse cannot change a frozen-state
MOTS. All binary stopped areas/charges above remain trial-surface readings.
See [t25-barriers-binary-N192.csv](t25-barriers-binary-N192.csv),
[t25-finder-stages.csv](t25-finder-stages.csv) and the raw retained shapes.

## Step 152, t=133 M_i

The sealed checkpoint supplies the numerical geometry. Track centres are
2239.7740195 and 2240.2259805 (separation 0.451961 M_i); midpoint 2240.
The 188 surfaces per N include per-hole radii from 3h12 to 0.1 M_i,
aspect ratios 0.75/1/1.5/2, and midpoint radii 0.25–8 M_i with aspect
ratios 1/1.5/2/3. All midpoint surfaces enclose both tracked punctures.
Two narrow oblate surfaces fail the three-cell semiaxis rule; 186 surfaces
remain eligible. There is **no sampled pointwise pair** in any eligible
family at either angular N. Ten mean-only crossings persist.

For centred individual spheres, r=0.047 has expansion range
[-0.509753, +6.84821], mean +0.069511 and negative area fraction 0.6134.
At r=0.08 the range is [-0.498515, +7.41647], mean +0.012320 and negative
fraction 0.6482. Even near the weaker mean crossing at r=0.09–0.1, the
surface retains both signs. These do not provide a resolved individual
horizon location.

| midpoint sphere r (M_i) | trial area | trial Q | signed min | signed max | mean |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 283.479 | 2.01855 | -0.399366 | +0.598692 | +0.057524 |
| 0.30 | 370.805 | 2.02678 | -0.481950 | +0.390505 | -0.144610 |
| 0.50 | 595.497 | 2.04383 | -0.639221 | +1.02770 | -0.305073 |
| 1.00 | 811.245 | 2.06636 | -0.573088 | +1.58527 | -0.315138 |
| 8.00 | 1252.782 | 2.09229 | -0.399125 | +0.583533 | +0.077799 |

The midpoint sphere mean changes sign at 0.25–0.3 and again at 6–8 M_i.
Prolate families have weaker mean crossings at a=4–6 (c/a=1.5),
a=3–4 (c/a=2), and a=0.25–0.3 and 2–3 (c/a=3). These are not pointwise
barriers or MOTS areas. No late finder search was launched. This finite
surface family cannot exclude a more general apparent horizon; it does
exclude treating these mean sign flips or coordinate-sphere areas as an
already measured individual/common horizon. The interval t=0 to t=133
does not determine when tracking became unsupported. The unchanged tool
can now be applied to the retained early checkpoints on a compute node.

![Signed expansion envelopes](t25-signed-expansion.png)

The shaded envelopes show angular minima/maxima and the line shows the
physical-area-weighted mean for centred spheres at N192. The initial
uniform-sign barriers are narrow, while the late families retain positive
and negative sectors across the mean crossings. No static solution is
plotted or subtracted.

## Reproduction and qualification packet

Build against a clean fork snapshot with the site's **serial** Chombo/HDF5
configuration. From `Tests/EMSNative`, the UCPH and Leonardo build line is
recorded in the C++ header and [t25-build-node.sh](t25-build-node.sh):

```sh
sh t25-build-node.sh ../.. "$CHOMBO_HOME" "$TMPDIR/t25-build" CXX=g++ OPT=HIGH
```

Use the loaded Leonardo GCC 11.5.0 `g++`; on UCPH use that site's configured
GCC/HDF5 compiler. The script defaults to MPI=FALSE, OPENMPCC=FALSE,
PRECISION=DOUBLE, USE_64=TRUE and -j1. The tool also enforces one MPI process
and one OpenMP thread. The build directory contains only the T25 main and
its own GNUmakefile, so Chombo's recursive make cannot pick up other Tests
mains. Native coverage access is enabled in that Tests translation unit
only. No production main or initial-data setter is linked. The site's
`Make.defs.local` supplies HDF5/link flags; serial Chombo libraries must
match this configuration. No remote build or cluster job was executed here.

The local isolated make build passed against the clean archive, using the
available OpenMP library configuration with runtime one thread. It peaked
at 0.539 GB RSS. Its maps and the strict direct-object build's maps are
bit-identical in every recorded value and angular row. This is local build
qualification, not an unperformed GCC11 cluster build. Exact local commands,
objects, library and author-source hashes are in [t25-build.json](t25-build.json).
The pending NaN-abort source changes were not used or changed.

Create a checkpoint parameter file from that run's own parameter file:
keep its domain, parity, transfer and coupling values; set `restart_file`
to the absolute checkpoint path and `hdf5_subpath = ""`, disable all live
tracking/extraction, and supply absent
static/CTT paths. Never increase its `max_level` beyond the stored hierarchy.
Surface CSV columns are `id,family,centre,a,c`; labels and paths must not
contain CSV commas. Each nested family keeps its centre/aspect ratio fixed.

```sh
./t25-expansion-map.ex checkpoint-params.txt \
  map_surfaces=surfaces.csv map_points=192 \
  map_output=map.csv map_angles=angles.csv map_punctures="x1 x2"
python3 t25-brackets.py map.csv brackets.csv --barriers barriers.csv
```

The C++ tool and standard-library bracket script are cluster-portable without
Python numerical packages. `map_angles` retains every signed angular sample,
weight, local level, spacing, lapse/chi and finder comparison. For a
time-resolved series, invoke the same program on each checkpoint with its
own contemporaneous numerical centres, then concatenate the per-checkpoint
CSVs. No time history is inferred from a late snapshot. The local
[t25-work.py](t25-work.py) records the exact initialization/map/seed/search
commands; its local paths are convenience defaults, not cluster input paths.

Optional `map_find=true` takes **one** seed surface and a matching
`map_finder_bracket` CSV. The tool rejects unmatched checkpoint/time/centre,
out-of-bracket axes, non-pointwise proof or unresolved seed before any
finder update. Parameters control angular N, chase, update cap, wall cap,
floor window and all three thresholds; the author algorithm is untouched.
The recorded searches use chase=quota=1, cap=100000, floor_window=100001,
235 s and thresholds 1e-7/1e-10/1e-12. A TIME_CAP or FLOOR is not a pass.

The measured peak RSS is **1.9999 GB for maps**, 1.932 GB for the binary
bounded finder and 1.268 GB for t=0 initialization. All run caps are met;
each job's own receipt and completion text were checked
([t25-receipts.csv](t25-receipts.csv)). Final-build single/late replays and
the isolated portable-build replay are bit-identical; positive and negative
bracket-guard cases and the unresolved exclusion check pass
([t25-regression.csv](t25-regression.csv)). No run is pending. Earlier local
setup failures and their receipts are retained rather than counted as
evidence. Nothing was evolved, submitted, SSH-accessed or committed.

Small surface/bracket/finder CSVs are in this packet. Angular dumps, native
receipts, regenerated checkpoints and stopped shapes remain under
`/private/tmp/ems-t25/`; their hashes and sizes are in `t25-manifest.txt`.
No bulk checkpoint is copied into the worktree. The required next
time-history reading is precisely the signed map through the retained early
checkpoints, not an extrapolation of a fixed-area radius.
