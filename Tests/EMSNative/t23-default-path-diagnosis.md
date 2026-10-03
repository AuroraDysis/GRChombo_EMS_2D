# T23 default-path diagnostic difference

Status: **QUALIFIED_UNWRITTEN_PLOT_DIAGNOSTIC_GHOSTS**. Exp-0025 admission remains blocked; no cluster retry, production fix or identity-check edit is performed. The cluster finding is verified against its own four native receipts (32 zero rank exits each) and failed identity receipt: 57,806 differences, six diagnostic fields, no evolved-field differences at plot0. Cluster HDF5 files were not collected, so the original per-cell locations cannot be reconstructed from the per-box receipt alone. Local controlled evidence supplies the cell classification, not invented cluster coordinates.

The source inspection finds a pre-existing hole in plot output. `Source/GRChomboCore/GRAMRLevel.cpp:907` allocates `plot_data` with three ghost cells and no initialization. With `write_plot_ghosts=true`, evolved fields use the boundary copier at lines 927–928, while diagnostic copyTo at lines 936–937 omits it. The final same-level exchange at line 967 cannot populate uncovered refinement ghosts or physical-boundary diagnostic ghosts. `Examples/EMS/EMSBH2DLevel.cpp:302–325` computes native constraints on valid cells only. This plot writer is unchanged from 5da576b through 758f0d2. A difference in heap reuse can therefore expose saved diagnostic values with no defined numerical meaning, even when evolved fields agree exactly. The controlled local receipts qualify this as an existing unwritten-output defect; they do not identify a deterministic first bad evolution commit.

The earlier T22 plot regression requested only 28 evolved components; it did not establish six diagnostic-ghost identities. The T23 matrix explicitly plots all 34 components, with plot ghosts enabled, in both the production main and dense-tag factory. It isolates 31176be, d507438, a pre-selector source with only the qualified T21 hook, and 758f0d2, plus every-unit access flags and the actual unused MovingGauge object linked into the old executable. Every native unit is rebuilt against its own source/parameter layout; cross-commit GRAMRLevel.o reuse would be invalid. Local GCC is 16 on macOS, whereas the cluster compiler is GCC 11.5 on Linux: local results test the source mechanism and do not certify that compiler/MPI combination.

The original matrix stopped at `build-candidate-native-flags`: own exit 1, no resource gate, because the candidate's checkpoint overrides call private base methods at `Examples/EMS/EMSBH2DLevel.cpp:712,729` (`GRAMRLevel.hpp:98,102`). The failed build is preserved under `/private/tmp/ems-t23/builds/candidate-native-flags.failed-private-base/`, with its own [failure record](t23-failed-build.csv). The Tests-only repair adds `-fno-access-control` to EMSBH2DLevel as well as the dense main; every other native unit retains normal access checks. Thus this control is explicitly **main plus the necessary level unit**, not a claim that the pinned candidate compiles with the old main-only scope. Production sources are not repaired. The five successful original builds were reused; all remaining cases completed under `/private/tmp/ems-t23/continuation/done.exit = 0`, each verified independently.

There is no deterministic first bad numerical commit in this matrix. The earliest source-stage contrast, 31176be versus 5da576b, already shows 3,144 undefined diagnostic-ghost differences at t=0; nevertheless the old/all-access executables and all their native objects have identical SHA-256 hashes. Running that **same executable content** twice produces 5,232 plot0 and 5,841 plot1 diagnostic-ghost differences in the production-main cases, with zero defined-state differences. The candidate all-access and narrow-access executables are also byte-identical. [Build identities](t23-build-identity.csv) therefore rule out an access-flag arithmetic change locally; source/build layout is not needed to generate the mismatch. The unused extra object can change undefined ghost output too, without changing any checked evolved or valid diagnostic value. The demonstrated cause is the unwritten buffer, while the precise allocator history of the original cluster process is not established.

All ordinary small-grid cases use levels 0–1, N=64x32, L=8 M, center=(4,0), h0=0.125 M, dt0=0.03125 M, one coarse step, the same E data, alpha_K, unchanged ExperimentalGauge/KO/point transfers and the positive-time static guard. If the ordinary small grid has no raw diagnostic difference, the matrix additionally initializes the real E-high max-level-14 hierarchy in both old/new dense variants with zero advances. It never launches an expensive max-level-14 coarse step. A separate private, environment-parametric buffer probe seeds only the newly allocated plot buffer with two distinct finite values. If these appear in saved diagnostic ghosts but nowhere in valid constraints/evolved data/checkpoints, that demonstrates the unwritten-output mechanism without subtracting or evaluating static data after initialization.

| comparison | file | evolved differences | valid diagnostic differences | diagnostic ghost differences |
| --- | --- | ---: | ---: | ---: |
| candidate-dense vs candidate-native-flags-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| candidate-dense vs candidate-native-flags-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| candidate-dense vs candidate-native-flags-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 2088 |
| candidate-dense vs candidate-native-flags-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 4380 |
| candidate-production vs candidate-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| candidate-production vs candidate-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| candidate-production vs candidate-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 0 |
| candidate-production vs candidate-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5019 |
| candidate-production vs candidate-native-flags-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| candidate-production vs candidate-native-flags-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| candidate-production vs candidate-native-flags-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 0 |
| candidate-production vs candidate-native-flags-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5784 |
| old-dense vs candidate-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs candidate-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs candidate-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs candidate-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs candidate-native-flags-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs candidate-native-flags-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs candidate-native-flags-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 2088 |
| old-dense vs candidate-native-flags-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 4380 |
| old-dense vs old-all-access-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs old-all-access-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs old-all-access-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 3144 |
| old-dense vs old-all-access-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 4380 |
| old-dense vs old-extra-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs old-extra-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs old-extra-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs old-extra-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 2187 |
| old-dense vs recording-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs recording-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs recording-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs recording-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 2807 |
| old-dense vs setter-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs setter-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs setter-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 3144 |
| old-dense vs setter-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 3110 |
| old-dense vs tracking-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs tracking-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs tracking-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 0 |
| old-dense vs tracking-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 2211 |
| old-production vs candidate-native-flags-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs candidate-native-flags-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs candidate-native-flags-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 3144 |
| old-production vs candidate-native-flags-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5841 |
| old-production vs candidate-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs candidate-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs candidate-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 3144 |
| old-production vs candidate-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 6177 |
| old-production vs old-all-access-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs old-all-access-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs old-all-access-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 5232 |
| old-production vs old-all-access-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5841 |
| old-production vs old-dense | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs old-dense | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs old-dense | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 3144 |
| old-production vs old-dense | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5841 |
| old-production vs old-extra-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs old-extra-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs old-extra-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs old-extra-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5841 |
| old-production vs poison-a | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs poison-a | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs poison-a | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 6480 |
| old-production vs poison-a | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 6480 |
| old-production vs recording-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs recording-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs recording-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 5232 |
| old-production vs recording-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5538 |
| old-production vs setter-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs setter-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs setter-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 3144 |
| old-production vs setter-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5841 |
| old-production vs tracking-production | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| old-production vs tracking-production | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| old-production vs tracking-production | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 5232 |
| old-production vs tracking-production | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 5538 |
| poison-a vs poison-b | EMS_000000.2d.hdf5 | 0 | 0 | 0 |
| poison-a vs poison-b | EMS_000001.2d.hdf5 | 0 | 0 | 0 |
| poison-a vs poison-b | EMS_Plot_000000.2d.hdf5 | 0 | 0 | 6480 |
| poison-a vs poison-b | EMS_Plot_000001.2d.hdf5 | 0 | 0 | 6480 |

Across completed comparisons: evolved bit mismatches **0**, valid diagnostic bit mismatches **0**, poison-dependent saved differences **12960**. The 20 matrix comparison cases cover 9,587,200 evolved Float64 value comparisons and 768,000 valid-diagnostic comparisons, with zero mismatches. These are array identities, not whole HDF5 file identities: the candidate intentionally adds gauge metadata. All 18 matrix native runs and two supplementary split-box poison runs produced both checkpoints and plots at t=0 and t=0.03125 M; [run validation](t23-run-validation.csv) seals the four outputs of each case. [t23-comparisons.csv](t23-comparisons.csv) records every field and region including zeros; [t23-location-summary.csv](t23-location-summary.csv) records differing boxes, levels, categories and maxima. The full per-cell ledger remains under `/private/tmp/ems-t23/ledgers/`; it is not committed. Physical-boundary/axis ghosts, uncovered refinement ghosts and covered same-level ghosts are distinguished using each saved hierarchy, not a chosen radial cut. Checkpoint0, plot0, checkpoint1 and plot1 are all compared; no early failure stops the later comparisons.

The [poison witness](t23-poison-witness.csv) is exact: at **each** plot clock all 6,480 changed values are 1.25 in one run and 2.5 in the other (difference exactly 1.25), in all six diagnostic fields. There are 612 coarse physical-boundary halo cells per field and, on level 1, 162 physical/axis halo cells plus 306 uncovered-refinement halo cells per field. The coarse valid box is i=0…63, j=0…31; the fine valid box is i=40…87, j=0…23 with physical faces x=±1.5 M, y=0…1.5 M. The poison occupies only the three-cell outside halos: coarse i=−3…66, j=−3…34 outside the valid box; fine axis j=−3…−1 and the lateral/top halos outside its valid box. No valid puncture cell changes. Ordinary old/new production plots have 3,144 ghost differences at t=0 and 6,177 at t=0.03125 M; the corresponding old/new dense plots happen to have zero differences at both clocks. Accidental equality of unwritten ghosts in one case is not a defined-output guarantee.

The smallest matrix grid has one box on each level and therefore has no internal same-level box interface. A supplementary, independently receipted poison pair uses the same private executable and grid with `max_box_size=32` to create internal interfaces. [Split-box comparisons](t23-split-comparisons.csv) and [qualification](t23-split-qualification.json) check **4032** covered same-level diagnostic ghost values with **zero** poison differences. Its **13824** changed values over both plots again consist exclusively of exact 1.25/2.5 pairs in uncovered physical/refinement halos; all **567488** defined-value comparisons are identical. Thus an internal box edge with neighbouring valid support is distinguished experimentally from an uncovered level face or physical/axis boundary. Its per-cell ledger is `/private/tmp/ems-t23/ledgers/split-cell-mismatches.csv`.

The production evolution RHS consumes evolved state, not the private plot buffer. An unwritten plot diagnostic ghost cannot feed it. This structural statement does not waive bitwise checks: any actual evolved or valid-diagnostic mismatch in the local matrix is reported as requiring further diagnosis. Undefined ghost output can occur at positive-time plots as well as t=0; whether it happens in the ordinary cases is a measured result in the tables. Defined valid constraints, their native formulas and term scales must be distinguished from these output halos.

For exp-0025, do not compare undefined diagnostic ghost values or interpolate them into the normalized magnitude gate. The registered masks, clocks and thresholds do not change. The exp-0023 reducer already strips output ghosts in `phase-3/audit.py:57–66`. Its existing phase-4 `EMSConstraintPieces.cpp` / `NativeConstraintPieces.impl.hpp` / `NativeGaussPieces.hpp`, with `production/pieces.py`, recomputes native valid-cell totals and individual term scales from saved numerical states without advances or static reads. `pieces.py:50–55` explicitly crops saved halos and checks evolved valid bits; lines 60–68 validate native totals against saved valid diagnostics and term sums with the existing 2048-epsilon scaled arithmetic allowance. Lines 70–92 use AMR-uncovered valid cells, the frozen masks, coordinate-volume weights, and the largest individual term's RMS on that same mask and clock. The T19 magnitudes consume these native tables; T21 supplies gauge-labelled checkpoint replay and current-state RHS/geometry, not a replacement normalized-constraint reducer by itself. One common, source/compiler-pinned native constraint tool can replay H and G at all 13 clocks, checking untouched evolved valid bits and the inherited normalization and Gauss shell integral. Existing H replay validation is evidence for H, not a waived validation for G. The reducer must fill private evolved stencils consistently and use valid/common-supported native constraints; no plot-ghost constraint is an input. Gauge metadata must be honoured when using a gauge-aware restart wrapper, although the native constraint formulas do not depend on the driver package. This preserves comparability without changing the screen or declaring the new gauge accurate.

The old writer is defective for saved diagnostic halos; the old and new valid constraints are the defined numerical quantities. Proposed controller action: retain strict bitwise evolved checkpoint/plot comparisons, including their defined ghosts, and strict bitwise comparison of all six diagnostic fields on every **valid** box cell and covered same-level ghost at both t=0 and one step. Check layouts, clocks, component names and gauge metadata separately, recognising the intentional new metadata. Report diagnostic ghosts outside same-level valid support as undefined rather than passing or failing them numerically. The cluster should classify its original differing indices before changing admission: the collected per-box receipt cannot establish those indices, so this local result does not itself lift that blocker. A future diagnostic-output repair would need explicit ghost definitions and separate qualification; copying uncomputed diagnostic source ghosts or silently zeroing them would not establish meaningful constraints. No production patch is applied, especially none to the in-flight exp-0024 fork d507438. [t23-plot-buffer-probe.patch](t23-plot-buffer-probe.patch) is private instrumentation only, not a proposed evolution or admission-check change.

Measured per-process RSS maximum: **610762752 bytes** (0.611 GB, below 6 GB). All builds/runs are serial with two OpenMP threads and single-threaded numerical libraries, a 6 GB per-process cap, 7 GB process-tree cap and 16 GB output gate. Every completed child needs its own exit/resource receipt and native output checks; a launcher exit is insufficient. The original `/private/tmp/ems-t23/pipeline/done.exit = 1` remains as the failed-build history; the successful continuation is `/private/tmp/ems-t23/continuation/done.exit = 0`. No SSH, cluster work or commit is performed. Regenerate numerical comparisons with `python Tests/EMSNative/t23-compare.py analyze`; once the continuation and qualification receipts are complete, `python Tests/EMSNative/t23-compare.py finalize` checks all own receipts, refreshes the report/README and seals the packet. Tests instrumentation and builds are confined to `/private/tmp/ems-t23/`.
