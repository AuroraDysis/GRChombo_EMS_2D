#!/bin/zsh
# All runs are frozen or t=0; deliberately sequential HDF5 restart access.
set -eu
cd "${0:A:h}/../.."
source Tests/EMSNP/config.sh
out=Tests/EMSNP/results-t3
logs=$EMS_NP_LOG_ROOT/t3
mkdir -p "$logs"
measure="$PWD/Tests/EMSNP/measure_t3.py"
python3 Tests/EMSNP/prepare_t3.py
python3 "$measure" make -C Tests/EMSWeylOffline all DIM=2 CXX="$CXX" -j 2 > "$logs/build-release.log" 2>&1
python3 Tests/EMSNP/syntax_check.py "$out" > "$logs/syntax-driver.log" 2>&1
env OMP_NUM_THREADS=1 "Tests/EMSWeylOffline/EMSNPTest${EMS_EXE_SUFFIX}" \
    "$out/rn-n256/params.txt" > "$logs/unit.log" 2>&1
for pair in 'np_static_metric static-cas' 'np_first_order first-order-cas' 'np_tracefree tracefree-cas'; do
    parts=(${=pair})
    python3 "$measure" "$GTIMEOUT" -k 5s 60s \
        "$WOLFRAM" -script "scripts/cas/$parts[1].wls" \
        > "$logs/$parts[2].log" 2>&1
    # Kernel parse errors can exit zero: require the actual closed witness.
    rg -q '^PROVED' "$logs/$parts[2].log"
    if rg -q 'NOT-CLOSED|Syntax::|Check::' "$logs/$parts[2].log"; then exit 1; fi
done
for spec in 'rn RN-identity 1 2 4' 'ems a0.9-e7 0.0512 0.2622 0.6069 1.0429'; do
    parts=(${=spec})
    python3 "$measure" "$JULIA" --project="$EMS_ROOT" \
        Tests/EMSNP/static_oracle.jl "$EMS_ROOT/artifacts/echo/$parts[2].trumpet" \
        "$out/$parts[1]-oracle.csv" $parts[3,-1] > "$logs/$parts[1]-oracle.log" 2>&1
done
for case in rn-n256 rn-n512 ems-n1024 ems-n2048 checkpoint cavity default; do
    exe=EMSNPStatic
    if [[ "$case" == checkpoint || "$case" == cavity || "$case" == default ]]; then exe=EMSWeylCheckpoint; fi
    (
        cd "$out/$case"
        python3 "$measure" "../../../EMSWeylOffline/${exe}${EMS_EXE_SUFFIX}" params.txt
    ) > "$logs/$case-release.log" 2>&1
done
dependency_status=0
env OMP_NUM_THREADS=1 "Tests/EMSWeylOffline/EMSNPStatic${EMS_EXE_SUFFIX}" \
    "$out/bad-activation.params.txt" > "$logs/activation-dependency.log" 2>&1 || dependency_status=$?
[[ "$dependency_status" == 2 ]]
rg -q 'ems_np_activate requires ems_radiation_activate' "$logs/activation-dependency.log"
python3 Tests/EMSNP/prepare_t3.py --verify-only > "$logs/input-final.log"
report_status=0
python3 Tests/EMSNP/report_t3.py > "$logs/report.log" 2>&1 || report_status=$?
cat "$logs/report.log"
print "RUN_COMPLETED validation_exit=$report_status"
exit $report_status
