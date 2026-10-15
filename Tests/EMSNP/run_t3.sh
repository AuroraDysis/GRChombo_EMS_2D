#!/bin/zsh
# All runs are frozen or t=0; deliberately sequential HDF5 restart access.
set -eu
cd "${0:A:h}/../.."
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
export OMP_NUM_THREADS=8
out=Tests/EMSNP/results-t3
measure="$PWD/Tests/EMSNP/measure_t3.py"
python3 Tests/EMSNP/prepare_t3.py
python3 "$measure" make -C Tests/EMSWeylOffline all DIM=2 -j 2 > "$out/build-release.log" 2>&1
python3 Tests/EMSNP/syntax_check.py "$out" > "$out/syntax-driver.log" 2>&1
env OMP_NUM_THREADS=1 Tests/EMSWeylOffline/EMSNPTest2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
    "$out/rn-n256/params.txt" > "$out/unit.log" 2>&1
for pair in 'np_static_metric static-cas' 'np_first_order first-order-cas' 'np_tracefree tracefree-cas'; do
    parts=(${=pair})
    python3 "$measure" /opt/homebrew/bin/gtimeout -k 5s 60s \
        /Applications/Wolfram.app/Contents/MacOS/wolfram -script "scripts/cas/$parts[1].wls" \
        > "$out/$parts[2].log" 2>&1
    # Kernel parse errors can exit zero: require the actual closed witness.
    rg -q '^PROVED' "$out/$parts[2].log"
    if rg -q 'NOT-CLOSED|Syntax::|Check::' "$out/$parts[2].log"; then exit 1; fi
done
for spec in 'rn RN-identity 1 2 4' 'ems a0.9-e7 0.0512 0.2622 0.6069 1.0429'; do
    parts=(${=spec})
    python3 "$measure" julia --project=/Users/auroradysis/Workspace/EMS \
        Tests/EMSNP/static_oracle.jl "/Users/auroradysis/Workspace/EMS/artifacts/echo/$parts[2].trumpet" \
        "$out/$parts[1]-oracle.csv" $parts[3,-1] > "$out/$parts[1]-oracle.log" 2>&1
done
for case in rn-n256 rn-n512 ems-n1024 ems-n2048 checkpoint cavity default; do
    exe=EMSNPStatic
    if [[ "$case" == checkpoint || "$case" == cavity || "$case" == default ]]; then exe=EMSWeylCheckpoint; fi
    (
        cd "$out/$case"
        python3 "$measure" "../../../EMSWeylOffline/${exe}2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex" params.txt
    ) > "$out/$case-release.log" 2>&1
done
dependency_status=0
env OMP_NUM_THREADS=1 Tests/EMSWeylOffline/EMSNPStatic2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex \
    "$out/bad-activation.params.txt" > "$out/activation-dependency.log" 2>&1 || dependency_status=$?
[[ "$dependency_status" == 2 ]]
rg -q 'ems_np_activate requires ems_radiation_activate' "$out/activation-dependency.log"
python3 Tests/EMSNP/prepare_t3.py --verify-only > "$out/input-final.log"
report_status=0
python3 Tests/EMSNP/report_t3.py > "$out/report.log" 2>&1 || report_status=$?
cat "$out/report.log"
print "RUN_COMPLETED validation_exit=$report_status"
exit $report_status
