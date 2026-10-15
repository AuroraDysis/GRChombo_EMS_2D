#!/bin/zsh
set -eu
cd "${0:A:h}/../.."
export CHOMBO_HOME=/Users/auroradysis/Workspace/EMS-deps/Chombo/lib
export OMP_NUM_THREADS=8
out=Tests/EMSNP/results-t4
mkdir -p "$out"
measure="$PWD/Tests/EMSNP/measure_t3.py"
python3 Tests/EMSNP/prepare_t4.py > "$out/prepare.log"
python3 "$measure" make -C Tests/EMSWeylOffline all DIM=2 -j 2 > "$out/build.log" 2>&1
python3 Tests/EMSNP/syntax_check.py "$out" > "$out/syntax.log" 2>&1
exe="$PWD/Tests/EMSWeylOffline/EMSWeylCheckpoint2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex"
negative_status=0
(cd "$out/insufficient-ghosts"; python3 "$measure" "$exe" params.txt) > "$out/insufficient-ghosts.log" 2>&1 || negative_status=$?
[[ "$negative_status" == 2 ]]
rg -q 'interpolation stencil outside FAB' "$out/insufficient-ghosts.log"
(cd "$out/study"; python3 "$measure" "$exe" params.txt) > "$out/study.log" 2>&1
python3 -c 'import sys,json;sys.path.insert(0,"Tests/EMSNP");from prepare_t4 import authenticate,OUT;(OUT/"inputs-final.json").write_text(json.dumps(authenticate(),indent=2)+"\n");print("PASS: final input authentication")' > "$out/input-final.log"
python3 Tests/EMSNP/report_t4.py > "$out/report.log"
cat "$out/report.log"
