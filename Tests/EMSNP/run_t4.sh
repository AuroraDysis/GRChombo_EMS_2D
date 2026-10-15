#!/bin/zsh
set -eu
cd "${0:A:h}/../.."
source Tests/EMSNP/config.sh
out=Tests/EMSNP/results-t4
logs=$EMS_NP_LOG_ROOT/t4
mkdir -p "$logs"
mkdir -p "$out"
measure="$PWD/Tests/EMSNP/measure_t3.py"
python3 Tests/EMSNP/prepare_t4.py > "$logs/prepare.log"
python3 "$measure" make -C Tests/EMSWeylOffline all DIM=2 CXX="$CXX" -j 2 > "$logs/build.log" 2>&1
python3 Tests/EMSNP/syntax_check.py "$out" > "$logs/syntax.log" 2>&1
exe="$PWD/Tests/EMSWeylOffline/EMSWeylCheckpoint${EMS_EXE_SUFFIX}"
negative_status=0
(cd "$out/insufficient-ghosts"; python3 "$measure" "$exe" params.txt) > "$logs/insufficient-ghosts.log" 2>&1 || negative_status=$?
[[ "$negative_status" == 2 ]]
rg -q 'interpolation stencil outside FAB' "$logs/insufficient-ghosts.log"
(cd "$out/study"; python3 "$measure" "$exe" params.txt) > "$logs/study.log" 2>&1
python3 -c 'import sys,json;sys.path.insert(0,"Tests/EMSNP");from prepare_t4 import authenticate,OUT;(OUT/"inputs-final.json").write_text(json.dumps(authenticate(),indent=2)+"\n");print("PASS: final input authentication")' > "$logs/input-final.log"
python3 Tests/EMSNP/report_t4.py > "$logs/report.log"
cat "$logs/report.log"
