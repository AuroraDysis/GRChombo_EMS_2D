#!/bin/zsh
# Detached continuation of the two registered high-resolution controls.
root=/Users/auroradysis/Workspace/EMS-deps/worktrees/wt-native-t1
runs=/private/tmp/ems-t3-ref
exe=$root/Examples/EMS/Main_EMSBH2DBH2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex
for flag in off on; do
    name=t3-ref-match-high-$flag
    dir=$runs/$name
    mkdir -p "$dir/plt" "$dir/chk"
    cd "$dir" || exit 98
    OMP_NUM_THREADS=8 "$exe" "$root/Tests/EMSNative/params/$name.txt" > "$dir/run.log" 2>&1
    code=$?
    print -r -- "$code" > "$dir/done.exit.tmp"
    mv "$dir/done.exit.tmp" "$dir/done.exit"
    if (( code != 0 )); then exit "$code"; fi
done
