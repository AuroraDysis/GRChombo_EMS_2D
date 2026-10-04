#!/bin/sh
# Serial tool build, no MPI launch, cluster submission, or initial-data reads.
# Arguments: clean-fork-root Chombo-lib build-directory [make overrides ...]
set -eu
test "$#" -ge 3
source_root=$1;chombo_lib=$2;build_dir=$3;shift 3
packet_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_root=$(CDPATH= cd -- "$source_root" && pwd)
chombo_lib=$(CDPATH= cd -- "$chombo_lib" && pwd)
mkdir -p "$build_dir"
cp "$packet_dir/t25-expansion-map.cpp" "$build_dir/"
cp "$packet_dir/t25.make" "$build_dir/GNUmakefile"
exec make -C "$build_dir" all SOURCE_ROOT="$source_root" CHOMBO_HOME="$chombo_lib" \
    DIM=2 MPI=FALSE OPENMPCC=FALSE PRECISION=DOUBLE USE_64=TRUE USE_HDF=TRUE -j1 "$@"
