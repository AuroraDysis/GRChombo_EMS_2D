# Machine defaults for the local serial harness; callers may override each value.
: ${EMS_ROOT:=/Users/auroradysis/Workspace/EMS}
: ${CHOMBO_HOME:=${EMS_ROOT:h}/EMS-deps/Chombo/lib}
: ${EMS_NP_LOG_ROOT:=$EMS_ROOT/.data/fork-np-logs}
: ${CXX:=/opt/homebrew/bin/g++-16}
: ${GTIMEOUT:=/opt/homebrew/bin/gtimeout}
: ${WOLFRAM:=/Applications/Wolfram.app/Contents/MacOS/wolfram}
: ${JULIA:=julia}
: ${EMS_EXE_SUFFIX:=2d.Darwin.64.g++-16.gfortran.OPTHIGH.OPENMPCC.ex}
: ${OMP_NUM_THREADS:=8}
export EMS_ROOT CHOMBO_HOME EMS_NP_LOG_ROOT CXX OMP_NUM_THREADS
