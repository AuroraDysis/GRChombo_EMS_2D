# Portable serial offline tool. Invoke through t25-build-node.sh so recursive
# Chombo make finds this file as GNUmakefile in an isolated build directory.
# Override SOURCE_ROOT to a clean fork snapshot; no initial-data setter is linked.
SOURCE_ROOT ?= ../..
.DEFAULT_GOAL := all
GRCHOMBO_SOURCE = $(SOURCE_ROOT)/Source
ebase := t25-expansion-map
# Make.test would otherwise compile every unrelated T1--T24 main in this directory.
# It scans this temporary build directory, which contains only the T25 main.
base_dir := .
LibNames := AMRTimeDependent AMRTools BoxTools
src_dirs := $(GRCHOMBO_SOURCE)/utils $(GRCHOMBO_SOURCE)/simd \
    $(GRCHOMBO_SOURCE)/BoxUtils $(GRCHOMBO_SOURCE)/CCZ4 \
    $(GRCHOMBO_SOURCE)/GRChomboCore $(GRCHOMBO_SOURCE)/AMRInterpolator
include $(CHOMBO_HOME)/mk/Make.test
CPPFLAGS += -I$(SOURCE_ROOT)/Examples/EMS -I$(GRCHOMBO_SOURCE)/RHFinder \
    -I$(GRCHOMBO_SOURCE)/Matter -I$(GRCHOMBO_SOURCE)/Cartoon \
    -I$(GRCHOMBO_SOURCE)/InitialConditions/EMSBH \
    -I$(GRCHOMBO_SOURCE)/InitialConditions/EMSBH/1D_SOL
CXXFLAGS += -fno-fast-math -ffp-contract=off
# Native coverage ownership is exposed in the Tests main only.
o/$(config)/t25-expansion-map.o: XTRACXXFLAGS += -fno-access-control
vpath %.cpp $(base_dir) $(src_dirs)
