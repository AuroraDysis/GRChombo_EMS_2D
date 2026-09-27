# Isolated t=0 and frozen checkpoint hierarchies, with production interpolation.
GRCHOMBO_SOURCE = ../../Source
ebase := EMSRHSurfaceTest EMSRHRestartTest EMSRHCheckpoint
base_dir := harness
LibNames := AMRTimeDependent AMRTools BoxTools
src_dirs := $(GRCHOMBO_SOURCE)/utils $(GRCHOMBO_SOURCE)/simd \
    $(GRCHOMBO_SOURCE)/BoxUtils $(GRCHOMBO_SOURCE)/CCZ4 \
    $(GRCHOMBO_SOURCE)/GRChomboCore $(GRCHOMBO_SOURCE)/AMRInterpolator
include $(CHOMBO_HOME)/mk/Make.test
CPPFLAGS += -I../../Examples/EMS \
    -I../../Source/RHFinder -I../../Source/Matter -I../../Source/Cartoon \
    -I../../Source/InitialConditions/EMSBH -I../../Source/InitialConditions/EMSBH/1D_SOL

# Make.test scans base_dir but its implicit compiler rule also needs this path.
vpath %.cpp $(base_dir) $(src_dirs)
