GRCHOMBO_SOURCE = ../../Source
ebase := EMSCTT2Readback EMSCTT2Hierarchy
LibNames := AMRTimeDependent AMRTools BoxTools
src_dirs := $(GRCHOMBO_SOURCE)/utils $(GRCHOMBO_SOURCE)/simd \
    $(GRCHOMBO_SOURCE)/BoxUtils $(GRCHOMBO_SOURCE)/CCZ4 \
    $(GRCHOMBO_SOURCE)/GRChomboCore $(GRCHOMBO_SOURCE)/AMRInterpolator
include $(CHOMBO_HOME)/mk/Make.test
CPPFLAGS += -I../../Examples/EMS \
    -I../../Source/utils -I../../Source/simd \
    -I../../Source/CCZ4 -I../../Source/Matter \
    -I../../Source/Cartoon -I../../Source/BoxUtils \
    -I../../Source/GRChomboCore \
    -I../../Source/InitialConditions/EMSBH \
    -I../../Source/InitialConditions/EMSBH/1D_SOL

CPPFLAGS += -I../../Source/RHFinder
vpath %.cpp $(src_dirs)
