# GRChombo
# Copyright 2012 The GRChombo collaboration.
# Please refer to LICENSE in GRChombo's root directory.

# -*- Mode: Makefile -*-
# GRChombo EMSTRUMPET 1 fixture and single-box tests.

ebase := EMSTrumpetFixtureTest EMSTrumpetGridConvergence EMSTrumpetBinaryGridConvergence
LibNames := BoxTools

include $(CHOMBO_HOME)/mk/Make.test

CPPFLAGS += -I../../Examples/EMS \
    -I../../Source/utils -I../../Source/simd \
    -I../../Source/CCZ4 -I../../Source/Matter \
    -I../../Source/Cartoon -I../../Source/BoxUtils \
    -I../../Source/GRChomboCore \
    -I../../Source/InitialConditions/EMSBH \
    -I../../Source/InitialConditions/EMSBH/1D_SOL
