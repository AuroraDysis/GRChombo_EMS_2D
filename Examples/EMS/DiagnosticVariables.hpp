/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_mod_F,
    c_phi_rad, // radiated scalar field
    c_Mscalar,
    c_Qscalar,

    c_Ham,

    c_Mom1,
    c_Mom2,

    c_Mom,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_Pheyl2_Re,
    c_Pheyl2_Im,

    c_rho, // stress tensor components
    c_rho_ADM,

    c_Sx,
    c_Sy,

    c_GaussE,
    c_GaussB,

    // c_Sxx,
    // c_Sxy,
    // c_Syy,

    // c_Sww,

    // c_S,

    // c_Sxx_TF,
    // c_Sxy_TF,
    // c_Syy_TF,
    // c_Sww_TF,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {

    "mod_F",
    "phi_rad",
    "Mscalar",
    "Qscalar",

    "Ham",
    "Mom1",
    "Mom2",
    "Mom",

    "Weyl4_Re",
    "Weyl4_Im",

    "Pheyl2_Re",
    "Pheyl2_Im",

    "rho",
    "rho_ADM",

    "Sx",
    "Sy",
    "GaussE",
    "GaussB"

  };
}

#endif /* DIAGNOSTICVARIABLES_HPP */
