/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
// Problem specific includes (robin):
#include "EMSBHParams.hpp"
#include "EMSCouplingFunction.hpp"

#ifdef USE_AHFINDER
#include "AHInitialGuess.hpp"
#endif

class SimulationParameters : public SimulationParametersBase
{
public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_threshold_A", regrid_threshold_A);
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        pp.load("G_Newton", m_G_Newton);

        // EMSBH initial data params
        pp.load("ems_data_path", emsbh_params.data_path); // the default should fail
        pp.load("ems_data_format", ems_data_format, std::string("legacy_dat"));
        pp.load("t2_guard_initial_data_after_t0", t2_guard_initial_data_after_t0,
                false);
        pp.load("ems_ctt_data_path", emsbh_params.ctt_data_path, std::string());
        pp.load("gridpoints", emsbh_params.gridpoints, 40000);
        pp.load("star_centre", emsbh_params.star_centre,
                {0.5 * L, 0.5 * L});
        pp.load("binary", emsbh_params.binary, false);
        pp.load("boosted", emsbh_params.boosted, false);
        pp.load("boost_rapidity", emsbh_params.rapidity, 0.0);
        pp.load("separation", emsbh_params.separation, 0.0);
        pp.load("bh_charge", emsbh_params.bh_charge, 0.0);
        pp.load("bh_mass", emsbh_params.bh_mass, 1.0);
        pp.load("G_Newton", emsbh_params.Newtons_constant, 1.0);

        // EMS vs RB
        pp.load("ems_not_rn", EMS_not_RN, false);
        // Coupling params
        pp.load("ems_alpha", coupling_function_params.alpha, 0.0);
        pp.load("ems_f0", coupling_function_params.f0, 0.0);
        pp.load("ems_f1", coupling_function_params.f1, 0.0);
        pp.load("ems_f2", coupling_function_params.f2, 0.0);

        // Perturbation shell params
        pp.load("Ylm_amplitude", emsbh_params.Ylm_amplitude, 0.);
        pp.load("Ylm_thickness", emsbh_params.Ylm_thickness, 3.);
        pp.load("Ylm_r0", emsbh_params.Ylm_r0 , 0.5 * L);


        // Do we want Weyl extraction, puncture tracking and constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.5 * L, 0.0});


        // Apparent Horizon stuff
        #ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5); // ~M/2
        pp.load("AH_num_horizons", AH_num_horizons, 0);
        pp.load("AH_level_to_run", AH_level_to_run, 0);
        pp.load("AH_expect_merger", AH_expect_merger, 0);
        pp.load("horizon_centre_1", horizon_centre_1,
                {0.5 * L, 0.0});
        pp.load("horizon_centre_2", horizon_centre_2,
                {0.5 * L, 0.0});
        #endif

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.5 * L, 0.0});

        // Weyl extraction
        pp.load("activate_gw_extraction", activate_weyl_extraction, 0);

        // Em extraction stuff
        pp.load("activate_em_extraction", activate_pheyl_extraction, 0);

        // real scalar extraction stuff
        pp.load("activate_rs_extraction", activate_rs_extraction, 0);

        // mass charge scalar extraction stuff
        pp.load("activate_mq_extraction", activate_mq_extraction, 0);


        // RH horizon finder
        pp.load("RH_activate", m_RH_activate, false);
        pp.load("RH_num_horizons", m_RH_num_horizons, 0);
        pp.load("RH_initial_radii",    m_RH_initial_radii,    m_RH_num_horizons, 1.0);
        pp.load("RH_initial_centre",   m_RH_initial_centre,   m_RH_num_horizons, 0.0);
        pp.load("RH_num_points",       m_RH_num_points,       m_RH_num_horizons, 100);
        pp.load("RH_level",            m_RH_level,            m_RH_num_horizons, 0);
        pp.load("RH_time_step_freq",   m_RH_time_step_freq,   m_RH_num_horizons, 1);
        pp.load("RH_newton_crit",      m_RH_newton_crit,      m_RH_num_horizons, 0.0);
        pp.load("RH_chase_speeds",     m_RH_chase_speeds,     m_RH_num_horizons, 1.0);
        pp.load("RH_start_times",      m_RH_start_times,      m_RH_num_horizons, 0.0);

        // Variables for outputting inf-norm
        pp.load("num_vars_inf_norm", num_vars_inf_norm, 0);
        pp.load("vars_inf_norm", vars_inf_norm, num_vars_inf_norm, 0);


        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);
    }

    void check_params()
    {
        if (ems_data_format != "legacy_dat" && ems_data_format != "emstrumpet1")
            MayDay::Error("unknown ems_data_format; use legacy_dat or emstrumpet1");
        if (ems_data_format == "emstrumpet1" && !EMS_not_RN)
            MayDay::Error("emstrumpet1 requires ems_not_rn = true");
        if (!emsbh_params.ctt_data_path.empty() &&
            (ems_data_format != "emstrumpet1" || !emsbh_params.binary))
            MayDay::Error("ems_ctt_data_path requires an emstrumpet1 binary");
    }

    // tagging
    bool activate_extraction;

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi, regrid_threshold_A;

    // extraction stuff
    int activate_weyl_extraction;
    int activate_pheyl_extraction;
    int activate_rs_extraction;
    int activate_mq_extraction;

    // // Layer regridding
    // bool m_do_layer_tagging;
    // double y_regrid_lim;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // EMS BH stuff
    EMSBH_params_t emsbh_params;
    std::string ems_data_format;
    bool t2_guard_initial_data_after_t0;
    CouplingFunction::params_t coupling_function_params;
    bool EMS_not_RN;

    double m_G_Newton;
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    // Vars for outputting inf-norms
    int num_vars_inf_norm;
    std::vector<int> vars_inf_norm;

    // RH horizon finder
    bool m_RH_activate;
    int m_RH_num_horizons;
    std::vector<double> m_RH_initial_radii;
    std::vector<double> m_RH_initial_centre; // x-coord per surface (y=0 by cartoon symmetry)
    std::vector<int>    m_RH_num_points;
    std::vector<int>    m_RH_level;
    std::vector<int>    m_RH_time_step_freq;
    std::vector<double> m_RH_newton_crit;
    std::vector<double> m_RH_chase_speeds;
    std::vector<double> m_RH_start_times;

#ifdef USE_AHFINDER
    double AH_initial_guess;
    int AH_num_horizons;
    int AH_expect_merger;
    int AH_level_to_run;
    std::array<double, CH_SPACEDIM> horizon_centre_1;
    std::array<double, CH_SPACEDIM> horizon_centre_2;
#endif

};
#endif /* SIMULATIONPARAMETERS_HPP_ */
