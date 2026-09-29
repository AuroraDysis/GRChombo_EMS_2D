/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERSBASE_HPP_
#define SIMULATIONPARAMETERSBASE_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "CCZ4RHS.hpp"
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
#include <limits>

#ifdef USE_AHFINDER
#include "AHFinder.hpp"
#endif

#include "SphericalExtraction.hpp"
// add this type alias here for backwards compatibility
using extraction_params_t = spherical_extraction_params_t;

class SimulationParametersBase : public ChomboParameters
{
  public:
    SimulationParametersBase(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
        check_params();
    }

  private:
    void read_params(GRParmParse &pp)
    {
        // Lapse evolution
        pp.load("lapse_advec_coeff", ccz4_params.lapse_advec_coeff, 1.0);
        pp.load("lapse_coeff", ccz4_params.lapse_coeff, 2.0);
        pp.load("lapse_power", ccz4_params.lapse_power, 1.0);

        // Shift Evolution
        pp.load("shift_advec_coeff", ccz4_params.shift_advec_coeff, 0.0);
        pp.load("shift_Gamma_coeff", ccz4_params.shift_Gamma_coeff, 0.75);
        pp.load("eta", ccz4_params.eta, 1.0);

        // CCZ4 parameters
        pp.load("formulation", formulation, 0);
        pp.load("kappa1", ccz4_base_params.kappa1, 0.1);
        pp.load("kappa2", ccz4_base_params.kappa2, 0.0);
        pp.load("kappa3", ccz4_base_params.kappa3, 1.0);
        pp.load("covariantZ4", ccz4_base_params.covariantZ4, true);
        ccz4_params.kappa1 = ccz4_base_params.kappa1;
        ccz4_params.kappa2 = ccz4_base_params.kappa2;
        ccz4_params.kappa3 = ccz4_base_params.kappa3;
        ccz4_params.covariantZ4 = ccz4_base_params.covariantZ4;

        // Dissipation
        pp.load("sigma", sigma, 0.1);
        pp.load("amr_transfer", amr_transfer, std::string("legacy"));
        if (amr_transfer != "legacy" && amr_transfer != "point")
            MayDay::Error("amr_transfer must be legacy or point");

        // Nan Check and min chi and lapse values
        pp.load("nan_check", nan_check, true);
        pp.load("min_chi", min_chi, 1e-4);
        pp.load("min_lapse", min_lapse, 1e-4);

        // directory to store data (extraction files, puncture data, constraint
        // norms)
        pp.load("data_subpath", data_path, std::string(""));
        if (!data_path.empty() && data_path.back() != '/')
            data_path += "/";
        if (output_path != "./" && !output_path.empty())
            data_path = output_path + data_path;

        // Extraction params
        pp.load("activate_extraction", activate_extraction, false);

        if (activate_extraction)
        {
            pp.load("num_extraction_radii",
                    extraction_params.num_extraction_radii, 1);

            // Check for multiple extraction radii, otherwise load single
            // radius/level (for backwards compatibility).
            if (pp.contains("extraction_levels"))
            {
                pp.load("extraction_levels",
                        extraction_params.extraction_levels,
                        extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("extraction_level", extraction_params.extraction_levels,
                        1, 0);
            }
            if (pp.contains("extraction_radii"))
            {
                pp.load("extraction_radii", extraction_params.extraction_radii,
                        extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("extraction_radius", extraction_params.extraction_radii,
                        1, 0.1);
            }

            pp.load("num_points_phi", extraction_params.num_points_phi, 2);
            pp.load("num_points_theta", extraction_params.num_points_theta, 5);
            if (extraction_params.num_points_theta % 2 == 0)
            {
                extraction_params.num_points_theta += 1;
                pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
                       << "rule so increased by 1.\n";
            }
            pp.load("extraction_center", extraction_params.center, center);

            if (pp.contains("modes"))
            {
                pp.load("num_modes", extraction_params.num_modes);
                std::vector<int> extraction_modes_vect(
                    2 * extraction_params.num_modes);
                pp.load("modes", extraction_modes_vect,
                        2 * extraction_params.num_modes);
                extraction_params.modes.resize(extraction_params.num_modes);
                for (int i = 0; i < extraction_params.num_modes; ++i)
                {
                    extraction_params.modes[i].first =
                        extraction_modes_vect[2 * i];
                    extraction_params.modes[i].second =
                        extraction_modes_vect[2 * i + 1];
                }
            }
            else
            {
                // by default extraction (l,m) = (2,0), (2,1) and (2,2)
                extraction_params.num_modes = 3;
                extraction_params.modes.resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    extraction_params.modes[i].first = 2;
                    extraction_params.modes[i].second = i;
                }
            }

            pp.load("write_extraction", extraction_params.write_extraction,
                    false);

            std::string extraction_path;
            if (pp.contains("extraction_subpath"))
            {
                pp.load("extraction_subpath", extraction_path);
                if (!extraction_path.empty() && extraction_path.back() != '/')
                    extraction_path += "/";
                if (output_path != "./" && !output_path.empty())
                    extraction_path = output_path + extraction_path;
            }
            else
                extraction_path = data_path;

            extraction_params.data_path = data_path;
            extraction_params.extraction_path = extraction_path;

            // default names to Weyl extraction
            pp.load("extraction_file_prefix",
                    extraction_params.extraction_file_prefix,
                    std::string("Weyl4_extraction_"));
            pp.load("integral_file_prefix",
                    extraction_params.integral_file_prefix,
                    std::string("Weyl4_mode_"));
        }

        ////////////
        // PHEYL2 PARAMS
        ////////////

        pp.load("activate_em_extraction", activate_em_extraction, false);

        if (activate_em_extraction)
        {
            pp.load("em_num_extraction_radii",
                    pheyl2_extraction_params.num_extraction_radii, 1);

            // Check for multiple extraction radii, otherwise load single
            // radius/level (for backwards compatibility).
            if (pp.contains("em_extraction_levels"))
            {
                pp.load("em_extraction_levels",
                        pheyl2_extraction_params.extraction_levels,
                        pheyl2_extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("em_extraction_level", pheyl2_extraction_params.extraction_levels,
                        1, 0);
            }
            if (pp.contains("em_extraction_radii"))
            {
                pp.load("em_extraction_radii", pheyl2_extraction_params.extraction_radii,
                        pheyl2_extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("em_extraction_radius", pheyl2_extraction_params.extraction_radii,
                        1, 0.1);
            }

            pp.load("em_num_points_phi", pheyl2_extraction_params.num_points_phi, 2);
            pp.load("em_num_points_theta", pheyl2_extraction_params.num_points_theta, 5);
            if (pheyl2_extraction_params.num_points_theta % 2 == 0)
            {
                pheyl2_extraction_params.num_points_theta += 1;
                pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
                       << "rule so increased by 1.\n";
            }
            pp.load("em_extraction_center", pheyl2_extraction_params.center, center);

            if (pp.contains("em_modes"))
            {
                pp.load("em_num_modes", pheyl2_extraction_params.num_modes);
                std::vector<int> pheyl2_extraction_modes_vect(
                    2 * pheyl2_extraction_params.num_modes);
                pp.load("em_modes", pheyl2_extraction_modes_vect,
                        2 * pheyl2_extraction_params.num_modes);
                pheyl2_extraction_params.modes.resize(pheyl2_extraction_params.num_modes);
                for (int i = 0; i < pheyl2_extraction_params.num_modes; ++i)
                {
                    pheyl2_extraction_params.modes[i].first =
                        pheyl2_extraction_modes_vect[2 * i];
                    pheyl2_extraction_params.modes[i].second =
                        pheyl2_extraction_modes_vect[2 * i + 1];
                }
            }
            else
            {
                // by default extraction (l,m) = (2,0), (2,1) and (2,2)
                pheyl2_extraction_params.num_modes = 3;
                pheyl2_extraction_params.modes.resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    pheyl2_extraction_params.modes[i].first = 2;
                    pheyl2_extraction_params.modes[i].second = i;
                }
            }

            pp.load("em_write_extraction", pheyl2_extraction_params.write_extraction,
                    false);

            std::string pheyl2_extraction_path;
            if (pp.contains("em_extraction_subpath"))
            {
                pp.load("em_extraction_subpath", pheyl2_extraction_path);
                if (!pheyl2_extraction_path.empty() && pheyl2_extraction_path.back() != '/')
                    pheyl2_extraction_path += "/";
                if (output_path != "./" && !output_path.empty())
                    pheyl2_extraction_path = output_path + pheyl2_extraction_path;
            }
            else
                pheyl2_extraction_path = data_path;

            pheyl2_extraction_params.data_path = data_path;
            pheyl2_extraction_params.extraction_path = pheyl2_extraction_path;

            // default names to Weyl extraction
            pp.load("em_extraction_file_prefix",
                    pheyl2_extraction_params.extraction_file_prefix,
                    std::string("Pheyl2_extraction_"));
            pp.load("em_integral_file_prefix",
                    pheyl2_extraction_params.integral_file_prefix,
                    std::string("Pheyl2_mode_"));
        }



        ////////////
        // REALSCALAR PARAMS
        ////////////

        pp.load("activate_rs_extraction", activate_rs_extraction, false);

        if (activate_rs_extraction)
        {
            pp.load("rs_num_extraction_radii",
                    rs_extraction_params.num_extraction_radii, 1);

            // Check for multiple extraction radii, otherwise load single
            // radius/level (for backwards compatibility).
            if (pp.contains("rs_extraction_levels"))
            {
                pp.load("rs_extraction_levels",
                        rs_extraction_params.extraction_levels,
                        rs_extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("rs_extraction_level", rs_extraction_params.extraction_levels,
                        1, 0);
            }
            if (pp.contains("rs_extraction_radii"))
            {
                pp.load("rs_extraction_radii", rs_extraction_params.extraction_radii,
                        rs_extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("rs_extraction_radius", rs_extraction_params.extraction_radii,
                        1, 0.1);
            }

            pp.load("rs_num_points_phi", rs_extraction_params.num_points_phi, 2);
            pp.load("rs_num_points_theta", rs_extraction_params.num_points_theta, 5);
            if (rs_extraction_params.num_points_theta % 2 == 0)
            {
                rs_extraction_params.num_points_theta += 1;
                pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
                       << "rule so increased by 1.\n";
            }
            pp.load("rs_extraction_center", rs_extraction_params.center, center);

            if (pp.contains("rs_modes"))
            {
                pp.load("rs_num_modes", rs_extraction_params.num_modes);
                std::vector<int> rs_extraction_modes_vect(
                    2 * rs_extraction_params.num_modes);
                pp.load("rs_modes", rs_extraction_modes_vect,
                        2 * rs_extraction_params.num_modes);
                rs_extraction_params.modes.resize(rs_extraction_params.num_modes);
                for (int i = 0; i < rs_extraction_params.num_modes; ++i)
                {
                    rs_extraction_params.modes[i].first =
                        rs_extraction_modes_vect[2 * i];
                    rs_extraction_params.modes[i].second =
                        rs_extraction_modes_vect[2 * i + 1];
                }
            }
            else
            {
                // by default extraction (l,m) = (2,0), (2,1) and (2,2)
                rs_extraction_params.num_modes = 3;
                rs_extraction_params.modes.resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    rs_extraction_params.modes[i].first = 2;
                    rs_extraction_params.modes[i].second = i;
                }
            }

            pp.load("rs_write_extraction", rs_extraction_params.write_extraction,
                    false);

            std::string rs_extraction_path;
            if (pp.contains("rs_extraction_subpath"))
            {
                pp.load("rs_extraction_subpath", rs_extraction_path);
                if (!rs_extraction_path.empty() && rs_extraction_path.back() != '/')
                    rs_extraction_path += "/";
                if (output_path != "./" && !output_path.empty())
                    rs_extraction_path = output_path + rs_extraction_path;
            }
            else
                rs_extraction_path = data_path;

            rs_extraction_params.data_path = data_path;
            rs_extraction_params.extraction_path = rs_extraction_path;

            // default names to Weyl extraction
            pp.load("rs_extraction_file_prefix",
                    rs_extraction_params.extraction_file_prefix,
                    std::string("RealScalar_extraction_"));
            pp.load("rs_integral_file_prefix",
                    rs_extraction_params.integral_file_prefix,
                    std::string("RealScalar_mode_"));
        }



        ////////////
        // MASSCHARGE PARAMS
        ////////////

        pp.load("activate_mq_extraction", activate_mq_extraction, false);

        if (activate_mq_extraction)
        {
            pp.load("mq_num_extraction_radii",
                    mq_extraction_params.num_extraction_radii, 1);

            // Check for multiple extraction radii, otherwise load single
            // radius/level (for backwards compatibility).
            if (pp.contains("mq_extraction_levels"))
            {
                pp.load("mq_extraction_levels",
                        mq_extraction_params.extraction_levels,
                        mq_extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("mq_extraction_level", mq_extraction_params.extraction_levels,
                        1, 0);
            }
            if (pp.contains("mq_extraction_radii"))
            {
                pp.load("mq_extraction_radii", mq_extraction_params.extraction_radii,
                        mq_extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("mq_extraction_radius", mq_extraction_params.extraction_radii,
                        1, 0.1);
            }

            pp.load("mq_num_points_phi", mq_extraction_params.num_points_phi, 2);
            pp.load("mq_num_points_theta", mq_extraction_params.num_points_theta, 5);
            if (mq_extraction_params.num_points_theta % 2 == 0)
            {
                mq_extraction_params.num_points_theta += 1;
                pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
                       << "rule so increased by 1.\n";
            }
            pp.load("mq_extraction_center", mq_extraction_params.center, center);


            if (pp.contains("mq_modes"))
            {
                pp.load("mq_num_modes", mq_extraction_params.num_modes);
                std::vector<int> mq_extraction_modes_vect(
                    2 * mq_extraction_params.num_modes);
                pp.load("mq_modes", mq_extraction_modes_vect,
                        2 * mq_extraction_params.num_modes);
                mq_extraction_params.modes.resize(mq_extraction_params.num_modes);
                for (int i = 0; i < mq_extraction_params.num_modes; ++i)
                {
                    mq_extraction_params.modes[i].first =
                        mq_extraction_modes_vect[2 * i];
                    mq_extraction_params.modes[i].second =
                        mq_extraction_modes_vect[2 * i + 1];
                }
            }
            else
            {
                // by default extraction (l,m) = (2,0), (2,1) and (2,2)
                mq_extraction_params.num_modes = 3;
                mq_extraction_params.modes.resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    mq_extraction_params.modes[i].first = 2;
                    mq_extraction_params.modes[i].second = i;
                }
            }

            pp.load("mq_write_extraction", mq_extraction_params.write_extraction,
                    false);

            std::string mq_extraction_path;
            if (pp.contains("mq_extraction_subpath"))
            {
                pp.load("mq_extraction_subpath", mq_extraction_path);
                if (!mq_extraction_path.empty() && mq_extraction_path.back() != '/')
                    mq_extraction_path += "/";
                if (output_path != "./" && !output_path.empty())
                    mq_extraction_path = output_path + mq_extraction_path;
            }
            else
                mq_extraction_path = data_path;

            mq_extraction_params.data_path = data_path;
            mq_extraction_params.extraction_path = mq_extraction_path;

            // default names to mq extraction
            pp.load("mq_extraction_file_prefix",
                    mq_extraction_params.extraction_file_prefix,
                    std::string("MQ_extraction_"));
            pp.load("mq_integral_file_prefix",
                    mq_extraction_params.integral_file_prefix,
                    std::string("MQ_mode_"));
        }

#ifdef USE_AHFINDER
        // Apparent horizon parameters
        pp.load("AH_activate", AH_activate, false);
        if (AH_activate)
            AH_params.read_params(pp, *this);
#endif
    }

    void check_params()
    {
        check_parameter("dt_multiplier", dt_multiplier, dt_multiplier < 1.0,
                        "must be < 1.0 for stability");
        warn_parameter("dt_multiplier", dt_multiplier, dt_multiplier <= 0.5,
                       "is unlikely to be stable for > 0.5");

        check_parameter("sigma", sigma,
                        (sigma >= 0.0) && (sigma <= 2.0 / dt_multiplier),
                        "must be >= 0.0 and <= 2 / dt_multiplier for stability "
                        "(see Alcubierre p344)");
        warn_parameter("nan_check", nan_check, nan_check,
                       "should not normally be disabled");
        // not sure these are necessary hence commented out
        // check_parameter("min_chi", min_chi, (min_chi >= 0.0), "must be >=
        // 0.0"); check_parameter("min_lapse", min_lapse, (min_lapse >= 0.0)
        // "must be >= 0.0");
        check_parameter("formulation", formulation,
                        (formulation == CCZ4RHS<>::USE_CCZ4) ||
                            (formulation == CCZ4RHS<>::USE_BSSN),
                        "must be 0 or 1");
        if (formulation == CCZ4RHS<>::USE_CCZ4)
        {
            warn_parameter(
                "kappa1", ccz4_params.kappa1, ccz4_params.kappa1 > 0.0,
                "should be greater than 0.0 to damp constraints (see "
                "arXiv:1106.2254).");
            warn_parameter("kappa2", ccz4_params.kappa2,
                           ccz4_params.kappa2 > -1.0,
                           "should be greater than -1.0 to damp constraints "
                           "(see arXiv:1106.2254)");
        }
        else if (formulation == CCZ4RHS<>::USE_BSSN)
        {
            // maybe we should just set these to zero and print a warning
            // in the BSSN case
            warn_parameter("kappa1", ccz4_params.kappa1,
                           ccz4_params.kappa1 == 0.0,
                           "setting to 0.0 as required for BSSN");
            warn_parameter("kappa2", ccz4_params.kappa2,
                           ccz4_params.kappa2 == 0.0,
                           "setting to 0.0 as required for BSSN");
            warn_parameter("kappa3", ccz4_params.kappa3,
                           ccz4_params.kappa3 == 0.0,
                           "setting to 0.0 as required for BSSN");
            // no warning necessary for ccz4_params.covariantZ4
            ccz4_params.kappa1 = 0.0;
            ccz4_params.kappa2 = 0.0;
            ccz4_params.kappa3 = 0.0;
        }

        // only warn for gauge parameters as there are legitimate cases you may
        // want to deviate from the norm
        warn_parameter("lapse_advec_coeff", ccz4_params.lapse_advec_coeff,
                       min(abs(ccz4_params.lapse_advec_coeff),
                           abs(ccz4_params.lapse_advec_coeff - 1.0)) <
                           std::numeric_limits<double>::epsilon(),
                       "usually set to 0.0 or 1.0");
        warn_parameter("lapse_power", ccz4_params.lapse_power,
                       abs(ccz4_params.lapse_power - 1.0) <
                           std::numeric_limits<double>::epsilon(),
                       "set to 1.0 for 1+log slicing");
        warn_parameter("lapse_coeff", ccz4_params.lapse_coeff,
                       abs(ccz4_params.lapse_coeff - 2.0) <
                           std::numeric_limits<double>::epsilon(),
                       "set to 2.0 for 1+log slicing");
        warn_parameter("shift_Gamma_coeff", ccz4_params.shift_Gamma_coeff,
                       abs(ccz4_params.shift_Gamma_coeff - 0.75) <
                           std::numeric_limits<double>::epsilon(),
                       "usually set to 0.75");
        warn_parameter("shift_advec_coeff", ccz4_params.shift_advec_coeff,
                       min(abs(ccz4_params.shift_advec_coeff),
                           abs(ccz4_params.shift_advec_coeff - 1.0)) <
                           std::numeric_limits<double>::epsilon(),
                       "usually set to 0.0 or 1.0");
        warn_parameter("eta", ccz4_params.eta,
                       ccz4_params.eta > 0.1 && ccz4_params.eta < 10,
                       "usually O(1/M_ADM) so typically O(1) in code units");

        // Now extraction parameters
        if (activate_extraction)
        {
            check_parameter(
                "num_extraction_radii", extraction_params.num_extraction_radii,
                extraction_params.num_extraction_radii > 0,
                "must be bigger than 0 when activate_extraction = 1");

            FOR(idir)
            {
                std::string center_name =
                    "extraction_center[" + std::to_string(idir) + "]";
                double center_in_dir = extraction_params.center[idir];
                check_parameter(
                    center_name, center_in_dir,
                    (center_in_dir >= reflective_domain_lo[idir]) &&
                        (center_in_dir <= reflective_domain_hi[idir]),
                    "must be in the computational domain after "
                    "applying reflective symmetry");
                for (int iradius = 0;
                     iradius < extraction_params.num_extraction_radii;
                     ++iradius)
                {
                    std::string radius_name =
                        "extraction_radii[" + std::to_string(iradius) + "]";
                    double radius = extraction_params.extraction_radii[iradius];
                    if (idir == 0)
                        check_parameter(radius_name, radius, radius >= 0.0,
                                        "must be >= 0.0");
                    check_parameter(
                        radius_name, radius,
                        (center_in_dir - radius >=
                         reflective_domain_lo[idir]) &&
                            (center_in_dir + radius <=
                             reflective_domain_hi[idir]),
                        "extraction sphere must lie within the computational "
                        "domain after applying reflective symmetry");
                }
            }
            for (int imode = 0; imode < extraction_params.num_modes; ++imode)
            {
                auto &mode = extraction_params.modes[imode];
                int l = mode.first;
                int m = mode.second;
                std::string mode_name = "modes[" + std::to_string(imode) + "]";
                std::string value_str = "(" + std::to_string(mode.first) +
                                        ", " + std::to_string(mode.second) +
                                        ")";
                check_parameter(
                    mode_name, value_str, (l >= 2) && (abs(m) <= l),
                    "l must be >= 2 and m must satisfy -l <= m <= l");
            }
        }
    }

  protected:
    // This is just the CCZ4 damping parameters in case you want to use
    // a different gauge (with different parameters)
    CCZ4_base_params_t ccz4_base_params;

  public:
    double sigma; // Kreiss-Oliger dissipation parameter
    std::string amr_transfer; // One opt-in switch for space AND RK-stage transfer.

    bool nan_check;

    double min_chi, min_lapse;

    int formulation; // Whether to use BSSN or CCZ4

    // Collection of parameters necessary for the CCZ4 RHS
    // Note the gauge parameters are specific to MovingPunctureGauge
    // If you are using a different gauge, you need to load your parameters
    // in your own SimulationParameters class.
    CCZ4_params_t<> ccz4_params;

    bool activate_extraction;
    bool activate_em_extraction;
    bool activate_rs_extraction;
    bool activate_mq_extraction;
    spherical_extraction_params_t extraction_params;
    spherical_extraction_params_t pheyl2_extraction_params;
    spherical_extraction_params_t rs_extraction_params;
    spherical_extraction_params_t mq_extraction_params;

#ifdef USE_AHFINDER
    bool AH_activate;
    AHParams_t<AHFunction> AH_params;
#endif

    std::string data_path;
};

#endif /* SIMULATIONPARAMETERSBASE_HPP_ */
