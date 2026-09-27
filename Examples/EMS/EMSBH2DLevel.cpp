/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#include <fstream>
#include <mutex>
#include <memory>

#include "EMSBH2DLevel.hpp"

// General headers
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"

// For tagging cells
#include "EMSExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "CCZ4Cartoon.hpp"
#include "ConstraintsCartoon.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"

// gauge
// #include "IntegratedMovingPunctureGauge.hpp"
// #include "MovingPunctureGauge.hpp"
#include "ExperimentalGauge.hpp"

// EMS includes
#include "EMSBH_read.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "EMSCouplingFunction.hpp"
#include "EMSCartoonGaussConstraints.hpp"
// #include "FixSuperposition_K.hpp"
#include "FixSuperposition_metric.hpp"
// RNBH test include initial data
#include "RNBH_read.hpp"

// For GW extraction
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

// For MQ extraction
#include "CrudeMassChargeExtraction.hpp"

// FOR EM extraction_level
#include "Pheyl2.hpp"
#include "PheylExtraction.hpp"

// For Scalar Field integrator
#include "RealScalarExtraction.hpp"

// For analysis outputs
#include "EMSCartoonLorentzScalars.hpp"
#include "SmallDataIO.hpp"

// For EM Tensor in 3D
// #include "EMTensor.hpp"

// other
#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"
#include "GammaCartoonCalculator.hpp"

void EMSBH2DLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}


// Initial data for field and metric variables
void EMSBH2DLevel::initialData()
{
    CH_TIME("EMSBH2DLevel::initialData");
    if (m_verbosity)
        pout() << "EMSBH2DLevel::initialData " << m_level << endl;

    // Read initial data for RN or EMS
    if (m_p.ems_data_format == "emstrumpet1")
    {
        EMSBH_trumpet_read emdbh(m_p.emsbh_params, m_p.coupling_function_params,
                                 m_p.m_G_Newton, m_dx, m_verbosity);
        emdbh.compute_1d_solution();
        if (m_verbosity)
            pout() << "EMSBH2DLevel::initialData - Interpolate to 3D grid "
                   << m_level << endl;
        BoxLoops::loop(make_compute_pack(SetValue(0.0), emdbh), m_state_new,
                       m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
        fillAllGhosts();
        BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
    }
    else if (m_p.EMS_not_RN) {
        // EMS data
        EMSBH_read emdbh(m_p.emsbh_params, m_p.coupling_function_params,
                            m_p.m_G_Newton, m_dx, m_verbosity);
        emdbh.compute_1d_solution();

        if (m_verbosity)
            pout() << "EMSBH2DLevel::initialData - Interpolate to 3D grid " << m_level << endl;
        // First set everything to zero ... we don't want undefined values in
        // constraints etc, then  initial conditions for EMSBH
        BoxLoops::loop(make_compute_pack(SetValue(0.0), emdbh), m_state_new,
                       m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

        if (m_p.emsbh_params.boosted)
        {
            // BoxLoops::loop(FixSuperposition_K(m_dx, m_p.emsbh_params.star_centre, m_p.emsbh_params.rapidity, m_p.emsbh_params.binary),
            //              m_state_new, m_state_new,
            //              INCLUDE_GHOST_CELLS, disable_simd());
            BoxLoops::loop(FixSuperposition_metric(m_dx, m_p.emsbh_params.star_centre, m_p.emsbh_params.rapidity, m_p.emsbh_params.binary),
                         m_state_new, m_state_new,
                         INCLUDE_GHOST_CELLS, disable_simd());
        }

        if (m_verbosity)
            pout() << "EMSBH2DLevel::initialData - GammaCalc " << m_level << endl;

        fillAllGhosts();
        BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());


    }
    // else not ems, algebraic rn
    else {
        // RN data
        RNBH_read emdbh(m_p.emsbh_params, m_p.coupling_function_params,
                            m_p.m_G_Newton, m_dx, m_verbosity);
        emdbh.compute_1d_solution();

        if (m_verbosity)
            pout() << "EMSBH2DLevel::initialData - Interpolate to 3D grid " << m_level << endl;
        // First set everything to zero ... we don't want undefined values in
        // constraints etc, then  initial conditions for EMSBH
        BoxLoops::loop(make_compute_pack(SetValue(0.0), emdbh), m_state_new,
                       m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
        if (m_verbosity)
            pout() << "EMSBH2DLevel::initialData - GammaCalc " << m_level << endl;

        fillAllGhosts();
        BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
    }

    fillAllGhosts();

    // Exp/Int gauge ONLY! to initialise B^i
    // auto my_gauge_conditions = IntegratedMovingPunctureGauge(m_p.ccz4_params);
    auto my_gauge_conditions = ExperimentalGauge(m_p.ccz4_params);
    BoxLoops::loop(my_gauge_conditions,
                      m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
}



// Things to do before a plot level - need to calculate the Weyl scalars
void EMSBH2DLevel::prePlotLevel()
{
    fillAllGhosts();
    // coupling not currently used in Constraint calcs
    CouplingFunction my_coupling(m_p.coupling_function_params);
    BoxLoops::loop(
        make_compute_pack(

            WeylOmScalar(m_p.extraction_params.center, m_dx),

            Constraints<CouplingFunction>(m_dx, my_coupling, m_p.m_G_Newton),
            EMSCartoonGaussConstraints(m_dx, m_p.coupling_function_params),

            EMSCartoonLorentzScalars<CouplingFunction>(m_dx,
                                        m_p.mq_extraction_params.center,
                                             m_p.coupling_function_params)

            //EMTensor<EinsteinMaxwellScalarFieldWithCoupling>(
                //emd_field, m_dx, c_rho,
                //Interval(c_s1, c_s3), Interval(c_s11, c_s33)),

            //Pheyl2(m_p.extraction_params.extraction_center,
                              //m_p.coupling_function_params, m_dx)
                                                                 ),

      m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void EMSBH2DLevel::specificEvalRHS(GRLevelData &a_soln,
                                         GRLevelData &a_rhs,
                                         const double a_time)
{
    ////////////////////////////////////////
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);


    ///////////////////////////////////
    // Coupling Function
    CouplingFunction coupling_function(m_p.coupling_function_params);


    ////////////////////////////
    // Integrated MPG
    // CCZ4Cartoon<IntegratedMovingPunctureGauge,
    //             FourthOrderDerivatives,
    //             CouplingFunction>
    // my_ccz4_cartoon(m_p.ccz4_params, m_dx, m_p.sigma, coupling_function,
    //                                        m_p.m_G_Newton, m_p.formulation);


    // //////////////////////////////
    // // MPG
    // CCZ4Cartoon<MovingPunctureGauge,
    //             FourthOrderDerivatives,
    //             CouplingFunction>
    // my_ccz4_cartoon(m_p.ccz4_params, m_dx, m_p.sigma, coupling_function,
    //                                       m_p.m_G_Newton, m_p.formulation);

    // //////////////////////////////
    // // XPG
    CCZ4Cartoon<ExperimentalGauge,
                FourthOrderDerivatives,
                CouplingFunction>
    my_ccz4_cartoon(m_p.ccz4_params, m_dx, m_p.sigma, coupling_function,
                                          m_p.m_G_Newton, m_p.formulation);


    ///////////////////////
    // zero diagnostic vars
    SetValue set_analysis_vars_zero(0.0, Interval(c_Xi + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_cartoon, set_analysis_vars_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void EMSBH2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}



void EMSBH2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state)
{
    BoxLoops::loop(EMSExtractionTaggingCriterion(
                     m_dx, m_level, m_p.mass_extraction_params,
                     m_p.regrid_threshold_A, m_p.regrid_threshold_phi,
                     m_p.regrid_threshold_chi), current_state,
                                                tagging_criterion);

}

void EMSBH2DLevel::specificPostTimeStep()
{
    CH_TIME("EMSBH2DLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main


    fillAllGhosts();
    CouplingFunction my_coupling(m_p.coupling_function_params);

    BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(
            Constraints<CouplingFunction>(m_dx, my_coupling, m_p.m_G_Newton),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(EMSCartoonGaussConstraints(m_dx, m_p.coupling_function_params),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(
            EMSCartoonLorentzScalars<CouplingFunction>(m_dx,
                                       m_p.mq_extraction_params.center,
                                            m_p.coupling_function_params),
                      m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

        //////////////////////////////////////////////
        // Horizon finding (if used)
        //////////////////////////////////////////////
#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_level_to_run)
    {
        pout() << "Started AHFINDER!" << std::endl;
        // VERY IMPORTANT: refresh interpolator before horizon finding
        fillAllGhosts();
        m_bh_amr.m_interpolator->refresh();


        // // Hack: if avg radius of found AH is negative, we reset initial guess
        // double current_avg_AH_radius =  m_bh_amr.m_ah_finder.get(0)->get_ave_F();
        // if (current_avg_AH_radius < 0.)
        // {
        //     m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
        //     pout() << "AHFinder: resetting initial guess as avg radius is "
        //               "negative."
        //            << endl;
        // }
        // else if (current_avg_AH_radius > 13.)
        // {
        //     m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
        //     pout() << "AHFinder: resetting initial guess as avg radius is "
        //               "too large."
        //            << endl;
        // }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
        pout() << "Finished AHFINDER!" << std::endl;

    }
#endif


    //////////////////////////////////////////////
    // RH Horizon Finder
    //////////////////////////////////////////////
    if (m_p.m_RH_activate)
    {
        if (m_bh_amr.m_rh_union.m_surfaces.empty())
        {
            m_bh_amr.m_rh_union.setup(m_p.m_RH_num_horizons,
                                       m_p.m_RH_initial_radii,
                                       m_p.m_RH_initial_centre,
                                       m_p.m_RH_num_points,
                                       m_p.m_RH_level,
                                       m_p.m_RH_time_step_freq,
                                       m_p.m_RH_newton_crit,
                                       m_p.m_RH_chase_speeds,
                                       m_p.m_RH_start_times,
                                       m_restart_time);
            m_bh_amr.m_rh_union.set_coupling_params(
                m_p.coupling_function_params.alpha,
                m_p.coupling_function_params.f0,
                m_p.coupling_function_params.f1,
                m_p.coupling_function_params.f2);
        }
        m_bh_amr.m_rh_union.update(m_time, m_level);
    }

    //////////////////////////////////////////////
    // Gravitational Wave Extraction
    //////////////////////////////////////////////
    if (m_p.activate_extraction == 1 &&
       at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
            {
                pout() << "EMSBH2DLevel::specificPostTimeStep:"
                          " Extracting gravitational waves." << endl;
            }


            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    //////////////////////////////////////////////
    // Constraints file and ADM mass
    //////////////////////////////////////////////
    if (m_level == 0)
    {
        bool first_step = (m_time == 0.);
        AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom2), 2, true);
        SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom});

        int adm_min_level = 0;
        bool calculate_adm = at_level_timestep_multiple(adm_min_level);
        if (calculate_adm)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
            double M_ADM = amr_reductions.sum(c_rho_ADM);
            SmallDataIO M_ADM_file(m_p.data_path + "M_ADM", m_dt, m_time,
                                m_restart_time, SmallDataIO::APPEND,
                                first_step);
            M_ADM_file.remove_duplicate_time_data();
            if (first_step)
            {
                M_ADM_file.write_header_line({"M_ADM"});
            }
            M_ADM_file.write_time_data_line({M_ADM});
        }

        //////////////////////////////////////////////
        // Do min chi extraction
        //////////////////////////////////////////////
        double min_chi = amr_reductions.min(c_chi);
        SmallDataIO min_chi_file(m_p.data_path + "min_chi",
                                        m_dt, m_time, m_restart_time,
                                        SmallDataIO::APPEND, first_step);
        min_chi_file.remove_duplicate_time_data();
        if (first_step)
        {
            min_chi_file.write_header_line({"min_chi"});
        }
        min_chi_file.write_time_data_line({min_chi});

        //////////////////////////////////////////////
        // Do max phi extraction
        //////////////////////////////////////////////
        // double max_phi = amr_reductions.max(c_phi);
        // SmallDataIO max_phi_file(m_p.data_path + "max_phi",
        //                                 m_dt, m_time, m_restart_time,
        //                                 SmallDataIO::APPEND, first_step);
        // max_phi_file.remove_duplicate_time_data();
        // if (first_step)
        // {
        //     max_phi_file.write_header_line({"max_phi"});
        // }
        // max_phi_file.write_time_data_line({max_phi});
      }


      //////////////////////////////////////////////
      // Do Mass Charge Integration
      //////////////////////////////////////////////
      if (m_p.activate_mq_extraction == 1 &&
          at_level_timestep_multiple(
              m_p.mq_extraction_params.min_extraction_level()))
      {
          CH_TIME("EMDBHLevel::doAnalysis::MassChargeExtraction");

          fillAllGhosts();

          // Do the extraction on the min extraction level
          if (m_level == m_p.mq_extraction_params.min_extraction_level())
          {
              if (m_verbosity)
              {
                  pout() << "BinaryBSLevel::specificPostTimeStep:"
                            " Extracting MassCharge integrals."
                         << endl;
              }

              // Refresh the interpolator and do the interpolation
              m_bh_amr.m_interpolator->refresh();
              CrudeMassChargeExtraction mq_extraction
                                                (m_p.mq_extraction_params,
                                        m_dt, m_time, first_step, m_restart_time);
              mq_extraction.execute_query(m_bh_amr.m_interpolator);
          }
      }


      ////////////////////////////////////////////////////
      // Do Real Scalar Spherical Harmonic Decomposition
      ////////////////////////////////////////////////////
      if (m_p.activate_rs_extraction == 1 &&
          at_level_timestep_multiple(
              m_p.rs_extraction_params.min_extraction_level()))
      {
          CH_TIME("EMDBHLevel::doAnalysis::RealScalarExtraction");

          fillAllGhosts();

          // Do the extraction on the min extraction level
          if (m_level == m_p.rs_extraction_params.min_extraction_level())
          {
              if (m_verbosity)
              {
                  pout() << "BinaryBSLevel::specificPostTimeStep:"
                            " Extracting RealScalar integrals."
                         << endl;
              }

              // Refresh the interpolator and do the interpolation
              m_bh_amr.m_interpolator->refresh();
              RealScalarExtraction rs_extraction
                                                (m_p.rs_extraction_params,
                                        m_dt, m_time, first_step, m_restart_time);
              rs_extraction.execute_query(m_bh_amr.m_interpolator);
          }
      }

      //////////////////////////////////////////////
      // Do Electromagnetic Radiation Integration
      //////////////////////////////////////////////
      if (m_p.activate_em_extraction == 1 &&
          at_level_timestep_multiple(
              m_p.pheyl2_extraction_params.min_extraction_level()))
      {
          CH_TIME("EMSBH2DLevel::doAnalysis::EMRadExtraction");

          fillAllGhosts();
          // CouplingFunction coupling_function(m_p.coupling_function_params);
          // EinsteinMaxwellScalarFieldWithCoupling emd_field(coupling_function);

          // fill grid with pheyl2 im and real components
          auto pheyl2_compute_pack = make_compute_pack(
              Pheyl2(m_p.pheyl2_extraction_params.extraction_center,
                                          m_p.coupling_function_params, m_dx));
          BoxLoops::loop(pheyl2_compute_pack, m_state_new, m_state_diagnostics,
                         EXCLUDE_GHOST_CELLS);


          // Do the extraction on the min extraction level
          if (m_level == m_p.pheyl2_extraction_params.min_extraction_level())
          {
              if (m_verbosity)
              {
                  pout() << "EMSBH2DLevel::specificPostTimeStep:"
                            " Extracting electromagnetic waves."
                         << endl;
              }

              // Refresh the interpolator and do the interpolation
              m_bh_amr.m_interpolator->refresh();
              PheylExtraction em_extraction(m_p.pheyl2_extraction_params,
                                           m_dt, m_time,
                                           first_step, m_restart_time);
              em_extraction.execute_query(m_bh_amr.m_interpolator);
          }
      }
}
