/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#include <fstream>
#include <mutex>
#include <memory>
#include <cstdint>
#include <typeinfo>

#include "EMSBH2DLevel.hpp"
#include "EMSRadiationExtraction.hpp"
#include "BoxIterator.H"

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
#include "MovingPunctureGauge.hpp"
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

void EMSBH2DLevel::t6_write_interface_strips() const
{
#ifdef CH_MPI
    MayDay::Error("t6_interface_diagnostics is a local serial diagnostic");
#endif
    const auto extent = [](const GRAMRLevel &level) {
        double a = 0.;
        const auto &layout = level.getLevelData().disjointBoxLayout();
        for (LayoutIterator it = layout.layoutIterator(); it.ok(); ++it)
            a = std::max(a, (layout[it()].bigEnd(1)+1)*level.get_dx());
        return a;
    };
    const auto levels = m_gr_amr.get_gramrlevels();
    const GRAMRLevel *finer = m_level+1 < int(levels.size()) ? levels[m_level+1] : nullptr;
    const double a = extent(*this), b = finer ? extent(*finer) : -1.;
    const double width = 2*m_dx;
    std::ofstream out(m_p.data_path + "strips/t6-level" + std::to_string(m_level) +
                      "-t" + std::to_string(m_time) + ".bin", std::ios::binary);
    if (!out) MayDay::Error("cannot open T6 interface strip");
    out.write("T6STRIP1", 8);
    const double meta[] = {m_time, m_dx, a, b};
    const std::uint32_t header[] = {std::uint32_t(m_level), NUM_VARS, NUM_VARS+8};
    out.write(reinterpret_cast<const char *>(meta), sizeof(meta));
    out.write(reinterpret_cast<const char *>(header), sizeof(header));
    const int constraints[] = {c_Ham, c_Mom1, c_Mom2, c_GaussE, c_GaussB};
    const auto near = [width](double x, double y, double face) {
        if (face < 0.) return false;
        x = std::abs(x);
        return (x <= face && y <= face)
            ? std::min(face-x, face-y) <= width
            : std::hypot(std::max(x-face, 0.), std::max(y-face, 0.)) <= width;
    };
    for (DataIterator it = m_state_new.dataIterator(); it.ok(); ++it)
        for (BoxIterator cell(m_state_new.disjointBoxLayout()[it()]); cell.ok(); ++cell)
        {
            const auto iv = cell();
            const double x = (iv[0]+.5)*m_dx-m_p.center[0], y = (iv[1]+.5)*m_dx;
            if (!near(x, y, a) && !near(x, y, b)) continue;
            double covered = 0.;
            if (finer)
            {
                const auto &layout = finer->getLevelData().disjointBoxLayout();
                for (LayoutIterator fine = layout.layoutIterator(); fine.ok(); ++fine)
                {
                    Box box = layout[fine()];
                    box.coarsen(2);
                    if (box.contains(iv)) { covered = 1.; break; }
                }
            }
            std::array<double, NUM_VARS+8> row;
            row[0] = x; row[1] = y; row[2] = covered;
            for (int c=0; c<NUM_VARS; ++c) row[3+c] = m_state_new[it()](iv, c);
            for (int c=0; c<5; ++c)
                row[3+NUM_VARS+c] = m_state_diagnostics[it()](iv, constraints[c]);
            out.write(reinterpret_cast<const char *>(row.data()), sizeof(row));
        }
    if (!out) MayDay::Error("T6 interface strip write failed");
}

void EMSBH2DLevel::specificAdvance()
{
    if (m_t7.enabled) t7_record(6,m_state_new,0);
    const bool capture=m_t13.stage_capture(m_level,m_p.max_level);
    if (capture) t13_record(6,m_state_new,m_time+m_dt);
    m_t13.floors(m_level,m_state_new,m_time+m_dt,"advance",m_p.min_chi,m_p.min_lapse,c_chi,c_lapse);
    // Enforce the trace free A_ij condition and positive chi and alpha
    if (capture)
    {
        BoxLoops::loop(TraceARemovalCartoon(),m_state_new,m_state_new,INCLUDE_GHOST_CELLS);
        t13_record(16,m_state_new,m_time+m_dt);
        BoxLoops::loop(PositiveChiAndAlpha(m_p.min_chi,m_p.min_lapse),m_state_new,m_state_new,INCLUDE_GHOST_CELLS);
        t13_record(17,m_state_new,m_time+m_dt);
    }
    else
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    if (m_t7.enabled) t7_record(7,m_state_new,0);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void EMSBH2DLevel::ems_t13_initial()
{
    if (!m_t13.enabled()) return;
    if (m_time!=0.) MayDay::Error("T13 launch audit requires a t=0 start");
    m_t13.floors(m_level,m_state_new,m_time,"initial",m_p.min_chi,m_p.min_lapse,c_chi,c_lapse);
    t13_snapshot();
}

void EMSBH2DLevel::ems_t7_initial()
{
    if (!m_t7.enabled || m_level<4 || m_level>6) return;
    // Initial diagnostics must also preserve the stored, unfilled ghost bits.
    std::vector<std::unique_ptr<FArrayBox>> saved;
    for (DataIterator it=m_state_new.dataIterator();it.ok();++it)
    {
        saved.emplace_back(new FArrayBox(m_state_new[it()].box(),NUM_VARS));
        saved.back()->copy(m_state_new[it()]);
    }
    prePlotLevel();int i=0;
    for (DataIterator it=m_state_new.dataIterator();it.ok();++it,++i)
        m_state_new[it()].copy(*saved[i]);
}

void EMSBH2DLevel::t7_snapshot()
{
    if (m_level<4 || m_level>6 || m_time==m_t7_snapshot_time ||
        std::abs(12*m_time-std::round(12*m_time))>1e-9) return;
    m_t7_snapshot_time=m_time;auto f=t7_faces();int source=0;
    for (DataIterator it=m_state_new.dataIterator();it.ok();++it,++source)
    {
        Box b=m_state_new.disjointBoxLayout()[it()];std::vector<IntVect> cells;
        FArrayBox data(b,NUM_VARS+5);
        for (BoxIterator bit(b);bit.ok();++bit)
        {
            double x=(bit()[0]+.5)*m_dx-m_p.center[0],y=(bit()[1]+.5)*m_dx,r=std::hypot(x,y);
            bool selected=(r>.5 && r<6.5 && (y<.05 || std::abs(x)<.05 || std::abs(y-std::abs(x))<.035));
            for (double a:{f.first,f.second}) if (a>0 && a<=4.+1e-10)
                selected|=std::abs(x)<=a+.1875 && y<=a+.1875 &&
                    (std::abs(std::abs(x)-a)<.1875 || std::abs(y-a)<.1875);
            if (!selected) continue;
            cells.push_back(bit());
            for (int c=0;c<NUM_VARS;++c) data(bit(),c)=m_state_new[it()](bit(),c);
            int k=NUM_VARS;
            for (int c:{c_Ham,c_Mom1,c_Mom2,c_GaussE,c_GaussB}) data(bit(),k++)=m_state_diagnostics[it()](bit(),c);
        }
        m_t7.frame(50,m_level,source,data,b,cells,0,NUM_VARS+5,m_time,m_dx,m_dt,f.first,f.second,true);
    }
}


// Initial data for field and metric variables
void EMSBH2DLevel::initialData()
{
    CH_TIME("EMSBH2DLevel::initialData");
    if ((m_p.amr_transfer == "point" || m_p.t2_guard_initial_data_after_t0) &&
        (m_time > 0. || m_gr_amr.get_current_time() > 0.))
        MayDay::Error("initialData called after t=0 (incomplete restart hierarchy)");
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

    if (m_ems_gauge.differential())
        BoxLoops::loop(SetValue(0.0, Interval(c_B1, c_B2)), m_state_new,
                       m_state_new, INCLUDE_GHOST_CELLS);
    else
    {
        auto my_gauge_conditions = ExperimentalGauge(m_p.ccz4_params);
        BoxLoops::loop(my_gauge_conditions,
                      m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    }
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
    if (m_t7.enabled) t7_snapshot();
}

#define EMS_RHS_METHOD specificEvalRHS
#define EMS_RHS_GAUGE ExperimentalGauge
#define EMS_RHS_DISPATCH_MOVING
#include "EMSBH2DRHS.impl.hpp"

void EMSBH2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    if (m_t7.enabled) t7_record(4,a_soln,0);
    const bool capture=m_t13.stage_capture(m_level,m_p.max_level);
    if (capture) t13_record(4,a_soln,m_time,0);
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
    if (capture) t13_record(5,a_soln,m_time,0);
    if (m_t7.enabled) t7_record(5,a_soln,0);
}



void EMSBH2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state)
{
    auto midpoint_params = m_p.mass_extraction_params;
    if (m_p.ems_binary_refinement)
        midpoint_params.num_extraction_radii = 0;
    BoxLoops::loop(EMSExtractionTaggingCriterion(
                     m_dx, m_level, midpoint_params,
                     m_p.regrid_threshold_A, m_p.regrid_threshold_phi,
                     m_p.regrid_threshold_chi), current_state,
                                                tagging_criterion);

    // Opt-in: reuse the author's criterion at both punctures. The outer
    // per-hole regions overlap to form the shared exterior hierarchy.
    if (m_p.ems_binary_refinement)
    {
        auto centres = m_bh_amr.m_puncture_tracker.get_puncture_coords();
        if (centres.empty())
        {
            auto left = m_p.emsbh_params.star_centre;
            auto right = left;
            left[0] -= .5 * m_p.emsbh_params.separation;
            right[0] += .5 * m_p.emsbh_params.separation;
            centres = {left, right};
        }
        for (const auto &centre : centres)
        {
            auto params = m_p.mass_extraction_params;
            params.extraction_center = centre;
            FArrayBox local(tagging_criterion.box(), 1);
            BoxLoops::loop(EMSExtractionTaggingCriterion(m_dx, m_level, params,
                m_p.regrid_threshold_A, m_p.regrid_threshold_phi,
                m_p.regrid_threshold_chi), current_state, local);
            for (BoxIterator bit(local.box()); bit.ok(); ++bit)
                tagging_criterion(bit(), 0) = std::max(
                    tagging_criterion(bit(), 0), local(bit(), 0));
        }
    }

    if (m_radiation.active)
        for (BoxIterator bit(tagging_criterion.box()); bit.ok(); ++bit)
        {
            const double x = (bit()[0] + .5) * m_dx - m_p.extraction_params.center[0];
            const double y = (bit()[1] + .5) * m_dx - m_p.extraction_params.center[1];
            const double r = std::hypot(x, y);
            if (m_level < m_radiation.wave_level && r < 1.2 * m_radiation.wave_radius)
                tagging_criterion(bit(), 0) = 100.;
            for (int i = 0; i < m_p.extraction_params.num_extraction_radii; ++i)
                if (m_level < m_p.extraction_params.extraction_levels[i] &&
                    r < 1.2 * m_p.extraction_params.extraction_radii[i])
                    tagging_criterion(bit(), 0) = 100.;
        }
}

void EMSBH2DLevel::specificPostTimeStep()
{
    if (m_p.ems_track_punctures && m_level == m_p.ems_puncture_tracking_level)
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
            m_dt, at_level_timestep_multiple(0));
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

    if (m_t6_interface_diagnostics && m_level >= 4 && at_level_timestep_multiple(0))
        t6_write_interface_strips();
    if (m_t7.enabled && at_level_timestep_multiple(0)) t7_snapshot();

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
    if (m_radiation.active &&
        m_level == m_p.extraction_params.min_extraction_level() &&
        at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
        ems_extract_radiation();
}

void EMSBH2DLevel::ems_prepare_radiation()
{
    fillAllGhosts();
    BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void EMSBH2DLevel::ems_extract_radiation()
{
    m_gr_amr.m_interpolator->refresh();
    EMSRadiationExtraction extraction(m_p.extraction_params, m_dt, m_time,
        m_time == 0., m_restart_time, m_p.coupling_function_params, m_radiation.phi_inf);
    extraction.execute_query(m_gr_amr.m_interpolator);
}

#ifdef CH_USE_HDF5
void EMSBH2DLevel::writeCheckpointHeader(HDF5Handle &handle) const
{
    GRAMRLevel::writeCheckpointHeader(handle);
    HDF5HeaderData header;
    header.m_string["ems_gauge"] = m_ems_gauge.name();
    header.m_string["ems_driver_semantics"] = m_ems_gauge.driver_semantics();
    header.writeToFile(handle);
}

void EMSBH2DLevel::readCheckpointHeader(HDF5Handle &handle)
{
    HDF5HeaderData header;
    header.readFromFile(handle);
    const auto found = header.m_string.find("ems_gauge");
    m_ems_gauge.check_checkpoint(found == header.m_string.end()
                                ? "experimental" : found->second);
    const auto driver = header.m_string.find("ems_driver_semantics");
    if (driver != header.m_string.end() && driver->second != m_ems_gauge.driver_semantics())
        MayDay::Error("EMS checkpoint driver semantics do not match selected gauge");
    GRAMRLevel::readCheckpointHeader(handle);
}
#endif
