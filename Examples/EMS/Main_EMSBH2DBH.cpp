/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CH_Timer.H"
#include "parstream.H" //Gives us pout()
#include <chrono>
#include <iostream>

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "EMSBH2DLevel.hpp"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);
    const EMSRadiationParameters radiation(pp, sim_params);
    const double rh_threshold = ems_rh_expansion_threshold(pp);

    if (sim_params.just_check_params)
        return 0;

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    BHAMR bh_amr;
    bh_amr.m_rh_union.m_thresh_super_low = rh_threshold;
    DefaultLevelFactory<EMSBH2DLevel> emdbh_level_fact(bh_amr, sim_params);
    setupAMRObject(bh_amr, emdbh_level_fact);

    // call this after amr object setup so grids known
    // and need it to stay in scope throughout run
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bh_amr.set_interpolator(&interpolator);
    if (radiation.active)
    {
        const auto levels = bh_amr.getAMRLevels();
        for (int i = 0; i < levels.size(); ++i)
            dynamic_cast<EMSBH2DLevel *>(levels[i])->ems_prepare_radiation();
        dynamic_cast<EMSBH2DLevel *>(levels[
            sim_params.extraction_params.min_extraction_level()])->ems_extract_radiation();
    }

    #ifdef USE_AHFINDER
        // one horizon
        if (sim_params.AH_activate && sim_params.AH_num_horizons==1)
        {
            AHSphericalGeometry sph1(sim_params.horizon_centre_1);
            bh_amr.m_ah_finder.add_ah(sph1, sim_params.AH_initial_guess,
                                      sim_params.AH_params);
        }
        // two horizons
        else if (sim_params.AH_activate && sim_params.AH_num_horizons==2)
        {
            AHSphericalGeometry sph1(sim_params.horizon_centre_1);
            AHSphericalGeometry sph2(sim_params.horizon_centre_2);
            bh_amr.m_ah_finder.add_ah(sph1, sim_params.AH_initial_guess,
                                      sim_params.AH_params);
            bh_amr.m_ah_finder.add_ah(sph2, sim_params.AH_initial_guess,
                                      sim_params.AH_params);
            if (sim_params.AH_expect_merger)
            {
                bh_amr.m_ah_finder.add_ah_merger(0, 1, sim_params.AH_params);
            }
        }
    #endif

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

    bh_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    bh_amr.conclude();

    CH_TIMER_REPORT(); // Report results when running with Chombo timers.

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
