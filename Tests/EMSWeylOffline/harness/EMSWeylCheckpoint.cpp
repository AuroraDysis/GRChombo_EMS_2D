// Weyl4 on a frozen production checkpoint, twice: the production vacuum form (WeylOmScalar) and the matter-complete form
// (EMSMatterWeylOmScalar), each followed by the production radiation extraction on the parameter file's spheres. No initial data,
// evolution or regrid. Writes vacuum/ems_radiation.csv and matter/ems_radiation.csv in the working directory.
// usage: EMSWeylCheckpoint<config>.ex <params with restart_file and the spheres>
#include "AMRInterpolator.hpp"
#include "BoxLoops.hpp"
#include "DefaultLevelFactory.hpp"
#include "EMSMatterWeylOmScalar.hpp"
#include "EMSRadiationExtraction.hpp"
#include "EMSRadiationParameters.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "Lagrange.hpp"
#include "SetupFunctions.hpp"
#include "WeylOmScalar.hpp"

class FrozenLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void initialData() override { throw std::runtime_error("restart_file required"); }
    void specificEvalRHS(GRLevelData &, GRLevelData &, double) override
    { throw std::runtime_error("offline snapshot must never advance"); }
    void computeTaggingCriterion(FArrayBox &, const FArrayBox &) override
    { throw std::runtime_error("offline snapshot must never regrid"); }

    // the production order: ghosts, then Weyl4 into the diagnostics on valid cells
    void weyl(bool matter)
    {
        fillAllGhosts();
        if (matter)
            BoxLoops::loop(EMSMatterWeylOmScalar(m_p.extraction_params.center, m_dx,
                                                 m_p.coupling_function_params, m_p.m_G_Newton),
                           m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        else
            BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
};

int run(int argc, char **argv)
{
    GRParmParse pp(argc - 2, argv + 2, nullptr, argv[1]);
    SimulationParameters p(pp);
    const EMSRadiationParameters radiation(pp, p);
    if (!p.restart_from_checkpoint || !radiation.active)
        throw std::runtime_error("restart_file and ems_radiation_activate = true required");
    GRAMR amr;
    DefaultLevelFactory<FrozenLevel> factory(amr, p);
    setupAMRObject(amr, factory);
    const double time = amr.getAMRLevels()[0]->time();
    const auto *coarse = amr.get_gramrlevels()[0];
    if (coarse->get_dx() != p.dx[0] || coarse->problemDomain().domainBox().bigEnd() != p.ivN)
        throw std::runtime_error("parameter grid does not match checkpoint");
    for (const auto *level : amr.get_gramrlevels())
        if (level->time() != time) throw std::runtime_error("checkpoint levels are not synchronized");
    AMRInterpolator<Lagrange<4>> interp(amr, p.origin, p.dx, p.boundary_params, 0);
    amr.set_interpolator(&interp);
    for (const bool matter : {false, true})
    {
        for (auto *level : amr.get_gramrlevels())
            dynamic_cast<FrozenLevel &>(*level).weyl(matter);
        interp.refresh();
        auto params = p.extraction_params;
        params.data_path = matter ? "matter/" : "vacuum/";
        EMSRadiationExtraction extraction(params, p.dx[0] * p.dt_multiplier, time, true, 0.,
                                          p.coupling_function_params, radiation.phi_inf);
        extraction.execute_query(&interp);
        pout() << "EMS_WEYL_OFFLINE t=" << time << " form=" << (matter ? "matter" : "vacuum") << std::endl;
    }
    return 0;
}

int main(int argc, char **argv)
{
    mainSetup(argc, argv);
    int status;
    try { status = run(argc, argv); }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; status = 2; }
    mainFinalize();
    return status;
}
