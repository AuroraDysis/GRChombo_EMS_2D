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
#include "../../EMSNP/CheckpointJet.hpp"
#include "../../EMSNP/ValidationExtraction.hpp"
#include "../../EMSNP/SixthWeyl.hpp"
#include "../../EMSNP/OrderStudy.hpp"

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

    void np_probe(double x, double y, bool matter, const std::string &path)
    {
        BoxLoops::loop(EMSNPCheckpointJet(m_p.extraction_params.center, m_dx,
            m_p.coupling_function_params, m_p.m_G_Newton, x, y, matter, path),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd());
    }

    void sixth_weyl(bool matter)
    {
        BoxLoops::loop(EMSNPSixthWeyl(m_p.extraction_params.center,m_dx,
            m_p.coupling_function_params,m_p.m_G_Newton,matter),m_state_new,
            m_state_diagnostics,EXCLUDE_GHOST_CELLS,disable_simd());
    }
};

int run(int argc, char **argv)
{
    GRParmParse pp(argc - 2, argv + 2, nullptr, argv[1]);
    SimulationParameters p(pp);
    const EMSRadiationParameters radiation(pp, p);
    int probe_level=-1;
    pp.query("np_probe_level", probe_level);
    double probe_x=0, probe_y=0;
    if (probe_level>=0)
    {
        pp.load("np_probe_x", probe_x);
        pp.load("np_probe_y", probe_y);
        if (probe_level>p.max_level || !std::isfinite(probe_x) || !std::isfinite(probe_y) || probe_y<=0)
            throw std::runtime_error("invalid NP probe point/level");
    }
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
    bool order_study=false;
    pp.query("np_order_study",order_study);
    if (order_study)
    {
        int repeats=3;
        pp.query("np_order_repeats",repeats);
        if (!radiation.np_active || repeats<1 || repeats>20)
            throw std::runtime_error("NP order study requires NP active and 1..20 repeats");
        std::ofstream cost("cost.csv");
        if (!cost) throw std::runtime_error("cannot write order-study cost.csv");
        cost << "order,matter,checked,route,sample,seconds\n" << std::setprecision(17);
        bool bounded=true;
        pp.query("np_order_bounded",bounded);
        if (!bounded)
        {
            for (auto *level:amr.get_gramrlevels()) dynamic_cast<FrozenLevel &>(*level).weyl(false);
            ems_np_order<6,true,false>(amr,p,radiation,false,cost,repeats);
            throw std::runtime_error("negative unbounded-stencil control unexpectedly returned");
        }
        for (bool matter:{false,true})
        {
            for (auto *level:amr.get_gramrlevels())
                dynamic_cast<FrozenLevel &>(*level).weyl(matter);
            ems_np_order<2>(amr,p,radiation,matter,cost,repeats);
            ems_np_order<4>(amr,p,radiation,matter,cost,repeats);
            ems_np_order<6>(amr,p,radiation,matter,cost,repeats);
            ems_np_order<4,false>(amr,p,radiation,matter,cost,repeats);
        }
        cost.close();
        if (!cost) throw std::runtime_error("failed writing order-study costs");
        return 0;
    }
    AMRInterpolator<Lagrange<4>> interp(amr, p.origin, p.dx, p.boundary_params, 0);
    amr.set_interpolator(&interp);
    for (const bool matter : {false, true})
    {
        for (auto *level : amr.get_gramrlevels())
            dynamic_cast<FrozenLevel &>(*level).weyl(matter);
        if (probe_level>=0)
        {
            if (probe_level>=int(amr.get_gramrlevels().size()))
                throw std::runtime_error("NP probe level absent from checkpoint");
            dynamic_cast<FrozenLevel &>(*amr.get_gramrlevels()[probe_level]).np_probe(
                probe_x, probe_y, matter, matter ? "matter/np_jet.txt" : "vacuum/np_jet.txt");
        }
        interp.refresh();
        auto params = p.extraction_params;
        params.data_path = matter ? "matter/" : "vacuum/";
        EMSRadiationExtraction extraction(params, p.dx[0] * p.dt_multiplier, time, true, 0.,
                                          p.coupling_function_params, radiation.phi_inf);
        extraction.execute_query(&interp);
        if (radiation.np_active)
        {
            EMSNPValidationExtraction np(params,p.dx[0]*p.dt_multiplier,time,true,0.,
                p.coupling_function_params,p.m_G_Newton,matter);
            np.execute_query(&interp);
        }
        pout() << "EMS_WEYL_OFFLINE t=" << time << " form=" << (matter ? "matter" : "vacuum") << std::endl;
    }
    if (radiation.np_active)
    {
        for (bool matter:{false,true})
        {
            for (auto *level:amr.get_gramrlevels())
                dynamic_cast<FrozenLevel &>(*level).sixth_weyl(matter);
            interp.refresh();
            auto params=p.extraction_params;
            params.data_path=matter ? "sixth-matter/" : "sixth-vacuum/";
            EMSRadiationExtraction radiation6(params,p.dx[0]*p.dt_multiplier,time,true,0.,
                p.coupling_function_params,radiation.phi_inf);
            radiation6.execute_query(&interp);
        }
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
