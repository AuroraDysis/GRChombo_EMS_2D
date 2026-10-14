// Real t=0 hierarchy from the production trumpet reader; no evolution or RH finder.
#include "AMRInterpolator.hpp"
#include "BoxLoops.hpp"
#include "DefaultLevelFactory.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "EMSMatterWeylOmScalar.hpp"
#include "EMSRadiationExtraction.hpp"
#include "EMSRadiationParameters.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "SetValue.hpp"
#include "SetupFunctions.hpp"
#include "../../EMSNP/ValidationExtraction.hpp"

class NPStaticLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void initialData() override
    {
        EMSBH_trumpet_read setter(m_p.emsbh_params,m_p.coupling_function_params,m_p.m_G_Newton,m_dx,0);
        setter.compute_1d_solution();
        BoxLoops::loop(make_compute_pack(SetValue(0.),setter),m_state_new,m_state_new,
            INCLUDE_GHOST_CELLS,disable_simd());
    }
    void specificEvalRHS(GRLevelData &,GRLevelData &,double) override
    { throw std::runtime_error("NP static control must never advance"); }
    void computeTaggingCriterion(FArrayBox &tags,const FArrayBox &) override { tags.setVal(0.); }
    void weyl()
    {
        fillAllGhosts();
        BoxLoops::loop(EMSMatterWeylOmScalar(m_p.extraction_params.center,m_dx,
            m_p.coupling_function_params,m_p.m_G_Newton),m_state_new,m_state_diagnostics,EXCLUDE_GHOST_CELLS);
    }
};
int run(int argc,char **argv)
{
    GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);
    SimulationParameters p(pp);
    EMSRadiationParameters radiation(pp,p);
    if (p.restart_from_checkpoint || !radiation.np_active)
        throw std::runtime_error("static control requires initial data and ems_np_activate");
    GRAMR amr;
    DefaultLevelFactory<NPStaticLevel> factory(amr,p);
    setupAMRObject(amr,factory);
    for (auto *level:amr.get_gramrlevels()) dynamic_cast<NPStaticLevel &>(*level).weyl();
    AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,0);
    amr.set_interpolator(&interp); interp.refresh();
    EMSNPValidationExtraction np(p.extraction_params,p.dx[0]*p.dt_multiplier,0,true,0.,
        p.coupling_function_params,p.m_G_Newton);
    np.execute_query(&interp);
    return 0;
}
int main(int argc,char **argv)
{
    mainSetup(argc,argv);
    int status=0;
    try { status=run(argc,argv); }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; status=2; }
    mainFinalize(); return status;
}
