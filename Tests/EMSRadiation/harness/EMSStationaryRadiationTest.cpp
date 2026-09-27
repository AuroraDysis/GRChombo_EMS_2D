// Test-only Killing-coordinate evolution. Production's precollapsed lapse and
// ExperimentalGauge intentionally remain unchanged in Examples/EMS*.
#include "../../../Examples/EMS/EMSBH2DLevel.cpp"
#include "SetupFunctions.hpp"

struct EMSKillingLapse
{
    const EMSBH_trumpet_read &setter;
    double dx,mass;
    std::array<double,CH_SPACEDIM> center;
    void compute(Cell<double> cell) const
    {
        const Coordinates<double> x(cell,dx,center);
        cell.store_vars(setter.compute_ems_adm_vars(x.x,x.y,0.,mass,0.,0.).lapse,c_lapse);
    }
};

class EMSStationaryRadiationLevel : public EMSBH2DLevel
{
  public:
    using EMSBH2DLevel::EMSBH2DLevel;
    void initialData() override
    {
        if (m_p.emsbh_params.binary || m_p.emsbh_params.boosted)
            throw std::runtime_error("stationary control requires one unboosted trumpet");
        EMSBH_trumpet_read setter(m_p.emsbh_params,m_p.coupling_function_params,m_p.m_G_Newton,m_dx,0);
        setter.compute_1d_solution();
        BoxLoops::loop(make_compute_pack(SetValue(0.),setter,
                      EMSKillingLapse{setter,m_dx,m_p.emsbh_params.bh_mass,m_p.center}),
                      m_state_new,m_state_new,INCLUDE_GHOST_CELLS,disable_simd());
        fillAllGhosts();
        BoxLoops::loop(GammaCartoonCalculator(m_dx),m_state_new,m_state_new,
                      EXCLUDE_GHOST_CELLS,disable_simd());
    }
    void specificEvalRHS(GRLevelData &sol,GRLevelData &rhs,double) override
    {
        BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(),
            PositiveChiAndAlpha(m_p.min_chi,m_p.min_lapse)),sol,sol,INCLUDE_GHOST_CELLS);
        CouplingFunction coupling(m_p.coupling_function_params);
        CCZ4Cartoon<ExperimentalGauge,FourthOrderDerivatives,CouplingFunction>
            equations(m_p.ccz4_params,m_dx,m_p.sigma,coupling,m_p.m_G_Newton,m_p.formulation);
        BoxLoops::loop(equations,sol,rhs,EXCLUDE_GHOST_CELLS);
        // Freeze gauge only, after dissipation as well. Every physical field
        // evolves under the same CCZ4/EMS RHS; this is not a frozen snapshot.
        BoxLoops::loop(SetValue(0.,Interval(c_lapse,c_B2)),rhs,rhs,EXCLUDE_GHOST_CELLS);
    }
};

int main(int argc,char **argv)
{
    mainSetup(argc,argv);
    int status=0;
    try
    {
        GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);
        SimulationParameters p(pp);
        BHAMR amr;
        DefaultLevelFactory<EMSStationaryRadiationLevel> factory(amr,p);
        setupAMRObject(amr,factory);
        AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,0);
        amr.set_interpolator(&interp);
        auto *level=dynamic_cast<EMSStationaryRadiationLevel *>(amr.getAMRLevels()[0]);
        level->ems_prepare_radiation(); level->ems_extract_radiation();
        amr.run(p.stop_time,p.max_steps); amr.conclude();
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; status=1; }
    mainFinalize();
    return status;
}
