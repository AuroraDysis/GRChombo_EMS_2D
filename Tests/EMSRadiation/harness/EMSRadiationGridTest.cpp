#include "GRAMR.hpp"
#include "SetupFunctions.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "BoxLoops.hpp"
#include "SetValue.hpp"
#include "EMSRadiationExtraction.hpp"

struct EMSFlatPulse
{
    double dx,time;
    std::array<double,CH_SPACEDIM> center;
    void compute(Cell<double> cell) const
    {
        const Coordinates<double> x(cell,dx,center);
        const double r=std::max(1.,std::hypot(x.x,x.y));
        const double u=time-r+1.5, sigma=.5;
        const double f=.03*std::exp(-u*u/(2*sigma*sigma))/std::sqrt(4*M_PI);
        cell.store_vars(f/r,c_phi);
        cell.store_vars(u*f/(sigma*sigma*r),c_Pi);
        for (int k:{c_chi,c_h11,c_h22,c_hww,c_lapse}) cell.store_vars(1.,k);
    }
};

class EMSFlatPulseLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void set_pulse(double t)
    {
        BoxLoops::loop(make_compute_pack(SetValue(0.),EMSFlatPulse{m_dx,t,m_p.center}),
                      m_state_new,m_state_new,INCLUDE_GHOST_CELLS,disable_simd());
        BoxLoops::loop(SetValue(0.),m_state_new,m_state_diagnostics,INCLUDE_GHOST_CELLS,disable_simd());
    }
    void initialData() override { set_pulse(0); }
    void specificEvalRHS(GRLevelData &,GRLevelData &,double) override
    { throw std::runtime_error("flat exact-solution extraction test has no evolution"); }
    void computeTaggingCriterion(FArrayBox &tags,const FArrayBox &) override
    { tags.setVal(0); }
};

int main(int argc,char **argv)
{
    mainSetup(argc,argv);
    int status=0;
    try
    {
        GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);
        SimulationParameters p(pp);
        GRAMR amr;
        DefaultLevelFactory<EMSFlatPulseLevel> factory(amr,p);
        setupAMRObject(amr,factory);
        AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,0);
        amr.set_interpolator(&interp);
        auto *level=dynamic_cast<EMSFlatPulseLevel *>(amr.getAMRLevels()[0]);
        for (int i=0;i<=240;++i)
        {
            const double t=i*.05;
            level->set_pulse(t); interp.refresh();
            EMSRadiationExtraction extraction(p.extraction_params,.05,t,i==0,0,{0,0,0,0},0);
            extraction.execute_query(&interp);
        }
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; status=1; }
    mainFinalize();
    return status;
}
