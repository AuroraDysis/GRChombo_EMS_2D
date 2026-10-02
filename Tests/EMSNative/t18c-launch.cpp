// Fixed native Chombo hierarchy; T16b's unchanged audit and launch observers.
#include "t18c-native.hpp"
#include "SetupFunctions.hpp"
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
int main(int argc,char **argv)
{
    mainSetup(argc,argv);GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);SimulationParameters p(pp);
    BHAMR amr;T18CFactory<EMSBH2DLevel> factory(amr,p);
    ProblemDomain domain(IntVect::Zero,p.ivN);
    for(int d=0;d<CH_SPACEDIM;++d)domain.setPeriodic(d,p.boundary_params.is_periodic[d]);
    amr.define(p.max_level,p.ref_ratios,domain,&factory);
    amr.gridBufferSize(p.grid_buffer_size);amr.blockFactor(p.block_factor);amr.maxGridSize(p.max_grid_size);
    amr.regridIntervals(p.regrid_interval);amr.fillRatio(p.fill_ratio);amr.verbosity(p.verbosity);
    amr.checkpointInterval(-1);amr.plotInterval(-1);
    Vector<Vector<Box>> grids(p.max_level+1);std::string path;pp.load("t18c_boxes",path);
    if(path.empty())amr.makeBaseLevelMesh(grids[0]);
    else
    {
        std::ifstream in(path);int l,x0,y0,x1,y1;
        while(in>>l>>x0>>y0>>x1>>y1)grids[l].push_back(Box(IntVect(x0,y0),IntVect(x1,y1)));
        if(!in.eof())MayDay::Error("T18c box-file read failed");
    }
    if(grids[0].size()==0)amr.makeBaseLevelMesh(grids[0]);
    for(int l=0;l<=p.max_level;++l)if(grids[l].size()==0)MayDay::Error("T18c missing level");
    amr.setupForFixedHierarchyRun(grids);
    AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,p.verbosity);amr.set_interpolator(&interp);
    t18c_initial_audit(amr,p);
    bool audit=false;pp.load("t18c_audit_only",audit,false);
    if(audit){pout()<<"T18c audit-only finished."<<std::endl;mainFinalize();return 0;}
    pout()<<"T18c evolution starts after completed t=0 audit"<<std::endl;
    bool stopped=false;try{amr.run(p.stop_time,p.max_steps);}catch(const T13LaunchStop &){stopped=true;}
    if(!stopped)MayDay::Error("T18c launch did not reach registered stop");
    pout()<<"T18c launch finished."<<std::endl;mainFinalize();return 0;
}
