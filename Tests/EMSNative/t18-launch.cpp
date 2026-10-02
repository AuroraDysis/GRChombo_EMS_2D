// Read-only launch capture on translated, exact T17 refined boxes.
#define main t18_window_main
#include "t18-native.cpp"
#undef main
#include <sstream>

class T18LaunchLevel : public T18Level
{
    T7OperationRecorder writer;
    double end;
  public:
    T18LaunchLevel(GRAMR &g,const SimulationParameters &p,int v)
        :T18Level(g,p,v),writer(p.data_path+"t18-launch-")
    {GRParmParse pp;pp.load("t18_capture_stop",end);}
    void snapshot()
    {
        if(m_level!=m_p.max_level)return;
        for(double sign:{-1.,1.})
        {
            const double centre=m_p.center[0]+sign*m_p.emsbh_params.separation/2;
            const auto windows=T13LaunchRecorder::windows(m_dx,centre);
            writer.record(50,m_level,m_state_new,windows,3,0,NUM_VARS,m_time,m_dx,m_dt,0.,0.);
        }
    }
    void postTimeStep() override
    {
        T18Level::postTimeStep();
        snapshot();
        if(m_level==m_p.max_level && m_time+m_dt>end+1e-14)throw T18Stop{};
    }
};
class T18LaunchFactory : public DefaultLevelFactory<EMSBH2DLevel>
{
  public:
    using DefaultLevelFactory<EMSBH2DLevel>::DefaultLevelFactory;
    AMRLevel *new_amrlevel() const override
    {auto *l=new T18LaunchLevel(m_gr_amr,m_p,m_p.verbosity);l->initialDtMultiplier(m_p.dt_multiplier);return l;}
};
int main(int argc,char **argv)
{
    mainSetup(argc,argv);GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);SimulationParameters p(pp);
    BHAMR amr;T18LaunchFactory factory(amr,p);
    ProblemDomain domain(IntVect::Zero,p.ivN);
    amr.define(p.max_level,p.ref_ratios,domain,&factory);
    amr.gridBufferSize(p.grid_buffer_size);amr.blockFactor(p.block_factor);amr.maxGridSize(p.max_grid_size);
    amr.regridIntervals(p.regrid_interval);amr.fillRatio(p.fill_ratio);amr.verbosity(p.verbosity);
    amr.checkpointInterval(-1);amr.plotInterval(-1);
    Vector<Vector<Box>> grids(p.max_level+1);amr.makeBaseLevelMesh(grids[0]);
    std::string path;pp.load("t18_capture_boxes",path);std::ifstream in(path);
    int level,x0,y0,x1,y1;
    while(in>>level>>x0>>y0>>x1>>y1)
        grids[level].push_back(Box(IntVect(x0,y0),IntVect(x1,y1)));
    if(!in.eof())MayDay::Error("T18 capture box read failed");
    for(int l=1;l<=p.max_level;++l)if(grids[l].size()==0)MayDay::Error("T18 missing refined level");
    amr.setupForFixedHierarchyRun(grids);
    AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,p.verbosity);amr.set_interpolator(&interp);
    const auto levels=amr.getAMRLevels();
    for(int l=0;l<levels.size();++l)dynamic_cast<T18LaunchLevel *>(levels[l])->snapshot();
    bool stopped=false;
    try{amr.run(p.stop_time,p.max_steps);}catch(const T18Stop &){stopped=true;}
    if(!stopped)MayDay::Error("T18 launch did not reach registered stop");
    std::ofstream out(p.data_path+"t18-capture-stop.txt");
    out<<std::setprecision(17)<<dynamic_cast<T18LaunchLevel *>(levels[levels.size()-1])->m_time<<'\n';
    if(!out)MayDay::Error("T18 capture stop write failed");
    mainFinalize();return 0;
}
