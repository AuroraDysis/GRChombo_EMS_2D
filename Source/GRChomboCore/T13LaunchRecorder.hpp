#ifndef T13LAUNCHRECORDER_HPP_
#define T13LAUNCHRECORDER_HPP_
#include "T7OperationRecorder.hpp"
#include "UserVariables.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>

// Opt-in, write-only local launch audit. Reuse the T7 lossless frame format.
struct T13LaunchStop { double time; };
class T13LaunchRecorder
{
    T7OperationRecorder m_writer;
    std::ofstream m_floors;
    std::string m_path;
  public:
    double stop=0.;
    int steps=0;
    explicit T13LaunchRecorder(const std::string &path)
        :m_writer(path+"t13-"),m_path(path)
    {
        GRParmParse pp; pp.load("t13_launch_stop_time",stop,0.);
        if (!std::isfinite(stop) || stop<0.)
            MayDay::Error("t13_launch_stop_time must be finite and nonnegative (0 disables)");
#ifdef CH_MPI
        if (stop>0.) MayDay::Error("T13 launch capture requires local serial Chombo");
#endif
#if CH_SPACEDIM != 2
        if (stop>0.) MayDay::Error("T13 launch capture requires DIM=2");
#endif
    }
    bool enabled() const { return stop>0.; }
    bool selected(int level,int finest) const
    { return enabled() && level>=finest-1 && level<=finest; }
    bool stage_capture(int level,int finest) const
    { return enabled() && level==finest && steps<4; }
    static std::vector<Box> windows(double h,double centre)
    {
        // Fixed physical region, includes puncture cells and every W stencil.
        const int c=std::llround(centre/h),k=std::ceil(.006/h);
        return {Box(IntVect(D_DECL(c-k,0,0)),IntVect(D_DECL(c+k-1,k-1,0)))};
    }
    void stage(int index,double time,double old,double next,double fraction)
    {
        m_writer.stage=index;m_writer.stage_time=time;
        m_writer.coarse_old=old;m_writer.coarse_new=next;
        m_writer.step_fraction=fraction;
    }
    void update(double dt) {m_writer.update_dt=dt;}
    void record(int phase,int level,const GRLevelData &data,double time,
                double h,double dt,double centre,int growth=3)
    {
        int source=0;
        for (DataIterator it=data.dataIterator();it.ok();++it,++source)
        {
            const Box valid=data.disjointBoxLayout()[it()];
            Box sample=windows(h,centre)[0]&valid;
            if (sample.isEmpty()) continue;
            sample.grow(growth);sample&=data[it()].box();
            std::vector<IntVect> cells;
            for (BoxIterator bit(sample);bit.ok();++bit) cells.push_back(bit());
            m_writer.frame(phase,level,source*16,data[it()],valid,cells,0,NUM_VARS,
                           time,h,dt,0.,0.);
        }
    }
    void parts(int level,int source,const FArrayBox &data,const Box &valid,
               double time,double h,double dt,double centre)
    {
        const Box sample=windows(h,centre)[0]&valid;
        std::vector<IntVect> cells;
        for (BoxIterator it(sample);it.ok();++it) cells.push_back(it());
        for (int part=0;part<2;++part)
            m_writer.frame(20+part,level,source*16,data,valid,cells,part*NUM_VARS,
                           NUM_VARS,time,h,dt,0.,0.);
    }
    void floors(int level,const GRLevelData &data,double time,const char *where,
                double min_chi,double min_lapse,int chi_component,int lapse_component)
    {
        if (!enabled()) return;
        if (!m_floors.is_open())
        {
            m_floors.open(m_path+"t13-floors-L"+std::to_string(level)+".csv");
            m_floors<<"time,step,stage,operation,valid_cells,ghost_cells,chi_activations,lapse_activations,chi_min,lapse_min,nonfinite\n";
        }
        long long valid_count=0,ghost_count=0,chi_count=0,lapse_count=0,nan_count=0;
        double chi_min=INFINITY,lapse_min=INFINITY;
        for (DataIterator it=data.dataIterator();it.ok();++it)
        {
            const auto &fab=data[it()];const Box valid=data.disjointBoxLayout()[it()];
            for (BoxIterator bit(fab.box());bit.ok();++bit)
            {
                (valid.contains(bit())?valid_count:ghost_count)++;
                double chi=fab(bit(),chi_component),lapse=fab(bit(),lapse_component);
                chi_count+=chi<min_chi;lapse_count+=lapse<min_lapse;
                chi_min=std::min(chi_min,chi);lapse_min=std::min(lapse_min,lapse);
                for (int c=0;c<NUM_VARS;++c) nan_count+=!std::isfinite(fab(bit(),c));
            }
        }
        m_floors<<std::setprecision(17)<<time<<','<<steps<<','<<m_writer.stage<<','<<where<<','
                <<valid_count<<','<<ghost_count<<','<<chi_count<<','<<lapse_count<<','
                <<chi_min<<','<<lapse_min<<','<<nan_count<<'\n';
        m_floors.flush();
        if (!m_floors) MayDay::Error("T13 floor audit write failed");
    }
};
#endif
