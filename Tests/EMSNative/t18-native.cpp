// Test-only observation around unmodified EMS evolution and point transfers.
#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "EMSBH2DLevel.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "BoxIterator.H"
#include <chrono>
#include <fstream>
#include <iomanip>

struct T18Stop {};
class T18Level : public EMSBH2DLevel
{
    bool monitor=false;
    double stop=0.;
    long long samples=0,bad=0,chi_floors=0,lapse_floors=0;
    double alpha_max=0.,chi_min=INFINITY,alpha_min=INFINITY;
    double light=0.,lapse=0.,gauge=0.,beta_max=0.,puncture_max=0.;
    GRLevelData initial;
    double gamma_delta_W=0.,lapse_delta_W=0.,shift_delta_W=0.;
    std::chrono::steady_clock::time_point started;
  public:
    T18Level(GRAMR &g,const SimulationParameters &p,int v):EMSBH2DLevel(g,p,v)
    {GRParmParse pp;pp.load("t18_monitor",monitor,false);pp.load("t18_stop_time",stop,0.);}
    void census(const GRLevelData &data)
    {
        const bool first=!initial.isDefined();
        if(first)initial.define(data.disjointBoxLayout(),5,IntVect::Zero);
        gamma_delta_W=lapse_delta_W=shift_delta_W=0.;
        ++samples;
        for(DataIterator it=data.dataIterator();it.ok();++it)
        {
            const auto &u=data[it()];const Box b=data.disjointBoxLayout()[it()];
            for(BoxIterator bit(b);bit.ok();++bit)
            {
                const auto iv=bit();const double c=u(iv,c_chi),a=u(iv,c_lapse);
                chi_min=std::min(chi_min,c);alpha_min=std::min(alpha_min,a);alpha_max=std::max(alpha_max,a);
                chi_floors+=c<m_p.min_chi;lapse_floors+=a<m_p.min_lapse;
                for(int n=0;n<NUM_VARS;++n)bad+=!std::isfinite(u(iv,n));
                const double h11=u(iv,c_h11),h12=u(iv,c_h12),h22=u(iv,c_h22),hww=u(iv,c_hww);
                const double det=h11*h22-h12*h12;
                if(!(det>0.)||!(h11>0.)||!(hww>0.)){++bad;continue;}
                const double inv=std::max(1./hww,(h11+h22+std::hypot(h11-h22,2*h12))/(2*det));
                const double beta=std::hypot(u(iv,c_shift1),u(iv,c_shift2));
                beta_max=std::max(beta_max,beta);
                light=std::max(light,beta+a*std::sqrt(c*inv));
                lapse=std::max(lapse,beta+std::sqrt(1.8*a*c*inv));
                gauge=std::max(gauge,beta+std::sqrt(inv));
                const double x=(iv[0]+.5)*m_dx-m_p.emsbh_params.star_centre[0];
                const double y=(iv[1]+.5)*m_dx;
                if(x*x+y*y<.08*.08)
                    for(int n:{c_Gamma1,c_Gamma2,c_K,c_Theta,c_Pi})
                        puncture_max=std::max(puncture_max,std::abs(u(iv,n)));
                int k=0;
                for(int n:{c_Gamma1,c_Gamma2,c_lapse,c_shift1,c_shift2})
                {if(first)initial[it()](iv,k)=u(iv,n);++k;}
                const double r2=x*x+y*y;
                if(r2>=.00075*.00075 && r2<=.0025*.0025)
                {
                    gamma_delta_W=std::max(gamma_delta_W,std::hypot(u(iv,c_Gamma1)-initial[it()](iv,0),u(iv,c_Gamma2)-initial[it()](iv,1)));
                    lapse_delta_W=std::max(lapse_delta_W,std::abs(a-initial[it()](iv,2)));
                    shift_delta_W=std::max(shift_delta_W,std::hypot(u(iv,c_shift1)-initial[it()](iv,3),u(iv,c_shift2)-initial[it()](iv,4)));
                }
            }
        }
    }
    void flush()
    {
        std::ofstream out(m_p.data_path+"t18-monitor-l"+std::to_string(m_level)+".csv");
        out<<"level,time_M,h,dt,samples,nonfinite_or_bad_metric,chi_floor_crossings,lapse_floor_crossings,chi_min,lapse_min,lapse_max,max_beta,max_light_speed,max_lapse_speed,max_shift_speed,puncture_field_max\n";
        out<<std::setprecision(17)<<m_level<<','<<m_time<<','<<m_dx<<','<<m_dt<<','<<samples<<','<<bad<<','
           <<chi_floors<<','<<lapse_floors<<','<<chi_min<<','<<alpha_min<<','<<alpha_max<<','
           <<beta_max<<','<<light<<','<<lapse<<','<<gauge<<','<<puncture_max<<'\n';
        if(!out)MayDay::Error("T18 monitor write failed");
    }
    Real advance() override
    {if(m_level==0)started=std::chrono::steady_clock::now();return GRAMRLevel::advance();}
    void specificEvalRHS(GRLevelData &u,GRLevelData &rhs,double time) override
    {if(monitor)census(u);EMSBH2DLevel::specificEvalRHS(u,rhs,time);}
    void postTimeStep() override
    {
        GRAMRLevel::postTimeStep();
        if(monitor && m_level==m_p.max_level)
        {
            census(m_state_new);
            std::ofstream history(m_p.data_path+"t18-launch-history.csv",std::ios::app);
            if(history.tellp()==0)history<<"time_M,delta_Gamma_W_max,delta_lapse_W_max,delta_shift_W_max,puncture_field_max\n";
            history<<std::setprecision(17)<<m_time<<','<<gamma_delta_W<<','<<lapse_delta_W<<','<<shift_delta_W<<','<<puncture_max<<'\n';
            if(!history)MayDay::Error("T18 launch history write failed");
            if(bad||chi_floors||lapse_floors||puncture_max>1e6||beta_max>2.||alpha_max>2.)
            {flush();MayDay::Error("T18 predeclared boundedness screen failed");}
            if(m_time+m_dt>stop+1e-14){flush();throw T18Stop{};}
            if(samples%160==0)flush();
        }
        if(m_level==0)
        {
            std::ofstream out(m_p.data_path+"t18-rates.csv",std::ios::app);
            if(out.tellp()==0)out<<"time_M,seconds\n";
            out<<std::setprecision(17)<<m_time<<','<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<'\n';
        }
    }
};
class T18Factory : public DefaultLevelFactory<EMSBH2DLevel>
{
  public:
    using DefaultLevelFactory<EMSBH2DLevel>::DefaultLevelFactory;
    AMRLevel *new_amrlevel() const override
    {auto *l=new T18Level(m_gr_amr,m_p,m_p.verbosity);l->initialDtMultiplier(m_p.dt_multiplier);return l;}
};

int main(int argc,char **argv)
{
    mainSetup(argc,argv);GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);SimulationParameters p(pp);
    bool fixed=false;pp.load("t18_fixed_hierarchy",fixed,false);
    BHAMR amr;T18Factory factory(amr,p);
    if(fixed)
    {
        // Native setup with its new-run call suppressed; then use Chombo's
        // supported fixed-hierarchy initializer with exactly matched cell sets.
        ProblemDomain domain(IntVect::Zero,p.ivN);
        for(int d=0;d<CH_SPACEDIM;++d)domain.setPeriodic(d,p.boundary_params.is_periodic[d]);
        amr.define(p.max_level,p.ref_ratios,domain,&factory);
        amr.gridBufferSize(4);amr.blockFactor(p.block_factor);amr.maxGridSize(p.max_grid_size);
        amr.regridIntervals(p.regrid_interval);amr.fillRatio(p.fill_ratio);amr.verbosity(p.verbosity);
        amr.checkpointInterval(p.checkpoint_interval);amr.checkpointPrefix(p.hdf5_path+p.checkpoint_prefix);
        amr.plotInterval(p.plot_interval);amr.plotPrefix(p.hdf5_path+p.plot_prefix);
        Vector<Vector<Box>> grids(p.max_level+1);amr.makeBaseLevelMesh(grids[0]);
        for(int l=1;l<=p.max_level;++l)
        {
            const double h=p.coarsest_dx/std::pow(2,l);
            for(double sign:{-1.,1.})
            {
                const double centre=p.emsbh_params.star_centre[0]+sign*p.emsbh_params.separation/2;
                const int x0=(int(std::floor(centre/h/8))-2)*8;
                for(int j=0;j<16;j+=p.max_grid_size)for(int i=0;i<32;i+=p.max_grid_size)
                    grids[l].push_back(Box(IntVect(x0+i,j),IntVect(x0+std::min(i+p.max_grid_size,32)-1,std::min(j+p.max_grid_size,16)-1)));
            }
        }
        amr.setupForFixedHierarchyRun(grids);
    }
    else setupAMRObject(amr,factory);
    AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,p.verbosity);
    amr.set_interpolator(&interp);
    const auto levels=amr.getAMRLevels();
    for(int i=0;i<levels.size();++i)
    {auto *l=dynamic_cast<T18Level *>(levels[i]);if(l->monitor){l->census(l->getLevelData());l->flush();}}
    bool stopped=false;
    try{amr.run(p.stop_time,p.max_steps);}catch(const T18Stop &){stopped=true;}
    if(stopped)
    {
        for(int i=0;i<levels.size();++i)dynamic_cast<T18Level *>(levels[i])->flush();
        auto *l=dynamic_cast<T18Level *>(levels[levels.size()-1]);
        std::ofstream out(p.data_path+"t18-stop.txt");out<<std::setprecision(17)<<l->m_time<<'\n';
    }
    else amr.conclude();
    mainFinalize();return 0;
}
