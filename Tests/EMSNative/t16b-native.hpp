// Tests-only fixed binary hierarchy and write-only diagnostics. Evolution calls
// are qualified calls to the unmodified native level implementation.
#include "T14DenseTags.hpp"
#include "t16-native.hpp"
#include "ConstraintsCartoon.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include "GammaCartoonCalculator.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "BHAMR.hpp"
#include <set>

class T16BLevel : public EMSBH2DLevel
{
  public:
    T16BLevel(GRAMR &g,const SimulationParameters &p,int v):EMSBH2DLevel(g,p,v) {}
    void tagCells(IntVectSet &tags) override
    {
        Box bounds;
        for(DataIterator it=m_state_new.dataIterator();it.ok();++it)
            bounds.minBox(m_state_new.disjointBoxLayout()[it()]);
        DenseIntVectSet bits(bounds);bits.makeEmptyBits();IntVectSet chosen(bits);
        const int k=std::llround(std::ldexp(128.,-m_level-1)/m_dx);
        for(double offset:{-16.,16.})
        {
            const int c=std::llround((m_p.center[0]+offset)/m_dx);
            Box target(IntVect(D_DECL(c-k,0,0)),IntVect(D_DECL(c+k-1,k-1,0)));
            target&=m_problem_domain;chosen|=target;
        }
        chosen.recalcMinBox();tags=chosen;
        if(!tags.isDense()) MayDay::Error("T16b fixed tags unexpectedly became a tree");
    }
    void twin_record(int phase,const GRLevelData &u,double time,int growth=3)
    {
        if(m_level!=m_p.max_level)return;
        for(double offset:{-16.,16.})
            m_t13.record(phase,m_level,u,time,m_dx,m_dt,m_p.center[0]+offset,growth);
    }
    void postTimeStep() override
    {
        try { GRAMRLevel::postTimeStep(); }
        catch(const T13LaunchStop &) { twin_record(50,m_state_new,m_time);throw; }
        if(m_time>0.)twin_record(50,m_state_new,m_time);
    }
    void specificEvalRHS(GRLevelData &u,GRLevelData &rhs,const double time) override
    {
        const bool capture=m_t13.stage_capture(m_level,m_p.max_level);
        if(capture)twin_record(2,u,time);
        EMSBH2DLevel::specificEvalRHS(u,rhs,time);
        if(capture){twin_record(15,u,time);twin_record(3,rhs,time,0);}
    }
};
template<class level_t> class T16BFactory : public DefaultLevelFactory<level_t>
{
  public:
    T16BFactory(GRAMR &g,SimulationParameters &p):DefaultLevelFactory<level_t>(g,p) {}
    AMRLevel *new_amrlevel() const override
    {
        auto *level=new T16BLevel(this->m_gr_amr,this->m_p,this->m_p.verbosity);
        level->initialDtMultiplier(this->m_p.dt_multiplier);return level;
    }
};

// Read-only selector, copied from the native consumer's map/coverage tests.
// The audit checks it against forced-panel evaluation on sampled actual cells.
inline int t16b_panel(const EMSCTTSolution_read &c,double x,double y)
{
    const double eps=8*std::numeric_limits<double>::epsilon();
    for(int k=0;k<c.m_panels.size();++k)
    {
        const auto &p=c.m_panels[k];const double X=p.stretch*(x-p.center);
        const double r=std::sqrt(X*X+y*y),mu=X/r;
        if(!(r>0.)||mu<p.mu_lo-eps||mu>p.mu_hi+eps)continue;
        double hi=p.hi;
        if(p.kind=="bridge")
        {
            if(p.boundary=="plane")hi=-p.center*p.stretch/mu;
            else{const double a=1-mu*mu+mu*mu/(p.stretch*p.stretch),b=p.center*mu/p.stretch;
                 hi=(-b+std::sqrt(b*b+a*(p.hi*p.hi-p.center*p.center)))/a;}
        }
        const double t=p.kind=="inverse"?(2/r-p.lo-p.hi)/(p.hi-p.lo):
            (2*std::log(r)-std::log(p.lo)-std::log(hi))/(std::log(hi)-std::log(p.lo));
        if(std::isfinite(t)&&t>=-1-eps&&t<=1+eps)return k+1;
    }
    MayDay::Error("T16b selector outside companion");return 0;
}
inline double t16b_coarse_h(double x,double y)
{
    double h=1.;
    for(int l=1;l<=13;++l)
    {
        const double f=std::ldexp(128.,-l);
        if(std::min(std::abs(x+16.),std::abs(x-16.))<f&&std::abs(y)<f)h=std::ldexp(1.,-l);
    }
    return h;
}

inline void t16b_initial_audit(BHAMR &amr,const SimulationParameters &p)
{
    auto iso=p.emsbh_params;iso.binary=false;iso.ctt_data_path.clear();
    iso.use_geometric_initial_lapse=true;iso.use_maximal_initial_lapse=false;
    EMSBH_trumpet_read geometry(iso,p.coupling_function_params,p.m_G_Newton,p.dx[0],0);
    geometry.compute_1d_solution();
    EMSCTTSolution_read companion;std::string reason;
    const double eta=iso.rapidity;
    EMSCTTSolution_read::binding_t binding{{1.,1.},{-16.,16.},{eta,-eta},
        {4*M_PI,0.,0.,-.9}};
    if(!companion.check_file(p.emsbh_params.ctt_data_path,geometry.m_1d_sol,binding,reason))
        MayDay::Error(reason.c_str());
    std::ofstream state("t16b-nonlapse.bin",std::ios::binary),native("t16b-t0-native.csv"),
        ranges("t16b-ranges.csv"),source("t16b-source.csv"),limits("t16b-limits.csv");
    std::ifstream reference;GRParmParse pp;std::string reference_path;
    pp.load("t16b_compare_nonlapse",reference_path,std::string());
    if(!reference_path.empty()){reference.open(reference_path,std::ios::binary);if(!reference)MayDay::Error("missing default non-lapse reference");}
    native<<"level,x,y,h,weight,Ham,Mom,GaussE,GaussB,C_Gamma,panel,sheets,amr_seam,native_crossing\n";
    ranges<<"level,box,x0,y0,x1,y1,h,valid_cells,a_left_min,a_left_max,a_right_min,a_right_max,alpha_bin_min,alpha_bin_max,native_lapse_min,native_lapse_max,geometric_native_max_error,nonlapse_values,selector_checks\n";
    source<<"hole,level,x,y,h,physical_Gamma_n,KO_Gamma_n,total_Gamma_n,physical_K,KO_K,lapse_rhs,driver_cancel,Base_bit_mismatches";
    for(int k=0;k<13;++k)source<<",Gamma_term_"<<k;source<<'\n';
    limits<<"hole,R,a_own,a_companion,alpha_bin,alpha_bin_over_a_own\n";
    native<<std::setprecision(17);ranges<<std::setprecision(17);source<<std::setprecision(17);limits<<std::setprecision(17);
    for(double center:{-16.,16.})for(int k=2;k<=10;++k)
    {
        double r=std::pow(10.,-k),x=center+(center<0?1:-1)*r/std::cosh(eta);
        double a=geometry.compute_ems_adm_vars(x,0.,0.,1.,center,center<0?eta:-eta).lapse;
        double b=geometry.compute_ems_adm_vars(x,0.,0.,1.,-center,center<0?-eta:eta).lapse;
        double lo=std::min(a,b),hi=std::max(a,b),ab=lo/std::sqrt(1+(lo/hi)*(lo/hi)-lo*lo);
        limits<<(center<0?"left":"right")<<','<<r<<','<<a<<','<<b<<','<<ab<<','<<ab/a<<'\n';
    }
    const auto levels=amr.getAMRLevels();long long bit_count=0,selector_count=0;
    for(int l=0;l<levels.size();++l)
    {
        auto *level=dynamic_cast<T16BLevel *>(levels[l]);auto &data=level->m_state_new;
        const double h=level->get_dx();int box_index=0;
        for(DataIterator it=data.dataIterator();it.ok();++it,++box_index)
        {
            const auto &u=data[it()];const Box valid=data.disjointBoxLayout()[it()];
            FArrayBox gamma(valid,NUM_VARS),con(valid,NUM_DIAGNOSTIC_VARS),rhs(valid,NUM_VARS);
            BoxLoops::loop(GammaCartoonCalculator(h),u,gamma,valid,disable_simd());
            BoxLoops::loop(make_compute_pack(Constraints<CouplingFunction>(h,CouplingFunction(p.coupling_function_params),p.m_G_Newton),
                EMSCartoonGaussConstraints(h,p.coupling_function_params)),u,con,valid,disable_simd());
            double amin[3]={INFINITY,INFINITY,INFINITY},amax[3]={0.,0.,0.},nmin=INFINITY,nmax=0.,error=0.;
            long long values=0,checks=0;
            for(BoxIterator cell(u.box());cell.ok();++cell)
            {
                const auto iv=cell();double x=(iv[0]+.5)*h-p.center[0],y=(iv[1]+.5)*h;
                double a=geometry.compute_ems_adm_vars(x,y,0.,1.,-16.,eta).lapse;
                double b=geometry.compute_ems_adm_vars(x,y,0.,1.,16.,-eta).lapse;
                double lo=std::min(a,b),hi=std::max(a,b),ab=lo/std::sqrt(1+(lo/hi)*(lo/hi)-lo*lo),all[]={a,b,ab};
                for(int k=0;k<3;++k){if(!(all[k]>0&&all[k]<=1&&std::isfinite(all[k])))MayDay::Error("T16b lapse bound; no clipping");
                    amin[k]=std::min(amin[k],all[k]);amax[k]=std::max(amax[k],all[k]);}
                nmin=std::min(nmin,u(iv,c_lapse));nmax=std::max(nmax,u(iv,c_lapse));
                if(!(u(iv,c_lapse)>0.&&u(iv,c_lapse)<=1.&&std::isfinite(u(iv,c_lapse))))
                    MayDay::Error("T16b native lapse bound; no clipping");
                if(p.emsbh_params.use_geometric_initial_lapse)error=std::max(error,std::abs(u(iv,c_lapse)-ab));
                for(int n=0;n<NUM_VARS;++n)if(n!=c_lapse)
                {
                    double value=u(iv,n);if(!std::isfinite(value))MayDay::Error("T16b nonfinite initial state");
                    if(!reference.is_open())state.write(reinterpret_cast<const char*>(&value),8);++values;++bit_count;
                    if(reference.is_open()){double old;reference.read(reinterpret_cast<char*>(&old),8);
                        if(!reference||std::memcmp(&old,&value,8))MayDay::Error("T16b non-lapse fields differ");}
                }
                if(!valid.contains(iv))continue;
                bool covered=false;
                if(l+1<levels.size())
                {
                    const auto &fine_layout=dynamic_cast<T16BLevel *>(levels[l+1])->getLevelData().disjointBoxLayout();
                    for(LayoutIterator finer=fine_layout.layoutIterator();finer.ok();++finer)
                    {Box f=fine_layout[finer()];f.coarsen(2);if(f.contains(iv)){covered=true;break;}}
                }
                if(covered)continue;
                // Outside every declared volume mask, the only possible sheet
                // beyond r=128 is the compact chart's x=0 angular split.
                if(std::hypot(x,y)>128.&&std::abs(x)>3.)continue;
                const int panel=t16b_panel(companion,x,y);const double coarse=t16b_coarse_h(x,y);
                std::set<std::pair<int,int>> sheets;bool seam=false,crossing=false;
                for(int j=-3;j<=3;++j)for(int i=-3;i<=3;++i)
                {
                    double X=x+i*coarse,Y=y+j*coarse;int other=t16b_panel(companion,X,Y);
                    if(other!=panel)sheets.insert(std::minmax(panel,other));
                    if(t16b_coarse_h(X,Y)!=coarse)seam=true;
                    if(t16b_panel(companion,x+i*h,y+j*h)!=panel)crossing=true;
                }
                if(iv[0]%128==0&&iv[1]%32==0)
                {
                    const auto q=companion.evaluate(x,y,0.),q1=companion.evaluate(x,y,0.,panel);
                    if(q.logpsi!=q1.logpsi||q.C!=q1.C)MayDay::Error("T16b panel selector disagreement");++checks;++selector_count;
                }
                std::string label;for(auto pair:sheets){if(!label.empty())label+=';';label+=std::to_string(pair.first)+":"+std::to_string(pair.second);}
                double cg=std::hypot(u(iv,c_Gamma1)-gamma(iv,c_Gamma1),u(iv,c_Gamma2)-gamma(iv,c_Gamma2));
                if(std::hypot(x,y)<=128.||!sheets.empty())native<<l<<','<<x<<','<<y<<','<<h<<','<<y*h*h<<','<<con(iv,c_Ham)<<','
                    <<std::hypot(con(iv,c_Mom1),con(iv,c_Mom2))<<','<<con(iv,c_GaussE)<<','<<con(iv,c_GaussB)<<','<<cg<<','
                    <<panel<<','<<label<<','<<seam<<','<<crossing<<'\n';
                if(l==p.max_level)for(double center:{-16.,16.})
                    if(std::abs(x-center)<h&&y<2*h)
                    {
                        Base::params_t gauge{};gauge.kappa1=.1;gauge.kappa2=0;gauge.kappa3=1;gauge.covariantZ4=true;
                        gauge.lapse_advec_coeff=1;gauge.shift_advec_coeff=1;gauge.shift_Gamma_coeff=.75;gauge.eta=1;
                        Audit kernel(gauge,h,1.,CouplingFunction(p.coupling_function_params),p.m_G_Newton);
                        BoxPointers pointers(u,rhs);Cell<double> current(iv,pointers);double result[Audit::columns];kernel.native(current,y,result);
                        double total[NUM_VARS];for(int n=0;n<NUM_VARS;++n)total[n]=rhs(iv,n);
                        kernel.Base::compute(current);int mismatch=0;
                        for(int n=0;n<NUM_VARS;++n){double z=rhs(iv,n);mismatch+=std::memcmp(&z,&total[n],8)!=0;}
                        if(mismatch)MayDay::Error("T16b source replay differs from native Base");
                        double sign=center<0?1.:-1.,physical=sign*result[NUM_VARS+c_Gamma1],ko=sign*result[2*NUM_VARS+c_Gamma1];
                        source<<(center<0?"left":"right")<<','<<l<<','<<x<<','<<y<<','<<h<<','<<physical<<','<<ko<<','<<physical+ko<<','
                            <<result[NUM_VARS+c_K]<<','<<result[2*NUM_VARS+c_K]<<','<<result[NUM_VARS+c_lapse]<<','
                            <<std::max(std::abs(result[138]),std::abs(result[139]))<<','<<mismatch;
                        for(int k=0;k<13;++k)source<<','<<sign*result[4*NUM_VARS+2*k];source<<'\n';
                    }
            }
            ranges<<l<<','<<box_index<<','<<valid.smallEnd(0)<<','<<valid.smallEnd(1)<<','<<valid.bigEnd(0)<<','<<valid.bigEnd(1)<<','<<h<<','<<valid.numPts();
            for(int k=0;k<3;++k)ranges<<','<<amin[k]<<','<<amax[k];
            ranges<<','<<nmin<<','<<nmax<<','<<error<<','<<values<<','<<checks<<'\n';
            pout()<<"T16b ranges level "<<l<<" box "<<box_index<<" a_left="<<amin[0]<<":"<<amax[0]
                <<" a_right="<<amin[1]<<":"<<amax[1]<<" alpha_bin="<<amin[2]<<":"<<amax[2]<<std::endl;
        }
        level->twin_record(50,data,0.);
    }
    if(reference.is_open()&&reference.peek()!=EOF)MayDay::Error("T16b extra non-lapse reference values");
    if(!state||!native||!ranges||!source||!limits)MayDay::Error("T16b initial audit write failed");
    pout()<<"T16b native initial audit complete; non-lapse values="<<bit_count<<" selector checks="<<selector_count<<std::endl;
}
