// Write-only extension: native tags, native advance, native enforcement/RHS.
#include "t17-dense.hpp"
#include "t16-native.hpp"
#include <chrono>
class T17Level : public T14DenseTagsLevel
{
    bool dense=false;
    std::chrono::steady_clock::time_point coarse_start;
  public:
    T17Level(GRAMR &g,const SimulationParameters &p,int v):T14DenseTagsLevel(g,p,v)
    { GRParmParse pp;pp.load("t14_dense_initial_tags",dense,false); }
    void tagCells(IntVectSet &tags) override
    { if(dense)T14DenseTagsLevel::tagCells(tags);else GRAMRLevel::tagCells(tags); }
    Real advance() override
    {
        if(m_level==0)coarse_start=std::chrono::steady_clock::now();
        return GRAMRLevel::advance();
    }
    void twin_record(int phase,const GRLevelData &u,double time,int growth=3)
    {
        if(!m_t13.enabled() || m_level!=m_p.max_level)return;
        if(m_p.emsbh_params.binary)
            for(double sign:{-1.,1.})m_t13.record(phase,m_level,u,time,m_dx,m_dt,
                m_p.emsbh_params.star_centre[0]+sign*m_p.emsbh_params.separation/2,growth);
        else m_t13.record(phase,m_level,u,time,m_dx,m_dt,m_p.emsbh_params.star_centre[0],growth);
    }
    void postTimeStep() override
    {
        try { GRAMRLevel::postTimeStep(); }
        catch(const T13LaunchStop &){twin_record(50,m_state_new,m_time);throw;}
        if(m_time>0.)twin_record(50,m_state_new,m_time);
        if(m_level==0&&m_time>0.)
        {
            const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-coarse_start).count();
            std::ofstream out(m_p.data_path+"t17-coarse-rates.csv",std::ios::app);
            if(out.tellp()==0)out<<"time_M,dt0_M,seconds\n";
            out<<std::setprecision(17)<<m_time<<','<<m_dt<<','<<wall<<'\n';
            if(!out)MayDay::Error("T17 rate write failed");
        }
    }
    void specificEvalRHS(GRLevelData &u,GRLevelData &rhs,double time) override
    {
        const bool capture=m_t13.stage_capture(m_level,m_p.max_level);
        if(capture)twin_record(2,u,time);
        EMSBH2DLevel::specificEvalRHS(u,rhs,time);
        if(capture){twin_record(15,u,time);twin_record(3,rhs,time,0);}
    }
};
template<class level_t> class T17Factory : public DefaultLevelFactory<level_t>
{
  public:
    T17Factory(GRAMR &g,SimulationParameters &p):DefaultLevelFactory<level_t>(g,p){}
    AMRLevel *new_amrlevel() const override
    {
        auto *level=new T17Level(this->m_gr_amr,this->m_p,this->m_p.verbosity);
        level->initialDtMultiplier(this->m_p.dt_multiplier);return level;
    }
};
inline void t17_initial_audit(BHAMR &amr,const SimulationParameters &p)
{
    std::ofstream source(p.data_path+"t17-source.csv"),boxes(p.data_path+"t17-boxes.csv");
    source<<"hole,level,x,y,h,physical_Gamma_n,KO_Gamma_n,total_Gamma_n,Base_bit_mismatches\n"<<std::setprecision(17);
    boxes<<"level,box,x0,y0,x1,y1,h,valid_cells\n"<<std::setprecision(17);
    for(int l=0;l<amr.getAMRLevels().size();++l)
    {
        auto *level=dynamic_cast<T17Level *>(amr.getAMRLevels()[l]);
        const auto &data=level->getLevelData();const double h=level->get_dx();int n=0;
        for(DataIterator it=data.dataIterator();it.ok();++it,++n)
        {
            const Box valid=data.disjointBoxLayout()[it()];const auto &u=data[it()];
            boxes<<l<<','<<n<<','<<valid.smallEnd(0)<<','<<valid.smallEnd(1)<<','<<valid.bigEnd(0)<<','<<valid.bigEnd(1)<<','<<h<<','<<valid.numPts()<<'\n';
            if(l!=p.max_level)continue;
            FArrayBox rhs(valid,NUM_VARS);
            Audit kernel(p.ccz4_params,h,p.sigma,CouplingFunction(p.coupling_function_params),p.m_G_Newton);
            for(double sign:{-1.,1.})
            {
                if(!p.emsbh_params.binary&&sign>0)continue;
                const double center=p.emsbh_params.star_centre[0]+(p.emsbh_params.binary?sign*p.emsbh_params.separation/2:0.);
                const int c=std::floor(center/h-.5);
                for(int j=0;j<2;++j)for(int i=c;i<=c+1;++i)
                {
                    IntVect iv(D_DECL(i,j,0));if(!valid.contains(iv))continue;
                    BoxPointers pointers(u,rhs);Cell<double> cell(iv,pointers);double result[Audit::columns];
                    kernel.native(cell,(j+.5)*h,result);double total[NUM_VARS];
                    for(int k=0;k<NUM_VARS;++k)total[k]=rhs(iv,k);
                    kernel.Base::compute(cell);int mismatch=0;
                    for(int k=0;k<NUM_VARS;++k){double z=rhs(iv,k);mismatch+=std::memcmp(&z,&total[k],8)!=0;}
                    if(mismatch)MayDay::Error("T17 audit differs from Base");
                    const double physical=-sign*result[NUM_VARS+c_Gamma1],ko=-sign*result[2*NUM_VARS+c_Gamma1];
                    source<<(sign<0?"left":"right")<<','<<l<<','<<(i+.5)*h-center<<','<<(j+.5)*h<<','<<h<<','<<physical<<','<<ko<<','<<physical+ko<<','<<mismatch<<'\n';
                }
            }
        }
        level->twin_record(50,data,0.);
    }
    if(!boxes||!source)MayDay::Error("T17 audit write failed");
    pout()<<"T17 native t=0 audit complete; evolution may now start"<<std::endl;
}
