// Tests-only factory: the qualified T14 dense representation, unchanged
// production initialData, and native point audits on its actual hierarchy.
#include "T14DenseTags.hpp"
#include "SetupFunctions.hpp"
#include "t22-audit.hpp"
#include "TraceARemovalCartoon.hpp"
#include "PositiveChiAndAlpha.hpp"
#include <vector>

class T22AuditLevel : public T14DenseTagsLevel
{
  public:
    using T14DenseTagsLevel::T14DenseTagsLevel;
    void dump()
    {
        GRParmParse pp;std::vector<int> levels;std::vector<double> radii;
        int puncture_cells;pp.load("t22_audit_levels",levels,pp.countval("t22_audit_levels"));
        pp.load("t22_audit_radii",radii,pp.countval("t22_audit_radii"));
        pp.load("t22_audit_puncture_cells",puncture_cells);
        std::ofstream census("t22-boxes-L"+std::to_string(m_level)+".csv");
        census<<"level,box,lo_i,lo_j,hi_i,hi_j,h_M,dt_M\n"<<std::setprecision(17);
        int box=0;long long finite=0,nonfinite=0,nonzero_KTheta=0,nonzero_B=0;
        double alpha_min=INFINITY,alpha_max=-INFINITY,beta_max=0.,light_max=0.,lapse_speed_max=0.,driver_speed_max=0.;
        for(DataIterator it=m_state_new.dataIterator();it.ok();++it,++box)
        {
            const auto b=m_state_new.disjointBoxLayout()[it()];
            census<<m_level<<','<<box<<','<<b.smallEnd(0)<<','<<b.smallEnd(1)<<','<<b.bigEnd(0)<<','<<b.bigEnd(1)<<','<<m_dx<<','<<m_dt<<'\n';
            const auto &u=m_state_new[it()];
            for(BoxIterator iv(b);iv.ok();++iv)
            {
                for(int c=0;c<NUM_VARS;++c) {finite+=std::isfinite(u(iv(),c));nonfinite+=!std::isfinite(u(iv(),c));}
                nonzero_KTheta+=u(iv(),c_K)!=0. || u(iv(),c_Theta)!=0.;
                nonzero_B+=u(iv(),c_B1)!=0. || u(iv(),c_B2)!=0.;
                const double a=u(iv(),c_lapse),chi=u(iv(),c_chi),h11=u(iv(),c_h11),h22=u(iv(),c_h22),h12=u(iv(),c_h12);
                alpha_min=std::min(alpha_min,a);alpha_max=std::max(alpha_max,a);
                double beta=std::hypot(u(iv(),c_shift1),u(iv(),c_shift2));beta_max=std::max(beta_max,beta);
                const double small=.5*(h11+h22-std::hypot(h11-h22,2*h12));
                if(!(small>0.) || !(a>0.) || !(chi>0.))throw std::runtime_error("T22 invalid initial geometry");
                const double inv=chi/small;
                light_max=std::max(light_max,beta+a*std::sqrt(inv));
                // Conservative frozen-coefficient envelopes, not a principal
                // symbol / strong-hyperbolicity proof on this curved state.
                lapse_speed_max=std::max(lapse_speed_max,beta+std::sqrt(m_p.ccz4_params.lapse_coeff*std::pow(a,m_p.ccz4_params.lapse_power)*inv));
                driver_speed_max=std::max(driver_speed_max,beta+std::sqrt((4./3.)*m_p.ccz4_params.shift_Gamma_coeff/small));
            }
        }
        std::ofstream totals("t22-state-L"+std::to_string(m_level)+".csv");
        totals<<"level,valid_values,nonfinite,nonzero_KTheta_cells,nonzero_B_cells,lapse_min,lapse_max,shift_max,light_envelope,lapse_envelope,driver_envelope,Courant_envelope\n"<<std::setprecision(17)
            <<m_level<<','<<finite<<','<<nonfinite<<','<<nonzero_KTheta<<','<<nonzero_B<<','<<alpha_min<<','<<alpha_max<<','<<beta_max<<','<<light_max<<','<<lapse_speed_max<<','<<driver_speed_max<<','<<m_dt/m_dx*std::max({light_max,lapse_speed_max,driver_speed_max})<<'\n';
        if(nonfinite || !census || !totals)throw std::runtime_error("T22 invalid initial census");
        if(std::find(levels.begin(),levels.end(),m_level)==levels.end())return;
        // Avoid a hierarchy-wide scratch duplicate: one native box at a time.
        std::ofstream native("t22-native-L"+std::to_string(m_level)+".bin",std::ios::binary);
        T22NativeAudit<MovingPunctureGauge> moving(m_p.ccz4_params,m_dx,m_p.sigma,CouplingFunction(m_p.coupling_function_params),m_p.m_G_Newton,m_p.formulation);
        T22NativeAudit<ExperimentalGauge> experimental(m_p.ccz4_params,m_dx,m_p.sigma,CouplingFunction(m_p.coupling_function_params),m_p.m_G_Newton,m_p.formulation);
        const int centre=std::llround(m_p.center[0]/m_dx);
        std::vector<IntVect> points;
        for(int j=0;j<puncture_cells;++j)for(int i=-puncture_cells;i<puncture_cells;++i)points.emplace_back(D_DECL(centre+i,j,0));
        for(double r:radii)for(double angle:std::array<double,3>{0.,M_PI/2,M_PI/4})
            points.emplace_back(D_DECL(std::llround((m_p.center[0]+r*std::cos(angle))/m_dx-.5),std::max(0,int(std::llround(r*std::sin(angle)/m_dx-.5))),0));
        std::sort(points.begin(),points.end());points.erase(std::unique(points.begin(),points.end()),points.end());
        for(DataIterator it=m_state_new.dataIterator();it.ok();++it)
        {
            const auto valid=m_state_new.disjointBoxLayout()[it()];
            bool needed=false;for(const auto &iv:points)needed|=valid.contains(iv);if(!needed)continue;
            FArrayBox u(m_state_new[it()].box(),NUM_VARS),rhs(u.box(),NUM_VARS);u.copy(m_state_new[it()]);
            BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(),PositiveChiAndAlpha(m_p.min_chi,m_p.min_lapse)),u,u,u.box(),disable_simd());
            BoxPointers pointers(u,rhs);
            for(const auto &iv:points)if(valid.contains(iv))
            {
                double meta[5]={double(m_level),double(iv[0]),double(iv[1]),m_dx,m_time};
                std::vector<double> out(T22NativeAudit<MovingPunctureGauge>::columns);
                Cell<double> cell(iv,pointers);
                if(m_ems_gauge.differential())moving.sample(cell,(iv[1]+.5)*m_dx,out.data());
                else experimental.sample(cell,(iv[1]+.5)*m_dx,out.data());
                native.write(reinterpret_cast<char*>(meta),sizeof(meta));
                native.write(reinterpret_cast<char*>(out.data()),out.size()*sizeof(double));
            }
        }
        if(!native)throw std::runtime_error("T22 native audit write failed");
    }
};

int main(int argc,char **argv)
{
    mainSetup(argc,argv);int status=0;
    try
    {
        GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);SimulationParameters p(pp);
        EMSGaugeSelection(pp).record(p);
        BHAMR amr;DefaultLevelFactory<T22AuditLevel> factory(amr,p);setupAMRObject(amr,factory);
        for(auto *l:amr.get_gramrlevels())l->fillAllEvolutionGhosts();
        for(auto *l:amr.get_gramrlevels())dynamic_cast<T22AuditLevel*>(l)->dump();
        if(!p.restart_from_checkpoint)amr.writeCheckpointFile();
        pout()<<"T22_NATIVE_INITIAL_AUDIT_COMPLETE; no advances"<<std::endl;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';status=2;}
    mainFinalize();return status;
}
