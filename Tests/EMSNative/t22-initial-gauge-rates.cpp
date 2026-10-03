// Offline reads of numerical t=0 blocks only. No initializer, static input or
// evolution. Native advection/gauge classes give all-valid-cell lapse rates;
// native complete RHS samples give the geometric/matter/KO driver split.
#include "t22-audit.hpp"
#include "TraceARemovalCartoon.hpp"
#include "PositiveChiAndAlpha.hpp"
#include <iostream>
#include <vector>
#include <set>

template<class data_t> struct T22LapseOnly
{
    data_t lapse;
    template<class mapper_t> void enum_mapping(mapper_t map)
    {VarsTools::define_enum_mapping(map,c_lapse,lapse);}
};
int main(int argc,char **argv)
{
    if(argc!=4)return 2;
    const bool moving=std::string(argv[1])=="moving_puncture";
    if(!moving && std::string(argv[1])!="experimental")return 2;
    std::ifstream in(argv[2],std::ios::binary);std::ofstream out(argv[3],std::ios::binary);
    std::ofstream defects(std::string(argv[3])+".mismatches.csv");
    defects << "level,box,i,j,component,reconstructed,native,reconstructed_bits,native_bits\n" << std::setprecision(17);
    std::ofstream aliases(std::string(argv[3])+".tensor-aliases.csv");
    aliases << "level,box,i,j,component,overwritten_upper,stored_lower\n" << std::setprecision(17);
    T22NativeAudit<MovingPunctureGauge>::params_t p{};
    p.kappa1=.1;p.kappa2=0.;p.kappa3=1.;p.covariantZ4=true;
    p.lapse_advec_coeff=moving?0.:1.;p.lapse_coeff=2.;p.lapse_power=1.;
    p.shift_advec_coeff=moving?0.:1.;p.shift_Gamma_coeff=moving?1.:.75;p.eta=1.;
    double info[9];
    while(in.read(reinterpret_cast<char*>(info),sizeof(info)))
    {
        const int l=int(info[0]),box=int(info[1]);const double h=info[2];
        const IntVect lo{D_DECL(int(info[3]),int(info[4]),0)},hi{D_DECL(int(info[5]),int(info[6]),0)};
        const double centre=info[7];const int count=int(info[8]);
        Box valid(lo,hi),allocated=valid;allocated.grow(3);
        FArrayBox u(allocated,NUM_VARS),rhs(allocated,NUM_VARS);
        if(!in.read(reinterpret_cast<char*>(u.dataPtr()),u.box().numPts()*NUM_VARS*sizeof(double)))return 3;
        std::vector<IntVect> points;
        for(int i=0;i<count;++i){int xy[2];if(!in.read(reinterpret_cast<char*>(xy),sizeof(xy)))return 3;points.emplace_back(D_DECL(xy[0],xy[1],0));}
        out.write(reinterpret_cast<char*>(info),sizeof(info));
        FourthOrderDerivatives deriv(h);
        ExperimentalGauge old(p);MovingPunctureGauge selected(p);
        // Rates are evaluated on the same saved state. The projection changes
        // A only; K/Theta/lapse/shift and their gauge rows are unchanged.
        BoxPointers pointers(u,rhs);
        for(BoxIterator bit(valid);bit.ok();++bit)
        {
            Cell<double> cell(bit(),pointers);
            auto v=cell.load_vars<CCZ4CartoonVars::VarsWithGauge>();
            auto al=deriv.advection<T22LapseOnly>(cell,v.shift);
            decltype(v) adv,rate;VarsTools::assign(adv,0.);VarsTools::assign(rate,0.);adv.lapse=al.lapse;
            CCZ4CartoonVars::VarsWithGauge<Tensor<1,double>> d1{};
            CCZ4CartoonVars::Diff2VarsWithGauge<Tensor<2,double>> d2{};
            if(moving)selected.rhs_gauge(rate,v,d1,d2,adv);else old.rhs_gauge(rate,v,d1,d2,adv);
            const double rows[2]={al.lapse,rate.lapse};
            out.write(reinterpret_cast<const char*>(rows),sizeof(rows));
        }
        BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(),PositiveChiAndAlpha(1e-12,1e-12)),u,u,u.box(),disable_simd());
        T22NativeAudit<ExperimentalGauge> aold(p,h,1.,CouplingFunction({4*M_PI,0.,0.,-.8}),1.);
        T22NativeAudit<MovingPunctureGauge> anew(p,h,1.,CouplingFunction({4*M_PI,0.,0.,-.8}),1.);
        for(const auto &iv:points)
        {
            std::vector<double> result(T22NativeAudit<MovingPunctureGauge>::columns);
            Cell<double> cell(iv,pointers);
            double alias[3]={-1.,0.,0.};
            if(moving)anew.sample(cell,(iv[1]+.5)*h,result.data(),alias);else aold.sample(cell,(iv[1]+.5)*h,result.data(),alias);
            if(alias[0]>=0)aliases << l << ',' << box << ',' << iv[0] << ',' << iv[1] << ',' << alias[0] << ',' << alias[1] << ',' << alias[2] << '\n';
            for(int c=0;c<NUM_VARS;++c)
            {
                const double reconstructed=result[3*NUM_VARS+c],actual=rhs(iv,c);
                std::uint64_t rb,ab;std::memcpy(&rb,&reconstructed,8);std::memcpy(&ab,&actual,8);
                if(rb!=ab)defects << l << ',' << box << ',' << iv[0] << ',' << iv[1] << ',' << c << ',' << reconstructed << ',' << actual << ',' << rb << ',' << ab << '\n';
            }
            out.write(reinterpret_cast<const char*>(result.data()),result.size()*sizeof(double));
        }
    }
    return in.eof() && out?0:3;
}
