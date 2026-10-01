// Current-state replay only. The T11 native algebra is shared, not re-derived.
#define main t11_t0_reader_main
#include "T11RHS.cpp"
#undef main
#include "GammaCartoonCalculator.hpp"
#include "ConstraintsCartoon.hpp"
#include "EMSCartoonGaussConstraints.hpp"

int main(int argc,char **argv)
{
    if (argc!=4) return 2;
    const std::string mode=argv[1];
    if (mode!="--stage" && mode!="--state") return 2;
    std::ifstream in(argv[2],std::ios::binary);
    std::ofstream out(argv[3],std::ios::binary);
    Base::params_t p{};p.kappa1=.1;p.kappa2=0;p.kappa3=1;p.covariantZ4=true;
    p.lapse_advec_coeff=1;p.shift_advec_coeff=1;p.shift_Gamma_coeff=.75;p.eta=1;
    CouplingFunction coupling({4*M_PI,0.,0.,-.9});
    std::vector<double> record(3+49*NUM_VARS);
    while (in.read(reinterpret_cast<char*>(record.data()),record.size()*8))
    {
        double h=record[0];int j=std::llround(record[1]/h-.5);
        IntVect c(D_DECL(0,j,0));Box core(c,c),b(c-3*IntVect::Unit,c+3*IntVect::Unit);
        FArrayBox u(b,NUM_VARS),rhs(core,NUM_VARS),gamma(core,NUM_VARS),con(core,NUM_DIAGNOSTIC_VARS);
        std::copy(record.begin()+3,record.end(),u.dataPtr());
        if (mode=="--stage")
            BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(),PositiveChiAndAlpha(1e-12,1e-12)),u,u,b,disable_simd());
        Audit kernel(p,h,1.,coupling,1.);double result[Audit::columns];
        BoxPointers pointers(u,rhs);Cell<double> cell(c,pointers);
        kernel.native(cell,record[1],result);
        double total[NUM_VARS];for (int n=0;n<NUM_VARS;++n) total[n]=rhs(c,n);
        kernel.Base::compute(cell);
        for (int n=0;n<NUM_VARS;++n)
        {double v=rhs(c,n);if (std::memcmp(&v,&total[n],8)!=0) result[Audit::columns-1]+=1.;}
        BoxLoops::loop(GammaCartoonCalculator(h),u,gamma,core,disable_simd());
        BoxLoops::loop(make_compute_pack(Constraints<CouplingFunction>(h,coupling,1.),
                       EMSCartoonGaussConstraints(h,{4*M_PI,0.,0.,-.9})),u,con,core,disable_simd());
        out.write(reinterpret_cast<char*>(result),sizeof(result));
        const double extra[]={gamma(c,c_Gamma1),gamma(c,c_Gamma2),con(c,c_Ham),
                 con(c,c_Mom1),con(c,c_Mom2),con(c,c_GaussE),con(c,c_GaussB)};
        out.write(reinterpret_cast<const char*>(extra),sizeof(extra));
        for (int n=0;n<NUM_VARS;++n)
        {double v=rhs(c,n);out.write(reinterpret_cast<const char*>(&v),8);}
    }
    return in.eof() && out?0:3;
}
