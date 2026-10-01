// Initialization only: native equations/KO and independent translation reference.
#include "t16-native.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "GammaCartoonCalculator.hpp"
#include "ConstraintsCartoon.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include <iomanip>

int main(int argc,char **argv)
{
    if(argc!=5) return 2;
    EMSBH_params_t p{};p.bh_mass=1.;p.star_centre={0.,0.};p.data_path=argv[1];
    p.boosted=true;p.rapidity=std::atanh(.0523381404947);
    p.use_geometric_initial_lapse=std::string(argv[2])=="geometric";
    EMSBH_trumpet_read setter(p,{4*M_PI,0.,0.,-.9},1.,1.,0);
    setter.compute_1d_solution();
    std::ifstream in(argv[3]);std::ofstream out(argv[4]);
    out<<"h,r_requested,ray,domain,component,x,y,state,physical_rhs,ko,translation,physical_residual,total_residual,constraint\n"<<std::setprecision(17);
    Base::params_t gp{};gp.kappa1=.1;gp.kappa2=0;gp.kappa3=1;gp.covariantZ4=true;
    gp.lapse_advec_coeff=1;gp.shift_advec_coeff=1;gp.shift_Gamma_coeff=.75;gp.eta=1;
    double h,r;int ray,domain;
    while(in>>h>>r>>ray>>domain)
    {
        const double x=ray?r/std::sqrt(2.):r;
        const int jy=ray?std::llround(r/std::sqrt(2.)/h-.5):0;
        const double y=(jy+.5)*h;
        const IntVect c(D_DECL(0,jy,0));const Box core(c,c),box(c-6*IntVect::Unit,c+6*IntVect::Unit);
        FArrayBox u(box,NUM_VARS),rhs(core,NUM_VARS),con(core,NUM_DIAGNOSTIC_VARS);
        for(BoxIterator it(box);it.ok();++it)
        {
            auto v=setter.compute_single_bh_vars(x+it()[0]*h,(it()[1]+.5)*h,1.,0.,p.rapidity);
            v.enum_mapping([&](int n,double z){u(it(),n)=z;});
        }
        const Box inner(c-4*IntVect::Unit,c+4*IntVect::Unit);
        BoxLoops::loop(GammaCartoonCalculator(h),u,u,inner,disable_simd());
        for(BoxIterator it(inner);it.ok();++it)for(int k=0;k<2;++k)
            u(it(),c_B1+k)=.75*u(it(),c_Gamma1+k)-u(it(),c_shift1+k);
        BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(),PositiveChiAndAlpha(1e-12,1e-12)),u,u,box,disable_simd());
        Audit kernel(gp,h,1.,CouplingFunction({4*M_PI,0.,0.,-.9}),1.);
        double result[Audit::columns];BoxPointers pointers(u,rhs);Cell<double> cell(c,pointers);
        kernel.native(cell,y,result);
        double total[NUM_VARS];for(int n=0;n<NUM_VARS;++n)total[n]=rhs(c,n);
        kernel.Base::compute(cell);
        for(int n=0;n<NUM_VARS;++n)
        {double actual=rhs(c,n);if(std::memcmp(&actual,&total[n],8)!=0)return 5;}
        BoxLoops::loop(make_compute_pack(Constraints<CouplingFunction>(h,CouplingFunction({4*M_PI,0.,0.,-.9}),1.),
            EMSCartoonGaussConstraints(h,{4*M_PI,0.,0.,-.9})),u,con,core,disable_simd());
        const double w[4]={4./5.,-1./5.,4./105.,-1./280.};
        for(int n=0;n<NUM_VARS;++n)
        {
            double dx=0.;for(int k=1;k<=4;++k)dx+=w[k-1]*(u(c+k*BASISV(0),n)-u(c-k*BASISV(0),n))/h;
            const double t=-.0523381404947*dx,physical=result[NUM_VARS+n],ko=result[2*NUM_VARS+n];
            double constraint=0.;
            if(n==c_K)constraint=con(c,c_Ham);
            if(n==c_A11)constraint=con(c,c_Mom1);
            if(n==c_A22)constraint=con(c,c_Mom2);
            if(n==c_Ex)constraint=con(c,c_GaussE);
            out<<h<<','<<r<<','<<ray<<','<<domain<<','<<n<<','<<x<<','<<y<<','<<u(c,n)<<','
               <<physical<<','<<ko<<','<<t<<','<<physical-t<<','<<physical+ko-t<<','<<constraint<<'\n';
            if(!std::isfinite(physical+ko-t))return 4;
        }
    }
    return in.eof()&&out?0:3;
}
