#include "t22-audit.hpp"
#include <random>
#include <iostream>
#include <vector>

template<class gauge_t> void generic(const std::string &path)
{
    using Audit=T22NativeAudit<gauge_t>;
    std::ofstream out(path,std::ios::binary);
    std::mt19937_64 random(22009);
    const auto uniform=[&] {return std::generate_canonical<double,53>(random);};
    Box b(IntVect(D_DECL(-3,7,0)),IntVect(D_DECL(3,13,0)));
    FArrayBox u(b,NUM_VARS),rhs(b,NUM_VARS);
    std::vector<double> result(Audit::columns);
    for(int s=0;s<64;++s)
    {
        typename Audit::params_t p{};
        p.kappa1=.1;p.kappa2=0.;p.kappa3=1.;p.covariantZ4=true;
        p.lapse_advec_coeff=s<32?0.:uniform();
        p.lapse_coeff=s<32?2.:1.2+uniform();p.lapse_power=s<32?1.:.5+uniform();
        p.shift_advec_coeff=s<32?0.:uniform();
        p.shift_Gamma_coeff=s<32?1.:.5+uniform();p.eta=s<32?1.:.3+uniform();
        const double h=.03+.01*uniform();
        for(BoxIterator cell(b);cell.ok();++cell)
        {
            const double x=cell()[0]*h,y=(cell()[1]+.5)*h;
            for(int c=0;c<NUM_VARS;++c)u(cell(),c)=.02*std::sin((c+1)*x+y)+.01*std::cos(x-(c+1)*y);
            u(cell(),c_chi)=.8+.02*x;
            u(cell(),c_h11)=1.03+.01*x;u(cell(),c_h22)=.97+.01*y;u(cell(),c_hww)=1.01+.01*y;
            u(cell(),c_h12)=.002*x;u(cell(),c_lapse)=.7+.02*y;
            u(cell(),c_K)=.02+.003*x;u(cell(),c_Theta)=.001+.0002*y;
            u(cell(),c_phi)=.05+.04*x+.03*y;u(cell(),c_Pi)=.17+.01*y;
            u(cell(),c_Ex)=.02+.01*x;u(cell(),c_Ey)=.03+.01*y;u(cell(),c_Ez)=.012;
            u(cell(),c_Bx)=.008;u(cell(),c_By)=.009;u(cell(),c_Bz)=.015;
        }
        Audit kernel(p,h,1.,CouplingFunction({4*M_PI,0.,0.,-.8}),1.);
        BoxPointers pointers(u,rhs);Cell<double> cell(IntVect(D_DECL(0,10,0)),pointers);
        kernel.sample(cell,10.5*h,result.data());
        out.write(reinterpret_cast<const char *>(result.data()),result.size()*sizeof(double));
    }
    if(!out) throw std::runtime_error("T22 generic output write failed");
}
int main(int argc,char **argv)
{
    if(argc!=3)return 2;
    if(std::string(argv[1])=="experimental")generic<ExperimentalGauge>(argv[2]);
    else if(std::string(argv[1])=="moving_puncture")generic<MovingPunctureGauge>(argv[2]);
    else return 2;
    return 0;
}
