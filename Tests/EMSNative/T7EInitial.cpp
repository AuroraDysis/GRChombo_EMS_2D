// Standalone t=0 diagnostic. No production path or evolution input is changed.
#include "EMSBH_trumpet_read.hpp"
#include <iostream>
#include <iomanip>
int main(int argc,char **argv)
{
    if (argc!=3 || (std::string(argv[2])!="--t0" && std::string(argv[2])!="--t0-fused")) return 2;
    EMSBH_params_t p{};p.data_path=argv[1];p.bh_mass=1.;p.star_centre={0.,0.};
    EMSBH_trumpet_read setter(p,{4*M_PI,0.,0.,-.8},1.,1.,0);
    auto *old=std::cout.rdbuf(std::cerr.rdbuf());
    setter.compute_1d_solution();std::cout.rdbuf(old);
    bool fused=std::string(argv[2])=="--t0-fused";double xy[3];
    while (std::cin.read(reinterpret_cast<char*>(xy),sizeof(double)*(fused?3:2)))
    {
        if(fused){const double h=xy[2];xy[0]=std::fma(xy[0]+.5,h,-224.);xy[1]=(xy[1]+.5)*h;}
        const double R=std::hypot(xy[0],xy[1]);const auto q=setter.m_1d_sol.compute_radial_vars(R);
        const double nu=setter.m_1d_sol.get_nu(),s=q.s,t=1-s,ks=q.compactification_factor;
        const double Y=setter.m_1d_sol.get_chebyshev_value(0,s),Y1=setter.m_1d_sol.get_chebyshev_value(0,s,1),Y2=setter.m_1d_sol.get_chebyshev_value(0,s,2);
        const double sr=nu*s*t/(R*ks),srr=sr*sr*(1/s-1/t-(nu-1)/ks)-sr/R;
        const double ls=2/(nu*s)-2*Y1/Y,lss=-2/(nu*s*s)-2*(Y2/Y-Y1*Y1/(Y*Y));
        const double chi=q.X*q.X,cr=chi*ls*sr,crr=chi*((ls*sr)*(ls*sr)+lss*sr*sr+ls*srr);
        const double phis=-setter.m_1d_sol.get_chebyshev_value(2,s)+t*setter.m_1d_sol.get_chebyshev_value(2,s,1),pr=phis*sr;
        const double F=std::exp(8*M_PI*.8*q.phi*q.phi);
        const double ricci=2*(crr+2*cr/R)-2.5*cr*cr/chi,aa=6*q.K_amplitude*q.K_amplitude;
        const double matter=16*M_PI*(q.Pi*q.Pi+chi*pr*pr+F*chi*q.E_R*q.E_R);
        auto v=setter.compute_single_bh_vars(xy[0],xy[1],1.,0.,0.);
        double result[35]={};
        v.enum_mapping([&](int c,double a){result[c]=a;});
        result[28]=ricci-aa-matter;result[29]=ricci;result[30]=aa;result[31]=matter;
        result[32]=cr;result[33]=crr;result[34]=q.s;
        std::fwrite(result,sizeof(double),35,stdout);
    }
    return std::cin.eof() ? 0 : 3;
}
