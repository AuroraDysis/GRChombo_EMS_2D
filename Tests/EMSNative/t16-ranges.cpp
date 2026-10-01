// Per-hole geometric lapse ranges on native grid footprints, t=0 only.
#include "EMSBH_trumpet_read.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char **argv)
{
    if(argc!=4)return 2;
    EMSBH_params_t p{};p.bh_mass=1.;p.star_centre={0.,0.};p.data_path=argv[1];
    p.boosted=true;p.rapidity=std::atanh(.0523381404947);p.use_geometric_initial_lapse=true;
    EMSBH_trumpet_read setter(p,{4*M_PI,0.,0.,-.9},1.,1.,0);setter.compute_1d_solution();
    const double c=std::cosh(p.rapidity),sh=std::sinh(p.rapidity);
    std::ifstream in(argv[2]);std::ofstream out(argv[3]);
    out<<"level,box,h,region,points,a_min,a_max,J_min\n"<<std::setprecision(17);
    int l,b,x0,y0,x1,y1;double h;
    while(in>>l>>b>>x0>>y0>>x1>>y1>>h)
    {
        // Each isolated grid footprint, and its per-hole value 32 M away.
        for(int side=0;side<3;++side)
        {
            double lo=INFINITY,hi=0.,jmin=INFINITY;long long count=0;
            for(int iy=y0-3;iy<=y1+3;++iy)for(int ix=x0-3;ix<=x1+3;++ix)
            {
                double x=(ix+.5)*h-336.+(side==1?32.:(side==2?-32.:0.));
                const double y=(iy+.5)*h,R=std::hypot(c*x,y);
                const auto q=setter.m_1d_sol.compute_radial_vars(R);
                const double w=c-(side==2?-sh:sh)*q.shift_R*c*x/R;
                const double J=w*w-sh*sh*q.lapse*q.lapse*q.X*q.X;
                const double a=q.lapse/std::sqrt(J);
                if(!(J>0&&w>0&&a>0&&a<=1.)||!std::isfinite(a))return 4;
                lo=std::min(lo,a);hi=std::max(hi,a);jmin=std::min(jmin,J);++count;
            }
            const char *label=side==0?"isolated-native":(side==1?"left-hole-near-right-footprint":"right-hole-near-left-footprint");
            out<<l<<','<<b<<','<<h<<','<<label<<','<<count<<','<<lo<<','<<hi<<','<<jmin<<'\n';
            std::cout<<"level="<<l<<" box="<<b<<" region="<<label<<" a_min="<<lo<<" a_max="<<hi<<'\n';
        }
    }
    return in.eof()&&out?0:3;
}
