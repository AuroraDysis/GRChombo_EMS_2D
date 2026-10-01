// T12: initialization-only reader/setter audit. Never used after t=0.
#include "EMSBH_trumpet_read.hpp"
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

static void hash_value(uint64_t &h, double x)
{
    unsigned char bytes[8]; std::memcpy(bytes, &x, 8);
    for (auto b : bytes) { h ^= b; h *= 1099511628211ULL; }
}
int main(int argc, char **argv)
{
    if (argc != 6) return 2;
    const std::string action=argv[1], mode=argv[3];
    EMSBH_params_t p{}; p.bh_mass=1.; p.star_centre={0.,0.}; p.data_path=argv[2];
#ifndef T12_BASELINE
    p.use_maximal_initial_lapse=(mode!="original");
    if(mode=="boost") p.rapidity=.1;
    if(mode=="boost_flag") p.boosted=true;
    if(mode=="binary") p.binary=true;
    if(mode=="separation") p.separation=1.;
    if(mode=="ctt") p.ctt_data_path="unsupported.ctt";
#else
    if(mode!="original") return 2;
#endif
    EMSBH_trumpet_read setter(p,{4*M_PI,0.,0.,-.8},1.,1.,0);
    auto *old=std::cout.rdbuf(std::cerr.rdbuf()); setter.compute_1d_solution(); std::cout.rdbuf(old);
    if(action=="--guards") {
        if(mode=="helper_boost") setter.compute_single_bh_vars(.1,.1,1.,0.,.1);
        if(mode=="helper_binary") setter.compute_binary_bh_vars(.1,.1,1.,1.,0.);
        return 0;
    }
    if(action=="--points") {
        std::ifstream in(argv[4],std::ios::binary); std::ofstream out(argv[5],std::ios::binary); double xy[2];
        while(in.read(reinterpret_cast<char*>(xy),16)) {
            auto v=setter.compute_single_bh_vars(xy[0],xy[1],1.,0.,0.);
            const auto q=setter.m_1d_sol.compute_radial_vars(std::sqrt(xy[0]*xy[0]+xy[1]*xy[1]));
            double a[32];v.enum_mapping([&](int c,double x){a[c]=x;});
            a[28]=q.lapse;a[29]=q.X;a[30]=setter.m_1d_sol.get_chebyshev_value(0,q.s);a[31]=q.dr_ds_numerator;
            out.write(reinterpret_cast<char*>(a),sizeof(a));
        }
        return in.eof() && out.good()?0:3;
    }
    if(action!="--census") return 2;
    std::ifstream in(argv[4]); std::ofstream out(argv[5]);
    out << "level,box,h_M,points,min_R,min_Y,min_G,min_alpha_K,max_alpha_K,min_chi,min_lapse,alpha_floor_margin,chi_floor_margin,driver_cancel_max,state_hash,nonlapse_hash,finite\n" << std::setprecision(17);
    int level,box,x0,y0,x1,y1;double h;
    while(in>>level>>box>>x0>>y0>>x1>>y1>>h) {
        double minR=1e300,minY=1e300,minG=1e300,minAK=1e300,maxAK=0,minChi=1e300,minL=1e300,cancel=0;
        uint64_t full=1469598103934665603ULL,nonlapse=full;long count=0;
        for(int j=y0-3;j<=y1+3;++j) for(int i=x0-3;i<=x1+3;++i) {
            const double x=(i+.5)*h-336.,y=(j+.5)*h,R=std::sqrt(x*x+y*y);
            const auto q=setter.m_1d_sol.compute_radial_vars(R);
            const double Y=setter.m_1d_sol.get_chebyshev_value(0,q.s);
            auto v=setter.compute_single_bh_vars(x,y,1.,0.,0.);
            // Gauge initialization after Gamma calculation has the same dependency in both options.
            for(int k=0;k<2;++k){v.B[k]=.75*v.Gamma[k]-v.shift[k];cancel=std::max(cancel,std::abs(.75*v.Gamma[k]-v.shift[k]-v.B[k]));}
            if(!(Y>0 && q.dr_ds_numerator>0 && q.lapse>0 && std::isfinite(q.lapse) && v.chi>1e-12 && v.lapse>1e-12)) return 4;
            bool finite=true;
            v.enum_mapping([&](int c,double z){finite &= std::isfinite(z);hash_value(full,z);if(c!=c_lapse)hash_value(nonlapse,z);});
            if(!finite) return 5;
            minR=std::min(minR,R);minY=std::min(minY,Y);minG=std::min(minG,q.dr_ds_numerator);
            minAK=std::min(minAK,q.lapse);maxAK=std::max(maxAK,q.lapse);minChi=std::min(minChi,v.chi);minL=std::min(minL,v.lapse);++count;
        }
        out<<level<<','<<box<<','<<h<<','<<count<<','<<minR<<','<<minY<<','<<minG<<','<<minAK<<','<<maxAK<<','<<minChi<<','<<minL<<','<<minL-1e-12<<','<<minChi-1e-12<<','<<cancel<<','<<full<<','<<nonlapse<<",1\n";
    }
    return in.eof() && out.good()?0:3;
}
