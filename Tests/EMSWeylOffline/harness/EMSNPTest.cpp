#include "EMSNPScalar.hpp"
#include "SetupFunctions.hpp"
#include <cassert>

int main(int argc,char **argv)
{
    mainSetup(argc,argv);
    EMSNPScalar compute({0.,0.},1.,CouplingFunction::params_t{},1.);
    EMSNPScalar::Vars<double> v{};
    EMSNPScalar::Vars<Tensor<1,double>> d1{};
    EMSNPScalar::Diff2Vars<Tensor<2,double>> d2{};
    VarsTools::assign(v,0.); VarsTools::assign(d1,0.); VarsTools::assign(d2,0.);
    v.chi=v.h[0][0]=v.h[1][1]=v.hww=v.lapse=1.;
    Coordinates<double> coords(IntVect::Zero,1.); coords.x=3.; coords.y=4.;
    auto g=compute.geometry(v,d1,d2,coords);
    for (const auto value:EMSNPScalar::evaluate(g,g.tetrad).values)
        if (std::abs(value)>1e-14) throw std::runtime_error("NP flat-space control failed");
    v.Ex=.3;v.Ey=.4;v.Ez=.2;v.Bx=.1;v.By=.2;v.Bz=.3;
    g=compute.geometry(v,d1,d2,coords);
    const auto s=EMSNPScalar::evaluate(g,g.tetrad);
    const std::array<EMSNPScalar::Complex,3> expected{{{.15,-.12},{.25,.11},{.15,-.08}}};
    for (int i=0;i<3;++i)
        if (std::abs(s.values[5+i]-expected[i])>1e-14)
            throw std::runtime_error("NP complex Maxwell normalization failed");
    bool rejected=false;
    v.chi=1e-5;
    try { compute.geometry(v,d1,d2,coords); } catch (const std::runtime_error &) { rejected=true; }
    if (!rejected) throw std::runtime_error("NP inconsistent metric floor was not rejected");
    pout() << "PASS NP flat, complex Maxwell signs/normalization, chi-floor rejection\n";
    mainFinalize(); return 0;
}
