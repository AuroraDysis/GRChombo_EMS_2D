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
    for (bool matter:{false,true})
    {
        const auto control=compute.geometry(v,d1,d2,coords,matter);
        const auto value=EMSNPScalar::evaluate(control,control.tetrad).values[4];
        if (std::abs(value.real()-control.legacy)>1e-13)
            throw std::runtime_error("NP vacuum/matter same-jet Weyl4 control failed");
    }
    // A general complex Weyl tensor exercises all ten degrees of freedom,
    // including components absent from the polar checkpoint sector.
    EMSNPScalar::Geometry algebraic;
    for (int i=0;i<4;++i) algebraic.g[i][i]=(i==0 ? -1 : 1);
    algebraic.tetrad=g.tetrad;
    const double frame[3][3]{{1,0,0},{0,1,0},{0,0,1}};
    const double e[3][3]{{.11,.12,.13},{.12,.22,.23},{.13,.23,-.33}};
    const double b[3][3]{{.31,.32,.33},{.32,.42,.43},{.33,.43,-.73}};
    EMSNPScalar::rebuild_weyl(algebraic,e,b,frame);
    for (int i=0;i<4;++i) for (int j=0;j<4;++j)
    {
        double trace=0;
        for (int k=0;k<4;++k) trace+=(k==0 ? -1 : 1)*algebraic.C[k][i][k][j];
        if (std::abs(trace)>1e-14) throw std::runtime_error("NP algebraic trace control failed");
        for (int k=0;k<4;++k) for (int l=0;l<4;++l)
            if (std::abs(algebraic.C[i][j][k][l]+algebraic.C[i][k][l][j]+
                         algebraic.C[i][l][j][k])>1e-14)
                throw std::runtime_error("NP algebraic Bianchi control failed");
    }
    bool rejected=false;
    v.chi=1e-5;
    try { compute.geometry(v,d1,d2,coords); } catch (const std::runtime_error &) { rejected=true; }
    if (!rejected) throw std::runtime_error("NP inconsistent metric floor was not rejected");
    pout() << "PASS NP flat, complex Maxwell signs/normalization, general E/B algebra, chi-floor rejection\n";
    mainFinalize(); return 0;
}
