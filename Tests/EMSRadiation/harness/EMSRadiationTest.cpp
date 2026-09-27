#include "EMSRadiationExtraction.hpp"
#include "SetupFunctions.hpp"
#include "WeylOmScalar.hpp"
#include <iostream>

struct EMSWeylProbe : WeylOmScalar
{
    using WeylOmScalar::WeylOmScalar;
    using WeylOmScalar::compute_Weyl_Om;
};

void weyl_normalization_check()
{
    // Local outgoing TT curvature, propagation along +y, h_plus''=0.13.
    // Supply curvature directly to isolate the author's null-tetrad factor.
    const double hdd=.13;
    Riemann_t<double> r{};
    r.riemann_0d0d[0][0]=-hdd/2;
    r.riemann_0ddd[0][0][1]=-hdd/2;
    r.riemann_dddd[1][0][1][0]=-hdd/2;
    r.riemann_0w0w=hdd/2;
    r.riemann_0wdw[1]=-hdd/2;
    r.riemann_dwdw[1][1]=hdd/2;
    CCZ4CartoonVars::VarsWithGauge<double> v{};
    v.chi=v.lapse=v.hww=v.h[0][0]=v.h[1][1]=1;
    CCZ4CartoonVars::VarsWithGauge<Tensor<1,double>> d1{};
    CCZ4CartoonVars::Diff2VarsNoGauge<Tensor<2,double>> d2{};
    Coordinates<double> x(IntVect::Zero,1); x.x=0; x.y=1;
    const auto out=EMSWeylProbe({0,0},1).compute_Weyl_Om(r,v,d1,d2,x);
    if (std::abs(out.Om22+hdd)>1e-12)
        throw std::runtime_error("WeylOmScalar tetrad normalization check failed");
    std::cout << "EMS GW tetrad check: Psi4=" << out.Om22 << " expected=-0.13\n";
}

void stress_check()
{
    using namespace EMSNativeRadiation;
    double worst=0;
    for (int sample=0;sample<24;++sample)
    {
        Fields f;
        f.lapse=.6+.02*sample; f.Pi=.2;
        f.g={{{1.2,.1,0},{.1,1.5,0},{0,0,.9}}};
        f.shift={.03,-.07,0}; f.grad={.3,-.2,0};
        f.E={.1,std::sin(sample+.3),.2}; f.B={.2,.3,std::cos(sample+.1)};
        const double det2=1.2*1.5-.01, rootdet=std::sqrt(det2*.9), coupling=1.7;
        const Mat inv{{{1.5/det2,-.1/det2,0},{-.1/det2,1.2/det2,0},{0,0,1/.9}}};
        double g[4][4]{}, gu[4][4]{}, field[4][4]{}, T[4][4]{}, Ts[4][4]{};
        g[0][0]=-f.lapse*f.lapse+dot(f.shift,mul(f.g,f.shift));
        gu[0][0]=-1/(f.lapse*f.lapse);
        for (int i=0;i<3;++i)
        {
            g[0][i+1]=g[i+1][0]=mul(f.g,f.shift)[i];
            gu[0][i+1]=gu[i+1][0]=f.shift[i]/(f.lapse*f.lapse);
            for (int j=0;j<3;++j)
            {
                g[i+1][j+1]=f.g[i][j];
                gu[i+1][j+1]=inv[i][j]-f.shift[i]*f.shift[j]/(f.lapse*f.lapse);
                Vec ei{},ej{}; ei[i]=1; ej[j]=1;
                field[i+1][j+1]=rootdet*dot(cross(ei,ej),mul(inv,f.B));
            }
        }
        for (int i=0;i<3;++i)
        {
            field[i+1][0]=f.lapse*f.E[i];
            for (int j=0;j<3;++j) field[i+1][0]+=field[i+1][j+1]*f.shift[j];
            field[0][i+1]=-field[i+1][0];
        }
        double grad[4]{-f.lapse*f.Pi+dot(f.shift,f.grad),f.grad[0],f.grad[1],0};
        double norm=0,ff=0;
        for (int a=0;a<4;++a) for (int b=0;b<4;++b)
        {
            norm+=gu[a][b]*grad[a]*grad[b];
            for (int c=0;c<4;++c) for (int d=0;d<4;++d)
                ff+=field[a][b]*gu[a][c]*gu[b][d]*field[c][d];
        }
        for (int a=0;a<4;++a) for (int b=0;b<4;++b)
        {
            Ts[a][b]=2*grad[a]*grad[b]-g[a][b]*norm;
            T[a][b]=-coupling*g[a][b]*ff/2;
            for (int c=0;c<4;++c) for (int d=0;d<4;++d)
                T[a][b]+=2*coupling*field[a][c]*gu[c][d]*field[b][d];
        }
        const double theta=.1+sample*.12;
        const Vec dr{std::cos(theta),std::sin(theta),0};
        double em=0,scalar=0;
        for (int i=0;i<3;++i) for (int a=0;a<4;++a)
        {
            em-=f.lapse*rootdet*dr[i]*gu[i+1][a]*T[a][0];
            scalar-=f.lapse*rootdet*dr[i]*gu[i+1][a]*Ts[a][0];
        }
        const auto v=evaluate(f,theta,coupling);
        worst=std::max(worst,std::max(std::abs(em-v.em_flux),std::abs(scalar-v.scalar_flux)));
    }
    if (worst>1e-12) throw std::runtime_error("four-tensor stress check failed");
    std::cout << "EMS stress four-tensor check: 24 curved metric/shift samples, max error " << worst << '\n';
}

// Exact surface fields use the same production evaluation, quadrature,
// harmonics and CSV writer. The separate hierarchy smoke tests interpolation.
class EMSAnalyticExtraction : public EMSRadiationExtraction
{
  public:
    using EMSRadiationExtraction::EMSRadiationExtraction;
    void fill(int mode,double time,int kind,double radius)
    {
        for (int j=0;j<m_params.num_points_u;++j)
        {
            const double theta=m_geom.u(j,m_params.num_points_u);
            const double y=SphericalHarmonics::spin_Y_lm(std::sin(theta),0.,
                                      std::cos(theta),0,mode,0).Real;
            auto put=[&](int k,double v){m_interp_data[k][j]=v;};
            for (int k=0;k<20;++k) put(k,0);
            for (int k:{2,5,6,8,9}) put(k,1);
            if (kind==0) // exact spherical solution phi=f(t-r)/r * Y00
            {
                const double u=time-4, amp=0.03, sigma=0.5;
                const double f=amp*std::exp(-u*u/(2*sigma*sigma));
                const double df=-u*f/(sigma*sigma);
                put(0,f*y/radius); put(1,-df*y/radius);
                const double grad=(-df/radius-f/(radius*radius))*y;
                put(16,grad*std::cos(theta)); put(17,grad*std::sin(theta));
            }
            else if (kind==1) { put(0,0.03*y/radius); put(1,-0.07*y/radius); }
            else // transverse outgoing / incoming dipole, F_inf=exp(-0.4)
            {
                const double et=0.02*std::sin(theta)/radius;
                put(10,-et*std::sin(theta)); put(11,et*std::cos(theta));
                put(15,(kind==2 ? 1 : -1)*et);
            }
        }
    }
};

int main(int argc,char **argv)
{
    mainSetup(argc,argv);
    int status=0;
    try
    {
        stress_check();
        weyl_normalization_check();
        if (argc!=5) throw std::runtime_error("usage: test output theta kind mode");
        spherical_extraction_params_t p;
        p.num_extraction_radii=1; p.extraction_radii={20}; p.extraction_levels={0};
        p.center={0,0}; p.num_points_theta=std::stoi(argv[2]); p.num_points_phi=2;
        p.num_modes=0; p.write_extraction=false; p.data_path=std::string(argv[1])+"/";
        const int kind=std::stoi(argv[3]), mode=std::stoi(argv[4]);
        const int steps=kind==0 ? 800 : 0;
        for (int i=0;i<=steps;++i)
        {
            const double t=i*0.01;
            CouplingFunction::params_t c{1,kind>=2 ? 0.2 : 0.,0,0};
            EMSAnalyticExtraction e(p,0.01,t,i==0,0,c,0);
            e.fill(mode,t,kind,20);
            e.write_modes();
        }
        std::cout << "EMS analytic extraction complete\n";
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; status=1; }
    mainFinalize();
    return status;
}
