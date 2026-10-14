#ifndef EMS_NP_VALIDATION_EXTRACTION_HPP
#define EMS_NP_VALIDATION_EXTRACTION_HPP
#include "EMSNPExtraction.hpp"
#include <limits>

// Test-only pointwise controls, evaluated from the production reconstructed
// tensors. FD rotates the tetrad itself; CAS coefficients are a separate path.
class EMSNPValidationExtraction : public EMSNPExtraction
{
  public:
    using EMSNPExtraction::EMSNPExtraction;
    void execute_query(AMRInterpolator<Lagrange<4>> *interp)
    {
        EMSNPExtraction::execute_query(interp);
        if (procID()!=0) return;
        std::ofstream file(m_params.data_path+"np_points.csv");
        if (!file) throw std::runtime_error("cannot write NP validation points");
        file << "r,theta,weight,area,legacy,trace";
        for (const char *name:EMSNPScalar::names) file << ',' << name << "_re," << name << "_im";
        for (const char *fixed:{"l","k"}) for (const char *x:{"out","in"})
            file << ',' << fixed << '_' << x << "_fd_re," << fixed << '_' << x << "_fd_im,"
                 << fixed << '_' << x << "_cas_re," << fixed << '_' << x << "_cas_im,"
                 << fixed << '_' << x << "_budget," << fixed << '_' << x << "_den";
        file << '\n' << std::setprecision(17);
        for (int ir=0;ir<m_params.num_surfaces;++ir) for (int j=1;j<m_params.num_points_u-1;++j)
        {
            const auto g=point(ir,j);
            const auto s=EMSNPScalar::evaluate(g,g.tetrad);
            const auto &v=s.values;
            const double theta=m_geom.u(j,m_params.num_points_u);
            const double weight=2*M_PI*m_du*std::sin(theta)*
                IntegrationMethod::simpson.weight(j,m_params.num_points_u,false);
            // Inverse metric from the supplied null tetrad: -lk-kl+m mb+mb m.
            double trace=0;
            for (int b=0;b<4;++b) for (int d=0;d<4;++d)
            {
                EMSNPScalar::Complex value=0;
                for (int a=0;a<4;++a) for (int c=0;c<4;++c)
                    value+=(-g.tetrad[0][a]*g.tetrad[1][c]-g.tetrad[1][a]*g.tetrad[0][c]
                        +g.tetrad[2][a]*g.tetrad[3][c]+g.tetrad[3][a]*g.tetrad[2][c])*g.C[a][b][c][d];
                trace=std::max(trace,std::abs(value));
            }
            file << m_params.surface_param_values[ir] << ',' << theta << ',' << weight
                 << ',' << g.area_factor << ',' << g.legacy << ',' << trace;
            for (const auto z:v) file << ',' << z.real() << ',' << z.imag();
            for (int fixed=0;fixed<2;++fixed)
            {
                // Degree <=4 in a, hence this four-point derivative cancels the
                // complete central-difference truncation error, even off shell.
                const double h=1./32;
                auto rotated=[&](double a) {
                    auto t=g.tetrad;
                    const int f=fixed==0 ? 0 : 1, other=1-f;
                    for (int i=0;i<4;++i)
                    {
                        t[other][i]+=a*(g.tetrad[2][i]+g.tetrad[3][i])+a*a*g.tetrad[f][i];
                        t[2][i]+=a*g.tetrad[f][i]; t[3][i]+=a*g.tetrad[f][i];
                    }
                    return EMSNPScalar::evaluate(g,t).values;
                };
                const auto p=rotated(h),m=rotated(-h),p2=rotated(2*h),m2=rotated(-2*h);
                const std::array<EMSNPScalar::Complex,2> cas=fixed==0 ?
                    std::array<EMSNPScalar::Complex,2>{2.*v[5]*v[3]-6.*v[7]*v[1],2.*v[6]*v[0]-4.*v[5]*v[1]} :
                    std::array<EMSNPScalar::Complex,2>{2.*v[6]*v[4]-4.*v[7]*v[3],2.*v[7]*v[1]-6.*v[5]*v[3]};
                const auto dp3=fixed==0 ? 3.*v[2] : v[4];
                const auto dp1=fixed==0 ? v[0] : 3.*v[2];
                const auto df2=fixed==0 ? 2.*v[6] : EMSNPScalar::Complex(0);
                const auto df0=fixed==0 ? EMSNPScalar::Complex(0) : 2.*v[6];
                const double den[2]{std::abs(2.*v[6]*dp3)+std::abs(3.*v[2]*df2),
                    std::abs(2.*v[6]*dp1)+std::abs(3.*v[2]*df0)};
                double cnorm=0,fnorm=0,tnorm=0;
                for (int a=0;a<4;++a) for (int b=0;b<4;++b)
                {
                    fnorm=std::max(fnorm,std::abs(g.F[a][b]));
                    for (int c=0;c<4;++c) for (int d=0;d<4;++d)
                        cnorm=std::max(cnorm,std::abs(g.C[a][b][c][d]));
                }
                for (const auto &t:g.tetrad)
                {
                    double sum=0; for (const auto z:t) sum+=std::abs(z);
                    tnorm=std::max(tnorm,sum);
                }
                const double budget=4096*std::numeric_limits<double>::epsilon()*cnorm*fnorm*
                    std::pow(tnorm*(1+2*h)*(1+2*h),6)/h;
                for (int x=0;x<2;++x)
                {
                    const auto fd=(8.*(p[8+x]-m[8+x])-(p2[8+x]-m2[8+x]))/(12*h);
                    file << ',' << fd.real() << ',' << fd.imag() << ',' << cas[x].real()
                         << ',' << cas[x].imag() << ',' << budget << ',' << den[x];
                }
            }
            file << '\n';
        }
        file.close();
        if (!file) throw std::runtime_error("failed writing NP validation points");
    }
};
#endif
