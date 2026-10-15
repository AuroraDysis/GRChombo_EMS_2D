#ifndef EMS_NP_VALIDATION_EXTRACTION_HPP
#define EMS_NP_VALIDATION_EXTRACTION_HPP
#include "EMSNPExtraction.hpp"
#include <limits>
#include "ConstraintsCartoon.hpp"

class EMSNPPointConstraints : public Constraints<CouplingFunction>
{
  public:
    using Constraints<CouplingFunction>::Constraints;
    using Constraints<CouplingFunction>::constraint_equations;
};

// Test-only pointwise controls, evaluated from the production reconstructed
// tensors. FD rotates the tetrad itself; CAS coefficients are a separate path.
class EMSNPValidationExtraction : public EMSNPExtraction
{
    EMSNPPointConstraints m_constraints;
    auto geometric_constraints(Jets jets) const
    {
        // ConstraintsCartoon normally consumes the independently evolved Gamma.
        // Rebuild Gamma and its derivative to also report the metric constraints.
        auto &v=jets.v; auto &d=jets.d1;
        const auto gu=TensorAlgebra::compute_inverse_sym(v.h);
        double dg[2][2][2]{};
        FOR(i,j,p,k,l) { dg[i][j][p]-=gu[i][k]*d.h[k][l][p]*gu[l][j]; }
        FOR(i)
        {
            v.Gamma[i]=(i==1 ? 1/(v.hww*jets.coords.y) : 0)-gu[i][1]/jets.coords.y;
            FOR(j,k,l) { v.Gamma[i]+=gu[i][l]*gu[j][k]*(d.h[l][k][j]-.5*d.h[j][k][l]); }
            FOR(j) { v.Gamma[i]-=.5*gu[i][j]*d.hww[j]/v.hww; }
            FOR(p)
            {
                d.Gamma[i][p]=-(i==1 ? d.hww[p]/(v.hww*v.hww*jets.coords.y)+
                    (p==1 ? 1/(v.hww*jets.coords.y*jets.coords.y) : 0) : 0)
                    -dg[i][1][p]/jets.coords.y+(p==1 ? gu[i][1]/(jets.coords.y*jets.coords.y) : 0);
                FOR(j,k,l)
                {
                    d.Gamma[i][p]+=(dg[i][l][p]*gu[j][k]+gu[i][l]*dg[j][k][p])*
                        (d.h[l][k][j]-.5*d.h[j][k][l])+gu[i][l]*gu[j][k]*
                        (jets.d2.h[l][k][j][p]-.5*jets.d2.h[j][k][l][p]);
                }
                FOR(j)
                {
                    d.Gamma[i][p]-=.5*((dg[i][j][p]*d.hww[j]+gu[i][j]*jets.d2.hww[j][p])/v.hww
                        -gu[i][j]*d.hww[j]*d.hww[p]/(v.hww*v.hww));
                }
            }
        }
        return m_constraints.constraint_equations(v,d,jets.d2,jets.coords.y);
    }
  public:
    EMSNPValidationExtraction(const spherical_extraction_params_t &p,double dt,double time,
        bool first,double restart,CouplingFunction::params_t coupling,double G,bool matter=true)
        : EMSNPExtraction(p,dt,time,first,restart,coupling,G,matter),
          m_constraints(1.,CouplingFunction(coupling),G)
    {
        for (int var:{c_Gamma1,c_Gamma2})
            for (const auto &d:{Derivative::dx,Derivative::dy}) request(var,d);
    }
    void execute_query(AMRInterpolator<Lagrange<4>> *interp)
    {
        EMSNPExtraction::execute_query(interp);
        if (procID()!=0) return;
        std::ofstream file(m_params.data_path+"np_points.csv");
        if (!file) throw std::runtime_error("cannot write NP validation points");
        file << "r,theta,weight,area,legacy,trace,E_trace,B_trace,B_axial0,B_axial1,B_axial2,"
                "Ham,Ham_metric,K_traceA,Mom0,Mom1,Mom2,psi0_change,psi4_change,projection_budget";
        for (const char *name:EMSNPScalar::names) file << ',' << name << "_re," << name << "_im";
        for (const char *fixed:{"l","k"}) for (const char *x:{"out","in"})
            file << ',' << fixed << '_' << x << "_fd_re," << fixed << '_' << x << "_fd_im,"
                 << fixed << '_' << x << "_cas_re," << fixed << '_' << x << "_cas_im,"
                 << fixed << '_' << x << "_budget," << fixed << '_' << x << "_den";
        file << '\n' << std::setprecision(17);
        for (int ir=0;ir<m_params.num_surfaces;++ir) for (int j=1;j<m_params.num_points_u-1;++j)
        {
            const auto g=point(ir,j);
            const auto jets=point_jets(ir,j);
            const auto constraints=m_constraints.constraint_equations(
                jets.v,jets.d1,jets.d2,jets.coords.y);
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
                 << ',' << g.area_factor << ',' << g.legacy << ',' << trace
                 << ',' << g.electric_trace << ',' << g.magnetic_trace;
            for (double x:g.magnetic_axial) file << ',' << x;
            const auto metric_constraints=geometric_constraints(jets);
            const auto inverse=TensorAlgebra::compute_inverse_sym(jets.v.h);
            double traceA=jets.v.Aww/jets.v.hww;
            FOR(i,k) { traceA+=inverse[i][k]*jets.v.A[i][k]; }
            file << ',' << constraints.Ham << ',' << metric_constraints.Ham
                 << ',' << jets.v.K*traceA/3;
            double frame[3][3]{};
            for (int i=0;i<3;++i)
            {
                frame[0][i]=((g.tetrad[0][i+1]-g.tetrad[1][i+1])/std::sqrt(2.)).real();
                frame[1][i]=((g.tetrad[2][i+1]+g.tetrad[3][i+1])/std::sqrt(2.)).real();
                frame[2][i]=((g.tetrad[2][i+1]-g.tetrad[3][i+1])/EMSNPScalar::Complex(0,std::sqrt(2.))).real();
            }
            for (int a=0;a<3;++a)
            {
                double mom=0; FOR(i) { mom+=frame[a][i]*constraints.Mom[i]; }
                file << ',' << mom;
            }
            // Independent historical field-Ricci subtraction, used only to
            // certify that the STF reconstruction preserves Psi0 and Psi4.
            auto old=g;
            double scalar=0;
            const auto hu=TensorAlgebra::compute_inverse_sym(jets.v.h);
            scalar=-g.ric[0][0];
            FOR(i,k) { scalar+=jets.v.chi*hu[i][k]*g.ric[i+1][k+1]; }
            scalar+=jets.v.chi*g.ric[3][3]/jets.v.hww;
            double norm=0,tnorm=0;
            for (int a=0;a<4;++a) for (int b=0;b<4;++b)
            for (int c=0;c<4;++c) for (int d=0;d<4;++d)
            {
                old.C[a][b][c][d]=g.R[a][b][c][d]
                    -(g.g[a][c]*g.ric[d][b]-g.g[a][d]*g.ric[c][b]
                     -g.g[b][c]*g.ric[d][a]+g.g[b][d]*g.ric[c][a])/2
                    +scalar*(g.g[a][c]*g.g[d][b]-g.g[a][d]*g.g[c][b])/6;
                norm=std::max(norm,std::abs(g.R[a][b][c][d])+std::abs(old.C[a][b][c][d]));
            }
            for (const auto &t:g.tetrad)
            {
                double sum=0; for (auto z:t) sum+=std::abs(z); tnorm=std::max(tnorm,sum);
            }
            const auto previous=EMSNPScalar::evaluate(old,g.tetrad);
            file << ',' << std::abs(previous.values[0]-v[0]) << ',' << std::abs(previous.values[4]-v[4])
                 << ',' << 4096*std::numeric_limits<double>::epsilon()*norm*std::pow(tnorm,4);
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
