#ifndef EMS_NP_CHECKPOINT_JET_HPP
#define EMS_NP_CHECKPOINT_JET_HPP

#include "EMSMatterWeylOmScalar.hpp"
#include <fstream>
#include <iomanip>

// Test-only, one valid cell: export the reconstruction before any NP contraction.
class EMSNPCheckpointJet : public EMSMatterWeylOmScalar
{
    double m_x, m_y;
    bool m_matter;
    std::string m_path;
  public:
    EMSNPCheckpointJet(const std::array<double, CH_SPACEDIM> &center,
        double dx, CouplingFunction::params_t coupling, double G,
        double x, double y, bool matter, const std::string &path)
        : EMSMatterWeylOmScalar(center, dx, coupling, G), m_x(x), m_y(y),
          m_matter(matter), m_path(path) {}

    void compute(Cell<double> cell) const
    {
        const Coordinates<double> coords(cell, m_deriv.m_dx, m_center);
        if (std::abs(coords.x-m_x) >= m_deriv.m_dx/2 ||
            std::abs(coords.y-m_y) >= m_deriv.m_dx/2) return;
        const auto vars = cell.load_vars<Vars>();
        const auto d1 = m_deriv.diff1<Vars>(cell);
        const auto d2 = m_deriv.diff2<Diff2Vars>(cell);
        if (vars.chi <= 1e-4)
            throw std::runtime_error("NP probe must be outside both Weyl chi floors");
        auto r = compute_riemann_tensor(vars, d1, d2, coords);
        if (m_matter) add_matter(r, vars, d1);
        const auto t = compute_null_tetrad(vars, coords);
        const auto legacy = compute_Weyl_Om(r, vars, d1, d2, coords);
        const auto hu = TensorAlgebra::compute_inverse_sym(vars.h);
        const auto chris = TensorAlgebra::compute_christoffel(d1.h, hu);
        const auto em = compute_EMS_EM_tensor(vars, d1, hu, 1./vars.hww,
            chris, GR_SPACEDIM-CH_SPACEDIM, m_coupling);
        double R[4][4][4][4]{}, ric[4][4]{}, g[4][4]{}, F[4][4]{};
        auto set = [&](int a, int b, int c, int d, double value) {
            R[a][b][c][d]=R[c][d][a][b]=R[b][a][d][c]=R[d][c][b][a]=value;
            R[b][a][c][d]=R[a][b][d][c]=R[d][c][a][b]=R[c][d][b][a]=-value;
        };
        FOR(i,j)
        {
            set(0,i+1,0,j+1,r.riemann_0d0d[i][j]);
            set(i+1,3,j+1,3,r.riemann_dwdw[i][j]);
            FOR(k)
            {
                set(0,i+1,j+1,k+1,r.riemann_0ddd[i][j][k]);
                FOR(l) { set(i+1,j+1,k+1,l+1,r.riemann_dddd[i][j][k][l]); }
            }
        }
        set(0,3,0,3,r.riemann_0w0w);
        FOR(i) { set(0,3,i+1,3,r.riemann_0wdw[i]); }
        g[0][0]=-1.;
        FOR(i,j) { g[i+1][j+1]=vars.h[i][j]/vars.chi; }
        g[3][3]=vars.hww/vars.chi;
        const double eightpiG=8*M_PI*m_G_Newton;
        ric[0][0]=eightpiG*(em.rho+em.S)/2;
        FOR(i)
        {
            // Pi=-n(phi), so the evolution's positive Si is -T_ni.
            ric[0][i+1]=ric[i+1][0]=-eightpiG*em.Si[i];
            FOR(j) { ric[i+1][j+1]=eightpiG*(em.Sij[i][j]-(em.S-em.rho)*g[i+1][j+1]/2); }
        }
        ric[3][3]=eightpiG*(em.Sww-(em.S-em.rho)*g[3][3]/2);
        const double E[3]{vars.Ex,vars.Ey,vars.Ez};
        const double B[3]{vars.Bx,vars.By,vars.Bz};
        double bu[3]{};
        FOR(i,j) { bu[i]+=vars.chi*hu[i][j]*B[j]; }
        bu[2]=vars.chi*B[2]/vars.hww;
        const double rootdet=std::sqrt((g[1][1]*g[2][2]-g[1][2]*g[2][1])*g[3][3]);
        for (int i=0;i<3;++i) { F[i+1][0]=E[i]; F[0][i+1]=-E[i]; }
        F[1][2]=rootdet*bu[2]; F[2][3]=rootdet*bu[0]; F[1][3]=-rootdet*bu[1];
        for (int i=1;i<4;++i) for (int j=i+1;j<4;++j) F[j][i]=-F[i][j];
        // The independent contraction reader receives n,u,v,w in that order.
        double frame[4][4]{{1,0,0,0},{0,t.m1[0],t.m1[1],0},
                            {0,t.m2[0],t.m2[1],0},{0,0,0,t.mw}};
        std::ofstream file(m_path);
        if (!file) throw std::runtime_error("cannot write NP checkpoint jet");
        file << std::setprecision(17) << "point " << coords.x << ' ' << coords.y
             << ' ' << m_deriv.m_dx << '\n' << "legacy " << legacy.Om22 << '\n';
        for (int a=0;a<4;++a) for (int b=0;b<4;++b)
        {
            file << "metric " << a << ' ' << b << ' ' << g[a][b] << '\n';
            file << "ricci " << a << ' ' << b << ' ' << ric[a][b] << '\n';
            file << "faraday " << a << ' ' << b << ' ' << F[a][b] << '\n';
            file << "frame " << a << ' ' << b << ' ' << frame[a][b] << '\n';
            for (int c=0;c<4;++c) for (int d=0;d<4;++d)
                file << "riemann " << a << ' ' << b << ' ' << c << ' ' << d << ' ' << R[a][b][c][d] << '\n';
        }
        if (!file) throw std::runtime_error("failed writing NP checkpoint jet");
    }
};
#endif
