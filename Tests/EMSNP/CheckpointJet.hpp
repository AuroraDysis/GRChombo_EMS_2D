#ifndef EMS_NP_CHECKPOINT_JET_HPP
#define EMS_NP_CHECKPOINT_JET_HPP

#include "EMSNPScalar.hpp"
#include <fstream>
#include <iomanip>

// Test-only, one valid cell: export the reconstruction before any NP contraction.
class EMSNPCheckpointJet : public EMSNPScalar
{
    double m_x, m_y;
    bool m_matter;
    std::string m_path;
  public:
    EMSNPCheckpointJet(const std::array<double, CH_SPACEDIM> &center,
        double dx, CouplingFunction::params_t coupling, double G,
        double x, double y, bool matter, const std::string &path)
        : EMSNPScalar(center, dx, coupling, G), m_x(x), m_y(y),
          m_matter(matter), m_path(path) {}

    void compute(Cell<double> cell) const
    {
        const Coordinates<double> coords(cell, m_deriv.m_dx, m_center);
        if (std::abs(coords.x-m_x) >= m_deriv.m_dx/2 ||
            std::abs(coords.y-m_y) >= m_deriv.m_dx/2) return;
        const auto vars = cell.load_vars<Vars>();
        const auto d1 = m_deriv.diff1<Vars>(cell);
        const auto d2 = m_deriv.diff2<Diff2Vars>(cell);
        const auto out=geometry(vars,d1,d2,coords,m_matter);
        const auto t=compute_null_tetrad(vars,coords);
        const auto &R=out.R;
        const auto &ric=out.ric;
        const auto &g=out.g;
        const auto &F=out.F;
        // The independent contraction reader receives n,u,v,w in that order.
        double frame[4][4]{{1,0,0,0},{0,t.m1[0],t.m1[1],0},
                            {0,t.m2[0],t.m2[1],0},{0,0,0,t.mw}};
        std::ofstream file(m_path);
        if (!file) throw std::runtime_error("cannot write NP checkpoint jet");
        file << std::setprecision(17) << "point " << coords.x << ' ' << coords.y
             << ' ' << m_deriv.m_dx << '\n' << "legacy " << out.legacy << '\n';
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
