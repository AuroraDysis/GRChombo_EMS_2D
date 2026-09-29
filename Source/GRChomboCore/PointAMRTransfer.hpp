/* Sixth-order transfers for cell-centred POINT values, ratio two, 2D.
 * No initial-data object or reference state is retained by this operator. */
#ifndef POINTAMRTRANSFER_HPP_
#define POINTAMRTRANSFER_HPP_

#include "BoundaryConditions.hpp"
#include "IntVectSet.H"
#include "LayoutData.H"
#include "TimeInterpolatorRK4.H"
#include <array>
#include <cmath>

// Keep Chombo's storage/dense polynomial, correcting its two diff12
// coefficients locally. Eq. 49 has a COARSE dt^2 denominator: the stage
// correction is dt_fine^3/dt_coarse^3, not the r^-2 used by this Chombo.
class PointRK4Interpolator : public TimeInterpolatorRK4
{
  public:
    void intermediate(LevelData<FArrayBox> &u, double theta, int stage,
                      const Interval &comps) const
    {
        TimeInterpolatorRK4::intermediate(u, theta, stage, comps);
        if (stage == 2 || stage == 3)
        {
            const double r = 1. / m_refineCoarse;
            const double correction = r * r * (r - 1.) * (stage == 2 ? -.25 : .5);
            for (DataIterator dit = u.dataIterator(); dit.ok(); ++dit)
                u[dit].plus(m_diff12[dit], correction, comps.begin(),
                            comps.begin(), comps.size());
        }
    }
};

class PointAMRTransfer
{
  public:
    struct Stencil
    {
        IntVect fine, first;
        std::array<std::array<double, 6>, 2> weights;
    };

    // Lagrange at a point on six consecutive cell centres: degree five exact.
    static std::array<double, 6> weights(double q, int first)
    {
        std::array<double, 6> w;
        for (int k = 0; k < 6; ++k)
        {
            w[k] = 1.;
            for (int j = 0; j < 6; ++j)
                if (j != k)
                    w[k] *= (q - first - j) / (k - j);
        }
        return w;
    }

    void define(const DisjointBoxLayout &fine,
                const DisjointBoxLayout &coarse, const ProblemDomain &domain,
                int ratio, int ghosts, double coarse_dx,
                const std::array<double, CH_SPACEDIM> &center,
                const BoundaryConditions::params_t &bdy)
    {
        if (CH_SPACEDIM != 2 || ratio != 2 || ghosts != 3)
            MayDay::Error("amr_transfer=point requires DIM=2, ratio=2, three ghosts");
        m_coarse_domain = coarsen(domain, 2);
        m_bdy = bdy;
        m_coarsened = DisjointBoxLayout();
        coarsen(m_coarsened, fine, 2);
        m_coarse_data.define(m_coarsened, NUM_VARS, 4 * IntVect::Unit);
        m_restricted.define(m_coarsened, NUM_VARS);
        m_copier.define(coarse, m_coarsened, m_coarse_domain,
                        4 * IntVect::Unit);
        m_coarse_boundaries.define(coarse_dx, center, bdy, m_coarse_domain, 4);
        time.define(fine, coarse, domain, 2, NUM_VARS, 4);
        m_ghost_stencils.define(fine);
        for (DataIterator dit = fine.dataIterator(); dit.ok(); ++dit)
        {
            Box ghost_domain = domain.domainBox();
            for (int dir = 0; dir < 2; ++dir)
                if (bdy.is_periodic[dir])
                    ghost_domain.grow(dir, 3);
            IntVectSet ghosts_to_fill(grow(fine[dit], 3) & ghost_domain);
            for (LayoutIterator lit = fine.layoutIterator(); lit.ok(); ++lit)
                for (int j = bdy.is_periodic[1] ? -1 : 0;
                     j <= (bdy.is_periodic[1] ? 1 : 0); ++j)
                    for (int i = bdy.is_periodic[0] ? -1 : 0;
                         i <= (bdy.is_periodic[0] ? 1 : 0); ++i)
                    {
                        Box image = fine[lit];
                        image.shift(IntVect(D_DECL(i * domain.domainBox().size(0),
                                                   j * domain.domainBox().size(1), 0)));
                        ghosts_to_fill -= image;
                    }
            auto &stencils = m_ghost_stencils[dit];
            for (IVSIterator iv(ghosts_to_fill); iv.ok(); ++iv)
                stencils.push_back(stencil(iv()));
        }
    }

    // Used by both normal ghost fills and regridding; exchange wins at seams.
    void fill(LevelData<FArrayBox> &fine, const LevelData<FArrayBox> &coarse,
              const Interval &comps, bool valid_too = false,
              VariableType var_type = VariableType::evolution)
    {
        coarse.copyTo(comps, m_coarse_data, comps, m_copier);
        fill_coarse_boundaries(comps, var_type);
        interpolate(fine, comps, valid_too);
    }

    // theta is the START of this fine step in the enclosing coarse step.
    // It is deliberately the same for all four stages (including the two
    // distinct half-time stages); Chombo supplies McCorquodale--Colella 43--45.
    void fill_stage(LevelData<FArrayBox> &fine, double theta, int stage)
    {
        const Interval comps(0, NUM_VARS - 1);
        time.intermediate(m_coarse_data, theta, stage, comps);
        fill_coarse_boundaries(comps);
        interpolate(fine, comps, false);
    }

    // Caller fills fine's same-level, coarse-fine and parity ghosts FIRST.
    // Six fine points straddle the coarse centre; degree five, O(h^6).
    void restrict_to_coarse(LevelData<FArrayBox> &coarse,
                            const LevelData<FArrayBox> &fine)
    {
        const auto w = weights(.5, -2);
        const DataIterator dit = fine.dataIterator();
#pragma omp parallel for schedule(static)
        for (int ibox = 0; ibox < dit.size(); ++ibox)
        {
            const DataIndex di = dit[ibox];
            for (BoxIterator bit(m_coarsened[di]); bit.ok(); ++bit)
                for (int comp = 0; comp < NUM_VARS; ++comp)
                {
                    const double base = fine[di](2 * bit(), comp);
                    double correction = 0.;
                    for (int j = 0; j < 6; ++j)
                        for (int i = 0; i < 6; ++i)
                            correction += w[i] * w[j] * (fine[di](
                                2 * bit() + IntVect(D_DECL(i - 2, j - 2, 0)), comp) - base);
                    m_restricted[di](bit(), comp) = base + correction;
                }
        }
        m_restricted.copyTo(coarse);
    }

    PointRK4Interpolator time;

  private:
    Stencil stencil(const IntVect &iv) const
    {
        Stencil s;
        s.fine = iv;
        for (int dir = 0; dir < 2; ++dir)
        {
            const double q = (iv[dir] + .5) / 2. - .5;
            int lo = m_coarse_domain.domainBox().smallEnd(dir);
            int hi = m_coarse_domain.domainBox().bigEnd(dir);
            // Reflective ghosts contain parity, so do not shift away from
            // the cartoon axis. Periodic copies also admit centred stencils.
            if (m_bdy.is_periodic[dir] ||
                m_bdy.lo_boundary[dir] == BoundaryConditions::REFLECTIVE_BC)
                lo -= 4;
            if (m_bdy.is_periodic[dir] ||
                m_bdy.hi_boundary[dir] == BoundaryConditions::REFLECTIVE_BC)
                hi += 4;
            if (hi - lo < 5)
                MayDay::Error("point transfer domain has fewer than six cells");
            s.first[dir] = std::max(lo, std::min(int(std::floor(q)) - 2, hi - 5));
            s.weights[dir] = weights(q, s.first[dir]);
        }
        return s;
    }

    void fill_coarse_boundaries(const Interval &comps,
                               VariableType var_type = VariableType::evolution)
    {
        if (var_type == VariableType::diagnostic)
        {
            m_coarse_boundaries.fill_diagnostic_boundaries(Side::Hi, m_coarse_data, comps);
            m_coarse_boundaries.fill_diagnostic_boundaries(Side::Lo, m_coarse_data, comps);
        }
        else
        {
            m_coarse_boundaries.fill_solution_boundaries(Side::Hi, m_coarse_data, comps);
            m_coarse_boundaries.fill_solution_boundaries(Side::Lo, m_coarse_data, comps);
        }
    }

    static void apply(FArrayBox &dst, const FArrayBox &src, const Stencil &s,
                      const Interval &comps)
    {
        CH_assert(src.box().contains(s.first));
        CH_assert(src.box().contains(s.first + IntVect(D_DECL(5, 5, 0))));
        for (int comp = comps.begin(); comp <= comps.end(); ++comp)
        {
            // Preserve constants exactly and reduce roundoff in second
            // derivatives. The anchor is a CURRENT source cell, not U_*.
            const double base = src(s.first + IntVect(D_DECL(2, 2, 0)), comp);
            double correction = 0.;
            for (int j = 0; j < 6; ++j)
                for (int i = 0; i < 6; ++i)
                    correction += s.weights[0][i] * s.weights[1][j] *
                                  (src(s.first + IntVect(D_DECL(i, j, 0)), comp) - base);
            dst(s.fine, comp) = base + correction;
        }
    }

    void interpolate(LevelData<FArrayBox> &fine, const Interval &comps,
                     bool valid_too)
    {
        const DataIterator dit = fine.dataIterator();
#pragma omp parallel for schedule(static)
        for (int ibox = 0; ibox < dit.size(); ++ibox)
        {
            const DataIndex di = dit[ibox];
            if (valid_too)
                for (BoxIterator bit(fine[di].box()); bit.ok(); ++bit)
                    apply(fine[di], m_coarse_data[di], stencil(bit()), comps);
            else
                for (const auto &s : m_ghost_stencils[di])
                    apply(fine[di], m_coarse_data[di], s, comps);
        }
    }

    DisjointBoxLayout m_coarsened;
    ProblemDomain m_coarse_domain;
    BoundaryConditions::params_t m_bdy;
    BoundaryConditions m_coarse_boundaries;
    GRLevelData m_coarse_data, m_restricted;
    Copier m_copier;
    LayoutData<std::vector<Stencil>> m_ghost_stencils;
};
#endif
