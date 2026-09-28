#ifndef EMSKS2_REFERENCE_CACHE_HPP_
#define EMSKS2_REFERENCE_CACHE_HPP_

#include "1D_SOL/EMSKS2Profile.hpp"
#include "ReferenceStationaryGauge.hpp"
#include "BoxIterator.H"
#include "FArrayBox.H"
#include <array>
#include <cmath>

class EMSKS2ReferenceCache
{
  public:
    static void fill(FArrayBox &out, const EMSKS2Profile &profile, double dx,
                     const std::array<double, CH_SPACEDIM> &centre)
    {
        for (BoxIterator bit(out.box()); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            const double x = (iv[0] + .5) * dx - centre[0];
            const double y = (iv[1] + .5) * dx - centre[1];
            if (x == 0 && y == 0)
                MayDay::Error("EMSKS 2 reference sampled at puncture centre");
            const auto o = profile.object(x, y);
            out(iv, ReferenceStationaryGauge::alpha_star) = o.s.alpha;
            out(iv, ReferenceStationaryGauge::q_star) = o.s.q_star;
            out(iv, ReferenceStationaryGauge::K_star) = o.s.K;
            out(iv, ReferenceStationaryGauge::beta_x) = o.beta[0];
            out(iv, ReferenceStationaryGauge::beta_y) = o.beta[1];
            out(iv, ReferenceStationaryGauge::taper) = o.s.taper;
        }
    }
};

#endif
