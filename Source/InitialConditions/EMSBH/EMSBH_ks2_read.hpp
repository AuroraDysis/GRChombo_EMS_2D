#ifndef EMSBH_KS2_READ_HPP_
#define EMSBH_KS2_READ_HPP_

#include "1D_SOL/EMSKS2Profile.hpp"
#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "EMSBHParams.hpp"
#include "EMSCouplingFunction.hpp"
#include "MayDay.H"
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include <atomic>
#include <cmath>
#include <memory>
#include <string>

class EMSBH_ks2_read
{
  public:
    EMSBH_ks2_read(const EMSBH_params_t &params,
                   const CouplingFunction::params_t &coupling, double dx,
                   double min_chi = 1e-8, double min_lapse = 1e-8)
        : params_(params), dx_(dx), min_chi_(min_chi), min_lapse_(min_lapse)
    {
        if (params.binary || params.boosted || params.rapidity != 0)
            MayDay::Error("emsks2 supports one unboosted object only");
        if (params.star_centre[1] != 0 ||
            std::abs(params.star_centre[0] / dx -
                     std::round(params.star_centre[0] / dx)) > 1e-10)
            MayDay::Error("emsks2 puncture must lie on a grid vertex on every level");
        std::string reason;
        if (!profile_.check_file(params.data_path, reason))
            MayDay::Error(("EMSKS 2: " + reason).c_str());
        if (profile_.get("ems_alpha") != coupling.alpha ||
            profile_.get("ems_f0") != coupling.f0 ||
            profile_.get("ems_f1") != coupling.f1 ||
            profile_.get("ems_f2") != coupling.f2)
            MayDay::Error("EMSKS 2 coupling differs from run parameters");
        if (std::abs(profile_.get("M") - params.bh_mass) >
            1e-12 * profile_.get("M"))
            MayDay::Error("EMSKS 2 mass differs from run parameter");
    }

    const EMSKS2Profile &profile() const { return profile_; }
    std::size_t chi_floor_count() const { return counts_->chi.load(); }
    std::size_t lapse_floor_count() const { return counts_->lapse.load(); }

    CCZ4CartoonVars::VarsWithGauge<double>
    evaluate(double x, double y, double z = 0) const
    {
        const auto o = profile_.object(x, y, z);
        const auto &s = o.s;
        CCZ4CartoonVars::VarsWithGauge<double> vars;
        VarsTools::assign(vars, 0.);
        vars.chi = s.chi;
        vars.K = s.K;
        FOR(i, j)
        { vars.h[i][j] = o.h[3 * i + j];
          vars.A[i][j] = o.A[3 * i + j]; }
        vars.hww = o.hww;
        vars.Aww = o.Aww;
        vars.Gamma[0] = o.Gamma[0];
        vars.Gamma[1] = o.Gamma[1];
        vars.lapse = s.alpha;
        vars.shift[0] = o.beta[0];
        vars.shift[1] = o.beta[1];
        vars.phi = s.phi;
        vars.Pi = s.Pi;
        vars.Ex = o.E[0];
        vars.Ey = o.E[1];
        vars.Ez = o.E[2];
        return vars;
    }

    template <class data_t> void compute(Cell<data_t> cell) const
    {
        const Coordinates<data_t> c(cell, dx_, params_.star_centre);
        auto vars = evaluate(c.x, c.y);
        counts_->chi += vars.chi < min_chi_;
        counts_->lapse += vars.lapse < min_lapse_;
        if (!(vars.chi >= 100 * min_chi_ &&
              vars.lapse >= 100 * min_lapse_))
            MayDay::Error("EMSKS 2 floor margin below factor 100 on active or ghost cell");
        cell.store_vars(vars);
    }

  private:
    EMSBH_params_t params_;
    EMSKS2Profile profile_;
    double dx_, min_chi_, min_lapse_;
    struct Counts { std::atomic<std::size_t> chi{0}, lapse{0}; };
    std::shared_ptr<Counts> counts_ = std::make_shared<Counts>();
};

#endif
