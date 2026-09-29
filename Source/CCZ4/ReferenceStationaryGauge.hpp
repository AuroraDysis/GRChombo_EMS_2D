#ifndef REFERENCE_STATIONARY_GAUGE_HPP_
#define REFERENCE_STATIONARY_GAUGE_HPP_

#include "FArrayBox.H"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"
#include <cmath>

class ReferenceStationaryGauge
{
  public:
    enum { alpha_star, q_star, K_star, beta_x, beta_y, taper, num_reference };
    struct params_t : MovingPunctureGauge::params_t
    {
        const FArrayBox *reference = nullptr;
        const FArrayBox *relative_lapse = nullptr;
        double mass = 1.;
        double onepluslog_n = 2.;
        bool onepluslog = false;
    };

    explicit ReferenceStationaryGauge(const params_t &params) : m_params(params) {}

    double reference(const IntVect &iv, int component) const
    {
        if (!m_params.reference)
            MayDay::Error("reference_stationary cache is absent");
        return (*m_params.reference)(iv, component);
    }

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void rhs_gauge(vars_t<data_t> &rhs, const vars_t<data_t> &,
                   const vars_t<Tensor<1, data_t>> &,
                   const diff2_vars_t<Tensor<2, data_t>> &,
                   const vars_t<data_t> &) const
    {
        rhs.lapse = 0.;
        FOR(i) { rhs.shift[i] = 0.; rhs.B[i] = 0.; }
    }

    template <class data_t, template <typename> class vars_t>
    void set_lapse_rhs(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                       const vars_t<data_t> &advec, const IntVect &iv) const
    {
        const double alpha_ref = reference(iv, alpha_star);
        const double principal = m_params.onepluslog ?
            m_params.onepluslog_n * vars.lapse : vars.lapse * vars.lapse;
        rhs.lapse = advec.lapse - principal *
            (vars.K - reference(iv, K_star) - 2 * vars.Theta) -
            reference(iv, q_star) * vars.lapse -
            vars.lapse / m_params.mass * std::log(vars.lapse / alpha_ref);
    }

    bool relative_lapse() const { return m_params.relative_lapse != nullptr; }

    template <class data_t, template <typename> class vars_t, class deriv_t>
    void set_relative_lapse_rhs(vars_t<data_t> &rhs,
                                const vars_t<data_t> &vars,
                                const IntVect &iv, const deriv_t &deriv,
                                double sigma) const
    {
        const auto terms = deriv.scalar_advection_dissipation(
            *m_params.relative_lapse, iv, vars.shift,
            sigma * reference(iv, taper));
        const double u = (*m_params.relative_lapse)(iv, 0);
        const double k = vars.K - reference(iv, K_star) - 2 * vars.Theta;
        rhs.lapse = vars.lapse *
            (terms.first + terms.second -
             (m_params.onepluslog ? m_params.onepluslog_n : vars.lapse) * k -
             u / m_params.mass);
    }

  private:
    params_t m_params;
};

#endif
