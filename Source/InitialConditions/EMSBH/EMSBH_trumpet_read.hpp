/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMSBH_TRUMPET_READ_HPP_
#define EMSBH_TRUMPET_READ_HPP_

#include "1D_SOL/EMSTrumpetSolution_read.hpp"
#include "1D_SOL/EMSCTTSolution_read.hpp"
#include <memory>
#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "EMSBHParams.hpp"
#include "EMSCouplingFunction.hpp"
#include "Tensor.hpp"

//! Sets EMSTRUMPET 1 initial data on the cartoon grid.
class EMSBH_trumpet_read
{
  public:
    template <class data_t> struct ems_adm_vars_t
    {
        Tensor<2, data_t, 3> gamma, K;
        Tensor<1, data_t, 3> shift = {0.}, E = {0.}, B = {0.};
        data_t phi, Pi, lapse, chi;
        data_t boost_metric_det;    //! J = w^2 - sinh(eta)^2 alpha_K^2 X^2
        data_t boost_normal_factor; //! w = cosh(eta) - sinh(eta) beta_x
    };

    EMSBH_trumpet_read(EMSBH_params_t a_params_EMSBH,
                       CouplingFunction::params_t a_params_coupling_function,
                       double a_G_Newton, double a_dx, int a_verbosity);
    //! Read the radial solution and verify the run coupling.
    void compute_1d_solution();
    //! Store physical gamma and K for FixSuperposition_metric.
    template <class data_t> void compute(Cell<data_t> a_current_cell) const;
    //! Physical fields for fixture comparisons.
    template <class data_t>
    ems_adm_vars_t<data_t>
    compute_ems_adm_vars(data_t a_x, data_t a_y, data_t a_z, double a_mass,
                         double a_center, double a_rapidity) const;
    template <class data_t>
    ems_adm_vars_t<data_t>
    compute_binary_ems_adm_vars(data_t a_x, data_t a_y, double a_mass,
                                double a_separation, double a_rapidity,
                                data_t a_z = data_t(0), int a_panel = 0) const;
    //! Conformal fields after the same physical metric finalization.
    template <class data_t>
    CCZ4CartoonVars::VarsWithGauge<data_t>
    compute_single_bh_vars(data_t a_x, data_t a_y, double a_mass,
                           double a_center, double a_rapidity) const;
    template <class data_t>
    CCZ4CartoonVars::VarsWithGauge<data_t>
    compute_binary_bh_vars(data_t a_x, data_t a_y, double a_mass,
                           double a_separation, double a_rapidity,
                           int a_panel = 0) const;

    const EMSCTTSolution_read &ctt_solution() const { return *m_ctt_sol; }

    EMSTrumpetSolution_read m_1d_sol;

  private:
    EMSBH_params_t m_params_EMSBH;
    CouplingFunction::params_t m_params_coupling_function;
    double m_G_Newton, m_dx;
    int m_verbosity;
    std::shared_ptr<const EMSCTTSolution_read> m_ctt_sol;

    static Tensor<2, double, 4>
    lorentz_transform(const Tensor<2, double, 4> &a_tensor,
                      const Tensor<2, double, 4> &a_lorentz_matrix);
    template <class data_t>
    static CCZ4CartoonVars::VarsWithGauge<data_t>
    conformal_decomposition(const ems_adm_vars_t<data_t> &a_adm_vars,
                            bool a_is_binary);
};

#include "EMSBH_trumpet_read.impl.hpp"

#endif /* EMSBH_TRUMPET_READ_HPP_ */
