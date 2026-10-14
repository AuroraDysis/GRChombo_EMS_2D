/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMSMATTERWEYLOMSCALAR_HPP_
#define EMSMATTERWEYLOMSCALAR_HPP_

#include "EMSCouplingFunction.hpp"
#include "CCZ4Cartoon.hpp" // compute_EMS_EM_tensor: the evolution's own EMS stress tensor
#include "WeylOmScalar.hpp"

//! WeylOmScalar with the matter term kept. WeylOmScalar takes the two normal-normal projections of the spacetime Riemann tensor from
//! the vacuum field equations; the contracted Gauss relation gives in general
//!     4R_{i n j n} = R_ij + K K_ij - K_ik K^k_j - 4R_ij,   4R_ij = 8 pi G [S_ij - gamma_ij (S - rho) / 2],
//! and the same for the ww component (the contracted Gauss relation, e.g. Gourgoulhon, 3+1 Formalism in General Relativity, ch. 2-4).
//! The Gauss and Codazzi parts are exact already.
//! In Weyl4_Re = Om22 - Omww the trace terms cancel, so the change against WeylOmScalar is -2 pi G (S_vv - S_ww / gamma_ww), v the unit
//! tetrad vector m2: second order in perturbations of a spherical scalar and radial electric background.
class EMSMatterWeylOmScalar : public WeylOmScalar
{
  public:
    EMSMatterWeylOmScalar(const std::array<double, CH_SPACEDIM> a_center,
                          const double a_dx,
                          const CouplingFunction::params_t &a_coupling,
                          const double a_G_Newton)
        : WeylOmScalar(a_center, a_dx), m_coupling(a_coupling),
          m_G_Newton(a_G_Newton)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
        const Coordinates<data_t> coords(current_cell, m_deriv.m_dx, m_center);
        Riemann_t<data_t> riemann = compute_riemann_tensor(vars, d1, d2, coords);
        add_matter(riemann, vars, d1);
        const NPScalarMatrix_t<data_t> out =
            compute_Weyl_Om(riemann, vars, d1, d2, coords);
        current_cell.store_vars(out.Om22, c_Weyl4_Re);
        current_cell.store_vars(0., c_Weyl4_Im);
    }

  protected:
    const CouplingFunction m_coupling;
    const double m_G_Newton;

    template <class data_t>
    void add_matter(Riemann_t<data_t> &riemann, const Vars<data_t> &vars,
                    const Vars<Tensor<1, data_t>> &d1) const
    {
        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);
        const data_t h_UU_ww = 1. / vars.hww;
        const emtensorCartoon_t<data_t> em = compute_EMS_EM_tensor(
            vars, d1, h_UU, h_UU_ww, chris, GR_SPACEDIM - CH_SPACEDIM,
            m_coupling);
        // the chi floor of compute_riemann_tensor, so both use one physical metric
        const data_t one_over_chi = 1. / simd_max(vars.chi, 1e-8);
        const data_t half_trace = 0.5 * (em.S - em.rho);
        const double eight_pi_G = 8. * M_PI * m_G_Newton;
        FOR(i, j)
        {
            riemann.riemann_0d0d[i][j] -=
                eight_pi_G *
                (em.Sij[i][j] - half_trace * vars.h[i][j] * one_over_chi);
        }
        riemann.riemann_0w0w -=
            eight_pi_G * (em.Sww - half_trace * vars.hww * one_over_chi);
    }
};

#endif /* EMSMATTERWEYLOMSCALAR_HPP_ */
