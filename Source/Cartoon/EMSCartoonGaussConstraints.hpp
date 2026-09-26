/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMSCARTOONGAUSSCONSTRAINTS_HPP_
#define EMSCARTOONGAUSSCONSTRAINTS_HPP_

#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "EMSCouplingFunction.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"

//! Signed D_i(F E^i) and D_i B^i on the cartoon plane.
class EMSCartoonGaussConstraints
{
    const FourthOrderDerivatives m_deriv;
    double m_dx;
    CouplingFunction::params_t m_coupling_params;

  public:
    EMSCartoonGaussConstraints(double a_dx,
                               CouplingFunction::params_t a_coupling_params)
        : m_deriv(a_dx), m_dx(a_dx), m_coupling_params(a_coupling_params)
    {
    }

    template <class data_t> void compute(Cell<data_t> a_current_cell) const
    {
        const auto vars =
            a_current_cell.template load_vars<CCZ4CartoonVars::VarsNoGauge>();
        const auto d1 = m_deriv.template diff1<CCZ4CartoonVars::VarsNoGauge>(
            a_current_cell);
        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        const data_t y = Coordinates<data_t>(a_current_cell, m_dx).y;
        Tensor<1, data_t> E = {vars.Ex, vars.Ey};
        Tensor<1, data_t> B = {vars.Bx, vars.By};
        Tensor<2, data_t> dE, dB;
        dE[0][0] = d1.Ex[0];
        dE[0][1] = d1.Ey[0];
        dE[1][0] = d1.Ex[1];
        dE[1][1] = d1.Ey[1];
        dB[0][0] = d1.Bx[0];
        dB[0][1] = d1.By[0];
        dB[1][0] = d1.Bx[1];
        dB[1][1] = d1.By[1];
        const auto compute_divergence =
            [&](const Tensor<1, data_t> &a_field,
                const Tensor<2, data_t> &a_derivative)
        {
            data_t divergence = 0;
            FOR(a)
            {
                data_t q = 0, dq = 0, tr = 0;
                FOR(j)
                {
                    q += h_UU[a][j] * a_field[j];
                    dq += h_UU[a][j] * a_derivative[a][j];
                    FOR(k)
                    {
                        tr += h_UU[j][k] * d1.h[j][k][a];
                        FOR(l)
                        dq -= h_UU[a][j] * d1.h[j][k][a] * h_UU[k][l] *
                              a_field[l];
                    }
                }
                const data_t up = vars.chi * q;
                const data_t log_volume = 0.5 * tr +
                                          0.5 * d1.hww[a] / vars.hww -
                                          1.5 * d1.chi[a] / vars.chi;
                divergence += vars.chi * dq + d1.chi[a] * q + up * log_volume;
            }
            data_t uy = 0;
            FOR(j) uy += vars.chi * h_UU[1][j] * a_field[j];
            return divergence + uy / y;
        };
        const data_t F =
            exp(-2 * m_coupling_params.alpha *
                (m_coupling_params.f0 + m_coupling_params.f1 * vars.phi +
                 m_coupling_params.f2 * vars.phi * vars.phi));
        data_t electric_phi_gradient = 0;
        FOR(a, j)
        electric_phi_gradient += vars.chi * h_UU[a][j] * E[j] * d1.phi[a];
        const data_t electric_gauss =
            F *
            (compute_divergence(E, dE) -
             2 * m_coupling_params.alpha *
                 (m_coupling_params.f1 + 2 * m_coupling_params.f2 * vars.phi) *
                 electric_phi_gradient);
        a_current_cell.store_vars(electric_gauss, c_GaussE);
        a_current_cell.store_vars(compute_divergence(B, dB), c_GaussB);
    }
};

#endif /* EMSCARTOONGAUSSCONSTRAINTS_HPP_ */
