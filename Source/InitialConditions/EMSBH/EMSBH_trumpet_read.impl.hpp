/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EMSBH_TRUMPET_READ_HPP_)
#error "This file should only be included through EMSBH_trumpet_read.hpp"
#endif
#ifndef EMSBH_TRUMPET_READ_IMPL_HPP_
#define EMSBH_TRUMPET_READ_IMPL_HPP_

#include "MayDay.H"
#include "TensorAlgebra.hpp"
#include "parstream.H"
#include <cmath>

inline EMSBH_trumpet_read::EMSBH_trumpet_read(
    EMSBH_params_t a_params_EMSBH,
    CouplingFunction::params_t a_params_coupling_function, double a_G_Newton,
    double a_dx, int a_verbosity)
    : m_params_EMSBH(a_params_EMSBH),
      m_params_coupling_function(a_params_coupling_function),
      m_G_Newton(a_G_Newton), m_dx(a_dx), m_verbosity(a_verbosity)
{
    if (m_params_EMSBH.star_centre[1] != 0)
        MayDay::Error("EMSTRUMPET object must lie on the cartoon axis");
}

inline void EMSBH_trumpet_read::compute_1d_solution()
{
    m_1d_sol.main(m_params_EMSBH.data_path);
    if (m_1d_sol.get_coupling_parameters() !=
        std::array<double, 4>{
            m_params_coupling_function.alpha, m_params_coupling_function.f0,
            m_params_coupling_function.f1, m_params_coupling_function.f2})
        MayDay::Error("EMSTRUMPET coupling differs from run parameters");
    if (m_verbosity)
        pout() << "EMSTRUMPET radial solution ready" << std::endl;
}

inline Tensor<2, double, 4> EMSBH_trumpet_read::lorentz_transform(
    const Tensor<2, double, 4> &a_tensor,
    const Tensor<2, double, 4> &a_lorentz_matrix)
{
    Tensor<2, double, 4> transformed_tensor = {0.};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int m = 0; m < 4; ++m)
                for (int n = 0; n < 4; ++n)
                    transformed_tensor[i][j] += a_lorentz_matrix[m][i] *
                                                a_tensor[m][n] *
                                                a_lorentz_matrix[n][j];
    return transformed_tensor;
}

template <class data_t>
inline typename EMSBH_trumpet_read::template ems_adm_vars_t<data_t>
EMSBH_trumpet_read::compute_ems_adm_vars(data_t a_x, data_t a_y, data_t a_z,
                                         double a_mass, double a_center,
                                         double a_rapidity) const
{
    if (!(a_mass > 0) || !std::isfinite(a_mass))
        MayDay::Error("EMSTRUMPET mass must be positive and finite");
    const double c = std::cosh(a_rapidity), sh = std::sinh(a_rapidity);
    const data_t xr = c * (a_x - a_center);
    const data_t R = std::sqrt(xr * xr + a_y * a_y + a_z * a_z) / a_mass;
    if (!(R > 0))
        MayDay::Error(
            "EMSTRUMPET object center is the limiting cylinder (R = 0)");
    const auto radial_vars = m_1d_sol.compute_radial_vars(R);
    const double P = 1 / radial_vars.X, alpha = radial_vars.lapse,
                 b = radial_vars.shift_R;
    Tensor<1, data_t, 3> l = {xr / (a_mass * R), a_y / (a_mass * R),
                              a_z / (a_mass * R)};
    const data_t bx = b * l[0], w = c - sh * bx;
    Tensor<2, double, 4> ghat = {0.};
    ghat[0][0] = b * b - (alpha * radial_vars.X) * (alpha * radial_vars.X);
    for (int i = 0; i < 3; ++i)
    {
        ghat[0][i + 1] = ghat[i + 1][0] = b * l[i];
        ghat[i + 1][i + 1] = 1;
    }
    Tensor<2, double, 4> lorentz_matrix = {0.};
    lorentz_matrix[0][0] = lorentz_matrix[1][1] = c;
    lorentz_matrix[0][1] = lorentz_matrix[1][0] = -sh;
    lorentz_matrix[2][2] = lorentz_matrix[3][3] = 1;
    const auto boosted_metric = lorentz_transform(ghat, lorentz_matrix);
    Tensor<2, data_t, 3> hatgamma = {0.};
    Tensor<1, data_t, 3> gt = {0.};
    for (int i = 0; i < 3; ++i)
    {
        gt[i] = boosted_metric[i + 1][0];
        for (int j = 0; j < 3; ++j)
            hatgamma[i][j] = boosted_metric[i + 1][j + 1];
    }
    const data_t J =
        w * w - sh * sh * alpha * alpha * radial_vars.X * radial_vars.X;
    if (!(J > 0 && w > 0))
        MayDay::Error("EMSTRUMPET boosted data loses positive normal margin");
    ems_adm_vars_t<data_t> adm_vars{};
    adm_vars.boost_metric_det = J;
    adm_vars.boost_normal_factor = w;
    const data_t det_hatgamma =
        TensorAlgebra::compute_determinant_sym(hatgamma);
    if (!(det_hatgamma > 0))
        MayDay::Error("EMSTRUMPET nonpositive spatial metric determinant");
    const auto inverse = TensorAlgebra::compute_inverse_sym(hatgamma);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            adm_vars.shift[i] += inverse[i][j] * gt[j];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            adm_vars.gamma[i][j] = P * P * hatgamma[i][j];
    Tensor<2, data_t, 3> K0 = {0.};
    Tensor<2, double, 4> G4 = {0.};
    const double alphar = radial_vars.dlapse_dR / a_mass,
                 br = radial_vars.dshift_dR / a_mass;
    const double pr = radial_vars.dP_dR / (a_mass * P),
                 k = radial_vars.K_amplitude / a_mass;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            const double dij = i == j ? 1 : 0;
            K0[i][j] = P * P * k * (3 * l[i] * l[j] - dij);
            const double Gammax =
                pr * ((i == 0 ? l[j] : 0) + (j == 0 ? l[i] : 0) - dij * l[0]);
            G4[i + 1][j + 1] = -w * K0[i][j] + sh * alpha * Gammax;
        }
    for (int i = 0; i < 3; ++i)
    {
        const double d0 = i == 0 ? 1 : 0;
        const double Dbx = br * l[i] * l[0] +
                           b / (a_mass * R) * (d0 - l[i] * l[0]) + b * pr * d0;
        const double Kxi = k * (3 * l[0] * l[i] - d0);
        G4[0][i + 1] = G4[i + 1][0] =
            w * (alphar * l[i] - 2 * b * P * P * k * l[i]) + sh * alpha * Dbx -
            sh * alpha * alpha * Kxi;
    }
    G4[0][0] =
        w * (b * alphar - 2 * b * b * P * P * k) +
        sh * alpha * alpha * radial_vars.X * radial_vars.X * alphar * l[0] +
        sh * alpha * (b * br * l[0] + b * b * pr * l[0]) -
        4 * sh * alpha * alpha * b * k * l[0];
    const auto boostedG = lorentz_transform(G4, lorentz_matrix);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            adm_vars.K[i][j] = -boostedG[i + 1][j + 1] / std::sqrt(J);
    const double phis =
        -m_1d_sol.get_chebyshev_value(2, radial_vars.s) +
        (1 - radial_vars.s) * m_1d_sol.get_chebyshev_value(2, radial_vars.s, 1);
    const data_t phix =
        phis *
        (m_1d_sol.get_nu() * radial_vars.s * (1 - radial_vars.s) /
         (R * radial_vars.compactification_factor)) *
        l[0] / a_mass;
    adm_vars.phi = radial_vars.phi;
    adm_vars.Pi = (w * radial_vars.Pi / a_mass +
                   alpha * sh * radial_vars.X * radial_vars.X * phix) /
                  std::sqrt(J);
    Tensor<2, double, 4> F = {0.};
    for (int i = 0; i < 3; ++i)
    {
        F[i + 1][0] = radial_vars.E_R * l[i] / a_mass;
        F[0][i + 1] = -F[i + 1][0];
    }
    const auto boosted_maxwell_tensor = lorentz_transform(F, lorentz_matrix);
    for (int i = 0; i < 3; ++i)
    {
        double ei = boosted_maxwell_tensor[i + 1][0];
        for (int j = 0; j < 3; ++j)
            ei -= adm_vars.shift[j] * boosted_maxwell_tensor[i + 1][j + 1];
        adm_vars.E[i] = std::sqrt(J) * ei;
    }
    Tensor<1, data_t, 3> H = {boosted_maxwell_tensor[2][3],
                              boosted_maxwell_tensor[3][1],
                              boosted_maxwell_tensor[1][2]};
    Tensor<1, data_t, 3> h = {0.};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            h[i] += hatgamma[i][j] * H[j];
    for (int i = 0; i < 3; ++i)
        adm_vars.B[i] = alpha * radial_vars.X / std::sqrt(J) * h[i];
    adm_vars.lapse = alpha / std::sqrt(J);
    adm_vars.chi = radial_vars.X * radial_vars.X / std::cbrt(J);
    return adm_vars;
}

template <class data_t>
inline typename EMSBH_trumpet_read::template ems_adm_vars_t<data_t>
EMSBH_trumpet_read::compute_binary_ems_adm_vars(data_t a_x, data_t a_y,
                                                double a_mass,
                                                double a_separation,
                                                double a_rapidity) const
{
    if (!(a_separation > 0))
        MayDay::Error("EMSTRUMPET separation must be positive");
    const auto left_bh = compute_ems_adm_vars(a_x, a_y, data_t(0), a_mass,
                                              -a_separation / 2, a_rapidity);
    const auto right_bh = compute_ems_adm_vars(a_x, a_y, data_t(0), a_mass,
                                               a_separation / 2, -a_rapidity);
    ems_adm_vars_t<data_t> superposed_vars{};
    Tensor<1, data_t, 3> electric_density = {0.}, magnetic_density = {0.};
    const double c = std::cosh(a_rapidity);
    const double charge = a_mass * m_1d_sol.get_native_charge();
    const ems_adm_vars_t<data_t> *holes[2] = {&left_bh, &right_bh};
    const double centers[2] = {-a_separation / 2, a_separation / 2};
    for (int hole = 0; hole < 2; ++hole)
    {
        const data_t X = a_x - centers[hole];
        const data_t R = std::sqrt(c * c * X * X + a_y * a_y);
        const Tensor<1, data_t, 3> position = {X, a_y, data_t(0)};
        const auto inverse =
            TensorAlgebra::compute_inverse_sym(holes[hole]->gamma);
        const data_t root_det = std::sqrt(
            TensorAlgebra::compute_determinant_sym(holes[hole]->gamma));
        for (int i = 0; i < 3; ++i)
        {
            electric_density[i] += charge * c * position[i] / (R * R * R);
            for (int j = 0; j < 3; ++j)
                magnetic_density[i] +=
                    root_det * inverse[i][j] * holes[hole]->B[j];
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        superposed_vars.shift[i] = left_bh.shift[i] + right_bh.shift[i];
        for (int j = 0; j < 3; ++j)
        {
            superposed_vars.gamma[i][j] =
                left_bh.gamma[i][j] + right_bh.gamma[i][j];
            superposed_vars.K[i][j] = left_bh.K[i][j] + right_bh.K[i][j];
        }
    }
    superposed_vars.phi = m_1d_sol.get_phi_inf() +
                          (left_bh.phi - m_1d_sol.get_phi_inf()) +
                          (right_bh.phi - m_1d_sol.get_phi_inf());
    superposed_vars.Pi = left_bh.Pi + right_bh.Pi;
    auto physical_gamma = superposed_vars.gamma;
    for (int i = 0; i < 3; ++i)
        physical_gamma[i][i] -= 1;
    const data_t det = TensorAlgebra::compute_determinant_sym(physical_gamma);
    if (!(det > 0))
        MayDay::Error("EMSTRUMPET nonpositive physical metric determinant");
    const data_t root_det = std::sqrt(det);
    const data_t phi = superposed_vars.phi;
    const data_t coupling = std::exp(
        -2 * m_params_coupling_function.alpha *
        (m_params_coupling_function.f0 + m_params_coupling_function.f1 * phi +
         m_params_coupling_function.f2 * phi * phi));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            superposed_vars.E[i] += physical_gamma[i][j] * electric_density[j] /
                                    (root_det * coupling);
            superposed_vars.B[i] +=
                physical_gamma[i][j] * magnetic_density[j] / root_det;
        }
    return superposed_vars;
}

template <class data_t>
inline CCZ4CartoonVars::VarsWithGauge<data_t>
EMSBH_trumpet_read::conformal_decomposition(
    const ems_adm_vars_t<data_t> &a_adm_vars, bool a_is_binary)
{
    Tensor<2, data_t, 3> gamma = a_adm_vars.gamma;
    if (a_is_binary)
        for (int i = 0; i < 3; ++i)
            gamma[i][i] -= 1;
    const data_t det = TensorAlgebra::compute_determinant_sym(gamma);
    if (!(det > 0))
        MayDay::Error("EMSTRUMPET nonpositive physical metric determinant");
    const data_t chi = 1 / std::cbrt(det);
    const auto inverse = TensorAlgebra::compute_inverse_sym(gamma);
    data_t K = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            K += inverse[i][j] * a_adm_vars.K[i][j];
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars, 0.);
    vars.chi = chi;
    vars.K = K;
    FOR(i, j)
    {
        vars.h[i][j] = chi * gamma[i][j];
        vars.A[i][j] = chi * (a_adm_vars.K[i][j] - K * gamma[i][j] / 3);
    }
    vars.hww = chi * gamma[2][2];
    vars.Aww = chi * (a_adm_vars.K[2][2] - K * gamma[2][2] / 3);
    vars.lapse = std::sqrt(chi);
    vars.shift[0] = a_adm_vars.shift[0];
    vars.shift[1] = a_adm_vars.shift[1];
    vars.phi = a_adm_vars.phi;
    vars.Pi = a_adm_vars.Pi;
    vars.Ex = a_adm_vars.E[0];
    vars.Ey = a_adm_vars.E[1];
    vars.Ez = a_adm_vars.E[2];
    vars.Bx = a_adm_vars.B[0];
    vars.By = a_adm_vars.B[1];
    vars.Bz = a_adm_vars.B[2];
    return vars;
}

template <class data_t>
inline CCZ4CartoonVars::VarsWithGauge<data_t>
EMSBH_trumpet_read::compute_single_bh_vars(data_t a_x, data_t a_y,
                                           double a_mass, double a_center,
                                           double a_rapidity) const
{
    return conformal_decomposition(
        compute_ems_adm_vars(a_x, a_y, data_t(0), a_mass, a_center, a_rapidity),
        false);
}

template <class data_t>
inline CCZ4CartoonVars::VarsWithGauge<data_t>
EMSBH_trumpet_read::compute_binary_bh_vars(data_t a_x, data_t a_y,
                                           double a_mass, double a_separation,
                                           double a_rapidity) const
{
    return conformal_decomposition(
        compute_binary_ems_adm_vars(a_x, a_y, a_mass, a_separation, a_rapidity),
        true);
}

template <class data_t>
inline void EMSBH_trumpet_read::compute(Cell<data_t> a_current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars{};
    a_current_cell.load_vars(vars);
    const Coordinates<data_t> coords(a_current_cell, m_dx,
                                     m_params_EMSBH.star_centre);
    const double eta = m_params_EMSBH.boosted ? m_params_EMSBH.rapidity : 0.0;
    const auto adm_vars =
        m_params_EMSBH.binary
            ? compute_binary_ems_adm_vars(coords.x, coords.y,
                                          m_params_EMSBH.bh_mass,
                                          m_params_EMSBH.separation, eta)
            : compute_ems_adm_vars(coords.x, coords.y, data_t(0),
                                   m_params_EMSBH.bh_mass, 0.0, eta);
    FOR(i, j)
    {
        vars.h[i][j] = adm_vars.gamma[i][j];
        vars.A[i][j] = adm_vars.K[i][j];
    }
    vars.hww = adm_vars.gamma[2][2];
    vars.Aww = adm_vars.K[2][2];
    vars.shift[0] = adm_vars.shift[0];
    vars.shift[1] = adm_vars.shift[1];
    vars.phi = adm_vars.phi;
    vars.Pi = adm_vars.Pi;
    vars.Ex = adm_vars.E[0];
    vars.Ey = adm_vars.E[1];
    vars.Ez = adm_vars.E[2];
    vars.Bx = adm_vars.B[0];
    vars.By = adm_vars.B[1];
    vars.Bz = adm_vars.B[2];
    vars.lapse = adm_vars.lapse;
    vars.Theta = 0;
    vars.Xi = 0;
    vars.Lambda = 0;
    a_current_cell.store_vars(vars);
}

#endif /* EMSBH_TRUMPET_READ_IMPL_HPP_ */
