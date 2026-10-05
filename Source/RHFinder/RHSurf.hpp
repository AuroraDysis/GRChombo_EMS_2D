#ifndef RHSURF_HPP_
#define RHSURF_HPP_

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "RHBandedMatrixInverter.hpp"

// Loop over the interior (non-ghost) array indices of a surface.
// ii runs from m_NG to m_n+m_NG-1 inclusive.
#define RH_FOR_INTERIOR(ii) for (int (ii) = m_NG; (ii) < m_n + m_NG; ++(ii))

// Star-shaped surface in the poloidal half-plane, parameterised as r = f(theta).
// Arrays have n + 2*NG entries: NG ghost cells on each side of the n interior points.
// theta[i] = (i + 0.5 - NG) * pi / n,  f[i] initialised to 1.
// data is cell centred
class RHSurf
{
  public:
    std::vector<double> m_centre; // [x, y] centre in grid coords
    int m_n;                      // number of interior angular points
    int m_NG;                     // ghost cells on each side
    int m_index;                  // surface index (order in params file)
    std::vector<double> m_theta;      // size n + 2*NG
    std::vector<double> m_f;          // size n + 2*NG, initialised to 1
    std::vector<double> m_Expansion;  // Theta_plus at each interior point, filled before chasing
    enum class SolverState { DORMANT, FAR, MEDIUM, CLOSE, FOUND };
    SolverState m_state = SolverState::DORMANT;
    bool m_dead         = false; // set on fatal error; surface is never updated again
    int m_level          = 0;    // AMR level this surface is updated on
    int m_time_step_freq = 1;    // chase iterations per update call
    double m_newton_crit = 0.0;  // switch to Newton below this expansion_error (0 = disabled)
    // Theta_plus uses df/d2f at +/-2; all field data/Cartesian derivatives
    // in DivS and KSSmK are point-local. Even polar reflection preserves this band.
    static constexpr int expansion_half_width = 2;
    bool m_use_newton = false;
    std::vector<double> m_cell_dx; // largest coordinate spacing at each interior point
    double m_rmin = 0.0001;
    double m_rmax = 10.0;
    double m_chase_speed = 1.0; // multiplier on the courant chase step
    double m_start_time  = 0.0; // simulation time before which the surface is dormant
    double m_d_theta;          // pi / m_n, angular spacing

    // field values at surface points (size n + 2*NG)
    std::vector<double> m_chi;
    std::vector<double> m_K;
    std::vector<double> m_A11, m_A12, m_A22, m_Aww;
    std::vector<double> m_h11, m_h12, m_h22, m_hww;

    // Cartesian spatial derivatives of chi and conformal metric
    std::vector<double> m_dx_chi, m_dy_chi;
    std::vector<double> m_dx_h11, m_dy_h11;
    std::vector<double> m_dx_h12, m_dy_h12;
    std::vector<double> m_dx_h22, m_dy_h22;
    std::vector<double> m_dx_hww, m_dy_hww;

    // EM and scalar field values at surface points
    std::vector<double> m_Ex, m_Ey, m_phi;
    // EMS coupling params: coupling = exp(-2*alpha*(f0+f1*phi+f2*phi^2))
    double m_alpha = 0., m_f0 = 0., m_f1 = 0., m_f2 = 0.;

    RHSurf(int a_n, int a_NG, std::vector<double> a_origin, int a_index = 0)
        : m_centre(a_origin), m_n(a_n), m_NG(a_NG), m_index(a_index),
          m_d_theta(M_PI / a_n),
          m_theta(a_n + 2 * a_NG, 0.0),
          m_f(a_n + 2 * a_NG, 1.0),
          m_Expansion(a_n + 2 * a_NG, 0.0),
          m_chi(a_n + 2 * a_NG, 0.0),
          m_K(a_n + 2 * a_NG, 0.0),
          m_A11(a_n + 2 * a_NG, 0.0), m_A12(a_n + 2 * a_NG, 0.0),
          m_A22(a_n + 2 * a_NG, 0.0), m_Aww(a_n + 2 * a_NG, 0.0),
          m_h11(a_n + 2 * a_NG, 0.0), m_h12(a_n + 2 * a_NG, 0.0),
          m_h22(a_n + 2 * a_NG, 0.0), m_hww(a_n + 2 * a_NG, 0.0),
          m_dx_chi(a_n + 2 * a_NG, 0.0), m_dy_chi(a_n + 2 * a_NG, 0.0),
          m_dx_h11(a_n + 2 * a_NG, 0.0), m_dy_h11(a_n + 2 * a_NG, 0.0),
          m_dx_h12(a_n + 2 * a_NG, 0.0), m_dy_h12(a_n + 2 * a_NG, 0.0),
          m_dx_h22(a_n + 2 * a_NG, 0.0), m_dy_h22(a_n + 2 * a_NG, 0.0),
          m_dx_hww(a_n + 2 * a_NG, 0.0), m_dy_hww(a_n + 2 * a_NG, 0.0),
          m_Ex(a_n + 2 * a_NG, 0.0), m_Ey(a_n + 2 * a_NG, 0.0),
          m_phi(a_n + 2 * a_NG, 0.0)
    {
        for (int i = 0; i < m_n + 2 * m_NG; ++i)
            m_theta[i] = (i + 0.5 - m_NG) * m_d_theta;
    }

    void reset(double a_r = 5.0)
    {
        std::fill(m_f.begin(), m_f.end(), a_r);
        m_state = SolverState::FAR;
    }

    double df(int ii) const
    {
        return (-m_f[ii+2] + 8.0*m_f[ii+1] - 8.0*m_f[ii-1] + m_f[ii-2])
               / (12.0 * m_d_theta);
    }

    double d2f(int ii) const
    {
        return (-m_f[ii+2] + 16.0*m_f[ii+1] - 30.0*m_f[ii]
                + 16.0*m_f[ii-1] - m_f[ii-2])
               / (12.0 * m_d_theta * m_d_theta);
    }

    // Logarithmic radial derivative of the surface: q = f'(theta) / f(theta).
    double q(int ii) const
    {
        return df(ii) / m_f[ii];
    }


    //////////////// Lengths ////////////

    // Proper poloidal arc-length element at point ii: sqrt(gamma_theta_theta) dtheta
    double dl(int ii) const
    {
        const double f   = m_f[ii];
        const double fp  = df(ii);
        const double th  = m_theta[ii];
        const double dxi_dth  = fp * std::cos(th) - f * std::sin(th);
        const double drho_dth = fp * std::sin(th) + f * std::cos(th);
        const double h_th_th = m_h11[ii] * dxi_dth  * dxi_dth
                             + 2.0 * m_h12[ii] * dxi_dth * drho_dth
                             + m_h22[ii] * drho_dth * drho_dth;
        return std::sqrt(h_th_th / m_chi[ii]) * m_d_theta;
    }

    // Proper azimuthal circumference at the equator (theta = pi/2)
    double equator_perim() const
    {
        const int i_eq = m_NG + m_n / 2;
        return 2.0 * M_PI * m_f[i_eq] * std::sqrt(m_hww[i_eq] / m_chi[i_eq]);
    }

    // Proper polar great-circle length (full loop through both poles)
    double polar_perim() const
    {
        double polarperim = 0.;
        RH_FOR_INTERIOR(ii)
        {
            polarperim += dl(ii);
        }
        return 2.0 * polarperim;
    }


    //////////////// Areas //////////////

    double dA(int ii) const
    {
        const double f   = m_f[ii];
        const double fp  = df(ii);
        const double th  = m_theta[ii];

        // tangent vector components along the surface in (xi, rho) poloidal space
        const double dxi_dth  = fp * std::cos(th) - f * std::sin(th);
        const double drho_dth = fp * std::sin(th) + f * std::cos(th);

        // conformal induced metric along theta direction: h_tilde_{theta theta}
        const double h_th_th = m_h11[ii] * dxi_dth  * dxi_dth
                             + 2.0 * m_h12[ii] * dxi_dth * drho_dth
                             + m_h22[ii] * drho_dth * drho_dth;

        // dA = 2pi * rho/chi * sqrt(h_th_th * hww) * dtheta
        // rho = f sin(theta), and gamma_phiphi = hww/chi * rho^2
        return 2.0 * M_PI * f * std::sin(th) / m_chi[ii]
               * std::sqrt(h_th_th * m_hww[ii]) * m_d_theta;
    }

    double Area() const
    {
        double sum = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            sum += dA(ii);
        }
        return sum;
    }

    double average_chi() const
    {
        double sum = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            sum += m_chi[ii];
        }
        return sum / m_n;
    }

    double average_f() const
    {
        double sum = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            sum += m_f[ii];
        }
        return sum / m_n;
    }

    // Area-averaged expansions
    double average_Theta_plus() const
    {
        double sum_Theta_dA = 0.0, sum_dA = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            const double dAii = dA(ii);
            sum_Theta_dA += Theta_plus(ii) * dAii;
            sum_dA       += dAii;
        }
        return sum_Theta_dA / sum_dA;
    }

    double average_Theta_minus() const
    {
        double sum_Theta_dA = 0.0, sum_dA = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            const double dAii = dA(ii);
            sum_Theta_dA += Theta_minus(ii) * dAii;
            sum_dA       += dAii;
        }
        return sum_Theta_dA / sum_dA;
    }

    double average_Theta() const { return average_Theta_plus(); }

    // Area-averaged |Theta+|^2 — used as convergence error in RHUnion.
    double expansion_error() const
    {
        double sum_sq = 0.0, sum_dA = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            const double dAii = dA(ii);
            sum_sq += Theta_plus(ii) * Theta_plus(ii) * dAii;
            sum_dA += dAii;
        }
        return sum_sq / sum_dA;
    }



    //////////////// Normal Vecs //////////////

    // Unit outward normal in spherical polar components (s^r_hat, s^theta_hat).
    // Follows from grad(r - f(theta)) / |grad(r - f(theta))| in the spherical polar metric.
    // std::array<double, 2> s_polar(int ii) const
    // {
    //     const double dtheta = M_PI / m_n;
    //     const double f  = m_f[ii];
    //     const double fp = (m_f[ii + 1] - m_f[ii - 1]) / (2.0 * dtheta);
    //     const double mag = std::sqrt(f * f + fp * fp);
    //     return { f / mag, -fp / mag };
    // }

    // un normalised curved space co-vector
    std::array<double, 3> s_heavy(int ii) const
    {
        // surf coords
        const double f = m_f[ii];
        const double th = m_theta[ii];

        // cart coords
        const double x = f * std::cos(th);
        const double y = f * std::sin(th);

        // conventional z axis normal, chombo x
        const double s_x = (x + q(ii) * y)/f;
        // conventional rho axis normal, chombo y
        const double s_y = (y - q(ii) * x)/f;

        return { s_x, s_y, 0 };
    }

    // un normalised second form
    // partial_i partial_j F 
    // F is the surface r-f(theta)
    std::array<std::array<double, 3>, 3> didjF(int ii) const
    {
        std::array<std::array<double, 3>, 3> out = {{ {0., 0., 0.},
                                                      {0., 0., 0.},
                                                      {0., 0., 0.} }};
        // surf coords
        const double f = m_f[ii];
        const double fpp = d2f(ii);
        const double f3 = f*f*f;
        const double th = m_theta[ii];

        // cart coords
        const double x = f * std::cos(th);
        const double y = f * std::sin(th);
        // z = 0

        // algebraic object needed later
        const double thing = 1. - f*fpp/(x*x+y*y);

        // evaluate d_ij F | i used mathematica
        // it is evalutaed on y=0, but calculated with dynamic z

        out[0][0] = y/f3  * (-2.*x*q(ii) + y*thing);
        out[1][1] = x/f3  * (2.*y*q(ii) + x*thing);
        out[0][1] = 1./f3 * ((x*x-y*y)*q(ii) - x*y*thing);
        out[1][0] = out[0][1]; // symmetry

        // cartoon term
        out[2][2] = 1./f  - x*(x*x+y*y)*q(ii)/(y*f3);

        return out;
    }

    // Inverse conformal metric h^{ij} at surface point ii.
    // Indices: 0=x, 1=y, 2=w (azimuthal).
    std::array<std::array<double, 3>, 3> h_inv(int ii) const
    {
        std::array<std::array<double, 3>, 3> out = {{ {0., 0., 0.},
                                                      {0., 0., 0.},
                                                      {0., 0., 0.} }};

        // const double det2 = m_h11[ii] * m_h22[ii] - m_h12[ii] * m_h12[ii];
        // note that det2 = 1/hww ...

        out[0][0] =  m_h22[ii] * m_hww[ii];
        out[1][1] =  m_h11[ii] * m_hww[ii];
        out[1][0] = -m_h12[ii] * m_hww[ii];
        out[0][1] = -m_h12[ii] * m_hww[ii];
        out[2][2] =  1.0 / m_hww[ii];
        return out;
    }

    // Derivative of inverse physical metric: out[i][j][k] = d_i gamma^{jk}.
    // gamma^{jk} = chi * h^{jk}  =>  d_i gamma^{jk} = (d_i chi) h^{jk} + chi d_i h^{jk}
    // with  d_i h^{jk} = -h^{jl} (d_i h_{lm}) h^{mk}  (matrix derivative identity).
    // Indices: 0=x, 1=y, 2=w. Azimuthal derivatives vanish by axisymmetry.
    std::array<std::array<std::array<double, 3>, 3>, 3> d_gamma_inv(int ii) const
    {
        std::array<std::array<std::array<double, 3>, 3>, 3> out{};
        for (auto &mat : out) for (auto &row : mat) row.fill(0.);

        const auto hU = h_inv(ii);

        // d_i h_{jk}: only i=0 (x) and i=1 (y) are non-zero
        std::array<std::array<std::array<double, 3>, 3>, 3> dh{};
        for (auto &mat : dh) for (auto &row : mat) row.fill(0.);

        // d_i h_jk

        // d/dx
        dh[0][0][0] = m_dx_h11[ii]; 
        dh[0][0][1] = dh[0][1][0] = m_dx_h12[ii];
        dh[0][1][1] = m_dx_h22[ii]; 
        dh[0][2][2] = m_dx_hww[ii];
        // d/dy
        dh[1][0][0] = m_dy_h11[ii]; 
        dh[1][0][1] = dh[1][1][0] = m_dy_h12[ii];
        dh[1][1][1] = m_dy_h22[ii]; 
        dh[1][2][2] = m_dy_hww[ii];

        // cartoon terms in z
        const double y_cyl = m_f[ii] * std::sin(m_theta[ii]);
        dh[2][1][2] = dh[2][2][1] = (m_h22[ii] - m_hww[ii]) / y_cyl;
        dh[2][0][2] = dh[2][2][0] = m_h12[ii] / y_cyl;

        const double chi = m_chi[ii];
        const std::array<double, 3> dchi = { m_dx_chi[ii], m_dy_chi[ii], 0. };

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                {
                    double dhUjk = 0.;
                    for (int l = 0; l < 3; ++l)
                        for (int m = 0; m < 3; ++m)
                            dhUjk -= hU[j][l] * dh[i][l][m] * hU[m][k];
                    out[i][j][k] = dchi[i] * hU[j][k] + chi * dhUjk;
                }
        return out;
    }

    // Unnormalized contravariant normal: raises s_heavy with g^{ij} = h^{ij} / chi.
    // Inverse conformal metric from det h = hww*(h11*h22 - h12^2) = 1:
    //   h^{11} = h22*hww,  h^{12} = -h12*hww,  h^{22} = h11*hww
    std::array<double, 3> s_heavy_vec(int ii) const
    {
        const auto   s  = s_heavy(ii);
        const auto hUU = h_inv(ii);
        const double sx = s[0], sy = s[1];

        return { (hUU[0][0] * sx + hUU[1][0] * sy) * m_chi[ii],
                 (hUU[0][1] * sx + hUU[1][1] * sy) * m_chi[ii],
                  0.};
    }

    // Physical length of the normal: sqrt(s_i s^i) = sqrt(s_heavy . s_heavy_vec).
    double s_length(int ii) const
    {
        const auto s  = s_heavy(ii);
        const auto sv = s_heavy_vec(ii);
        return std::sqrt(s[0] * sv[0] + s[1] * sv[1] + s[2] * sv[2] );
    }

    // Covariant unit normal.
    std::array<double, 3> s_unit(int ii) const
    {
        const auto   s   = s_heavy(ii);
        const double len = s_length(ii);
        return { s[0] / len, s[1] / len, s[2] / len,};
    }

    // Contravariant unit normal.
    std::array<double, 3> s_unit_vec(int ii) const
    {
        const auto   sv  = s_heavy_vec(ii);
        const double len = s_length(ii);
        return { sv[0] / len, sv[1] / len, sv[2] / len };
    }

    // K_ij s^i s^j - K using GRChombo BSSN variables and the physical unit normal.
    // K_ij = (A_ij + K/3 * h_ij) / chi; s_w = 0 since the normal is purely poloidal.
    double KSSmK(int ii) const
    {
        // curved space outward unit normal 
        const auto   shat = s_unit_vec(ii); 
        const double sx = shat[0], sy = shat[1];

        const double Ass    = m_A11[ii] * sx * sx 
                            + 2.0 * m_A12[ii] * sx * sy           
                            + m_A22[ii] * sy * sy;
        return Ass/m_chi[ii] - (2./3.) * m_K[ii];
    }

    double Theta_plus(int ii)  const { return  DivS(ii) + KSSmK(ii); }
    double Theta_minus(int ii) const { return -DivS(ii) + KSSmK(ii); }

    double K0() const {
        double integral_k = 0.;
        RH_FOR_INTERIOR(ii)
        {
            integral_k += DivS(ii) * dA(ii);
        }
        return integral_k/Area();
    }

    double K2() const {
        double integral_k2 = 0.;
        const double k_0 = K0();
        RH_FOR_INTERIOR(ii)
        {
            integral_k2 += pow(DivS(ii) - k_0, 2) * dA(ii);
        }
        return integral_k2/Area();
    }

    double K4() const {
        double integral_k4 = 0.;
        const double k_0 = K0();
        const double k_2 = K2();
        double integrand;
        RH_FOR_INTERIOR(ii)
        {
            integrand = k_2 * k_2;
            integrand -= pow(DivS(ii) - k_0, 2);
            integral_k4 += pow(integrand, 2) * dA(ii);
        }
        return integral_k4/Area();
    }

    // Area-weighted L=2 Legendre projection of the mean curvature:
    // C_2 = integral( K * P_2(cos theta) * dA ) / integral( dA )
    // Zero for a round surface; non-zero encodes pole-vs-equator curvature asymmetry.
    double KP2() const {
        double num = 0., denom = 0.;
        RH_FOR_INTERIOR(ii)
        {
            const double x   = std::cos(m_theta[ii]);
            const double P2  = 0.5 * (3.*x*x - 1.);
            const double dAi = dA(ii);
            num   += DivS(ii) * P2 * dAi;
            denom += dAi;
        }
        return num / denom;
    }

    // Area-weighted L=4 Legendre projection of the mean curvature:
    // C_4 = integral( K * P_4(cos theta) * dA ) / integral( dA )
    double KP4() const {
        double num = 0., denom = 0.;
        RH_FOR_INTERIOR(ii)
        {
            const double x   = std::cos(m_theta[ii]);
            const double x2  = x * x;
            const double P4  = (35.*x2*x2 - 30.*x2 + 3.) / 8.;
            const double dAi = dA(ii);
            num   += DivS(ii) * P4 * dAi;
            denom += dAi;
        }
        return num / denom;
    }

    // D_i s^i = (1/sqrt(gamma)) partial_i (sqrt(gamma) s^i), summed over i = x, y, w.
    // The w-term does NOT vanish: careful of the cartoon condition!
    double DivS(int ii) const
    {

        // vars
        const double chi = m_chi[ii];
        const double L = s_length(ii);
        const double L3 = L*L*L;
        const auto sh = s_heavy(ii); // co-vec, non normalised
        const auto d_ij_F = didjF(ii);
        const auto dgi = d_gamma_inv(ii);

        // inverse metric gamma^{ij} = chi * h^{ij}
        auto ginv = h_inv(ii); 
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                ginv[a][b] *= chi;

        // Div s = root-g^-1 partial_i (root-g s^i)
        double out = 0;

        // -(3/2)/chi * n^i d_i chi = -(3/2)/(chi*L) * gamma^{ij} s_j d_i chi
        // partial_z chi vanishes, so dchi[2] = 0
        const double dchi[3] = {m_dx_chi[ii], m_dy_chi[ii], 0.0};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                out += -(3.0 / (2.0 * chi * L)) * ginv[i][j] * sh[j] * dchi[i];

        // now the remainder is 
        // partial_i s_unit_vec[i]
        // s-unit  is s_j h^ij chi / sqrt( s_m s_n h^mn chi)
        //    =    -> s_j g^ij / sqrt( s_m s_n g^mn)
          
        // d_i g^jk terms
        for (size_t i = 0; i < 3; i++)
        {
            for (size_t j = 0; j < 3; j++)
            {
                // metric derivs
                out += dgi[i][i][j]*sh[j]/L;

                // s derivs 
                out += ginv[i][j]*d_ij_F[i][j]/L;
                
                for (size_t k = 0; k < 3; k++)
                {
                    for (size_t l = 0; l < 3; l++)
                    {
                        // metric derivs 
                        out += - ginv[i][j] * dgi[i][k][l]
                                            * sh[j] * sh[k] * sh[l] 
                                            / (2.*L3);                 
                        // s derivs
                        out += - ginv[i][j] * ginv[k][l] * d_ij_F[i][k]
                                            * sh[j] * sh[l] / L3;
                    }
                }
            }
        }


        return out;
    }


    // Irreducible (horizon) mass: M_irr = sqrt(A / 16pi), from A = 16pi M^2 for Schwarzschild.
    double M_irr() const
    {
        return std::sqrt(Area() / (16.0 * M_PI));
    }

    // Total BH mass via Christodoulou-Smarr: M = M_irr + Q^2 / (4 M_irr).
    // Exact for RN; approximate for EMS (ignores scalar hair contribution).
    double M_total() const
    {
        const double Mirr = M_irr();
        const double Q    = Q_charge();
        return Mirr + Q * Q / (4. * Mirr);
    }

        // Electric charge on the same scale as bh_charge / MQ extraction.
    // E in the code is normalised as E_code = E_physical / sqrt(8pi),
    // so the Gauss flux integral (1/4pi) * integral gives Q_physical / sqrt(8pi).
    // Multiplying by sqrt(8pi) / (4pi) = 1/sqrt(2pi) recovers the physical charge
    // that reaches 1 at the extremal RN limit for M_ADM = 1.
    double Q_charge() const
    {
        double sum = 0.;
        RH_FOR_INTERIOR(ii)
        {
            const double phi = m_phi[ii];
            const double f   = m_f0 + m_f1 * phi + m_f2 * phi * phi;
            const double coupling = std::exp(-2. * m_alpha * f);
            const auto sv = s_unit_vec(ii);
            sum += coupling * (m_Ex[ii] * sv[0] + m_Ey[ii] * sv[1]) * dA(ii);
        }
        return sum / std::sqrt(2. * M_PI);
    }

    // ADM linear momentum in x: (1/16pi) * integral of (K_xj - K*gamma_xj) s^j dA
    // Using BSSN vars: integrand = [(A11 - 2K/3 h11)*sx + (A12 - 2K/3 h12)*sy] / chi
    // s^j is the curved-space unit normal; G=1 assumed.
    double PX() const
    {
        double sum = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            const auto   sc = s_unit_vec(ii);
            const double sx = sc[0], sy = sc[1];
            const double integrand =
                ( (m_A11[ii] - (2.0 / 3.0) * m_K[ii] * m_h11[ii]) * sx
                + (m_A12[ii] - (2.0 / 3.0) * m_K[ii] * m_h12[ii]) * sy )
                / m_chi[ii];
            sum += integrand * dA(ii);
        }
        return sum / (8.0 * M_PI);
    }

    // Shift the x-origin to the area-weighted centroid of the surface.
    void re_centre()
    {
        double sum_x_dA = 0.0, sum_dA = 0.0;
        RH_FOR_INTERIOR(ii)
        {
            double x_ii = m_centre[0] + m_f[ii] * std::cos(m_theta[ii]);
            sum_x_dA += x_ii * dA(ii);
            sum_dA   += dA(ii);
        }
        m_centre[0] = sum_x_dA / sum_dA;
    }

    const char *state_str() const
    {
        if (m_dead) return "dead  ";
        switch (m_state)
        {
            case SolverState::FOUND:   return "found ";
            case SolverState::CLOSE:   return "close ";
            case SolverState::MEDIUM:  return "medium";
            case SolverState::DORMANT: return "dormnt";
            default:                   return "far   ";
        }
    }

    void print_diagnostics(std::ostream &os, double a_time) const
    {
        os << std::scientific << std::setprecision(8)
           << std::setw(16) << a_time
           << std::setw(16) << m_index
           << std::setw(16) << m_centre[0]
           << std::setw(16) << average_f()
           << std::setw(16) << Area()
           << std::setw(16) << M_irr()
           << std::setw(16) << PX()
           << std::setw(16) << average_Theta_plus()
           << std::setw(16) << average_Theta_minus()
           << std::setw(16) << expansion_error()
           << std::setw(16) << K0()
           << std::setw(16) << K2()
           << std::setw(16) << K4()
           << std::setw(16) << KP2()
           << std::setw(16) << KP4()
           << std::setw(16) << equator_perim()
           << std::setw(16) << polar_perim()
           << std::setw(16) << Q_charge()
           << std::setw(16) << M_total()
           << std::setw(8)  << state_str()
           << std::endl;
    }

    // Write one row of f values: t  f[0]  f[1]  ...  f[n-1]  (interior points only)
    void print_f(std::ostream &os, double a_time) const
    {
        os << std::scientific << std::setprecision(8) << std::setw(16) << a_time;
        for (int i = 0; i < m_n; ++i)
            os << std::setw(16) << m_f[m_NG + i];
        os << std::endl;
    }

    // Remove all rows with t > a_restart_time from this surface's output files.
    // Called before re-opening in append mode on a Chombo restart.
    void prune(double a_restart_time) const
    {
        auto prune_file = [&](const std::string &fname)
        {
            std::ifstream in(fname);
            if (!in.is_open()) return; // file doesn't exist yet
            std::vector<std::string> kept;
            std::string line;
            while (std::getline(in, line))
            {
                if (line.empty() || line.front() == '#') { kept.push_back(line); continue; }
                std::istringstream iss(line);
                double t;
                if (!(iss >> t) || t <= a_restart_time) kept.push_back(line);
            }
            in.close();
            std::ofstream out(fname, std::ios::out | std::ios::trunc);
            for (const auto &l : kept)
                out << l << '\n';
        };
        prune_file("rh_surf_" + std::to_string(m_index) + ".dat");
        prune_file("rh_f"     + std::to_string(m_index) + ".dat");
    }

    void hello() const
    {
        pout() << std::fixed << std::setprecision(6)
               << "  surf = "         << m_index
               << "  x = "            << std::setw(10) << m_centre[0]
               << "  r = "            << std::setw(10) << average_f()
               << "  n = "            << std::setw(6)  << m_n
               << "  level = "        << std::setw(4)  << m_level
               << "  chase_speed = "  << std::setw(8)  << m_chase_speed
               << std::endl;
    }

    void fancy_hello() const
    {
        const std::string bar(52, '-');
        const char *mode = state_str();
        pout() << "\n  +" << bar << "+\n"
               << "  |  Horizon " << m_index
               << "   @  x = " << std::fixed << std::setprecision(4) << m_centre[0]
               << "  n= " << m_n
               << std::string(52 - 35 - std::to_string(m_index).size()-2, ' ') << "|\n"
               << "  +" << bar << "+\n"
               << std::setprecision(6)
               << "  |  <r>     = " << std::setw(12) << average_f()
               << "     Area    = " << std::setw(12) << Area()      << "  |\n"
               << "  |  M_irr   = " << std::setw(12) << M_irr()
               << "     P_x     = " << std::setw(12) << PX()        << "  |\n"
               << "  |  <Th+>   = " << std::setw(12) << average_Theta_plus()
               << "     e       = " << std::setw(12) << 1.0 - polar_perim() / equator_perim() << "  |\n"
               << "  |  err     = " << std::setw(12) << expansion_error()
               << "     mode    =       " << mode               << "  |\n"
               << "  |  K0      = " << std::setw(12) << K0()
               << "     KP2     = " << std::setw(12) << KP2()   << "  |\n"
               << "  |  Q       = " << std::setw(12) << Q_charge()
               << "     M_tot   = " << std::setw(12) << M_total()  << "  |\n"
               << "  +" << bar << "+\n";
    }

    // Move the surface one step toward Theta = 0.
    // Theta = D_i s^i + K_ij s^i s^j - K is computed at each interior point.
    // f is updated by -a_speed * Theta, then ghosts are refreshed.
    void compute_expansion()
    {
        RH_FOR_INTERIOR(ii)
        {
            m_Expansion[ii] = Theta_plus(ii);
        }
    }

    void chase_step(double a_courant = 0.25)
    {
        compute_expansion();
        RH_FOR_INTERIOR(ii)
        {
            if (!std::isfinite(m_Expansion[ii]))
                throw std::runtime_error("NaN/Inf in expansion at chase_step");
            const double speed = a_courant * m_d_theta * m_d_theta * m_f[ii] * m_f[ii];
            const double EP = m_Expansion[ii];
            const double sign = 0.01 * EP/sqrt(EP*EP + 0.0000000000001);
            const double size = abs(EP);
            m_f[ii] -= speed * (sign * sqrt(size) + EP);
            m_f[ii] = std::max(m_rmin, std::min(m_rmax, m_f[ii]));
        }
        fill_ghost_even(m_f);
    }

    // Uniformly shift every interior f value by a_delta and ghost-fill.
    // Only m_f changes; field arrays are unaffected (re-interpolated separately).
    void shift_f(double a_delta)
    {
        RH_FOR_INTERIOR(ii)
            m_f[ii] += a_delta;
        fill_ghost_even(m_f);
    }

    // Called by RHUnion::newton_step AFTER:
    //   1. compute_expansion() stored Theta_0 in m_Expansion
    //   2. shift_f(+a_delta_f) moved the surface outward
    //   3. fields were re-interpolated at the shifted surface (fills all ghosts)
    // Computes Theta_1 at the shifted surface, estimates dTheta/df from the
    // two-point difference, and jumps f to the Newton root from the original position.
    void newton_hop(double a_delta_f, double a_max_hop)
    {
        RH_FOR_INTERIOR(ii)
        {
            const double theta0 = m_Expansion[ii]; // Theta at f (before shift)
            const double theta1 = Theta_plus(ii);  // Theta at f + delta_f (current)

            // Newton root from original f:
            //   hop = -theta0 / ((theta1-theta0)/delta_f)
            // Current m_f = f + delta_f, so apply as:
            //   m_f <- m_f - delta_f + hop
            const double dT = theta1 - theta0;
            double hop;
            if (std::abs(dT) > 1.0e-12)
                hop = -theta0 * a_delta_f / dT;
            else
                hop = 0.0; // zero gradient: undo the probe shift, don't move

            // clip to Courant-like max step so a bad gradient can't explode f
            hop = std::max(-a_max_hop, std::min(a_max_hop, hop));

            m_f[ii] = m_f[ii] - a_delta_f + hop;
            m_f[ii] = std::max(m_rmin, std::min(m_rmax, m_f[ii]));
        }
        // m_f is even; field-array ghosts were already set by fill_all_ghosts()
        // inside the preceding interpolate_fields() call
        fill_ghost_even(m_f);
    }

    void fill_all_ghosts()
    {
        for (auto *v : even_arrays())
        {
            fill_ghost_even(*v);
        }
        for (auto *v : odd_arrays())
        {
            fill_ghost_odd(*v);
        }
    }

  private:
    // Even fields: Q(x, -y) = +Q(x, y).
    // Includes f, scalars, diagonal metric/K components, and x-derivatives of
    // even fields / y-derivatives of odd fields.
    std::vector<std::vector<double> *> even_arrays()
    {
        return {&m_f,
                &m_chi,    &m_K,
                &m_A11,    &m_A22,    &m_Aww,
                &m_h11,    &m_h22,    &m_hww,
                &m_dx_chi,
                &m_dx_h11, &m_dy_h12,
                &m_dx_h22, &m_dx_hww,
                &m_phi,    &m_Ex};     // phi: scalar (even); Ex: even at y=0 (parity ODD_X)
    }

    // Odd fields: Q(x, -y) = -Q(x, y).
    // Off-diagonal (xy) metric/K components and y-derivatives of even fields /
    // x-derivatives of odd fields.
    std::vector<std::vector<double> *> odd_arrays()
    {
        return {&m_h12,    &m_A12,
                &m_dy_chi,
                &m_dy_h11, &m_dx_h12,
                &m_dy_h22, &m_dy_hww,
                &m_Ey};               // Ey: odd at y=0 (parity ODD_Y)
    }

    // Even (zero-gradient) ghost fill: ghost = +nearest interior.
    void fill_ghost_even(std::vector<double> &v)
    {
        for (int k = 0; k < m_NG; ++k)
        {
            v[m_NG - 1 - k]   = v[m_NG + k];
            v[m_NG + m_n + k] = v[m_NG + m_n - 1 - k];
        }
    }

    // Odd (antisymmetric) ghost fill: ghost = -nearest interior.
    void fill_ghost_odd(std::vector<double> &v)
    {
        for (int k = 0; k < m_NG; ++k)
        {
            v[m_NG - 1 - k]   = -v[m_NG + k];
            v[m_NG + m_n + k] = -v[m_NG + m_n - 1 - k];
        }
    }
};

#endif /* RHSURF_HPP_ */
