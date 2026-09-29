/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4CARTOON_HPP_)
#error "This file should only be included through CCZ4Cartoon.hpp"
#endif

#ifndef CCZ4CARTOON_IMPL_HPP_
#define CCZ4CARTOON_IMPL_HPP_

#define COVARIANTZ4
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"
#include <iostream>
#include <stdio.h>

template <class gauge_t, class deriv_t, class coupling_t>
inline CCZ4Cartoon<gauge_t, deriv_t, coupling_t>::CCZ4Cartoon(
        params_t params, double dx, double sigma,
        coupling_t a_coupling, double a_G_Newton, int formulation,
        double cosmological_constant)
        : CCZ4RHS<gauge_t, deriv_t>(params, dx, sigma, formulation,
                                    cosmological_constant),
          m_coupling(a_coupling), m_G_Newton(a_G_Newton)
{
}

template <class data_t, template <typename> class vars_t>
emtensorCartoon_t<data_t>
compute_EMS_EM_tensor(const vars_t<data_t> &vars,
                     const vars_t<Tensor<1, data_t>> &d1,
const Tensor<2, data_t> &h_UU, const data_t &h_UU_ww,
const chris_t<data_t> &chris, const int &nS, const CouplingFunction &a_coupling)
{
emtensorCartoon_t<data_t> out;

// set coupling functions
data_t coupling = 0., f_coupling = 0., fprime = 0.;
a_coupling.compute_coupling(f_coupling, fprime, coupling, vars);

// 3D B and E fields
Tensor<1, data_t, 3> E, B;
E[0] = vars.Ex;
E[1] = vars.Ey;
E[2] = vars.Ez;
B[0] = vars.Bx;
B[1] = vars.By;
B[2] = vars.Bz;

// epsilon tensor
Tensor<3, data_t, 3> epsLLL = {0.};
data_t eps012 = pow(vars.chi, -1.5);
epsLLL[0][1][2] = eps012;
epsLLL[0][2][1] = -eps012;
epsLLL[1][2][0] = eps012;
epsLLL[1][0][2] = -eps012;
epsLLL[2][0][1] = eps012;
epsLLL[2][1][0] = -eps012;

// 3d spatial inverse metric
Tensor<2, data_t, 3> gammaUU_3d = {0.};
FOR2(i,j)
{
    gammaUU_3d[i][j] = h_UU[i][j] * vars.chi;
}
gammaUU_3d[2][2] = h_UU_ww * vars.chi;

// eps_{ijk} E^j B^k
Tensor<1, data_t, 3> epsEB = {0.};
for (int i=0; i<3; i++) {
for (int j=0; j<3; j++) {
for (int k=0; k<3; k++) {
for (int m=0; m<3; m++) {
for (int n=0; n<3; n++) {
    epsEB[i] += epsLLL[i][j][k] * gammaUU_3d[j][m] * E[m]
                                * gammaUU_3d[k][n] * B[n];
}}}}}

// some bits we need
data_t Dphi_sqr = 0.;
data_t BB = 0.;
data_t EE = 0.;
FOR2(i,j)
{
    Dphi_sqr += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    BB += vars.chi * h_UU[i][j] * B[i] * B[j];
    EE += vars.chi * h_UU[i][j] * E[i] * E[j];
}
// cartoon terms
BB += vars.chi * h_UU_ww * B[2] * B[2];
EE += vars.chi * h_UU_ww * E[2] * E[2];

// stress tensor components
out.rho = vars.Pi * vars.Pi + Dphi_sqr + coupling * (BB + EE);

FOR(i)
{
    out.Si[i] = 2. * (vars.Pi * d1.phi[i] + coupling * epsEB[i]);
}

FOR(i, j)
{
    out.Sij[i][j] = 2. * d1.phi[i] * d1.phi[j]
                       - vars.h[i][j] * (Dphi_sqr - vars.Pi * vars.Pi)/vars.chi
                    + coupling * (-2. * (E[i] * E[j] + B[i] * B[j])
                                      + vars.h[i][j] * (BB + EE)/vars.chi );
}



// Sww
out.Sww = - (vars.hww/vars.chi) * (Dphi_sqr - vars.Pi * vars.Pi)
                + coupling * (-2. * (vars.Ez*vars.Ez + vars.Bz*vars.Bz)
                                  + (vars.hww/vars.chi) * (BB + EE) );

// S
out.S = TensorAlgebra::compute_trace(out.Sij, h_UU);
//cartoon term
out.S += h_UU_ww * out.Sww;
// mult by chi as gamma_UU = chi * h_UU, used in trace
out.S *= vars.chi;

return out;
}

template <class gauge_t, class deriv_t, class coupling_t>
template <class data_t>
void CCZ4Cartoon<gauge_t, deriv_t, coupling_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
            this->m_deriv.template advection<Vars>(current_cell, vars.shift);

    Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);

    Vars<data_t> rhs;
    CCZ4Cartoon<gauge_t, deriv_t>::rhs_equation(rhs, vars, d1, d2, advec,
                                                coords.y);

    this->m_deriv.add_dissipation(rhs, current_cell, this->m_sigma);

    current_cell.store_vars(rhs); // Write the rhs into the output FArrayBox
}

template <class gauge_t, class deriv_t, class coupling_t>
template <class data_t, template <typename> class vars_t,
        template <typename> class diff2_vars_t>
void CCZ4Cartoon<gauge_t, deriv_t, coupling_t>::rhs_equation(
        vars_t<data_t> &rhs, const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1,
const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
const double &cartoon_coord) const
{
using namespace TensorAlgebra;
const int dI = CH_SPACEDIM - 1;
const int nS =
        GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
const double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
const double two_over_gr_spacedim = 2. * one_over_gr_spacedim;

auto h_UU = compute_inverse_sym(vars.h);
auto h_UU_ww = 1. / vars.hww;
auto chris = compute_christoffel(d1.h, h_UU);
const double one_over_cartoon_coord = 1. / cartoon_coord;
const double one_over_cartoon_coord2 = 1. / (cartoon_coord * cartoon_coord);

Tensor<1, data_t> chris_ww;
FOR(i)
        {
                chris_ww[i] =
                        one_over_cartoon_coord * (delta(i, dI) - h_UU[i][dI] * vars.hww);
        FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
        }
Tensor<1, data_t> chris_contracted; //!< includes the higher D contributions
FOR(i)
        {
                chris_contracted[i] = chris.contracted[i] + nS * h_UU_ww * chris_ww[i];
        }

Tensor<1, data_t> Z_over_chi;
Tensor<1, data_t> Z;
if (this->m_formulation == CCZ4::USE_BSSN)
{
FOR(i) Z_over_chi[i] = 0.0;
}
else
{
FOR(i) Z_over_chi[i] = 0.5 * (vars.Gamma[i] - chris_contracted[i]);
}
FOR(i) Z[i] = vars.chi * Z_over_chi[i];

auto ricci = CCZ4CartoonGeometry::compute_ricci_Z(
        vars, d1, d2, h_UU, h_UU_ww, chris, Z_over_chi, cartoon_coord);

data_t divshift = compute_trace(d1.shift) +
                  nS * vars.shift[dI] /
                  cartoon_coord; //!< includes the higher D contribution
data_t divshift_w = compute_trace(d1.shift) -
                    CH_SPACEDIM * vars.shift[dI] /
                    cartoon_coord; //!< combination that appears in the
//!< rhs of the ww components
data_t Z_dot_d1lapse = compute_dot_product(Z, d1.lapse);
data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);
data_t dlapse_dot_dhww = compute_dot_product(d1.lapse, d1.hww, h_UU);

data_t hUU_dlapse_cartoon = 0;
FOR(i)
        {
                hUU_dlapse_cartoon +=
                        one_over_cartoon_coord * h_UU[dI][i] * d1.lapse[i];
        }

Tensor<1, data_t> reg_03, reg_04;
FOR(i)
        {
                reg_03[i] = one_over_cartoon_coord * d1.shift[i][dI] -
                            one_over_cartoon_coord2 * delta(i, dI) * vars.shift[dI];
        reg_04[i] = one_over_cartoon_coord * d1.shift[dI][i] -
        one_over_cartoon_coord2 * delta(i, dI) * vars.shift[dI];
        }

Tensor<2, data_t> covdtilde2lapse;
Tensor<2, data_t> covd2lapse; // NOTE: we compute chi * D_i D_j lapse
FOR(k, l)
{
covdtilde2lapse[k][l] = d2.lapse[k][l];
FOR(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
covd2lapse[k][l] =
vars.chi * covdtilde2lapse[k][l] +
0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
vars.h[k][l] * dlapse_dot_dchi);
}
data_t covd2lapse_ww =
        vars.chi * (0.5 * dlapse_dot_dhww + vars.hww * hUU_dlapse_cartoon) -
        0.5 * vars.hww * dlapse_dot_dchi;

/*data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
FOR(i)
{
    tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
    FOR(j)
    {
        tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                       d1.lapse[i] * d1.chi[j]);
    }
}*/
data_t tr_covd2lapse = compute_trace(covd2lapse, h_UU);
tr_covd2lapse += nS * h_UU_ww * covd2lapse_ww;

// Make A_{ij} trace free
/*data_t trA = nS * h_UU_ww * vars.Aww;
FOR(i,j)
{
    trA += h_UU[i][j] * vars.A[i][j];
}*/
data_t trA = compute_trace(vars.A, h_UU);
trA += nS * h_UU_ww * vars.Aww;

Tensor<2, data_t> A_TF;
FOR(i, j)
{
A_TF[i][j] = vars.A[i][j] - trA * vars.h[i][j] * one_over_gr_spacedim;
}
data_t Aww_TF;
Aww_TF = vars.Aww - trA * vars.hww * one_over_gr_spacedim;

Tensor<2, data_t> A_UU = raise_all(A_TF, h_UU);
data_t A_UU_ww = h_UU_ww * h_UU_ww * Aww_TF;
// A^{ij} A_{ij}. - Note the abuse of the compute trace function.
data_t tr_A2 = compute_trace(A_TF, A_UU) + nS * A_UU_ww * Aww_TF;

rhs.chi = advec.chi + two_over_gr_spacedim * vars.chi *
                      (vars.lapse * vars.K - divshift);

FOR(i, j)
{
rhs.h[i][j] = advec.h[i][j] - 2.0 * vars.lapse * A_TF[i][j] -
two_over_gr_spacedim * vars.h[i][j] * divshift;
FOR(k)
        {
                rhs.h[i][j] +=
                        vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
        }
}

rhs.hww = advec.hww - 2.0 * vars.lapse * Aww_TF -
          two_over_gr_spacedim * vars.hww * divshift_w;

Tensor<2, data_t> Adot_TF;
FOR(i, j)
{
Adot_TF[i][j] =
-covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
}
data_t Adot_TF_ww = -covd2lapse_ww + vars.chi * vars.lapse * ricci.LLww;

data_t trAdot = compute_trace(Adot_TF, h_UU) + nS * h_UU_ww * Adot_TF_ww;
// make_trace_free(Adot_TF, vars.h, h_UU);
FOR(i, j) { Adot_TF[i][j] -= one_over_gr_spacedim * vars.h[i][j] * trAdot; }
Adot_TF_ww -= one_over_gr_spacedim * vars.hww * trAdot;

FOR(i, j)
{
rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] +
vars.A[i][j] * (vars.lapse * (vars.K - 2 * vars.Theta) -
two_over_gr_spacedim * divshift);
FOR(k)
        {
                rhs.A[i][j] +=
                        A_TF[k][i] * d1.shift[k][j] + A_TF[k][j] * d1.shift[k][i];
        FOR(l)
        {
            rhs.A[i][j] -=
                    2 * vars.lapse * h_UU[k][l] * A_TF[i][k] * A_TF[l][j];
        }
        }
}

rhs.Aww = advec.Aww + Adot_TF_ww +
          Aww_TF * (vars.lapse *
                    ((vars.K - 2 * vars.Theta) - 2 * h_UU_ww * Aww_TF) -
                    two_over_gr_spacedim * divshift_w);

// DANGER, ROBINS HACKED CCZ4!
data_t f_of_alpha = vars.lapse / (0.005 + vars.lapse);
#ifdef COVARIANTZ4
data_t kappa1_lapse = this->m_params.kappa1 * f_of_alpha;
#else
data_t kappa1_lapse = this->m_params.kappa1 * vars.lapse;
#endif

rhs.Theta =
advec.Theta +
0.5 * vars.lapse *
(ricci.scalar - tr_A2 +
 ((GR_SPACEDIM - 1.0) * one_over_gr_spacedim) * vars.K * vars.K -
 2 * vars.Theta * vars.K) -
0.5 * vars.Theta * kappa1_lapse *
((GR_SPACEDIM + 1) + this->m_params.kappa2 * (GR_SPACEDIM - 1)) -
Z_dot_d1lapse;

rhs.Theta += -vars.lapse * this->m_cosmological_constant;

rhs.K =
advec.K +
vars.lapse * (ricci.scalar + vars.K * (vars.K - 2 * vars.Theta)) -
kappa1_lapse * GR_SPACEDIM * (1 + this->m_params.kappa2) * vars.Theta -
tr_covd2lapse;
rhs.K += -2 * vars.lapse * GR_SPACEDIM / ((double)GR_SPACEDIM - 1.) *
this->m_cosmological_constant;

Tensor<1, data_t> Gammadot;
FOR(i)
        {
                Gammadot[i] =
                        two_over_gr_spacedim *
                        (divshift * (chris_contracted[i] +
                                     2 * this->m_params.kappa3 * Z_over_chi[i]) -
                         2 * vars.lapse * vars.K * Z_over_chi[i]) -
                        2 * kappa1_lapse * Z_over_chi[i] +
                        nS *
                        (2. * vars.lapse * chris_ww[i] * A_UU_ww + h_UU_ww * reg_03[i]);
        FOR(j)
        {
            Gammadot[i] +=
                    2 * h_UU[i][j] *
                    (vars.lapse * d1.Theta[j] - vars.Theta * d1.lapse[j]) -
                    2 * A_UU[i][j] * d1.lapse[j] -
                    vars.lapse * (two_over_gr_spacedim * (GR_SPACEDIM - 1.0) *
                                  h_UU[i][j] * d1.K[j] +
                                  GR_SPACEDIM * A_UU[i][j] * d1.chi[j] / vars.chi) -
                    (chris_contracted[j] +
                     2 * this->m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift[i][j] +
                    (nS * (GR_SPACEDIM - 2.0) * one_over_gr_spacedim) * h_UU[i][j] *
                    reg_04[j];
            FOR(k)
            {
                Gammadot[i] +=
                        2 * vars.lapse * chris.ULL[i][j][k] * A_UU[j][k] +
                        h_UU[j][k] * d2.shift[i][j][k] +
                        ((GR_SPACEDIM - 2.0) * one_over_gr_spacedim) * h_UU[i][j] *
                        d2.shift[k][j][k];
            }
        }
        rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
        }

this->m_gauge.rhs_gauge(rhs, vars, d1, d2, advec);

////////////////////////
// Stress Tensor Calc
///////////////////////

// calculate coupling
data_t coupling =0., f_coupling=0., fprime=0.;

m_coupling.compute_coupling(f_coupling, fprime, coupling, vars);

emtensorCartoon_t<data_t> emtensor =
        compute_EMS_EM_tensor(vars, d1, h_UU, h_UU_ww, chris, nS, m_coupling);

Tensor<2, data_t> Sij_TF;
data_t Sww_TF;

FOR(i, j)
{
    Sij_TF[i][j] = emtensor.Sij[i][j]
                  - one_over_gr_spacedim * emtensor.S * vars.h[i][j] / vars.chi;
}
Sww_TF = emtensor.Sww - one_over_gr_spacedim * emtensor.S * vars.hww / vars.chi;

FOR(i, j)
{
rhs.A[i][j] +=
-8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sij_TF[i][j];
}

rhs.Aww += -8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sww_TF;

rhs.K +=
4.0 * M_PI * m_G_Newton * vars.lapse * (emtensor.S - 3 * emtensor.rho);
rhs.Theta += -8.0 * M_PI * m_G_Newton * vars.lapse * emtensor.rho;

FOR(i)
{
    data_t matter_term_Gamma = 0.0;
    FOR(j)
    {
        matter_term_Gamma += -16.0 * M_PI * m_G_Newton * vars.lapse *
                             h_UU[i][j] * emtensor.Si[j];
    }

    rhs.Gamma[i] += matter_term_Gamma;

    // TERM needed for MPG, but not IMPG
    // rhs.B[i] += matter_term_Gamma;
}

////////////////////////////////////////
// Matter evolution  - Robin's code below (also stress tensor is mine)
///////////////////////////////////////

// some conveniences and damping parameters
//data_t ooy = sqrt(1./(pow(1./one_over_cartoon_coord,2)+10e-20)); // 1/y
data_t ooy = one_over_cartoon_coord; // 1/y
const int iy = dI; // index of y coord
const int ix = 0;
const int iz = 2;

///////////////////////
// // some 3d variables
///////////////////

// fields and simple derivs
Tensor<1, data_t, 3> E;
Tensor<1, data_t, 3> B;
Tensor<1, data_t, 3> Dphi = {0.};
Tensor<1, data_t, 3> Dchi = {0.};
Tensor<1, data_t, 3> Dlapse = {0.};
// E_i
E[ix] = vars.Ex;
E[iy] = vars.Ey;
E[iz] = vars.Ez;
// B_i
B[ix] = vars.Bx;
B[iy] = vars.By;
B[iz] = vars.Bz;

FOR(i) {
    Dphi[i] = d1.phi[i];
    Dchi[i] = d1.chi[i];
    Dlapse[i] = d1.lapse[i];
}
Dphi[2] = 0.;
Dchi[2] = 0.;
Dlapse[2] = 0.;

// h^ij 3d
Tensor<2, data_t, 3> h_UU_3d = {0.};
FOR2(i,j)
{
    h_UU_3d[i][j] = h_UU[i][j];
}
h_UU_3d[2][2] = h_UU_ww;



////

// physical metric gamma_ij
Tensor<2, data_t, 3> gammaLL_3d = {0.};
FOR2(i,j)
{
    gammaLL_3d[i][j] = vars.h[i][j]/vars.chi;
}
gammaLL_3d[iz][iz] = vars.hww/vars.chi;

// extrinsic curvture - K_ij
Tensor<2, data_t, 3> Kij = {0.}; // K_{ij}
data_t Kzz = (vars.Aww + vars.K*vars.hww/3.)/vars.chi;
FOR2(i,j)
{
    Kij[i][j] = (vars.A[i][j] + vars.K*vars.h[i][j]/3.)/vars.chi;
}
Kij[iz][iz] = Kzz;

// epsilon tensor eps_ijk
Tensor<3, data_t, 3> epsLLL = {0.};
data_t eps012 = pow(vars.chi, -1.5);
epsLLL[0][1][2] = eps012;
epsLLL[0][2][1] = -eps012;
epsLLL[1][2][0] = eps012;
epsLLL[1][0][2] = -eps012;
epsLLL[2][0][1] = eps012;
epsLLL[2][1][0] = -eps012;


/////////////////////////
// some derivatives required for calculations
Tensor<2, data_t, 3> dE = {0.}; // d[i] E[j] = dE[i][j]
Tensor<2, data_t, 3> dB = {0.}; // d[i] B[j] = dB[i][j]
Tensor<2, data_t, 3> dij_phi = {0.}; // d[i] B[j] = dB[i][j]
FOR(i)
{
    dE[i][0] = d1.Ex[i];
    dE[i][1] = d1.Ey[i];
    dE[i][2] = d1.Ez[i];
    dB[i][0] = d1.Bx[i];
    dB[i][1] = d1.By[i];
    dB[i][2] = d1.Bz[i];
}
//cartoon terms for z derivatives
dE[2][0] = 0.;
dE[2][1] = - ooy * vars.Ez;
dE[2][2] = ooy * vars.Ey;
dB[2][0] = 0.;
dB[2][1] = - ooy * vars.Bz;
dB[2][2] = ooy * vars.By;

FOR2(i,j)
{
    dij_phi[i][j] = d2.phi[i][j];
}
dij_phi[iz][iz] = ooy * Dphi[iy];


//////////////////////////
// derivs of various 3d metrics
//
Tensor<3, data_t, 3> dh_3d = {0.};
Tensor<3, data_t, 3> dgamma = {0.};
Tensor<3, data_t, 3> dgammainv = {0.};

//d_i h_jk
FOR3(i,j,k)
{
    dh_3d[i][j][k] = d1.h[j][k][i];
}
FOR(i)
{
    dh_3d[i][2][2] = d1.hww[i];
}
dh_3d[2][0][2] = ooy * vars.h[0][1];
dh_3d[2][2][0] = ooy * vars.h[0][1];
dh_3d[2][1][2] = ooy * (vars.h[1][1]-vars.hww);
dh_3d[2][2][1] = ooy * (vars.h[1][1]-vars.hww);

for (int i=0; i<3; i++){
for (int j=0; j<3; j++){
for (int k=0; k<3; k++){
    dgamma[i][j][k] = (dh_3d[i][j][k] - gammaLL_3d[j][k]*Dchi[i])/vars.chi;
}}}
// deriv of inverse metric
for (int i=0; i<3; i++){
for (int a=0; a<3; a++){
for (int b=0; b<3; b++){
for (int c=0; c<3; c++){
for (int d=0; d<3; d++){
    dgammainv[i][a][b] += - vars.chi * vars.chi
                              * h_UU_3d[a][c] * h_UU_3d[b][d]
                              * dgamma[i][c][d];
}}}}}

// maybe beter for well posedness ..
// d_i (h^ij chi) = - chi Gamma^j + h^ij d_i chi
// Gamma^j ~ -d_i h^ij + CCZ4 stuff
// FOR2(i,j)
// {
//     dgammainv[i][i][j] = - vars.chi * vars.Gamma[j] + h_UU[i][j] * d1.chi[i];
// }


///////////////////////////////
// connection symbols - all physical
// Tensor<3, data_t, 3> my_chris_LLL = {0.};
// Tensor<3, data_t, 3> my_chris = {0.};
// Tensor<1, data_t, 3> my_chris_contracted = {0.};
//
// // LLL christoffel symbol
// for (int i=0; i<3; i++){
// for (int j=0; j<3; j++){
// for (int k=0; k<3; k++){
//     my_chris_LLL[i][j][k] = 0.5 * ( - dgamma[i][j][k]
//                                     + dgamma[j][k][i]
//                                     + dgamma[k][i][j]);
// }}}
//
// // 3d christoffel symbol
// for (int i=0; i<3; i++){
// for (int j=0; j<3; j++){
// for (int k=0; k<3; k++){
// for (int l=0; l<3; l++){
//     my_chris[i][j][k] += h_UU_3d[i][l]*vars.chi*my_chris_LLL[l][j][k];
// }}}}
//
// // WARNING, CHECK THIS IS CORRECT FOR USE LATER BY CONVENTIONS
// // contracted physical 3d christoffel symbol
// for (int i=0; i<3; i++){
// for (int j=0; j<3; j++){
// for (int k=0; k<3; k++){
//     my_chris_contracted[i] += my_chris[i][j][k]*h_UU_3d[j][k]*vars.chi;
// }}}



/////////////////////////////
// some intermediate objects


Tensor<2, data_t, 3> the_E_term = {0.};
Tensor<2, data_t, 3> the_B_term = {0.};
for (int k = 0; k <3; k++){
for (int j = 0; j <3; j++){
    the_E_term[j][k] = B[k] * Dlapse[j] + vars.lapse * dB[j][k]
                            - 2. * vars.lapse * fprime * B[k] * Dphi[j];
    the_B_term[j][k] = - E[k] * Dlapse[j] - vars.lapse * dE[j][k];
}}

// Make some objects we need to keep the equations readable
data_t FF = 0.;
data_t EdotDphi = 0.;
data_t Dphi_dot_Dlapse = 0.;
data_t divE = 0.;
data_t divB = 0.;
data_t lap_phi = 0.; // laplacian
Tensor<1, data_t, 3> KijBj = {0.};
Tensor<1, data_t, 3> KijEj = {0.};

for (int i = 0; i <3; i++){
for (int j = 0; j <3; j++){
for (int k = 0; k <3; k++){
      KijEj[i] += E[j] * Kij[k][i] * vars.chi * h_UU_3d[j][k];
      KijBj[i] += B[j] * Kij[k][i] * vars.chi * h_UU_3d[j][k];
}}}

for (int i = 0; i <3; i++){
for (int j = 0; j <3; j++){

    EdotDphi += vars.chi * h_UU_3d[i][j] * E[i] * Dphi[j];
    Dphi_dot_Dlapse += vars.chi * h_UU_3d[i][j] * Dphi[i] * Dlapse[j];

    lap_phi += dij_phi[i][j] * vars.chi * h_UU_3d[i][j];
    lap_phi += -1.5 * h_UU_3d[i][j] * Dphi[j] * Dchi[i];
    lap_phi += Dphi[j] * dgammainv[i][i][j];

    divE += vars.chi * h_UU_3d[i][j] * dE[i][j];
    divE += -1.5 * h_UU_3d[i][j] * E[j] * Dchi[i];
    divE += E[j] * dgammainv[i][i][j];

    divB += vars.chi * h_UU_3d[i][j] * dB[i][j];
    divB += -1.5 * h_UU_3d[i][j] * B[j] * Dchi[i];
    divB += B[j] * dgammainv[i][i][j];

    FF += 2. * vars.chi * h_UU_3d[i][j] * (B[i] * B[j] - E[i] * E[j]);
}}




////////////////////////////////
// Lie Derivatives
////////////////////////////////

// Lie derivatives, uses d_i beta_j = d1.shift[j][i]
// L_beta w_i = advec w_i + w_j d_i beta^j
// e
data_t LieBetaEx = advec.Ex
                    + vars.Ex * d1.shift[0][0] + vars.Ey * d1.shift[1][0];
data_t LieBetaEy = advec.Ey
                    + vars.Ex * d1.shift[0][1] + vars.Ey * d1.shift[1][1];
data_t LieBetaEz = advec.Ez + ooy * vars.Ez * vars.shift[1];
// b
data_t LieBetaBx = advec.Bx
                    + vars.Bx * d1.shift[0][0] + vars.By * d1.shift[1][0];
data_t LieBetaBy = advec.By
                    + vars.Bx * d1.shift[0][1] + vars.By * d1.shift[1][1];
data_t LieBetaBz = advec.Bz + ooy * vars.Bz * vars.shift[iy];

////////////////////////////////
// Fill RHS of evolution variables
///////////////////////////////

/////////////////////////////////
// scalar field vars
rhs.phi = advec.phi - vars.lapse * vars.Pi;
rhs.Pi = advec.Pi + vars.lapse * (vars.K * vars.Pi - lap_phi
                                  - 0.5 * fprime * coupling * FF)
                  - Dphi_dot_Dlapse;


////////////////////////////
// electromagnetism

// electric fields - ansatz 1 of Lehner paper for Xi (more detail below)
rhs.Ex = LieBetaEx + vars.lapse * (vars.K * vars.Ex - 2. * KijEj[ix]
                                   - 2. * fprime * vars.Pi * vars.Ex
                                                        - d1.Xi[ix] );

rhs.Ey = LieBetaEy + vars.lapse * (vars.K * vars.Ey - 2. * KijEj[iy]
                                   - 2. * fprime * vars.Pi * vars.Ey
                                                        - d1.Xi[iy] );

rhs.Ez = LieBetaEz + vars.lapse * (vars.K * vars.Ez - 2. * KijEj[iz]
                                   - 2. * fprime * vars.Pi * vars.Ez);


// start filling in Magnetic fields
rhs.Bx = LieBetaBx + vars.lapse * (vars.K * vars.Bx - 2. * KijBj[ix] + d1.Lambda[ix]);
rhs.By = LieBetaBy + vars.lapse * (vars.K * vars.By - 2. * KijBj[iy] + d1.Lambda[iy]);
rhs.Bz = LieBetaBz + vars.lapse * (vars.K * vars.Bz - 2. * KijBj[iz]);

// explicit loops here as we want to loop over z cartoon terms included here
for (int j=0; j<3; j++){
for (int k=0; k<3; k++){
for (int m=0; m<3; m++){
for (int n=0; n<3; n++){

    rhs.Ex += epsLLL[ix][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_E_term[m][n];

    rhs.Ey += epsLLL[iy][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_E_term[m][n];

    rhs.Ez += epsLLL[iz][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_E_term[m][n];

    rhs.Bx += epsLLL[ix][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_B_term[m][n];

    rhs.By += epsLLL[iy][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_B_term[m][n];

    rhs.Bz += epsLLL[iz][j][k] * vars.chi * vars.chi
                              * h_UU_3d[j][m] * h_UU_3d[k][n]
                              * the_B_term[m][n];
}}}}


//////////////////////////////////////////////
// Constraint damping stuff
// https://arxiv.org/pdf/1706.09875

// decay parameters
data_t kappa_B = 1.0, kappa_E = 1.0;

// maxwell constraint damper vars

// magnetic constraint damper is Lambda
rhs.Lambda = advec.Lambda + vars.lapse * ( divB - kappa_B * vars.Lambda);

// electric constraint damper is Xi

// ansats 1 (Lehner paper)
rhs.Xi = advec.Xi - vars.lapse * (divE - 2. * EdotDphi * fprime
                                       + kappa_E * vars.Xi);


// were done, phew!

}

#endif /* CCZ4CARTOON_IMPL_HPP_ */



// old

// /////////////////////////////////
// // scalar field vars
// rhs.phi = advec.phi - vars.lapse * vars.Pi;
// rhs.Pi = advec.Pi + vars.lapse * vars.K * vars.Pi
//                   - ooy * vars.lapse * h_UU_ww * vars.chi * d1.phi[iy]
//                   -0.5 * vars.lapse * fprime * coupling * FF;
// // strictly the Gamma here absorbs the Z vector automatically.
// // chris_contracted[i] doesnt, it is -partial_i h^{ij}
// FOR (i)
// {
//     rhs.Pi += (1.-latch) * vars.chi * vars.lapse * chris_contracted[i] * d1.phi[i];
//     rhs.Pi += latch * vars.chi * vars.lapse * vars.Gamma[i] * d1.phi[i];
// }
//
// FOR2  (i,j)
// {
//     rhs.Pi += h_UU[i][j]*( - vars.chi * vars.lapse * d2.phi[i][j]
//                            + 0.5 * d1.chi[i] * vars.lapse * d1.phi[j]
//                            - vars.chi * d1.lapse[i] * d1.phi[j]);
// }
