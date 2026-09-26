/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXSUPERPOSITION_METRIC_HPP_
#define FIXSUPERPOSITION_METRIC_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

class FixSuperposition_metric
{
    // Only variables needed are metric
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;
        Tensor<2, data_t> A;
        Tensor<1, data_t> shift;
        data_t hww;
        data_t Aww;
        data_t lapse;
        data_t K;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_h11, c_h22>(), h);

            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_A11, c_A22>(), A);

            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_shift1, c_shift2>(), shift);

            VarsTools::define_enum_mapping(mapping_function, c_hww, hww);

            VarsTools::define_enum_mapping(mapping_function, c_Aww, Aww);

            VarsTools::define_enum_mapping(mapping_function, c_lapse, lapse);

            VarsTools::define_enum_mapping(mapping_function, c_K, K);
        }
    };

  protected:
    const FourthOrderDerivatives m_deriv;
    const double m_dx;
    std::array<double, CH_SPACEDIM> m_center;
    const double m_rapidity;
    const bool m_binary;

  public:
    FixSuperposition_metric(
        double a_dx, std::array<double, CH_SPACEDIM> a_center, double a_rapidity, bool a_binary)
        : m_deriv(a_dx), m_dx(a_dx), m_center(a_center), m_rapidity(a_rapidity), m_binary(a_binary)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        //////////////////////////////
        // my stuff
        //////////////////////////////

        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const auto vars = current_cell.template load_vars<Vars>();

        using namespace TensorAlgebra;

        double r = sqrt(coords.x*coords.x + coords.y*coords.y);
        double ramp_r = 1.;
        double ramp_frac = 0.2;
        // if (r < 0.5 * ramp_frac)
        // {
        //     ramp_r = pow(sin(M_PI * r / ramp_frac),2);
        // }

        // INPUT H IS ACTAULLY GAMMA HERE!!!!!!!
        // fix superposition error in physical metric and lapse
        double one_if_binary = 0.;
        if (m_binary)
        {
            one_if_binary = 1.;
        }
        double lapse = vars.lapse - one_if_binary;

        // load space metric
        Tensor<2, data_t, 3> gamma = {0.};
        FOR2(i,j) gamma[i][j] = vars.h[i][j];
        FOR1(i) gamma[i][i] -= one_if_binary; // superposition fix
        gamma[2][2] = vars.hww - one_if_binary;

        // load Kij (stored in Aij previously)
        Tensor<2, data_t, 3> Kij = {0.};
        Tensor<2, data_t, 3> Aij = {0.};
        FOR2(i,j) Kij[i][j] = vars.A[i][j];
        Kij[2][2] = vars.Aww;

        // inverse physical metric
        const auto gamma_UU = compute_inverse_sym(gamma);


        // conformally decompose metric
        double det_gamma = compute_determinant_sym(gamma);
        double chi = pow(det_gamma,-1./3.);

        // calculate trace
        double K = 0.;
        for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          K += gamma_UU[i][j] * Kij[i][j];
        }}

        // definition of A_ij
        for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          Aij[i][j] = chi * (Kij[i][j] - K * gamma[i][j]/3.);
        }}


        // // assign values in the output FArrayBox

        // metric
        current_cell.store_vars(chi*gamma[2][2],c_hww);
        current_cell.store_vars(chi*gamma[0][0],c_h11);
        current_cell.store_vars(chi*gamma[0][1],c_h12);
        current_cell.store_vars(chi*gamma[1][1],c_h22);
        //
        // curvature
        current_cell.store_vars(Aij[2][2],c_Aww);
        current_cell.store_vars(Aij[0][0],c_A11);
        current_cell.store_vars(Aij[0][1],c_A12);
        current_cell.store_vars(Aij[1][1],c_A22);

        // trace K
        current_cell.store_vars(K, c_K);

        // // precollapsed lapse
        current_cell.store_vars(sqrt(chi),c_lapse);
        // // raw lapse
        // current_cell.store_vars(lapse, c_lapse);

        // conformal factor
        current_cell.store_vars(chi, c_chi);



    }
};

#endif /* FIXSUPERPOSITION_METRIC_HPP_ */
