/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMSTRUMPETSOLUTION_READ_HPP_
#define EMSTRUMPETSOLUTION_READ_HPP_

#include <array>
#include <map>
#include <string>
#include <vector>

//! Reads the EMSTRUMPET 1 Chebyshev solution and reconstructs radial data.
class EMSTrumpetSolution_read
{
  public:
    struct ems_radial_vars_t
    {
        double s;                       //! s: compactified coordinate
        double r;                       //! r: areal radius
        double delta;                   //! delta: radial metric exponent
        double phi;                     //! phi: scalar field
        double X;                       //! X: radial metric factor
        double lapse;                   //! alpha_K: trumpet lapse
        double shift_R;                 //! beta_R: radial shift
        double K_amplitude;             //! k: trace-free curvature amplitude
        double Pi;                      //! Pi: scalar momentum
        double E_R;                     //! E_R: radial electric field
        double dr_ds_numerator;         //! G = Y + (1-s) Y_s
        double compactification_factor; //! k_s = 1 + (nu-1) s
        double dP_dR;                   //! P_R, where P = 1/X
        double dlapse_dR;               //! alpha_R
        double dshift_dR;               //! beta_R derivative
    };

    //! Abort on an invalid profile, as for other initial data readers.
    void main(const std::string &a_data_path);
    //! Validate and load a profile without aborting (used by file
    //! tests).
    bool check_file(const std::string &a_data_path, std::string &a_reason);
    bool read_from_file(std::string &a_reason);
    double get_chebyshev_value(int a_field, double a_s, int a_order = 0) const;
    double get_phi_inf() const { return get_metadata_value("phi_inf"); }
    double get_native_charge() const { return get_metadata_value("q_native"); }
    double get_nu() const { return get_metadata_value("nu"); }
    std::array<double, 4> get_coupling_parameters() const;
    double get_compactified_coordinate(double a_R) const;
    ems_radial_vars_t compute_radial_vars(double a_R) const;

  private:
    using coefficients_t = std::array<std::vector<double>, 4>;
    std::array<coefficients_t, 3> m_coefficients;
    std::map<std::string, double> m_metadata_values;
    std::string m_data_path;

    static bool parse_finite_double(const std::string &a_text, double &a_value);
    static bool parse_degree(const std::string &a_text, int &a_value);
    static std::vector<double>
    compute_derivative_coefficients(const std::vector<double> &a_coefficients);
    static double evaluate_clenshaw(const std::vector<double> &a_coefficients,
                                    double a_s);
    bool invert_compactified_coordinate(double a_R, double &a_s) const;
    double get_metadata_value(const std::string &a_key) const
    {
        return m_metadata_values.at(a_key);
    }
    double compute_coupling(double a_phi) const;
};

#include "EMSTrumpetSolution_read.impl.hpp"

#endif /* EMSTRUMPETSOLUTION_READ_HPP_ */
