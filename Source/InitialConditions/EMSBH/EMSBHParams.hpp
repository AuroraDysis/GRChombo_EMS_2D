/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef EMSBHPARAMS_HPP_
#define EMSBHPARAMS_HPP_

//! A structure for the input params for the emdbh
struct EMSBH_params_t
{
    int gridpoints;
    double separation;
    double bh_charge;
    double bh_mass;
    bool binary;
    bool boosted;
    double rapidity;
    double Newtons_constant;
    double Ylm_amplitude;
    double Ylm_thickness;
    double Ylm_r0;
    std::string data_path;
    std::string ctt_data_path; //!< optional EMSCTT 1 companion
    bool use_maximal_initial_lapse = false; //!< t=0 unboosted single hole only
    bool use_geometric_initial_lapse = false; //!< t=0 full boosted/binary geometry
    std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre
};

#endif /* EMSBHPARAMS_HPP_ */
