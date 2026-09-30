/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4CARTOON_HPP_
#define CCZ4CARTOON_HPP_

#include "CCZ4CartoonGeometry.hpp"
#include "CCZ4CartoonVars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"
//#include "DefaultEMDCouplingFunction.hpp"

#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

#include <array>

/// Compute class to calculate the CCZ4 right hand side with cartoon
/**
 * This compute class implements the CCZ4 cartoon right hand side equations. Use
 *it by handing it to a loop in the BoxLoops namespace. CCZ4 cartoon inherits
 *from CCZ4 and adds the cartoon to the standar CCZ4 right
 *hand side.
 **/
template <class gauge_t = MovingPunctureGauge,
          class deriv_t = FourthOrderDerivatives,
          class coupling_t = CouplingFunction>
class CCZ4Cartoon : public CCZ4RHS<gauge_t>
{
  public:
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;

    using params_t = typename CCZ4RHS<gauge_t, deriv_t>::params_t;

    /// CCZ4 cartoon variables
    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;

    /// CCZ4 cartoon variables
    template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsWithGauge<data_t>;

    /// Constructor
    CCZ4Cartoon(
        params_t params, //!< The CCZ4 cartoon parameters
        double dx,       //!< The grid spacing
        double sigma,    //!< Kreiss-Oliger dissipation coefficient
        coupling_t a_coupling, // EMS coupling function
        double a_G_Newton,
        int formulation = CCZ4::USE_CCZ4, //!< Switches between CCZ4, BSSN,...
        double cosmological_constant = 0  //!< Value of the cosmological const.

    );

    // ScalarBubble_2D::params_t m_bubble_params;

    /// Compute function
    /** This function orchestrates the calculation of the rhs for one specific
     * grid cell. This function is called by the BoxLoops::loop for each grid
     * cell; there should rarely be a need to call it directly.
     */
    template <class data_t> void compute(Cell<data_t> current_cell) const;
    void t7_capture(FArrayBox *out,const std::vector<Box> *regions)
    { m_t7_capture=out;m_t7_regions=regions; }

  protected:
    /// Calculates the rhs for CCZ4 with cartoon terms
    /** Calculates the right hand side for CCZ4 with slicing \f$- n \alpha^m (K
     *- 2\Theta)\f$ and Gamma-Driver shift condition. The variables (the
     *template argument vars_t) must contain at least the members: chi,
     *h[i][j], Gamma[i], A[i][j], Theta, lapse and shift[i].
     **/
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void rhs_equation(
        vars_t<data_t> &rhs, //!< Reference to the variables into which the
        //! output right hand side is written
        const vars_t<data_t> &vars, //!< The values of the current variables
        const vars_t<Tensor<1, data_t>>
            &d1, //!< First derivative of the variables
        const diff2_vars_t<Tensor<2, data_t>>
            &d2, //!< The second derivative the variables
        const vars_t<data_t>
            &advec, //!< The advection derivatives of the variables
        const double &cartoon_coord) const; //!< cartoon radial coordinate

    const double m_G_Newton;
    const coupling_t m_coupling;
    FArrayBox *m_t7_capture=nullptr;
    const std::vector<Box> *m_t7_regions=nullptr;
};

#include "CCZ4Cartoon.impl.hpp"

#endif /* CCZ4CARTOON_HPP_ */
