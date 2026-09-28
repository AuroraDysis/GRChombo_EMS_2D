/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMSBH2DLEVEL_HPP_
#define EMSBH2DLEVEL_HPP_

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include <memory>
// Problem specific includes
#include "EMSCouplingFunction.hpp"
#include "EMSRadiationParameters.hpp"
// #include "EinsteinMaxwellDilatonField.hpp" //shouldnt need

class EMSBH2DLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<EMSBH2DLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);
    EMSRadiationParameters m_radiation = [this] {
        GRParmParse pp;
        return EMSRadiationParameters(pp, m_p);
    }();

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    virtual void postRegrid(int a_base_level) override;

#ifdef CH_USE_HDF5
    virtual void postRestart() override;
#endif

    /// Any actions that should happen just before checkpointing
    virtual void prePlotLevel() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    void eval_reference_rhs(GRLevelData &a_soln, GRLevelData &a_rhs,
                            const CouplingFunction &coupling_function);
    void rebuild_reference();
    void impose_reference_shift(GRLevelData &state);
    void check_reference_floors(const GRLevelData &state) const;
    void write_reference_diagnostics() const;

    std::unique_ptr<LevelData<FArrayBox>> m_reference;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Identify and tag the cells that need higher resolution
    virtual void
    computeTaggingCriterion(FArrayBox &tagging_criterion,
                            const FArrayBox &current_state) override;

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;

  public:
    // Prepare every level before the initial extraction after interpolator setup.
    void ems_prepare_radiation();
    void ems_extract_radiation();
};

#endif /* EMSBHLEVEL_HPP_ */
