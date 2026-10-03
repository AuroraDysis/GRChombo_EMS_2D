#ifndef EMS_GAUGE_SELECTION_HPP_
#define EMS_GAUGE_SELECTION_HPP_

#include "GRParmParse.hpp"
#include "MayDay.H"
#include "SPMD.H"
#include "parstream.H"
#include <fstream>
#include <iomanip>
#include <string>

// A gauge package is part of checkpoint provenance, not a runtime switch on
// an evolved state. Untagged checkpoints belong to the historical package.
class EMSGaugeSelection
{
    std::string m_name;
  public:
    explicit EMSGaugeSelection(GRParmParse &pp)
    {
        pp.load("ems_gauge", m_name, std::string("experimental"));
        if (m_name != "experimental" && m_name != "moving_puncture")
            MayDay::Error("ems_gauge must be experimental or moving_puncture");
        bool capture = false;
        if (pp.contains("t21_rhs_capture")) pp.load("t21_rhs_capture", capture);
        if (capture)
        {
            std::string label, driver;
            pp.load("t21_gauge_label", label);
            pp.load("t21_driver_semantics", driver);
            if (label != name() || driver != driver_semantics())
                MayDay::Error("T21 recording labels do not match selected EMS gauge package");
        }
    }
    const std::string &name() const { return m_name; }
    bool differential() const { return m_name == "moving_puncture"; }
    std::string driver_semantics() const
    {
        return differential() ? "differential_gamma_driver"
                              : "integrated_offset_decaying_0p1";
    }
    void check_checkpoint(const std::string &stored) const
    {
        if (stored != m_name)
            MayDay::Error("EMS checkpoint gauge does not match ems_gauge; start fresh initial data");
    }
    template <class parameters_t>
    void record(const parameters_t &p) const
    {
        pout() << "EMS gauge package: " << m_name << "; B: "
               << driver_semantics() << std::endl;
        if (procID() != 0) return;
        std::ofstream out(p.data_path + "ems-gauge.txt");
        out << std::setprecision(17) << "ems_gauge = " << m_name
            << "\nB_semantics = " << driver_semantics()
            << "\nlapse_advec_coeff = " << p.ccz4_params.lapse_advec_coeff
            << "\nlapse_coeff = " << p.ccz4_params.lapse_coeff
            << "\nlapse_power = " << p.ccz4_params.lapse_power
            << "\nshift_advec_coeff = " << p.ccz4_params.shift_advec_coeff
            << "\nshift_Gamma_coeff = " << p.ccz4_params.shift_Gamma_coeff
            << "\neta = " << p.ccz4_params.eta
            << "\ncoefficient_semantics = "
            << (differential() ? "all six coefficients used by MovingPunctureGauge"
                              : "declared parameters; legacy ExperimentalGauge retains its own fixed lapse/advection/decay coefficients")
            << '\n';
        if (!out) MayDay::Error("cannot write EMS gauge run record");
    }
};
#endif
