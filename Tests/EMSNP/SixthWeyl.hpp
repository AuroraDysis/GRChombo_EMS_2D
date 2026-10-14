#ifndef EMS_NP_SIXTH_WEYL_HPP
#define EMS_NP_SIXTH_WEYL_HPP
#include "EMSMatterWeylOmScalar.hpp"
#include "SixthOrderDerivatives.hpp"

// Test-only independent order estimate for the existing cell-then-interpolate W.
class EMSNPSixthWeyl : public EMSMatterWeylOmScalar
{
    bool m_matter;
  public:
    EMSNPSixthWeyl(std::array<double,CH_SPACEDIM> center,double dx,
        CouplingFunction::params_t coupling,double G,bool matter)
        : EMSMatterWeylOmScalar(center,dx,coupling,G),m_matter(matter) {}
    void compute(Cell<double> cell) const
    {
        const auto vars=cell.load_vars<Vars>();
        const SixthOrderDerivatives deriv(m_deriv.m_dx);
        const auto d1=deriv.diff1<Vars>(cell);
        const auto d2=deriv.diff2<Diff2Vars>(cell);
        const Coordinates<double> coords(cell,m_deriv.m_dx,m_center);
        auto r=compute_riemann_tensor(vars,d1,d2,coords);
        if (m_matter) add_matter(r,vars,d1);
        cell.store_vars(compute_Weyl_Om(r,vars,d1,d2,coords).Om22,c_Weyl4_Re);
        cell.store_vars(0.,c_Weyl4_Im);
    }
};
#endif
