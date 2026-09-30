#ifndef T7RHS_HPP_
#define T7RHS_HPP_
#include "CCZ4Cartoon.hpp"
#include "ExperimentalGauge.hpp"
#include "EMSCouplingFunction.hpp"

// Read-only replay of the production kernel; never stores into evolution data.
class T7RHS : public CCZ4Cartoon<ExperimentalGauge, FourthOrderDerivatives, CouplingFunction>
{
    using Base = CCZ4Cartoon<ExperimentalGauge, FourthOrderDerivatives, CouplingFunction>;
  public:
    using Base::Base;
    static constexpr int components = 3*NUM_VARS+4;
    void compute(Cell<double> cell) const
    {
        auto u = cell.load_vars<Vars>();
        auto d1 = this->m_deriv.diff1<Vars>(cell);
        auto d2 = this->m_deriv.diff2<Diff2Vars>(cell);
        auto adv = this->m_deriv.advection<Vars>(cell, u.shift);
        Vars<double> physical, ko, clean;
        this->rhs_equation(physical, u, d1, d2, adv, Coordinates<double>(cell,this->m_deriv.m_dx).y);
        VarsTools::assign(ko,0.); VarsTools::assign(clean,0.);
        this->m_deriv.add_dissipation(ko,cell,this->m_sigma);
        clean.Ex=-u.lapse*d1.Xi[0]; clean.Ey=-u.lapse*d1.Xi[1];
        clean.Bx=u.lapse*d1.Lambda[0]; clean.By=u.lapse*d1.Lambda[1];
        physical.enum_mapping([&](const int &c,double &v){cell.store_vars(v,c);});
        ko.enum_mapping([&](const int &c,double &v){cell.store_vars(v,NUM_VARS+c);});
        clean.enum_mapping([&](const int &c,double &v){cell.store_vars(v,2*NUM_VARS+c);});
        cell.store_vars((adv.Xi-physical.Xi)/u.lapse-u.Xi,3*NUM_VARS);
        cell.store_vars((physical.Lambda-adv.Lambda)/u.lapse+u.Lambda,3*NUM_VARS+1);
        cell.store_vars(adv.Xi,3*NUM_VARS+2); cell.store_vars(adv.Lambda,3*NUM_VARS+3);
    }
};
#endif
