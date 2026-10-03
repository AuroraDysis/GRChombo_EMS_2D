#ifndef EMS_T22_AUDIT_HPP
#define EMS_T22_AUDIT_HPP
#include <cmath>
#include "EMSCouplingFunction.hpp"
#include "ExperimentalGauge.hpp"
#include "CCZ4Cartoon.hpp"
#include "BoxLoops.hpp"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>

// Tests only. Native derivatives, native complete RHS, and an independent
// native stress-tensor route to the Gamma matter source. No static reference.
template <class gauge_t>
class T22NativeAudit : public CCZ4Cartoon<gauge_t,FourthOrderDerivatives,CouplingFunction>
{
    using Base=CCZ4Cartoon<gauge_t,FourthOrderDerivatives,CouplingFunction>;
  public:
    using typename Base::params_t;
    static constexpr int columns=4*NUM_VARS+12;
    void sample(Cell<double> cell, double y, double *out, double *alias_out=nullptr) const
    {
        using namespace TensorAlgebra;
        auto u=cell.template load_vars<Base::template Vars>();
        auto d1=this->m_deriv.template diff1<Base::template Vars>(cell);
        auto d2=this->m_deriv.template diff2<Base::template Diff2Vars>(cell);
        auto adv=this->m_deriv.template advection<Base::template Vars>(cell,u.shift);
        typename Base::template Vars<double> pre,ko,total,geo;
        this->rhs_equation(pre,u,d1,d2,adv,y);
        VarsTools::assign(ko,0.);
        this->m_deriv.add_dissipation(ko,cell,this->m_sigma);
        total=pre;
        this->m_deriv.add_dissipation(total,cell,this->m_sigma);
        const auto store=[&](auto v,int offset) {
            v.enum_mapping([&](int c,double x){out[offset+c]=x;});
        };
        store(u,0);store(pre,NUM_VARS);store(ko,2*NUM_VARS);store(total,3*NUM_VARS);
        const auto inv=compute_inverse_sym(u.h);
        const auto ch=compute_christoffel(d1.h,inv);
        const auto em=compute_EMS_EM_tensor(u,d1,inv,1/u.hww,ch,1,this->m_coupling);
        T22NativeAudit<gauge_t> geometric(this->m_params,this->m_deriv_dx(),
                                         this->m_sigma,this->m_coupling,0.,this->m_formulation);
        geometric.rhs_equation(geo,u,d1,d2,adv,y);
        FOR(i)
        {
            double matter=0.;
            FOR(j) matter+=-16.*M_PI*this->m_G_Newton*u.lapse*inv[i][j]*em.Si[j];
            out[4*NUM_VARS+i]=geo.Gamma[i];
            out[4*NUM_VARS+2+i]=matter;
            out[4*NUM_VARS+4+i]=em.Si[i];
            const auto &p=this->m_params;
            out[4*NUM_VARS+6+i]=p.shift_advec_coeff*adv.B[i]
                -p.shift_advec_coeff*adv.Gamma[i]+pre.Gamma[i]-p.eta*u.B[i];
            out[4*NUM_VARS+8+i]=p.shift_advec_coeff*adv.shift[i]+p.shift_Gamma_coeff*u.B[i];
        }
        const auto &p=this->m_params;
        out[4*NUM_VARS+10]=p.lapse_advec_coeff*adv.lapse
            -p.lapse_coeff*std::pow(u.lapse,p.lapse_power)*(u.K-2*u.Theta);
        // Independent orchestration identity: production kernel, including KO.
        this->Base::compute(cell);
        long long mismatches=0;
        // Symmetric enum_mapping visits both A12 and A21 under the same
        // stored component. Native storage keeps the last visit. Compare
        // exactly those 28 stored doubles, not the overwritten tensor alias.
        if(alias_out)total.enum_mapping([&](int c,double x) {
            const auto &pointers=cell.get_box_pointers();
            double actual=pointers.m_out_ptr[c][pointers.get_out_index(cell)];
            if(std::memcmp(&actual,&x,sizeof(double))!=0)
            {alias_out[0]=c;alias_out[1]=x;alias_out[2]=actual;}
        });
        for(int c=0;c<NUM_VARS;++c)
        {
            const auto &pointers=cell.get_box_pointers();
            double actual=pointers.m_out_ptr[c][pointers.get_out_index(cell)];
            mismatches+=std::memcmp(&actual,out+3*NUM_VARS+c,sizeof(double))!=0;
        }
        out[4*NUM_VARS+11]=double(mismatches);
    }
  private:
    // FourthOrderDerivatives stores inverse dx, but the audit's constructor
    // retains its own input for the zero-G matter-isolation instance.
    double m_dx_audit;
    double m_deriv_dx() const { return m_dx_audit; }
  public:
    T22NativeAudit(params_t p,double dx,double sigma,CouplingFunction c,double g,int formulation=0)
        :Base(p,dx,sigma,c,g,formulation),m_dx_audit(dx) {}
};
#endif
