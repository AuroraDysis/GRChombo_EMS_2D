// Offline t=0 only. Native kernels, saved current-field ghosts; no evolution.
#include <cmath>
#include "EMSCouplingFunction.hpp"
#include "ExperimentalGauge.hpp"
#include "CCZ4Cartoon.hpp"
#include "TraceARemovalCartoon.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "BoxLoops.hpp"
#include "EMSTrumpetSolution_read.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <cstring>

using Base = CCZ4Cartoon<ExperimentalGauge,FourthOrderDerivatives,CouplingFunction>;
class Audit : public Base
{
public:
    using Base::Base;
    static constexpr int groups=13, columns=4*NUM_VARS+2*groups+7;
    void evaluate(const Vars<double>& u,const Vars<Tensor<1,double>>& d1,
                  const Diff2Vars<Tensor<2,double>>& d2,const Vars<double>& adv,
                  const Vars<double>& ko,double y,double *out) const
    {
        Vars<double> rhs; rhs_equation(rhs,u,d1,d2,adv,y);
        auto store=[&](auto v,int start){v.enum_mapping([&](int c,double x){out[start+c]=x;});};
        store(u,0);store(rhs,NUM_VARS);store(ko,2*NUM_VARS);store(adv,3*NUM_VARS);
        using namespace TensorAlgebra;
        const auto inv=compute_inverse_sym(u.h); const double iw=1/u.hww;
        const auto ch=compute_christoffel(d1.h,inv);
        const double tr=compute_trace(u.A,inv)+u.Aww*iw;
        auto at=u.A; FOR(i,j) at[i][j]-=tr*u.h[i][j]*(1./3.);
        const auto au=raise_all(at,inv); const double aw=iw*iw*(u.Aww-tr*u.hww*(1./3.));
        const double iy=1/y,iy2=1/(y*y);
        Tensor<1,double> ww,ct,z;
        FOR(i){ww[i]=iy*((i==1?1.:0.)-inv[i][1]*u.hww);
            FOR(j) ww[i]-=.5*inv[i][j]*d1.hww[j];
            ct[i]=ch.contracted[i]+iw*ww[i]; z[i]=.5*(u.Gamma[i]-ct[i]);}
        const double div=compute_trace(d1.shift)+u.shift[1]/y;
        const double kap=m_params.kappa1*u.lapse/(.005+u.lapse);
        const auto em=compute_EMS_EM_tensor(u,d1,inv,iw,ch,1,m_coupling);
        double terms[groups][2]={};
        FOR(i){
            terms[0][i]=adv.Gamma[i];
            terms[8][i]=2*u.lapse*ww[i]*aw;
            terms[9][i]=iw*(iy*d1.shift[i][1]-iy2*(i==1?1.:0.)*u.shift[1]);
            terms[10][i]=(2./3.)*(div*(ct[i]+2*m_params.kappa3*z[i])-2*u.lapse*u.K*z[i])-2*kap*z[i];
            FOR(j){
                terms[3][i]-=2*au[i][j]*d1.lapse[j];
                terms[4][i]-=3*u.lapse*au[i][j]*d1.chi[j]/u.chi;
                terms[5][i]-=(4./3.)*u.lapse*inv[i][j]*d1.K[j];
                terms[6][i]+=2*inv[i][j]*(u.lapse*d1.Theta[j]-u.Theta*d1.lapse[j]);
                terms[9][i]+=(1./3.)*inv[i][j]*(iy*d1.shift[1][j]-iy2*(j==1?1.:0.)*u.shift[1]);
                terms[10][i]-=(ct[j]+2*m_params.kappa3*z[j])*d1.shift[i][j];
                terms[11][i]-=16*M_PI*m_G_Newton*u.lapse*inv[i][j]*em.Si[j];
                FOR(k){terms[1][i]+=inv[j][k]*d2.shift[i][j][k];
                    terms[2][i]+=inv[i][j]*d2.shift[k][j][k]/3;
                    terms[7][i]+=2*u.lapse*ch.ULL[i][j][k]*au[j][k];}
            }
            terms[12][i]=ko.Gamma[i];
        }
        double defect=0;
        FOR(i){double sum=0,scale=0;
            for(int g=0;g<groups;++g){out[4*NUM_VARS+2*g+i]=terms[g][i];if(g<12){sum+=terms[g][i];scale+=std::abs(terms[g][i]);}}
            defect=std::max(defect,std::abs(sum-rhs.Gamma[i])/(std::numeric_limits<double>::epsilon()*std::max(1.,scale)));
        }
        int start=4*NUM_VARS+2*groups;
        FOR(i){out[start+i]=m_params.shift_Gamma_coeff*u.Gamma[i]-m_params.eta*u.shift[i]-u.B[i];out[start+3+i]=-.1*u.B[i];}
        out[start+2]=-1.8*u.lapse*(u.K-2*u.Theta);out[start+5]=defect;
        out[start+6]=0.;
    }
    void native(Cell<double> cell,double y,double *out) const
    {
        auto u=cell.load_vars<Vars>();auto d1=m_deriv.diff1<Vars>(cell);auto d2=m_deriv.diff2<Diff2Vars>(cell);
        auto adv=m_deriv.advection<Vars>(cell,u.shift);Vars<double> ko;VarsTools::assign(ko,0.);
        m_deriv.add_dissipation(ko,cell,m_sigma);evaluate(u,d1,d2,adv,ko,y,out);
        Vars<double> total;rhs_equation(total,u,d1,d2,adv,y);m_deriv.add_dissipation(total,cell,m_sigma);
        cell.store_vars(total);
    }
};

