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

int main(int argc,char **argv)
{
    if(argc!=4)return 2;
    const std::string mode=argv[1];
    if(mode=="--t0-cheb"){
        EMSTrumpetSolution_read reader; auto *old=std::cout.rdbuf(std::cerr.rdbuf()); reader.main(argv[2]);std::cout.rdbuf(old);
        std::ifstream in(argv[3],std::ios::binary);double r;
        while(in.read(reinterpret_cast<char*>(&r),8)){
            const auto q=reader.compute_radial_vars(r);double v[20]={r,q.s,reader.get_nu(),q.X,q.lapse,q.shift_R,q.K_amplitude,q.Pi};
            for(int f=0;f<3;++f)for(int d=0;d<4;++d)v[8+4*f+d]=reader.get_chebyshev_value(f,q.s,d);
            std::fwrite(v,8,20,stdout);
        }return in.eof()?0:3;
    }
    const bool analytic=mode=="--t0-analytic";
    if(!analytic && mode!="--t0-native")return 2;
    std::ifstream in(argv[2],std::ios::binary);std::ofstream out(argv[3],std::ios::binary);
    Base::params_t p{};p.kappa1=.1;p.kappa2=0;p.kappa3=1;p.covariantZ4=true;
    p.lapse_advec_coeff=1;p.shift_advec_coeff=1;p.shift_Gamma_coeff=.75;p.eta=1;
    std::vector<double> record(analytic?3+8*NUM_VARS:3+49*NUM_VARS);
    Box b(IntVect(D_DECL(-3,-3,0)),IntVect(D_DECL(3,3,0)));FArrayBox u(b,NUM_VARS),dummy(b,NUM_VARS);
    double result[Audit::columns];double last_h=-1;std::unique_ptr<Audit> kernel;
    while(in.read(reinterpret_cast<char*>(record.data()),record.size()*8)){
        const double h=analytic?1.:record[0];if(h!=last_h){kernel=std::make_unique<Audit>(p,h,1.,CouplingFunction({4*M_PI,0.,0.,-.8}),1.);last_h=h;}
        if(analytic){
            Base::Vars<double> v,adv,ko;Base::Vars<Tensor<1,double>> d1;Base::Diff2Vars<Tensor<2,double>> d2;
            v.enum_mapping([&](int c,double& x){x=record[3+c];});
            d1.enum_mapping([&](int c,Tensor<1,double>& x){FOR(j)x[j]=record[3+NUM_VARS+2*c+j];});
            d2.enum_mapping([&](int c,Tensor<2,double>& x){FOR(j,k)x[j][k]=record[3+3*NUM_VARS+4*c+2*j+k];});
            adv.enum_mapping([&](int c,double& x){x=record[3+7*NUM_VARS+c];});VarsTools::assign(ko,0.);
            kernel->evaluate(v,d1,d2,adv,ko,record[1],result);
        }else{
            const int j=std::llround(record[1]/h-.5);const IntVect c(D_DECL(0,j,0));
            u.shift(c-u.box().smallEnd()-3*IntVect::Unit);dummy.shift(c-dummy.box().smallEnd()-3*IntVect::Unit);
            std::copy(record.begin()+3,record.end(),u.dataPtr());
            // Production specificEvalRHS projects and floors the entire input, including ghosts.
            BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(),PositiveChiAndAlpha(1e-12,1e-12)),u,u,u.box(),disable_simd());
            BoxPointers pointers(u,dummy);Cell<double> cell(c,pointers);
            kernel->native(cell,record[1],result);
            double total[NUM_VARS];for(int v=0;v<NUM_VARS;++v)total[v]=dummy(c,v);
            kernel->Base::compute(cell);
            for(int v=0;v<NUM_VARS;++v){double actual=dummy(c,v);if(std::memcmp(&actual,&total[v],8)!=0)result[Audit::columns-1]+=1.;}
        }
        out.write(reinterpret_cast<char*>(result),sizeof(result));
    }
    return in.eof() && out?0:3;
}
