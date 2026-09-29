// Offline sensitivity audit of saved Float64 fields; no evolution or reader.
#include <cmath>
#include <fstream>
#include <limits>

// First-order absolute error budget: one epsilon per input and operation.
// This is a conservative sensitivity scale, not an interval certificate.
struct Budget
{
    double v = 0., e = 0.;
    static constexpr double eps = std::numeric_limits<double>::epsilon();
    Budget() = default;
    Budget(double value) : v(value), e(eps * std::abs(value)) {}
    Budget(double value, double error) : v(value), e(error) {}
    friend Budget operator+(Budget a, Budget b)
    { double v = a.v + b.v; return {v, a.e + b.e + eps * std::abs(v)}; }
    friend Budget operator-(Budget a, Budget b)
    { double v = a.v - b.v; return {v, a.e + b.e + eps * std::abs(v)}; }
    friend Budget operator-(Budget a) { return {-a.v, a.e}; }
    friend Budget operator*(Budget a, Budget b)
    { double v = a.v * b.v; return {v, std::abs(b.v)*a.e + std::abs(a.v)*b.e + eps*std::abs(v)}; }
    friend Budget operator/(Budget a, Budget b)
    { double v = a.v / b.v; return {v, (a.e + std::abs(v)*b.e)/std::abs(b.v) + eps*std::abs(v)}; }
    Budget &operator+=(Budget b) { return *this = *this + b; }
    Budget &operator-=(Budget b) { return *this = *this - b; }
    Budget &operator*=(Budget b) { return *this = *this * b; }
    Budget &operator/=(Budget b) { return *this = *this / b; }
    friend bool operator>(Budget a, Budget b) { return a.v > b.v; }
    friend bool operator<=(Budget a, Budget b) { return a.v <= b.v; }
};
Budget sqrt(Budget a) { double v = std::sqrt(a.v); return {v, a.e/(2*v) + Budget::eps*v}; }
Budget exp(Budget a) { double v = std::exp(a.v); return {v, v*a.e + Budget::eps*v}; }
Budget pow(Budget a, double b)
{ double v = std::pow(a.v,b); return {v, std::abs(b*v/a.v)*a.e + Budget::eps*std::abs(v)}; }
Budget simd_max(double a, Budget b) { return b.v > a ? b : Budget(a); }

#include "EMSCouplingFunction.hpp"
#include "CCZ4Cartoon.hpp"
#include "ConstraintsCartoon.hpp"
#include <iostream>

struct Audit : Constraints<CouplingFunction>
{
    Audit() : Constraints(1., CouplingFunction({4*M_PI,0.,0.,-.8}), 1.) {}
    using Constraints::constraint_equations;
};

using V = CCZ4CartoonVars::VarsNoGauge<Budget>;
using D1 = CCZ4CartoonVars::VarsNoGauge<Tensor<1,Budget>>;
using D2 = CCZ4CartoonVars::Diff2VarsNoGauge<Tensor<2,Budget>>;

// Same contractions as GammaCartoonCalculator; the initialization's derivative
// rounding must be propagated when the stored Gamma is differentiated again.
Tensor<1,Budget> gamma_budget(const V &v, const D1 &d, double y)
{
    using namespace TensorAlgebra;
    auto inv = compute_inverse_sym(v.h);
    auto c = compute_christoffel(d.h, inv);
    Tensor<1,Budget> g;
    FOR(i)
    {
        Budget cw = (delta(i,1) - inv[i][1]*v.hww)/y;
        FOR(j) cw -= .5*inv[i][j]*d.hww[j];
        g[i] = c.contracted[i] + cw/v.hww;
    }
    return g;
}

Budget gauss(const V &v, const D1 &d, double y)
{
    auto inv = TensorAlgebra::compute_inverse_sym(v.h);
    Budget E[2] = {v.Ex,v.Ey};
    Budget de[2][2] = {{d.Ex[0],d.Ey[0]},{d.Ex[1],d.Ey[1]}};
    Budget div = 0.;
    FOR(a)
    {
        Budget q=0.,dq=0.,tr=0.;
        FOR(j)
        {
            q += inv[a][j]*E[j];
            dq += inv[a][j]*de[a][j];
            FOR(k)
            {
                tr += inv[j][k]*d.h[j][k][a];
                FOR(l) dq -= inv[a][j]*d.h[j][k][a]*inv[k][l]*E[l];
            }
        }
        Budget up = v.chi*q;
        Budget lv = .5*tr + .5*d.hww[a]/v.hww - 1.5*d.chi[a]/v.chi;
        div += v.chi*dq + d.chi[a]*q + up*lv;
    }
    Budget uy=0.,pg=0.;
    FOR(j) uy += v.chi*inv[1][j]*E[j];
    FOR(a,j) pg += v.chi*inv[a][j]*E[j]*d.phi[a];
    return exp(-2.*(4*M_PI)*(-.8*v.phi*v.phi)) *
        (div + uy/y - 2.*(4*M_PI)*(-1.6*v.phi)*pg);
}

int main(int argc, char **argv)
{
    if (argc != 3) return 2;
    std::ifstream in(argv[1],std::ios::binary);
    std::ofstream out(argv[2],std::ios::binary);
    Audit audit;
    double record[3+NUM_VARS*25];
    size_t cells=0;
    while (in.read(reinterpret_cast<char*>(record),sizeof(record)))
    {
        double h=record[0],y=record[1],xc=record[2];
        auto value=[&](int n,int x,int z) { return record[3+n*25+(z+2)*5+x+2]; };
        double initialization_error[NUM_VARS];
        for(int n=0;n<NUM_VARS;++n)
        {
            auto derivative=[&](bool x)
            {
                auto u=[&](int k) { return x ? value(n,k,0) : value(n,0,k); };
                return ((1./12.)*u(-2)-(2./3.)*u(-1)+(2./3.)*u(1)-(1./12.)*u(2))/h;
            };
            double ux=derivative(true),uy=derivative(false);
            // Initial reader terminates its log-radius inversion at 100 eps.
            // Slopes are measured from current fields, not from the data file.
            // Native -O3 on this ARM build fuses (i+.5)*h-center (fnmsub);
            // coordinate rounding scales with the local result, not 224 M.
            initialization_error[n]=Budget::eps*(
                100.*std::abs(xc*ux+y*uy) +std::abs(xc*ux)+std::abs(y*uy));
        }
        double result[13];
        for(int mode=0;mode<2;++mode)
        {
        auto s = [&](int n,int x,int z)
        { double u=value(n,x,z); return Budget(u,Budget::eps*std::abs(u)+(mode ? initialization_error[n] : 0.)); };
        auto diff = [&](int n,int dir)
        {
            auto u = [&](int k) { return dir==0 ? s(n,k,0) : s(n,0,k); };
            return ((1./12.)*u(-2) -(2./3.)*u(-1) +(2./3.)*u(1) -(1./12.)*u(2))/h;
        };
        V v; D1 d; D2 dd;
        v.enum_mapping([&](int n,Budget &a) { a=s(n,0,0); });
        d.enum_mapping([&](int n,Tensor<1,Budget> &a) { FOR(i) a[i]=diff(n,i); });
        // Use the local metric magnitudes as the envelope for Gamma's stored
        // initialization error at the neighbouring cells. Two-cell variation
        // is separately checked in the Python census; no solution is supplied.
        auto gb=gamma_budget(v,d,y);
        for (int a=0;a<2;++a)
        {
            int n=c_Gamma1+a;
            v.Gamma[a].e += gb[a].e;
            for (int i=0;i<2;++i) d.Gamma[a][i].e += 1.5*gb[a].e/h;
        }
        dd.enum_mapping([&](int n,Tensor<2,Budget> &a)
        {
            FOR(i)
            {
                auto u = [&](int k) { return i==0 ? s(n,k,0) : s(n,0,k); };
                a[i][i]=(-(1./12.)*u(-2) +(4./3.)*u(-1) -2.5*u(0)
                          +(4./3.)*u(1) -(1./12.)*u(2))/(h*h);
            }
            Budget mix=0.;
            double w[5]={1./12.,-2./3.,0.,2./3.,-1./12.};
            for (int x=-2;x<=2;++x) for(int z=-2;z<=2;++z)
                if (x && z) mix += w[x+2]*w[z+2]*s(n,x,z);
            a[0][1]=a[1][0]=mix/(h*h);
        });
        auto c=audit.constraint_equations(v,d,dd,y);
        auto g=gauss(v,d,y);
        // Algebraic term separation, using only the saved current metric and
        // Gamma. Setting chi derivatives to zero isolates the metric part of
        // native Ricci; it is never used as a field, target or evolution input.
        auto md=d; auto mdd=dd;
        FOR(i) { md.chi[i]=0.; FOR(j) mdd.chi[i][j]=0.; }
        auto inv=TensorAlgebra::compute_inverse_sym(v.h);
        auto chris=TensorAlgebra::compute_christoffel(md.h,inv);
        auto mr=CCZ4CartoonGeometry::compute_ricci(v,md,mdd,inv,1./v.hww,chris,y);
        double metric_span=0.;
        for (int n : {c_h11,c_h22,c_hww})
        {
            double lo=s(n,0,0).v,hi=lo;
            for(int x=-2;x<=2;++x) for(int z=-2;z<=2;++z)
            { lo=std::min(lo,s(n,x,z).v); hi=std::max(hi,s(n,x,z).v); }
            metric_span=std::max(metric_span,hi-lo);
        }
        double baseline[10]={c.Ham.v,std::hypot(c.Mom[0].v,c.Mom[1].v),g.v,
                          c.Ham.e,std::hypot(c.Mom[0].e,c.Mom[1].e),g.e,
                          mr.scalar.v,c.Ham.v-mr.scalar.v,metric_span,
                          std::hypot(v.Gamma[0].v-gb[0].v,v.Gamma[1].v-gb[1].v)};
        if (!mode) std::copy(baseline,baseline+10,result);
        else { result[10]=c.Ham.e; result[11]=std::hypot(c.Mom[0].e,c.Mom[1].e); result[12]=g.e; }
        }
        out.write(reinterpret_cast<char*>(result),sizeof(result));
        ++cells;
    }
    if (!in.eof() || !out || cells==0) return 3;
    std::cout << cells << " current-state cells audited\n";
}
