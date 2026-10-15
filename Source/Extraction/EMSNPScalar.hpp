#ifndef EMS_NP_SCALAR_HPP
#define EMS_NP_SCALAR_HPP

#include <cmath>
#include "EMSMatterWeylOmScalar.hpp"
#include <complex>
#include <stdexcept>

// One convention for every scalar: (-+++), l=(n+u)/sqrt2,
// k=(n-u)/sqrt2, m=(v+i w)/sqrt2, positive Weyl contractions.
class EMSNPScalar : public EMSMatterWeylOmScalar
{
  public:
    using Complex = std::complex<double>;
    using Vector = std::array<Complex,4>;
    using Tetrad = std::array<Vector,4>; // l,k,m,mbar
    struct Geometry
    {
        double R[4][4][4][4]{}, C[4][4][4][4]{};
        double g[4][4]{}, ric[4][4]{}, F[4][4]{};
        Tetrad tetrad;
        double area_factor=NAN, legacy=NAN;
        bool valid=true;
        double electric_trace{}, magnetic_trace{}, magnetic_axial[3]{};
    };
    struct Sample
    {
        std::array<Complex,10> values; // Psi0..4,Phi0..2,Xout,Xin
    };
    static constexpr std::array<int,10> spins{{2,1,0,-1,-2,1,0,-1,-1,1}};
    static constexpr std::array<const char *,10> names{{
        "Psi0","Psi1","Psi2","Psi3","Psi4","Phi0","Phi1","Phi2","X_out","X_in"}};
    using EMSMatterWeylOmScalar::EMSMatterWeylOmScalar;

    Geometry geometry(const Vars<double> &vars,
        const Vars<Tensor<1,double>> &d1,
        const Diff2Vars<Tensor<2,double>> &d2,
        const Coordinates<double> &coords, bool matter=true) const
    {
        // Below this value the inherited tetrad and Riemann use different chi
        // floors, so they cannot define the requested normalized common frame.
        const double det2=vars.h[0][0]*vars.h[1][1]-vars.h[0][1]*vars.h[0][1];
        if (!(vars.chi>=1e-4 && vars.h[0][0]>0 && det2>0 && vars.hww>0 && coords.y>0))
        {
            Geometry invalid; invalid.valid=false;
            return invalid;
        }
        auto r=compute_riemann_tensor(vars,d1,d2,coords);
        if (matter) add_matter(r,vars,d1);
        const auto t=compute_null_tetrad(vars,coords);
        const auto hu=TensorAlgebra::compute_inverse_sym(vars.h);
        const auto chris=TensorAlgebra::compute_christoffel(d1.h,hu);
        const auto em=compute_EMS_EM_tensor(vars,d1,hu,1./vars.hww,chris,
            GR_SPACEDIM-CH_SPACEDIM,m_coupling);
        Geometry out;
        auto set=[&](int a,int b,int c,int d,double value) {
            if (a==b || c==d) return;
            out.R[a][b][c][d]=out.R[c][d][a][b]=out.R[b][a][d][c]=out.R[d][c][b][a]=value;
            out.R[b][a][c][d]=out.R[a][b][d][c]=out.R[d][c][a][b]=out.R[c][d][b][a]=-value;
        };
        FOR(i,j)
        {
            set(0,i+1,0,j+1,r.riemann_0d0d[i][j]);
            set(i+1,3,j+1,3,r.riemann_dwdw[i][j]);
            FOR(k)
            {
                set(0,i+1,j+1,k+1,r.riemann_0ddd[i][j][k]);
                FOR(l) { set(i+1,j+1,k+1,l+1,r.riemann_dddd[i][j][k][l]); }
            }
        }
        set(0,3,0,3,r.riemann_0w0w);
        FOR(i) { set(0,3,i+1,3,r.riemann_0wdw[i]); }
        auto &g=out.g; auto &ric=out.ric;
        g[0][0]=-1.;
        FOR(i,j) { g[i+1][j+1]=vars.h[i][j]/vars.chi; }
        g[3][3]=vars.hww/vars.chi;
        const double eightpiG=8*M_PI*m_G_Newton, trace=em.S-em.rho;
        ric[0][0]=eightpiG*(em.rho+em.S)/2;
        FOR(i)
        {
            // Pi=-n(phi); the evolution's Si is -T_ni.
            ric[0][i+1]=ric[i+1][0]=-eightpiG*em.Si[i];
            FOR(j) { ric[i+1][j+1]=eightpiG*(em.Sij[i][j]-trace*g[i+1][j+1]/2); }
        }
        ric[3][3]=eightpiG*(em.Sww-trace*g[3][3]/2);
        const double scalar=-eightpiG*trace;
        // The vacuum control also omits Ricci here, matching its normal block Q.
        const double ricci_weight=matter ? 1. : 0.;
        // Electric Weyl before TF: its anisotropic part is Q_ij-4piG S_ij,
        // Q=3Ric+K K_ij-K_ik K^k_j. Raw R_injn alone would give -8piG S_ij.
        const double frame[3][3]{{t.m1[0],t.m1[1],0},
                                 {t.m2[0],t.m2[1],0},{0,0,t.mw}};
        double electric[3][3]{}, magnetic[3][3]{};
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        for (int a=0;a<3;++a) for (int b=0;b<3;++b)
        {
            electric[i][j]+=frame[i][a]*frame[j][b]*(out.R[0][a+1][0][b+1]
                +ricci_weight*(ric[a+1][b+1]/2-g[a+1][b+1]*(ric[0][0]/2+scalar/6)));
            for (int k=0;k<3;++k) for (int l=0;l<3;++l) for (int c=0;c<3;++c)
                magnetic[i][j]+=.5*epsilon(j,k,l)*frame[i][a]*frame[k][b]*frame[l][c]*
                    (out.R[0][a+1][b+1][c+1]+ricci_weight*(g[a+1][b+1]*ric[0][c+1]/2
                     -g[a+1][c+1]*ric[0][b+1]/2));
        }
        rebuild_weyl(out,electric,magnetic,frame);
        const double E[3]{vars.Ex,vars.Ey,vars.Ez}, B[3]{vars.Bx,vars.By,vars.Bz};
        double bu[3]{};
        FOR(i,j) { bu[i]+=vars.chi*hu[i][j]*B[j]; }
        bu[2]=vars.chi*B[2]/vars.hww;
        const double rootdet=std::sqrt((g[1][1]*g[2][2]-g[1][2]*g[2][1])*g[3][3]);
        for (int i=0;i<3;++i) { out.F[i+1][0]=E[i]; out.F[0][i+1]=-E[i]; }
        out.F[1][2]=rootdet*bu[2]; out.F[2][3]=rootdet*bu[0]; out.F[1][3]=-rootdet*bu[1];
        for (int i=1;i<4;++i) for (int j=i+1;j<4;++j) out.F[j][i]=-out.F[i][j];
        const double s=1/std::sqrt(2.);
        out.tetrad={Vector{s,s*t.m1[0],s*t.m1[1],0},Vector{s,-s*t.m1[0],-s*t.m1[1],0},
            Vector{0,s*t.m2[0],s*t.m2[1],Complex(0,s*t.mw)},
            Vector{0,s*t.m2[0],s*t.m2[1],Complex(0,-s*t.mw)}};
        const double radius=std::hypot(coords.x,coords.y);
        const double vx=-coords.y/radius, vy=coords.x/radius;
        out.area_factor=std::sqrt((g[1][1]*vx*vx+2*g[1][2]*vx*vy+g[2][2]*vy*vy)*g[3][3]);
        out.legacy=compute_Weyl_Om(r,vars,d1,d2,coords).Om22;
        return out;
    }

    static int epsilon(int i,int j,int k)
    { return (i-j)*(j-k)*(k-i)/2; }

    // Oriented orthonormal components; the STF operation uses the full x,y,w
    // spatial metric through the frame, never the reduced 2x2 trace.
    static void rebuild_weyl(Geometry &out,const double (&e)[3][3],
        const double (&b)[3][3],const double (&frame)[3][3])
    {
        out.electric_trace=out.magnetic_trace=0;
        for (double &x:out.magnetic_axial) x=0;
        for (int i=0;i<3;++i)
        {
            out.electric_trace+=e[i][i]; out.magnetic_trace+=b[i][i];
            for (int j=0;j<3;++j) for (int k=0;k<3;++k)
                out.magnetic_axial[i]+=epsilon(i,j,k)*b[j][k];
        }
        double coframe[3][3]{}, ec[3][3]{}, bc[3][3]{}, gu[3][3]{};
        for (int i=0;i<3;++i) for (int a=0;a<3;++a) for (int c=0;c<3;++c)
        {
            coframe[i][a]+=out.g[a+1][c+1]*frame[i][c];
            gu[a][c]+=frame[i][a]*frame[i][c];
        }
        for (int a=0;a<3;++a) for (int c=0;c<3;++c)
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        {
            ec[a][c]+=coframe[i][a]*coframe[j][c]*
                ((e[i][j]+e[j][i])/2-(i==j ? out.electric_trace/3 : 0));
            bc[a][c]+=coframe[i][a]*coframe[j][c]*
                ((b[i][j]+b[j][i])/2-(i==j ? out.magnetic_trace/3 : 0));
        }
        const double rootdet=1/(frame[0][0]*frame[1][1]*frame[2][2]
                              -frame[0][1]*frame[1][0]*frame[2][2]);
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        {
            out.C[0][i+1][0][j+1]=out.C[i+1][0][j+1][0]=ec[i][j];
            out.C[i+1][0][0][j+1]=out.C[0][i+1][j+1][0]=-ec[i][j];
            for (int k=0;k<3;++k)
            {
                double mixed=0;
                for (int l=0;l<3;++l) for (int m=0;m<3;++m)
                    mixed+=rootdet*epsilon(j,k,l)*gu[l][m]*bc[i][m];
                out.C[0][i+1][j+1][k+1]=out.C[j+1][k+1][0][i+1]=mixed;
                out.C[i+1][0][j+1][k+1]=out.C[j+1][k+1][i+1][0]=-mixed;
                for (int l=0;l<3;++l)
                    out.C[i+1][j+1][k+1][l+1]=out.g[i+1][k+1]*ec[j][l]
                        -out.g[i+1][l+1]*ec[j][k]-out.g[j+1][k+1]*ec[i][l]
                        +out.g[j+1][l+1]*ec[i][k];
            }
        }
    }

    static Complex contract(const double (&C)[4][4][4][4],
        const Vector &u,const Vector &v,const Vector &w,const Vector &z)
    {
        Complex value=0;
        for (int a=0;a<4;++a) for (int b=0;b<4;++b)
        for (int c=0;c<4;++c) for (int d=0;d<4;++d)
            value+=C[a][b][c][d]*u[a]*v[b]*w[c]*z[d];
        return value;
    }
    static Complex contract(const double (&F)[4][4],const Vector &u,const Vector &v)
    {
        Complex value=0;
        for (int a=0;a<4;++a) for (int b=0;b<4;++b) value+=F[a][b]*u[a]*v[b];
        return value;
    }
    static Sample evaluate(const Geometry &g,const Tetrad &t)
    {
        if (!g.valid) return invalid_sample();
        const auto &l=t[0], &k=t[1], &m=t[2], &mb=t[3];
        Sample out;
        auto &v=out.values;
        v[0]=contract(g.C,l,m,l,m); v[1]=contract(g.C,l,k,l,m);
        v[2]=contract(g.C,l,m,mb,k); v[3]=contract(g.C,l,k,mb,k);
        v[4]=contract(g.C,k,mb,k,mb);
        v[5]=contract(g.F,l,m); v[6]=(contract(g.F,l,k)+contract(g.F,mb,m))/2.;
        v[7]=contract(g.F,mb,k);
        v[8]=2.*v[6]*v[3]-3.*v[2]*v[7]; v[9]=2.*v[6]*v[1]-3.*v[2]*v[5];
        for (const auto z:v)
            if (!std::isfinite(z.real()) || !std::isfinite(z.imag()))
                return invalid_sample();
        return out;
    }
    static Sample invalid_sample()
    {
        Sample out;
        out.values.fill(Complex(NAN,NAN));
        return out;
    }
};
#endif
