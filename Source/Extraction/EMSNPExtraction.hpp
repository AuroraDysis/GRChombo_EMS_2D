#ifndef EMS_NP_EXTRACTION_HPP
#define EMS_NP_EXTRACTION_HPP

#include "EMSNPScalar.hpp"
#include "SphericalExtraction.hpp"
#include <fstream>
#include <iomanip>
#include <map>

// Reconstruct on interpolated fourth-order jets, avoiding any extra permanent
// diagnostic variables or changes to checkpoint/plot component layouts.
class EMSNPExtraction : public SphericalExtraction
{
    EMSNPScalar m_compute;
    bool m_matter;
    std::map<std::pair<int,Derivative>,int> m_indices;
  protected:
    void request(int var,const Derivative &deriv=Derivative::LOCAL)
    {
        const auto key=std::make_pair(var,deriv);
        if (m_indices.count(key)) return;
        m_indices[key]=m_vars.size();
        add_var(var,VariableType::evolution,deriv);
    }
  public:
    EMSNPExtraction(const spherical_extraction_params_t &p,double dt,double time,
        bool first,double restart,CouplingFunction::params_t coupling,double G,bool matter=true)
        : SphericalExtraction(p,dt,time,first,restart),
          m_compute(p.center,1.,coupling,G),m_matter(matter)
    {
        EMSNPScalar::Vars<double> vars{};
        vars.enum_mapping([&](int var,double &) { request(var); });
        for (int var:{c_chi,c_h11,c_h12,c_h22,c_hww,c_A11,c_A12,c_A22,c_Aww,c_K,c_phi})
            for (const auto &deriv:{Derivative::dx,Derivative::dy}) request(var,deriv);
        for (int var:{c_chi,c_h11,c_h12,c_h22,c_hww})
            for (const auto &deriv:{Derivative::dxdx,Derivative::dxdy,Derivative::dydy}) request(var,deriv);
    }

    struct Jets
    {
        EMSNPScalar::Vars<double> v{};
        EMSNPScalar::Vars<Tensor<1,double>> d1{};
        EMSNPScalar::Diff2Vars<Tensor<2,double>> d2{};
        Coordinates<double> coords{IntVect::Zero,1.};
    };
    Jets point_jets(int ir,int j) const
    {
        const int q=index(ir,j);
        auto at=[&](int var,const Derivative &deriv=Derivative::LOCAL) {
            const auto it=m_indices.find(std::make_pair(var,deriv));
            return it==m_indices.end() ? 0. : m_interp_data[it->second][q];
        };
        Jets jets;
        auto &v=jets.v; auto &d1=jets.d1; auto &d2=jets.d2;
        v.enum_mapping([&](int var,double &value) { value=at(var); });
        d1.enum_mapping([&](int var,Tensor<1,double> &value) {
            value[0]=at(var,Derivative::dx); value[1]=at(var,Derivative::dy); });
        d2.enum_mapping([&](int var,Tensor<2,double> &value) {
            value[0][0]=at(var,Derivative::dxdx); value[1][1]=at(var,Derivative::dydy);
            value[0][1]=value[1][0]=at(var,Derivative::dxdy); });
        auto &coords=jets.coords;
        coords.x=m_interp_coords[0][q]-m_center[0];
        coords.y=m_interp_coords[1][q]-m_center[1];
        return jets;
    }
    EMSNPScalar::Geometry point(int ir,int j) const
    {
        const auto jets=point_jets(ir,j);
        return m_compute.geometry(jets.v,jets.d1,jets.d2,jets.coords,m_matter);
    }

    template <class InterpAlgo>
    void execute_query(AMRInterpolator<InterpAlgo> *interp)
    { extract(interp); write_modes(); }

    void write_modes()
    {
        if (procID()!=0) return;
        const std::string path=m_params.data_path+"ems_np.csv";
        std::ofstream file(path,m_first_step && m_restart_time==0 ? std::ios::trunc : std::ios::app);
        if (!file) { pout()<<"EMS NP warning: cannot write "<<path<<" at t="<<m_time<<'\n'; return; }
        if (file.tellp()==0)
        {
            file << "t,r,R,l,m";
            for (const char *name:EMSNPScalar::names) file << ',' << name << "_re," << name << "_im";
            file << '\n';
        }
        file << std::setprecision(17);
        for (int ir=0;ir<m_params.num_surfaces;++ir)
        {
            double area=0;
            bool valid=true;
            std::array<std::array<EMSNPScalar::Complex,10>,7> modes{};
            const double r=m_params.surface_param_values[ir];
            // Poles have zero quadrature measure and the polar frame is singular.
            for (int j=1;j<m_params.num_points_u-1;++j)
            {
                const double theta=m_geom.u(j,m_params.num_points_u);
                const double weight=2*M_PI*m_du*std::sin(theta)*
                    IntegrationMethod::simpson.weight(j,m_params.num_points_u,false);
                const auto g=point(ir,j);
                if (!g.valid || !std::isfinite(g.area_factor)) { valid=false; break; }
                const auto value=EMSNPScalar::evaluate(g,g.tetrad);
                for (auto z:value.values)
                    valid=valid && std::isfinite(z.real()) && std::isfinite(z.imag());
                if (!valid) break;
                area+=weight*r*r*g.area_factor;
                for (int l=0;l<=6;++l) for (int n=0;n<10;++n)
                    if (l>=std::abs(EMSNPScalar::spins[n]))
                        modes[l][n]+=weight*SphericalHarmonics::spin_Y_lm(std::sin(theta),0.,
                            std::cos(theta),EMSNPScalar::spins[n],l,0).Real*value.values[n];
            }
            double R=std::sqrt(area/(4*M_PI));
            valid=valid && R>0 && std::isfinite(R);
            for (const auto &mode:modes) for (auto z:mode)
                valid=valid && std::isfinite(R*z.real()) && std::isfinite(R*z.imag());
            if (!valid)
            {
                R=NAN;
                for (auto &mode:modes) mode.fill(EMSNPScalar::Complex(NAN,NAN));
                static std::map<std::pair<std::string,double>,double> warned;
                const auto key=std::make_pair(path,r);
                if (!warned.count(key) || warned[key]!=m_time)
                {
                    pout()<<"EMS NP warning: invalid sphere r="<<r<<" at t="<<m_time
                          <<"; writing NaN modes\n";
                    warned[key]=m_time;
                }
            }
            for (int l=0;l<=6;++l)
            {
                file << m_time << ',' << r << ',' << R << ',' << l << ",0";
                for (const auto value:modes[l]) file << ',' << R*value.real() << ',' << R*value.imag();
                file << '\n';
            }
        }
        file.close();
        if (!file) pout()<<"EMS NP warning: failed writing "<<path<<" at t="<<m_time<<'\n';
    }
};
#endif
