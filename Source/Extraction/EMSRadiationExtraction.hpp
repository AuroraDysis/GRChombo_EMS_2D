#ifndef EMS_RADIATION_EXTRACTION_HPP
#define EMS_RADIATION_EXTRACTION_HPP

#include "EMSNativeRadiation.hpp"
#include "EMSCouplingFunction.hpp"
#include "SphericalExtraction.hpp"
#include <fstream>
#include <iomanip>

// Same spheres/interpolator/harmonics as the author's extraction. Axisymmetry
// makes all m != 0 vanish; record every allowed l <= 6, including scalar l=0.
class EMSRadiationExtraction : public SphericalExtraction
{
    CouplingFunction::params_t m_c;
    double m_phi_inf;
  public:
    EMSRadiationExtraction(const spherical_extraction_params_t &p, double dt,
        double time, bool first, double restart, CouplingFunction::params_t c,
        double phi_inf)
        : SphericalExtraction(p,dt,time,first,restart), m_c(c), m_phi_inf(phi_inf)
    {
        add_evolution_vars({c_phi,c_Pi,c_lapse,c_shift1,c_shift2,c_chi,
                            c_h11,c_h12,c_h22,c_hww,c_Ex,c_Ey,c_Ez,c_Bx,c_By,c_Bz});
        add_var(c_phi,VariableType::evolution,Derivative::dx);
        add_var(c_phi,VariableType::evolution,Derivative::dy);
        add_diagnostic_vars({c_Weyl4_Re,c_Weyl4_Im});
    }

    void execute_query(AMRInterpolator<Lagrange<4>> *interp)
    { extract(interp); write_modes(); }

    void write_modes()
    {
        if (procID() != 0) return;
        const std::string path = m_params.data_path + "ems_radiation.csv";
        // Append-only across restarts: the offline reader keeps the last
        // complete (time,r,l) record, just like the native small-data readers.
        std::ofstream file(path, m_first_step && m_restart_time == 0 ?
                           std::ios::trunc : std::ios::app);
        if (!file) throw std::runtime_error("cannot write " + path);
        if (file.tellp() == 0)
            file << "t,r,R,l,m,S_re,S_im,D_re,D_im,P_re,P_im,W_re,W_im,"
                    "L_scalar_stress,L_EM_stress,alpha2,F_inf,Pi_rms\n";
        file << std::setprecision(17);
        const double finf = std::exp(-2*m_c.alpha*(m_c.f0+m_c.f1*m_phi_inf
                                                  +m_c.f2*m_phi_inf*m_phi_inf));
        for (int ir=0; ir<m_params.num_surfaces; ++ir)
        {
            double area=0, alpha2=0, ls=0, le=0, pi2=0;
            std::array<double,7> sm{}, dm{};
            std::array<std::complex<double>,7> pm{}, wm{};
            const double r = m_params.surface_param_values[ir];
            for (int j=0; j<m_params.num_points_u; ++j)
            {
                const int q = index(ir,j);
                auto at = [&](int k) { return m_interp_data[k][q]; };
                const double theta = m_geom.u(j,m_params.num_points_u);
                const double weight = 2*M_PI*m_du*std::sin(theta)*
                    IntegrationMethod::simpson.weight(j,m_params.num_points_u,false);
                EMSNativeRadiation::Fields f;
                f.phi=at(0); f.Pi=at(1); f.lapse=at(2);
                f.shift={at(3),at(4),0}; f.grad={at(16),at(17),0};
                f.g={{{at(6)/at(5),at(7)/at(5),0},
                       {at(7)/at(5),at(8)/at(5),0},{0,0,at(9)/at(5)}}};
                f.E={at(10),at(11),at(12)}; f.B={at(13),at(14),at(15)};
                const double coupling=std::exp(-2*m_c.alpha*(m_c.f0+
                                        m_c.f1*f.phi+m_c.f2*f.phi*f.phi));
                const auto v=EMSNativeRadiation::evaluate(f,theta,coupling);
                area += weight*r*r*v.area_factor;
                alpha2 += weight*r*r*v.area_factor*f.lapse*f.lapse;
                pi2 += weight*r*r*v.area_factor*f.Pi*f.Pi;
                ls += weight*r*r*v.scalar_flux;
                le += weight*r*r*v.em_flux;
                for (int l=0; l<=6; ++l)
                {
                    auto harmonic = [&](int spin) {
                        return SphericalHarmonics::spin_Y_lm(std::sin(theta),
                                      0.,std::cos(theta),spin,l,0).Real; };
                    sm[l] += weight*harmonic(0)*(f.phi-m_phi_inf);
                    dm[l] += weight*harmonic(0)*v.dtphi;
                    if (l>=1) pm[l] += weight*harmonic(-1)*v.phi2;
                    if (l>=2) wm[l] += weight*harmonic(-2)*
                                              std::complex<double>(at(18),at(19));
                }
            }
            const double R=std::sqrt(area/(4*M_PI));
            for (int l=0; l<=6; ++l)
                file << m_time << ',' << r << ',' << R << ',' << l << ",0,"
                     << R*sm[l] << ",0," << R*dm[l] << ",0,"
                     << R*pm[l].real() << ',' << R*pm[l].imag() << ','
                     << R*wm[l].real() << ',' << R*wm[l].imag() << ','
                     << ls << ',' << le << ',' << alpha2/area << ',' << finf
                     << ',' << std::sqrt(pi2/area) << '\n';
        }
        if (!file) throw std::runtime_error("failed writing " + path);
    }
};
#endif
