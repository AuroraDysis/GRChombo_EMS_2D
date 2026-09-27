#ifndef EMS_NATIVE_RADIATION_HPP
#define EMS_NATIVE_RADIATION_HPP

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

// Native fields at y >= 0, z = 0. E and B are covectors. Polar axis: x.
namespace EMSNativeRadiation
{
using Vec = std::array<double, 3>;
using Mat = std::array<Vec, 3>;
inline double dot(const Vec &a, const Vec &b)
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline Vec mul(const Mat &a, const Vec &b)
{ return {dot(a[0], b), dot(a[1], b), dot(a[2], b)}; }
inline Vec cross(const Vec &a, const Vec &b)
{ return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}; }

struct Fields
{
    double phi = 0, Pi = 0, lapse = 1;
    Vec shift{}, grad{}, E{}, B{};
    Mat g{{{1,0,0}, {0,1,0}, {0,0,1}}};
};
struct Sample
{
    double dtphi, scalar_flux, em_flux, area_factor;
    std::complex<double> phi2;
};

inline Sample evaluate(const Fields &f, double theta, double coupling)
{
    const double det2 = f.g[0][0]*f.g[1][1]-f.g[0][1]*f.g[0][1];
    const double det = det2*f.g[2][2];
    if (!(f.lapse > 0 && f.g[0][0] > 0 && det2 > 0 && det > 0 && coupling > 0))
        throw std::runtime_error("EMS radiation: invalid surface metric/lapse/coupling");
    Mat inv{{{f.g[1][1]/det2,-f.g[0][1]/det2,0},
             {-f.g[0][1]/det2,f.g[0][0]/det2,0}, {0,0,1/f.g[2][2]}}};
    const Vec dr{std::cos(theta), std::sin(theta), 0};
    Vec et{-std::sin(theta), std::cos(theta), 0};
    const double gt = dot(et, mul(f.g,et));
    for (auto &x : et) x /= std::sqrt(gt);
    const Vec ep{0,0,1/std::sqrt(f.g[2][2])};
    // s = gamma^{-1} dr / |dr|; (s,e_theta,e_phi) is right handed.
    // k=(n-s)/sqrt(2), m=(e_theta+i e_phi)/sqrt(2).
    // F_ij=epsilon_ijk B^k, E_i=F_ia n^a. No sqrt(coupling) in Phi2.
    Sample out{};
    out.phi2 = {0.5*(dot(f.E,et)+dot(f.B,ep)),
                0.5*(dot(f.B,et)-dot(f.E,ep))};
    out.dtphi = -f.lapse*f.Pi + dot(f.shift,f.grad);
    auto grad_up = mul(inv,f.grad);
    Vec js{};
    for (int i=0; i<3; ++i)
        js[i] = -2*out.dtphi*(grad_up[i]-f.shift[i]*f.Pi/f.lapse);
    const auto eu = mul(inv,f.E), bu = mul(inv,f.B);
    const double rho = coupling*(dot(f.E,eu)+dot(f.B,bu));
    auto su = cross(f.E,f.B);
    for (auto &x : su) x *= 2*coupling/std::sqrt(det);
    const double beta_s = dot(f.shift,mul(f.g,su));
    Vec je{};
    for (int i=0; i<3; ++i)
    {
        const double stress_beta = rho*f.shift[i]
            -2*coupling*(eu[i]*dot(f.E,f.shift)+bu[i]*dot(f.B,f.shift));
        // -T^i_t = alpha S^i - beta^i rho - S^i_j beta^j
        //            + beta^i beta^j S_j/alpha.
        je[i] = f.lapse*su[i]-f.shift[i]*rho-stress_beta
                +f.shift[i]*beta_s/f.lapse;
    }
    // Multiply by r^2 dOmega in the extractor, not by the proper area again.
    out.scalar_flux = f.lapse*std::sqrt(det)*dot(dr,js);
    out.em_flux = f.lapse*std::sqrt(det)*dot(dr,je);
    out.area_factor = std::sqrt(gt*f.g[2][2]);
    return out;
}
}
#endif
