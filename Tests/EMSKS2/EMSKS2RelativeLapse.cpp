#include "EMSKS2ReferenceCache.hpp"
#include "FourthOrderDerivatives.hpp"
#include "BoxIterator.H"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>

template <typename T> struct GaugeVars
{
    T lapse{}, K{}, Theta{};
    Tensor<1, T> shift{};
};

int main(int argc, char **argv)
{
    if (argc != 3) return 2;
    const auto start = std::chrono::steady_clock::now();
    const int N = std::stoi(argv[1]);
    EMSKS2Profile profile;
    profile.load(argv[2]);
    const double M = profile.get("M"), rh = profile.get("r_h"), dx = M / N;
    const double ra = profile.rho_a(), rm = profile.get("r_m");
    const int half = int(std::ceil(2.5 * rh / dx));
    const Box valid(IntVect(0, 0), IntVect(2 * half - 1, half - 1));
    Box ghost = valid;
    ghost.grow(3);
    FArrayBox ref(ghost, ReferenceStationaryGauge::num_reference);
    EMSKS2ReferenceCache::fill(ref, profile, dx, {half * dx, 0.});
    FArrayBox alpha(ghost, 1), u(ghost, 1), alpha_ref(ghost, 1);
    for (BoxIterator bit(ghost); bit.ok(); ++bit)
        alpha_ref(bit(), 0) = ref(bit(), ReferenceStationaryGauge::alpha_star);
    FourthOrderDerivatives deriv(dx);
    ReferenceStationaryGauge::params_t gauge_params;
    gauge_params.reference = &ref;
    gauge_params.relative_lapse = &u;
    gauge_params.mass = M;
    ReferenceStationaryGauge gauge(gauge_params);
    std::cout << std::setprecision(17)
              << "N,amplitude,mask,count,L2,Linf,tau_L2,tau_Linf,seconds\n";
    for (double amplitude : {1e-3, 1e-1})
    {
        for (BoxIterator bit(ghost); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            const double x = (iv[0] + .5 - half) * dx;
            const double y = (iv[1] + .5) * dx;
            const double z = (x * x + y * y - rh * rh) / (.9 * rh * rh);
            const double um = amplitude * std::exp(-z * z) * (1 + .2 * x / rh);
            alpha(iv, 0) = ref(iv, ReferenceStationaryGauge::alpha_star) *
                           std::exp(um);
            u(iv, 0) = std::log(alpha(iv, 0) /
                                 ref(iv, ReferenceStationaryGauge::alpha_star));
        }
        double squares[4]{}, maxima[4]{}, tau_squares[4]{}, tau_maxima[4]{};
        int counts[4]{};
        for (BoxIterator bit(valid); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            const double x = (iv[0] + .5 - half) * dx;
            const double y = (iv[1] + .5) * dx;
            const double rho = std::hypot(x, y);
            if (rho - 4 * dx < ra / 2 || rho > 2.3 * rh) continue;
            const int mask = rho < ra ? -1 : rho < rm ? 0 :
                             rho < rh ? 1 : rho < 2 * rh ? 2 : 3;
            if (mask < 0) continue;
            Tensor<1, double> beta;
            beta[0] = ref(iv, ReferenceStationaryGauge::beta_x);
            beta[1] = ref(iv, ReferenceStationaryGauge::beta_y);
            GaugeVars<double> vars, rhs;
            vars.lapse = alpha(iv, 0);
            const double z = (rho * rho - rh * rh) / (.9 * rh * rh);
            const double bump = amplitude * std::exp(-z * z);
            vars.K = ref(iv, ReferenceStationaryGauge::K_star) + .03 * bump;
            vars.Theta = .01 * bump;
            vars.shift = beta;
            gauge.set_relative_lapse_rhs(rhs, vars, iv, deriv, 0.);
            const double factor = 1 + .2 * x / rh;
            const double dzdx = 2 * x / (.9 * rh * rh);
            const double dzdy = 2 * y / (.9 * rh * rh);
            const double dudx = bump * (-2 * z * dzdx * factor + .2 / rh);
            const double dudy = bump * (-2 * z * dzdy * factor);
            const double continuum = alpha(iv, 0) *
                (beta[0] * dudx + beta[1] * dudy -
                 alpha(iv, 0) * (vars.K -
                     ref(iv, ReferenceStationaryGauge::K_star) -
                     2 * vars.Theta) - u(iv, 0) / M);
            const double error = rhs.lapse - continuum;
            squares[mask] += error * error;
            maxima[mask] = std::max(maxima[mask], std::abs(error));
            // Unperturbed reference, with tapered KO at sigma=1.
            const auto r = deriv.scalar_advection_dissipation(
                alpha_ref, iv, beta,
                ref(iv, ReferenceStationaryGauge::taper));
            const double t = (r.first + r.second) /
                ref(iv, ReferenceStationaryGauge::alpha_star) -
                ref(iv, ReferenceStationaryGauge::q_star);
            tau_squares[mask] += t * t;
            tau_maxima[mask] = std::max(tau_maxima[mask], std::abs(t));
            ++counts[mask];
        }
        const char *names[] = {"join", "KS_collar", "exterior", "far"};
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        struct rusage usage{};
        getrusage(RUSAGE_SELF, &usage);
        for (int m = 0; m < 4; ++m)
        {
            if (!counts[m] || !std::isfinite(maxima[m]) ||
                maxima[m] > 1e-6 || !std::isfinite(tau_maxima[m]))
                throw std::runtime_error("relative lapse manufactured check failed");
            std::cout << N << ',' << amplitude << ',' << names[m] << ','
                      << counts[m] << ',' << std::sqrt(squares[m] / counts[m])
                      << ',' << maxima[m] << ','
                      << std::sqrt(tau_squares[m] / counts[m]) << ','
                      << tau_maxima[m] << ',' << elapsed << '\n';
        }
        std::cerr << "stage=amplitude elapsed=" << elapsed
                  << " peak_RSS_bytes=" << usage.ru_maxrss
                  << " items=" << (amplitude == 1e-3 ? 1 : 2) << "/2\n";
    }
}
