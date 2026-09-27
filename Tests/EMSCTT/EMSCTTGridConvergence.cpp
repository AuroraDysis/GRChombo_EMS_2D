// Extends EMSTrumpetBinaryGridConvergence's t=0 production-operator path.
#include <cmath>
#include "EMSCouplingFunction.hpp"
#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ConstraintsCartoon.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include "ExperimentalGauge.hpp"
#include "FixSuperposition_metric.hpp"
#include "GammaCartoonCalculator.hpp"
#include "SetValue.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

struct region_t
{
    int count = 0, chi_floor = 0, lapse_floor = 0;
    double error = 0, ham_antisym_sum = 0;
    std::array<double, 6> sum{}, maximum{};
    void add(const FArrayBox &state, const FArrayBox &diag, const IntVect &iv,
             const CCZ4CartoonVars::VarsWithGauge<double> &v, int nx)
    {
        const double antisym =
            (diag(iv, c_Ham) - diag(IntVect(nx - 1 - iv[0], iv[1]), c_Ham)) / 2;
        ham_antisym_sum += antisym * antisym;
        const double values[] = {
            diag(iv, c_Ham),    diag(iv, c_Mom1),
            diag(iv, c_Mom2),   std::hypot(diag(iv, c_Mom1), diag(iv, c_Mom2)),
            diag(iv, c_GaussE), diag(iv, c_GaussB)};
        for (int k = 0; k < 6; ++k)
        {
            if (!std::isfinite(values[k]))
                throw std::runtime_error("nonfinite constraint");
            sum[k] += values[k] * values[k];
            maximum[k] = std::max(maximum[k], std::abs(values[k]));
        }
        const std::pair<int, double> fields[] = {
            {c_chi, v.chi},     {c_h11, v.h[0][0]},     {c_h12, v.h[0][1]},
            {c_h22, v.h[1][1]}, {c_hww, v.hww},         {c_K, v.K},
            {c_A11, v.A[0][0]}, {c_A12, v.A[0][1]},     {c_A22, v.A[1][1]},
            {c_Aww, v.Aww},     {c_phi, v.phi},         {c_Pi, v.Pi},
            {c_Ex, v.Ex},       {c_Ey, v.Ey},           {c_Ez, v.Ez},
            {c_Bx, v.Bx},       {c_By, v.By},           {c_Bz, v.Bz},
            {c_lapse, v.lapse}, {c_shift1, v.shift[0]}, {c_shift2, v.shift[1]}};
        for (const auto &f : fields)
        {
            const double e = std::abs(state(iv, f.first) - f.second) /
                             (1 + std::abs(f.second));
            if (!std::isfinite(e))
                throw std::runtime_error("nonfinite staging");
            error = std::max(error, e);
        }
        chi_floor += state(iv, c_chi) < 1e-6; // actual ConstraintsCartoon floor
        lapse_floor += state(iv, c_lapse) < 1e-12;
        ++count;
    }
    void print(const std::string &mode, double level, double dx,
               const char *name) const
    {
        if (!count || chi_floor || lapse_floor || error > 5e-13)
            throw std::runtime_error("empty/floored mask or incorrect staging");
        std::cout << mode << ',' << level << ',' << dx << ',' << name << ','
                  << count << ',' << chi_floor << ',' << lapse_floor << ','
                  << error;
        for (int k = 0; k < 6; ++k)
            std::cout << ',' << std::sqrt(sum[k] / count) << ',' << maximum[k];
        std::cout << ',' << std::sqrt(ham_antisym_sum / count) << std::endl;
    }
};
static void box_run(double level, bool outer, const std::string &profile,
                    const std::string &companion)
{
    const double dx = outer ? 4.0 / level : 1.0 / (32 * level);
    const double xmax = outer ? 96 : 18, ymax = outer ? 96 : 2;
    const int nx = std::lround(2 * xmax / dx), ny = std::lround(ymax / dx);
    CouplingFunction::params_t coupling{12.566370614359172, 0, 0, -20};
    EMSBH_params_t params{};
    params.bh_mass = 1;
    params.star_centre = {xmax, 0};
    params.boosted = params.binary = true;
    params.rapidity = .20273255;
    params.separation = 32;
    params.data_path = profile;
    params.ctt_data_path = companion;
    EMSBH_trumpet_read reader(params, coupling, 1., dx, 0);
    reader.compute_1d_solution();
    const Box interior(IntVect(0, 0), IntVect(nx - 1, ny - 1));
    Box ghost = interior;
    ghost.grow(5);
    FArrayBox state(ghost, NUM_VARS), diagnostic(interior, NUM_DIAGNOSTIC_VARS);
    BoxLoops::loop(make_compute_pack(SetValue(0.), reader), state, state,
                   disable_simd());
    BoxLoops::loop(
        FixSuperposition_metric(dx, params.star_centre, params.rapidity, true),
        state, state, disable_simd());
    Box gamma = interior;
    gamma.grow(2);
    BoxLoops::loop(GammaCartoonCalculator(dx), state, state, gamma,
                   disable_simd());
    ExperimentalGauge::params_t gauge{};
    gauge.shift_Gamma_coeff = .75;
    gauge.eta = 1;
    BoxLoops::loop(ExperimentalGauge(gauge), state, state, interior,
                   disable_simd());
    BoxLoops::loop(
        Constraints<CouplingFunction>(dx, CouplingFunction(coupling), 1.),
        state, diagnostic, interior, disable_simd());
    BoxLoops::loop(EMSCartoonGaussConstraints(dx, coupling), state, diagnostic,
                   interior, disable_simd());
    std::array<region_t, 7> regions;
    const double c = std::cosh(params.rapidity);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const double x = (i + .5) * dx - xmax, y = (j + .5) * dx;
            int index = -1;
            if (outer)
            {
                const double r = std::hypot(x, y);
                if (r >= 64 && r <= 96)
                    index = 6;
            }
            else
            {
                const double rl = std::hypot(c * (x + 16), y),
                             rr = std::hypot(c * (x - 16), y);
                // The old binary harness starts at .5, inside R_h=.63594.
                // Reuse the evolution harness's .75 exterior collar instead.
                if (rl >= .75 && rl <= 1.5)
                    index = 0;
                else if (rr >= .75 && rr <= 1.5)
                    index = 1;
                else if (std::abs(x) <= 14 && y <= 1.5)
                    index = 2;
            }
            if (index < 0)
                continue;
            const IntVect iv(i, j);
            const auto reference =
                reader.compute_binary_bh_vars(x, y, 1., 32., params.rapidity);
            regions[index].add(state, diagnostic, iv, reference, nx);
            if (!outer && y <= .125)
                regions[index + 3].add(state, diagnostic, iv, reference, nx);
        }
    const char *names[] = {"left",       "right",       "middle", "left_axis",
                           "right_axis", "middle_axis", "outer"};
    for (int k = outer ? 6 : 0; k < (outer ? 7 : 6); ++k)
        regions[k].print(companion.empty() ? "density" : "CTT", level, dx,
                         names[k]);
}
int main(int argc, char **argv)
{
    try
    {
        if (argc != 4)
            throw std::runtime_error("usage: EMSCTTGridConvergence refinement "
                                     "profile companion-or-density");
        const double level = std::stod(argv[1]);
        if (level != .5 && level != 1 && level != 2 && level != 4)
            throw std::runtime_error("refinement must be 0.5, 1, 2 or 4");
        const std::string companion =
            std::string(argv[3]) == "density" ? "" : argv[3];
        const auto start = std::chrono::steady_clock::now();
        std::cout << std::setprecision(17);
        std::cout
            << "mode,level,dx,region,count,chi_floor,lapse_floor,field_error,H_"
               "L2,H_Linf,Mx_L2,Mx_Linf,My_L2,My_Linf,M_L2,M_Linf,GaussE_L2,"
               "GaussE_Linf,GaussB_L2,GaussB_Linf,H_antisym_L2\n";
        box_run(level, false, argv[2], companion);
        box_run(level, true, argv[2], companion);
        std::cout << "seconds="
                  << std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - start)
                         .count()
                  << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
