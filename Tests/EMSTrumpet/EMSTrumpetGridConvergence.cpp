/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <cmath>

#include "EMSCouplingFunction.hpp"

#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ConstraintsCartoon.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include "ExperimentalGauge.hpp"
#include "GammaCartoonCalculator.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <utility>

struct norm_t
{
    double sum = 0, max = 0;
    void accumulate(double a_value)
    {
        sum += a_value * a_value;
        max = std::max(max, std::abs(a_value));
    }
};
struct region_t
{
    std::array<norm_t, 5> norms{};
    int count = 0, chi_floor = 0, lapse_floor = 0;
    double field_error = 0;
    void accumulate(const FArrayBox &a_state, const FArrayBox &a_diagnostics,
                    const IntVect &a_iv,
                    const CCZ4CartoonVars::VarsWithGauge<double> &a_vars)
    {
        const int vars[5] = {c_Ham, c_Mom1, c_Mom2, c_GaussE, c_GaussB};
        for (int k = 0; k < 5; ++k)
            norms[k].accumulate(a_diagnostics(a_iv, vars[k]));
        const std::pair<int, double> fields[] = {{c_chi, a_vars.chi},
                                                 {c_h11, a_vars.h[0][0]},
                                                 {c_h12, a_vars.h[0][1]},
                                                 {c_h22, a_vars.h[1][1]},
                                                 {c_hww, a_vars.hww},
                                                 {c_K, a_vars.K},
                                                 {c_A11, a_vars.A[0][0]},
                                                 {c_A12, a_vars.A[0][1]},
                                                 {c_A22, a_vars.A[1][1]},
                                                 {c_Aww, a_vars.Aww},
                                                 {c_lapse, a_vars.lapse},
                                                 {c_shift1, a_vars.shift[0]},
                                                 {c_shift2, a_vars.shift[1]},
                                                 {c_phi, a_vars.phi},
                                                 {c_Pi, a_vars.Pi},
                                                 {c_Ex, a_vars.Ex},
                                                 {c_Ey, a_vars.Ey},
                                                 {c_Ez, a_vars.Ez},
                                                 {c_Bx, a_vars.Bx},
                                                 {c_By, a_vars.By},
                                                 {c_Bz, a_vars.Bz}};
        for (const auto &f : fields)
            field_error = std::max(field_error,
                                   std::abs(a_state(a_iv, f.first) - f.second) /
                                       (1 + std::abs(f.second)));
        chi_floor += a_state(a_iv, c_chi) < 1e-12;
        lapse_floor += a_state(a_iv, c_lapse) < 1e-12;
        ++count;
    }
    void print_summary(const char *a_name) const
    {
        std::cout << a_name << ',' << count << ',' << chi_floor << ','
                  << lapse_floor << ',' << field_error;
        for (const auto &n : norms)
            std::cout << ',' << std::sqrt(n.sum / count) << ',' << n.max;
        std::cout << '\n';
    }
};
int main(int argc, char **argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "usage: EMSTrumpetGridConvergence N eta profile.trumpet "
                     "[echo]\n";
        return 2;
    }
    const bool echo = argc == 5 && std::string(argv[4]) == "echo";
    if (argc == 5 && !echo)
        return 2;
    const int N = std::stoi(argv[1]);
    const double eta = std::stod(argv[2]), dx = 1.0 / N;
    if (echo ? (N != 128 && N != 256 && N != 512)
             : (N != 32 && N != 64 && N != 128))
    {
        std::cerr << "unsupported N\n";
        return 2;
    }
    const auto start = std::chrono::steady_clock::now();
    CouplingFunction::params_t coupling{12.566370614359172, 0, 0, -20};
    if (echo)
    {
        EMSTrumpetSolution_read profile;
        std::string reason;
        if (!profile.check_file(argv[3], reason))
        {
            std::cerr << "invalid trumpet: " << reason << '\n';
            return 2;
        }
        const auto c = profile.get_coupling_parameters();
        coupling = {c[0], c[1], c[2], c[3]};
    }

    EMSBH_params_t params{};
    params.bh_mass = 1;
    params.star_centre = {echo ? 0.5 : 2.0, 0};
    params.data_path = argv[3];
    params.boosted = eta != 0;
    params.rapidity = eta;
    EMSBH_trumpet_read reader(params, coupling, 1.0, dx, 0);
    reader.compute_1d_solution();
    const Box interior(IntVect(0, 0), IntVect((echo ? N : 4 * N) - 1,
                                              (echo ? N / 2 : 2 * N) - 1));
    Box ghost = interior;
    ghost.grow(5);
    FArrayBox state(ghost, NUM_VARS), diagnostic(interior, NUM_DIAGNOSTIC_VARS);
    BoxLoops::loop(make_compute_pack(SetValue(0.0), reader), state, state,
                   disable_simd());
    // Analytic ghost cells give the reflected Gamma values needed by d1.Gamma
    // at the first two axis rows, as the level's second ghost fill does.
    Box gamma_box = interior;
    gamma_box.grow(2);
    BoxLoops::loop(GammaCartoonCalculator(dx), state, state, gamma_box,
                   disable_simd());
    ExperimentalGauge::params_t gauge{};
    gauge.shift_Gamma_coeff = 0.75;
    gauge.eta = 1.0;
    BoxLoops::loop(ExperimentalGauge(gauge), state, state, interior,
                   disable_simd());
    BoxLoops::loop(
        Constraints<CouplingFunction>(dx, CouplingFunction(coupling), 1.0),
        state, diagnostic, interior, disable_simd());
    BoxLoops::loop(EMSCartoonGaussConstraints(dx, coupling), state, diagnostic,
                   interior, disable_simd());
    region_t full, axis;
    const double c = std::cosh(eta);
    for (int j = 0; j < (echo ? N / 2 : 2 * N); ++j)
        for (int i = 0; i < (echo ? N : 4 * N); ++i)
        {
            const double x = (i + 0.5) * dx - (echo ? 0.5 : 2.0),
                         y = (j + 0.5) * dx;
            const double rest = std::hypot(c * x, y);
            if (rest < (echo ? 0.01 : 0.25) || rest > (echo ? 0.5 : 1.5))
                continue;
            const IntVect iv(i, j);
            const auto reference =
                reader.compute_single_bh_vars(x, y, 1.0, 0.0, eta);
            full.accumulate(state, diagnostic, iv, reference);
            if (j < 2)
                axis.accumulate(state, diagnostic, iv, reference);
        }
    if (full.count == 0 || axis.count == 0)
        return 3;
    std::cout << std::setprecision(17);
    std::cout << "N=" << N << " eta=" << eta << " seconds="
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
                     .count()
              << '\n';
    std::cout << "region,count,chi_floor,lapse_floor,field_error,H_L2,H_Linf,"
                 "Mx_L2,Mx_Linf,My_L2,My_Linf,GaussE_L2,GaussE_Linf,GaussB_L2,"
                 "GaussB_Linf\n";
    full.print_summary("full");
    axis.print_summary("axis");
}
