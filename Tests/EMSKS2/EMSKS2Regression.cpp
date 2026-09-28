#include <cmath>
#include "EMSCouplingFunction.hpp"
#include "BoxLoops.hpp"
#include "EMSBHParams.hpp"
#include "CCZ4CartoonVars.hpp"
#include "EMSBH_read.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "ExperimentalGauge.hpp"
#include "GammaCartoonCalculator.hpp"
#include "SetValue.hpp"
#include "SHA256.hpp"
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

// Run this identical source against 34f2ef0 and the EMSKS2 branch. Hash every
// evolved variable at every interior and ghost cell after the old gauge pass.
int main(int argc, char **argv)
{
    if (argc != 6)
    { std::cerr << "usage: EMSKS2Regression legacy.dat reference.trumpet B.trumpet companion.ctt case\n"; return 2; }
    const std::string which = argv[5];
    const bool legacy = which == "legacy";
    const bool binary = which == "ctt";
    const bool echo = which == "echo-B";
    if (!legacy && !binary && !echo && which != "reference") return 2;
    const double dx = binary ? .5 : echo ? 1. / 128 : 1. / 32;
    const int nx = binary ? 128 : echo ? 128 : 128;
    const int ny = binary ? 8 : echo ? 64 : 64;
    EMSBH_params_t params{};
    params.bh_mass = 1;
    params.star_centre = {binary ? 32. : echo ? .5 : 2., 0};
    params.binary = binary;
    params.boosted = binary;
    params.rapidity = binary ? .20273255 : 0;
    params.separation = binary ? 32 : 0;
    params.data_path = legacy ? argv[1] : echo ? argv[3] : argv[2];
    params.ctt_data_path = binary ? argv[4] : "";
    CouplingFunction::params_t coupling{12.566370614359172, 0, 0,
                                        echo ? -.6 : -20.};
    const Box interior(IntVect(0, 0), IntVect(nx - 1, ny - 1));
    Box ghost = interior;
    ghost.grow(5);
    FArrayBox state(ghost, NUM_VARS);
    if (legacy)
    {
        EMSBH_read reader(params, coupling, 1., dx, 0);
        reader.compute_1d_solution();
        BoxLoops::loop(make_compute_pack(SetValue(0.), reader), state, state,
                       disable_simd());
    }
    else
    {
        EMSBH_trumpet_read reader(params, coupling, 1., dx, 0);
        reader.compute_1d_solution();
        BoxLoops::loop(make_compute_pack(SetValue(0.), reader), state, state,
                       disable_simd());
    }
    Box gamma = interior;
    gamma.grow(2);
    BoxLoops::loop(GammaCartoonCalculator(dx), state, state, gamma,
                   disable_simd());
    ExperimentalGauge::params_t gauge{};
    gauge.shift_Gamma_coeff = .75;
    gauge.eta = 1.;
    BoxLoops::loop(ExperimentalGauge(gauge), state, state, interior,
                   disable_simd());
    std::string bytes;
    bytes.reserve(static_cast<std::size_t>(ghost.numPts()) * NUM_VARS * 8);
    for (int j = ghost.smallEnd(1); j <= ghost.bigEnd(1); ++j)
        for (int i = ghost.smallEnd(0); i <= ghost.bigEnd(0); ++i)
            for (int k = 0; k < NUM_VARS; ++k)
            {
                const double value = state(IntVect(i, j), k);
                std::uint64_t bits;
                std::memcpy(&bits, &value, 8);
                for (int byte = 0; byte < 8; ++byte)
                    bytes.push_back(static_cast<char>(bits >> (8 * byte)));
            }
    std::cout << which << ',' << ghost.numPts() << ',' << NUM_VARS << ','
              << SHA256::digest(bytes) << '\n';
}
