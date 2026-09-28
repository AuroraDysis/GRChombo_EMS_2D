#include "EMSKS2Profile.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    if (argc != 4) return 2;
    EMSKS2Profile profile;
    profile.load(argv[3]);
    const int N = std::stoi(argv[2]);
    const double h = profile.get("M") / N;
    const auto sample = profile.object(h / 2, h / 2).s;
    constexpr double min_lapse = 1e-8, min_chi = 1e-8;
    std::cout << std::setprecision(17) << argv[1] << ',' << N << ','
              << h << ',' << h / std::sqrt(2.) << ',' << sample.alpha << ','
              << sample.alpha / min_lapse << ',' << sample.chi << ','
              << sample.chi / min_chi << ','
              << (sample.alpha >= 100 * min_lapse &&
                  sample.chi >= 100 * min_chi ? "pass" : "reject") << '\n';
}
