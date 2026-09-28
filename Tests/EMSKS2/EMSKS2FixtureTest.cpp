#include "EMSKS2Profile.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static std::vector<std::string> split(const std::string &line)
{
    std::vector<std::string> out;
    std::istringstream in(line);
    for (std::string s; std::getline(in, s, '\t');) out.push_back(s);
    return out;
}
static double field(const std::string &key, const EMSKS2Profile::Cartesian &o,
                    double x, double y, double z)
{
    const auto &s = o.s;
    if (key == "rho") return s.rho;
    if (key == "r") return s.r;
    if (key == "x") return x;
    if (key == "y") return y;
    if (key == "z") return z;
    if (key == "N") return s.N;
    if (key == "H") return s.H;
    if (key == "sigma") return s.sigma;
    if (key == "g") return s.g;
    if (key == "lambda") return s.lambda;
    if (key == "w") return s.w;
    if (key == "alpha" || key == "alpha_star") return s.alpha;
    if (key == "alpha_r") return s.alpha_r;
    if (key == "S") return s.S;
    if (key == "S_r") return s.S_r;
    if (key == "k_r") return s.k_r;
    if (key == "k_t") return s.k_t;
    if (key == "K" || key == "K_star") return s.K;
    if (key == "phi") return s.phi;
    if (key == "phi_r") return s.phi_r;
    if (key == "Pi") return s.Pi;
    if (key == "E_rho") return s.E_rho;
    if (key == "rho_a") return s.rho_a;
    if (key == "chi") return s.chi;
    if (key == "gtilde_ww") return o.hww;
    if (key == "Atilde_ww") return o.Aww;
    if (key == "alpha_star_rho") return s.alpha_rho;
    if (key == "beta_star_rho") return s.beta_rho;
    if (key == "q_star") return s.q_star;
    if (key == "ko_taper") return s.taper;
    if (key.rfind("gtilde_", 0) == 0 && key.size() == 9)
        return o.h[3 * (key[7] - '1') + key[8] - '1'];
    if (key.rfind("Atilde_", 0) == 0 && key.size() == 9)
        return o.A[3 * (key[7] - '1') + key[8] - '1'];
    if (key.rfind("Gamma_", 0) == 0) return o.Gamma.at(key[6] - '1');
    if (key.rfind("beta_", 0) == 0) return o.beta.at(key[5] - '1');
    if (key.rfind("E_", 0) == 0) return o.E.at(key[2] - '1');
    throw std::runtime_error("unhandled fixture column " + key);
}
static void compare(const std::string &name, double actual, double expected,
                    double tolerance, double &maximum, int &count)
{
    const double error = std::abs(actual - expected) / (1 + std::abs(expected));
    maximum = std::max(maximum, error);
    ++count;
    if (!(error <= tolerance))
    {
        std::ostringstream out;
        out << std::setprecision(17) << name << ": " << actual << " versus "
            << expected << ", normalized error " << error;
        throw std::runtime_error(out.str());
    }
}
static void test_member(const std::string &dir, const std::string &name)
{
    EMSKS2Profile p;
    p.load(dir + "/" + name + ".ks2");
    double maximum = 0;
    int count = 0;
    std::ifstream table(dir + "/" + name + ".tsv");
    std::string line;
    if (!std::getline(table, line)) throw std::runtime_error("missing fixture table");
    const auto keys = split(line);
    while (std::getline(table, line))
    {
        const auto row = split(line);
        if (row.size() != keys.size()) throw std::runtime_error("bad fixture row");
        double x = 0, y = 0, z = 0;
        for (std::size_t k = 0; k < keys.size(); ++k)
        { if (keys[k] == "x") x = std::stod(row[k]);
          if (keys[k] == "y") y = std::stod(row[k]);
          if (keys[k] == "z") z = std::stod(row[k]); }
        const auto o = p.object(x, y, z);
        for (std::size_t k = 1; k < keys.size(); ++k)
        {
            const bool first = keys[k] == "alpha_r" || keys[k] == "S_r" ||
                keys[k] == "phi_r" || keys[k] == "alpha_star_rho";
            compare(name + "/" + row[0] + "/" + keys[k],
                    field(keys[k], o, x, y, z), std::stod(row[k]),
                    first ? 1e-10 : 1e-12, maximum, count);
        }
    }
    std::ifstream elements(dir + "/" + name + "-elements.tsv");
    if (!std::getline(elements, line)) throw std::runtime_error("missing element table");
    const auto ek = split(line);
    while (std::getline(elements, line))
    {
        const auto row = split(line);
        if (row.size() != ek.size()) throw std::runtime_error("bad element row");
        const bool interior = row[0] == "interior";
        const double u = std::stod(row[1]), r = std::stod(row[2]);
        const auto v = p.source(r, interior), du = p.source(r, interior, 1);
        const double h = p.get("exterior_h"), beta = p.get("exterior_beta");
        const double t = h == 0 ? u : (1 + h) * u / (h + u);
        const double mapped = interior ? u : beta == 0 ? t :
            .5 + std::asinh((2 * t - 1) * std::sinh(beta)) / (2 * beta);
        const double lo = interior ? p.get("interior_lo") : 0;
        const double hi = interior ? p.get("interior_hi") : 1;
        const double z = (2 * mapped - lo - hi) / (hi - lo);
        for (const auto &kv : {std::pair<const char *, double>{"u", p.get("r_h") / r},
                               {"r", p.get("r_h") / u},
                               {"mapped_y", mapped}, {"cheb_z", z}})
        {
            const auto it = std::find(ek.begin(), ek.end(), kv.first);
            compare(name + "/" + row[0] + "/" + kv.first, kv.second,
                    std::stod(row[std::distance(ek.begin(), it)]),
                    1e-12, maximum, count);
        }
        const char *names[5] = {"mu", "phi", "p", "delta", "V"};
        for (int j = 0; j < 5; ++j)
        {
            for (int d = 0; d < 2; ++d)
            {
                const std::string key = std::string(names[j]) + (d ? "_u" : "");
                const auto it = std::find(ek.begin(), ek.end(), key);
                if (it == ek.end()) throw std::runtime_error("missing source column");
                const auto k = std::distance(ek.begin(), it);
                compare(name + "/" + row[0] + "/" + key,
                        d ? du[j] : v[j], std::stod(row[k]),
                        d ? 1e-10 : 1e-12, maximum, count);
            }
        }
    }
    std::cout << name << "," << count << "," << std::setprecision(12)
              << maximum << "," << p.fingerprint() << '\n';
}
int main(int argc, char **argv)
{
    if (argc != 2) { std::cerr << "usage: EMSKS2FixtureTest fixtures-dir\n"; return 2; }
    try
    {
        std::cout << "member,comparisons,max_normalized_error,sha256\n";
        for (const char *name : {"B", "E", "reference"})
            test_member(argv[1], name);
        std::ifstream malformed(std::string(argv[1]) + "/malformed.tsv");
        std::string line;
        std::getline(malformed, line);
        int rejected = 0;
        while (std::getline(malformed, line))
        {
            const auto row = split(line);
            EMSKS2Profile p;
            std::string reason;
            if (row.size() != 2 ||
                p.check_file(std::string(argv[1]) + "/" + row[0], reason))
                throw std::runtime_error("accepted malformed file " + row[0]);
            ++rejected;
            std::cout << "reject," << row[0] << "," << reason << '\n';
        }
        if (rejected != 7) throw std::runtime_error("missing malformed files");
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
