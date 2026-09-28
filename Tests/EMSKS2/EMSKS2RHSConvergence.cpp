#include "BoxLoops.hpp"
#include "EMSCouplingFunction.hpp"
#include "CCZ4Cartoon.hpp"
#include "EMSBH_ks2_read.hpp"
#include "EMSKS2ReferenceCache.hpp"
#include "GammaCartoonCalculator.hpp"
#include "ReferenceStationaryGauge.hpp"
#include "SetValue.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sys/resource.h>

static long peak_rss()
{
    struct rusage r{};
    getrusage(RUSAGE_SELF, &r);
    return r.ru_maxrss;
}

struct Norm
{
    double square = 0., maximum = 0.;
    void add(double x)
    {
        square += x * x;
        maximum = std::max(maximum, std::abs(x));
    }
};

struct Mask
{
    std::array<Norm, 14> norms{};
    int count = 0;
};

static const char *group_names[] = {
    "chi", "gamma_tilde", "K", "A_tilde", "Theta", "Gamma_tilde",
    "alpha", "beta", "phi", "Pi", "E", "B_magnetic", "Xi", "Lambda"};

static const std::array<std::array<int, 5>, 14> group_vars = {{
    {{c_chi, -1, -1, -1, -1}},
    {{c_h11, c_h12, c_h22, c_hww, -1}},
    {{c_K, -1, -1, -1, -1}},
    {{c_A11, c_A12, c_A22, c_Aww, -1}},
    {{c_Theta, -1, -1, -1, -1}},
    {{c_Gamma1, c_Gamma2, -1, -1, -1}},
    {{c_lapse, -1, -1, -1, -1}},
    {{c_shift1, c_shift2, -1, -1, -1}},
    {{c_phi, -1, -1, -1, -1}},
    {{c_Pi, -1, -1, -1, -1}},
    {{c_Ex, c_Ey, c_Ez, -1, -1}},
    {{c_Bx, c_By, c_Bz, -1, -1}},
    {{c_Xi, -1, -1, -1, -1}},
    {{c_Lambda, -1, -1, -1, -1}}
}};

static void check_reference(const FArrayBox &ref, const EMSKS2Profile &p,
                            const Box &box, double dx, double cx)
{
    const int stride = std::max(1, static_cast<int>(box.numPts() / 1024));
    int n = 0;
    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
        if (n++ % stride != 0) continue;
        const IntVect iv = bit();
        const auto o = p.object((iv[0] + .5) * dx - cx, (iv[1] + .5) * dx);
        const double expected[6] = {o.s.alpha, o.s.q_star, o.s.K,
                                    o.beta[0], o.beta[1], o.s.taper};
        for (int k = 0; k < 6; ++k)
            if (ref(iv, k) != expected[k])
                throw std::runtime_error("reference cache differs from direct file evaluation");
    }
}

int main(int argc, char **argv)
{
    if (argc != 4 && argc != 5 && argc != 6)
    {
        std::cerr << "usage: EMSKS2RHSConvergence member M/h profile.ks2 [fine-floor-file|- [join-profile.csv]]\n";
        return 2;
    }
    const auto start = std::chrono::steady_clock::now();
    auto stage = [&](const char *name, std::size_t done, std::size_t total)
    {
        std::cerr << "stage=" << name << " member=" << argv[1] << " N=" << argv[2]
                  << " elapsed=" << std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - start).count()
                  << " peak_RSS_bytes=" << peak_rss() << " items=" << done
                  << '/' << total << '\n';
    };
    const std::string member = argv[1];
    const int N = std::stoi(argv[2]);
    double floor = 1e-8;
    if (argc >= 5 && std::string(argv[4]) != "-")
    {
        std::ifstream input(argv[4]);
        std::string key, equals;
        double value;
        bool lapse = false, chi = false;
        while (input >> key >> equals >> value)
        {
            if (equals != "=" || value != 1e-9)
                throw std::runtime_error("invalid fine-grid floor file");
            if (key == "min_lapse" && !lapse) lapse = true;
            else if (key == "min_chi" && !chi) chi = true;
            else throw std::runtime_error("invalid fine-grid floor key");
        }
        if (!input.eof() || !lapse || !chi)
            throw std::runtime_error("incomplete fine-grid floor file");
        floor = 1e-9;
    }
    EMSKS2Profile profile;
    profile.load(argv[3]);
    const double M = profile.get("M"), rh = profile.get("r_h"), dx = M / N;
    const int half = static_cast<int>(std::ceil(2.5 * rh / dx));
    const int nx = 2 * half, ny = half;
    const Box interior(IntVect(0, 0), IntVect(nx - 1, ny - 1));
    Box ghost = interior;
    ghost.grow(5);
    FArrayBox state(ghost, NUM_VARS);
    EMSBH_params_t parameters{};
    parameters.bh_mass = M;
    parameters.data_path = argv[3];
    parameters.star_centre = {half * dx, 0.};
    CouplingFunction::params_t coupling{profile.get("ems_alpha"),
        profile.get("ems_f0"), profile.get("ems_f1"), profile.get("ems_f2")};
    EMSBH_ks2_read reader(parameters, coupling, dx, floor, floor);
    stage("allocate", 0, ghost.numPts());
    BoxLoops::loop(make_compute_pack(SetValue(0.), reader), state, state,
                   disable_simd());
    double min_alpha_margin = std::numeric_limits<double>::infinity();
    double min_chi_margin = min_alpha_margin;
    for (BoxIterator bit(ghost); bit.ok(); ++bit)
    {
        const IntVect iv = bit();
        min_alpha_margin = std::min(min_alpha_margin, state(iv, c_lapse) / floor);
        min_chi_margin = std::min(min_chi_margin, state(iv, c_chi) / floor);
    }
    if (min_alpha_margin < 100 || min_chi_margin < 100 ||
        reader.chi_floor_count() || reader.lapse_floor_count())
        throw std::runtime_error("floor margin or activation failed");
    stage("setter", ghost.numPts(), ghost.numPts());
    Box gamma = interior;
    gamma.grow(2);
    BoxLoops::loop(GammaCartoonCalculator(dx), state, state, gamma,
                   disable_simd());
    FArrayBox ref(ghost, ReferenceStationaryGauge::num_reference);
    EMSKS2ReferenceCache::fill(ref, profile, dx, parameters.star_centre);
    stage("cache_fill", ghost.numPts(), ghost.numPts());
    check_reference(ref, profile, ghost, dx, half * dx);
    // A new patch and a fresh reader model regrid and restart cache rebuilds.
    const Box new_box(IntVect(half + 7, 3), IntVect(half + 15, 11));
    FArrayBox new_ref(new_box, ReferenceStationaryGauge::num_reference);
    EMSKS2Profile reloaded;
    reloaded.load(argv[3]);
    EMSKS2ReferenceCache::fill(new_ref, reloaded, dx, parameters.star_centre);
    check_reference(new_ref, profile, new_box, dx, half * dx);
    stage("reference", ghost.numPts(), ghost.numPts());

    for (int f = 0; f < 2; ++f)
    {
        ReferenceStationaryGauge::params_t gp;
        gp.reference = &ref;
        gp.mass = M;
        gp.onepluslog = f != 0;
        gp.onepluslog_n = 2.;
        ReferenceStationaryGauge gauge(gp);
        const IntVect iv(half + 7, 4);
        CCZ4CartoonVars::VarsWithGauge<double> vars{}, rhs{}, advec{};
        vars.lapse = ref(iv, ReferenceStationaryGauge::alpha_star);
        vars.K = ref(iv, ReferenceStationaryGauge::K_star);
        vars.Theta = 0.;
        advec.lapse = ref(iv, ReferenceStationaryGauge::q_star) * vars.lapse;
        gauge.set_lapse_rhs(rhs, vars, advec, iv);
        if (std::abs(rhs.lapse) > 1e-15 * std::max(1., vars.lapse))
            throw std::runtime_error("reference lapse fixed point failed");
    }

    const double ra = profile.rho_a(), rm = profile.get("r_m");
    const char *mask_names[] = {"join", "KS_collar", "exterior", "far"};
    std::cout << std::setprecision(17)
              << "member,N,KO,mask,group,count,L2,Linf,chi_floor,lapse_floor,min_alpha_margin,min_chi_margin,seconds\n";
    for (int ko = 0; ko < 2; ++ko)
    {
        FArrayBox rhs(interior, NUM_VARS);
        CCZ4_params_t<ReferenceStationaryGauge::params_t> cp{};
        cp.kappa1 = .1; cp.kappa2 = 0.; cp.kappa3 = 1.; cp.covariantZ4 = true;
        cp.reference = &ref; cp.mass = M; cp.onepluslog = false;
        CCZ4Cartoon<ReferenceStationaryGauge, FourthOrderDerivatives,
                    CouplingFunction> op(cp, dx, ko ? 1. : 0.,
                                          CouplingFunction(coupling), 1., 0);
        BoxLoops::loop(op, state, rhs, interior, disable_simd());
        stage(ko ? "rhs_KO_on" : "rhs_KO_off", interior.numPts(), interior.numPts());
        std::array<Mask, 4> masks{};
        std::array<std::array<Norm, 5>, 40> join_profile{};
        std::array<int, 40> join_counts{};
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const double x = (i + .5 - half) * dx, y = (j + .5) * dx;
                const double rho = std::hypot(x, y);
                if (rho - 4 * dx < ra / 2 || rho > 2.3 * rh) continue;
                const int m = rho < ra ? -1 : rho < rm ? 0 :
                              rho < rh ? 1 : rho < 2 * rh ? 2 : 3;
                if (m < 0) continue;
                auto &mask = masks[m];
                const IntVect iv(i, j);
                for (int k = 0; k < 14; ++k)
                {
                    double magnitude2 = 0.;
                    for (int c : group_vars[k])
                        if (c >= 0)
                        {
                            const double v = rhs(iv, c);
                            if (!std::isfinite(v))
                                throw std::runtime_error("nonfinite RHS on mask");
                            magnitude2 += v * v;
                        }
                    mask.norms[k].add(std::sqrt(magnitude2));
                    if (ko && m == 0 && argc == 6 &&
                        (k == 2 || k == 3 || k == 5 || k == 4 || k == 6))
                    {
                        const int bin = std::min(39, static_cast<int>(
                            40 * (rho - ra) / (rm - ra)));
                        const int group = k == 2 ? 0 : k == 3 ? 1 :
                                          k == 5 ? 2 : k == 4 ? 3 : 4;
                        join_profile[bin][group].add(std::sqrt(magnitude2) /
                                                       std::pow(dx, 4));
                    }
                }
                if (ko && m == 0 && argc == 6)
                    ++join_counts[std::min(39, static_cast<int>(
                        40 * (rho - ra) / (rm - ra)))];
                if (rhs(iv, c_shift1) != 0 || rhs(iv, c_shift2) != 0 ||
                    rhs(iv, c_B1) != 0 || rhs(iv, c_B2) != 0)
                    throw std::runtime_error("prescribed shift or gauge B has nonzero RHS");
                ++mask.count;
            }
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        for (int m = 0; m < 4; ++m)
        {
            if (!masks[m].count) throw std::runtime_error("empty mask");
            for (int k = 0; k < 14; ++k)
                std::cout << member << ',' << N << ',' << ko << ',' << mask_names[m]
                          << ',' << group_names[k] << ',' << masks[m].count << ','
                          << std::sqrt(masks[m].norms[k].square / masks[m].count)
                          << ',' << masks[m].norms[k].maximum << ','
                          << reader.chi_floor_count() << ',' << reader.lapse_floor_count()
                          << ',' << min_alpha_margin << ',' << min_chi_margin
                          << ',' << seconds << '\n';
        }
        if (ko && argc == 6)
        {
            std::ofstream profile_out(argv[5]);
            if (!profile_out) throw std::runtime_error("cannot write join RHS profile");
            profile_out << "member,N,bin,rho_mid,r_mid,z_mid,count,group,RMS_h4,max_h4\n"
                        << std::setprecision(17);
            const char *groups[] = {"K", "A_tilde", "Gamma_tilde", "Theta", "alpha"};
            for (int bin = 0; bin < 40; ++bin)
            {
                const double rho = ra + (rm - ra) * (bin + .5) / 40;
                const double r = profile.sample(rho).r;
                for (int group = 0; group < 5; ++group)
                    profile_out << member << ',' << N << ',' << bin << ','
                                << rho << ',' << r << ','
                                << (r - profile.get("r_a")) /
                                   (rm - profile.get("r_a")) << ','
                                << join_counts[bin] << ',' << groups[group] << ','
                                << (join_counts[bin] ? std::sqrt(
                                    join_profile[bin][group].square / join_counts[bin]) : 0.)
                                << ',' << join_profile[bin][group].maximum << '\n';
            }
        }
        const IntVect first(half, 0);
        std::cerr << "audit member=" << member << " N=" << N << " KO=" << ko
                  << " rho=" << dx / std::sqrt(2.)
                  << " chi_rhs=" << rhs(first, c_chi)
                  << " K_rhs=" << rhs(first, c_K)
                  << " alpha_rhs=" << rhs(first, c_lapse) << '\n';
    }
    stage("report", interior.numPts() * 2, interior.numPts() * 2);
}
