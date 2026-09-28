#include <cmath>
#include "EMSCouplingFunction.hpp"
#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ConstraintsCartoon.hpp"
#include "EMSBH_ks2_read.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include "GammaCartoonCalculator.hpp"
#include "SetValue.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <sys/resource.h>
#include <vector>

static long peak_rss()
{ struct rusage r{}; getrusage(RUSAGE_SELF, &r); return r.ru_maxrss; }

struct Norm
{
    double square = 0, maximum = 0;
    void add(double x) { square += x * x; maximum = std::max(maximum, std::abs(x)); }
};
struct Mask
{
    std::array<Norm, 6> norm{};
    int count = 0;
    void add(const FArrayBox &state, const FArrayBox &diag,
             const IntVect &iv, const EMSKS2Profile::Cartesian &analytic)
    {
        const int vars[5] = {c_Ham, c_Mom1, c_Mom2, c_GaussE, c_GaussB};
        for (int k = 0; k < 5; ++k) norm[k].add(diag(iv, vars[k]));
        norm[5].add(std::hypot(state(iv, c_Gamma1) - analytic.Gamma[0],
                               state(iv, c_Gamma2) - analytic.Gamma[1]));
        ++count;
    }
};

// Test-only decomposition, following Constraints::constraint_equations.
class ConstraintTerms : public Constraints<CouplingFunction>
{
  public:
    ConstraintTerms(double dx, const CouplingFunction::params_t &coupling)
        : Constraints(dx, CouplingFunction(coupling), 1.) {}

    void compute(Cell<double> cell) const
    {
        using namespace TensorAlgebra;
        const auto vars = cell.load_vars<Vars>();
        const auto d1 = m_deriv.diff1<Vars>(cell);
        const auto d2 = m_deriv.diff2<Diff2Vars>(cell);
        const double y = Coordinates<double>(cell, m_dx).y;
        const int dI = CH_SPACEDIM - 1, nS = GR_SPACEDIM - CH_SPACEDIM;
        const double chi_reg = std::max(1e-6, vars.chi);
        const auto h_UU = compute_inverse_sym(vars.h);
        const double h_UU_ww = 1. / vars.hww;
        const auto chris = compute_christoffel(d1.h, h_UU);
        const auto ricci = CCZ4CartoonGeometry::compute_ricci(
            vars, d1, d2, h_UU, h_UU_ww, chris, y);
        const auto A_UU = raise_all(vars.A, h_UU);
        const double tr_A2 = compute_trace(vars.A, A_UU) +
            nS * vars.Aww * vars.Aww * h_UU_ww * h_UU_ww;
        const auto em = compute_EMS_EM_tensor(vars, d1, h_UU, h_UU_ww,
                                               chris, nS, m_coupling);
        const double H_kin = (GR_SPACEDIM - 1.) * vars.K * vars.K /
                             GR_SPACEDIM;
        const double H_matter = -16. * M_PI * em.rho;
        cell.store_vars(std::abs(ricci.scalar) + std::abs(H_kin) +
                        std::abs(tr_A2) + std::abs(H_matter), 0);
        cell.store_vars(ricci.scalar + H_kin - tr_A2 + H_matter, 3);

        Tensor<1, double> chris_ww, reg_07;
        FOR(i)
        {
            chris_ww[i] = (delta(i, dI) - h_UU[i][dI] * vars.hww) / y;
            FOR(j) chris_ww[i] -= .5 * h_UU[i][j] * d1.hww[j];
            reg_07[i] = (vars.A[i][dI] - delta(i, dI) * vars.Aww) / y;
        }
        Tensor<2, double> covd_A[CH_SPACEDIM];
        FOR(i, j, k)
        {
            covd_A[i][j][k] = d1.A[j][k][i];
            FOR(l) covd_A[i][j][k] -= chris.ULL[l][i][j] * vars.A[l][k] +
                                      chris.ULL[l][i][k] * vars.A[l][j];
        }
        FOR(i)
        {
            double first = -(GR_SPACEDIM - 1.) * d1.K[i] / GR_SPACEDIM;
            double cartoon = nS * h_UU_ww *
                (reg_07[i] - .5 * h_UU_ww * vars.Aww * d1.hww[i]);
            double sum = first + cartoon;
            double den = std::abs(first) + std::abs(cartoon);
            FOR(j)
            {
                const double ww = -nS * h_UU_ww * chris_ww[j] * vars.A[i][j];
                sum += ww;
                den += std::abs(ww);
                FOR(k)
                {
                    const double div = h_UU[j][k] *
                        (covd_A[k][j][i] - GR_SPACEDIM * vars.A[i][j] *
                         d1.chi[k] / (2 * chi_reg));
                    sum += div;
                    den += std::abs(div);
                }
            }
            const double matter = -8. * M_PI * em.Si[i];
            cell.store_vars(den + std::abs(matter), i + 1);
            cell.store_vars(sum + matter, i + 4);
        }
    }
};

struct Bin
{
    double H2 = 0, M2 = 0, HD2 = 0, MD2 = 0;
    double Mx2 = 0, My2 = 0, MxD2 = 0, MyD2 = 0;
    double Hmax = 0, Mmax = 0, rsum = 0, zsum = 0;
    int count = 0;
};
int main(int argc, char **argv)
{
    if (argc < 4 || argc > 6) {
        std::cerr << "usage: EMSKS2GridConvergence member M/h profile.ks2 [profile.csv | fine-params.txt profile.csv]\n";
        return 2;
    }
    const auto start = std::chrono::steady_clock::now();
    const std::string member = argv[1];
    const int N = std::stoi(argv[2]);
    double min_lapse = 1e-8, min_chi = 1e-8;
    if (argc == 6)
    {
        std::ifstream file(argv[4]);
        std::string name, equals;
        bool got_lapse = false, got_chi = false;
        double value;
        while (file >> name >> equals >> value)
        {
            if (equals != "=" || value != 1e-9) throw std::runtime_error("bad fine-grid floors");
            if (name == "min_lapse" && !got_lapse) { min_lapse = value; got_lapse = true; }
            else if (name == "min_chi" && !got_chi) { min_chi = value; got_chi = true; }
            else throw std::runtime_error("bad fine-grid floor key");
        }
        if (!file.eof() || !got_lapse || !got_chi)
            throw std::runtime_error("missing or malformed fine-grid floors");
    }
    const char *profile_path = argc >= 5 ? argv[argc - 1] : nullptr;
    EMSKS2Profile p;
    p.load(argv[3]);
    const double M = p.get("M"), rh = p.get("r_h"), dx = M / N;
    const int half = static_cast<int>(std::ceil(2.5 * rh / dx));
    const int nx = 2 * half, ny = half;
    std::cerr << "stage=allocate member=" << member << " N=" << N
              << " elapsed=0 peak_RSS_bytes=" << peak_rss()
              << " items=0/" << nx * ny << '\n';
    EMSBH_params_t params{};
    params.bh_mass = M;
    params.data_path = argv[3];
    params.star_centre = {half * dx, 0};
    CouplingFunction::params_t coupling{p.get("ems_alpha"), p.get("ems_f0"),
                                        p.get("ems_f1"), p.get("ems_f2")};
    EMSBH_ks2_read reader(params, coupling, dx, min_chi, min_lapse);
    const Box interior(IntVect(0, 0), IntVect(nx - 1, ny - 1));
    Box ghost = interior;
    ghost.grow(5);
    FArrayBox state(ghost, NUM_VARS), diag(interior, NUM_DIAGNOSTIC_VARS);
    BoxLoops::loop(make_compute_pack(SetValue(0.), reader), state, state,
                   disable_simd());
    double alpha_margin = std::numeric_limits<double>::infinity();
    double chi_margin = std::numeric_limits<double>::infinity();
    for (int j = -5; j < ny + 5; ++j)
        for (int i = -5; i < nx + 5; ++i)
        {
            const IntVect iv(i, j);
            alpha_margin = std::min(alpha_margin, state(iv, c_lapse) / min_lapse);
            chi_margin = std::min(chi_margin, state(iv, c_chi) / min_chi);
        }
    if (alpha_margin < 100 || chi_margin < 100 ||
        reader.chi_floor_count() || reader.lapse_floor_count())
        throw std::runtime_error("fine-grid floor margin or count failed");
    std::cerr << "stage=setter member=" << member << " N=" << N
              << " elapsed=" << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start).count()
              << " peak_RSS_bytes=" << peak_rss()
              << " items=" << ghost.numPts() << "/" << ghost.numPts()
              << " min_alpha_margin=" << alpha_margin
              << " min_chi_margin=" << chi_margin << '\n';
    Box gamma = interior;
    gamma.grow(2);
    BoxLoops::loop(GammaCartoonCalculator(dx), state, state, gamma,
                   disable_simd());
    BoxLoops::loop(Constraints<CouplingFunction>(dx, CouplingFunction(coupling), 1.),
                   state, diag, interior, disable_simd());
    BoxLoops::loop(EMSCartoonGaussConstraints(dx, coupling), state, diag,
                   interior, disable_simd());
    std::cerr << "stage=constraints member=" << member << " N=" << N
              << " elapsed=" << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start).count()
              << " peak_RSS_bytes=" << peak_rss()
              << " items=" << nx * ny << "/" << nx * ny << '\n';
    const double ra = p.rho_a(), rm = p.get("r_m");
    Box probe_box = interior;
    FArrayBox *terms = nullptr;
    if (profile_path)
    {
        const int extent = static_cast<int>(std::ceil(rm / dx)) + 1;
        probe_box = Box(IntVect(std::max(0, half - extent), 0),
                        IntVect(std::min(nx - 1, half + extent),
                                std::min(ny - 1, extent)));
        terms = new FArrayBox(probe_box, 6);
        BoxLoops::loop(ConstraintTerms(dx, coupling), state, *terms,
                       probe_box, disable_simd());
        std::cerr << "stage=terms member=" << member << " N=" << N
                  << " elapsed=" << std::chrono::duration<double>(
                      std::chrono::steady_clock::now() - start).count()
                  << " peak_RSS_bytes=" << peak_rss() << " items="
                  << probe_box.numPts() << "/" << probe_box.numPts() << '\n';
    }
    Mask masks[4];
    std::array<Bin, 40> bins{};
    const double inner = ra / 2;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const double x = (i + .5 - half) * dx, y = (j + .5) * dx;
            const double rho = std::hypot(x, y);
            if (rho - 4 * dx < inner || rho > 2.3 * rh) continue;
            const int mask = rho < ra ? -1 : rho < rm ? 0 :
                             rho < rh ? 1 : rho < 2 * rh ? 2 : 3;
            if (mask < 0) continue;
            const IntVect iv(i, j);
            const auto analytic = p.object(x, y);
            masks[mask].add(state, diag, iv, analytic);
            if (mask == 0 && terms)
            {
                const int b = std::min(39, int(40 * (rho - ra) / (rm - ra)));
                auto &bin = bins[b];
                const double H = diag(iv, c_Ham), Mx = diag(iv, c_Mom1);
                const double My = diag(iv, c_Mom2), M = std::hypot(Mx, My);
                const double HD = (*terms)(iv, 0), MxD = (*terms)(iv, 1);
                const double MyD = (*terms)(iv, 2), MD = std::hypot(MxD, MyD);
                if (std::abs(H - (*terms)(iv, 3)) > 1e-9 * std::max(1., HD) ||
                    std::abs(Mx - (*terms)(iv, 4)) > 1e-9 * std::max(1., MxD) ||
                    std::abs(My - (*terms)(iv, 5)) > 1e-9 * std::max(1., MyD))
                    throw std::runtime_error("constraint term reconstruction differs");
                bin.H2 += H * H; bin.M2 += M * M;
                bin.Mx2 += Mx * Mx; bin.My2 += My * My;
                bin.HD2 += HD * HD; bin.MD2 += MD * MD;
                bin.MxD2 += MxD * MxD; bin.MyD2 += MyD * MyD;
                bin.Hmax = std::max(bin.Hmax, std::abs(H));
                bin.Mmax = std::max(bin.Mmax, M);
                bin.rsum += analytic.s.r;
                bin.zsum += (analytic.s.r - p.get("r_a")) /
                            (p.get("r_m") - p.get("r_a"));
                ++bin.count;
            }
        }
    if (terms)
    {
        std::ofstream out(profile_path);
        if (!out) throw std::runtime_error("cannot write radial profile");
        out << std::setprecision(17)
            << "member,N,bin,rho_lo,rho_hi,rho_mid,r_mean,z_mean,count,H_max_h4,H_rms_h4,M_max_h4,M_rms_h4,H_ratio_rms,M_ratio_rms,Mx_ratio_rms,My_ratio_rms\n";
        const double h4 = std::pow(dx, 4);
        for (int b = 0; b < 40; ++b)
        {
            const auto &v = bins[b];
            if (!v.count) throw std::runtime_error("empty join profile bin");
            const double lo = ra + (rm - ra) * b / 40.;
            const double hi = ra + (rm - ra) * (b + 1) / 40.;
            out << member << ',' << N << ',' << b << ',' << lo << ','
                << hi << ',' << (lo + hi) / 2 << ',' << v.rsum / v.count
                << ',' << v.zsum / v.count << ',' << v.count << ','
                << v.Hmax / h4 << ',' << std::sqrt(v.H2 / v.count) / h4
                << ',' << v.Mmax / h4 << ',' << std::sqrt(v.M2 / v.count) / h4
                << ',' << std::sqrt(v.H2 / v.HD2)
                << ',' << std::sqrt(v.M2 / v.MD2)
                << ',' << std::sqrt(v.Mx2 / v.MxD2)
                << ',' << std::sqrt(v.My2 / v.MyD2) << '\n';
        }
        delete terms;
    }
    std::cout << std::setprecision(17);
    std::cout << "member,N,mask,residual,count,L2,Linf,chi_floor,lapse_floor,seconds\n";
    const char *names[4] = {"join", "KS_collar", "exterior", "far"};
    const char *residuals[6] = {"H", "Mx", "My", "GaussE", "GaussB", "Gamma_error"};
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    for (int m = 0; m < 4; ++m)
    {
        if (!masks[m].count) throw std::runtime_error("empty KS2 mask");
        for (int k = 0; k < 6; ++k)
            std::cout << member << ',' << N << ',' << names[m] << ','
                      << residuals[k] << ',' << masks[m].count << ','
                      << std::sqrt(masks[m].norm[k].square / masks[m].count)
                      << ',' << masks[m].norm[k].maximum << ','
                      << reader.chi_floor_count() << ','
                      << reader.lapse_floor_count() << ',' << seconds << '\n';
    }
    const double x = dx / 2, y = dx / 2, rho = std::hypot(x, y);
    const auto nearest = p.object(x, y);
    const double F = std::exp(-2 * coupling.alpha *
        (coupling.f0 + coupling.f1 * nearest.s.phi +
         coupling.f2 * nearest.s.phi * nearest.s.phi));
    const double E2 = nearest.s.E_rho * nearest.s.E_rho /
                      (nearest.s.lambda * std::pow(nearest.s.r / rho, 2));
    const auto core = p.source(p.get("r_0"), true);
    const double F0 = std::exp(-2 * coupling.alpha *
        (coupling.f0 + coupling.f1 * core[1] + coupling.f2 * core[1] * core[1]));
    const double q = p.get("Q") / std::sqrt(8 * std::acos(-1.));
    const double E2limit = q * q /
        (F0 * F0 * std::pow(p.get("r_0"), 4));
    std::cout << "audit," << member << ',' << N << ",rho," << rho
              << ",alpha," << nearest.s.alpha << ",alpha_over_rho_nu,"
              << nearest.s.alpha / std::pow(rho, p.get("nu"))
              << ",chi," << nearest.s.chi << ",chi_power_ratio,"
              << nearest.s.chi * std::pow(p.get("r_0") / rho, 2)
              << ",E2," << E2 << ",E2_limit_ratio," << E2 / E2limit << '\n';
    std::cerr << "stage=report member=" << member << " N=" << N
              << " elapsed=" << seconds << " peak_RSS_bytes="
              << peak_rss() << " items=" << nx * ny << "/" << nx * ny << '\n';
}
