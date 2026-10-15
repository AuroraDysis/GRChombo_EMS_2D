// Isolated t=0 hierarchy; no evolution RHS exists in this executable.
#include <cmath>

#include "EMSCouplingFunction.hpp"

#include "BoxIterator.H"
#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ConstraintsCartoon.hpp"
#include "DefaultLevelFactory.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "GammaCartoonCalculator.hpp"
#include "SetValue.hpp"
#include "SetupFunctions.hpp"
#include <chrono>
#include <iomanip>

struct Norms
{
    long count = 0, floored = 0;
    double volume = 0, staging = 0;
    std::array<double, 4> sum{}, max{};
    void print(const char *name) const
    {
        if (!count || floored || staging > 5e-13)
            throw std::runtime_error("empty/floored mask or staging failure");
        std::cout << "NORM," << name << ',' << count << ',' << floored << ','
                  << staging;
        for (int k = 0; k < 4; ++k)
            std::cout << ',' << std::sqrt(sum[k] / volume) << ',' << max[k];
        std::cout << '\n';
    }
};
class CTT2Level : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void initialData() override
    {
        EMSBH_trumpet_read setter(m_p.emsbh_params,
                                  m_p.coupling_function_params, m_p.m_G_Newton,
                                  m_dx, 0);
        setter.compute_1d_solution();
        BoxLoops::loop(make_compute_pack(SetValue(0.), setter), m_state_new,
                       m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
    }
    void specificEvalRHS(GRLevelData &, GRLevelData &, double) override
    {
        throw std::runtime_error("EMSCTT2 hierarchy must never advance");
    }
    void computeTaggingCriterion(FArrayBox &tags, const FArrayBox &) override
    {
        // Fixed geometric refinement across the resolution ladder, at t=0 only.
        for (BoxIterator bit(tags.box()); bit.ok(); ++bit)
        {
            const auto iv = bit();
            const double x = (iv[0] + .5) * m_dx -
                             m_p.emsbh_params.star_centre[0],
                         y = (iv[1] + .5) * m_dx;
            tags(iv, 0) = std::hypot(x, y) < .4 * m_p.L ? 100. : 0.;
        }
    }
    void measure(const std::vector<Box> &covered, std::array<Norms, 8> &norms,
                 double radius_scale, bool analytic_ghosts,
                 const std::array<double, 4> &interfaces, double half_width)
    {
        // Restore analytic t=0 values, then exercise the production AMR ghost
        // fill.
        initialData();
        if (!analytic_ghosts)
            fillAllGhosts();
        for (auto dit = m_state_new.dataIterator(); dit.ok(); ++dit)
        {
            Box region = m_grids[dit];
            region.grow(m_num_ghosts >= 4 ? 2 : 0);
            BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new[dit],
                           m_state_new[dit], region, disable_simd());
        }
        if (m_num_ghosts < 4)
            fillAllGhosts(VariableType::evolution, Interval(c_Gamma1, c_Gamma2));
        BoxLoops::loop(Constraints<CouplingFunction>(
                           m_dx, CouplingFunction(m_p.coupling_function_params),
                           m_p.m_G_Newton),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS,
                       disable_simd());
        BoxLoops::loop(
            EMSCartoonGaussConstraints(m_dx, m_p.coupling_function_params),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS,
            disable_simd());
        EMSBH_trumpet_read setter(m_p.emsbh_params,
                                  m_p.coupling_function_params, m_p.m_G_Newton,
                                  m_dx, 0);
        setter.compute_1d_solution();
        const double w = m_dx * m_dx;
        auto dit = m_state_new.dataIterator();
        long valid = 0, skipped = 0;
        for (dit.begin(); dit.ok(); ++dit)
        {
            const auto &a = m_state_new[dit];
            const auto &d = m_state_diagnostics[dit];
            for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit)
            {
                const auto iv = bit();
                if (std::any_of(covered.begin(), covered.end(),
                                [&](const Box &b) { return b.contains(iv); }))
                {
                    ++skipped;
                    continue;
                }
                ++valid;
                const double x = (iv[0] + .5) * m_dx -
                                 m_p.emsbh_params.star_centre[0],
                             y = (iv[1] + .5) * m_dx, r = std::hypot(x, y);
                const double R = r / radius_scale, Y = y / radius_scale;
                const bool masks[] = {R >= 4 && R <= 8, R >= 25 && R <= 55,
                                      R >= 68 && R <= 76,
                                      R >= 25 && R <= 55 && Y <= 2,
                                      std::abs(R-interfaces[0]) <= half_width,
                                      std::abs(R-interfaces[1]) <= half_width,
                                      std::abs(R-interfaces[2]) <= half_width,
                                      std::abs(R-interfaces[3]) <= half_width};
                if (!std::any_of(std::begin(masks), std::end(masks),
                                 [](bool v) { return v; }))
                    continue;
                const auto v = setter.compute_single_bh_vars(x, y, 1., 0., 0.);
                const std::pair<int, double> fields[] = {
                    {c_chi, v.chi},         {c_h11, v.h[0][0]},
                    {c_h12, v.h[0][1]},     {c_h22, v.h[1][1]},
                    {c_hww, v.hww},         {c_K, v.K},
                    {c_A11, v.A[0][0]},     {c_A12, v.A[0][1]},
                    {c_A22, v.A[1][1]},     {c_Aww, v.Aww},
                    {c_phi, v.phi},         {c_Pi, v.Pi},
                    {c_Ex, v.Ex},           {c_Ey, v.Ey},
                    {c_Ez, v.Ez},           {c_lapse, v.lapse},
                    {c_shift1, v.shift[0]}, {c_shift2, v.shift[1]}};
                double err = 0;
                for (const auto &f : fields)
                    err = std::max(err, std::abs(a(iv, f.first) - f.second) /
                                            (1 + std::abs(f.second)));
                const double values[] = {
                    d(iv, c_Ham), std::hypot(d(iv, c_Mom1), d(iv, c_Mom2)),
                    d(iv, c_GaussE), d(iv, c_GaussB)};
                for (int j = 0; j < 8; ++j)
                    if (masks[j])
                    {
                        auto &n = norms[j];
                        ++n.count;
                        n.volume += w;
                        n.staging = std::max(n.staging, err);
                        n.floored +=
                            a(iv, c_chi) < 1e-6 || a(iv, c_lapse) < 1e-12;
                        for (int k = 0; k < 4; ++k)
                        {
                            if (!std::isfinite(values[k]))
                                throw std::runtime_error(
                                    "nonfinite constraints");
                            n.sum[k] += w * values[k] * values[k];
                            n.max[k] = std::max(n.max[k], std::abs(values[k]));
                        }
                    }
            }
        }
        for (LayoutIterator it=m_grids.layoutIterator(); it.ok(); ++it)
        {
            const Box &b=m_grids[it];
            std::cout<<"BOX,"<<m_level<<','<<b.smallEnd(0)<<','<<b.smallEnd(1)
                     <<','<<b.bigEnd(0)<<','<<b.bigEnd(1)<<'\n';
        }
        std::cout << "LEVEL," << m_level << ',' << m_dx << ',' << valid << ','
                  << skipped << '\n';
    }
};
int run(int argc, char **argv)
{
    GRParmParse pp(argc - 2, argv + 2, nullptr, argv[1]);
    SimulationParameters p(pp);
    bool check_params_only;
    pp.load("ctt2_test_check_params_only", check_params_only, false);
    if (check_params_only)
    {
        std::cout << "startup parameters accepted\n";
        return 0;
    }
    if (p.restart_from_checkpoint || p.max_steps != 0 || p.stop_time != 0 ||
        p.max_level != 1)
        throw std::runtime_error("requires a fresh two-level t=0 hierarchy");
    double radius_scale;
    pp.load("ctt2_test_radius_scale", radius_scale, 1.);
    if (!std::isfinite(radius_scale) || radius_scale <= 0.)
        throw std::runtime_error("positive finite test radius scale required");
    bool analytic_ghosts;
    pp.load("ctt2_test_analytic_ghosts", analytic_ghosts, false);
    std::cout<<"CONFIG,"<<p.amr_transfer<<','<<p.num_ghosts<<','<<analytic_ghosts<<'\n';
    GRAMR amr;
    DefaultLevelFactory<CTT2Level> factory(amr, p);
    setupAMRObject(amr, factory);
    std::array<Norms, 8> norms;
    std::array<double, 4> interfaces{-1., -1., -1., -1.};
    const bool interface_masks = pp.contains("ctt2_test_interfaces");
    double half_width = 0.;
    if (interface_masks)
    {
        pp.load("ctt2_test_interfaces", interfaces);
        pp.load("ctt2_test_interface_half_width", half_width);
        if (!(half_width > 0.) || !std::isfinite(half_width) ||
            !std::all_of(interfaces.begin(), interfaces.end(),
                         [](double r) { return r > 0. && std::isfinite(r); }))
            throw std::runtime_error("invalid interface masks");
    }
    const auto levels = amr.get_gramrlevels();
    if (levels.size() != 2)
        throw std::runtime_error("two initialized AMR levels required");
    for (size_t l = 0; l < levels.size(); ++l)
    {
        std::vector<Box> covered;
        if (l + 1 < levels.size())
        {
            const auto &layout =
                levels[l + 1]->getLevelData().disjointBoxLayout();
            for (LayoutIterator it = layout.layoutIterator(); it.ok(); ++it)
            {
                Box b = layout[it];
                b.coarsen(2);
                covered.push_back(b);
            }
        }
        dynamic_cast<CTT2Level &>(*levels[l])
            .measure(covered, norms, radius_scale, analytic_ghosts, interfaces, half_width);
    }
    const char *names[] = {"inner", "pulse", "outer", "pulse_axis",
                           "interface_1", "interface_2", "interface_3", "interface_4"};
    for (int j = 0; j < (interface_masks ? 8 : 4); ++j)
        norms[j].print(names[j]);
    return 0;
}
int main(int argc, char **argv)
{
    mainSetup(argc, argv);
    int status = 0;
    const auto start = std::chrono::steady_clock::now();
    std::cout << std::setprecision(17);
    try
    {
        status = run(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        status = 2;
    }
    std::cout << "seconds="
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
                     .count()
              << '\n';
    mainFinalize();
    return status;
}
