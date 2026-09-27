#include "GRAMR.hpp"
#include "SetupFunctions.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "BoxLoops.hpp"
#include "EMSBH_trumpet_read.hpp"
#include "SetValue.hpp"
#include "RHUnion.hpp"
#include "../../../Source/Extraction/EMSRadiationParameters.hpp"


// Like AMRInterpolatorTest: a real hierarchy and production interpolation,
// with initial data only. There is deliberately no evolution RHS.
class EMSIsolatedRHLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void initialData() override
    {
        EMSBH_trumpet_read setter(m_p.emsbh_params, m_p.coupling_function_params,
                                  m_p.m_G_Newton, m_dx, 0);
        setter.compute_1d_solution();
        BoxLoops::loop(make_compute_pack(SetValue(0.), setter), m_state_new,
                       m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
    }
    void specificEvalRHS(GRLevelData &, GRLevelData &, double) override
    { throw std::runtime_error("EMSIsolatedRHTest is t=0 only"); }
    void computeTaggingCriterion(FArrayBox &tags, const FArrayBox &) override
    { tags.setVal(0.); }
};

int run(int argc, char **argv)
{
    GRParmParse pp(argc-2, argv+2, nullptr, argv[1]);
    SimulationParameters p(pp);
    struct PlotAMR : GRAMR
    {
        using AMR::writePlotFile;
        using AMR::writeCheckpointFile;
    };
    PlotAMR amr;
    DefaultLevelFactory<EMSIsolatedRHLevel> factory(amr, p);
    setupAMRObject(amr, factory);
    AMRInterpolator<Lagrange<4>> interp(amr, p.origin, p.dx, p.boundary_params, 0);
    amr.set_interpolator(&interp);
    RHUnion geometry;
    geometry.setup(p.m_RH_num_horizons, p.m_RH_initial_radii,
                   p.m_RH_initial_centre, p.m_RH_num_points, p.m_RH_level,
                   p.m_RH_time_step_freq, p.m_RH_newton_crit,
                   p.m_RH_chase_speeds, p.m_RH_start_times);
    const auto c = p.coupling_function_params;
    geometry.set_coupling_params(c.alpha, c.f0, c.f1, c.f2);
    geometry.set_interpolator(&interp);
    std::string initial_shape;
    pp.load("ems_test_initial_shape", initial_shape, std::string());
    if (!initial_shape.empty())
    {
        if (geometry.m_surfaces.size() != 1)
            throw std::runtime_error("shape restart is a single-surface test control");
        std::ifstream input(initial_shape);
        if (!input) throw std::runtime_error("cannot open initial shape");
        std::string line, last;
        while (std::getline(input, line))
        {
            const auto first = line.find_first_not_of(" \t");
            if (first != std::string::npos && line[first] != '#') last = line;
        }
        std::istringstream values(last);
        double time;
        if (!(values >> time) || time != 0)
            throw std::runtime_error("initial shape must be at t=0");
        auto &s = geometry.m_surfaces.front();
        for (int j=s.m_NG; j<s.m_NG+s.m_n; ++j)
            if (!(values >> s.m_f[j]) || !std::isfinite(s.m_f[j]) || s.m_f[j] <= 0)
                throw std::runtime_error("invalid initial shape");
        std::string extra;
        if (values >> extra) throw std::runtime_error("extra initial shape values");
        s.fill_all_ghosts();
        std::cout << "EMSRH_TEST_INITIAL_SHAPE " << initial_shape << '\n';
    }
    double threshold;
    pp.load("ems_test_expansion_squared", threshold, 1e-7);
    if (pp.contains("ems_rh_expansion_threshold"))
        threshold = ems_rh_expansion_threshold(pp);
    if (!std::isfinite(threshold) || threshold <= 0 || threshold > 1e-7)
        throw std::runtime_error("invalid test expansion threshold");
    geometry.m_thresh_super_low = threshold;
    bool verify_only;
    pp.load("ems_test_verify_only", verify_only, false);
    if (verify_only)
    {
        if (initial_shape.empty()) throw std::runtime_error("verification requires saved shape");
        // Test-only replay: no chase or re-centring after this fresh interpolation.
        geometry.interpolate_fields();
        const auto &s = geometry.m_surfaces.front();
        const double error = s.expansion_error();
        std::cout << std::setprecision(17) << "EMSRH_FRESH " << s.Area() << ' '
                  << s.Q_charge() << ' ' << s.average_Theta_plus() << ' '
                  << s.average_Theta_minus() << ' ' << error << '\n';
        return std::isfinite(error) && error <= threshold ? 0 : 1;
    }
    int max_updates;
    pp.load("ems_test_max_updates", max_updates, 60);
    amr.writePlotFile();
    bool checkpoint;
    pp.load("ems_test_checkpoint", checkpoint, false);
    if (checkpoint) amr.writeCheckpointFile();
    bool done = false;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < max_updates; ++i)
    {
        geometry.update(0.,0);
        done = std::all_of(geometry.m_surfaces.begin(), geometry.m_surfaces.end(),
                         [](const RHSurf &s) { return !s.m_dead && s.m_state == RHSurf::SolverState::FOUND; });
        if (done) break;
        if (std::any_of(geometry.m_surfaces.begin(), geometry.m_surfaces.end(),
                        [](const RHSurf &s) { return s.m_dead; })) break;
    }
    std::cout << "\nEMSRH_T0 " << (done ? "FOUND" : "FAILED")
              << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
              << '\n';
    if (done)
    {
        // Preserve full-precision geometry for testing serialization error only.
        // The author's rh_surf/rh_f outputs above remain untouched.
        geometry.interpolate_fields();
        for (size_t k=0; k<geometry.m_surfaces.size(); ++k)
        {
            const auto &s=geometry.m_surfaces[k];
            std::cout << std::setprecision(17) << "EMSRH_FINAL_FRESH " << k << ' '
                      << s.expansion_error() << '\n';
            std::ofstream f("ems_test_shape"+std::to_string(k)+".dat");
            f << std::setprecision(17) << "0";
            for (int j=s.m_NG; j<s.m_NG+s.m_n; ++j) f << ' ' << s.m_f[j];
            f << '\n';
        }
    }
    return done ? 0 : 1;
}
int main(int argc, char **argv)
{
    mainSetup(argc, argv);
    int status = 0;
    try { status = run(argc, argv); }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; status = 2; }
    mainFinalize();
    return status;
}
