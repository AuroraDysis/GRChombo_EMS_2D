// Frozen production checkpoint; no initial-data setter or evolution path.
#include "GRAMR.hpp"
#include "SetupFunctions.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "RHUnion.hpp"
#include "T7RHS.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include "BoxLoops.hpp"
#include <cstdint>

class FrozenLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void initialData() override { throw std::runtime_error("restart_file required"); }
    void specificEvalRHS(GRLevelData &, GRLevelData &, double) override
    { throw std::runtime_error("offline snapshot must never advance"); }
    void computeTaggingCriterion(FArrayBox &, const FArrayBox &) override
    { throw std::runtime_error("offline snapshot must never regrid"); }
    void t7_prepare() { fillAllEvolutionGhosts(); }
    void t7_dump(const SimulationParameters &p)
    {
        std::ofstream out("maxwell.bin",std::ios::binary);
        if (!out) throw std::runtime_error("T7 Maxwell output open failed");
        out.write("T7MAX001",8);
        const double meta[]={m_time,m_dx,p.center[0]};
        const std::uint32_t header[]={std::uint32_t(m_level),NUM_VARS,T7RHS::components};
        out.write(reinterpret_cast<const char*>(meta),sizeof(meta));
        out.write(reinterpret_cast<const char*>(header),sizeof(header));
        T7RHS kernel(p.ccz4_params,m_dx,p.sigma,CouplingFunction(p.coupling_function_params),p.m_G_Newton,p.formulation);
        for (DataIterator it=m_state_new.dataIterator();it.ok();++it)
        {
            const Box b=m_state_new.disjointBoxLayout()[it()];
            FArrayBox rhs(b,T7RHS::components),gauss(b,NUM_DIAGNOSTIC_VARS);
            BoxLoops::loop(kernel,m_state_new[it()],rhs,b,disable_simd());
            BoxLoops::loop(EMSCartoonGaussConstraints(m_dx,p.coupling_function_params),m_state_new[it()],gauss,b,disable_simd());
            const std::int32_t bounds[]={b.smallEnd(0),b.smallEnd(1),b.bigEnd(0),b.bigEnd(1)};
            out.write(reinterpret_cast<const char*>(bounds),sizeof(bounds));
            for (BoxIterator cell(b);cell.ok();++cell)
            {
                for (int c=0;c<NUM_VARS;++c) {double v=m_state_new[it()](cell(),c);out.write(reinterpret_cast<const char*>(&v),8);}
                for (int c=0;c<T7RHS::components;++c) {double v=rhs(cell(),c);out.write(reinterpret_cast<const char*>(&v),8);}
                for (int c:{c_GaussE,c_GaussB}) {double v=gauss(cell(),c);out.write(reinterpret_cast<const char*>(&v),8);}
            }
        }
        if (!out) throw std::runtime_error("T7 Maxwell output write failed");
    }
};

// Seed interpolation only. The unchanged finder subsequently solves at this N.
RHSurf resample(const RHSurf &old, int n)
{
    RHSurf s(n, old.m_NG, old.m_centre, old.m_index);
    for (int i = s.m_NG; i < n + s.m_NG; ++i)
    {
        const double x = s.m_theta[i] / old.m_d_theta - .5 + old.m_NG;
        const int j = std::floor(x);
        s.m_f[i] = old.m_f[j] + (x-j) * (old.m_f[j+1]-old.m_f[j]);
    }
    s.fill_all_ghosts();
    s.m_state = RHSurf::SolverState::FAR;
    s.m_time_step_freq = 400;
    s.m_chase_speed = .125;
    return s;
}

void record(std::ofstream &out, const RHSurf &s, double time, int stage,
            double threshold, const char *status, int updates, double seconds,
            const RHUnion &rh)
{
    if (procID() != 0) return;
    const double area = s.Area(), charge = s.Q_charge(), error = s.expansion_error();
    double mean = 0., variance = 0.;
    for (int i=s.m_NG; i<s.m_NG+s.m_n; ++i) mean += s.dA(i)*s.m_phi[i];
    mean /= area;
    for (int i=s.m_NG; i<s.m_NG+s.m_n; ++i)
        variance += s.dA(i)*std::pow(s.m_phi[i]-mean, 2);
    if (!std::isfinite(area + charge + error + mean + variance) || area <= 0. || variance < 0.)
        throw std::runtime_error("invalid offline surface measurement");
    out << std::setprecision(17) << time << ',' << s.m_index << ',' << s.m_n << ','
        << stage << ',' << threshold << ',' << status << ',' << updates << ',' << seconds
        << ',' << s.m_centre[0] << ',' << area << ',' << charge << ',' << mean << ','
        << std::sqrt(variance/area) << ',' << error << ',' << s.average_Theta_minus()
        << ',' << s.M_total();
    double cells = std::numeric_limits<double>::infinity();
    for (int j = 0; j < s.m_n; ++j)
        cells = std::min(cells, s.m_f[j+s.m_NG] / s.m_cell_dx.at(j));
    out << ',' << rh.m_interpolation_calls << ',' << rh.m_newton_iterations
        << ',' << rh.m_newton_failures << ',' << rh.m_flow_steps << ',' << cells << std::endl;
    std::ofstream shape("shape-"+std::to_string(s.m_index)+"-"+std::to_string(stage)+".dat");
    shape << std::setprecision(17) << time << ' ' << s.m_centre[0];
    for (int i=s.m_NG; i<s.m_NG+s.m_n; ++i) shape << ' ' << s.m_f[i];
    shape << '\n';
}

int run(int argc, char **argv)
{
    GRParmParse pp(argc-2, argv+2, nullptr, argv[1]);
    SimulationParameters p(pp);
    int n, max_updates, floor_window;
    double limit;
    bool t6_checkpoint_diagnostics;
    bool t7_diagnostics;
    std::string solver;
    bool check_jacobian;
    double fd_relative, step_cells;
    int backtracks;
    pp.load("offline_solver", solver, std::string("flow"));
    pp.load("offline_check_jacobian", check_jacobian, false);
    pp.load("offline_newton_fd_relative", fd_relative, std::sqrt(std::numeric_limits<double>::epsilon()));
    pp.load("offline_newton_max_step_cells", step_cells, 1.);
    pp.load("offline_newton_backtracks", backtracks, 20);
    pp.load("t7_diagnostics",t7_diagnostics,false);
    pp.load("t6_checkpoint_diagnostics", t6_checkpoint_diagnostics, false);
    pp.load("offline_points", n, 96);
    pp.load("offline_max_updates", max_updates, 2000);
    pp.load("offline_floor_window", floor_window, 64);
    pp.load("offline_seconds", limit, 235.);
    if (!p.restart_from_checkpoint || p.m_RH_num_horizons < 1 || n < 8 || n > 4096 ||
        max_updates < 1 || floor_window < 4 || !std::isfinite(limit) || limit <= 0. || limit > 1795. ||
        (solver != "flow" && solver != "newton") || !std::isfinite(fd_relative) || fd_relative <= 0. ||
        !std::isfinite(step_cells) || step_cells <= 0. || backtracks < 1 || backtracks > 100)
        throw std::runtime_error("invalid offline controls");
    GRAMR amr;
    DefaultLevelFactory<FrozenLevel> factory(amr, p);
    setupAMRObject(amr, factory);
    const double time = amr.getAMRLevels()[0]->time();
    const auto *coarse = amr.get_gramrlevels()[0];
    if (coarse->get_dx() != p.dx[0] ||
        coarse->problemDomain().domainBox().bigEnd() != p.ivN)
        throw std::runtime_error("parameter grid does not match checkpoint");
    for (const auto *level : amr.get_gramrlevels())
        if (level->time() != time &&
            (!t6_checkpoint_diagnostics || !std::isfinite(level->time()) ||
             std::abs(level->time()-time) > 64*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(time))))
            throw std::runtime_error("checkpoint levels are not synchronized");
    if (t7_diagnostics)
    {
        for (auto *level:amr.get_gramrlevels()) dynamic_cast<FrozenLevel*>(level)->t7_prepare();
        dynamic_cast<FrozenLevel*>(amr.get_gramrlevels().back())->t7_dump(p);
        return 0;
    }
    AMRInterpolator<Lagrange<4>> interp(amr, p.origin, p.dx, p.boundary_params, 0);
    amr.set_interpolator(&interp);
    RHUnion rh;
    rh.m_measure_cell_sizes = true;
    rh.m_use_newton = solver == "newton";
    rh.m_newton_fd_relative = fd_relative;
    rh.m_newton_max_step_cells = step_cells;
    rh.m_newton_backtracks = backtracks;
    // setup's existing reader is gated on t>0. nextafter also restores t=0 rows.
    rh.setup(p.m_RH_num_horizons, p.m_RH_initial_radii, p.m_RH_initial_centre,
             p.m_RH_num_points, p.m_RH_level, p.m_RH_time_step_freq,
             p.m_RH_newton_crit, p.m_RH_chase_speeds, p.m_RH_start_times,
             time == 0. ? std::nextafter(0., 1.) : time);
    std::ofstream out, trace;
    if (procID() == 0)
    {
        out.open("surfaces.csv");
        trace.open("progress.csv");
        trace << "stage,update,search_index,expansion_squared,A,Q,centre\n";
        out << "time,search_index,N_theta,stage,threshold,status,updates,seconds,centre,A,Q,phi_mean,phi_rms,expansion_squared,theta_minus,M_RN_legacy,interpolation_calls,newton_iterations,newton_failures,flow_steps,minimum_local_cells\n";
        std::ofstream skipped("skipped.csv");
        skipped << "search_index,status\n";
        for (const auto &s : rh.m_surfaces)
            if (s.m_dead || time < s.m_start_time)
                skipped << s.m_index << ',' << (s.m_dead ? "DEAD" : "DORMANT") << '\n';
    }
    for (int k=int(rh.m_surfaces.size())-1; k>=0; --k)
        if (rh.m_surfaces[k].m_dead || time < rh.m_surfaces[k].m_start_time)
        {
            rh.m_surfaces.erase(rh.m_surfaces.begin()+k);
            rh.m_outfiles.erase(rh.m_outfiles.begin()+k);
            rh.m_ffiles.erase(rh.m_ffiles.begin()+k);
        }
    if (rh.m_surfaces.empty()) return 0;
    const auto c = p.coupling_function_params;
    rh.set_coupling_params(c.alpha, c.f0, c.f1, c.f2);
    rh.set_interpolator(&interp);
    interp.refresh();
    rh.m_skip_interpolator_refresh = true;
    rh.interpolate_fields();
    for (auto &s : rh.m_surfaces)
    {
        record(out, s, time, -1, 0., "SEED_REPLAY", 0, 0., rh);
        s = resample(s, n);
        s.m_use_newton = rh.m_use_newton;
        // Frozen E snapshots have a strongly collapsed conformal factor.
        // Only this opt-in diagnostic changes the author's chase multiplier.
        if (t6_checkpoint_diagnostics) s.m_chase_speed = 2.;
    }
    if (t6_checkpoint_diagnostics)
    {
        // Fixed coordinate spheres use exactly RHSurf's EMS displacement flux.
        // They do not chase, re-centre, or supply an evolution boundary value.
        RHUnion spheres;
        spheres.set_interpolator(&interp);
        const double radii[] = {.02, .05, .1};
        for (int i=0; i<3; ++i)
        {
            spheres.m_surfaces.emplace_back(n, RHUnion::NG,
                std::vector<double>{p.center[0], 0.}, i);
            auto &s = spheres.m_surfaces.back();
            std::fill(s.m_f.begin(), s.m_f.end(), radii[i]);
        }
        spheres.set_coupling_params(c.alpha, c.f0, c.f1, c.f2);
        spheres.interpolate_fields();
        if (procID() == 0)
        {
            std::ofstream fixed("fixed-spheres.csv");
            fixed << "time,N_theta,radius,centre,A,Q\n";
            for (const auto &s : spheres.m_surfaces)
                fixed << std::setprecision(17) << time << ',' << n << ','
                      << radii[s.m_index] << ',' << s.m_centre[0] << ','
                      << s.Area() << ',' << s.Q_charge() << '\n';
        }
    }
    if (procID() == 0)
        for (size_t k=0; k<rh.m_surfaces.size(); ++k)
        {
            const auto id = std::to_string(rh.m_surfaces[k].m_index);
            rh.m_outfiles[k].close();
            rh.m_outfiles[k].open("offline_rh_surf_"+id+".dat");
            rh.m_ffiles[k].close();
            rh.m_ffiles[k].open("offline_rh_f"+id+".dat");
        }
    rh.set_coupling_params(c.alpha, c.f0, c.f1, c.f2);
    const auto flow_seed = rh.m_surfaces;
    if (check_jacobian)
    {
        rh.interpolate_fields();
        std::ofstream check("jacobian-check.csv");
        check << "surface,half_width,colors,max_absolute,max_relative,outside_band\n";
        for (auto &s : rh.m_surfaces)
        {
            const auto colored = rh.newton_jacobian(s);
            const auto columns = rh.newton_jacobian(s, false);
            double difference = 0., scale = 0., outside = 0.;
            for (int i = 0; i < s.m_n; ++i)
                for (int j = 0; j < s.m_n; ++j)
                {
                    difference = std::max(difference, std::abs(colored[i][j]-columns[i][j]));
                    scale = std::max(scale, std::abs(columns[i][j]));
                    if (std::abs(i-j) > RHSurf::expansion_half_width)
                        outside = std::max(outside, std::abs(columns[i][j]));
                }
            check << std::setprecision(17) << s.m_index << ',' << RHSurf::expansion_half_width
                  << ',' << 2*RHSurf::expansion_half_width+1 << ',' << difference
                  << ',' << difference/scale << ',' << outside << '\n';
            if (!std::isfinite(difference+scale+outside) || scale == 0. ||
                difference > 1e-10*scale || outside > 1e-10*scale)
                throw std::runtime_error("colored Newton Jacobian check failed");
        }
    }
    const auto start = std::chrono::steady_clock::now();
    auto seconds = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); };
    const double schedule[] = {1e-7, 1e-10, 1e-12};
    for (int stage=0; stage<3; ++stage)
    {
        rh.m_thresh_super_low = schedule[stage];
        std::vector<double> errors;
        double previous_min = 0., previous_max = 0.;
        int updates = 0;
        const char *stop = "UPDATE_CAP";
        bool fallback = false;
        for (; updates<max_updates;)
        {
            rh.update(time, 0);
            ++updates;
            if (rh.m_use_newton && rh.m_newton_failures)
            {
                fallback = true;
                break;
            }
            double worst = 0.;
            for (const auto &s : rh.m_surfaces)
            {
                const double e = s.expansion_error();
                if (s.m_dead || !std::isfinite(e)) throw std::runtime_error("offline surface died");
                worst = std::max(worst, e);
                if (procID() == 0)
                    trace << std::setprecision(17) << stage << ',' << updates << ',' << s.m_index
                          << ',' << e << ',' << s.Area() << ',' << s.Q_charge() << ','
                          << s.m_centre[0] << std::endl;
            }
            if (worst <= schedule[stage]) { stop = "FOUND"; break; }
            if (seconds() >= limit) { stop = "TIME_CAP"; break; }
            errors.push_back(worst);
            if (errors.size() == size_t(floor_window))
            {
                const auto range = std::minmax_element(errors.begin(), errors.end());
                const double low = *range.first, high = *range.second;
                // A transient can rise before relaxing. Compare consecutive
                // windows, not time since the best value at the original seed.
                // ponytail: a worst-residual plateau stops the group; use
                // per-surface windows if mixed convergence needs independent stops.
                if (std::abs(low-previous_min) < .01*low &&
                    std::abs(high-previous_max) < .01*high)
                {
                    pout() << "OFFLINE_FLOOR min=" << std::setprecision(17) << low
                           << " max=" << high << '\n';
                    stop = "FLOOR"; break;
                }
                previous_min = low; previous_max = high;
                errors.clear();
            }
        }
        if (rh.m_use_newton && (fallback || std::string(stop) != "FOUND"))
        {
            // A failed Newton search is a transaction: replay the exact original
            // resampled seed with the unchanged flow controls and row time cap.
            // Probe time remains charged to that cap.
            pout() << "OFFLINE_NEWTON_ROLLBACK seconds=" << std::setprecision(17)
                   << seconds() << " reason=" << (fallback ? "STEP_FAILED" : stop) << '\n';
            rh.m_surfaces = flow_seed;
            rh.m_use_newton = false;
            for (auto &s : rh.m_surfaces) s.m_use_newton = false;
            stage = -1;
            continue;
        }
        for (const auto &s : rh.m_surfaces)
            record(out, s, time, stage, schedule[stage],
                   s.expansion_error() <= schedule[stage] ? "FOUND" : stop, updates, seconds(), rh);
        if (std::string(stop) != "FOUND") return 1;
    }
    return 0;
}

int main(int argc, char **argv)
{
    mainSetup(argc, argv);
    int status;
    try { status = run(argc, argv); }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; status = 2; }
    mainFinalize();
    return status;
}
