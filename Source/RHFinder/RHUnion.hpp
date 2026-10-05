#ifndef RHUNION_HPP_
#define RHUNION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "RHSurf.hpp"
#include "UserVariables.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Owns all tracked horizons.
class RHUnion
{
  public:
    enum : int { NG = 2 }; // ghost cells on each side of each surface array
    int    m_num_stale_repeats = 5;       // stale repeats when err > thresh_high
    double m_courant           = 0.25;     // normal chase step size
    double m_thresh_high       = 0.0001;    // above: slow chase + stale repeats
    double m_thresh_close     = 0.0005;   // below: switch from fast chase to Newton
    double m_thresh_super_low  = 0.0000001; // below: surface converged, stop
    // Frozen hierarchy only: caller must refresh once before opting in.
    bool m_skip_interpolator_refresh = false;
    bool m_use_newton = false;
    bool m_measure_cell_sizes = false;
    double m_newton_fd_relative = std::sqrt(std::numeric_limits<double>::epsilon());
    double m_newton_max_step_cells = 1.;
    int m_newton_backtracks = 20;
    unsigned long long m_interpolation_calls = 0, m_newton_iterations = 0;
    unsigned long long m_newton_failures = 0, m_flow_steps = 0;
    std::vector<RHSurf> m_surfaces;
    std::vector<std::ofstream> m_outfiles;
    std::vector<std::ofstream> m_ffiles;

    RHUnion() {}

    void setup(int a_num,
               const std::vector<double> &a_radii,
               const std::vector<double> &a_centres_x,
               const std::vector<int>    &a_num_points,
               const std::vector<int>    &a_levels,
               const std::vector<int>    &a_time_step_freqs,
               const std::vector<double> &a_newton_crits,
               const std::vector<double> &a_chase_speeds,
               const std::vector<double> &a_start_times,
               double                     a_restart_time = 0.)
    {
        m_surfaces.clear();
        m_outfiles.clear();
        m_outfiles.resize(a_num);
        m_ffiles.clear();
        m_ffiles.resize(a_num);
        for (int i = 0; i < a_num; ++i)
        {
            m_surfaces.emplace_back(a_num_points[i], NG, std::vector<double>{a_centres_x[i], 0.0}, i);
            m_surfaces.back().m_level          = a_levels[i];
            m_surfaces.back().m_time_step_freq = a_time_step_freqs[i];
            m_surfaces.back().m_newton_crit    = a_newton_crits[i];
            m_surfaces.back().m_chase_speed    = a_chase_speeds[i];
            m_surfaces.back().m_start_time     = a_start_times[i];
            std::fill(m_surfaces.back().m_f.begin(),
                      m_surfaces.back().m_f.end(), a_radii[i]);
            if (a_restart_time > 0.)
                restore_surface(m_surfaces.back(), a_restart_time);
            if (procID() == 0)
            {
                const bool is_restart = (a_restart_time > 0.);
                if (is_restart)
                    m_surfaces.back().prune(a_restart_time);
                const auto mode = is_restart
                    ? std::ios::out | std::ios::app
                    : std::ios::out | std::ios::trunc;

                m_outfiles[i].open("rh_surf_" + std::to_string(i) + ".dat", mode);
                if (!is_restart)
                    m_outfiles[i] << std::setw(16) << "#t"
                                  << std::setw(16) << "surf"
                                  << std::setw(16) << "x"
                                  << std::setw(16) << "<r>"
                                  << std::setw(16) << "Area"
                                  << std::setw(16) << "M_irr"
                                  << std::setw(16) << "P_x"
                                  << std::setw(16) << "<Theta+>"
                                  << std::setw(16) << "<Theta->"
                                  << std::setw(16) << "err"
                                  << std::setw(16) << "K0"
                                  << std::setw(16) << "K2"
                                  << std::setw(16) << "K4"
                                  << std::setw(16) << "KP2"
                                  << std::setw(16) << "KP4"
                                  << std::setw(16) << "eq_perim"
                                  << std::setw(16) << "pol_perim"
                                  << std::setw(16) << "Q_charge"
                                  << std::setw(16) << "M_total"
                                  << std::setw(8)  << "mode"
                                  << std::endl;

                m_ffiles[i].open("rh_f" + std::to_string(i) + ".dat", mode);
                if (!is_restart)
                {
                    m_ffiles[i] << std::setw(16) << "#t";
                    for (int j = 0; j < a_num_points[i]; ++j)
                        m_ffiles[i] << std::setw(16) << ("f[" + std::to_string(j) + "]");
                    m_ffiles[i] << std::endl;
                }
            }
        }

        if (procID() == 0)
        {
            pout() << "\nRHUnion::Construction  Initial surfaces are:" << std::endl;
            for (const auto &surf : m_surfaces)
                surf.hello();
        }
    }

    // Set EMS coupling parameters on all surfaces so Q_charge() uses the right coupling.
    void set_coupling_params(double alpha, double f0, double f1, double f2)
    {
        for (auto &surf : m_surfaces)
        {
            surf.m_alpha = alpha;
            surf.m_f0    = f0;
            surf.m_f1    = f1;
            surf.m_f2    = f2;
        }
    }

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    // Interpolate all needed fields at interior surface points in one combined
    // query, then store results in each surface's arrays and fill ghost cells.
    void interpolate_fields()
    {
        int total_pts = 0;
        for (const auto &surf : m_surfaces) {
            total_pts += surf.m_n;
        }

        std::vector<double> qx(total_pts), qy(total_pts);
        int offset = 0;
        for (const auto &surf : m_surfaces)
        {
            for (int i = 0; i < surf.m_n; ++i)
            {
                int ii = i + surf.m_NG; // non-ghost array index
                qx[offset + i] = surf.m_centre[0] + surf.m_f[ii] * std::cos(surf.m_theta[ii]);
                qy[offset + i] = surf.m_centre[1] + surf.m_f[ii] * std::sin(surf.m_theta[ii]);
            }
            offset += surf.m_n;
        }

        // One list defines the query, packed broadcast and surface destinations.
        const struct
        {
            int component;
            Derivative derivative;
            std::vector<double> RHSurf::*destination;
        } fields[] = {
            {c_chi, Derivative::LOCAL, &RHSurf::m_chi},
            {c_K, Derivative::LOCAL, &RHSurf::m_K},
            {c_A11, Derivative::LOCAL, &RHSurf::m_A11},
            {c_A12, Derivative::LOCAL, &RHSurf::m_A12},
            {c_A22, Derivative::LOCAL, &RHSurf::m_A22},
            {c_Aww, Derivative::LOCAL, &RHSurf::m_Aww},
            {c_h11, Derivative::LOCAL, &RHSurf::m_h11},
            {c_h12, Derivative::LOCAL, &RHSurf::m_h12},
            {c_h22, Derivative::LOCAL, &RHSurf::m_h22},
            {c_hww, Derivative::LOCAL, &RHSurf::m_hww},
            {c_chi, Derivative::dx, &RHSurf::m_dx_chi},
            {c_chi, Derivative::dy, &RHSurf::m_dy_chi},
            {c_h11, Derivative::dx, &RHSurf::m_dx_h11},
            {c_h11, Derivative::dy, &RHSurf::m_dy_h11},
            {c_h12, Derivative::dx, &RHSurf::m_dx_h12},
            {c_h12, Derivative::dy, &RHSurf::m_dy_h12},
            {c_h22, Derivative::dx, &RHSurf::m_dx_h22},
            {c_h22, Derivative::dy, &RHSurf::m_dy_h22},
            {c_hww, Derivative::dx, &RHSurf::m_dx_hww},
            {c_hww, Derivative::dy, &RHSurf::m_dy_hww},
            {c_Ex, Derivative::LOCAL, &RHSurf::m_Ex},
            {c_Ey, Derivative::LOCAL, &RHSurf::m_Ey},
            {c_phi, Derivative::LOCAL, &RHSurf::m_phi},
        };
        std::vector<double> values(std::size(fields) * total_pts);
        InterpolationQuery query(procID() == 0 ? total_pts : 0);
        query.setCoords(0, qx.data()).setCoords(1, qy.data());
        for (size_t c = 0; c < std::size(fields); ++c)
            query.addComp(fields[c].component, values.data() + c * total_pts,
                          fields[c].derivative, VariableType::evolution);

        // Every rank participates: non-root ranks answer queries for their boxes.
        m_interpolator->interp(query);
        ++m_interpolation_calls;
        std::vector<double> cell_dx;
        if (m_use_newton || m_measure_cell_sizes)
        {
            cell_dx.resize(total_pts);
            if (procID() == 0)
                for (int i = 0; i < total_pts; ++i)
                {
                    const auto &dx = m_interpolator->get_query_dx(i);
                    cell_dx[i] = *std::max_element(dx.begin(), dx.end());
                }
#ifdef CH_MPI
            MPI_Bcast(cell_dx.data(), cell_dx.size(), MPI_DOUBLE, 0, Chombo_MPI::comm);
#endif
        }
#ifdef CH_MPI
        MPI_Bcast(values.data(), values.size(), MPI_DOUBLE, 0, Chombo_MPI::comm);
#endif

        offset = 0;
        for (auto &surf : m_surfaces)
        {
            for (size_t c = 0; c < std::size(fields); ++c)
                std::copy_n(values.data() + c * total_pts + offset, surf.m_n,
                            (surf.*fields[c].destination).begin() + surf.m_NG);
            surf.fill_all_ghosts();
            if (!cell_dx.empty())
                surf.m_cell_dx.assign(cell_dx.begin() + offset,
                                      cell_dx.begin() + offset + surf.m_n);
            offset += surf.m_n;
        }
    }

    // Forward differences: exactly one fresh field query per color (or column).
    // Each reflected row stencil has at most one independent point of each color.
    std::vector<std::vector<double>> newton_jacobian(RHSurf &surf, bool colored = true)
    {
        constexpr int b = RHSurf::expansion_half_width;
        const RHSurf saved = surf;
        std::vector<double> base(surf.m_n), h(surf.m_n);
        for (int j = 0; j < surf.m_n; ++j)
        {
            const double radius = saved.m_f[j + saved.m_NG];
            const double length = std::max(radius, saved.m_cell_dx.at(j));
            // Account for rounding in centre + f*cos(theta), not just f itself.
            const double coordinate = std::max(length, std::max(
                std::abs(saved.m_centre[0]), std::abs(saved.m_centre[1])));
            h[j] = m_newton_fd_relative * std::sqrt(length * coordinate);
            base[j] = saved.Theta_plus(j + saved.m_NG);
        }
        std::vector<std::vector<double>> J(surf.m_n, std::vector<double>(surf.m_n));
        const int stride = colored ? 2*b + 1 : surf.m_n;
        try
        {
            for (int color = 0; color < stride; ++color)
            {
                surf = saved;
                for (int j = color; j < surf.m_n; j += stride)
                    surf.m_f[j + surf.m_NG] += h[j];
                surf.fill_all_ghosts();
                interpolate_fields();
                for (int j = color; j < surf.m_n; j += stride)
                    for (int i = colored ? std::max(0, j-b) : 0;
                         i < (colored ? std::min(surf.m_n, j+b+1) : surf.m_n); ++i)
                        J[i][j] = (surf.Theta_plus(i + surf.m_NG) - base[i]) / h[j];
            }
        }
        catch (...)
        {
            surf = saved;
            throw;
        }
        surf = saved; // restore geometry AND its fields, including mirrored ghosts
        return J;
    }

    // Returns false with the complete pre-step state restored; caller uses flow.
    bool newton_step(RHSurf &surf)
    {
        static_assert(RHSurf::expansion_half_width == 2,
                      "Expansion stencil changed: replace invert_banded5 with banded LU");
        const RHSurf saved = surf;
        ++m_newton_iterations;
        try
        {
            const double error = saved.expansion_error();
            const auto J = newton_jacobian(surf);
            const auto inverse = invert_banded5(J);
            std::vector<double> step(surf.m_n), residual(surf.m_n);
            double alpha = 1.;
            for (int i = 0; i < surf.m_n; ++i)
                residual[i] = saved.Theta_plus(i + saved.m_NG);
            for (int i = 0; i < surf.m_n; ++i)
            {
                for (int j = 0; j < surf.m_n; ++j) step[i] -= inverse[i][j]*residual[j];
                if (!std::isfinite(step[i])) throw std::runtime_error("invalid Newton step");
                if (step[i] != 0.) alpha = std::min(alpha,
                    m_newton_max_step_cells * saved.m_cell_dx[i] / std::abs(step[i]));
            }
            // Guard the unpivoted existing band inverse with its actual solve residual.
            double defect = 0., scale = 0.;
            for (int i = 0; i < surf.m_n; ++i)
            {
                double r = residual[i];
                for (int j = 0; j < surf.m_n; ++j) r += J[i][j]*step[j];
                defect = std::max(defect, std::abs(r));
                scale = std::max(scale, std::abs(residual[i]));
            }
            if (!std::isfinite(defect) || defect > 1e-8*scale)
                throw std::runtime_error("Newton band solve residual failed");
            for (int trial = 0; trial < m_newton_backtracks; ++trial, alpha *= .5)
            {
                surf = saved;
                bool valid = true;
                for (int j = 0; j < surf.m_n; ++j)
                {
                    const int ii = j + surf.m_NG;
                    surf.m_f[ii] += alpha*step[j];
                    valid &= surf.m_f[ii] > surf.m_rmin && surf.m_f[ii] < surf.m_rmax;
                }
                if (!valid) continue;
                surf.fill_all_ghosts();
                interpolate_fields();
                // Also bound against any finer cells entered by this trial.
                for (int j = 0; j < surf.m_n; ++j)
                    valid &= std::abs(alpha*step[j]) <= m_newton_max_step_cells *
                        std::min(saved.m_cell_dx[j], surf.m_cell_dx[j]);
                const double next = surf.expansion_error();
                if (valid && std::isfinite(next) && next <= (1.-1e-4*alpha)*error)
                    return true;
            }
        }
        catch (const std::exception &e)
        {
            pout() << "RHFinder Newton fallback: " << e.what() << '\n';
        }
        surf = saved;
        ++m_newton_failures;
        return false;
    }

    // One finder step for all surfaces assigned to a_level.
    // Flow: dormancy check → interpolate → classify → chase loop → output.
    void update(double a_time, int a_level)
    {
        // skip entirely if no surface runs on this AMR level
        bool any = false;
        for (const auto &s : m_surfaces)
            if (s.m_level == a_level) { any = true; break; }
        if (!any) return;

        const auto t_start = std::chrono::steady_clock::now();

        // reset any surface whose mean radius has drifted out of bounds
        for (auto &surf : m_surfaces)
        {
            if (surf.m_dead || surf.m_state == RHSurf::SolverState::DORMANT) continue;
            const double r = surf.average_f();
            if (r <= surf.m_rmin || r >= surf.m_rmax)
            {
                pout() << "RHFinder: surface " << surf.m_index
                       << " out of bounds (mean_r=" << r << "), resetting to r=5\n";
                surf.reset(5.0);
            }
        }

        // update dormancy: surfaces before their start time stay/become dormant;
        // surfaces that have just passed their start time wake to FAR for classification
        for (auto &surf : m_surfaces)
        {
            if (surf.m_dead) continue;
            if (a_time < surf.m_start_time)
                surf.m_state = RHSurf::SolverState::DORMANT;
            else if (surf.m_state == RHSurf::SolverState::DORMANT)
                surf.m_state = RHSurf::SolverState::FAR;
        }

        // pull fresh grid data for all surfaces in one MPI round-trip
        if (!m_skip_interpolator_refresh) m_interpolator->refresh();
        interpolate_fields();

        // re-centre all surfaces using the fresh field data
        for (auto &surf : m_surfaces)
            surf.re_centre();
        // Newton must differentiate data at the current centre, including its
        // point-local position channel. Default flow retains its original path.
        if (m_use_newton) interpolate_fields();

        const int n = (int)m_surfaces.size();

        // measure pre-chase expansion error and assign solver regime;
        // this regime is held fixed for the entire chase below
        std::vector<double> errs(n);
        for (int k = 0; k < n; ++k)
        {
            auto &surf = m_surfaces[k];
            if (surf.m_state == RHSurf::SolverState::DORMANT) { errs[k] = 0.; continue; }
            errs[k] = surf.expansion_error();
            if      (errs[k] <= m_thresh_super_low) surf.m_state = RHSurf::SolverState::FOUND;
            else if (errs[k] <= m_thresh_high)      surf.m_state = RHSurf::SolverState::CLOSE;
            else                                    surf.m_state = RHSurf::SolverState::FAR;
            if (m_use_newton && surf.m_use_newton && !surf.m_dead &&
                surf.m_level == a_level && errs[k] > m_thresh_super_low)
            {
                if (newton_step(surf)) errs[k] = 0.; // one iteration per update
                else surf.m_use_newton = false; // stop probing; resume the existing flow
            }
        }

        // precompute per-surface courants; max_iters is the longest active surface's quota
        std::vector<double> courants(n);
        int max_iters = 0;
        for (int k = 0; k < n; ++k)
        {
            courants[k] = m_courant * m_surfaces[k].m_chase_speed;
            if (!m_surfaces[k].m_dead &&
                m_surfaces[k].m_state != RHSurf::SolverState::DORMANT &&
                m_surfaces[k].m_level == a_level &&
                errs[k] > m_thresh_super_low)
                max_iters = std::max(max_iters, m_surfaces[k].m_time_step_freq);
        }

        // chase — step-first so all active surfaces advance in lockstep and share
        // the same interpolated field snapshot at every iteration.
        // FAR surfaces take extra stale steps to close the gap faster.
        for (int i = 0; i < max_iters; ++i)
        {
            for (int k = 0; k < n; ++k)
            {
                auto &surf = m_surfaces[k];
                if (surf.m_level != a_level || surf.m_dead) continue;
                if (surf.m_state == RHSurf::SolverState::DORMANT) continue;
                if (errs[k] <= m_thresh_super_low) continue; // FOUND: nothing to do
                if (i >= surf.m_time_step_freq) continue;    // surface finished its quota

                try
                {
                    surf.chase_step(courants[k]);
                    ++m_flow_steps;
                    if (errs[k] > m_thresh_high) // FAR: stale repeats to close gap faster
                        for (int j = 0; j < m_num_stale_repeats; ++j)
                        {
                            surf.chase_step(courants[k]);
                            ++m_flow_steps;
                        }
                }
                catch (const std::exception &e)
                {
                    surf.m_dead = true;
                    pout() << "\nRHFinder: surface " << k << " marked dead: "
                           << e.what() << "\n";
                }
            }
            interpolate_fields(); // one batched MPI call refreshes all surfaces

            // early exit if every active surface on this level is converged
            bool all_found = true;
            for (int k = 0; k < n; ++k)
            {
                const auto &surf = m_surfaces[k];
                if (surf.m_level != a_level || surf.m_dead) continue;
                if (surf.m_state == RHSurf::SolverState::DORMANT) continue;
                if (surf.expansion_error() > m_thresh_super_low) { all_found = false; break; }
            }
            if (all_found) break;
        }

        // re-evaluate states post-chase so output and display reflect the final error
        for (int k = 0; k < n; ++k)
        {
            auto &surf = m_surfaces[k];
            if (surf.m_dead) continue;
            if (surf.m_state == RHSurf::SolverState::DORMANT) continue;
            const double err = surf.expansion_error();
            if      (err <= m_thresh_super_low) surf.m_state = RHSurf::SolverState::FOUND;
            else if (err <= m_thresh_high)      surf.m_state = RHSurf::SolverState::CLOSE;
            else                                surf.m_state = RHSurf::SolverState::FAR;
        }

        // write diagnostics (rhout) and full f(theta) profile (fout) for all surfaces
        rhout(a_time);
        fout(a_time);

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_start).count();

        if (procID() == 0)
        {
            pout() << "\nRHFinder::update  t = " << std::fixed << std::setprecision(4) << a_time
                   << "  (" << std::setprecision(3) << elapsed << " s)";
            for (const auto &surf : m_surfaces)
                surf.fancy_hello();
        }
    }

    void rhout(double a_time)
    {
        if (procID() != 0) return;
        for (int i = 0; i < (int)m_surfaces.size(); ++i)
        {
            if (m_surfaces[i].m_state == RHSurf::SolverState::DORMANT) continue;
            m_surfaces[i].print_diagnostics(m_outfiles[i], a_time);
        }
    }

    void fout(double a_time)
    {
        if (procID() != 0) return;
        for (int i = 0; i < (int)m_surfaces.size(); ++i)
        {
            if (m_surfaces[i].m_state == RHSurf::SolverState::DORMANT) continue;
            m_surfaces[i].print_f(m_ffiles[i], a_time);
        }
    }

    // Remove rows with t > a_restart_time from all output files.
    // Call this after setup() (so surface indices exist), then re-call setup()
    // with a_append = true so subsequent writes continue from the pruned tail.
    void prune_all(double a_restart_time)
    {
        if (procID() != 0) return;
        for (const auto &surf : m_surfaces)
            surf.prune(a_restart_time);
    }

    void hello() const
    {
        for (const auto &surf : m_surfaces)
        {
            surf.hello();
        }
    }

  private:
    // Read before pruning. Only root touches history; all ranks receive exactly
    // the same centre/shape, including the saved dead flag.
    static void restore_surface(RHSurf &surf, double restart_time)
    {
        std::vector<double> saved(surf.m_n + 2);
        saved[0] = surf.m_centre[0];
        std::copy_n(surf.m_f.begin() + surf.m_NG, surf.m_n, saved.begin() + 2);
        if (procID() == 0)
        {
            auto last_row = [restart_time](const std::string &filename)
            {
                std::ifstream input(filename);
                std::string line, last;
                while (std::getline(input, line))
                {
                    std::istringstream row(line);
                    double time;
                    if (row >> time && time <= restart_time) last = line;
                }
                return last;
            };
            const auto index = std::to_string(surf.m_index);
            const auto diagnostics = last_row("rh_surf_" + index + ".dat");
            const auto shape = last_row("rh_f" + index + ".dat");
            if (diagnostics.empty() && shape.empty() && restart_time <= surf.m_start_time)
            {
                // A surface that has not started has no output history yet.
            }
            else
            {
                std::istringstream d(diagnostics), f(shape);
                double td, tf;
                int id;
                if (!(d >> td >> id >> saved[0]) || id != surf.m_index ||
                    !std::isfinite(saved[0]) || !(f >> tf) || td != tf)
                    MayDay::Error("RHFinder restart requires matching rh_surf/rh_f history");
                for (int j = 0; j < surf.m_n; ++j)
                    if (!(f >> saved[j + 2]) || !std::isfinite(saved[j + 2]) ||
                        saved[j + 2] <= 0.)
                        MayDay::Error("RHFinder restart has an invalid or short shape");
                std::string token, mode;
                if (f >> token)
                    MayDay::Error("RHFinder restart shape point count changed");
                while (d >> token) mode = token;
                saved[1] = (mode == "dead");
                pout() << "RHFinder: surface " << surf.m_index
                       << " restored centre and " << surf.m_n
                       << " radii from t = " << td << "\n";
            }
        }
#ifdef CH_MPI
        MPI_Bcast(saved.data(), saved.size(), MPI_DOUBLE, 0, Chombo_MPI::comm);
#endif
        surf.m_centre[0] = saved[0];
        surf.m_dead = saved[1] != 0.;
        std::copy_n(saved.begin() + 2, surf.m_n, surf.m_f.begin() + surf.m_NG);
        surf.fill_all_ghosts();
        surf.m_state = restart_time < surf.m_start_time
                           ? RHSurf::SolverState::DORMANT : RHSurf::SolverState::FAR;
    }

    AMRInterpolator<Lagrange<4>> *m_interpolator = nullptr;
};

#endif /* RHUNION_HPP_ */
