#include "EMSBH2DLevel.hpp"
#include "BoxIterator.H"
#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ComputePack.hpp"
#include "EMSKS2ReferenceCache.hpp"
#include "ReferenceStationaryGauge.hpp"
#include "SetValue.hpp"
#ifdef CH_MPI
#include "SPMD.H"
#endif
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

void EMSBH2DLevel::rebuild_reference()
{
    if (m_p.gauge_type != "reference_stationary") return;
    const auto start = std::chrono::steady_clock::now();
    EMSKS2Profile profile;
    profile.load(m_p.emsbh_params.data_path);
    auto next = std::make_unique<LevelData<FArrayBox>>();
    next->define(m_grids, ReferenceStationaryGauge::num_reference,
                 m_num_ghosts * IntVect::Unit);
    std::size_t cells = 0;
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
        auto &fab = (*next)[dit];
        EMSKS2ReferenceCache::fill(fab, profile, m_dx,
                                   m_p.emsbh_params.star_centre);
        cells += fab.box().numPts();
    }
    m_reference = std::move(next);
    pout() << "EMSKS2 reference cache level=" << m_level << " cells=" << cells
           << " seconds=" << std::chrono::duration<double>(
               std::chrono::steady_clock::now() - start).count() << endl;
}

void EMSBH2DLevel::impose_reference_shift(GRLevelData &state)
{
    if (m_p.gauge_type != "reference_stationary") return;
    if (!m_reference) rebuild_reference();
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
        FArrayBox &fab = state[dit];
        const FArrayBox &ref = (*m_reference)[dit];
        for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            fab(iv, c_shift1) = ref(iv, ReferenceStationaryGauge::beta_x);
            fab(iv, c_shift2) = ref(iv, ReferenceStationaryGauge::beta_y);
            fab(iv, c_B1) = 0.;
            fab(iv, c_B2) = 0.;
        }
    }
}

void EMSBH2DLevel::check_reference_floors(const GRLevelData &state) const
{
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
        const auto &fab = state[dit];
        for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            if (!(fab(iv, c_lapse) >= m_p.min_lapse) ||
                !(fab(iv, c_chi) >= m_p.min_chi))
                MayDay::Error("reference_stationary floor activation or nonfinite lapse/chi");
        }
    }
}

void EMSBH2DLevel::write_reference_diagnostics() const
{
    if (m_p.reference_diagnostics_path.empty() || m_level != m_p.max_level)
        return;
    const auto due = [this](double interval) {
        if (interval == 0.) return true;
        const double eps = 1e-12 * std::max(interval, m_dt);
        return std::floor((m_time + eps) / interval) >
               std::floor((m_time - m_dt + eps) / interval);
    };
    if (!due(m_p.reference_diagnostics_interval)) return;
    const bool write_radial = due(m_p.reference_radial_interval);
    EMSKS2Profile p;
    p.load(m_p.emsbh_params.data_path);
    const double ra = p.rho_a(), rm = p.get("r_m"), rh = p.get("r_h");
    struct D
    {
        double H2 = 0., M2 = 0., GE2 = 0., alpha_drift = 0., rate = 0.;
        int count = 0;
    };
    std::array<D, 4> masks{};
    struct Radial
    {
        std::array<double, 8> squares{}, maxima{};
        int count = 0;
    };
    constexpr int radial_bins = 60;
    std::array<Radial, radial_bins> radial{};
    double barrier_max = -std::numeric_limits<double>::infinity();
    int barrier_count = 0, nonfinite = 0;
    double min_lapse_margin = std::numeric_limits<double>::infinity();
    double min_chi_margin = min_lapse_margin;
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
        const auto &state = m_state_new[dit];
        const auto &diag = m_state_diagnostics[dit];
        const auto &old = m_state_old[dit];
        const auto &ref = (*m_reference)[dit];
        for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            const double x = (iv[0] + .5) * m_dx - m_p.emsbh_params.star_centre[0];
            const double y = (iv[1] + .5) * m_dx - m_p.emsbh_params.star_centre[1];
            const double rho = std::hypot(x, y);
            min_lapse_margin = std::min(min_lapse_margin,
                                       state(iv, c_lapse) / m_p.min_lapse);
            min_chi_margin = std::min(min_chi_margin,
                                     state(iv, c_chi) / m_p.min_chi);
            for (int v = 0; v < NUM_VARS; ++v)
                nonfinite += !std::isfinite(state(iv, v));
            if (std::abs(rho - .95 * rh) <= m_dx)
            {
                const double nx = x / rho, ny = y / rho;
                const double hr = nx * nx * state(iv, c_h11) +
                    2 * nx * ny * state(iv, c_h12) + ny * ny * state(iv, c_h22);
                const double beta = nx * state(iv, c_shift1) +
                                    ny * state(iv, c_shift2);
                const double speed = -beta + state(iv, c_lapse) *
                    std::sqrt(state(iv, c_chi) / hr);
                barrier_max = std::max(barrier_max, speed);
                ++barrier_count;
            }
            if (rho >= ra / 2 && rho <= 3 * rh && rho - 4 * m_dx >= ra / 2)
            {
                const int bin = std::min(radial_bins - 1, static_cast<int>(
                    radial_bins * (rho - ra / 2) / (3 * rh - ra / 2)));
                auto &b = radial[bin];
                const double values[8] = {
                    diag(iv, c_Ham),
                    std::hypot(diag(iv, c_Mom1), diag(iv, c_Mom2)),
                    diag(iv, c_GaussE), state(iv, c_Theta),
                    state(iv, c_Xi), state(iv, c_Lambda),
                    state(iv, c_K) - ref(iv, ReferenceStationaryGauge::K_star),
                    state(iv, c_lapse) /
                        ref(iv, ReferenceStationaryGauge::alpha_star) - 1.};
                for (int k = 0; k < 8; ++k)
                {
                    b.squares[k] += values[k] * values[k];
                    b.maxima[k] = std::max(b.maxima[k], std::abs(values[k]));
                }
                ++b.count;
            }
            if (rho - 4 * m_dx < ra / 2 || rho > 2.3 * rh) continue;
            const int m = rho < ra ? -1 : rho < rm ? 0 :
                          rho < rh ? 1 : rho < 2 * rh ? 2 : 3;
            if (m < 0) continue;
            auto &d = masks[m];
            const double H = diag(iv, c_Ham), Mx = diag(iv, c_Mom1);
            const double My = diag(iv, c_Mom2), GE = diag(iv, c_GaussE);
            d.H2 += H * H;
            d.M2 += Mx * Mx + My * My;
            d.GE2 += GE * GE;
            d.alpha_drift = std::max(d.alpha_drift,
                std::abs(state(iv, c_lapse) /
                         ref(iv, ReferenceStationaryGauge::alpha_star) - 1.));
            if (m_time > 0)
                for (int v = 0; v < NUM_VARS; ++v)
                    d.rate = std::max(d.rate,
                        std::abs((state(iv, v) - old(iv, v)) / m_dt));
            ++d.count;
        }
    }
#ifdef CH_MPI
    std::array<double, 4 * 3 + radial_bins * 8> sums{}, sums_global{};
    std::array<double, 4 * 2 + radial_bins * 8 + 1> maxima{}, maxima_global{};
    std::array<double, 2> minima{min_lapse_margin, min_chi_margin}, minima_global{};
    std::array<int, 4 + radial_bins + 2> counts{}, counts_global{};
    for (int m = 0; m < 4; ++m)
    {
        sums[3 * m] = masks[m].H2;
        sums[3 * m + 1] = masks[m].M2;
        sums[3 * m + 2] = masks[m].GE2;
        maxima[2 * m] = masks[m].alpha_drift;
        maxima[2 * m + 1] = masks[m].rate;
        counts[m] = masks[m].count;
    }
    for (int bin = 0; bin < radial_bins; ++bin)
    {
        for (int k = 0; k < 8; ++k)
        {
            sums[12 + 8 * bin + k] = radial[bin].squares[k];
            maxima[8 + 8 * bin + k] = radial[bin].maxima[k];
        }
        counts[4 + bin] = radial[bin].count;
    }
    maxima.back() = barrier_max;
    counts[4 + radial_bins] = barrier_count;
    counts[5 + radial_bins] = nonfinite;
    MPI_Allreduce(sums.data(), sums_global.data(), sums.size(), MPI_DOUBLE,
                  MPI_SUM, Chombo_MPI::comm);
    MPI_Allreduce(maxima.data(), maxima_global.data(), maxima.size(),
                  MPI_DOUBLE, MPI_MAX, Chombo_MPI::comm);
    MPI_Allreduce(minima.data(), minima_global.data(), minima.size(),
                  MPI_DOUBLE, MPI_MIN, Chombo_MPI::comm);
    MPI_Allreduce(counts.data(), counts_global.data(), counts.size(), MPI_INT,
                  MPI_SUM, Chombo_MPI::comm);
    if (procID() != 0) return;
    for (int m = 0; m < 4; ++m)
    {
        masks[m].H2 = sums_global[3 * m];
        masks[m].M2 = sums_global[3 * m + 1];
        masks[m].GE2 = sums_global[3 * m + 2];
        masks[m].alpha_drift = maxima_global[2 * m];
        masks[m].rate = maxima_global[2 * m + 1];
        masks[m].count = counts_global[m];
    }
    for (int bin = 0; bin < radial_bins; ++bin)
    {
        for (int k = 0; k < 8; ++k)
        {
            radial[bin].squares[k] = sums_global[12 + 8 * bin + k];
            radial[bin].maxima[k] = maxima_global[8 + 8 * bin + k];
        }
        radial[bin].count = counts_global[4 + bin];
    }
    barrier_max = maxima_global.back();
    barrier_count = counts_global[4 + radial_bins];
    nonfinite = counts_global[5 + radial_bins];
    min_lapse_margin = minima_global[0];
    min_chi_margin = minima_global[1];
#endif
    const bool first = m_restart_time == 0. &&
        (m_p.reference_diagnostics_interval == 0. ?
         m_time <= m_dt * (1. + 1e-12) :
         m_time - m_dt < m_p.reference_diagnostics_interval);
    std::ofstream out(m_p.reference_diagnostics_path,
        first ? std::ios::trunc : std::ios::app);
    if (!out) MayDay::Error("cannot write reference diagnostics");
    if (first)
        out << "time,time_over_M,level,mask,count,H_L2,M_L2,GaussE_L2,max_alpha_drift,max_rate,barrier_cplus_max,barrier_samples,min_lapse_margin,min_chi_margin,nonfinite,floor_count\n";
    const char *names[] = {"join", "KS_collar", "exterior", "far"};
    out << std::setprecision(17);
    for (int m = 0; m < 4; ++m)
    {
        const auto &d = masks[m];
        if (!d.count) MayDay::Error("reference diagnostic mask has no cells");
        out << m_time << ',' << m_time / p.get("M") << ',' << m_level << ','
            << names[m] << ',' << d.count << ',' << std::sqrt(d.H2 / d.count)
            << ',' << std::sqrt(d.M2 / d.count) << ','
            << std::sqrt(d.GE2 / d.count) << ',' << d.alpha_drift << ','
            << d.rate << ',' << barrier_max << ',' << barrier_count << ','
            << min_lapse_margin << ',' << min_chi_margin << ','
            << nonfinite << ",0\n";
    }
    if (!write_radial) return;
    const bool first_radial = m_restart_time == 0. &&
        (m_p.reference_radial_interval == 0. ? first :
         m_time - m_dt < m_p.reference_radial_interval);
    std::ofstream radial_out(m_p.reference_diagnostics_path + ".radial.csv",
        first_radial ? std::ios::trunc : std::ios::app);
    if (!radial_out) MayDay::Error("cannot write reference radial diagnostics");
    if (first_radial)
        radial_out << "time,time_over_M,level,bin,rho_lo,rho_hi,rho_mid,r_mid,z_mid,count,H_RMS,H_max,M_RMS,M_max,GaussE_RMS,GaussE_max,Theta_RMS,Theta_max,Xi_RMS,Xi_max,Lambda_RMS,Lambda_max,K_minus_Kstar_RMS,K_minus_Kstar_max,alpha_drift_RMS,alpha_drift_max\n";
    radial_out << std::setprecision(17);
    for (int bin = 0; bin < radial_bins; ++bin)
    {
        const double lo = ra / 2 + (3 * rh - ra / 2) * bin / radial_bins;
        const double hi = ra / 2 + (3 * rh - ra / 2) * (bin + 1) / radial_bins;
        const double mid = (lo + hi) / 2;
        const double r = p.sample(mid).r;
        radial_out << m_time << ',' << m_time / p.get("M") << ',' << m_level
                   << ',' << bin << ',' << lo << ',' << hi << ',' << mid
                   << ',' << r << ',';
        if (mid >= ra && mid <= rm)
            radial_out << (r - p.get("r_a")) / (rm - p.get("r_a"));
        radial_out << ',' << radial[bin].count;
        for (int k = 0; k < 8; ++k)
            radial_out << ',' << (radial[bin].count ?
                std::sqrt(radial[bin].squares[k] / radial[bin].count) : 0.)
                       << ',' << radial[bin].maxima[k];
        radial_out << '\n';
    }
}

void EMSBH2DLevel::postRegrid(int base_level)
{
    GRAMRLevel::postRegrid(base_level);
    if (m_p.gauge_type == "reference_stationary")
    {
        rebuild_reference();
        impose_reference_shift(m_state_new);
    }
}

#ifdef CH_USE_HDF5
void EMSBH2DLevel::postRestart()
{
    if (m_p.gauge_type == "reference_stationary")
    {
        rebuild_reference();
        impose_reference_shift(m_state_new);
    }
}
#endif



void EMSBH2DLevel::eval_reference_rhs(GRLevelData &a_soln,
                                        GRLevelData &a_rhs,
                                        const CouplingFunction &coupling_function)
{
        SetValue zero_diagnostics(0.0, Interval(c_Xi + 1, NUM_VARS - 1));
        for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
            CCZ4_params_t<ReferenceStationaryGauge::params_t> params;
            static_cast<CCZ4_base_params_t &>(params) =
                static_cast<const CCZ4_base_params_t &>(m_p.ccz4_params);
            auto &gauge = static_cast<ReferenceStationaryGauge::params_t &>(params);
            static_cast<MovingPunctureGauge::params_t &>(gauge) =
                static_cast<const MovingPunctureGauge::params_t &>(m_p.ccz4_params);
            gauge.reference = &(*m_reference)[dit];
            gauge.mass = m_p.emsbh_params.bh_mass;
            gauge.onepluslog = m_p.reference_f == "onepluslog";
            gauge.onepluslog_n = m_p.reference_onepluslog_n;
            CCZ4Cartoon<ReferenceStationaryGauge,
                        FourthOrderDerivatives, CouplingFunction>
                rhs(params, m_dx, m_p.sigma, coupling_function,
                    m_p.m_G_Newton, m_p.formulation);
            BoxLoops::loop(make_compute_pack(rhs, zero_diagnostics), a_soln[dit],
                           a_rhs[dit], m_grids[dit], disable_simd());
        }
}
