// Frozen production checkpoint; no initial-data setter or evolution path.
#include "GRAMR.hpp"
#include "SetupFunctions.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "RHUnion.hpp"
#include "EMSKS2Profile.hpp"

// One row per coordinate or fixed-areal sphere. The local R field comes from
// gamma_phiphi; its full Cartesian gradient enters the spatial mass formula.
void exterior(RHUnion &rh, AMRInterpolator<Lagrange<4>> &interp,
              const SimulationParameters &p, double time, int n)
{
    EMSKS2Profile profile;
    profile.load(p.emsbh_params.data_path);
    const double rhorizon = profile.get("r_h"), cx = p.emsbh_params.star_centre[0];
    auto optical_slope = [&](double r) {
        const double h=1e-4*r;
        const auto a=profile.sample(r+h), b=profile.sample(r-h);
        return (a.N*a.sigma*a.sigma/((r+h)*(r+h))-
                b.N*b.sigma*b.sigma/((r-h)*(r-h)))/(2*h);
    };
    double ring = NAN, previous = rhorizon*1.001;
    for (int j=1; j<=512; ++j)
    {
        const double next=rhorizon*(1.001+29.*j/512.);
        if (optical_slope(previous)>0 && optical_slope(next)<0)
        {
            double lo=previous, hi=next;
            for (int k=0;k<40;++k)
            {
                const double mid=(lo+hi)/2;
                if (optical_slope(mid)>0) lo=mid; else hi=mid;
            }
            ring=(lo+hi)/2;
        }
        previous=next;
    }
    if (!std::isfinite(ring)) throw std::runtime_error("outer light ring not found");
    const double radii[] = {1.2,1.5,2.,3.,5.,8.,12.,ring/rhorizon,20.,30.};
    std::ofstream out;
    if (procID() == 0)
    {
        out.open("exterior.csv");
        out << "time,N_theta,kind,target_over_rh,coordinate_mean,R_mean,R_max_departure,m_mean,m_max_departure,phi_mean,phi_max_departure,lapse_mean,lapse_max_departure,K_mean,K_max_departure,electric_E2_mean,electric_E2_max_departure,Q_S,charge_flux_per_steradian_mean,charge_flux_per_steradian_max_departure,m_static,phi_static,Q_static,lapse_static,K_static\n";
    }
    auto measure = [&](const char *kind, double target, bool areal)
    {
        RHUnion sample;
        sample.set_interpolator(&interp);
        RHSurf s(n, RHUnion::NG, {cx,0.}, 0);
        s.m_alpha = p.coupling_function_params.alpha;
        s.m_f0 = p.coupling_function_params.f0;
        s.m_f1 = p.coupling_function_params.f1;
        s.m_f2 = p.coupling_function_params.f2;
        for (int i=0; i<n; ++i) s.m_f[i+s.m_NG] = target;
        s.fill_all_ghosts();
        sample.m_surfaces.push_back(s);
        for (int iteration=0; iteration<(areal ? 16 : 1); ++iteration)
        {
            sample.interpolate_fields();
            if (!areal) break;
            double worst = 0.;
            auto &v = sample.m_surfaces[0];
            for (int i=v.m_NG; i<v.m_NG+n; ++i)
            {
                const double R = v.m_f[i]*std::sqrt(v.m_hww[i]/v.m_chi[i]);
                const double change = target/R;
                if (!std::isfinite(change) || change < .5 || change > 2.)
                    throw std::runtime_error("fixed-areal root lost bracket");
                worst = std::max(worst, std::abs(R-target)/target);
                v.m_f[i] *= change;
            }
            v.fill_all_ghosts();
            if (worst < 1e-12) { sample.interpolate_fields(); break; }
            if (iteration == 15) throw std::runtime_error("fixed-areal root did not converge");
        }
        const auto &v = sample.m_surfaces[0];
        std::vector<double> qx(n), qy(n), lapse(n);
        for (int j=0; j<n; ++j)
        {
            const int i=j+v.m_NG;
            qx[j]=cx+v.m_f[i]*std::cos(v.m_theta[i]);
            qy[j]=v.m_f[i]*std::sin(v.m_theta[i]);
        }
        InterpolationQuery query(procID()==0 ? n : 0);
        query.setCoords(0,qx.data()).setCoords(1,qy.data());
        query.addComp(c_lapse,lapse.data(),Derivative::LOCAL,VariableType::evolution);
        interp.interp(query);
#ifdef CH_MPI
        MPI_Bcast(lapse.data(), n, MPI_DOUBLE, 0, Chombo_MPI::comm);
#endif
        double mean[6]{}, deviation[6]{};
        double weight=0., coordinate=0., solid_angle=0., flux_mean=0., flux_departure=0.;
        std::vector<std::array<double,6>> values(n);
        std::vector<double> flux(n);
        for (int j=0; j<n; ++j)
        {
            const int i=j+v.m_NG;
            const double rho=v.m_f[i], R=rho*std::sqrt(v.m_hww[i]/v.m_chi[i]);
            const double x=std::cos(v.m_theta[i]), y=std::sin(v.m_theta[i]);
            const double W=v.m_hww[i]/v.m_chi[i], sq=std::sqrt(W);
            const double Wx=(v.m_dx_hww[i]-W*v.m_dx_chi[i])/v.m_chi[i];
            const double Wy=(v.m_dy_hww[i]-W*v.m_dy_chi[i])/v.m_chi[i];
            const double gx=x*sq+rho*Wx/(2*sq), gy=y*sq+rho*Wy/(2*sq);
            const auto hi=v.h_inv(i);
            const double grad2=v.m_chi[i]*(hi[0][0]*gx*gx+2*hi[0][1]*gx*gy+hi[1][1]*gy*gy);
            const double E2=v.m_chi[i]*(hi[0][0]*v.m_Ex[i]*v.m_Ex[i]
                +2*hi[0][1]*v.m_Ex[i]*v.m_Ey[i]+hi[1][1]*v.m_Ey[i]*v.m_Ey[i]);
            values[j]={R,R*(1-grad2)/2,v.m_phi[i],lapse[j],v.m_K[i],E2};
            const double w=v.dA(i);
            weight+=w; coordinate+=w*rho;
            for (int k=0;k<6;++k) mean[k]+=w*values[j][k];
            const double domega=2*M_PI*std::sin(v.m_theta[i])*v.m_d_theta;
            const double phi=v.m_phi[i];
            const double F=std::exp(-2*v.m_alpha*(v.m_f0+v.m_f1*phi+v.m_f2*phi*phi));
            const auto normal=v.s_unit_vec(i);
            flux[j]=F*(v.m_Ex[i]*normal[0]+v.m_Ey[i]*normal[1])*w/
                (std::sqrt(2*M_PI)*domega);
            flux_mean+=flux[j]*domega;
            solid_angle+=domega;
        }
        for (double &x:mean) x/=weight;
        coordinate/=weight;
        flux_mean/=solid_angle;
        for (const auto &row:values)
            for (int k=0;k<6;++k) deviation[k]=std::max(deviation[k],std::abs(row[k]-mean[k]));
        for (double value:flux) flux_departure=std::max(flux_departure,std::abs(value-flux_mean));
        const auto ref=profile.sample(areal ? target : coordinate);
        const double mref=ref.r*(1-1/(ref.g*ref.g*ref.lambda))/2;
        if (!(weight > 0.) || !std::isfinite(weight+coordinate+flux_mean+flux_departure+mref+v.Q_charge()))
            throw std::runtime_error("nonfinite exterior measurement");
        for (int k=0;k<6;++k)
            if (!std::isfinite(mean[k]+deviation[k]))
                throw std::runtime_error("nonfinite exterior invariant");
        if (procID()==0)
        {
            out << std::setprecision(17) << time << ',' << n << ',' << kind << ','
                << target/rhorizon << ',' << coordinate;
            for (int k=0;k<6;++k) out << ',' << mean[k] << ',' << deviation[k];
            out << ',' << v.Q_charge() << ',' << flux_mean << ',' << flux_departure
                << ',' << mref << ',' << ref.phi << ','
                << profile.get("Q") << ',' << ref.alpha << ',' << ref.K << '\n';
        }
    };
    for (double radius:radii) measure("coordinate",radius*rhorizon,false);
    for (double radius:radii) measure("areal",radius*rhorizon,true);
}

class FrozenLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void initialData() override { throw std::runtime_error("restart_file required"); }
    void specificEvalRHS(GRLevelData &, GRLevelData &, double) override
    { throw std::runtime_error("offline snapshot must never advance"); }
    void computeTaggingCriterion(FArrayBox &, const FArrayBox &) override
    { throw std::runtime_error("offline snapshot must never regrid"); }
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

void record(std::ofstream &out, const RHSurf &s, AMRInterpolator<Lagrange<4>> &interp,
            double dx, double time, int stage, double threshold,
            const char *status, int updates, double seconds)
{
    const double area = s.Area(), charge = s.Q_charge(), error = s.expansion_error();
    // Radial perturbations use the same interpolator and the same angular shape.
    RHUnion probe;
    probe.m_surfaces.push_back(s);
    probe.set_interpolator(&interp);
    const double step = std::min(dx/8., *std::min_element(s.m_f.begin()+s.m_NG,
                                               s.m_f.begin()+s.m_NG+s.m_n)/8.);
    auto displaced = [&](double shift) {
        auto &v = probe.m_surfaces[0];
        v = s;
        for (int i=0; i<v.m_n; ++i) v.m_f[i+v.m_NG] += shift;
        v.fill_all_ghosts();
        probe.interpolate_fields();
        return v;
    };
    const RHSurf plus = displaced(step), minus = displaced(-step);
    double theta_max = 0., position_max = 0.;
    for (int i=s.m_NG; i<s.m_NG+s.m_n; ++i)
    {
        const double deriv = (plus.Theta_plus(i)-minus.Theta_plus(i))/(2*step);
        const double residual = std::abs(s.Theta_plus(i));
        theta_max = std::max(theta_max, residual);
        position_max = std::max(position_max, deriv == 0. ? INFINITY : residual/std::abs(deriv));
    }
    const double delta_area = std::abs(plus.Area()-minus.Area())/(2*step)*position_max;
    const double delta_charge = std::abs(plus.Q_charge()-minus.Q_charge())/(2*step)*position_max;
    double mean = 0., variance = 0.;
    for (int i=s.m_NG; i<s.m_NG+s.m_n; ++i) mean += s.dA(i)*s.m_phi[i];
    mean /= area;
    for (int i=s.m_NG; i<s.m_NG+s.m_n; ++i)
        variance += s.dA(i)*std::pow(s.m_phi[i]-mean, 2);
    if (!std::isfinite(area + charge + error + mean + variance + theta_max + position_max
                       + delta_area + delta_charge) || area <= 0. || variance < 0.)
        throw std::runtime_error("invalid offline surface measurement");
    if (procID() != 0) return;
    out << std::setprecision(17) << time << ',' << s.m_index << ',' << s.m_n << ','
        << stage << ',' << threshold << ',' << status << ',' << updates << ',' << seconds
        << ',' << s.m_centre[0] << ',' << area << ',' << charge << ',' << mean << ','
        << std::sqrt(variance/area) << ',' << error << ',' << s.average_Theta_minus()
        << ',' << s.M_total() << ',' << std::sqrt(area/(4*M_PI)) << ','
        << theta_max << ',' << position_max << ',' << delta_area << ',' << delta_charge << std::endl;
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
    int sphere_seed, skip_mass;
    pp.load("offline_points", n, 96);
    pp.load("offline_max_updates", max_updates, 2000);
    pp.load("offline_floor_window", floor_window, 64);
    pp.load("offline_seconds", limit, 235.);
    pp.load("offline_sphere_seed", sphere_seed, 0);
    pp.load("offline_skip_mass", skip_mass, 0);
    if (!p.restart_from_checkpoint || p.m_RH_num_horizons < 1 || n < 8 || n > 4096 ||
        max_updates < 1 || floor_window < 4 || !std::isfinite(limit) || limit <= 0. || limit > 1795.)
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
        if (level->time() != time) throw std::runtime_error("checkpoint levels are not synchronized");
    AMRInterpolator<Lagrange<4>> interp(amr, p.origin, p.dx, p.boundary_params, 0);
    amr.set_interpolator(&interp);
    RHUnion rh;
    // setup's existing reader is gated on t>0. nextafter also restores t=0 rows.
    rh.setup(p.m_RH_num_horizons, p.m_RH_initial_radii, p.m_RH_initial_centre,
             p.m_RH_num_points, p.m_RH_level, p.m_RH_time_step_freq,
             p.m_RH_newton_crit, p.m_RH_chase_speeds, p.m_RH_start_times,
             sphere_seed ? 0. : (time == 0. ? std::nextafter(0., 1.) : time));
    std::ofstream out, trace;
    if (procID() == 0)
    {
        out.open("surfaces.csv");
        trace.open("progress.csv");
        trace << "stage,update,search_index,expansion_squared,A,Q,centre\n";
        out << "time,search_index,N_theta,stage,threshold,status,updates,seconds,centre,A,Q,phi_mean,phi_rms,expansion_squared,theta_minus,M_RN_legacy,R_areal,theta_max,position_delta,area_delta,charge_delta\n";
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
    if (skip_mass) exterior(rh, interp, p, time, n);
    rh.interpolate_fields();
    for (auto &s : rh.m_surfaces)
    {
        record(out, s, interp, coarse->get_dx(), time, -1, 0., "SEED_REPLAY", 0, 0.);
        s = resample(s, n);
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
    rh.interpolate_fields();
    const auto start = std::chrono::steady_clock::now();
    auto seconds = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); };
    const double schedule[] = {1e-7, 1e-10, skip_mass ? 1e-14 : 1e-12};
    for (int stage=0; stage<3; ++stage)
    {
        rh.m_thresh_super_low = schedule[stage];
        std::vector<double> errors;
        double previous_min = 0., previous_max = 0.;
        int updates = 0;
        const char *stop = "UPDATE_CAP";
        auto best = rh.m_surfaces;
        double best_error = 0.;
        for (const auto &s : best) best_error = std::max(best_error, s.expansion_error());
        if (best_error <= schedule[stage]) stop = "FOUND";
        for (; updates<max_updates && std::string(stop) != "FOUND";)
        {
            rh.update(time, 0);
            ++updates;
            // update re-centres; acceptance always uses fields at the final points.
            rh.interpolate_fields();
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
            if (worst < best_error) { best_error = worst; best = rh.m_surfaces; }
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
        if (std::string(stop) != "FOUND")
        {
            rh.m_surfaces = best;
            rh.interpolate_fields();
        }
        for (const auto &s : rh.m_surfaces)
            record(out, s, interp, coarse->get_dx(), time, stage, schedule[stage],
                   s.expansion_error() <= schedule[stage] ? "FOUND" : stop, updates, seconds());
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
