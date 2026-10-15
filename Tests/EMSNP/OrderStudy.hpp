#ifndef EMS_NP_ORDER_STUDY_HPP
#define EMS_NP_ORDER_STUDY_HPP
#include "EMSNPExtraction.hpp"
#include <chrono>
#include <filesystem>
#include <type_traits>
#include <cstdio>
#include <cstdlib>

// Preserve the checkpoint's three-ghost transfer; shift wide stencils to fit.
template <int Order,bool Bounded=true> class EMSNPCheckedLagrange : public Lagrange<Order>
{
    const InterpSource<> &m_source;
    std::array<int,CH_SPACEDIM> m_deriv;
    std::array<double,CH_SPACEDIM> m_eval;
    struct BoundedSource : public InterpSource<>
    {
        const InterpSource<> &source;
        const Box &box;
        BoundedSource(const InterpSource<> &s,const Box &b):source(s),box(b) {}
        const LevelData<FArrayBox> &getLevelData(VariableType type) const override
        { return source.getLevelData(type); }
        double get_dx(int i) const override { return source.get_dx(i); }
        bool contains(const std::array<double,CH_SPACEDIM> &point) const override
        {
            for (int i=0;i<CH_SPACEDIM;++i)
                if (point[i]<box.smallEnd(i) || point[i]>box.bigEnd(i)) return false;
            return source.contains(point);
        }
    };
    struct CheckedFab
    {
        const FArrayBox &fab;
        double get(const IntVect &iv,int comp) const
        {
            if (!fab.box().contains(iv))
            {
                // Queries run in OpenMP workers: stop before the invalid read.
                std::fputs("NP order study: interpolation stencil outside FAB\n",stderr);
                std::fflush(stderr);
                std::_Exit(2);
            }
            return fab.get(iv,comp);
        }
    };
  public:
    explicit EMSNPCheckedLagrange(const InterpSource<> &source)
        :Lagrange<Order>(source),m_source(source) {}
    void setup(const std::array<int,CH_SPACEDIM> &deriv,
        const std::array<double,CH_SPACEDIM> &eval)
    {
        m_deriv=deriv; m_eval=eval;
        if constexpr (!Bounded) Lagrange<Order>::setup(deriv,eval);
    }
    double interpData(const FArrayBox &fab,int comp)
    {
        if constexpr (Bounded)
        {
            BoundedSource source(m_source,fab.box());
            Lagrange<Order> algorithm(source);
            algorithm.setup(m_deriv,m_eval);
            return algorithm.interpData(CheckedFab{fab},comp);
        }
        else return Lagrange<Order>::interpData(CheckedFab{fab},comp);
    }
};

template <int Order,bool Checked=true,bool Bounded=true>
void ems_np_order(GRAMR &amr,const SimulationParameters &p,
    const EMSRadiationParameters &radiation,bool matter,std::ofstream &cost,int repeats)
{
    const double time=amr.getAMRLevels()[0]->time();
    auto params=p.extraction_params;
    params.data_path=(Checked ? "order"+std::to_string(Order) : "production-cost")+
        (matter ? "/matter/" : "/vacuum/");
    std::filesystem::create_directories(params.data_path);
    using Algorithm=std::conditional_t<Checked,EMSNPCheckedLagrange<Order,Bounded>,Lagrange<Order>>;
    AMRInterpolator<Algorithm> interp(amr,p.origin,p.dx,p.boundary_params,0);
    using Clock=std::chrono::steady_clock;
    auto start=Clock::now();
    interp.refresh();
    cost << Order << ',' << matter << ',' << Checked << ",refresh,-1,"
         << std::chrono::duration<double>(Clock::now()-start).count() << '\n';
    // First pass warms allocation/interpolator caches; timed calls include CSV IO.
    for (int sample=-1;sample<repeats;++sample)
    {
        start=Clock::now();
        EMSRadiationExtraction old(params,p.dx[0]*p.dt_multiplier,time,true,0.,
            p.coupling_function_params,radiation.phi_inf);
        old.execute_query(&interp);
        cost << Order << ',' << matter << ',' << Checked << ",production," << sample << ','
             << std::chrono::duration<double>(Clock::now()-start).count() << '\n';
        start=Clock::now();
        EMSNPExtraction np(params,p.dx[0]*p.dt_multiplier,time,true,0.,
            p.coupling_function_params,p.m_G_Newton,matter);
        np.execute_query(&interp);
        cost << Order << ',' << matter << ',' << Checked << ",np," << sample << ','
             << std::chrono::duration<double>(Clock::now()-start).count() << '\n';
    }
    pout() << "NP_ORDER_STUDY order=" << Order << " matter=" << matter
           << " bounds_checked=" << Checked << '\n';
}
#endif
