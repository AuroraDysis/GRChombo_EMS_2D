// Exercises the production point transfer and GRAMRLevel's RK4/subcycling.
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SetupFunctions.hpp"
#include <fstream>
#include <iomanip>

static void require(bool value, const char *message)
{
    if (!value)
        MayDay::Error(message);
}

static DisjointBoxLayout layout(std::initializer_list<Box> input,
                                const ProblemDomain &domain)
{
    Vector<Box> boxes;
    for (const auto &box : input)
        boxes.push_back(box);
    return DisjointBoxLayout(boxes, Vector<int>(boxes.size(), 0), domain);
}

static double polynomial(const IntVect &iv, double h, int a, int b)
{
    return std::pow((iv[0] + .5) * h, a) * std::pow((iv[1] + .5) * h, b);
}

static double smooth(const IntVect &iv, double h, bool odd)
{
    const double x = (iv[0] + .5) * h, y = (iv[1] + .5) * h;
    return std::exp(2. * x) * (odd ? std::sin(3. * y) : std::cos(3. * y));
}

static void operator_checks(const BoundaryConditions::params_t &input)
{
    auto bdy = input;
    bdy.set_is_periodic({false, false});
    bdy.set_lo_boundary({BoundaryConditions::SOMMERFELD_BC,
                         BoundaryConditions::REFLECTIVE_BC});
    bdy.set_hi_boundary({BoundaryConditions::SOMMERFELD_BC,
                         BoundaryConditions::SOMMERFELD_BC});
    bdy.vars_parity.fill(BoundaryConditions::EVEN);
    bdy.vars_parity[1] = BoundaryConditions::ODD_Y;
    std::ofstream out("operator.csv");
    out << std::unitbuf;
    out << "test,N,region,error,order\n" << std::setprecision(17);
    double previous[2][2][2] = {};
    for (int n : {32, 64, 128})
    {
        const double h = 1. / n;
        const ProblemDomain cd(Box(IntVect::Zero, IntVect(n - 1, n - 1)));
        const ProblemDomain fd = refine(cd, 2);
        auto cl = layout({cd.domainBox()}, cd);
        // Two boxes test exchange seams as well as the axis and patch edge.
        auto fl = layout({Box(IntVect(n / 2, 0), IntVect(n - 1, n - 1)),
                          Box(IntVect(n, 0), IntVect(3 * n / 2 - 1, n - 1))}, fd);
        GRLevelData coarse, fine, restricted;
        coarse.define(cl, NUM_VARS, 3 * IntVect::Unit);
        fine.define(fl, NUM_VARS, 3 * IntVect::Unit);
        restricted.define(cl, NUM_VARS);
        PointAMRTransfer point;
        point.define(fl, cl, fd, 2, 3, h, {0., 0.}, bdy);
        if (n == 32)
        {
            double previous_stage_error = 0.;
            for (double H : {.125, .0625, .03125})
            {
                coarse.setVal(1.);
                point.time.setDt(H);
                point.time.saveInitialSoln(coarse);
                GRLevelData rhs;
                rhs.define(cl, NUM_VARS, 3 * IntVect::Unit);
                const double lambda = .7;
                double k[4];
                k[0] = lambda;
                k[1] = lambda * (1. + H * k[0] / 2.);
                k[2] = lambda * (1. + H * k[1] / 2.);
                k[3] = lambda * (1. + H * k[2]);
                for (double value : k)
                {
                    rhs.setVal(value);
                    point.time.saveRHS(rhs);
                }
                double error = 0.;
                for (double theta : {0., .5})
                {
                    const double start = std::exp(lambda * theta * H);
                    const double expected[4] = {start,
                        start * (1. + lambda * H / 4.),
                        start * (1. + lambda * H / 4. + lambda * lambda * H * H / 16.),
                        start * (1. + lambda * H / 2. + lambda * lambda * H * H / 8. +
                                 lambda * lambda * lambda * H * H * H / 32.)};
                    for (int stage = 0; stage < 4; ++stage)
                    {
                        point.fill_stage(fine, theta, stage);
                        for (DataIterator dit = fine.dataIterator(); dit.ok(); ++dit)
                        {
                            const IntVect iv(fl[dit].smallEnd(0) == n / 2 ?
                                                 fl[dit].smallEnd(0) - 1 :
                                                 fl[dit].bigEnd(0) + 1, 0);
                            error = std::max(error, std::abs(fine[dit](iv, 0) - expected[stage]));
                        }
                    }
                }
                const double order = previous_stage_error > 0.
                                         ? std::log2(previous_stage_error / error) : 0.;
                require(error < 2. * std::pow(H, 4), "RK stage ghost differs from fine RK4 jet");
                if (previous_stage_error > 0.)
                    require(order > 3.8, "RK stage ghost residual is below fourth order");
                previous_stage_error = error;
                out << "stage_jet," << 1./H << ",all," << error << ',' << order << '\n';
            }
        }
        double pmax = 0., rmax = 0.;
        for (int a = 0; a <= 5; ++a)
            for (int b = 0; b <= 5; ++b)
            {
                const int comp = b % 2;
                coarse.setVal(0.);
                fine.setVal(0.);
                for (DataIterator dit = coarse.dataIterator(); dit.ok(); ++dit)
                    for (BoxIterator bit(coarse[dit].box()); bit.ok(); ++bit)
                        coarse[dit](bit(), comp) = polynomial(bit(), h, a, b);
                point.fill(fine, coarse, Interval(comp, comp), true);
                for (DataIterator dit = fine.dataIterator(); dit.ok(); ++dit)
                    for (BoxIterator bit(fine[dit].box() & fd); bit.ok(); ++bit)
                        pmax = std::max(pmax, std::abs(fine[dit](bit(), comp) -
                                                    polynomial(bit(), h / 2, a, b)));
                // Independent fine samples exercise restriction, not P then R.
                for (DataIterator dit = fine.dataIterator(); dit.ok(); ++dit)
                    for (BoxIterator bit(fine[dit].box()); bit.ok(); ++bit)
                        fine[dit](bit(), comp) = polynomial(bit(), h / 2, a, b);
                point.restrict_to_coarse(restricted, fine);
                for (DataIterator dit = restricted.dataIterator(); dit.ok(); ++dit)
                    for (BoxIterator bit(Box(IntVect(n / 4, 0),
                                             IntVect(3 * n / 4 - 1, n / 2 - 1)));
                         bit.ok(); ++bit)
                        rmax = std::max(rmax, std::abs(restricted[dit](bit(), comp) -
                                                    polynomial(bit(), h, a, b)));
            }
        require(pmax < 2e-11 && rmax < 2e-11, "degree-five polynomial transfer failed");
        out << "polynomial_prolong," << n << ",all," << pmax << ",\n";
        out << "polynomial_restrict," << n << ",all," << rmax << ",\n";
        for (int odd = 0; odd < 2; ++odd)
        {
            coarse.setVal(0.);
            fine.setVal(0.);
            for (DataIterator dit = coarse.dataIterator(); dit.ok(); ++dit)
                for (BoxIterator bit(coarse[dit].box()); bit.ok(); ++bit)
                    coarse[dit](bit(), odd) = smooth(bit(), h, odd);
            point.fill(fine, coarse, Interval(odd, odd), true);
            double errors[2][2] = {};
            for (DataIterator dit = fine.dataIterator(); dit.ok(); ++dit)
                for (BoxIterator bit(fine[dit].box() & fd); bit.ok(); ++bit)
                {
                    const double error = std::abs(fine[dit](bit(), odd) -
                                                  smooth(bit(), h / 2, odd));
                    errors[0][0] = std::max(errors[0][0], error);
                    if (bit()[1] < 2)
                        errors[0][1] = std::max(errors[0][1], error);
                }
            for (DataIterator dit = fine.dataIterator(); dit.ok(); ++dit)
                for (BoxIterator bit(fine[dit].box()); bit.ok(); ++bit)
                    fine[dit](bit(), odd) = smooth(bit(), h / 2, odd);
            point.restrict_to_coarse(restricted, fine);
            for (DataIterator dit = restricted.dataIterator(); dit.ok(); ++dit)
                for (BoxIterator bit(Box(IntVect(n / 4, 0),
                                         IntVect(3 * n / 4 - 1, n / 2 - 1)));
                     bit.ok(); ++bit)
                {
                    const double error = std::abs(restricted[dit](bit(), odd) -
                                                  smooth(bit(), h, odd));
                    errors[1][0] = std::max(errors[1][0], error);
                    if (bit()[1] < 2)
                        errors[1][1] = std::max(errors[1][1], error);
                }
            for (int op = 0; op < 2; ++op)
                for (int axis = 0; axis < 2; ++axis)
                {
                    const double order = previous[odd][op][axis] > 0.
                                             ? std::log2(previous[odd][op][axis] / errors[op][axis])
                                             : 0.;
                    // Last grid is close to double roundoff; the first pair
                    // is the registered rate check, the last is reported.
                    if (n == 64)
                        require(order > 5.5, "sixth-order smooth transfer failed");
                    previous[odd][op][axis] = errors[op][axis];
                    out << (op ? "smooth_restrict_" : "smooth_prolong_")
                        << (odd ? "odd" : "even") << ',' << n << ','
                        << (axis ? "axis" : "all") << ',' << errors[op][axis]
                        << ',' << order << '\n';
                }
        }
    }
}

// Test-only flat-space scalar wave, u_t=v, v_t=laplacian(u).
// No analytic state is used by the RHS, ghosts, restriction or regrid.
class WaveLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;
    void initialData() override
    {
        m_state_new.setVal(0.);
        for (DataIterator dit = m_state_new.dataIterator(); dit.ok(); ++dit)
            for (BoxIterator bit(m_state_new[dit].box()); bit.ok(); ++bit)
            {
                const double x = (bit()[0] + .5) * m_dx;
                double u = 0., v = 0.;
                for (int image = -1; image <= 1; ++image)
                {
                    const double s = (x - 1.5 + 8. * image) / .4;
                    const double pulse = std::exp(-s * s);
                    u += pulse;
                    v += 2. * s * pulse / .4;
                }
                m_state_new[dit](bit(), 0) = u;
                m_state_new[dit](bit(), 1) = v;
            }
    }
    void specificEvalRHS(GRLevelData &u, GRLevelData &rhs, double) override
    {
        rhs.setVal(0.);
        const FourthOrderDerivatives diff(m_dx);
        for (DataIterator dit = u.dataIterator(); dit.ok(); ++dit)
        {
            const BoxPointers pointers(u[dit], rhs[dit]);
            for (BoxIterator bit(u.disjointBoxLayout()[dit]); bit.ok(); ++bit)
            {
                const IntVect iv = bit();
                const Cell<double> cell(iv, pointers);
                double values[2] = {u[dit](iv, 1), 0.};
                for (int dir = 0; dir < 2; ++dir)
                {
                    values[1] += diff.diff2<double>(pointers.m_in_ptr[0],
                        cell.get_in_index(), pointers.m_in_stride[dir]);
                    diff.add_dissipation(values, cell, m_p.sigma, dir);
                }
                rhs[dit](iv, 0) = values[0];
                rhs[dit](iv, 1) = values[1];
            }
        }
    }
    void prePlotLevel() override { fillAllGhosts(); }
};

int main(int argc, char **argv)
{
    mainSetup(argc, argv);
    GRParmParse pp(argc - 2, argv + 2, nullptr, argv[1]);
    SimulationParameters params(pp);
    bool operators;
    pp.load("operator_checks", operators, false);
    if (operators)
        operator_checks(params.boundary_params);
    else
    {
        GRAMR amr;
        DefaultLevelFactory<WaveLevel> factory(amr, params);
        ProblemDomain domain(Box(IntVect::Zero, params.ivN));
        domain.setPeriodic(0, true);
        amr.define(1, params.ref_ratios, domain, &factory);
        amr.verbosity(0);
        amr.regridIntervals(Vector<int>(2, 0));
        amr.plotInterval(params.plot_interval);
        amr.plotPrefix(params.plot_prefix);
        amr.checkpointInterval(-1);
        amr.maxGridSize(params.max_grid_size);
        amr.blockFactor(params.block_factor);
        Vector<Vector<Box>> grids(2);
        grids[0].push_back(domain.domainBox());
        const int n = params.ivN[0] + 1;
        grids[1].push_back(Box(IntVect(n / 2, 0), IntVect(n - 1, n - 1)));
        grids[1].push_back(Box(IntVect(n, 0), IntVect(3 * n / 2 - 1, n - 1)));
        amr.setupForFixedHierarchyRun(grids);
        amr.run(params.stop_time, params.max_steps);
        amr.conclude();
    }
    mainFinalize();
    return 0;
}
