/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRAMRLevel.hpp"

GRAMRLevel::GRAMRLevel(GRAMR &gr_amr, const SimulationParameters &a_p,
                       int a_verbosity)
    : m_gr_amr(gr_amr), m_p(a_p), m_verbosity(a_verbosity),
      m_t7(a_p.data_path),m_t13(a_p.data_path),m_num_ghosts(a_p.num_ghosts)
{
    if (m_verbosity)
        pout() << "GRAMRLevel constructor" << endl;
}

GRAMRLevel::~GRAMRLevel() {}

void GRAMRLevel::t13_record(int phase,const GRLevelData &data,double time,int growth)
{
    m_t13.record(phase,m_level,data,time,m_dx,m_dt,m_p.center[0],growth);
}

void GRAMRLevel::t13_snapshot()
{
    if (!m_t13.selected(m_level,m_p.max_level)) return;
    // Fill a COPY: the recorder never changes the evolved state or its ghosts.
    GRLevelData sample;
    sample.define(m_state_new.disjointBoxLayout(),NUM_VARS,m_state_new.ghostVect());
    for (DataIterator it=sample.dataIterator();it.ok();++it)
        sample[it()].copy(m_state_new[it()]);
    sample.exchange(m_exchange_copier);
    fillBdyGhosts(sample);
    t13_record(50,sample,m_time);
}

std::pair<double,double> GRAMRLevel::t7_faces() const
{
    auto extent=[](const GRAMRLevel *p) {
        if (!p) return -1.;
        double a=0.;const auto &g=p->m_state_new.disjointBoxLayout();
        for (LayoutIterator it=g.layoutIterator();it.ok();++it) a=std::max(a,(g[it()].bigEnd(1)+1)*p->m_dx);
        return a;
    };
    return {extent(this),extent(m_finer_level_ptr?gr_cast(m_finer_level_ptr):nullptr)};
}
void GRAMRLevel::t7_record(int phase,const GRLevelData &data,int growth,int first,int n,double h)
{
    if (!m_t7.enabled || m_level<4 || m_level>6) return;
    auto f=t7_faces();if (h==0.) h=m_dx;
    if (phase==40 || phase==10)
    {
        const bool parent=phase==40;
        auto regions=T7OperationRecorder::windows(parent?m_dx:2*m_dx,m_p.center[0],f.first,-1.);
        int source=0;
        for (DataIterator it=data.dataIterator();it.ok();++it,++source)
        {
            int region=0;Box valid=data.disjointBoxLayout()[it()];
            for (const auto &window:regions)
            {
                IntVectSet cells;
                if (parent)
                {
                    auto dep=T7OperationRecorder::dependencies(window&m_state_new.disjointBoxLayout()[it()]);
                    IntVectSet mirrored=dep;
                    for (IVSIterator iv(dep);iv.ok();++iv) if (iv()[1]<0)
                    { auto p=iv();p[1]=-p[1]-1;mirrored|=p; }
                    for (const auto &s:m_point_transfer.t7_ghost_stencils()[it()]) if (mirrored.contains(s.fine))
                        for (int y=0;y<6;++y) for (int x=0;x<6;++x)
                            cells|=s.first+IntVect(D_DECL(x,y,0));
                }
                else
                {
                    Box coarse=coarsen(valid,2)&window;
                    for (BoxIterator bit(coarse);bit.ok();++bit)
                        for (int y=-2;y<=3;++y) for (int x=-2;x<=3;++x)
                            cells|=2*bit()+IntVect(D_DECL(x,y,0));
                }
                std::vector<IntVect> selected;
                for (IVSIterator iv(cells);iv.ok();++iv)
                { if (!data[it()].box().contains(iv())) MayDay::Error("T7 incomplete operator support");selected.push_back(iv()); }
                m_t7.frame(phase,m_level,source*16+region,data[it()],valid,selected,first,n,m_time,h,m_dt,f.first,-1.);
                ++region;
            }
        }
        return;
    }
    if (phase==40) f.second=-1.; // only this level's own ghosts use the parent
    auto regions=T7OperationRecorder::windows(h,m_p.center[0],f.first,f.second);
    m_t7.record(phase,m_level,data,regions,growth,first,n,m_time+((phase==6 || phase==7)?m_dt:0.),h,m_dt,f.first,f.second);
}

void GRAMRLevel::define(AMRLevel *a_coarser_level_ptr,
                        const Box &a_problem_domain, int a_level,
                        int a_ref_ratio)
{
    ProblemDomain physdomain(a_problem_domain);
    define(a_coarser_level_ptr, physdomain, a_level, a_ref_ratio);

    // Also define the boundaries
    m_boundaries.define(m_dx, m_p.center, m_p.boundary_params, physdomain,
                        m_num_ghosts);
}

void GRAMRLevel::define(AMRLevel *a_coarser_level_ptr,
                        const ProblemDomain &a_problem_domain, int a_level,
                        int a_ref_ratio)
{
    if (m_verbosity)
        pout() << "GRAMRLevel::define " << a_level << endl;

    AMRLevel::define(a_coarser_level_ptr, a_problem_domain, a_level,
                     a_ref_ratio);

    if (a_coarser_level_ptr)
    {
        GRAMRLevel *coarser_level_ptr = gr_cast(a_coarser_level_ptr);
        m_dx = coarser_level_ptr->m_dx / Real(coarser_level_ptr->m_ref_ratio);
    }
    else
    {
        m_dx = m_p.L / (a_problem_domain.domainBox().longside());
    }

    // Also define the boundaries
    m_boundaries.define(m_dx, m_p.center, m_p.boundary_params, a_problem_domain,
                        m_num_ghosts);
}

/// Do casting from AMRLevel to GRAMRLevel and stop if this isn't possible
const GRAMRLevel *GRAMRLevel::gr_cast(const AMRLevel *const amr_level_ptr)
{
    const GRAMRLevel *gr_amr_level_ptr =
        dynamic_cast<const GRAMRLevel *>(amr_level_ptr);
    if (gr_amr_level_ptr == nullptr)
    {
        MayDay::Error("in GRAMRLevel::gr_cast: amr_level_ptr is not castable "
                      "to GRAMRLevel*");
    }
    return gr_amr_level_ptr;
}

/// Do casting from AMRLevel to GRAMRLevel and stop if this isn't possible
GRAMRLevel *GRAMRLevel::gr_cast(AMRLevel *const amr_level_ptr)
{
    return const_cast<GRAMRLevel *>(
        gr_cast(static_cast<const AMRLevel *const>(amr_level_ptr)));
}

const GRLevelData &GRAMRLevel::getLevelData(VariableType var_type) const
{
    if (var_type == VariableType::evolution)
        return m_state_new;
    else
        return m_state_diagnostics;
}

bool GRAMRLevel::contains(const std::array<double, CH_SPACEDIM> &point) const
{
    const Box &domainBox = problemDomain().domainBox();
    for (int i = 0; i < CH_SPACEDIM; ++i)
    {
        if (point[i] < domainBox.smallEnd(i) - m_num_ghosts ||
            point[i] > domainBox.bigEnd(i) + m_num_ghosts)
        {
            return false;
        }
    }
    return true;
}

// advance by one timestep
Real GRAMRLevel::advance()
{
    CH_TIME("GRAMRLevel::advance");

    // if 'print_progress_only_to_rank_0', still print if it's level 0 or
    // t=restart_time
    if (!m_p.print_progress_only_to_rank_0 || (procID() == 0) ||
        m_time == m_restart_time)
        printProgress("GRAMRLevel::advance");

    // copy soln to old state to save it
    m_state_new.copyTo(m_state_new.interval(), m_state_old,
                       m_state_old.interval());

    // Specifically copy boundary cells
    copyBdyGhosts(m_state_new, m_state_old);

    // The level classes take flux-register parameters, use dummy ones here
    LevelFluxRegister *coarser_fr = nullptr;
    LevelFluxRegister *finer_fr = nullptr;
    // Default coarser level pointers to an empty LevelData
    // (Chombo usually checks isDefined() rather than != nullptr so we shouldn't
    // use the nullptr)
    const GRLevelData null_gr_level_data;
    const GRLevelData *coarser_data_old = &null_gr_level_data;
    const GRLevelData *coarser_data_new = &null_gr_level_data;

    Real t_coarser_old = 0.0;
    Real t_coarser_new = 0.0;

    // A coarser level exists
    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        coarser_data_old = &coarser_gr_amr_level_ptr->m_state_old;
        coarser_data_new = &coarser_gr_amr_level_ptr->m_state_new;

        t_coarser_new = coarser_gr_amr_level_ptr->m_time;
        t_coarser_old = t_coarser_new - coarser_gr_amr_level_ptr->m_dt;
    }

    m_rk_stage = 0;
    if (m_finer_level_ptr != nullptr)
    {
        GRAMRLevel *fine_gr_amr_level_ptr = gr_cast(m_finer_level_ptr);
        RK4LevelAdvance(m_state_new, m_state_old,
                        m_p.amr_transfer == "point"
                            ? fine_gr_amr_level_ptr->m_point_transfer.time
                            : fine_gr_amr_level_ptr->m_patcher.getTimeInterpolator(),
                        *coarser_data_old, t_coarser_old, *coarser_data_new,
                        t_coarser_new, *coarser_fr, *finer_fr, m_time, m_dt,
                        *this);
    }
    else
    {
        RK4LevelAdvance(m_state_new, m_state_old, *coarser_data_old,
                        t_coarser_old, *coarser_data_new, t_coarser_new,
                        *coarser_fr, *finer_fr, m_time, m_dt, *this);
    }

    specificAdvance();
    // enforce solution BCs - in case of updates in specificAdvance
    fillBdyGhosts(m_state_new);

    m_time += m_dt;
    return m_dt;
}

// things to do after a timestep
void GRAMRLevel::postTimeStep()
{
    if (m_verbosity)
        pout() << "GRAMRLevel::postTimeStep " << m_level << endl;

    if (m_finer_level_ptr != nullptr)
    {
        GRAMRLevel *finer_gr_amr_level_ptr = gr_cast(m_finer_level_ptr);
        if (m_p.amr_transfer == "point")
        {
            if (m_t7.enabled) {t7_record(9,m_state_new,0);finer_gr_amr_level_ptr->t7_record(8,finer_gr_amr_level_ptr->m_state_new,0);}
            finer_gr_amr_level_ptr->fillAllEvolutionGhosts();
            if (m_t7.enabled) finer_gr_amr_level_ptr->t7_record(10,finer_gr_amr_level_ptr->m_state_new,4);
            finer_gr_amr_level_ptr->m_point_transfer.restrict_to_coarse(
                m_state_new, finer_gr_amr_level_ptr->m_state_new);
            if (m_t7.enabled) t7_record(11,m_state_new,0);
        }
        else
            finer_gr_amr_level_ptr->m_coarse_average.averageToCoarse(
                m_state_new, finer_gr_amr_level_ptr->m_state_new);
        // Synchronise times to avoid floating point errors for finer levels
        finer_gr_amr_level_ptr->time(m_time);
    }

    specificPostTimeStep();
    //doAnalysis(); // some old codes by miren/robin use this

    // enforce solution BCs - this is required after the averaging
    // and postentially after specificPostTimeStep actions
    fillBdyGhosts(m_state_new);

    if (m_t13.enabled() && m_level==m_p.max_level)
    {
        t13_snapshot();
        if (m_coarser_level_ptr) gr_cast(m_coarser_level_ptr)->t13_snapshot();
        ++m_t13.steps;
        // Keep the frozen timestep: finish at the last native step inside the
        // requested window, rather than taking a shortened, inconsistent step.
        if (m_time+m_dt>m_t13.stop+32*std::numeric_limits<double>::epsilon()*m_t13.stop)
            throw T13LaunchStop{m_time};
    }

    if (m_verbosity)
        pout() << "GRAMRLevel::postTimeStep " << m_level << " finished" << endl;
}

// for examples that don't implement a computeTaggingCriterion with diagnostic
// variables
void GRAMRLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    computeTaggingCriterion(tagging_criterion, current_state);
}

// things to do before tagging cells
void GRAMRLevel::preTagCells()
{
    CH_TIME("GRAMRLevel::preTagCells");
    fillAllEvolutionGhosts(); // We need filled ghost cells to calculate
                              // gradients etc
}

// create tags
void GRAMRLevel::tagCells(IntVectSet &a_tags)
{
    CH_TIME("GRAMRLevel::tagCells");
    if (m_verbosity)
        pout() << "GRAMRLevel::tagCells " << m_level << endl;

    preTagCells();

    IntVectSet local_tags;

    const DisjointBoxLayout &level_domain = m_state_new.disjointBoxLayout();
    DataIterator dit0 = level_domain.dataIterator();
    int nbox = dit0.size();
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        const Box &b = level_domain[di];
        const FArrayBox &state_fab = m_state_new[di];
        FArrayBox invalid_box; // no array memory allocated
        const FArrayBox *diagnostics_fab = &invalid_box;
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            diagnostics_fab = &m_state_diagnostics[di];
        }

        // FAB to store value of criterion
        FArrayBox tagging_criterion(b, 1);
        computeTaggingCriterion(tagging_criterion, state_fab, *diagnostics_fab);

        const IntVect &smallEnd = b.smallEnd();
        const IntVect &bigEnd = b.bigEnd();

        D_TERM(const int xmin = smallEnd[0];, const int ymin = smallEnd[1];
               , const int zmin = smallEnd[2];)

        D_TERM(const int xmax = bigEnd[0];, const int ymax = bigEnd[1];
               , const int zmax = bigEnd[2];)

#pragma omp parallel for collapse(CH_SPACEDIM) schedule(static) default(shared)
        D_INVTERM(for (int ix = xmin; ix <= xmax; ++ix),
                  for (int iy = ymin; iy <= ymax; ++iy),
                  for (int iz = zmin; iz <= zmax; ++iz))
        {
            IntVect iv(D_DECL(ix, iy, iz));
            if (tagging_criterion(iv, 0) >= m_p.regrid_thresholds[m_level])
            {
// local_tags |= is not thread safe.
#pragma omp critical
                {
                    local_tags |= iv;
                }
            }
        }
    }

    local_tags.grow(m_p.tag_buffer_size);

    // Need to do this in two steps unless a IntVectSet::operator &=
    // (ProblemDomain) operator is defined
    Box local_tags_box = local_tags.minBox();
    local_tags_box &= m_problem_domain;
    local_tags &= local_tags_box;

    a_tags = local_tags;
}

// create tags at initialization
void GRAMRLevel::tagCellsInit(IntVectSet &a_tags)
{
    // the default is to use the standard tagging function
    tagCells(a_tags);
}

// regrid
void GRAMRLevel::regrid(const Vector<Box> &a_new_grids)
{
    CH_TIME("GRAMRLevel::regrid");

    if (m_verbosity)
        pout() << "GRAMRLevel::regrid " << m_level << endl;

    m_level_grids = a_new_grids;

    mortonOrdering(m_level_grids);
    const DisjointBoxLayout level_domain = m_grids = loadBalance(a_new_grids);

    // save data for later copy, including boundary cells
    m_state_new.copyTo(m_state_new.interval(), m_state_old,
                       m_state_old.interval());

    // Specifically copy boundary cells
    copyBdyGhosts(m_state_new, m_state_old);

    // reshape state with new grids
    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_state_new.define(level_domain, NUM_VARS, iv_ghosts);

    // maintain interlevel stuff
    defineExchangeCopier(level_domain);
    m_coarse_average.define(level_domain, NUM_VARS, m_ref_ratio);
    m_fine_interp.define(level_domain, NUM_VARS, m_ref_ratio, m_problem_domain);

    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        if (m_p.amr_transfer == "point")
            m_point_transfer.define(level_domain, coarser_gr_amr_level_ptr->m_grids,
                                    m_problem_domain, m_ref_ratio, m_num_ghosts,
                                    coarser_gr_amr_level_ptr->m_dx, m_p.center,
                                    m_p.boundary_params);
        else
            m_patcher.define(level_domain, coarser_gr_amr_level_ptr->m_grids,
                             NUM_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                             m_ref_ratio, m_num_ghosts);
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            m_patcher_diagnostics.define(
                level_domain, coarser_gr_amr_level_ptr->m_grids,
                NUM_DIAGNOSTIC_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                m_ref_ratio, m_num_ghosts);
        }

        // interpolate from coarser level
        if (m_p.amr_transfer == "point")
            m_point_transfer.fill(m_state_new,
                                  coarser_gr_amr_level_ptr->m_state_new,
                                  Interval(0, NUM_VARS - 1), true);
        else
            m_fine_interp.interpToFine(m_state_new,
                                       coarser_gr_amr_level_ptr->m_state_new);

        // also interpolate fine boundary cells
        if (m_p.amr_transfer != "point" &&
            m_p.boundary_params.nonperiodic_boundaries_exist)
        {
            m_boundaries.interp_boundaries(
                m_state_new, coarser_gr_amr_level_ptr->m_state_new, Side::Hi);
            m_boundaries.interp_boundaries(
                m_state_new, coarser_gr_amr_level_ptr->m_state_new, Side::Lo);
        }
    }

    // copy from old state
    m_state_old.copyTo(m_state_old.interval(), m_state_new,
                       m_state_new.interval());

    // Specifically copy boundary cells (only if same layout, otherwise
    // interpolated)
    copyBdyGhosts(m_state_old, m_state_new);

    // enforce solution BCs (overwriting any interpolation)
    fillBdyGhosts(m_state_new);

    m_state_old.define(level_domain, NUM_VARS, iv_ghosts);
    if (NUM_DIAGNOSTIC_VARS > 0)
    {
        m_state_diagnostics.define(level_domain, NUM_DIAGNOSTIC_VARS,
                                   iv_ghosts);
    }

    // if 'print_progress_only_to_rank_0', print progress only on regrids
    // (except for rank 0, which kept doing prints)
    // print here instead of 'postRegrid' to avoid prints in reverse level order
    if (m_p.print_progress_only_to_rank_0 && (procID() != 0))
        printProgress("GRAMRLevel::regrid");
}

/// things to do after regridding
void GRAMRLevel::postRegrid(int a_base_level)
{
    // set m_restart_time to same as the coarser level
    if (m_level > a_base_level && m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        m_restart_time = coarser_gr_amr_level_ptr->m_restart_time;
    }
}

// initialize grid
void GRAMRLevel::initialGrid(const Vector<Box> &a_new_grids)
{
    CH_TIME("GRAMRLevel::initialGrid");

    if (m_verbosity)
        pout() << "GRAMRLevel::initialGrid " << m_level << endl;

    m_level_grids = a_new_grids;

    const DisjointBoxLayout level_domain = m_grids = loadBalance(a_new_grids);

    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_state_new.define(level_domain, NUM_VARS, iv_ghosts);
    m_state_old.define(level_domain, NUM_VARS, iv_ghosts);
    if (NUM_DIAGNOSTIC_VARS > 0)
    {
        m_state_diagnostics.define(level_domain, NUM_DIAGNOSTIC_VARS,
                                   iv_ghosts);
    }

    defineExchangeCopier(level_domain);
    m_coarse_average.define(level_domain, NUM_VARS, m_ref_ratio);
    m_fine_interp.define(level_domain, NUM_VARS, m_ref_ratio, m_problem_domain);

    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        if (m_p.amr_transfer == "point")
            m_point_transfer.define(level_domain, coarser_gr_amr_level_ptr->m_grids,
                                    m_problem_domain, m_ref_ratio, m_num_ghosts,
                                    coarser_gr_amr_level_ptr->m_dx, m_p.center,
                                    m_p.boundary_params);
        else
            m_patcher.define(level_domain, coarser_gr_amr_level_ptr->m_grids,
                             NUM_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                             m_ref_ratio, m_num_ghosts);
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            m_patcher_diagnostics.define(
                level_domain, coarser_gr_amr_level_ptr->m_grids,
                NUM_DIAGNOSTIC_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                m_ref_ratio, m_num_ghosts);
        }
    }
}

// things to do after initialization
void GRAMRLevel::postInitialize() { m_restart_time = 0.; }

// compute dt
Real GRAMRLevel::computeDt()
{
    if (m_verbosity)
        pout() << "GRAMRLevel::computeDt " << m_level << endl;
    return m_dt;
}

// compute dt with initial data
Real GRAMRLevel::computeInitialDt()
{
    if (m_verbosity)
        pout() << "GRAMRLevel::computeInitialDt " << m_level << endl;

    m_dt = m_initial_dt_multiplier * m_dx;
    return m_dt;
}

DisjointBoxLayout GRAMRLevel::loadBalance(const Vector<Box> &a_grids)
{
    CH_TIME("GRAMRLevel::loadBalance");

    // load balance and create boxlayout
    Vector<int> procMap;

    // appears to be faster for all procs to do the loadbalance (ndk)
    LoadBalance(procMap, a_grids);

    if (m_verbosity == 1)
    {
        pout() << "GRAMRLevel::::loadBalance" << endl;
    }
    else if (m_verbosity > 1)
    {
        pout() << "GRAMRLevel::::loadBalance: procesor map: " << endl;
        for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
            pout() << igrid << ": " << procMap[igrid] << "  " << endl;
        }
        pout() << endl;
    }

    DisjointBoxLayout dbl(a_grids, procMap, m_problem_domain);
    dbl.close();

    return dbl;
}

// write checkpoint header
#ifdef CH_USE_HDF5
void GRAMRLevel::writeCheckpointHeader(HDF5Handle &a_handle) const
{
    if (m_verbosity)
        pout() << "GRAMRLevel::writeCheckpointHeader" << endl;

    HDF5HeaderData header;
    header.m_int["num_components"] = NUM_VARS;
    char comp_str[30];
    for (int comp = 0; comp < NUM_VARS; ++comp)
    {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = UserVariables::variable_names[comp];
    }
    header.writeToFile(a_handle);

    if (m_verbosity)
        pout() << header << endl;
}

void GRAMRLevel::writeCheckpointLevel(HDF5Handle &a_handle) const
{
    CH_TIME("GRAMRLevel::writeCheckpointLevel");

    if (m_verbosity)
        pout() << "GRAMRLevel::writeCheckpointLevel" << endl;

    char level_str[20];
    sprintf(level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;

    a_handle.setGroup(label);

    HDF5HeaderData header;

    header.m_int["ref_ratio"] = m_ref_ratio;
    header.m_int["tag_buffer_size"] = m_p.tag_buffer_size;
    header.m_real["dx"] = m_dx;
    header.m_real["dt"] = m_dt;
    header.m_real["time"] = m_time;
    header.m_box["prob_domain"] = m_problem_domain.domainBox();

    // Setup the periodicity info
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        char dir_str[20];
        sprintf(dir_str, "%d", dir);
        const std::string periodic_label =
            std::string("is_periodic_") + dir_str;
        header.m_int[periodic_label] = m_problem_domain.isPeriodic(dir);
    }

    header.writeToFile(a_handle);

    if (m_verbosity)
        pout() << header << endl;

    write(a_handle, m_state_new.boxLayout());

    // only need to write ghosts when non periodic BCs exist
    IntVect ghost_vector = IntVect::Zero;
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        ghost_vector = m_num_ghosts * IntVect::Unit;
    }
    write(a_handle, m_state_new, "data", ghost_vector);
}

void GRAMRLevel::readCheckpointHeader(HDF5Handle &a_handle)
{
    CH_TIME("GRAMRLevel::readCheckpointHeader");

    if (m_verbosity)
        pout() << "GRAMRLevel::readCheckpointHeader" << endl;

    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if (m_verbosity)
        pout() << "hdf5 header data:" << endl;
    if (m_verbosity)
        pout() << header << endl;

    // read number of components
    if (header.m_int.find("num_components") == header.m_int.end())
    {
        MayDay::Error("GRAMRLevel::readCheckpointHeader: checkpoint file does "
                      "not have num_components");
    }
    int num_comps = header.m_int["num_components"];
    if (num_comps != NUM_VARS)
    {
        MayDay::Error("GRAMRLevel::readCheckpointHeader: num_components in "
                      "checkpoint file does not match solver");
    }

    // read component names
    std::string state_name;
    char comp_str[60];
    for (int comp = 0; comp < NUM_VARS; ++comp)
    {
        sprintf(comp_str, "component_%d", comp);
        if (header.m_string.find(comp_str) == header.m_string.end())
        {
            MayDay::Error("GRAMRLevel::readCheckpointHeader: checkpoint file "
                          "does not have enough component names");
        }
        state_name = header.m_string[comp_str];
        if (state_name != UserVariables::variable_names[comp])
        {
            if (m_p.ignore_checkpoint_name_mismatch)
            {
                MayDay::Warning("GRAMRLevel::readCheckpointHeader: state_name "
                                "mismatch error silenced by user.");
            }
            else
                MayDay::Error("GRAMRLevel::readCheckpointHeader: state_name in "
                              "checkpoint does not match solver");
        }
    }
}

void GRAMRLevel::readCheckpointLevel(HDF5Handle &a_handle)
{
    CH_TIME("GRAMRLevel::readCheckpointLevel");
    if (m_verbosity)
        pout() << "GRAMRLevel::readCheckpointLevel" << endl;

    char level_str[20];
    sprintf(level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;

    a_handle.setGroup(label);

    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if (m_verbosity)
        pout() << "hdf5 header data:" << endl;
    if (m_verbosity)
        pout() << header << endl;

    // read refinement ratio
    if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
        MayDay::Error(
            "GRAMRLevel::readCheckpointLevel: file does not contain ref_ratio");
    }
    m_ref_ratio = header.m_int["ref_ratio"];

    if (m_verbosity)
        pout() << "read ref_ratio = " << m_ref_ratio << endl;

    // read dx
    if (header.m_real.find("dx") == header.m_real.end())
    {
        MayDay::Error(
            "GRAMRLevel::readCheckpointLevel: file does not contain dx");
    }
    m_dx = header.m_real["dx"];

    if (m_verbosity)
        pout() << "read dx = " << m_dx << endl;

    // Since we have fixed time steping it is better to take dt from the
    // parameter file
    computeInitialDt();
    if (m_verbosity)
        pout() << "dt = " << m_dt << endl;

    // read time
    if (header.m_real.find("time") == header.m_real.end())
    {
        MayDay::Error(
            "GRAMRLevel::readCheckpointLevel: file does not contain time");
    }
    m_time = header.m_real["time"];
    m_restart_time = m_time;
    if (m_verbosity)
        pout() << "read time = " << m_time << endl;

    // read problem domain
    if (header.m_box.find("prob_domain") == header.m_box.end())
    {
        MayDay::Error("GRAMRLevel::readCheckpointLevel: file does not contain "
                      "prob_domain");
    }
    Box domainBox = header.m_box["prob_domain"];

    // Get the periodicity info
    bool isPeriodic[SpaceDim] = {
        false}; // default to false unless other info is available
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        char dir_str[20];
        sprintf(dir_str, "%d", dir);
        const std::string periodic_label =
            std::string("is_periodic_") + dir_str;
        if (!(header.m_int.find(periodic_label) == header.m_int.end()))
        {
            isPeriodic[dir] = (header.m_int[periodic_label] == true);
        }
    }
    m_problem_domain = ProblemDomain(domainBox, isPeriodic);

    // read grids
    Vector<Box> grids;
    const int grid_status = read(a_handle, grids);
    if (grid_status != 0)
    {
        MayDay::Error("GRAMRLevel::readCheckpointLevel: file does not contain "
                      "a Vector<Box>");
    }

    // create level domain
    const DisjointBoxLayout level_domain = m_grids = loadBalance(grids);

    if (m_verbosity)
        pout() << "read level domain: " << endl;
    LayoutIterator lit = level_domain.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
        const Box &b = level_domain[lit()];
        if (m_verbosity)
            pout() << lit().intCode() << ": " << b << endl;
        m_level_grids.push_back(b);
    }
    if (m_verbosity)
        pout() << endl;

    // maintain interlevel stuff
    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;

    defineExchangeCopier(level_domain);
    m_coarse_average.define(level_domain, NUM_VARS, m_ref_ratio);
    m_fine_interp.define(level_domain, NUM_VARS, m_ref_ratio, m_problem_domain);

    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        if (m_p.amr_transfer == "point")
            m_point_transfer.define(level_domain, coarser_gr_amr_level_ptr->m_grids,
                                    m_problem_domain, m_ref_ratio, m_num_ghosts,
                                    coarser_gr_amr_level_ptr->m_dx, m_p.center,
                                    m_p.boundary_params);
        else
            m_patcher.define(level_domain, coarser_gr_amr_level_ptr->m_grids,
                             NUM_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                             m_ref_ratio, m_num_ghosts);
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            m_patcher_diagnostics.define(
                level_domain, coarser_gr_amr_level_ptr->m_grids,
                NUM_DIAGNOSTIC_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                m_ref_ratio, m_num_ghosts);
        }
    }

    // reshape state with new grids
    m_state_new.define(level_domain, NUM_VARS, iv_ghosts);
    bool redefine_data = false;
    Interval comps(0, NUM_VARS - 1);
    const int data_status = read<FArrayBox>(a_handle, m_state_new, "data",
                                            level_domain, comps, redefine_data);
    if (data_status != 0)
    {
        MayDay::Error("GRAMRLevel::readCheckpointLevel: file does not contain "
                      "state data");
    }
    m_state_old.define(level_domain, NUM_VARS, iv_ghosts);
    if (NUM_DIAGNOSTIC_VARS > 0)
    {
        m_state_diagnostics.define(level_domain, NUM_DIAGNOSTIC_VARS,
                                   iv_ghosts);
    }
}

void GRAMRLevel::writePlotLevel(HDF5Handle &a_handle) const
{
    if (m_verbosity)
        pout() << "GRAMRLevel::writePlotLevel" << endl;

    // number and index of states to print
    const std::vector<std::pair<int, VariableType>> &plot_states =
        m_p.plot_vars;
    int num_states = plot_states.size();

    if (num_states > 0)
    {
        // Setup the level string
        char levelStr[20];
        sprintf(levelStr, "%d", m_level);
        const std::string label = std::string("level_") + levelStr;

        a_handle.setGroup(label);

        // Setup the level header information
        HDF5HeaderData header;

        header.m_int["ref_ratio"] = m_ref_ratio;
        header.m_int["tag_buffer_size"] = m_p.tag_buffer_size;
        header.m_real["dx"] = m_dx;
        header.m_real["dt"] = m_dt;
        header.m_real["time"] = m_time;
        header.m_box["prob_domain"] = m_problem_domain.domainBox();

        // Setup the periodicity info
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "%d", dir);
            const std::string periodic_label =
                std::string("is_periodic_") + dir_str;
            header.m_int[periodic_label] = m_problem_domain.isPeriodic(dir);
        }

        // Write the header for this level
        header.writeToFile(a_handle);

        if (m_verbosity)
            pout() << header << endl;

        const DisjointBoxLayout &levelGrids = m_state_new.getBoxes();
        IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
        LevelData<FArrayBox> plot_data(levelGrids, num_states, iv_ghosts);

        // only need to write ghosts when non periodic BCs exist
        IntVect ghost_vector = IntVect::Zero;
        if (m_p.write_plot_ghosts)
        {
            ghost_vector = m_num_ghosts * IntVect::Unit;
            Box grown_domain_box = m_problem_domain.domainBox();
            grown_domain_box.grow(ghost_vector);
            Copier boundary_copier;
            boundary_copier.ghostDefine(
                m_state_new.disjointBoxLayout(), plot_data.disjointBoxLayout(),
                grown_domain_box, ghost_vector, ghost_vector);
            for (int comp = 0; comp < num_states; comp++)
            {
                Interval currentComp(comp, comp);
                if (plot_states[comp].second == VariableType::evolution)
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    m_state_new.copyTo(plotComps, plot_data, currentComp,
                                       boundary_copier);
                }
                else
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    if (NUM_DIAGNOSTIC_VARS > 0)
                    {
                        m_state_diagnostics.copyTo(plotComps, plot_data,
                                                   currentComp);
                    }
                }
            }
        }

        else
        {
            for (int comp = 0; comp < num_states; comp++)
            {
                Interval currentComp(comp, comp);
                if (plot_states[comp].second == VariableType::evolution)
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    m_state_new.copyTo(plotComps, plot_data, currentComp);
                }
                else
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    if (NUM_DIAGNOSTIC_VARS > 0)
                    {
                        m_state_diagnostics.copyTo(plotComps, plot_data,
                                                   currentComp);
                    }
                }
            }
        }

        plot_data.exchange(plot_data.interval());

        // Write the data for this level
        write(a_handle, levelGrids);
        write(a_handle, plot_data, "data", ghost_vector);
    }
}

void GRAMRLevel::writePlotHeader(HDF5Handle &a_handle) const
{
    if (m_verbosity)
        pout() << "GRAMRLevel::writePlotHeader" << endl;

    // number and index of states to print.
    const std::vector<std::pair<int, VariableType>> &plot_states =
        m_p.plot_vars;
    int num_states = plot_states.size();

    if (num_states > 0)
    {
        // Setup the number of components
        HDF5HeaderData header;
        header.m_int["num_components"] = num_states;

        // Setup the component names
        char compStr[30];
        for (int comp = 0; comp < num_states; ++comp)
        {
            sprintf(compStr, "component_%d", comp);
            if (plot_states[comp].second == VariableType::evolution)
            {
                header.m_string[compStr] =
                    UserVariables::variable_names[plot_states[comp].first];
            }
            else
            {
                header.m_string[compStr] =
                    DiagnosticVariables::variable_names[plot_states[comp]
                                                            .first];
            }
        }

        // Write the header
        header.writeToFile(a_handle);

        if (m_verbosity)
            pout() << header << endl;
    }
    else
    {
        MayDay::Warning("GRAMRLevel::writePlotLevel: A plot interval is "
                        "provided but no components are selected for plotting. "
                        "Plot files will be empty.");
    }
}
#endif /*ifdef CH_USE_HDF5*/

void GRAMRLevel::evalRHS(GRLevelData &rhs, GRLevelData &soln,
                         LevelFluxRegister &fineFR, LevelFluxRegister &crseFR,
                         const GRLevelData &oldCrseSoln, Real oldCrseTime,
                         const GRLevelData &newCrseSoln, Real newCrseTime,
                         Real time, Real fluxWeight)
{
    CH_TIME("GRAMRLevel::evalRHS");
    if (m_t13.stage_capture(m_level,m_p.max_level))
        m_t13.stage(m_rk_stage,time,oldCrseTime,newCrseTime,
                   oldCrseSoln.isDefined()?(m_time-oldCrseTime)/(newCrseTime-oldCrseTime):0.);
    if (m_t7.enabled)
    {
        m_t7.stage=m_rk_stage;m_t7.stage_time=time;
        m_t7.coarse_old=oldCrseTime;m_t7.coarse_new=newCrseTime;
        m_t7.step_fraction=oldCrseSoln.isDefined()?(m_time-oldCrseTime)/(newCrseTime-oldCrseTime):0.;
    }
    if (m_verbosity)
        pout() << "GRAMRLevel::evalRHS" << endl;

    soln.exchange(m_exchange_copier);

    if (oldCrseSoln.isDefined())
    {
        // "time" falls between the old and the new coarse times
        Real alpha = (time - oldCrseTime) / (newCrseTime - oldCrseTime);

        // Assuming RK4, we know that there can only be 5 different alpha so fix
        // them with tolerance to prevent floating point problems
        Real eps = 0.01;
        if (abs(alpha) < eps)
            alpha = 0.0;
        else if (abs(alpha - 0.25) < eps)
            alpha = 0.25;
        else if (abs(alpha - 0.5) < eps)
            alpha = 0.5;
        else if (abs(alpha - 0.75) < eps)
            alpha = 0.75;
        else if (abs(alpha - 1.) < eps)
            alpha = 1.0;
        else
        {
            pout() << "alpha: " << alpha << endl;
            MayDay::Error(
                "Time interpolation coefficient is incompatible with RK4.");
        }

        // Interpolate ghost cells from next coarser level in space and time
        if (m_p.amr_transfer == "point")
        {
            // Chombo's intermediate expects the fine-step START, not the
            // stage time. Stages 1 and 2 are distinct at identical time.
            const double theta = std::round(
                2. * (m_time - oldCrseTime) / (newCrseTime - oldCrseTime)) / 2.;
            if (m_t7.enabled) m_t7.step_fraction=theta;
            m_point_transfer.fill_stage(soln, theta, m_rk_stage);
            if (m_t7.enabled) t7_record(40,m_point_transfer.t7_coarse_support(),4,0,NUM_VARS,2*m_dx);
            soln.exchange(m_exchange_copier);
        }
        else
            m_patcher.fillInterp(soln, alpha, 0, 0, NUM_VARS);
    }

    if (m_p.amr_transfer == "point")
        ++m_rk_stage;

    fillBdyGhosts(soln);

    specificEvalRHS(soln, rhs, time); // Call the problem specific rhs

    // evolution of the boundaries according to conditions
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.fill_rhs_boundaries(Side::Lo, soln, rhs);
        m_boundaries.fill_rhs_boundaries(Side::Hi, soln, rhs);
    }
    if (m_t13.stage_capture(m_level,m_p.max_level)) t13_record(22,rhs,m_time,0);
}

// implements soln += dt*rhs
void GRAMRLevel::updateODE(GRLevelData &soln, const GRLevelData &rhs, Real dt)
{
    CH_TIME("GRAMRLevel::updateODE");
    if (m_t13.stage_capture(m_level,m_p.max_level))
    {m_t13.update(dt);t13_record(12,soln,m_time,0);}
    if (m_t7.enabled) {m_t7.update_dt=dt;t7_record(12,soln,0);}
    // m_grown_grids will include outer boundary ghosts in the case of
    // nonperiodic BCs but will just be the problem domain otherwise.
    soln.plus(rhs, dt, m_grown_grids);
    if (m_t13.stage_capture(m_level,m_p.max_level)) t13_record(13,soln,m_time,0);
    if (m_t7.enabled) t7_record(13,soln,0);

    specificUpdateODE(soln, rhs, dt);
    fillBdyGhosts(soln);
}

// define data holder newSoln based on existingSoln,
// including ghost cell specification
void GRAMRLevel::defineSolnData(GRLevelData &newSoln,
                                const GRLevelData &existingSoln)
{
    newSoln.define(existingSoln.disjointBoxLayout(), existingSoln.nComp(),
                   existingSoln.ghostVect());
}

// define data holder for RHS based on existingSoln including ghost cell
// specification (which in most cases is no ghost cells)
void GRAMRLevel::defineRHSData(GRLevelData &newRHS,
                               const GRLevelData &existingSoln)
{
    // only need ghosts for non periodic boundary case
    IntVect ghost_vector = IntVect::Zero;
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        ghost_vector = m_num_ghosts * IntVect::Unit;
    }
    newRHS.define(existingSoln.disjointBoxLayout(), existingSoln.nComp(),
                  ghost_vector);
}

/// copy data in src into dest
void GRAMRLevel::copySolnData(GRLevelData &dest, const GRLevelData &src)
{
    src.copyTo(src.interval(), dest, dest.interval());

    // Specifically copy boundary cells if non periodic as
    // cells outside the domain are not copied by default
    copyBdyGhosts(src, dest);
}

bool GRAMRLevel::at_level_timestep_multiple(int a_level) const
{
    double target_dt = m_p.coarsest_dt;
    for (int ilevel = 0; ilevel < a_level; ++ilevel)
    {
        target_dt /= m_p.ref_ratios[ilevel];
    }
    // get difference to nearest multiple of target_dt
    const double time_remainder = remainder(m_time, target_dt);
    return (abs(time_remainder) < m_gr_amr.timeEps() * m_p.coarsest_dt);
}

void GRAMRLevel::fillAllGhosts(const VariableType var_type,
                               const Interval &a_comps)
{
    if (var_type == VariableType::evolution)
    {
        Interval comps(a_comps.begin(),
                       std::min<int>(NUM_VARS - 1, a_comps.end()));
        fillAllEvolutionGhosts(comps);
    }
    else if (var_type == VariableType::diagnostic)
    {
        Interval comps(a_comps.begin(),
                       std::min<int>(NUM_DIAGNOSTIC_VARS - 1, a_comps.end()));
        fillAllDiagnosticsGhosts(comps);
    }
}

void GRAMRLevel::fillAllEvolutionGhosts(const Interval &a_comps)
{
    CH_TIME("GRAMRLevel::fillAllEvolutionGhosts()");
    if (m_verbosity)
        pout() << "GRAMRLevel::fillAllEvolutionGhosts" << endl;

    // If there is a coarser level then interpolate undefined ghost cells
    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        if (m_p.amr_transfer == "point")
            m_point_transfer.fill(m_state_new,
                                  coarser_gr_amr_level_ptr->m_state_new, a_comps);
        else
            m_patcher.fillInterp(m_state_new, coarser_gr_amr_level_ptr->m_state_new,
                                 a_comps.begin(), a_comps.begin(), a_comps.size());
    }
    fillIntralevelGhosts(a_comps);
}

void GRAMRLevel::fillAllDiagnosticsGhosts(const Interval &a_comps)
{
    CH_TIME("GRAMRLevel::fillAllDiagnosticsGhosts");
    if (m_verbosity)
        pout() << "GRAMRLevel::fillAllDiagnosticsGhosts" << endl;

    // If there is a coarser level then interpolate undefined ghost cells
    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        if (m_p.amr_transfer == "point")
            m_point_transfer.fill(m_state_diagnostics,
                                  coarser_gr_amr_level_ptr->m_state_diagnostics,
                                  a_comps, false, VariableType::diagnostic);
        else
            m_patcher_diagnostics.fillInterp(
                m_state_diagnostics, coarser_gr_amr_level_ptr->m_state_diagnostics,
                a_comps.begin(), a_comps.begin(), a_comps.size());
    }
    m_state_diagnostics.exchange(a_comps, m_exchange_copier);

    // We should always fill the boundary ghosts to avoid nans
    // if we have non periodic directions
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.fill_diagnostic_boundaries(Side::Hi, m_state_diagnostics,
                                                a_comps);
        m_boundaries.fill_diagnostic_boundaries(Side::Lo, m_state_diagnostics,
                                                a_comps);
    }
}

void GRAMRLevel::fillIntralevelGhosts(const Interval &a_comps)
{
    m_state_new.exchange(a_comps, m_exchange_copier);
    fillBdyGhosts(m_state_new);
}

void GRAMRLevel::fillBdyGhosts(GRLevelData &a_state, const Interval &a_comps)
{
    // enforce solution BCs after filling ghosts
    if (m_p.boundary_params.boundary_solution_enforced)
    {
        m_boundaries.fill_solution_boundaries(Side::Hi, a_state, a_comps);
        m_boundaries.fill_solution_boundaries(Side::Lo, a_state, a_comps);
    }
}

void GRAMRLevel::copyBdyGhosts(const GRLevelData &a_src, GRLevelData &a_dest)
{
    // Specifically copy boundary cells if non periodic as
    // cells outside the domain are not copied by default
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.copy_boundary_cells(Side::Hi, a_src, a_dest);
        m_boundaries.copy_boundary_cells(Side::Lo, a_src, a_dest);
    }
}

void GRAMRLevel::defineExchangeCopier(const DisjointBoxLayout &a_level_grids)
{
    // if there are Sommerfeld BCs, expand boxes along those sides
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.expand_grids_to_boundaries(m_grown_grids, a_level_grids);
    }
    else
    { // nothing to do if periodic BCs
        m_grown_grids = a_level_grids;
    }

    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_exchange_copier.exchangeDefine(m_grown_grids, iv_ghosts);
}

void GRAMRLevel::printProgress(const std::string &from) const
{
    // Work out roughly how fast the evolution is going since restart
    double speed = (m_time - m_restart_time) / m_gr_amr.get_walltime();

    // Get information on number of boxes on this level (helps with better
    // load balancing)
    const DisjointBoxLayout &level_domain = m_state_new.disjointBoxLayout();
    int nbox = level_domain.dataIterator().size();
    int total_nbox = level_domain.size();

    pout() << from << " level " << m_level << " at time " << m_time << " ("
           << speed << " M/hr)"
           << ". Boxes on this rank: " << nbox << " / " << total_nbox << endl;
}
