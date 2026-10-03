// One parametric RHS implementation, instantiated in separate translation
// units so the opt-in template cannot change legacy SIMD code generation.
// EMS_RHS_METHOD and EMS_RHS_GAUGE are supplied by the two callers.
void EMSBH2DLevel::EMS_RHS_METHOD(GRLevelData &a_soln,
                                         GRLevelData &a_rhs,
                                         const double a_time)
{
#ifdef EMS_RHS_DISPATCH_MOVING
    if (m_ems_gauge.differential())
    {
        specificEvalRHSMoving(a_soln, a_rhs, a_time);
        return;
    }
#endif

    if (m_t21.selected(m_level,m_time,m_dt))
        m_t21.before_projection(m_level,m_rk_stage,a_time,m_dx,
                               {m_p.center[0],m_p.center[1]},a_soln,m_p.min_chi,m_p.min_lapse);
    if (m_t7.enabled) t7_record(2,a_soln);
    const bool t13=m_t13.stage_capture(m_level,m_p.max_level);
    if (t13) t13_record(2,a_soln,m_time);
    m_t13.floors(m_level,a_soln,a_time,"stage",m_p.min_chi,m_p.min_lapse,c_chi,c_lapse);
    ////////////////////////////////////////
    // Enforce positive chi and alpha and trace free A
    if (t13)
    {
        BoxLoops::loop(TraceARemovalCartoon(),a_soln,a_soln,INCLUDE_GHOST_CELLS);
        t13_record(14,a_soln,m_time);
        BoxLoops::loop(PositiveChiAndAlpha(m_p.min_chi,m_p.min_lapse),a_soln,a_soln,INCLUDE_GHOST_CELLS);
        t13_record(15,a_soln,m_time);
    }
    else BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);
    if (t13) t13_record(3,a_soln,m_time);
    if (m_t7.enabled) t7_record(3,a_soln);


    ///////////////////////////////////
    // Coupling Function
    CouplingFunction coupling_function(m_p.coupling_function_params);


    ////////////////////////////
    // Integrated MPG
    // CCZ4Cartoon<IntegratedMovingPunctureGauge,
    //             FourthOrderDerivatives,
    //             CouplingFunction>
    // my_ccz4_cartoon(m_p.ccz4_params, m_dx, m_p.sigma, coupling_function,
    //                                        m_p.m_G_Newton, m_p.formulation);


    // //////////////////////////////
    // // MPG
    // CCZ4Cartoon<MovingPunctureGauge,
    //             FourthOrderDerivatives,
    //             CouplingFunction>
    // my_ccz4_cartoon(m_p.ccz4_params, m_dx, m_p.sigma, coupling_function,
    //                                       m_p.m_G_Newton, m_p.formulation);

    // //////////////////////////////
    // // XPG
    CCZ4Cartoon<EMS_RHS_GAUGE,
                FourthOrderDerivatives,
                CouplingFunction>
    my_ccz4_cartoon(m_p.ccz4_params, m_dx, m_p.sigma, coupling_function,
                                          m_p.m_G_Newton, m_p.formulation);


    ///////////////////////
    // zero diagnostic vars
    SetValue set_analysis_vars_zero(0.0, Interval(c_Xi + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_cartoon, set_analysis_vars_zero);
    const bool t21=m_t21.selected(m_level,m_time,m_dt);
    if (t21 && (t13 || m_t7.enabled))
        MayDay::Error("T21 production recorder cannot be combined with T7/T13 capture");
    if (t21)
    {
        m_t21.next_call();
        m_t21.set_kernel_type(typeid(my_ccz4_cartoon).name());
        for (DataIterator it=a_soln.dataIterator();it.ok();++it)
        {
            const Box valid=a_soln.disjointBoxLayout()[it()];
            FArrayBox capture(a_rhs[it()].box(),2*NUM_VARS);
            std::vector<Box> regions{valid};
            my_ccz4_cartoon.t7_capture(&capture,&regions);
            BoxLoops::loop(make_compute_pack(my_ccz4_cartoon,set_analysis_vars_zero),a_soln[it()],a_rhs[it()],valid);
            const auto &record_state=m_t21.checkpoint_source?(*m_t21.checkpoint_source)[it()]:a_soln[it()];
            m_t21.write(m_level,m_rk_stage,a_time,m_time,m_dx,m_dt,
                        {m_p.center[0],m_p.center[1]},valid,record_state,capture,a_rhs[it()]);
        }
    }
    else if (t13 || (m_t7.enabled && m_level>=4 && m_level<=6))
    {
        auto f=t7_faces();auto windows=t13?T13LaunchRecorder::windows(m_dx,m_p.center[0]):T7OperationRecorder::windows(m_dx,m_p.center[0],f.first,f.second);
        int source=0;
        for (DataIterator it=a_soln.dataIterator();it.ok();++it,++source)
        {
            FArrayBox capture(a_rhs[it()].box(),2*NUM_VARS);
            my_ccz4_cartoon.t7_capture(&capture,&windows);
            BoxLoops::loop(make_compute_pack(my_ccz4_cartoon,set_analysis_vars_zero),a_soln[it()],a_rhs[it()],a_soln.disjointBoxLayout()[it()]);
            if (t13)
            {
                m_t13.parts(m_level,source,capture,a_soln.disjointBoxLayout()[it()],m_time,m_dx,m_dt,m_p.center[0]);
                continue;
            }
            int region=0;
            for (const auto &window:windows)
            {
                Box b=window&a_soln.disjointBoxLayout()[it()];std::vector<IntVect> cells;
                for (BoxIterator bit(b);bit.ok();++bit) cells.push_back(bit());
                for (int part=0;part<2;++part)
                    m_t7.frame(20+part,m_level,source*16+region,capture,a_soln.disjointBoxLayout()[it()],cells,
                               part*NUM_VARS,NUM_VARS,m_time,m_dx,m_dt,f.first,f.second);
                ++region;
            }
        }
        t7_record(22,a_rhs,0);
    }
    else BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

#undef EMS_RHS_METHOD
#undef EMS_RHS_GAUGE
#ifdef EMS_RHS_DISPATCH_MOVING
#undef EMS_RHS_DISPATCH_MOVING
#endif
