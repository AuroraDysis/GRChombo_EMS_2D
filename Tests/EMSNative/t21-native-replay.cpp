// Frozen numerical checkpoint only; uses whichever gauge the production level selects.
#include "EMSBH2DLevel.hpp"
#include "SetupFunctions.hpp"
#include "DefaultLevelFactory.hpp"
#include "BoxIterator.H"
#include <cstring>
#include <fstream>

class T21FrozenLevel:public EMSBH2DLevel
{
  public:
    T21FrozenLevel(GRAMR &g,const SimulationParameters &p,int v):EMSBH2DLevel(g,p,v) {}
    void initialData() override { throw std::runtime_error("T21 replay requires a complete numerical checkpoint"); }
    void dump()
    {
        m_t21.snapshot=true;
        m_t21.checkpoint_source=&m_state_new;
        if (!m_t21.selected(m_level,m_time,m_dt)) return;
        GRLevelData u,rhs;
        u.define(m_state_new.disjointBoxLayout(),NUM_VARS,m_state_new.ghostVect());
        rhs.define(m_state_new.disjointBoxLayout(),NUM_VARS,m_state_new.ghostVect());
        std::vector<std::unique_ptr<FArrayBox>> saved;
        for (DataIterator it=u.dataIterator();it.ok();++it) {
            u[it()].copy(m_state_new[it()]);
            saved.emplace_back(new FArrayBox(m_state_new.disjointBoxLayout()[it()],NUM_VARS));
            saved.back()->copy(m_state_new[it()]);
        }
        specificEvalRHS(u,rhs,m_time);
        long long checked=0,mismatches=0,nonfinite=0;int box=0;
        // RHS projection acts on a private copy. The checkpoint's valid bits are preserved.
        for (DataIterator it=u.dataIterator();it.ok();++it,++box)
            for (BoxIterator b(m_state_new.disjointBoxLayout()[it()]);b.ok();++b)
                for(int c=0;c<NUM_VARS;++c) {
                    ++checked;
                    double before=(*saved[box])(b(),c),after=m_state_new[it()](b(),c);
                    mismatches+=std::memcmp(&before,&after,8)!=0;
                    nonfinite+=!std::isfinite(after);
                }
        std::ofstream out("replay-valid-L"+std::to_string(m_level)+"-r"+std::to_string(procID())+".csv");
        out<<"level,valid_values,bit_mismatches,nonfinite,advances,static_reads\n"<<m_level<<','<<checked<<','<<mismatches<<','<<nonfinite<<",0,0\n";
        if (!out || mismatches || nonfinite) throw std::runtime_error("T21 invalid replay state");
    }
};

int main(int argc,char **argv)
{
    mainSetup(argc,argv);int status=0;
    try {
        GRParmParse pp(argc-2,argv+2,nullptr,argv[1]);SimulationParameters p(pp);
        bool capture=false;pp.load("t21_rhs_capture",capture,false);
        if (!p.restart_from_checkpoint || !p.t2_guard_initial_data_after_t0 || !capture)
            throw std::runtime_error("T21 requires guarded checkpoint restart and recording config");
        BHAMR amr;DefaultLevelFactory<T21FrozenLevel> factory(amr,p);setupAMRObject(amr,factory);
        auto all=amr.get_gramrlevels();double time=all[0]->time();
        for (auto *l:all) {
            if (std::abs(l->time()-time)>1e-9) throw std::runtime_error("T21 unsynchronized hierarchy");
            l->fillAllEvolutionGhosts();
        }
        for(auto *l:all) dynamic_cast<T21FrozenLevel*>(l)->dump();
    } catch(const std::exception &e) { std::cerr<<e.what()<<'\n';status=2; }
    mainFinalize();return status;
}
