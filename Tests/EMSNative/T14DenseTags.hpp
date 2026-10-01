#ifndef T14_DENSE_TAGS_HPP_
#define T14_DENSE_TAGS_HPP_
#include "EMSBH2DLevel.hpp"
#include "DefaultLevelFactory.hpp"
#include "DenseIntVectSet.H"
#include "BoxIterator.H"

// Tests-only representation workaround for TreeIntVectSet's fixed 24-entry
// traversal stacks at large global indices. The selected cells and dilation
// are identical; all levels below 13 call the original tagging method.
// Compile this test harness with -fno-access-control to call the existing
// private constructor/tagging methods, without editing production headers.
class T14DenseTagsLevel : public EMSBH2DLevel
{
  public:
    T14DenseTagsLevel(GRAMR &g,const SimulationParameters &p,int v)
        :EMSBH2DLevel(g,p,v) {}
    void tagCells(IntVectSet &tags) override
    {
        if(m_level<13) { GRAMRLevel::tagCells(tags);return; }
        fillAllEvolutionGhosts();
        const auto &layout=m_state_new.disjointBoxLayout();
        Box bounds;
        for(DataIterator it=m_state_new.dataIterator();it.ok();++it) bounds.minBox(layout[it()]);
        DenseIntVectSet bits(bounds);bits.makeEmptyBits();
        IntVectSet local(bits);
        for(DataIterator it=m_state_new.dataIterator();it.ok();++it)
        {
            const auto &b=layout[it()];FArrayBox criterion(b,1),invalid;
            const auto &diagnostics=NUM_DIAGNOSTIC_VARS>0?m_state_diagnostics[it()]:invalid;
            GRAMRLevel::computeTaggingCriterion(criterion,m_state_new[it()],diagnostics);
            for(BoxIterator cell(b);cell.ok();++cell)
                if(criterion(cell(),0)>=m_p.regrid_thresholds[m_level]) local|=cell();
        }
        local.grow(m_p.tag_buffer_size);
        local.recalcMinBox();
        Box clipped=local.minBox();clipped&=m_problem_domain;local&=clipped;
        if(!local.isDense()) MayDay::Error("T14 initialization tags unexpectedly became a tree");
        tags=local;
        pout()<<"T14 dense initialization tags level "<<m_level<<" cells "<<local.numPts()<<std::endl;
    }
};

template<class level_t> class T14LevelFactory : public DefaultLevelFactory<level_t>
{
    bool enabled=false;
  public:
    T14LevelFactory(GRAMR &g,SimulationParameters &p):DefaultLevelFactory<level_t>(g,p)
    { GRParmParse pp;pp.load("t14_dense_initial_tags",enabled,false); }
    AMRLevel *new_amrlevel() const override
    {
        if(!enabled) return DefaultLevelFactory<level_t>::new_amrlevel();
        auto *level=new T14DenseTagsLevel(this->m_gr_amr,this->m_p,this->m_p.verbosity);
        level->initialDtMultiplier(this->m_p.dt_multiplier);return level;
    }
};
#endif
