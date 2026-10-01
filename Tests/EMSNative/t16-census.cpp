// Initialization only, using the frozen production level factory and objects.
// No advance, extra static profile evaluation, or production-code changes.
#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "EMSBH2DLevel.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
#include "T14DenseTags.hpp"
#include "BoxIterator.H"
#include <fstream>
#include <iomanip>
#include <cmath>

int main(int argc,char **argv)
{
    if(argc!=3) return 2;
    mainSetup(argc,argv);
    GRParmParse pp(0,nullptr,nullptr,argv[1]);
    SimulationParameters p(pp);
    BHAMR amr;
    T14LevelFactory<EMSBH2DLevel> factory(amr,p);
    setupAMRObject(amr,factory);
    AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,p.verbosity);
    amr.set_interpolator(&interp);
    std::ofstream out(argv[2]);
    out<<"level,box,x0,y0,x1,y1,h_M,valid_cells,ghost_cells,chi_min,lapse_min,chi_activations,lapse_activations,nonfinite,max_shift_speed,max_light_speed,max_lapse_speed\n"<<std::setprecision(17);
    for(int l=0;l<amr.getAMRLevels().size();++l)
    {
        auto *lev=dynamic_cast<EMSBH2DLevel *>(amr.getAMRLevels()[l]);
        // Exactly the production t=0 recorder call, including its floor census.
        lev->ems_t13_initial();
        const auto &data=lev->getLevelData();int source=0;
        for(DataIterator it=data.dataIterator();it.ok();++it,++source)
        {
            const Box valid=data.disjointBoxLayout()[it()];const auto &a=data[it()];
            double chi=INFINITY,alpha=INFINITY,shift=0.,light=0.,lapse=0.;
            long long bad=0,cf=0,af=0;
            for(BoxIterator k(a.box());k.ok();++k)
            {
                const auto iv=k();const double c=a(iv,0),al=a(iv,13);
                chi=std::min(chi,c);alpha=std::min(alpha,al);
                cf+=c<p.min_chi;af+=al<p.min_lapse;
                for(int n=0;n<NUM_VARS;++n) bad+=!std::isfinite(a(iv,n));
                const double d=a(iv,1)*a(iv,3)-a(iv,2)*a(iv,2);
                const double inverse=std::max(1./a(iv,4),
                    (a(iv,1)+a(iv,3)+std::hypot(a(iv,1)-a(iv,3),2*a(iv,2)))/(2*d));
                const double beta=std::hypot(a(iv,14),a(iv,15));
                shift=std::max(shift,beta+std::sqrt(inverse));
                light=std::max(light,beta+al*std::sqrt(c*inverse));
                lapse=std::max(lapse,beta+std::sqrt(1.8*al*c*inverse));
            }
            out<<l<<','<<source<<','<<valid.smallEnd(0)<<','<<valid.smallEnd(1)<<','
               <<valid.bigEnd(0)<<','<<valid.bigEnd(1)<<','<<lev->get_dx()<<','
               <<valid.numPts()<<','<<a.box().numPts()-valid.numPts()<<','<<chi<<','
               <<alpha<<','<<cf<<','<<af<<','<<bad<<','<<shift<<','<<light<<','<<lapse<<'\n';
            if(bad||cf||af||!(chi>0.)||!(alpha>0.)) return 4;
        }
    }
    if(!out) return 5;
    mainFinalize();return 0;
}
