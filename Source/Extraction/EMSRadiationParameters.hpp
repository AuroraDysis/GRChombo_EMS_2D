#ifndef EMS_RADIATION_PARAMETERS_HPP
#define EMS_RADIATION_PARAMETERS_HPP
#include "GRParmParse.hpp"
#include "SimulationParameters.hpp"
#include <cmath>
#include <stdexcept>

inline double ems_rh_expansion_threshold(GRParmParse &pp)
{
    double value=1e-7;
    pp.query("ems_rh_expansion_threshold",value);
    if (!std::isfinite(value) || value<=0 || value>1e-7)
        throw std::runtime_error("ems_rh_expansion_threshold must be in (0,1e-7]");
    return value;
}

struct EMSRadiationParameters
{
    bool active=false;
    double phi_inf=0;
    int wave_level=0;
    double wave_radius=0;
    EMSRadiationParameters(GRParmParse &pp,const SimulationParameters &p)
    {
        pp.query("ems_radiation_activate",active);
        if (!active) return;
        pp.query("ems_radiation_phi_inf",phi_inf);
        pp.query("ems_radiation_wave_level",wave_level);
        pp.query("ems_radiation_wave_radius",wave_radius);
        if (!p.activate_extraction)
            throw std::runtime_error("EMS radiation requires activate_extraction=1 (shared spheres)");
        if (!std::isfinite(phi_inf) || wave_level<0 || wave_level>p.max_level ||
            !std::isfinite(wave_radius) || wave_radius<0 ||
            (wave_level>0 && wave_radius==0) || p.extraction_params.center[1]!=0 ||
            p.extraction_params.num_points_theta<17 || p.extraction_params.num_points_theta%2==0)
            throw std::runtime_error("invalid EMS radiation parameters (axis y=0, odd theta>=17)");
        for (int l:p.extraction_params.extraction_levels)
            if (l<0 || l>p.max_level)
                throw std::runtime_error("EMS radiation extraction level outside hierarchy");
        for (double r:p.extraction_params.extraction_radii)
            if (!std::isfinite(r) || r<=0)
                throw std::runtime_error("EMS radiation radii must be finite and positive");
        const auto c=p.coupling_function_params;
        const double finf=std::exp(-2*c.alpha*(c.f0+c.f1*phi_inf+c.f2*phi_inf*phi_inf));
        if (!std::isfinite(finf) || finf<=0)
            throw std::runtime_error("EMS radiation F_inf must be finite and positive");
    }
};
#endif
