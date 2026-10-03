#ifndef EMS_T21_RECORDER_HPP
#define EMS_T21_RECORDER_HPP
// Diagnostic storage only. No gauge type, static profile or evolution callback.
#include "GRParmParse.hpp"
#include "BoxIterator.H"
#include "SPMD.H"
#include "FourthOrderDerivatives.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>

class T21RHSRecorder
{
    bool enabled=false;
    int first_steps=0, first_cells=0, strip_cells=0;
    double radius=0., profile_radius=0.;
    std::vector<int> levels, fields, parities;
    std::vector<double> directions;
    std::string prefix, gauge, semantics, kernel_type;
    std::uint32_t calls=0;
    bool metadata_written=false;

    bool selected_cell(const IntVect &iv,double h,const std::array<double,2> &centre) const
    {
        double x=(iv[0]+.5)*h-centre[0],y=(iv[1]+.5)*h-centre[1],r=std::hypot(x,y);
        if (!snapshot) return r<=radius || (std::abs(x)<=first_cells*h && std::abs(y)<=first_cells*h);
        if (r>profile_radius+strip_cells*h) return false;
        for (std::size_t i=0;i<directions.size();i+=2)
            if (x*directions[i]+y*directions[i+1]>=-strip_cells*h
                && std::abs(-x*directions[i+1]+y*directions[i])<=strip_cells*h) return true;
        return false;
    }

    static double areal_proxy(const FArrayBox &u, const IntVect &iv,
                              double h, const std::array<double,2> &centre)
    {
        double x=(iv[0]+.5)*h-centre[0], y=(iv[1]+.5)*h-centre[1];
        double r=std::hypot(x,y);
        if (!(r>0.) || !(u(iv,c_chi)>0.)) return NAN;
        double tx=-y/r, ty=x/r;
        double ht=tx*tx*u(iv,c_h11)+2*tx*ty*u(iv,c_h12)+ty*ty*u(iv,c_h22);
        double prod=ht*u(iv,c_hww);
        return ht>0. && u(iv,c_hww)>0. ? r*std::sqrt(std::sqrt(prod))/std::sqrt(u(iv,c_chi)) : NAN;
    }

    static std::array<double,11> geometry(const FArrayBox &u,
                 const FArrayBox &parts, const FArrayBox &rhs,
                 const IntVect &iv, double h, const std::array<double,2> &centre)
    {
        std::array<double,11> a;
        a.fill(NAN); a[10]=0.;
        double x=(iv[0]+.5)*h-centre[0], y=(iv[1]+.5)*h-centre[1];
        double r=std::hypot(x,y), R=areal_proxy(u,iv,h,centre);
        if (!(r>0.) || !std::isfinite(R)) return a;
        double nx=x/r,ny=y/r,tx=-ny,ty=nx;
        double ht=tx*tx*u(iv,c_h11)+2*tx*ty*u(iv,c_h12)+ty*ty*u(iv,c_h22);
        double det=u(iv,c_h11)*u(iv,c_h22)-u(iv,c_h12)*u(iv,c_h12);
        if (!(det>0.) || !(u(iv,c_h11)>0.) || !(u(iv,c_lapse)>0.)) return a;
        auto time_rate=[&](const FArrayBox &f) {
            double dh=tx*tx*f(iv,c_h11)+2*tx*ty*f(iv,c_h12)+ty*ty*f(iv,c_h22);
            return R*(.25*(dh/ht+f(iv,c_hww)/u(iv,c_hww))-.5*f(iv,c_chi)/u(iv,c_chi));
        };
        a[0]=R; a[1]=time_rate(rhs); a[2]=time_rate(parts);
        FourthOrderDerivatives deriv(h);
        for (int d=0;d<2;++d) {
            std::array<double,5> stencil;
            for (int offset=-2;offset<=2;++offset)
                stencil[offset+2]=areal_proxy(u,iv+offset*BASISV(d),h,centre);
            a[3+d]=deriv.diff1<double>(stencil.data(),2,1);
        }
        double spatial=u(iv,c_chi)*(u(iv,c_h22)*a[3]*a[3]
                     -2*u(iv,c_h12)*a[3]*a[4]+u(iv,c_h11)*a[4]*a[4])/det;
        double advection=u(iv,c_shift1)*a[3]+u(iv,c_shift2)*a[4];
        for (int k=0;k<2;++k) {
            double normal=(a[1+k]-advection)/u(iv,c_lapse);
            a[5+k]=spatial-normal*normal;
        }
        a[7]=.5*R*(1-a[5]);
        a[8]=(ht-u(iv,c_hww))/(ht+u(iv,c_hww));
        a[9]=(nx*tx*u(iv,c_h11)+(nx*ty+ny*tx)*u(iv,c_h12)+ny*ty*u(iv,c_h22))/ht;
        a[10]=std::all_of(a.begin(),a.begin()+10,[](double v){return std::isfinite(v);})?1.:0.;
        return a;
    }

  public:
    // Set only by the checkpoint-only diagnostic executable, never by a run parameter.
    bool snapshot=false;
    const GRLevelData *checkpoint_source=nullptr;
    void set_kernel_type(const std::string &s) { kernel_type=s; }
    T21RHSRecorder()
    {
        GRParmParse pp;
        if (pp.contains("t21_rhs_capture")) pp.load("t21_rhs_capture",enabled);
        if (!enabled) return;
        pp.load("t21_gauge_label",gauge);
        pp.load("t21_driver_semantics",semantics);
        pp.load("t21_rhs_prefix",prefix);
        pp.load("t21_initial_rk_steps",first_steps);
        pp.load("t21_puncture_cells",first_cells);
        pp.load("t21_puncture_radius",radius);
        pp.load("t21_profile_radius",profile_radius);
        pp.load("t21_strip_half_cells",strip_cells);
        pp.load("t21_levels",levels,pp.countval("t21_levels"));
        pp.load("t21_ray_directions",directions,pp.countval("t21_ray_directions"));
        std::vector<std::string> names;
        pp.load("t21_fields",names,pp.countval("t21_fields"));
        pp.load("vars_parity",parities,NUM_VARS);
        for (const auto &s:names) {
            auto it=std::find(UserVariables::variable_names.begin(),UserVariables::variable_names.end(),s);
            if (it==UserVariables::variable_names.end()) throw std::runtime_error("T21 unknown field "+s);
            fields.push_back(it-UserVariables::variable_names.begin());
        }
        if (names.empty() || levels.empty() || directions.empty() || directions.size()%2
            || first_steps<0 || first_cells<0 || radius<0 || !(profile_radius>0.) || strip_cells<6)
            throw std::runtime_error("T21 invalid recorder parameters");
        for (const auto &s:{gauge,semantics})
            if (s.empty() || s.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)
                throw std::runtime_error("T21 metadata label must be a nonempty identifier");
        for (std::size_t i=0;i<directions.size();i+=2)
            if (std::abs(std::hypot(directions[i],directions[i+1])-1.)>1e-12)
                throw std::runtime_error("T21 ray must be a unit vector");
    }
    bool selected(int level,double level_time,double dt) const
    {
        return enabled && std::find(levels.begin(),levels.end(),level)!=levels.end()
            && (snapshot || level_time < first_steps*dt-16*std::numeric_limits<double>::epsilon()*std::abs(level_time));
    }
    void next_call() { ++calls; }
    void before_projection(int level,int stage,double time,double h,
             const std::array<double,2> &centre,const GRLevelData &u,
             double chi_floor,double lapse_floor)
    {
        long long cells=0,cf=0,af=0;double chi=INFINITY,alpha=INFINITY;
        for (DataIterator it=u.dataIterator();it.ok();++it)
            for (BoxIterator b(u.disjointBoxLayout()[it()]);b.ok();++b)
                if (selected_cell(b(),h,centre)) {
                    ++cells;double c=u[it()](b(),c_chi),a=u[it()](b(),c_lapse);
                    chi=std::min(chi,c);alpha=std::min(alpha,a);cf+=c<chi_floor;af+=a<lapse_floor;
                }
        if (!cells) return;
        std::ofstream out(prefix+"-L"+std::to_string(level)+"-r"+std::to_string(procID())+"-floors.csv",std::ios::app);
        if (out.tellp()==0) out<<"level,rank,stage,call,time_M,cells,chi_min_before,lapse_min_before,chi_activations,lapse_activations,chi_floor,lapse_floor\n";
        out<<std::setprecision(17)<<level<<','<<procID()<<','<<stage<<','<<calls+1<<','<<time<<','<<cells<<','<<chi<<','<<alpha<<','<<cf<<','<<af<<','<<chi_floor<<','<<lapse_floor<<'\n';
        if (!out) throw std::runtime_error("T21 preprojection floor census write failed");
    }
    void write(int level,int stage,double time,double level_time,double h,double dt,
               const std::array<double,2> &centre,const Box &valid,
               const FArrayBox &u,const FArrayBox &parts,const FArrayBox &rhs)
    {
        std::vector<IntVect> cells;
        for (BoxIterator b(valid);b.ok();++b)
            if (selected_cell(b(),h,centre)) cells.push_back(b());
        if (cells.empty()) return;
        auto stem=prefix+"-L"+std::to_string(level)+"-r"+std::to_string(procID());
        if (!metadata_written) {
            std::ofstream meta(stem+".json");
            meta<<"{\"gauge\":\""<<gauge<<"\",\"driver_semantics\":\""<<semantics
                <<"\",\"native_kernel_type\":\""<<kernel_type<<"\",\"fields\":[";
            for (std::size_t i=0;i<fields.size();++i) {
                if (i) meta<<',';
                meta<<'"'<<UserVariables::variable_names[fields[i]]<<'"';
            }
            meta<<"],\"parities\":[";
            for (std::size_t i=0;i<fields.size();++i) { if(i)meta<<',';meta<<parities[fields[i]]; }
            meta<<"],\"state_semantics\":\""<<(snapshot?"checkpoint_before_rhs_projection":"RK_input_after_native_projection")
                <<"\",\"geometry\":[\"R_area_density\",\"Rdot_total\",\"Rdot_preKO\",\"dRdx\",\"dRdy\",\"X_total\",\"X_preKO\",\"m_MS_proxy\",\"angular_anisotropy\",\"radial_tangent_cross\",\"geometry_valid\"],\"format\":\"T21RHS01-native-endian\"}\n";
            if (!meta) throw std::runtime_error("T21 metadata write failed");
            metadata_written=true;
        }
        std::ofstream out(stem+".bin",std::ios::binary|std::ios::app);
        out.write("T21RHS01",8);
        const std::uint32_t header[]={std::uint32_t(level),std::uint32_t(procID()),std::uint32_t(stage),calls,
            std::uint32_t(fields.size()),std::uint32_t(cells.size()),std::uint32_t(snapshot)};
        const double clocks[]={time,level_time,h,dt,centre[0],centre[1]};
        out.write(reinterpret_cast<const char*>(header),sizeof(header));
        out.write(reinterpret_cast<const char*>(clocks),sizeof(clocks));
        for (const auto &iv:cells) {
            const std::int32_t key[]={iv[0],iv[1]};
            out.write(reinterpret_cast<const char*>(key),sizeof(key));
            for (int part=0;part<4;++part) for (int c:fields) {
                double v=part==0?u(iv,c):part==1?parts(iv,c):part==2?parts(iv,NUM_VARS+c):rhs(iv,c);
                out.write(reinterpret_cast<const char*>(&v),8);
            }
            auto a=geometry(u,parts,rhs,iv,h,centre);
            out.write(reinterpret_cast<const char*>(a.data()),a.size()*8);
        }
        if (!out) throw std::runtime_error("T21 frame write failed");
    }
};
#endif
