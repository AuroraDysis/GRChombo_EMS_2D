// T25: frozen, serial native checkpoint expansion map. No initial-data/evolution calls.
// CLI: ./t25-expansion-map.ex checkpoint-params.txt map_surfaces=surfaces.csv
//      map_output=map.csv map_angles=angles.csv map_punctures='x1 x2' map_points=96
// Surface CSV (no commas in labels): id,family,centre,a,c. Family identifies nested
// surfaces at ONE centre/aspect ratio; a is cylindrical and c is axial semiaxis.
// Cluster build (from Tests/EMSNative, with SERIAL Chombo/HDF5 libraries):
// UCPH: sh t25-build-node.sh ../.. "$CHOMBO_HOME" "$TMPDIR/t25-build" CXX=g++ OPT=HIGH
// Leonardo (GCC 11.5.0 loaded): same command, CXX=g++. Site HDF5 flags in
// Chombo's Make.defs.local. MPI=FALSE, OPENMPCC=FALSE, Float64, -j1 are forced.
// No fast-math; -ffp-contract=off. Run the executable alone, not under mpirun.
// Then: python3 t25-brackets.py map.csv brackets.csv (standard library only).
// See t25-build.json for the exact local command and dependency hashes.
#include "GRAMR.hpp"
#include "SetupFunctions.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "RHUnion.hpp"
#include <chrono>
#include <limits>
#include <map>
#include <stdexcept>

class FrozenLevel : public GRAMRLevel
{
  public:
    static bool finder_mode;
    using GRAMRLevel::GRAMRLevel;
    void initialData() override { throw std::runtime_error("checkpoint required"); }
    void specificEvalRHS(GRLevelData &, GRLevelData &, double) override
    { throw std::runtime_error("T25 must never evolve"); }
    void computeTaggingCriterion(FArrayBox &, const FArrayBox &) override
    { throw std::runtime_error("T25 must never regrid"); }
    void readCheckpointLevel(HDF5Handle &h) override
    {
        GRAMRLevel::readCheckpointLevel(h);
        // Frozen interpolation uses only m_state_new, including coarse support.
        // Release unused evolution/diagnostic buffers immediately, per level.
        m_state_old.clear();
        if(!finder_mode) m_state_diagnostics.clear();
    }
};
bool FrozenLevel::finder_mode=false;

struct Surface { std::string id, family; double centre, a, c; };
std::vector<Surface> load_surfaces(const std::string &path)
{
    std::ifstream in(path); if (!in) throw std::runtime_error("surface CSV missing");
    std::vector<Surface> out; std::string line;
    std::getline(in,line);
    if(!line.empty() && line.back()=='\r') line.pop_back();
    if (line != "id,family,centre,a,c") throw std::runtime_error("surface CSV header");
    while (std::getline(in,line))
    {
        if (line.empty()) continue;
        std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line);
        Surface s; if (!(row>>s.id>>s.family>>s.centre>>s.a>>s.c) ||
            !std::isfinite(s.centre+s.a+s.c) || s.a<=0 || s.c<=0)
            throw std::runtime_error("invalid surface row");
        out.push_back(s);
    }
    if (out.empty()) throw std::runtime_error("empty surfaces");
    return out;
}

// Outward unit normal from F=(x-xc)^2/c^2+(rho^2+w^2)/a^2-1.
// Actual physical inverse, not the unit-determinant cofactor approximation.
struct Point { double theta, da, lapse, chi, det, finder, dx; int level; };
Point geometry(const RHSurf &s, int ii, const Surface &v, double lapse, double dx, int level)
{
    const double t=s.m_theta[ii], r=s.m_f[ii], rho=r*std::sin(t);
    const double h[3][3]={{s.m_h11[ii],s.m_h12[ii],0},
        {s.m_h12[ii],s.m_h22[ii],0},{0,0,s.m_hww[ii]}};
    const double chi=s.m_chi[ii], d=h[0][0]*h[1][1]-h[0][1]*h[0][1];
    if (!(chi>0 && h[0][0]>0 && d>0 && h[2][2]>0))
        return {NAN,NAN,lapse,chi,d*h[2][2],NAN,dx,level};
    const double gU[3][3]={{chi*h[1][1]/d,-chi*h[0][1]/d,0},
        {-chi*h[0][1]/d,chi*h[0][0]/d,0},{0,0,chi/h[2][2]}};
    double dh[3][3][3]={};
    dh[0][0][0]=s.m_dx_h11[ii]; dh[0][0][1]=dh[0][1][0]=s.m_dx_h12[ii];
    dh[0][1][1]=s.m_dx_h22[ii]; dh[0][2][2]=s.m_dx_hww[ii];
    dh[1][0][0]=s.m_dy_h11[ii]; dh[1][0][1]=dh[1][1][0]=s.m_dy_h12[ii];
    dh[1][1][1]=s.m_dy_h22[ii]; dh[1][2][2]=s.m_dy_hww[ii];
    dh[2][0][2]=dh[2][2][0]=h[0][1]/rho;
    dh[2][1][2]=dh[2][2][1]=(h[1][1]-h[2][2])/rho;
    const double dc[3]={s.m_dx_chi[ii],s.m_dy_chi[ii],0};
    double dg[3][3][3]={};
    for (int k=0;k<3;++k) for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        dg[k][i][j]=dh[k][i][j]/chi-h[i][j]*dc[k]/(chi*chi);
    const double n[3]={2*r*std::cos(t)/(v.c*v.c),2*rho/(v.a*v.a),0};
    const double f2[3]={2/(v.c*v.c),2/(v.a*v.a),2/(v.a*v.a)};
    double norm2=0, su[3]={};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) norm2+=gU[i][j]*n[i]*n[j];
    const double norm=std::sqrt(norm2);
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) su[i]+=gU[i][j]*n[j]/norm;
    double div=0;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
    {
        double hf=i==j ? f2[i] : 0;
        for (int k=0;k<3;++k) for (int l=0;l<3;++l)
            hf-=.5*gU[k][l]*(dg[i][l][j]+dg[j][l][i]-dg[l][i][j])*n[k];
        div+=(gU[i][j]-su[i]*su[j])*hf/norm;
    }
    const double ass=(s.m_A11[ii]*su[0]*su[0]+2*s.m_A12[ii]*su[0]*su[1]+
        s.m_A22[ii]*su[1]*su[1])/chi;
    // Exact analytic tangent of the radial graph of this spheroid.
    const double denom=std::cos(t)*std::cos(t)/(v.c*v.c)+std::sin(t)*std::sin(t)/(v.a*v.a);
    const double rp=-r*std::sin(t)*std::cos(t)*(1/(v.a*v.a)-1/(v.c*v.c))/denom;
    const double tx=rp*std::cos(t)-r*std::sin(t), ty=rp*std::sin(t)+r*std::cos(t);
    const double da=2*M_PI*rho/chi*std::sqrt((h[0][0]*tx*tx+2*h[0][1]*tx*ty+h[1][1]*ty*ty)*h[2][2])*s.m_d_theta;
    return {div+ass-2*s.m_K[ii]/3,da,lapse,chi,d*h[2][2],s.Theta_plus(ii),dx,level};
}

// One bounded search, seeded by an externally recorded resolved map bracket.
// RHUnion::update and RHSurf chase/recentring/Newton are unchanged.
int find_from_seed(AMRInterpolator<Lagrange<4>> &interp, const SimulationParameters &p,
                  const Surface &v, int n, double time, GRParmParse &pp)
{
    std::string bracket; pp.load("map_finder_bracket",bracket);
    std::ifstream proof(bracket);
    if(!proof) throw std::runtime_error("finder requires recorded bracket CSV");
    std::string header,row;std::getline(proof,header);std::getline(proof,row);
    if(header.find("status")==std::string::npos || row.find("POINTWISE")==std::string::npos)
        throw std::runtime_error("finder requires sampled pointwise barrier (not an average root)");
    std::map<std::string,std::string> values;std::istringstream hs(header),rs(row);std::string key,value;
    while(std::getline(hs,key,','))
    {
        if(!std::getline(rs,value,',')) throw std::runtime_error("malformed bracket proof");
        values[key]=value;
    }
    if(values.at("checkpoint")!=p.restart_file || std::stod(values.at("time"))!=time ||
       std::stod(values.at("centre"))!=v.centre ||
       v.a<std::stod(values.at("a_inner")) || v.a>std::stod(values.at("a_outer")) ||
       v.c<std::stod(values.at("c_inner")) || v.c>std::stod(values.at("c_outer")) ||
       !(std::stod(values.at("theta_inner_max"))<0 && std::stod(values.at("theta_outer_min"))>0))
        throw std::runtime_error("bracket does not match checkpoint/time/centre/seed or signs");
    int max_updates,window; double seconds_cap,chase;std::vector<double> thresholds;
    pp.load("map_find_max_updates",max_updates,100000);
    pp.load("map_find_floor_window",window,100001);
    pp.load("map_find_seconds",seconds_cap,235.);
    pp.load("map_find_chase",chase,1.);
    pp.load("map_find_thresholds",thresholds,3,std::vector<double>{1e-7,1e-10,1e-12});
    if(max_updates<1 || window<4 || seconds_cap<=0 || seconds_cap>1795 || chase<=0)
        throw std::runtime_error("invalid finder limits");
    RHUnion rh; rh.set_interpolator(&interp);
    rh.setup(1,{v.a},{v.centre},{n},{0},{1},{0.},{chase},{0.},0.);
    auto &s=rh.m_surfaces[0];
    for(int i=0;i<n+2*s.m_NG;++i)
        s.m_f[i]=1/std::sqrt(std::pow(std::cos(s.m_theta[i])/v.c,2)+std::pow(std::sin(s.m_theta[i])/v.a,2));
    const auto cp=p.coupling_function_params;
    rh.set_coupling_params(cp.alpha,cp.f0,cp.f1,cp.f2);
    std::ofstream out("finder.csv"),trace("finder-progress.csv");
    out<<"time,N,stage,threshold,status,updates,seconds,expansion_squared,A,Q,centre,r_min,r_max\n"<<std::setprecision(17);
    trace<<"stage,update,seconds,expansion_squared\n"<<std::setprecision(17);
    const auto start=std::chrono::steady_clock::now();
    auto seconds=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    interp.refresh();rh.interpolate_fields();
    for(size_t stage=0;stage<thresholds.size();++stage)
    {
        rh.m_thresh_super_low=thresholds[stage];std::vector<double> errors;
        double previous_min=0,previous_max=0;int updates=0;std::string status="UPDATE_CAP";
        for(;updates<max_updates;)
        {
            rh.update(time,0);++updates;rh.interpolate_fields();
            const double error=s.expansion_error();
            if(s.m_dead || !std::isfinite(error)) throw std::runtime_error("finder surface died");
            if(updates==1 || updates%100==0) trace<<stage<<','<<updates<<','<<seconds()<<','<<error<<std::endl;
            if(error<=thresholds[stage]) {status="FOUND";break;}
            if(seconds()>=seconds_cap) {status="TIME_CAP";break;}
            errors.push_back(error);
            if(errors.size()==size_t(window))
            {
                const auto range=std::minmax_element(errors.begin(),errors.end());
                if(std::abs(*range.first-previous_min)<.01**range.first &&
                   std::abs(*range.second-previous_max)<.01**range.second)
                {status="FLOOR";break;}
                previous_min=*range.first;previous_max=*range.second;errors.clear();
            }
        }
        const auto begin=s.m_f.begin()+s.m_NG,end=begin+n;
        const auto range=std::minmax_element(begin,end);
        out<<time<<','<<n<<','<<stage<<','<<thresholds[stage]<<','<<status<<','<<updates<<','<<seconds()<<','<<s.expansion_error()<<','<<s.Area()<<','<<s.Q_charge()<<','<<s.m_centre[0]<<','<<*range.first<<','<<*range.second<<std::endl;
        std::ofstream shape("finder-shape-"+std::to_string(stage)+".csv");
        shape<<"theta,radius\n"<<std::setprecision(17);
        for(int i=0;i<n;++i) shape<<s.m_theta[i+s.m_NG]<<','<<s.m_f[i+s.m_NG]<<'\n';
        if(status!="FOUND") {pout()<<"T25_FINDER_COMPLETE "<<status<<"; no advances\n";return 1;}
    }
    pout()<<"T25_FINDER_COMPLETE FOUND_ALL_STAGES; no advances\n";return 0;
}

int run(int argc,char **argv)
{
    GRParmParse pp(argc-2,argv+2,nullptr,argv[1]); SimulationParameters p(pp);
    pp.load("map_find",FrozenLevel::finder_mode,false);
    if (numProc()!=1 || !p.restart_from_checkpoint) throw std::runtime_error("serial checkpoint only");
    std::string input, output, angles;
    int n; std::vector<double> punctures;
    pp.load("map_surfaces",input); pp.load("map_output",output,std::string("map.csv"));
    pp.load("map_angles",angles,std::string("angles.csv")); pp.load("map_points",n,96);
    if(pp.contains("map_punctures")) pp.load("map_punctures",punctures,pp.countval("map_punctures"));
    if(n<8 || n>4096) throw std::runtime_error("map_points outside 8..4096");
    const auto specs=load_surfaces(input);
    GRAMR amr; DefaultLevelFactory<FrozenLevel> factory(amr,p); setupAMRObject(amr,factory);
    const auto levels=amr.get_gramrlevels(); const double time=levels[0]->time();
    if(levels[0]->get_dx()!=p.dx[0] || levels[0]->problemDomain().domainBox().bigEnd()!=p.ivN)
        throw std::runtime_error("parameter grid differs from checkpoint");
    for(const auto *level:levels) if(std::abs(level->time()-time)>1e-10*std::max(1.,std::abs(time)))
        throw std::runtime_error("unsynchronized checkpoint");
    AMRInterpolator<Lagrange<4>> interp(amr,p.origin,p.dx,p.boundary_params,0);
    amr.set_interpolator(&interp);
    interp.refresh(false);
    interp.fill_multilevel_ghosts(VariableType::evolution);
    std::ofstream out(output), detail(angles);
    if(!out || !detail) throw std::runtime_error("output open");
    out<<"checkpoint,time,id,family,centre,a,c,N,status,A,Q,x_min,x_max,rho_max,dx_min,dx_max,min_extent_over_dx,theta_min,theta_max,theta_mean,theta_rms,negative_area_fraction,lapse_mean,chi_mean,det_h_max_error,finder_theta_max_difference";
    for(size_t k=0;k<punctures.size();++k) out<<",encloses_"<<k;
    out<<'\n'<<std::setprecision(17);
    detail<<"id,theta,x,rho,level,dx,expansion,dA,lapse,chi,det_h,finder_expansion\n"<<std::setprecision(17);
    const auto cp=p.coupling_function_params;
    bool all_resolved=true;
    for(size_t k=0;k<specs.size();++k)
    {
        const auto &v=specs[k]; RHUnion rh; rh.set_interpolator(&interp);
        rh.m_surfaces.emplace_back(n,RHUnion::NG,std::vector<double>{v.centre,0},0);
        auto &s=rh.m_surfaces[0];
        for(int i=0;i<n+2*s.m_NG;++i)
            s.m_f[i]=1/std::sqrt(std::pow(std::cos(s.m_theta[i])/v.c,2)+std::pow(std::sin(s.m_theta[i])/v.a,2));
        rh.set_coupling_params(cp.alpha,cp.f0,cp.f1,cp.f2); rh.interpolate_fields();
        std::vector<double> x(n),y(n),lapse(n);
        for(int i=0;i<n;++i) {int ii=i+s.m_NG;x[i]=v.centre+s.m_f[ii]*std::cos(s.m_theta[ii]); y[i]=s.m_f[ii]*std::sin(s.m_theta[ii]);}
        InterpolationQuery query(n); query.setCoords(0,x.data()).setCoords(1,y.data()).addComp(c_lapse,lapse.data());
        interp.interp(query);
        // Native valid-box ownership, solely to report actual AMR coverage.
        // -fno-access-control is Tests-only. Values/derivatives above are native polynomial interpolation.
        const auto layout=interp.findBoxes(query);
        double A=0,Q=0,mean=0,sq=0,negative=0,ml=0,mc=0;
        double low=INFINITY,high=-INFINITY,dmin=INFINITY,dmax=0,deterr=0,fdiff=0;
        bool valid=true;
        for(int i=0;i<n;++i)
        {
            const int ii=i+s.m_NG, lev=layout.level_idx[i];
            if(lev<0 || lev>=int(levels.size())) throw std::runtime_error("native coverage unavailable");
            const auto z=geometry(s,ii,v,lapse[i],levels[lev]->get_dx(),lev);
            detail<<v.id<<','<<s.m_theta[ii]<<','<<x[i]<<','<<y[i]<<','<<lev<<','<<z.dx<<','<<z.theta<<','<<z.da<<','<<z.lapse<<','<<z.chi<<','<<z.det<<','<<z.finder<<'\n';
            valid &= std::isfinite(z.theta+z.da+z.lapse+z.chi+z.det+z.finder);
            low=std::min(low,z.theta);high=std::max(high,z.theta);dmin=std::min(dmin,z.dx);dmax=std::max(dmax,z.dx);
            deterr=std::max(deterr,std::abs(z.det-1));fdiff=std::max(fdiff,std::abs(z.theta-z.finder));
            A+=z.da; mean+=z.theta*z.da;sq+=z.theta*z.theta*z.da;
            if(z.theta<0) negative+=z.da;
            ml+=z.lapse*z.da;mc+=z.chi*z.da;
            // Stored E_i; contract with the physical contravariant unit normal.
            const double d=s.m_h11[ii]*s.m_h22[ii]-s.m_h12[ii]*s.m_h12[ii];
            const double nx=2*(x[i]-v.centre)/(v.c*v.c),ny=2*y[i]/(v.a*v.a);
            const double ux=s.m_chi[ii]*(s.m_h22[ii]*nx-s.m_h12[ii]*ny)/d;
            const double uy=s.m_chi[ii]*(s.m_h11[ii]*ny-s.m_h12[ii]*nx)/d;
            const double norm=std::sqrt(nx*ux+ny*uy);
            const double F=std::exp(-2*cp.alpha*(cp.f0+cp.f1*s.m_phi[ii]+cp.f2*s.m_phi[ii]*s.m_phi[ii]));
            Q+=F*(s.m_Ex[ii]*ux+s.m_Ey[ii]*uy)*z.da/norm/std::sqrt(2*M_PI);
        }
        const double cells=std::min(v.a,v.c)/dmax;
        const std::string status=!valid?"INVALID_GEOMETRY":(cells<3?"UNRESOLVED":"RESOLVED");
        all_resolved &= status=="RESOLVED";
        out<<p.restart_file<<','<<time<<','<<v.id<<','<<v.family<<','<<v.centre<<','<<v.a<<','<<v.c<<','<<n<<','<<status<<','<<A<<','<<Q<<','<<v.centre-v.c<<','<<v.centre+v.c<<','<<v.a<<','<<dmin<<','<<dmax<<','<<cells<<','<<low<<','<<high<<','<<mean/A<<','<<std::sqrt(sq/A)<<','<<negative/A<<','<<ml/A<<','<<mc/A<<','<<deterr<<','<<fdiff;
        for(double puncture:punctures) out<<','<<(std::abs(puncture-v.centre)<v.c?1:0);
        out<<std::endl;
        pout()<<"T25_SURFACE "<<v.id<<" "<<status<<" theta=["<<low<<","<<high<<"]\n";
    }
    if(!out || !detail) throw std::runtime_error("output write");
    pout()<<"T25_EXPANSION_MAP_COMPLETE; no advances\n";
    if(FrozenLevel::finder_mode)
    {
        if(specs.size()!=1) throw std::runtime_error("finder takes exactly one bracket-derived seed");
        if(!all_resolved) throw std::runtime_error("finder seed is unresolved or invalid");
        return find_from_seed(interp,p,specs[0],n,time,pp);
    }
    return 0;
}
int main(int argc,char **argv)
{
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    mainSetup(argc,argv); int code=0;
    try {code=run(argc,argv);} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';code=2;}
    mainFinalize(); return code;
}
