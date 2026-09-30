#ifndef T7OPERATIONRECORDER_HPP_
#define T7OPERATIONRECORDER_HPP_
#include "GRLevelData.hpp"
#include "GRParmParse.hpp"
#include "BoxIterator.H"
#include <cstdio>
#include <map>
#include <vector>
#include <cstring>
#include <cstdint>
#include "IntVectSet.H"

// Serial diagnostics, lossless Float64 bit-XOR between matching sample frames.
class T7OperationRecorder
{
    std::string m_path;
    std::map<int,FILE *> m_files;
    std::map<std::vector<int>,std::vector<std::uint64_t>> m_previous;
    void bytes(FILE *f,const void *p,std::size_t n)
    { if (std::fwrite(p,1,n,f)!=n) MayDay::Error("T7 diagnostic write failed"); }
  public:
    bool enabled=false;
    double stage_time=0.,coarse_old=0.,coarse_new=0.,step_fraction=0.,update_dt=0.;
    int stage=-1;
    explicit T7OperationRecorder(const std::string &path):m_path(path)
    { GRParmParse pp;pp.load("t7_diagnostics",enabled,false); }
    ~T7OperationRecorder()
    { for (auto &p:m_files) if (::pclose(p.second)!=0) MayDay::Error("T7 diagnostic compressor failed"); }
    static std::vector<Box> windows(double h,double centre,double a,double b)
    {
        std::vector<Box> out;
        for (double f:{a,b})
        {
            if (f<=0. || f>4.+1e-10) continue;
            int c=std::llround(centre/h),k=std::llround(f/h);
            out.emplace_back(IntVect(D_DECL(c+k-1,0,0)),IntVect(D_DECL(c+k,0,0)));
            out.emplace_back(IntVect(D_DECL(c-1,k-1,0)),IntVect(D_DECL(c-1,k,0)));
            out.emplace_back(IntVect(D_DECL(c+k-1,k-1,0)),IntVect(D_DECL(c+k,k,0)));
        }
        return out;
    }
    static IntVectSet dependencies(const Box &core)
    {
        IntVectSet cells;
        for (BoxIterator bit(core);bit.ok();++bit)
        {
            // All mixed second derivatives plus axial d1/d2, upwind and KO.
            for (int y=-2;y<=2;++y) for (int x=-2;x<=2;++x)
                cells|=bit()+IntVect(D_DECL(x,y,0));
            for (int k:{-3,3}) for (int dir=0;dir<2;++dir)
                cells|=bit()+k*BASISV(dir);
        }
        return cells;
    }
    void frame(int phase,int level,int source,const FArrayBox &fab,const Box &valid,
               const std::vector<IntVect> &cells,int first,int n,double time,double h,double dt,
               double a,double b,bool snapshot=false)
    {
        if (cells.empty()) return;
#ifdef CH_MPI
        MayDay::Error("t7_diagnostics requires local serial Chombo");
#endif
        int stream=level+(snapshot?100:0);
        if (!m_files.count(stream))
        {
            const auto path=m_path+(snapshot?"t7-snapshot-L":"t7-stage-L")+std::to_string(level)+".xz";
            std::string quoted="'";
            for (char c:path) quoted+=(c=='\''?"'\\''":std::string(1,c));
            quoted+="'";
            // Installed xz, one thread per stream; lossless, opt-in serial only.
            auto f=::popen(("xz -c -3 -T1 > "+quoted).c_str(),"w");
            if (!f) MayDay::Error("T7 diagnostic open failed");
            m_files[stream]=f;bytes(f,"T7OP0002",8);
        }
        auto f=m_files[stream];
        std::vector<int> key{stream,source,n,phase==40?40:phase>=20&&phase<40?phase:0};
        for (const auto &iv:cells) {key.push_back(iv[0]);key.push_back(iv[1]);}
        std::vector<std::uint64_t> bits(cells.size()*n),encoded(bits.size());
        for (size_t i=0;i<cells.size();++i) for (int c=0;c<n;++c)
        { double v=fab(cells[i],first+c);std::memcpy(&bits[i*n+c],&v,8); }
        auto &previous=m_previous[key];bool delta=previous.size()==bits.size();
        for (size_t i=0;i<bits.size();++i) encoded[i]=bits[i]^(delta?previous[i]:0);
        previous=bits;
        const std::int32_t header[]={phase,level,source,stage,n,int(cells.size()),int(delta),
            valid.smallEnd(0),valid.smallEnd(1),valid.bigEnd(0),valid.bigEnd(1)};
        const double meta[]={time,h,dt,stage_time,coarse_old,coarse_new,step_fraction,update_dt,a,b};
        bytes(f,header,sizeof(header));bytes(f,meta,sizeof(meta));
        for (const auto &iv:cells) {const std::int32_t p[]={iv[0],iv[1]};bytes(f,p,8);}
        std::vector<unsigned char> shuffled(encoded.size()*8);
        for (size_t i=0;i<encoded.size();++i) for (int j=0;j<8;++j)
            shuffled[j*encoded.size()+i]=reinterpret_cast<unsigned char*>(&encoded[i])[j];
        bytes(f,shuffled.data(),shuffled.size());
    }
    void record(int phase,int level,const GRLevelData &data,const std::vector<Box> &regions,
                int grow_cells,int first,int n,double time,double h,double dt,double a,double b)
    {
        int source=0;
        for (DataIterator it=data.dataIterator();it.ok();++it,++source)
        {
            const Box valid=data.disjointBoxLayout()[it()];int region=0;
            for (const auto &window:regions)
            {
                Box sample=window&valid;
                if (!sample.isEmpty())
                {
                    std::vector<IntVect> cells;
                    if (phase==2 || phase==3)
                    {
                        auto needed=dependencies(sample);
                        for (IVSIterator iv(needed);iv.ok();++iv) cells.push_back(iv());
                    }
                    else
                    {
                        sample.grow(grow_cells);sample&=data[it()].box();
                        for (BoxIterator bit(sample);bit.ok();++bit) cells.push_back(bit());
                    }
                    frame(phase,level,source*16+region,data[it()],valid,cells,first,n,time,h,dt,a,b);
                }
                ++region;
            }
        }
    }
};
#endif
