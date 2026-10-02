// Geometry-only replay of native radius tags and the actual Chombo mesher/LB.
// No evolved state is allocated; output includes every box and rank critical path.
#include "BRMeshRefine.H"
#include "DenseIntVectSet.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include <cmath>
#include <fstream>
#include <iostream>
#include "UsingNamespace.H"

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    std::ofstream boxes(std::string(argv[1])+"-boxes.csv");
    std::ofstream loads(std::string(argv[1])+"-loads.csv");
    boxes<<"max_box_size,level,box,x0,y0,x1,y1,valid_cells,halo_cells\n";
    loads<<"max_box_size,level,ranks,boxes,valid_cells,halo_cells,max_rank_valid,max_rank_halo,max_rank_weight_quarter,active_ranks\n";
    const int finest=12,block=8;
    const Box domain(IntVect(0,0),IntVect(2559,1279));
    Vector<int> ratios(finest+1,2);
    for (int size : {24})
    {
        Vector<Vector<Box>> old(finest+1), grid;
        // AMR::makeBaseLevelMesh: split the block-coarsened domain evenly.
        const int nx=(320+size/block-1)/(size/block);
        const int ny=(160+size/block-1)/(size/block);
        const int sx=(320+nx-1)/nx,sy=(160+ny-1)/ny;
        for(int j=0;j<ny;++j)for(int i=0;i<nx;++i)
            old[0].push_back(refine(Box(IntVect(i*sx,j*sy),
                IntVect(std::min((i+1)*sx,320)-1,std::min((j+1)*sy,160)-1)),block));
        Vector<IntVectSet> tags(finest);
        for(int l=0;l<finest;++l)
        {
            const double h=1.75/std::pow(2,l);
            const double face=l==0?224.:l==1?144.:224./std::pow(2,l+1);
            // Preserve the native 1.2 multiplication of the parameter radius.
            const double radius=1.2*(face/1.2);
            Box bounds(IntVect(std::floor((2232-radius)/h-.5)-4,0),
                       IntVect(std::ceil((2248+radius)/h-.5)+4,std::ceil(radius/h)+4));
            DenseIntVectSet bits(bounds);bits.makeEmptyBits();
            IntVectSet set(bits);
            for(double centre:{2232.,2248.})
            {
                const int lo=std::floor((centre-radius)/h-.5);
                const int hi=std::ceil((centre+radius)/h-.5);
                for(int j=0;j<=std::ceil(radius/h);++j)for(int i=lo;i<=hi;++i)
                    if(std::hypot((i+.5)*h-centre,(j+.5)*h)<radius)
                        set|=IntVect(i,j);
            }
            set.grow(3);set.recalcMinBox();set&=refine(domain,1<<l);
            tags[l]=set;
        }
        BRMeshRefine mesher(domain,ratios,.7,block,8,size);
        const int achieved=mesher.regrid(grid,tags,0,finest-1,old);
        if(achieved!=finest)return 3;
        grid[0]=old[0];
        for(int l=0;l<=finest;++l)
        {
            long long cells=0,halo=0;
            for(int k=0;k<grid[l].size();++k)
            {
                const Box &b=grid[l][k];const auto v=b.numPts(),g=grow(b,3).numPts();
                cells+=v;halo+=g;
                boxes<<size<<','<<l<<','<<k<<','<<b.smallEnd(0)<<','<<b.smallEnd(1)<<','
                     <<b.bigEnd(0)<<','<<b.bigEnd(1)<<','<<v<<','<<g<<'\n';
            }
            for(int ranks:{32,64,128,256,512})
            {
                Vector<int> assignment;
                if(LoadBalance(assignment,grid[l],ranks))return 4;
                Vector<long long> valid(ranks,0),ghost(ranks,0);
                for(int k=0;k<grid[l].size();++k)
                {valid[assignment[k]]+=grid[l][k].numPts();ghost[assignment[k]]+=grow(grid[l][k],3).numPts();}
                int active=0;long long mv=0,mg=0,mq=0;
                for(int p=0;p<ranks;++p){active+=valid[p]>0;mv=std::max(mv,valid[p]);mg=std::max(mg,ghost[p]);mq=std::max(mq,3*valid[p]+ghost[p]);}
                loads<<size<<','<<l<<','<<ranks<<','<<grid[l].size()<<','<<cells<<','<<halo<<','<<mv<<','<<mg<<','<<mq/4.<<','<<active<<'\n';
            }
        }
        std::cout<<"completed max_box_size="<<size<<std::endl;
    }
    return (!boxes||!loads)?5:0;
}
