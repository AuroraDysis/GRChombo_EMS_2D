#include "BoxLoops.hpp"
#include "SafeNanAbort.hpp"
#include <cstring>
#include <limits>
#include <omp.h>
int main(int argc,char **argv)
{
#ifdef CH_MPI
    MPI_Init(&argc,&argv);
#endif
    if (argc!=4) return 2;
    omp_set_dynamic(0); omp_set_num_threads(2);
    const bool fail=std::string(argv[1])=="nan";
    const int bad_rank=std::stoi(argv[3]);
    const Box box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(1,1,0)));
    FArrayBox data(box,NUM_VARS),before(box,NUM_VARS);
    data.setVal(0.);
    if (fail && procID()==bad_rank)
        data(IntVect(D_DECL(0,1,0)),c_K)=std::numeric_limits<double>::quiet_NaN();
    before.copy(data);
    SafeNanAbort check(1e20);
    BoxLoops::loop(check.collector(),data,data,disable_simd());
    const bool identical=std::memcmp(before.dataPtr(),data.dataPtr(),NUM_VARS*box.numPts()*sizeof(double))==0;
    std::ofstream bits(std::string(argv[2])+".bits.rank"+std::to_string(procID())+".json");
    bits<<"{\"bit_identical\":"<<(identical?"true":"false")
        <<",\"values\":"<<NUM_VARS*box.numPts()<<",\"failed\":"<<(check.failed()?"true":"false")<<"}\n";
    bits.close();
    if (!identical || check.failed()!=(fail && procID()==bad_rank)) return 3;
    check.terminate_if_failed(argv[2],12,135.174,.00042724609375,.000213623046875,86);
#ifdef CH_MPI
    // The healthy rank is deliberately blocked here when its peer aborts.
    MPI_Barrier(Chombo_MPI::comm);
    MPI_Finalize();
#endif
    return 0;
}
