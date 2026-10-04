// Tests-only witness of the pinned NanCheck failure; no production edit.
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include <atomic>
#include <cstring>
#include <fstream>
#include <limits>
#include <omp.h>
#include <signal.h>
#include <unistd.h>

static void timeout(int)
{
    const char msg[] = "T24_PROBE_TIMEOUT_AFTER_NANCHECK_DUMP\n";
    write(STDERR_FILENO, msg, sizeof(msg)-1);
    _Exit(124);
}

struct Collect
{
    std::atomic<bool> *bad;
    void compute(Cell<double> cell) const
    {
        for (int n=0; n<cell.get_num_in_vars(); ++n)
        {
            const double v = cell.load_vars(n);
            if (std::isnan(v) || std::abs(v)>1e20)
                bad->store(true, std::memory_order_relaxed);
        }
    }
};

int main(int argc, char **argv)
{
    if (argc!=3) return 2;
    const std::string mode=argv[1];
    omp_set_dynamic(0);
    omp_set_num_threads(2);
    const Box box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(1,1,0)));
    FArrayBox state(box, NUM_VARS);
    state.setVal(0.);
    state.setVal(1., c_chi);
    state.setVal(1., c_lapse);
    state.setVal(1., c_h11);
    state.setVal(1., c_h22);
    state.setVal(1., c_hww);
    // With static two-row work sharing, y=1 belongs to the non-master worker.
    if (mode!="finite") state(IntVect(D_DECL(0,mode=="legacy-master"?0:1,0)),c_K)=std::numeric_limits<double>::quiet_NaN();
    FArrayBox before(box,NUM_VARS);
    before.copy(state);
    const std::array<double,2> centre{0.,0.};
    signal(SIGALRM,timeout);
    alarm(3);
    if (mode=="legacy-master" || mode=="legacy-worker")
        BoxLoops::loop(NanCheck(1.,centre,"T24_ACTUAL_NANCHECK"),state,state,disable_simd());
    else
    {
        std::atomic<bool> bad{false};
        BoxLoops::loop(Collect{&bad},state,state,disable_simd());
        // All workers have joined here: serial witness/abort, no conditional barrier.
        const bool identical=std::memcmp(state.dataPtr(),before.dataPtr(),
            box.numPts()*NUM_VARS*sizeof(double))==0;
        std::ofstream record(argv[2]);
        record<<"{\"mode\":\""<<mode<<"\",\"threads\":2,\"worker_joined\":true,"
              <<"\"state_bit_identical\":"<<(identical?"true":"false")
              <<",\"nonfinite\":"<<(bad.load()?"true":"false")<<"}\n";
        record.close();
        if (!identical) return 2;
        if (bad.load()) { std::cerr<<"T24_TEST_NONFINITE_ABORT\n"; return 1; }
    }
    alarm(0);
    return 0;
}
