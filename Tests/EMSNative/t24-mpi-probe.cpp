#include <mpi.h>
#include <omp.h>
#include <cstdio>

// Tests only: qualify the local launcher before the checkpoint replay.
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int sum = 0;
    MPI_Allreduce(&rank, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    std::printf("T24_MPI_PROBE rank=%d size=%d sum=%d threads=%d\n",
                rank, size, sum, omp_get_max_threads());
    const int ok = sum == size * (size - 1) / 2;
    MPI_Finalize();
    return ok ? 0 : 1;
}
