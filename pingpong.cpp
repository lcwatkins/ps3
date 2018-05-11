#include <omp.h>
#include <mpi.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <assert.h>


int main(int argc, char *argv[]){
    /////////////// Initialize MPI /////////////////
    int nprocs, stat, myrank, ngrid, subgridlen;
    MPI_Init(&argc, &argv);
    stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert(stat == MPI_SUCCESS);
    stat = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    assert(stat == MPI_SUCCESS);


    ///////////////   Ping-pong   /////////////////
    int j;
    int npings = 100000;
    double t0, t1, totaltime, localtime;
    MPI_Status status;
    t0 = omp_get_wtime();
    for (int i=0;i<npings;i++) {
        if (myrank == 0){
            //printf("sending %d\n",myrank);
            MPI_Send(&i, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            //printf("receiving %d\n",myrank);
            MPI_Recv(&j, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
            //printf("send %d, recv %d\n",i,j);
        } else {
            //printf("receiving %d\n",myrank);
            MPI_Recv(&j, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            //printf("sending %d\n",myrank);
            MPI_Send(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    t1 = omp_get_wtime();
    localtime = t1 - t0;
    MPI_Reduce(&localtime, &totaltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Finalize();

    if (myrank == 0){
        printf("Average time (ms): %f\n", (totaltime/npings)*1000000);
    }
    return 0;
}

