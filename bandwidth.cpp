
#include <omp.h>
#include <mpi.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
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

    std::string testsize;
    int datasize;
    testsize = argv[1];
    if (testsize.compare("KB") == 0){
        datasize = 1024/4;
    } else if (testsize.compare("MB") == 0){
        datasize = 1024*1024/4;
    } else if (testsize.compare("GB") == 0){
        datasize = 1024*1024*1024/4;
    }

    ///////////////   Ping-pong   /////////////////
    int npings = 10;
    double t0, t1, totaltime, localtime;
    MPI_Status status;
    int* gb;
    gb = new int[datasize];
    for (int j=0;j<datasize;j++){
        gb[j] = 9;
    }
    t0 = omp_get_wtime();
    for (int i=0;i<npings;i++) {
        if (myrank == 0){
            //printf("sending %d\n",myrank);
            MPI_Send(&(gb[0]), datasize, MPI_INT, 1, 0, MPI_COMM_WORLD);
            //printf("receiving %d\n",myrank);
        } else {
            //printf("receiving %d\n",myrank);
            stat = MPI_Recv(&(gb[0]), datasize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            assert(stat == MPI_SUCCESS);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = omp_get_wtime();
    localtime = t1 - t0;
    MPI_Reduce(&localtime, &totaltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double nmb = (datasize*4.0)/(1024.0*1024.0);
    if (myrank == 0){
        std::cout << "Average speed for 1 " << testsize << ": " ;
        std::cout << ((nmb/totaltime)/npings) << " (MB/sec)" << std::endl;
    }
    MPI_Finalize();

    return 0;
}










