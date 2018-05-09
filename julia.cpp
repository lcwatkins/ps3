
#include <omp.h>
#include <mpi.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <assert.h>
#include <fstream>

double **alloc2darray(int size);
int **alloc2darrayInts(int size);
void printToFile(int N, int **A);
void printToBinFile(int N, int **A);
void printToBinFile_mpi(int N, int *A, int nprocs, int myrank, int mysize);
void getStartCoord(int myproc, int npts, int ex, int N, int &xinit, int &yinit);
int calcJulia(int maxiter, double cx, double cy, double dx, double dy, int i, int j);

int main(int argc, char *argv[]){
    int N, iter, maxiter;
    double cx, cy, dx, dy, zx, zy, temp;

    if (argc != 2){
        printf("wrong num args\n");
        return 1;
    }
    N = std::stoi(argv[1]);
    maxiter = 1000;
    cx = -0.7;
    cy = 0.26;
    dx = 3.0/N;
    dy = 2.0/N;

    /////////////// Initialize MPI /////////////////
    int nprocs, stat, myrank;
    MPI_Init(&argc, &argv);
    stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert(stat == MPI_SUCCESS);
    stat = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    assert(stat == MPI_SUCCESS);

    //////////// MPI ///////////////
    int i, j, nPtsPerRank, nPtsMyRank, exPts, iters;
    nPtsPerRank = (N*N)/nprocs;
    exPts = (N*N)%nprocs;
    if (myrank == 0){
        printf("per rank: %d, extra: %d\n",nPtsPerRank,exPts);
    }
    nPtsMyRank = nPtsPerRank;
    if (myrank < exPts) {
        nPtsMyRank += 1;
    }
    getStartCoord(myrank, nPtsPerRank, exPts, N, i, j);
    printf("proc %d: (%d,%d)\n",myrank,i,j);
    int* P = new int[nPtsMyRank];
    iters = calcJulia(maxiter, cx, cy, dx, dy, i, j);
    P[0] = iters;
    for (int pt=1;pt<nPtsMyRank;pt++){
        j += 1;
        if (j == N) {
            i += 1;
            j = 0;
        }
        iters = calcJulia(maxiter, cx, cy, dx, dy, i, j);
        P[pt] = iters;
    }
    printf("%d end: (%d,%d)\n",myrank,i,j);
    printToBinFile_mpi(N, P, nprocs, myrank, nPtsMyRank);

//    //// SERIAL ////
//    int **P;
//    P = alloc2darrayInts(N);
//    for (int i=0;i<N;i++){
//        for (int j=0;j<N;j++){
//            zx = -1.5 + dx*i;
//            zy = -1.0 + dy*j;
//            iter = 0;
//            while (((zx*zx + zy*zy) < 4.0) && (iter <= maxiter)) {
//                temp = zx*zx - zy*zy;
//                zy = 2*zx*zy + cy ;
//                zx = temp + cx ;
//                iter += 1;
//            }
//            P[i][j] = iter ;
//        }
//    }
//    printToBinFile(N,P);
//    /////////////////
    
    MPI_Finalize();
    return 0;
}


int calcJulia(int maxiter, double cx, double cy, double dx, double dy, int i, int j){
    double zx, zy, temp;
    int iter;
    zx = -1.5 + dx*i;
    zy = -1.0 + dy*j;
    iter = 0;
    while (((zx*zx + zy*zy) < 4.0) && (iter <= maxiter)) {
        temp = zx*zx - zy*zy;
        zy = 2*zx*zy + cy ;
        zx = temp + cx ;
        iter += 1;
    }
    return iter;
}

// a == index of first point in flattened array, then find coordinate (i,j) of 
// first point on this proc from the entire NxN grid
void getStartCoord(int myproc, int npts, int ex, int N, int &xinit, int &yinit){
    int a;
    a = npts*myproc;
    if (myproc < ex) {
        for (int i=0;i<myproc;i++){
            a += 1;
        }
    } else {
        a += ex;
    }
    xinit = 0; 
    yinit = 0;
    while (a >= N){
        xinit += 1;
        a -= N;
    }
    yinit = a;
}

// allocate a 2d array stored contiguously in memory 
double **alloc2darray(int size) {
    double *data = (double *)malloc(size*size*sizeof(double));
    double **array = (double **)malloc(size*sizeof(double*));
    for (int i=0;i<size;i++){
        array[i] = &(data[size*i]);
        for (int j=0;j<size;j++){
            array[i][j] = 0.0;
        }
    }
    return array;
}

int *alloc1darrayInts(int size){
    int *array = (int *)malloc(size*sizeof(int*));
    for (int i=0;i<size;i++){
        array[i] = 0;
    }
    return array;
}

// allocate a 2d array of ints stored contiguously in memory 
int **alloc2darrayInts(int size) {
    int *data = (int *)malloc(size*size*sizeof(int));
    int **array = (int **)malloc(size*sizeof(int*));
    for (int i=0;i<size;i++){
        array[i] = &(data[size*i]);
        for (int j=0;j<size;j++){
            array[i][j] = 0;
        }
    }
    return array;
}

void printToFile(int N, int **A){
    FILE* file;
    char outname [40];
    int n = sprintf(outname, "out-%d.dat", N);
    file = fopen(outname, "w");
    
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            fprintf(file, "%d ", A[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    return;
}

void printToBinFile(int N, int **A){
    FILE* file;
    char outname [40];
    int n = sprintf(outname, "out-%d.bin", N);
    std::ofstream outf(outname, std::ios::out | std::ios::binary);
    if (!outf){
        printf("cannot open output file\n");
        return;
    }
    outf.write((char *) &N, sizeof(int));
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            outf.write((char *) (&A[i][j]), sizeof(int));
        }
    }
    outf.close();
    return;
}

void printToBinFile_mpi(int N, int *A, int nprocs, int myrank, int mysize){
    // proc 0 is "master" to collect and print all info
    if (myrank == 0){
        MPI_Status status;
        int recvsize;
        FILE* file;
        char outname [40];
        int stat;
        int n = sprintf(outname, "out-%d.bin", N);
        std::ofstream outf(outname, std::ios::out | std::ios::binary);
        //std::ofstream outf(outname, std::ios::out );
        if (!outf){
            printf("cannot open output file\n");
            return;
        }

        outf.write((char *) &N, sizeof(int));
        // first print p0
        for (int i=0;i<mysize;i++){
            outf.write((char *) (&A[i]), sizeof(int));
        }

        int recvbuf[mysize];
        // receive from all other procs
        for (int i=1; i<nprocs; i++){
            stat = MPI_Recv(&recvsize, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            assert(stat == MPI_SUCCESS);
            printf("size %d from %d\n",recvsize,i);
            //int *recvbuf = new int[recvsize];
            stat = MPI_Recv(&recvbuf, recvsize, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            assert(stat == MPI_SUCCESS);
            for (int j=0; j<recvsize; j++){
                outf.write((char *) (&recvbuf[j]), sizeof(int));
            }
        }
        outf.close();
    } else {
        // have each other proc send its array one row at a time
        MPI_Send(&mysize, 1, MPI_INT, 0, 99, MPI_COMM_WORLD);
        MPI_Send(&A[0], mysize, MPI_INT, 0, 99, MPI_COMM_WORLD);
    }
    return;
}
