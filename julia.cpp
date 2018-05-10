
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
void appendToBinFile(int N, int *A, int size, int firstflag);
void printToBinFile_mpi(int N, int *A, int nprocs, int myrank, int mysize);
void getStartCoord(int myproc, int npts, int ex, int N, int &xinit, int &yinit);
void getCoordFromPt(int pt, int N, int &x, int &y);
int calcJulia(int maxiter, double cx, double cy, double dx, double dy, int i, int j);

int main(int argc, char *argv[]){
    int N, iter, maxiter, chunksize;
    double cx, cy, dx, dy, zx, zy, temp;
    std::string mpimode;

    if (argc < 3){
        printf("wrong num args\n");
        return 1;
    }
    N = std::stoi(argv[1]);
    mpimode = argv[2];
    if (mpimode.compare("dynamic") == 0){
        if (argc != 4){
            printf("wrong num args\n");
            return 1;
        }
        chunksize = std::stoi(argv[3]);
    }


    maxiter = 1000;
    cx = -0.7;
    cy = 0.26;
    dx = 3.0/N;
    dy = 2.0/N;

    ///////////////////////////////////////    INTIALIZE MPI    ///////////////////////////////////////////
    int nprocs, stat, myrank;
    MPI_Init(&argc, &argv);
    stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert(stat == MPI_SUCCESS);
    stat = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    assert(stat == MPI_SUCCESS);


    ///////////////////////////////////////     MPI: STATIC     ///////////////////////////////////////////
    if (mpimode.compare("static") == 0){
        int i, j, nPtsPerRank, nPtsMyRank, exPts, iters;
        nPtsPerRank = (N*N)/nprocs;
        exPts = (N*N)%nprocs;
        //if (myrank == 0){
        //    printf("per rank: %d, extra: %d\n",nPtsPerRank,exPts);
        //}
        nPtsMyRank = nPtsPerRank;
        if (myrank < exPts) {
            nPtsMyRank += 1;
        }
        getStartCoord(myrank, nPtsPerRank, exPts, N, i, j);
        //printf("proc %d: (%d,%d)\n",myrank,i,j);
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
        //printf("%d end: (%d,%d)\n",myrank,i,j);
        MPI_Barrier(MPI_COMM_WORLD);
        printToBinFile_mpi(N, P, nprocs, myrank, nPtsMyRank);
        MPI_Finalize();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////     MPI: DYNAMIC     ///////////////////////////////////////////
    ////  proc 0 is the master: it waits for signal from a proc that it is ready, the message
    ////  it sends is its proc #. Then p0 sends first the next chunksize, then the next starting point 
    ////  (for pt = 0...(N*N-1)) and the slave proc figures out what (i,j) that corresponds to, and 
    ////  calculates iterations for that point through point+chunksize.
    //// 
    ////  Now, each proc just stores all results then sends them one by one at the end, but the var
    ////  "nextSend" keeps track of the next chunk to send, which could be used to send chunks
    ////  whenever, and the command to receive a "send data" signal is in place if the master asked
    ////  for data to print periodically.
    //// 
    ////  MPI TAGS:
    ////    999 = "All finished signal"
    ////    555 = "Send/sending data"
    ////    0/1 = "Chunk done, sending next"
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if (mpimode.compare("dynamic") == 0){
        int thischunk, pt, source, ndone, ntotchunks;
        //chunksize = ((N*N)/(nprocs*2));
        thischunk = chunksize;
        ntotchunks = 0;
        pt = 0;
        MPI_Status status;

        if (myrank == 0) {
            // master
            //printf("%d points over %d procs, chunksize %d\n", (N*N), (nprocs-1), chunksize);
            // resultProcs stores the proc # for each distributed chunk in order
            // allChunkSizes stores the chunksize sent to each proc (all the same except the last one
            // might be smaller -- this framework would easily extend to changing chunk sizes though)
            int *resultProcs = new int[((N*N)/(chunksize)+1)];
            int *allChunkSizes = new int[((N*N)/(chunksize)+1)];
            // print out N to file
            appendToBinFile(N, resultProcs, 0, 1);
            while (pt < (N*N)) {
                if ((N*N)-pt < chunksize) {
                    // and we're done after this one
                    // last chunk is whatever points remain
                    thischunk = (N*N)-pt;
                }
                // wait for message saying 'done/ready'
                // then send more to next proc
                stat = MPI_Recv(&source, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
                assert(stat == MPI_SUCCESS);
                resultProcs[ntotchunks] = source;
                allChunkSizes[ntotchunks] = thischunk;
                MPI_Send(&pt, 1, MPI_INT, source, 0, MPI_COMM_WORLD);
                MPI_Send(&thischunk, 1, MPI_INT, source, 0, MPI_COMM_WORLD);
                pt += thischunk;
                ntotchunks += 1;
            }
            // send done signal
            // nprocs-1 b/c not counting this proc
            while (ndone < (nprocs-1)){
                stat = MPI_Recv(&source, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
                assert(stat == MPI_SUCCESS);
                MPI_Send(&myrank, 1, MPI_INT, source, 999, MPI_COMM_WORLD);
                ndone += 1;
            }
            // get any data stored on procs still
            int *recvbuf = new int[chunksize];
            for (int k=0;k<ntotchunks;k++){
                //printf("receiving chunk %d from proc %d, size %d\n",k,resultProcs[k], allChunkSizes[k]);
                stat = MPI_Recv(&(recvbuf[0]), allChunkSizes[k], MPI_INT, resultProcs[k], 555, MPI_COMM_WORLD, &status);
                assert(stat == MPI_SUCCESS);
                appendToBinFile(N, recvbuf, allChunkSizes[k], 0);
            }
        } else {
            // slave
            // keep A to store pointers to where results are stored -- if divided evenly, this 
            // would be (N*N)/(chunksize*nprocs). Estimate as ... ?? 
            
            int i,j, nchunks, work, mypt, iters, nextSend;
            int estcount = (2*sqrt(N)*N*N)/(chunksize*nprocs);
            int **A = new int*[estcount];
            int *S = new int[estcount];

            nextSend = 0;
            nchunks = 0;
            work = 1;

            while (work){
                MPI_Send(&myrank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                stat = MPI_Recv(&mypt, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                assert(stat == MPI_SUCCESS);
                // check for done signal
                 if (status.MPI_TAG == 999) {
                    work = 0;
                 } else if (status.MPI_TAG == 555) {
                     // if master asks for data to print periodically -- *not* finished implementing!
                     MPI_Send(&(A[nextSend]), S[nextSend], MPI_INT, 0, 555, MPI_COMM_WORLD);
                     delete[] A[nextSend];
                     nextSend += 1;
                 } else {
                    MPI_Recv(&thischunk, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
                    getCoordFromPt(mypt, N, i, j);

                    A[nchunks] = new int[thischunk];
                    S[nchunks] = thischunk;
                    for (int pt=0;pt<thischunk;pt++){
                        if (pt > 0){
                            j += 1;
                            if (j == N) {
                                i += 1;
                                j = 0;
                            }
                        }
                        iters = calcJulia(maxiter, cx, cy, dx, dy, i, j);
                        A[nchunks][pt] = iters;
                    }
                    nchunks += 1;
                }
            }
            // send all remaining data
            for (int k=nextSend;k<nchunks;k++){
                //printf("proc %d sending chunk %d of %d, size %d\n",myrank,k,nchunks,S[k]);
                MPI_Send(&(A[k][0]), S[k], MPI_INT, 0, 555, MPI_COMM_WORLD);
                delete[] A[k];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }


    //////////////////////////////////////       SERIAL       ///////////////////////////////////////////
    if (mpimode.compare("serial") == 0){
        int **P;
        P = alloc2darrayInts(N);
        for (int i=0;i<N;i++){
            for (int j=0;j<N;j++){
                zx = -1.5 + dx*i;
                zy = -1.0 + dy*j;
                iter = 0;
                while (((zx*zx + zy*zy) < 4.0) && (iter <= maxiter)) {
                    temp = zx*zx - zy*zy;
                    zy = 2*zx*zy + cy ;
                    zx = temp + cx ;
                    iter += 1;
                }
                P[i][j] = iter ;
            }
        }
        printToBinFile(N,P);
    }
    
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

// check this
void getCoordFromPt(int pt, int N, int &x, int &y){
    x = 0;
    y = 0;
    while (pt >= N){
        x += 1;
        pt -= N;
    }
    y = pt;
    return;
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

void appendToBinFile(int N, int *A, int size, int firstflag){
    FILE* file;
    char outname [40];
    int n = sprintf(outname, "out-%d.bin", N);
    std::ofstream outf;
    if (firstflag) {
        outf.open (outname, std::ios::out | std::ios::binary);
        outf.write((char *) &N, sizeof(int));
    } else {
        outf.open (outname, std::ios::out | std::ios::app | std::ios::binary);
    }
    if (!outf){
        printf("cannot open output file\n");
        return;
    }
    for (int i=0; i<size; i++){
        outf.write((char *) (&A[i]), sizeof(int));
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
