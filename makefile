
CC=mpiCC
CFLAGS=-fopenmp -std=c++11

bandwidth.o: bandwidth.cpp 
	$(CC) $(CFLAGS) bandwidth.cpp -o bandwidth.o

