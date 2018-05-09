
CC=mpiCC
CFLAGS=-fopenmp -std=c++11

julia.o: julia.cpp 
	$(CC) $(CFLAGS) julia.cpp -o julia.o

