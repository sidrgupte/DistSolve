#pragma once

#include <mpi.h>

struct GridInfo{
    // interior grid:
    int Nx;
    int Ny;

    // MPI process grid:
    int Rx;
    int Ry;

    // rank coord:
    int rx;
    int ry;

    // local owned grid dims:
    int nx;
    int ny;

    // starting coord of rank:
    int sx;
    int sy;

    // neighbour ranks:
    int up;
    int down;
    int left;
    int right;
};

int local_size(int global_size, int num_procs, int proc_coord);
int start_index(int global_size, int num_procs, int proc_coord);
GridInfo create_grid(int Nx, int Ny, MPI_Comm world, MPI_Comm& cart_comm);
