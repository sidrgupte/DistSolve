#include <iostream>
#include <algorithm>

#include "grid_info.hh"


int local_size(int global_size, int num_procs, int proc_coord){
    /*
    args:
        global_size = number of rows/cols (x/y)
        num_procs   = number of processes
        proc_coord  = processes' coordinates
    returns:
        local grid dimension for the process
    */
    if(proc_coord < (global_size % num_procs)){
        return (global_size / num_procs) + 1;
    }
    return global_size / num_procs;
}

int start_index(int global_size, int num_procs, int proc_coord){
    /*
    args:
        global_size = number of rows/cols (x/y)
        num_procs   = number of processes
        proc_coord  = processes' coordinates
    returns:
        starting global coordinate for the process
    */
    return proc_coord * (global_size / num_procs) + std::min(proc_coord, (global_size % num_procs));
}

GridInfo create_grid(int Nx, int Ny, MPI_Comm world, MPI_Comm& cart_comm){
    /*
    args:
        Nx          = num of global interior rows (x)
        Ny          = num of global interior cols (y)
        world       = input MPI comm
        cart_comm   = output cartesian MPI comm
    returns:
        GridInfo obj containing global grid's dims, process-grid's dims, 
        this rank's coords, local grid dims, starting global indices, nbrng ranks
    */

    int rank, size;
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &size);

    constexpr int ndims{2};
    int dims[ndims]{0, 0};                  // num of processes along rows/cols (x/y)
    MPI_Dims_create(size, 2, dims);

    if((Nx < dims[0]) || (Ny < dims[1])){
        if(rank==0) std::cerr << "Error: process grid exceeds interior grid dims.\n";
        MPI_Abort(world, 1);
    }

    int periods[ndims] = {0, 0};
    int reorder = 0;
    MPI_Cart_create(world, ndims, dims, periods, reorder, &cart_comm);

    int cart_rank;
    MPI_Comm_rank(cart_comm, &cart_rank);

    int coords[2];                                          // rank's process coord. row/col (x/y)
    MPI_Cart_coords(cart_comm, cart_rank, ndims, coords);

    // neighbours:
    int nbr_ranks_x[2];             // nbrs: top, bottom  
    MPI_Cart_shift(cart_comm, 0, 1, &nbr_ranks_x[0], &nbr_ranks_x[1]);

    int nbr_ranks_y[2];             // nbrs: left, right 
    MPI_Cart_shift(cart_comm, 1, 1, &nbr_ranks_y[0], &nbr_ranks_y[1]);

    GridInfo temp_grid;
    temp_grid.Nx = Nx;
    temp_grid.Ny = Ny;

    temp_grid.Rx = dims[0];
    temp_grid.Ry = dims[1];

    temp_grid.rx = coords[0];
    temp_grid.ry = coords[1];

    temp_grid.nx = local_size(Nx, dims[0], coords[0]);
    temp_grid.ny = local_size(Ny,  dims[1], coords[1]);

    temp_grid.sx = start_index(Nx, dims[0], coords[0]);
    temp_grid.sy = start_index(Ny, dims[1], coords[1]);

    temp_grid.up = nbr_ranks_x[0];
    temp_grid.down = nbr_ranks_x[1];
    temp_grid.left = nbr_ranks_y[0];
    temp_grid.right = nbr_ranks_y[1];

    return temp_grid;
}