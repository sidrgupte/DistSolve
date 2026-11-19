#include <iostream>
#include <mpi.h>
#include "grid_info.hh"

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);

    MPI_Comm cart_comm = MPI_COMM_NULL;
    GridInfo grid = create_grid(4, 4, MPI_COMM_WORLD, cart_comm);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::cout << "Rank " << rank
            << ": coords=(" << grid.rx << "," << grid.ry << ")"
            << ", local=(" << grid.nx << "," << grid.ny << ")"
            << ", start=(" << grid.sx << "," << grid.sy << ")\n";

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}