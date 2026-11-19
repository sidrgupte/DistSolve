#include <mpi.h>
#include <iostream>

#include "grid_info.hh"
#include "local_field.hh"
#include "halo_exchange.hh"

#include <cassert>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm cart_comm = MPI_COMM_NULL;
    GridInfo grid = create_grid(4, 4, MPI_COMM_WORLD, cart_comm);

    LocalField u(grid.nx, grid.ny);

    // Fill owned cells with the rank number.
    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            u.at(i, j) = static_cast<double>(rank + 1);
        }
    }

    exchange_halos(u, grid, cart_comm);

    auto expected = [](int neighbor) -> double {
        return neighbor == MPI_PROC_NULL
            ? 0.0
            : static_cast<double>(neighbor + 1);
    };

    for (int j = 1; j <= grid.ny; ++j) {
        assert(u.at(0, j) == expected(grid.up));
        assert(u.at(grid.nx + 1, j) == expected(grid.down));
    }

    for (int i = 1; i <= grid.nx; ++i) {
        assert(u.at(i, 0) == expected(grid.left));
        assert(u.at(i, grid.ny + 1) == expected(grid.right));
    }

    std::cout << "Rank " << rank << ": halo test PASSED\n";

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
}