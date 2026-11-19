#include <mpi.h>
#include <iostream>

#include "grid_info.hh"
#include "local_field.hh"
#include "parallel_solvers.hh"
#include "serial_solvers.hh"

#include <cassert>
#include <cmath>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm cart_comm = MPI_COMM_NULL;
    GridInfo grid = create_grid(4, 4, MPI_COMM_WORLD, cart_comm);

    LocalField u(grid.nx, grid.ny);
    LocalField b(grid.nx, grid.ny);
    LocalField next(grid.nx, grid.ny);

    // Initialize the RHS using global column indices.
    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            int global_j = grid.sy + j - 1;
            b.at(i, j) = (global_j < 2) ? 1.0 : 0.0;
        }
    }

    for (int k = 0; k < 10; ++k) {
        parallel_jacobi_step(u, b, next, grid, cart_comm);
        u.data.swap(next.data);
    }

    const double reference[4][4] = {
        {0.63616371, 0.72081661, 0.31264687, 0.11484528},
        {0.85647583, 0.98622894, 0.46491051, 0.17698765},
        {0.85647583, 0.98622894, 0.46491051, 0.17698765},
        {0.63616371, 0.72081661, 0.31264687, 0.11484528}
    };

    std::vector<double> serial_u(16, 0.0);
    std::vector<double> serial_next(16, 0.0);
    std::vector<double> serial_b(16, 0.0);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            serial_b[i * 4 + j] = (j < 2) ? 1.0 : 0.0;
        }
    }

    for (int k = 0; k < 10; ++k) {
        jacobi_step(serial_u, serial_b, serial_next, 4, 4);
        serial_u.swap(serial_next);
    }

    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            const int global_i = grid.sx + i - 1;
            const int global_j = grid.sy + j - 1;

            const double expected = serial_u[global_i * 4 + global_j];

            assert(std::abs(u.at(i, j) - expected) < 1e-12);
        }
    }

    std::cout << "Rank " << rank << ": 10-iteration MPI Jacobi comparison PASSED\n";

    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            std::cout << "Rank " << rank
                    << " global=(" << grid.sx + i - 1
                    << "," << grid.sy + j - 1
                    << ") value=" << u.at(i, j)
                    << '\n';
        }
    }

    std::cout << "Rank " << rank << ": parallel Jacobi test PASSED\n";

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
}