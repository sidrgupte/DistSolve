#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "grid_info.hh"
#include "local_field.hh"
#include "parallel_solvers.hh"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm cart_comm = MPI_COMM_NULL;
    GridInfo grid = create_grid(4, 4, MPI_COMM_WORLD, cart_comm);

    LocalField u(grid.nx, grid.ny);
    LocalField b(grid.nx, grid.ny);

    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            int global_j = grid.sy + j - 1;
            b.at(i, j) = (global_j < 2) ? 1.0 : 0.0;
        }
    }

    for (int k = 0; k < 10; ++k) {
        parallel_rbgs_step(u, b, grid, cart_comm);
    }

    std::vector<double> serial_u(16, 0.0);
    std::vector<double> serial_b(16, 0.0);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            serial_b[i * 4 + j] = (j < 2) ? 1.0 : 0.0;
        }
    }

    // Independent serial red-black Gauss-Seidel reference.
    for (int k = 0; k < 10; ++k) {
        for (int color = 0; color < 2; ++color) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    if ((i + j) % 2 != color)
                        continue;

                    const int idx = i * 4 + j;

                    const double up = (i > 0) ? serial_u[idx - 4] : 0.0;
                    const double down = (i < 3) ? serial_u[idx + 4] : 0.0;
                    const double left = (j > 0) ? serial_u[idx - 1] : 0.0;
                    const double right = (j < 3) ? serial_u[idx + 1] : 0.0;

                    serial_u[idx] = (
                        serial_b[idx] + up + down + left + right
                    ) / 4.0;
                }
            }
        }
    }

    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            const int global_i = grid.sx + i - 1;
            const int global_j = grid.sy + j - 1;

            const double expected = serial_u[global_i * 4 + global_j];

            assert(std::abs(u.at(i, j) - expected) < 1e-12);
        }
    }

    std::cout << "Rank " << rank << ": 10-iteration MPI RBGS comparison PASSED\n";

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
}