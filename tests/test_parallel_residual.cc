#include <mpi.h>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "grid_info.hh"
#include "local_field.hh"
#include "parallel_solvers.hh"
#include "parallel_residual.hh"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    constexpr int Nx = 128;
    constexpr int Ny = 128;
    constexpr int iterations = 1000;

    MPI_Comm cart_comm = MPI_COMM_NULL;
    GridInfo grid = create_grid(Nx, Ny, MPI_COMM_WORLD, cart_comm);

    int rank;
    MPI_Comm_rank(cart_comm, &rank);

    LocalField u(grid.nx, grid.ny);
    LocalField b(grid.nx, grid.ny);
    LocalField next(grid.nx, grid.ny);

    // Same right-hand side as the benchmark.
    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            const int global_j = grid.sy + j - 1;
            b.at(i, j) = (global_j < Ny / 2) ? 1.0 : 0.0;
        }
    }

    const double initial = parallel_residual_norm(
        u, b, grid, cart_comm
    );

    // Jacobi
    for (int k = 0; k < iterations; ++k) {
        parallel_jacobi_step(u, b, next, grid, cart_comm);
        u.data.swap(next.data);
    }

    const double jacobi_final = parallel_residual_norm(
        u, b, grid, cart_comm
    );

    // Reset solution, including halos.
    std::fill(u.data.begin(), u.data.end(), 0.0);

    // Red-black Gauss-Seidel
    for (int k = 0; k < iterations; ++k) {
        parallel_rbgs_step(u, b, grid, cart_comm);
    }

    const double rbgs_final = parallel_residual_norm(
        u, b, grid, cart_comm
    );

    if (rank == 0) {
        std::cout << std::setprecision(15);
        std::cout << "Processes: " << grid.Rx * grid.Ry << '\n';
        std::cout << "Initial residual: " << initial << '\n';
        std::cout << "Jacobi final residual: " << jacobi_final << '\n';
        std::cout << "RBGS final residual: " << rbgs_final << '\n';
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
}