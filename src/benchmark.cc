#include <mpi.h>
#include <iostream>
#include <algorithm>

#include "grid_info.hh"
#include "local_field.hh"
#include "parallel_solvers.hh"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    constexpr int Nx = 128;
    constexpr int Ny = 128;
    constexpr int iterations = 1000;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm cart_comm = MPI_COMM_NULL;
    GridInfo grid = create_grid(Nx, Ny, MPI_COMM_WORLD, cart_comm);

    LocalField u(grid.nx, grid.ny);
    LocalField b(grid.nx, grid.ny);
    LocalField next(grid.nx, grid.ny);

    // Initialize the right-hand side.
    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            const int global_j = grid.sy + j - 1;
            b.at(i, j) = (global_j < Ny / 2) ? 1.0 : 0.0;
        }
    }

    // Synchronize all ranks before timing.
    MPI_Barrier(cart_comm);
    const double start = MPI_Wtime();

    for (int k = 0; k < iterations; ++k) {
        parallel_jacobi_step(u, b, next, grid, cart_comm);
        u.data.swap(next.data);
    }

    const double elapsed = MPI_Wtime() - start;

    // The slowest rank determines total parallel execution time.
    double max_elapsed = 0.0;

    MPI_Reduce(
        &elapsed,
        &max_elapsed,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        0,
        cart_comm
    );

    if (rank == 0) {
        std::cout << "Jacobi"
                << ", processes=" << grid.Rx * grid.Ry
                << ", grid=" << Nx << "x" << Ny
                << ", iterations=" << iterations
                << ", total_seconds=" << max_elapsed
                << ", seconds_per_iteration="
                << max_elapsed / iterations
                << '\n';
    }
    std::fill(u.data.begin(), u.data.end(), 0.0);

    MPI_Barrier(cart_comm);
    const double rbgs_start = MPI_Wtime();

    for (int k = 0; k < iterations; ++k) {
        parallel_rbgs_step(u, b, grid, cart_comm);
    }

    const double rbgs_elapsed = MPI_Wtime() - rbgs_start;

    double rbgs_max_elapsed = 0.0;

    MPI_Reduce(
        &rbgs_elapsed,
        &rbgs_max_elapsed,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        0,
        cart_comm
    );

    if (rank == 0) {
        std::cout << "RBGS"
                << ", processes=" << grid.Rx * grid.Ry
                << ", grid=" << Nx << "x" << Ny
                << ", iterations=" << iterations
                << ", total_seconds=" << rbgs_max_elapsed
                << ", seconds_per_iteration="
                << rbgs_max_elapsed / iterations
                << '\n';
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
}