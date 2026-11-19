#include "parallel_solvers.hh"
#include "halo_exchange.hh"

void parallel_jacobi_step(
    LocalField& u,
    const LocalField& b,
    LocalField& next,
    const GridInfo& grid,
    MPI_Comm cart_comm
) {
    exchange_halos(u, grid, cart_comm);

    const int nx = grid.nx;
    const int ny = grid.ny;

    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            next.at(i, j) = (
                b.at(i, j)
                + u.at(i - 1, j)
                + u.at(i + 1, j)
                + u.at(i, j - 1)
                + u.at(i, j + 1)
            ) / 4.0;
        }
    }
}


void parallel_rbgs_step(
    LocalField& u,
    const LocalField& b,
    const GridInfo& grid,
    MPI_Comm cart_comm
) {
    // Phase 1: Update red cells.
    exchange_halos(u, grid, cart_comm);

    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            const int global_i = grid.sx + i - 1;
            const int global_j = grid.sy + j - 1;

            if ((global_i + global_j) % 2 == 0) {
                u.at(i, j) = (
                    b.at(i, j)
                    + u.at(i - 1, j)
                    + u.at(i + 1, j)
                    + u.at(i, j - 1)
                    + u.at(i, j + 1)
                ) / 4.0;
            }
        }
    }

    // Phase 2: Exchange updated red values.
    exchange_halos(u, grid, cart_comm);

    // Phase 3: Update black cells.
    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {
            const int global_i = grid.sx + i - 1;
            const int global_j = grid.sy + j - 1;

            if ((global_i + global_j) % 2 == 1) {
                u.at(i, j) = (
                    b.at(i, j)
                    + u.at(i - 1, j)
                    + u.at(i + 1, j)
                    + u.at(i, j - 1)
                    + u.at(i, j + 1)
                ) / 4.0;
            }
        }
    }
}
