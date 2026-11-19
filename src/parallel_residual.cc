#include "parallel_residual.hh"
#include "halo_exchange.hh"

#include <cmath>

double parallel_residual_norm(
    LocalField& u,
    const LocalField& b,
    const GridInfo& grid,
    MPI_Comm cart_comm
) {
    // Ensure neighboring solution values are current.
    exchange_halos(u, grid, cart_comm);

    double local_sum = 0.0;

    for (int i = 1; i <= grid.nx; ++i) {
        for (int j = 1; j <= grid.ny; ++j) {

            const double Au =
                4.0 * u.at(i, j)
                - u.at(i - 1, j)
                - u.at(i + 1, j)
                - u.at(i, j - 1)
                - u.at(i, j + 1);

            const double r = b.at(i, j) - Au;

            local_sum += r * r;
        }
    }

    double global_sum = 0.0;

    MPI_Allreduce(
        &local_sum,
        &global_sum,
        1,
        MPI_DOUBLE,
        MPI_SUM,
        cart_comm
    );

    return std::sqrt(global_sum);
}