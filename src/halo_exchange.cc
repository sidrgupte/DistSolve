#include "halo_exchange.hh"

void exchange_halos(
    LocalField& u,
    const GridInfo& grid,
    MPI_Comm cart_comm
) {
    const int nx = grid.nx;
    const int ny = grid.ny;

    MPI_Sendrecv(
        &u.at(1, 1), ny, MPI_DOUBLE, grid.up, 0,
        &u.at(nx + 1, 1), ny, MPI_DOUBLE, grid.down, 0,
        cart_comm, MPI_STATUS_IGNORE
    );

    MPI_Sendrecv(
        &u.at(nx, 1), ny, MPI_DOUBLE, grid.down, 1,
        &u.at(0, 1), ny, MPI_DOUBLE, grid.up, 1,
        cart_comm, MPI_STATUS_IGNORE
    );

    MPI_Datatype column_type;

    MPI_Type_vector(
        nx,       // num elements
        1,        // 1 per row
        ny + 2,   // dist b/w consecutive rows
        MPI_DOUBLE,
        &column_type
    );

    MPI_Type_commit(&column_type);

    MPI_Sendrecv(
        &u.at(1, 1), 1, column_type, grid.left, 2,
        &u.at(1, ny + 1), 1, column_type, grid.right, 2,
        cart_comm, MPI_STATUS_IGNORE
    );

    MPI_Sendrecv(
        &u.at(1, ny), 1, column_type, grid.right, 3,
        &u.at(1, 0), 1, column_type, grid.left, 3,
        cart_comm, MPI_STATUS_IGNORE
    );

    MPI_Type_free(&column_type);
}