#pragma once

#include <mpi.h>

#include "grid_info.hh"
#include "local_field.hh"

void parallel_jacobi_step(
    LocalField& u,
    const LocalField& b,
    LocalField& next,
    const GridInfo& grid,
    MPI_Comm cart_comm
);

void parallel_rbgs_step(
    LocalField& u,
    const LocalField& b,
    const GridInfo& grid,
    MPI_Comm cart_comm
);