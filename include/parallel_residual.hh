#pragma once

#include <mpi.h>

#include "grid_info.hh"
#include "local_field.hh"

double parallel_residual_norm(
    LocalField& u,
    const LocalField& b,
    const GridInfo& grid,
    MPI_Comm cart_comm
);