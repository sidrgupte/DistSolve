#pragma once

#include <mpi.h>

#include "grid_info.hh"
#include "local_field.hh"

void exchange_halos(
    LocalField& u,
    const GridInfo& grid,
    MPI_Comm cart_comm
);