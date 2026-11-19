#pragma once

#include <vector>

void jacobi_step(
    const std::vector<double>& u,
    const std::vector<double>& b,
    std::vector<double>& next,
    int Nx,
    int Ny
);

void gauss_seidel_step(
    std::vector<double>& u,
    const std::vector<double>& b,
    int Nx,
    int Ny
);

double residual_norm(
    const std::vector<double>& u,
    const std::vector<double>& b,
    int Nx,
    int Ny
);