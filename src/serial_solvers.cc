#include "serial_solvers.hh"
#include <iostream>
#include <stdexcept>
#include <cmath>

void jacobi_step(
    const std::vector<double>& u,
    const std::vector<double>& b,
    std::vector<double>& next,
    int Nx,
    int Ny
){
    if (Nx <= 0 || Ny <= 0) {
        throw std::invalid_argument("Grid dimensions must be positive");
    }

    const auto N = static_cast<std::size_t>(Nx) * Ny;
    if ((u.size() != N) || (b.size() != N) || (next.size() != N)) {
        throw std::invalid_argument("Invalid vector sizes");
    }

    for(int i{0}; i<Nx; ++i){
        for(int j{0}; j<Ny; ++j){

            int idx = i*Ny + j;

            double up{0.0};
            double bottom{0.0};
            double left{0.0};
            double right{0.0};

            if (i > 0) up = u[(i - 1) * Ny + j];
            if (i < Nx - 1) bottom = u[(i + 1) * Ny + j];
            if (j > 0) left = u[i * Ny + j - 1];
            if (j < Ny - 1) right = u[i * Ny + j + 1];

            next[idx] = (b[idx] + up + bottom + left + right) / 4.0;
        }
    }
}

void gauss_seidel_step(
    std::vector<double>& u,
    const std::vector<double>& b,
    int Nx,
    int Ny
){
    if (Nx <= 0 || Ny <= 0) {
        throw std::invalid_argument("Grid dimensions must be positive");
    }

    const auto N = static_cast<std::size_t>(Nx) * Ny;
    if ((u.size() != N) || (b.size() != N)) {
        throw std::invalid_argument("Invalid vector sizes");
    }

    for(int i{0}; i<Nx; ++i){
        for(int j{0}; j<Ny; ++j){

            int idx = i*Ny + j;

            double up{0.0};
            double bottom{0.0};
            double left{0.0};
            double right{0.0};

            if (i > 0) up = u[(i - 1) * Ny + j];
            if (i < Nx - 1) bottom = u[(i + 1) * Ny + j];
            if (j > 0) left = u[i * Ny + j - 1];
            if (j < Ny - 1) right = u[i * Ny + j + 1];

            u[idx] = (b[idx] + up + bottom + left + right) / 4.0;
        }
    }
}

double residual_norm(
    const std::vector<double>& u,
    const std::vector<double>& b,
    int Nx,
    int Ny
){
    if (Nx <= 0 || Ny <= 0) {
        throw std::invalid_argument("Grid dimensions must be positive");
    }

    const auto N = static_cast<std::size_t>(Nx) * Ny;
    if ((u.size() != N) || (b.size() != N)) {
        throw std::invalid_argument("Invalid vector sizes");
    }

    double sum_sq{0.0};
    for(int i{0}; i<Nx; ++i){
        for(int j{0}; j<Ny; ++j){
            int idx = i*Ny + j;

            double up{0.0};
            double bottom{0.0};
            double left{0.0};
            double right{0.0};

            if (i > 0) up = u[(i - 1) * Ny + j];
            if (i < Nx - 1) bottom = u[(i + 1) * Ny + j];
            if (j > 0) left = u[i * Ny + j - 1];
            if (j < Ny - 1) right = u[i * Ny + j + 1];

            double Au = 4.0 *u[idx] - (up + bottom + left + right);
            double r = Au - b[idx];

            sum_sq += r*r;
        }
    }
    return std::sqrt(sum_sq);   // 2 norm
}