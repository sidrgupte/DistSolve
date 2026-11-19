#include <iostream>
#include <iomanip>

#include <vector>
#include "serial_solvers.hh"

int main(){
    constexpr int Nx{4};
    constexpr int Ny{4};

    std::vector<double> b{
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 0, 0
    };

    std::vector<double> u(Nx * Ny, 0.0);
    std::vector<double> next(Nx * Ny, 0.0);
    std::vector<double> gs(Nx * Ny, 0.0);

    std::cout << std::setprecision(10);
    std::cout << "Iteration\tJacobi\t\tGauss-Seidel\n";
    std::cout << 0 << '\t' << residual_norm(u, b, Nx, Ny) << '\t' << residual_norm(gs, b, Nx, Ny) << '\n';

    for (int k = 1; k <= 10; ++k){
        jacobi_step(u, b, next, Nx, Ny);
        u.swap(next);

        gauss_seidel_step(gs, b, Nx, Ny);

        std::cout << k << '\t' << residual_norm(u, b, Nx, Ny) << '\t' << residual_norm(gs, b, Nx, Ny) << '\n';
    }

    return 0;
}