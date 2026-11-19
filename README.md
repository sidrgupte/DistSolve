
# DistSolve — Parallel Poisson Solvers

Parallel Computing | Stony Brook University  
Project 3, Problem 3.4

## Overview

C++17 and MPI implementation of Jacobi and Gauss–Seidel methods for the 2D Poisson equation using a five-point finite-difference stencil.

The matrix is applied implicitly without storing the full sparse matrix.

## Methods

- Serial Jacobi and Gauss–Seidel
- MPI Jacobi and red–black Gauss–Seidel
- 2D Cartesian domain decomposition
- Halo exchange using `MPI_Sendrecv`
- Global residual calculation using `MPI_Allreduce`

## Requirements

- C++17 compiler
- CMake
- MPI (Open MPI or MPICH)
- Python 3, NumPy, Matplotlib (for plotting)

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

Serial residual table:

```bash
./build/test_serial
```

Parallel correctness tests:

```bash
mpirun -np 4 ./build/test_parallel_jacobi
mpirun -np 4 ./build/test_parallel_rbgs
mpirun -np 4 ./build/test_parallel_residual
```

Performance benchmarks:

```bash
mpirun -np 1 ./build/benchmark
mpirun -np 2 ./build/benchmark
mpirun -np 4 ./build/benchmark
mpirun --use-hwthread-cpus -np 16 ./build/benchmark
```

## Results

The benchmark uses a 128 × 128 interior grid and 1,000 iterations per solver.

| MPI ranks | Jacobi speedup | RBGS speedup |
|---:|---:|---:|
| 1 | 1.00× | 1.00× |
| 2 | 1.46× | 1.71× |
| 4 | 1.60× | 2.39× |
| 16 | 0.95× | 1.50× |

Speedups use median wall-clock times from five runs.

The 16-process benchmark uses hardware-thread allocation on a single machine.

## Report

See `results/report.pdf` for the numerical results, speedup plot, and performance discussion.
