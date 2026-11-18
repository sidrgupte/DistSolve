cmake_minimum_required(VERSION 3.20)

project(DistSolve VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED COMPONENTS CXX)

add_executable(mpi_smoke_test tests/mpi_smoke_test.cc)
target_link_libraries(mpi_smoke_test PRIVATE MPI::MPI_CXX)
target_compile_options(mpi_smoke_test PRIVATE -Wall -Wextra -Wpedantic)