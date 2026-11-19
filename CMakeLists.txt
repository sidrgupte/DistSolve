cmake_minimum_required(VERSION 3.20)

project(DistSolve VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED COMPONENTS CXX)


# TESTS:
add_executable(mpi_smoke_test tests/mpi_smoke_test.cc)
target_link_libraries(mpi_smoke_test PRIVATE MPI::MPI_CXX)
target_compile_options(mpi_smoke_test PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_decomp tests/test_decomposition.cc)
target_link_libraries(test_decomp PRIVATE MPI::MPI_CXX grid_info)
target_include_directories(test_decomp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(test_decomp PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_serial tests/test_serial.cc)
target_link_libraries(test_serial PRIVATE MPI::MPI_CXX serial_solvers)
target_include_directories(test_serial PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(test_serial PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_halo tests/test_halo.cc)
target_link_libraries(test_halo PRIVATE halo_exchange)
target_compile_options(test_halo PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_parallel_jacobi tests/test_parallel_jacobi.cc)
target_link_libraries(test_parallel_jacobi PRIVATE parallel_solvers serial_solvers)
target_compile_options(test_parallel_jacobi PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_parallel_rbgs tests/test_parallel_rbgs.cc)
target_link_libraries(test_parallel_rbgs PRIVATE parallel_solvers)
target_compile_options(test_parallel_rbgs PRIVATE -Wall -Wextra -Wpedantic)

add_executable(benchmark src/benchmark.cc)
target_link_libraries(benchmark PRIVATE parallel_solvers)
target_compile_options(benchmark PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_parallel_residual tests/test_parallel_residual.cc)
target_link_libraries(test_parallel_residual PRIVATE parallel_solvers parallel_residual)
target_compile_options(test_parallel_residual PRIVATE -Wall -Wextra -Wpedantic)


# LIBS:
add_library(grid_info STATIC src/grid_info.cc)
target_link_libraries(grid_info PRIVATE MPI::MPI_CXX)
target_include_directories(grid_info PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(grid_info PRIVATE -Wall -Wextra -Wpedantic)

add_library(serial_solvers STATIC src/serial_solvers.cc)
target_link_libraries(serial_solvers PRIVATE MPI::MPI_CXX)
target_include_directories(serial_solvers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(serial_solvers PRIVATE -Wall -Wextra -Wpedantic)

add_library(halo_exchange STATIC src/halo_exchange.cc)
target_link_libraries(halo_exchange PUBLIC MPI::MPI_CXX grid_info)
target_include_directories(halo_exchange PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(halo_exchange PRIVATE -Wall -Wextra -Wpedantic)

add_library(parallel_solvers STATIC src/parallel_solvers.cc)
target_link_libraries(parallel_solvers PUBLIC halo_exchange)
target_include_directories(parallel_solvers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(parallel_solvers PRIVATE -Wall -Wextra -Wpedantic)

add_library(parallel_residual src/parallel_residual.cc)
target_include_directories(parallel_residual PUBLIC include)
target_link_libraries(parallel_residual PUBLIC halo_exchange)
target_compile_options(parallel_residual PRIVATE -Wall -Wextra -Wpedantic)
