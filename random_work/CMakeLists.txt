cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(random_walk)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(random_walk mpi_random_walk.cc)

