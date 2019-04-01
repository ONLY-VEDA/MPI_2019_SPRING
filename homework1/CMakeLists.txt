cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(homework1)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(all_gather problem1/all_gather.c)
add_executable(matrix problem2/matrix.cc)

