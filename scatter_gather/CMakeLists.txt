cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(scatter_gather)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(all_avg all_avg.c)
add_executable(all_avg_compare all_avg_compare.c)
add_executable(avg avg.c)
add_executable(avg_compare avg_compare.c)

