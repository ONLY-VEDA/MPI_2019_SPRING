cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(dynamci_receive)

set(CMAKE_C_COMPILER "mpicc")

add_executable(status mpi_status.c)
add_executable(probe mpi_probe.c)

