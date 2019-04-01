cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(send_recv)

set(CMAKE_C_COMPILER "mpicc")

add_executable(send_recive mpi_send_recive.c)
add_executable(ring mpi_ring.c)
add_executable(ping_pong mpi_ping_pong.c)

