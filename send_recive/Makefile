EXECS=mpi_send_recive mpi_ping_pong mpi_ring
MPICC?=mpicc

all: ${EXECS}

mpi_send_recive: mpi_send_recive.c
	${MPICC} -o mpi_send_recive mpi_send_recive.c

mpi_ring: mpi_ring.c
	${MPICC} -o mpi_ring mpi_ring.c

mpi_ping_pong: mpi_ping_pong.c
	${MPICC} -o mpi_ping_pong mpi_ping_pong.c

clean:
	rm ${EXECS}
