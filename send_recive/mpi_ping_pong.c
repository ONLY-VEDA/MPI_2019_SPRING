#define PING_PONG_LIMMIT 3

#include<mpi.h>
#include<stdio.h>
#include<time.h>

int main(){
	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int ping_pong_count = 0;
	int parterner_rank = (world_rank + 1) % 2;

	struct timespec spec;
	long ms;
	time_t s;
	while(ping_pong_count < PING_PONG_LIMMIT){
		if(world_rank == 0){
			ping_pong_count++;
		
			clock_gettime(CLOCK_REALTIME, &spec);
			printf("Before process %d is sending number %d to %d, time: %ld \n", world_rank, ping_pong_count, parterner_rank,  spec.tv_nsec);
			MPI_Send(&ping_pong_count, 1, MPI_INT, parterner_rank, 0, MPI_COMM_WORLD);

			clock_gettime(CLOCK_REALTIME, &spec);
			printf("After process %d is sending number %d to %d, time: %ld \n", world_rank, ping_pong_count, parterner_rank,  spec.tv_nsec);
		
		}else if(world_rank == 1){
			clock_gettime(CLOCK_REALTIME, &spec);
			printf("Before process %d is recive number %d to %d, time: %ld \n", world_rank, ping_pong_count, parterner_rank,  spec.tv_nsec);
			
			MPI_Recv(&ping_pong_count, 1, MPI_INT, parterner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			clock_gettime(CLOCK_REALTIME, &spec);
			printf("After process %d is recive number %d to %d, time: %ld \n", world_rank, ping_pong_count, parterner_rank,  spec.tv_nsec);
		
		}
	}
	MPI_Finalize();
}
