// Author: Wes Kendall
// Copyright 2012 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// Program that computes the average of an array of elements in parallel using
// MPI_Scatter and MPI_Allgather
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>

void NAIVE_Bcast(void* data, int count, MPI_Datatype datatype, int root,
              MPI_Comm communicator) {
  int world_rank;
  MPI_Comm_rank(communicator, &world_rank);
  int world_size;
  MPI_Comm_size(communicator, &world_size);

  if (world_rank == root) {
    int i;
    for (i = 0; i < world_size; i++) {
      if (i != world_rank) {
        MPI_Send(data, count, datatype, i, 0, communicator);
      }
    }
  } else {
    MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
  }
}

void NAIVE_Gather(void* send_data,
                int send_count,
                MPI_Datatype send_datatype,
                void* recv_data,
                int recv_count,
                MPI_Datatype recv_datatype,
                int root,
                MPI_Comm communicator) {
  int world_rank;
  MPI_Comm_rank(communicator, &world_rank);
  int world_size;
  MPI_Comm_size(communicator, &world_size);

  if (root == world_rank) {
      int i;
      int type_size;
      MPI_Type_size(recv_datatype, &type_size);
      for (i = 0; i < world_size; i++) {
          if( i == world_rank){
            memcpy( (recv_data + recv_count * type_size * i), send_data, send_count * type_size );
          } else {
            MPI_Recv( (recv_data + recv_count * type_size * i), recv_count, recv_datatype, i, 0, communicator, MPI_STATUS_IGNORE);
          }
      }
  } else {
      MPI_Send(send_data , send_count, send_datatype, root, 0, communicator);
  }
}

void NAIVE_Allgather(void* send_data,
    int send_count,
    MPI_Datatype send_datatype,
    void* recv_data,
    int recv_count,
    MPI_Datatype recv_datatype,
    MPI_Comm communicator) {

  int world_rank;
  MPI_Comm_rank(communicator, &world_rank);
  int world_size;
  MPI_Comm_size(communicator, &world_size);
  
  // root rank 0
  NAIVE_Gather(send_data, send_count, send_datatype, recv_data, recv_count, recv_datatype, 0, MPI_COMM_WORLD);
  NAIVE_Bcast(recv_data, recv_count * world_size, recv_datatype, 0, MPI_COMM_WORLD);
}

// Creates an array of random numbers. Each number has a value from 0 - 1
float *create_rand_nums(int num_elements) {
  float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (float)RAND_MAX);
  }
  return rand_nums;
}

// Computes the average of an array of numbers
float compute_avg(float *array, int num_elements) {
  float sum = 0.f;
  int i;
  for (i = 0; i < num_elements; i++) {
    sum += array[i];
  }
  return sum / num_elements;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: avg num_elements_per_proc num_trials\n");
    exit(1);
  }

  int num_elements_per_proc = atoi(argv[1]);
  int num_trials = atoi(argv[2]);
  // Seed the random number generator to get different results each time
  srand(time(NULL));

  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Create a random array of elements on the root process. Its total
  // size will be the number of elements per process times the number
  // of processes
  float *rand_nums = NULL;
  if (world_rank == 0) {
    rand_nums = create_rand_nums(num_elements_per_proc * world_size);
  }

  // For each process, create a buffer that will hold a subset of the entire
  // array
  float *sub_rand_nums = (float *)malloc(sizeof(float) * num_elements_per_proc);
  assert(sub_rand_nums != NULL);

  // Scatter the random numbers from the root process to all processes in
  // the MPI world
  MPI_Scatter(rand_nums, num_elements_per_proc, MPI_FLOAT, sub_rand_nums,
              num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  // Compute the average of your subset
  float sub_avg = compute_avg(sub_rand_nums, num_elements_per_proc);

  // Gather all partial averages down to all the processes
  float *sub_avgs = (float *)malloc(sizeof(float) * world_size);
  assert(sub_avgs != NULL);

  //double total_mpi_all_gather_time;

  //MPI_Barrier(MPI_COMM_WORLD);
  //total_mpi_all_gather_time -= MPI_Wtime();
  //for (int i = 0; i < num_trials; i++) {
  //  // Time my_bcast
  //  // Synchronize before starting timing
  //  //MPI_Barrier(MPI_COMM_WORLD);
  //  //total_my_bcast_time -= MPI_Wtime();
  //  //my_bcast(data, num_elements, MPI_INT, 0, MPI_COMM_WORLD);
  //  //// Synchronize again before obtaining final time
  //  //MPI_Barrier(MPI_COMM_WORLD);
  //  //total_my_bcast_time += MPI_Wtime();
  //
  //  // Time MPI_Bcast
  //  //MPI_Barrier(MPI_COMM_WORLD);
  //  //total_mpi_all_gather_time -= MPI_Wtime();
  //  MPI_Allgather(&sub_avg, 1, MPI_FLOAT, sub_avgs, 1, MPI_FLOAT, MPI_COMM_WORLD);
  //  //MPI_Barrier(MPI_COMM_WORLD);
  //  //total_mpi_all_gather_time += MPI_Wtime();
  //  //printf("Total mpi all gather time : %f at %d\n", total_mpi_all_gather_time, i);
  //}
  //MPI_Barrier(MPI_COMM_WORLD);
  //total_mpi_all_gather_time += MPI_Wtime();
  //printf("Total mpi all gather time : %f at process %d\n", total_mpi_all_gather_time, world_rank);
  //MPI_Allgather(&sub_avg, 1, MPI_FLOAT, sub_avgs, 1, MPI_FLOAT, MPI_COMM_WORLD);
  NAIVE_Allgather(&sub_avg, 1, MPI_FLOAT, sub_avgs, 1, MPI_FLOAT, MPI_COMM_WORLD);

  // Now that we have all of the partial averages, compute the
  // total average of all numbers. Since we are assuming each process computed
  // an average across an equal amount of elements, this computation will
  // produce the correct answer.
  float avg = compute_avg(sub_avgs, world_size);
  printf("Avg of all elements from proc %d is %f\n", world_rank, avg);

  // Clean up
  if (world_rank == 0) {
    free(rand_nums);
  }
  free(sub_avgs);
  free(sub_rand_nums);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
