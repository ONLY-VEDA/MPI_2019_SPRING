// Author: Wes Kendall
// Copyright 2012 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// Program that computes the average of an array of elements in parallel using
// MPI_Scatter and MPI_Allgather
//
#include <mpi.h>
#include <cassert>
#include <random>
#include <iostream>
#include <boost/mpi/datatype.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

using namespace boost::numeric::ublas;

void decompose_domain(int domain_size, int world_rank, int world_size,
		int *subdomain_start, int *subdomain_size){
	if( world_size > domain_size){
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	*subdomain_start = domain_size / world_size * world_rank;
	*subdomain_size = domain_size / world_size;
	if( world_rank == world_size - 1){
		*subdomain_size = domain_size % world_size;
	}
}

template <typename T>
matrix<T>* gen_matrix(size_t rows = 1024, size_t cols = 1024){
  std::random_device seeder;
  boost::random::mt19937 gen(seeder());
  boost::random::uniform_int_distribution<> dist(0, 255);

  matrix<T>* m = new matrix<T> (rows, cols);
  for (unsigned i = 0; i < m->size1 (); ++ i)
    for (unsigned j = 0; j < m->size2 (); ++ j)
        m->at_element(i, j) = dist(gen);;
  return m;
}

template <typename T>
matrix<T> matrix_mul_row_wise(matrix<T> m1, matrix<T> m2){
  assert(m1.size2() == m2.size1());

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  auto mpi_type = boost::mpi::get_mpi_datatype(m1(0,0));

  // row of first matrix
  int domain_size = m1.size1(); 
  int subdomain_start, subdomain_size;
  decompose_domain(domain_size, world_rank, world_size, &subdomain_start, &subdomain_size);
  // std::cout << subdomain_start << subdomain_size << std::endl;

  int recv_count =  subdomain_size * m1.size2();
  T *sub_vectors = new T[recv_count];
  assert(sub_vectors != NULL);

  if( world_rank == 0){
    // MPI_Scatter( m1.data().begin(), recv_count, mpi_type, sub_vectors, recv_count, mpi_type, 0, MPI_COMM_WORLD);
    MPI_Bcast(m2.data().begin(), m2.size1() * m2.size2(), mpi_type, 0, MPI_COMM_WORLD);
  }

  for(auto i = 0; i < m2.size1() * m2.size2(); i++){
    // std::cout << sub_vectors[i] << std::endl;
    std::cout << m2.data().begin()[i] << std::endl;
    
  }

  delete sub_vectors;

  return m1;
}

int main(int argc, char** argv) {
  // if (argc != 3) {
  //   fprintf(stderr, "Usage: avg num_elements_per_proc num_trials\n");
  //   exit(1);
  // }

  // int num_elements_per_proc = atoi(argv[1]);
  // int num_trials = atoi(argv[2]);


  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  matrix<int>* m1 = NULL;
  matrix<int>* m2 = NULL; 
  if (world_rank == 0) {
    m1 = gen_matrix<int>(5, 5);
    m2 = gen_matrix<int>(5, 5);
    std::cout << *m1 << std::endl;
    std::cout << *m2 << std::endl;
  }

  matrix_mul_row_wise(*m1, *m2);
  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank == 0) {
    delete m1;
    delete m2;
  }
  MPI_Finalize();
}
