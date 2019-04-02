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
#include <vector>
#include <random>
#include <iostream>
#include <boost/mpi/datatype.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

using namespace boost::numeric::ublas;

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
void vm_mul_row_wise(std::vector<T>* sub_vectors, matrix<T>* m, std::vector<T>* sub_out){
  const int m_rows = m->size1();
  const int m_cols = m->size2();
  const int v_cols = m->size1();
  const int v_rows = sub_vectors->size() / v_cols;

  for(int i=0; i < v_rows; i++){
    for(int j=0; j < m_cols; j++){
      T temp = 0;
      for(int k=0; k < m_rows; k++){
        temp += sub_vectors->at(i * v_cols + k) * m->at_element(k, j);
      }
      sub_out->at(i*v_cols + j) = temp;
    }
  }
}

template <typename T>
bool check_equality(matrix<T>* m1, matrix<T>* m2){
  for(int i=0; i < m1->size1(); i++){
    for(int j=0; j < m1->size2(); j++){
      if(m1->at_element(i,j) != m2->at_element(i,j)){
      	return false;
      }
    }
  }
  return true;
}

void matrix_mul(const int rows=1024, const int cols=1024){

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto mpi_type = MPI_INT;
  const int root = 0;

  matrix<int>* m1 = NULL;
  matrix<int>* m2 = NULL; 
  matrix<int> m_origin_out; 
  if (world_rank == root) {
    m1 = gen_matrix<int>(rows, cols);
    m2 = gen_matrix<int>(cols, rows);
    m_origin_out = matrix<int>(prod(*m1, *m2));
    MPI_Bcast(m2->data().begin(), rows * cols, mpi_type, root, MPI_COMM_WORLD);
    MPI_Bcast(m_origin_out.data().begin(), rows * cols, mpi_type, root, MPI_COMM_WORLD);
  }else{
    m2 = new matrix<int>(rows, cols);
    //m_origin_out = new matrix<int>(rows, cols);
    m_origin_out = matrix<int>(rows, cols);
    MPI_Bcast(m2->data().begin(), rows * cols, mpi_type, root, MPI_COMM_WORLD);
    MPI_Bcast(m_origin_out.data().begin(), rows * rows, mpi_type, root, MPI_COMM_WORLD);
  }
  
  if( rows % world_size != 0){
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  int recv_count =  rows / world_size * cols;
  
  std::vector<int>* sub_vectors = new std::vector<int>(recv_count);
  assert(sub_vectors != NULL);

  if(world_rank == root){
    MPI_Scatter( m1->data().begin(), recv_count, mpi_type, sub_vectors->data(), recv_count, mpi_type, root, MPI_COMM_WORLD);
  } else{
    MPI_Scatter( NULL, recv_count, mpi_type, sub_vectors->data(), recv_count, mpi_type, root, MPI_COMM_WORLD);
  }

  // v_rows * cols ,  cols * rows,  v_rows * rows
  std::vector<int>* sub_out = new std::vector<int>(recv_count);
  matrix<int>* m_out = new matrix<int>(rows, rows);
  assert(m_out != NULL);

  vm_mul_row_wise<int>(sub_vectors, m2, sub_out);

  MPI_Allgather(sub_out->data(), recv_count, mpi_type, m_out->data().begin(), recv_count, mpi_type, MPI_COMM_WORLD); 

  auto is_same = check_equality(m_out, &m_origin_out);
  std::cout <<  "Rank:" << world_rank  << ",The two results are same:" << is_same << std::endl; 
  MPI_Barrier(MPI_COMM_WORLD);
  //if(world_rank == 0){
  //  std::cout << "Rank " << world_rank << " Out: " << *m_out <<  std::endl
  //          << m_origin_out << std::endl;
  //}
  delete sub_out;
  delete m_out;
  delete m2;
  if (world_rank == root) {
    delete m1;
  }
}

void matrix_pool(const int rows=1024, const int cols=1024){
}

int main(int argc, char** argv) {

  MPI_Init(NULL, NULL);

  matrix_mul(16, 16);

  MPI_Finalize();
}
