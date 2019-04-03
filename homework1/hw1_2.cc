#define MIN(a,b) a < b ? a : b  
#define MAX(a,b) a > b ? a : b  

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
matrix<T>* gen_matrix(size_t rows = 1024, size_t cols = 1024,int min = 0, int max = 3){
  std::random_device seeder;
  boost::random::mt19937 gen(seeder());
  boost::random::uniform_int_distribution<> dist(min, max);

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
  }else{
    m2 = new matrix<int>(rows, cols);
    m_origin_out = matrix<int>(rows, cols);
  }
  MPI_Bcast(m2->data().begin(), rows * cols, mpi_type, root, MPI_COMM_WORLD);
  MPI_Bcast(m_origin_out.data().begin(), rows * rows, mpi_type, root, MPI_COMM_WORLD);
  
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
  // MPI_Barrier(MPI_COMM_WORLD);

  delete sub_out;
  delete m_out;
  delete m2;
  if (world_rank == root) {
    delete m1;
  }
}

template <typename T>
int get_max(vector<T>* v){
  T max = v->at(0);
  for(int i=0; i<v->size(); i++){
    max = (max < v->at(i)) ? v->at(i) : max;
  }
  return max;
}

// (W−F+2P)/S+1
// here just implement max pooling
void matrix_pool(const int rows=1024, const int cols=1024, const int k_h = 4, const int k_w = 4,
		const int p_h = 0, const int p_w = 0, const int stride_h = 4, const int stride_w = 4){
 
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto mpi_type = MPI_INT;
  const int root = 0;

  matrix<int>* M = NULL;
  if( world_rank == root){
    M = gen_matrix<int>(rows, cols);
  }else{
    M = new matrix<int>(rows, cols);
  }
  MPI_Bcast(M->data().begin(), rows*cols, mpi_type, root, MPI_COMM_WORLD);

  if( rows % k_h != 0 || cols % k_w != 0){
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  const int new_rows = (rows - k_h + 2 * p_h) / stride_h + 1;
  const int new_cols = (cols - k_w + 2 * p_w) / stride_w + 1;

  matrix<int>* N = new matrix<int>(new_rows, new_cols);
  int num_elements_per_proc = new_rows * new_cols / world_size;
  std::vector<int>* sub_vec = new std::vector<int>(num_elements_per_proc);
  MPI_Scatter(N->data().begin(), num_elements_per_proc, mpi_type, sub_vec->data(), num_elements_per_proc, mpi_type, root, MPI_COMM_WORLD);
  
  for(int i=0; i < num_elements_per_proc; i++){
    int pool_index = i*world_size + world_rank;

    int nr = pool_index / new_cols;
    int nc = pool_index % new_cols;
    
    int hstart = nr * stride_h - p_h;
    int wstart = nc * stride_w - p_w;
    int hend = MIN(hstart + k_h, rows);
    int wend = MIN(wstart + k_w, cols);
    hstart = MAX(hstart, 0);
    wstart = MAX(wstart, 0); 

    int temp_max = INT_MIN;
    for(int h=hstart; h < hend; h++){
      for(int w=wstart; w < wend; w++){
        int index = h * cols + w;
        int bottom_data = M->at_element(h, w);
        if(bottom_data > temp_max ){
          temp_max = bottom_data;
        }
      }
    }
    sub_vec->at(i) = temp_max;
  }

  MPI_Allgather(sub_vec->data(), num_elements_per_proc, mpi_type, N->data().begin(), num_elements_per_proc,mpi_type, MPI_COMM_WORLD);
  std::cout<< *M << std::endl;
  std::cout<< *N << std::endl;
}

// (W−F+2P)/S+1
// just implement naive conv
// TODO: im2col  
void matrix_conv(const int rows=1024, const int cols=1024, const int k_h = 4, const int k_w = 4,
		const int p_h = 0, const int p_w = 0, const int stride_h = 4, const int stride_w = 4){
 
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto mpi_type = MPI_INT;
  const int root = 0;

  matrix<int>* M = NULL;
  matrix<int>* K = NULL;
  if( world_rank == root){
    M = gen_matrix<int>(rows, cols);
    K = gen_matrix<int>(k_h, k_w);
  }else{
    M = new matrix<int>(rows, cols);
    K = new matrix<int>(k_h, k_w);
  }
  MPI_Bcast(M->data().begin(), rows*cols, mpi_type, root, MPI_COMM_WORLD);
  MPI_Bcast(K->data().begin(), k_h * k_w, mpi_type, root, MPI_COMM_WORLD);

  if( rows % k_h != 0 || cols % k_w != 0){
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  const int new_rows = (rows - k_h + 2 * p_h) / stride_h + 1;
  const int new_cols = (cols - k_w + 2 * p_w) / stride_w + 1;

  matrix<int>* N = new matrix<int>(new_rows, new_cols);
  int num_elements_per_proc = new_rows * new_cols / world_size;
  std::vector<int>* sub_vec = new std::vector<int>(num_elements_per_proc);
  MPI_Scatter(N->data().begin(), num_elements_per_proc, mpi_type, sub_vec->data(), num_elements_per_proc, mpi_type, root, MPI_COMM_WORLD);
  
  for(int i=0; i < num_elements_per_proc; i++){
    int pool_index = i*world_size + world_rank;

    int nr = pool_index / new_cols;
    int nc = pool_index % new_cols;
    
    int hstart = nr * stride_h - p_h;
    int wstart = nc * stride_w - p_w;
    int hend = MIN(hstart + k_h, rows);
    int wend = MIN(wstart + k_w, cols);
    hstart = MAX(hstart, 0);
    wstart = MAX(wstart, 0); 

    int temp_prod = 0;
    for(int h=hstart; h < hend; h++){
      for(int w=wstart; w < wend; w++){
        int index = h * cols + w;
        int bottom_data = M->at_element(h, w);
        int weight = K->at_element(h-hstart, w-wstart);
        temp_prod += bottom_data * weight;
      }
    }
    sub_vec->at(i) = temp_prod;
  }

  MPI_Allgather(sub_vec->data(), num_elements_per_proc, mpi_type, N->data().begin(), num_elements_per_proc,mpi_type, MPI_COMM_WORLD);
  std::cout<< *M << std::endl;
  std::cout<< *K << std::endl;
  std::cout<< *N << std::endl;

}

int main(int argc, char** argv) {

  MPI_Init(NULL, NULL);

  matrix_mul(1024, 1024);
  matrix_pool(1024,1024);
  matrix_conv(1024, 1024);

  MPI_Finalize();
}
