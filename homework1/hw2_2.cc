#include <omp.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include <boost/random.hpp>
#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#define N 10000000

using namespace boost::random;
using namespace boost::program_options;


int* seq_nums = new int[N];
int* par_nums = new int[N];

int random_gen(){
	std::random_device seeder;
	boost::random::mt19937 gen(seeder());
	boost::random::uniform_int_distribution<> dist(1,10000);
	
	int rand_num = dist(gen);
	return rand_num;
}

void init(){
	#pragma omp parallel for shared(seq_nums)
	for(int i = 0; i< N; i++){
		seq_nums[i] = random_gen();
	}
	std::cout << std::endl << "Finished gen data" << std::endl;
	memcpy(par_nums, seq_nums, sizeof(int)*N);
}

void swap(int *a, int *b){
	int temp = *a;
	*a = *b;
	*b = temp;
}

int partition(int *nums, int low, int high){
	int i = low;
	int j = high;
	int pivot = nums[low];
	while( i < j){
		if(nums[j] < pivot){
			swap(&nums[i], &nums[j]);
			i++;
			
		}
		else if(nums[i] > pivot){
			swap(&nums[i], &nums[j]);
			j--;
		}
		else{
			i++;
		}
	}
	nums[i] = pivot;
	return i;
		
}


void seq_qsort(int *nums, int low, int high){
	if(low < high){
		int pivot_idx = partition(nums, low, high);
		
		seq_qsort(nums, low, pivot_idx-1 );
		seq_qsort(nums, pivot_idx+1, high );
	}
}


void seq_qsort_wrap(){
	int count;
	int i;
	int n = N;
	
	seq_qsort(seq_nums, 0, n-1);

	std::cout << std::endl;
}


void par_qsort(int *nums, int low, int high){
	if(low < high){
		int pivot_idx;
		pivot_idx = partition(nums, low, high);

		#pragma omp task 
		par_qsort(nums, low, pivot_idx-1);
		#pragma omp task 
		par_qsort(nums, pivot_idx+1, high);
		
	}
}

void par_qsort_wrap(){
	int count;
	int i;
	int n = N;
	#pragma omp parallel default(none) shared(par_nums, n)
	{
	#pragma omp single nowait
	{par_qsort(par_nums, 0, n-1);}
	}

}


std::string time_cal(void (*func)()){
	
	struct timespec start,end;

	clock_gettime(CLOCK_REALTIME, &start);
	(*func)();
	clock_gettime(CLOCK_REALTIME, &end);
	std::stringstream ss;
	ss <<  end.tv_sec - start.tv_sec; 
	ss <<  " ";
	ss <<  end.tv_nsec - start.tv_nsec; 
	
	return ss.str();
}

void print_nums(){
	for(int i = 0; i< N; i++){
		std::cout << seq_nums[i] << " ";
	}
	std::cout << std::endl;


	for(int i = 0; i< N; i++){
		std::cout << par_nums[i] << " ";
	}
	std::cout << std::endl;
}

bool check_same(){
	for(int i=0;i <N; i++){
		if(seq_nums[i] != par_nums[i]){
			return false;
		}
	}
	return true;
}


int main(int argc, char *argv[]){

	init();
	auto seq_time = time_cal(seq_qsort_wrap);
	auto par_time = time_cal(par_qsort_wrap);
	std::cout << "Same:" << check_same() << std::endl;

	std::cout << seq_time << std::endl;
	std::cout << par_time << std::endl;
}
