#include <iostream>
#include <ctime>
#include <boost/random.hpp>
#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#define N 100000

using namespace boost::random;
using namespace boost::multiprecision;
using namespace boost::program_options;

cpp_bin_float_50 random_gen(){
	std::random_device seeder;
	boost::random::mt19937 gen(seeder());

	cpp_bin_float_50 rand_num = generate_canonical<cpp_bin_float_50, std::numeric_limits<cpp_bin_float_50>::digits>(gen);
	return rand_num;
}

void seq_monte_carlo(){
	int count;
	int i;
	//int n = vm["n"].as<int>();
	int n = N;

	cpp_bin_float_50 pi;
	cpp_bin_float_50 x,y;
	cpp_bin_float_50 PI = boost::math::constants::pi<cpp_bin_float_50>();

	count = 0;
	for(i = 0; i < n; i++){
		x = random_gen();
		y = random_gen();
		if(x*x + y*y <= 1.0){
			count++;
		}
	}
	pi = 4.0 * (double) count / (double) n;
	
	cpp_bin_float_50 error = (pi - PI) / PI;

	std::cout << std::setprecision(20); 
	std::cout << pi << std::endl;
	std::cout << PI << std::endl;
	std::cout << error << std::endl;


}

void par_monte_carlo(){
	int count_private[N] = {0};
	int count = 0;
	int i;
	//int n = vm["n"].as<int>();
	int n = N;

	cpp_bin_float_50 pi;
	cpp_bin_float_50 x,y;
	cpp_bin_float_50 PI = boost::math::constants::pi<cpp_bin_float_50>();

	count = 0;
	#pragma omp parallel for \
	   private(i,x,y,count_private) \
	   schedule(static,count_private)
	for(i = 0; i < n; i++){
		x = random_gen();
		y = random_gen();
		if(x*x + y*y <= 1.0){
			//count++;
			count_private[i]++;
		}
	}
	#pragma omp parallel for reduction(+:count)
	{
		for(i = 0; i < n; i++){
			count += count_private[i];
		}
	}
	pi = 4.0 * (double) count / (double) n;
	
	cpp_bin_float_50 error = (pi - PI) / PI;

	std::cout << std::setprecision(20); 
	std::cout << pi << std::endl;
	std::cout << PI << std::endl;
	std::cout << error << std::endl;
}


double time_cal(void (*func)()){
	clock_t begin = clock();
	(*func)();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	return elapsed_secs;
}


int main(int argc, char *argv[]){

//	options_description desc{"Options"};
//	desc.add_options()
//		("help,h", "Help screen")
//    		("n", value<int>()->default_value(10000), "Number of random counts")
//	;
//	variables_map vm;
//	store(parse_command_line(argc, argv, desc), vm);
//	notify(vm);
//
//	std::random_device seeder;
//	boost::random::mt19937 gen(seeder());

	auto seq_time = time_cal(seq_monte_carlo);
	auto par_time = time_cal(par_monte_carlo);
	std::cout << seq_time << std::endl;
	std::cout << par_time << std::endl;
}
