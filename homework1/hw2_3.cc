#include <omp.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <functional>
#include <cstring>
#include <sys/time.h>
#include <boost/random.hpp>

#define N_NODES 1024000

using namespace std;
using namespace boost::random;

int random_gen(int low_lim, int high_lim);

typedef struct Edge{
	int src, dest;
}Edge;

class Graph{
	public:
	vector<vector<int>>* adjList;
	Graph(){
		this->adjList = new vector<vector<int>>(N_NODES, vector<int>(10));
	};
	Graph(const Graph &g){
		this->adjList = new vector<vector<int>>(N_NODES);
		*adjList = *(g.adjList);
	};
	Graph(vector<vector<int>> *adjList){
		this->adjList =  new vector<vector<int>>();
		this->adjList->resize(adjList->size());
		for(int i=0; i< adjList->size(); i++){
			for(auto v:adjList->at(i)){
				(this->adjList)->at(i).push_back(v);
			}
		}
		*(this->adjList) = *adjList;
	};
	Graph(vector<Edge> const &edges, int n_nodes){
		this->adjList =  new vector<vector<int>>();
		adjList->resize(n_nodes);
		for(auto &edge:edges){
			(*adjList)[edge.src].push_back(edge.dest);
		}
	};
};

class PageRank{
	public:
	vector<float>* PR;
	Graph* g;
	PageRank(Graph *g){
		this->PR = new vector<float>(N_NODES, 1);
		this->g = new Graph(*g);
	};
	void Update(int num_iter){
		vector<float>* PR_new = new vector<float>(N_NODES, 0.0);
		while(num_iter--){
			PR_new->resize(N_NODES, 0.0);
			for(int i = 0; i < this->g->adjList->size(); i++){
				int adjLen = this->g->adjList->at(i).size();
				for(int j = 0; j < adjLen; j++){
					int dest = this->g->adjList->at(i).at(j);
					PR_new->at(dest) += this->PR->at(i) / adjLen;
				}
			}
			this->PR = PR_new;
		}
	};
	void UpdatePar(int num_iter){
		vector<float>* PR_new = new vector<float>(N_NODES, 0.0);
		omp_lock_t lock[N_NODES];
		for (int i=0; i<N_NODES; i++){
			omp_init_lock(&(lock[i]));
		}
		while(num_iter--){
			PR_new->resize(N_NODES, 0.0);
			int i,j;
			#pragma omp parallel for shared(PR_new) private(i,j)
			for(i = 0; i < this->g->adjList->size(); i++){
				int adjLen = this->g->adjList->at(i).size();
				for(j = 0; j < adjLen; j++){
					int dest = this->g->adjList->at(i).at(j);
					//#pragma omp critical
					omp_set_lock(&(lock[dest]));
					PR_new->at(dest) += this->PR->at(i) / adjLen;

					omp_unset_lock(&(lock[dest]));
				}
			}
			this->PR = PR_new;
		}
		for (int i=0; i<N_NODES; i++){
			omp_destroy_lock(&(lock[i]));
		}
	};
	void Print(){
		for(auto pr:*PR){
			std::cout << pr << " ";
		}
		std::cout <<  std::endl;
	};
	bool operator ==(const PageRank& rt){
		if(this->PR->size() != rt.PR->size()){
			return false;
		}
		for(int i=0; i < this->PR->size(); i++){
			//std::cout <<"i: "<< i << " " << this->PR->at(i) << " " << rt.PR->at(i) << std::endl;
			if( fabs(this->PR->at(i) - rt.PR->at(i)) > 0.001 ){
				return false;
			}
		};
		return true;
	}
};

int random_gen(int low_lim, int high_lim){
	std::random_device seeder;
	boost::random::mt19937 gen(seeder());
	boost::random::uniform_int_distribution<> dist(low_lim,high_lim);
	
	int rand_num = dist(gen);
	return rand_num;
}

Graph *g;

void init(){
	vector<vector<int>>* adjList = new vector<vector<int>>(N_NODES);
	#pragma omp parallel for shared(adjList)
	for(int i = 0; i< N_NODES; i++){
		int n_edge_this_node = random_gen(1,10);
		vector<int> edges;
		while(n_edge_this_node > 0){
			bool edge_exists = false;
			int dest = random_gen(0, N_NODES - 1);
			if(dest == i){
				continue;
			}
			for(auto it = edges.begin(); it != edges.end(); ++it){
				if(*it == dest){
					edge_exists = true;
					break;
				}
			}
			if(edge_exists){
				continue;
			}
			n_edge_this_node--;
			edges.push_back(dest);
		}
		//adjList->at(i) = *(new vector<int>(edges));
		adjList->at(i) = edges;
	}
	std::cout << "Finished Init Graph" << std::endl;
	g = new Graph(adjList);
}

std::string time_cal(std::function<void()> func){
	
	struct timespec start,end;

	clock_gettime(CLOCK_REALTIME, &start);
	func();
	clock_gettime(CLOCK_REALTIME, &end);
	std::stringstream ss;
	ss <<  end.tv_sec - start.tv_sec; 
	ss <<  " ";
	ss <<  end.tv_nsec - start.tv_nsec; 
	
	return ss.str();
}

void print_graph(Graph* g){
	for(int i=0; i < g->adjList->size(); i++){
		std:cout << "Edges: src at " << i << "Len dests:" << g->adjList->at(i).size() << std::endl;
		for(auto v:g->adjList->at(i)){
			std::cout << v <<" ";
		}
		std::cout << std::endl;
	}
}

int main(int argc, char *argv[]){

	init();
	PageRank* pr = new PageRank(g);
	PageRank* pr_par = new PageRank(g);
	bool is_same;
	int num_iter = 100;

	auto up_seq = std::bind(&PageRank::Update, pr, num_iter);
	auto up_par = std::bind(&PageRank::UpdatePar, pr_par, num_iter);
	auto seq_time = time_cal(up_seq);
	auto par_time = time_cal(up_par);

	is_same = (*pr) == (*pr_par);
	std::cout << is_same << std::endl;
	std::cout << seq_time << std::endl;
	std::cout << par_time << std::endl;
}
