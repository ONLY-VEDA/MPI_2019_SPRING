#include <mpi.h>
#include <cassert>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <boost/range/iterator_range.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

namespace fs = boost::filesystem;


std::vector<std::string>* scan_files_recursive(const char* root_path){
    fs::path* full_path = new fs::path(root_path);
    std::vector<std::string>* paths = new std::vector<std::string>();

    if(fs::is_directory(*full_path)) {
        // std::cout << *full_path << " is a directory containing:\n";
        for(auto& entry : boost::make_iterator_range(fs::directory_iterator(*full_path), {})){
            paths->push_back(entry.path().string());
        }
    }
    return paths;
}

bool check_path_length_same(const std::vector<std::string>* paths){
    int len = paths->at(0).length();
    for(auto it = paths->begin(); it!=paths->end(); it++){
        if(it->length() != len){
            return false;
        }
    }
    return true;
}


void word_count_big(){

}

void word_count_small(){
    const int root = 0;

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int path_size;
    int num_elements_per_proc;
    std::vector<std::string>* paths = new std::vector<std::string>();
    if(world_rank == root){
        paths = scan_files_recursive("../data/small_file");
        if(!check_path_length_same(paths)){
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        path_size = paths->size();
        num_elements_per_proc = path_size / world_size;
        if(path_size % world_size){
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
        // std::cout << sizeof(paths->at(0).length())  << std::endl;
    }
    MPI_Bcast(&path_size, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&num_elements_per_proc, 1, MPI_INT, root, MPI_COMM_WORLD);
    if(world_rank != root)
        paths = new std::vector<std::string>(path_size);
    // std::cout << world_rank << path_size*sizeof(std::string) << std::endl;
    MPI_Bcast(paths->data(), path_size*sizeof(std::string), MPI_BYTE, root, MPI_COMM_WORLD);
    std::cout << world_rank << " after " << *(paths->data()) << std::endl;

    std::vector<std::string>* sub_paths = new std::vector<std::string>(num_elements_per_proc);
    int scatter_data_length = num_elements_per_proc * sizeof(std::string);

    // std::cout << "rank: " << world_rank << ",path:" << sub_paths->size() << std::endl;

    MPI_Scatter(paths->data(), scatter_data_length , MPI_BYTE, sub_paths->data(), scatter_data_length, MPI_BYTE, root, MPI_COMM_WORLD);
    // for(auto it = paths->data(); it != (paths->data())+paths->size(); it++){
    //     std::cout << *it << std::endl;
    // }
    // std::cout << world_rank << " after " << *(paths->data()+25) << std::endl;
    // for(auto it = sub_paths->begin(); it!=sub_paths->end(); it++){
    //     std::cout << *it << std::endl;
    // }
    // std::cout << num_elements_per_proc << " " << sizeof(fs::path) << " " << scatter_data_length 
    //     << " " << sub_paths->size() <<  std::endl;

    // std::vector<std::string>* res_paths = new std::vector<std::string>(path_size);
    // MPI_Allgather(sub_paths->data(), scatter_data_length, MPI_BYTE, res_paths->data(), scatter_data_length,MPI_BYTE, MPI_COMM_WORLD);
    
    // for(auto it = res_paths->begin(); it!=res_paths->end(); it++){
        // std::cout << *it << std::endl;
    // }
    // std::cout << "rank: " << world_rank << ",path:" << sub_paths->size() << ",out:" << res_paths->size() << std::endl;
    delete paths;
    delete sub_paths;
}

int main(int argc, char** argv){
    MPI_Init(NULL,NULL);

    word_count_small();

    MPI_Finalize();
}