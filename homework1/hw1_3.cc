#include <mpi.h>
#include <cassert>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
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

void split_small(std::vector<std::string>* paths, std::vector<std::string>* sub_paths){
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int num_per_process = paths->size() / world_size;
    int offset = world_rank * num_per_process;
    for(auto it = paths->begin() + offset; it != paths->begin()+offset+num_per_process ;it++){
        sub_paths->push_back(*it);
    }
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

void word_counter_file(std::string path, std::string word, int* count ){
    std::ifstream f(path.c_str());
    std::stringstream buffer;
    buffer << f.rdbuf();
    std::string full_str = buffer.str();

    (*count) = 0;
    std::string::size_type pos = 0;
    while ( (pos = full_str.find(word, pos)) != std::string::npos) {
          (*count)++;
          pos += word.length();
   }
}


void word_count_big(){
    const int root = 0;
    int world_size;
    int world_rank;
    auto comm = MPI_COMM_WORLD;
    const char* path = "../data/big_file/big_100.txt"; 
    const std::string word("cursus");

    MPI_File fh;
    MPI_Status status;
    MPI_Offset size;
    MPI_File_open( comm, path, MPI_MODE_RDWR, MPI_INFO_NULL, &fh );
    MPI_Comm_size( comm, &world_size );
    MPI_Comm_rank( comm, &world_rank );

    MPI_File_get_size(fh, &size);

    long long num_per_process = size / world_size;

    MPI_Offset offset = num_per_process * world_rank;

    char* buffer = new char[num_per_process];
    MPI_File_read_at(fh, offset, buffer, num_per_process, MPI_BYTE, &status ) ;
        
    std::stringstream full_str_buf;
    // buffer << f.rdbuf();
    while(*(buffer) != '\0'){
        full_str_buf << *buffer;
        buffer++;
    }
    auto full_str = full_str_buf.str();

    int count = 0;
    std::string::size_type pos = 0;
    while ( (pos = full_str.find(word, pos)) != std::string::npos) {
            count++;
            pos += word.length();
    }    

    int all_count;
    MPI_Reduce(&count, &all_count, 1, MPI_INT, MPI_SUM, root, comm);

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "rank: " << world_rank << " count: " << count << std::endl; 
    if(world_rank == root){
        std::cout << "big all count:" << all_count << std::endl;
    }
}

void word_count_small(){
    const int root = 0;
    const std::string word("where");

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int path_size;
    int num_elements_per_proc;
    std::vector<std::string>* paths = new std::vector<std::string>();
    paths = scan_files_recursive("../data/small_file");
    if(!check_path_length_same(paths)){
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    path_size = paths->size();
    num_elements_per_proc = path_size / world_size;
    if(path_size % world_size){
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    std::vector<std::string>* sub_paths = new std::vector<std::string>();

    split_small(paths, sub_paths);

    int sub_count_sum = 0;
    int sub_count[num_elements_per_proc];
    int i=0;
    for(auto it = sub_paths->begin(); it!=sub_paths->end(); it++){
        int temp_count;
        word_counter_file(*it, word, &temp_count);
        sub_count[i++] = temp_count;
        sub_count_sum += temp_count;
    }


    // int all_count[path_size];
    // MPI_Allgather(sub_count, num_elements_per_proc, MPI_INT, all_count, num_elements_per_proc, MPI_INT, MPI_COMM_WORLD);

    int all_count;
    MPI_Reduce(&sub_count_sum, &all_count, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "rank: " << world_rank << " sub count sum: " << sub_count_sum << std::endl; 
    if(world_rank == root){
        std::cout << "small all count:" << all_count << std::endl;
    }

    delete paths;
    delete sub_paths;
}

int main(int argc, char** argv){
    MPI_Init(NULL,NULL);

    word_count_small();
    word_count_big();

    MPI_Finalize();
}