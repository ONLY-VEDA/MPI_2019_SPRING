cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(homework1)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(Boost COMPONENTS program_options system filesystem REQUIRED)
find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

add_executable(hw1_1 hw1_1.c)
add_executable(hw1_2 hw1_2.cc)
add_executable(hw1_3 hw1_3.cc)
add_executable(hw2_1 hw2_1.cc)
add_executable(monte_carlo monte_carlo.cc)

target_link_libraries(hw1_3 
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
)

target_link_libraries(monte_carlo
  ${Boost_LIBRARY}
)

target_compile_options(hw2_1 PRIVATE
  ${OpenMP_CXX_FLAGS}
)  

target_link_libraries(hw2_1 PRIVATE
  ${OpenMP_CXX_FLAGS}
)

