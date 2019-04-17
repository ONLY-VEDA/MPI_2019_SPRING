cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(homework2)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(IS_DEBUG TRUE)

if(IS_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
endif()


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

find_package(Boost COMPONENTS program_options system filesystem random REQUIRED)
find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${OpenMP_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${SYSTEM_INCLUDE_DIRS})

add_executable(hw2_1 hw2_1.cc)
add_executable(hw2_2 hw2_2.cc)
add_executable(hw2_3 hw2_3.cc)

target_compile_options(hw2_1 PRIVATE
  ${OpenMP_CXX_FLAGS}
)  

target_link_libraries(hw2_1 PRIVATE
  ${Boost_LIBRARIES}
  ${OpenMP_CXX_FLAGS}
)

target_compile_options(hw2_2 PRIVATE
  ${OpenMP_CXX_FLAGS}
)  

target_link_libraries(hw2_2 PRIVATE
  ${Boost_LIBRARIES}
  ${OpenMP_CXX_FLAGS}
)

target_compile_options(hw2_3 PRIVATE
  ${OpenMP_CXX_FLAGS}
)  

target_link_libraries(hw2_3 PRIVATE
  ${Boost_LIBRARIES}
  ${OpenMP_CXX_FLAGS}
)

