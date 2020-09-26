cmake_minimum_required(VERSION 2.8.3)
project(simplex_method)
ADD_COMPILE_OPTIONS(-std=c++14 )
find_package(Eigen3 3.3.0 REQUIRED
             PATHS ~/3rd-party/eigen-3.3.90/install/ # look here
             NO_DEFAULT_PATH)
add_library(simplex_method SHARED
    simplex_method.cpp)
add_executable(exe main.cpp)
target_link_libraries(exe simplex_method)
