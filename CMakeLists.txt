cmake_minimum_required(VERSION 3.12)
project(gatk)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS " -mavx -mavx2 -lomp -O3 -lboost_program_options ")

add_executable(gatk src/main.cpp)
