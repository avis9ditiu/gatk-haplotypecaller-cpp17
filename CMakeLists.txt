cmake_minimum_required(VERSION 3.12)
project(Cigar)

set(CMAKE_CXX_STANDARD 17)

add_executable(Cigar main.cpp cigar.hpp smithwaterman.hpp)

target_link_libraries(Cigar gtest)