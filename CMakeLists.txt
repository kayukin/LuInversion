cmake_minimum_required(VERSION 3.3)
project(LU_Inversion)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp Matrix.cpp Matrix.h Generator.cpp Generator.h Vector.cpp Vector.h LuDecomposition.cpp LuDecomposition.h)
add_executable(LU_Inversion ${SOURCE_FILES})