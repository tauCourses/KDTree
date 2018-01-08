cmake_minimum_required(VERSION 3.8)
project(KDTree)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(SOURCE_FILES main.cpp knn.h)
add_executable(KDTree ${SOURCE_FILES})