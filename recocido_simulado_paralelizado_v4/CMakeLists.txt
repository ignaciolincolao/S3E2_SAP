cmake_minimum_required(VERSION 3.15)
project(paralelizado LANGUAGES CUDA CXX)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_VERBOSE_MAKEFILE TRUE)



add_executable(paralelizado main.cu)

target_link_libraries(paralelizado  ${CMAKE_SOURCE_DIR}/UTM.cpp)