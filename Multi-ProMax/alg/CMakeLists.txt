cmake_minimum_required(VERSION 3.10)
project(alg)

set(CMAKE_CXX_STANDARD 14)

add_executable(alg pm.cpp company.h backup/infgraph2_v1.h backup/infgraph,h_v2_gamma.h backup/pm.h backup/update-time-tag.h)
