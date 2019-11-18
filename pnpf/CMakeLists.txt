cmake_minimum_required(VERSION 3.0) 

set(SOURCE_LIB src/pnpf.cpp include/pnpf/solvers.hpp)

add_library(pnpf STATIC ${SOURCE_LIB})
target_include_directories(pnpf PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(pnpf Eigen3::Eigen)
