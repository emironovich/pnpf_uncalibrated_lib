cmake_minimum_required(VERSION 3.0) 
project(pnpf)

set(SOURCE_LIB pnpf.cpp)

add_library(pnpf STATIC ${SOURCE_LIB})
