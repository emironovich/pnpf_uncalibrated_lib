cmake_minimum_required(VERSION 3.0) 
project(p35p_double)

set(SOURCE_LIB  rtGetNaN.cpp
				rt_nonfinite.cpp
				p35p_solver.cpp
				rtGetInf.cpp

)
	
#include_directories(/home/elizaveta/Matlab/extern/include)

add_library(p35p_solver_double STATIC ${SOURCE_LIB})
