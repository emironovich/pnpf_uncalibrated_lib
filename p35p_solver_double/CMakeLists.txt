cmake_minimum_required(VERSION 3.0) 
project(p35p_double)

set(SOURCE_LIB  sqrt.cpp
				mult_poly42.cpp
				p35p_solver_initialize.cpp
				xzggbak.cpp
				xdlanv2.cpp
				p35p_solver_terminate.cpp
				mldivide.cpp
				xgerc.cpp
				tmp/CMakeFiles/3.15.3/CompilerIdCXX/CMakeCXXCompilerId.cpp
				xzhgeqz.cpp
				xzlartg.cpp
				rtGetNaN.cpp
				p35p_solver_data.cpp
				xrot.cpp
				xdhseqr.cpp
				mult_poly22.cpp
				rt_nonfinite.cpp
				xzggev.cpp
				xzgghrd.cpp
				cat.cpp
				qr.cpp
				equations_for_groebner.cpp
				xnrm2.cpp
				xzlarf.cpp
				mult_for_groebner.cpp
				eig.cpp
				xzggbal.cpp
				p35p_solver.cpp
				find_f.cpp
				xztgevc.cpp
				rtGetInf.cpp
				schur.cpp
				quadruple_constraint.cpp
				p35p_solver_rtwutil.cpp
)
	
include_directories(/home/elizaveta/Matlab/extern/include)

add_library(p35p_solver_double STATIC ${SOURCE_LIB})
