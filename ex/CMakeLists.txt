cmake_minimum_required (VERSION 3.0)
project (myproject)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(/home/elizaveta/Matlab/extern/include)

set(SOURCE_EXE main.cpp)

include_directories(p35p_solver_double)	
add_executable(main ${SOURCE_EXE})

add_subdirectory(p35p_solver_double)

target_link_libraries(main p35p_solver_double Eigen3::Eigen)	