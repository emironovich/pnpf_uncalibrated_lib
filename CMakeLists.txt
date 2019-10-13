cmake_minimum_required(VERSION 3.0)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
project(pnp_test)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Now simply link against gtest or gtest_main as needed. Eg
include_directories(datagen)
add_subdirectory(datagen)

include_directories(p35p_solver_double)
add_subdirectory(p35p_solver_double)

include_directories(p35p_solver_single)
add_subdirectory(p35p_solver_single)

include_directories(solve_P4Pf_double)
add_subdirectory(solve_P4Pf_double)

include_directories(solve_P4Pf_mixed)
add_subdirectory(solve_P4Pf_mixed)

 add_executable(pnpTest test.cpp accuracyTest.cpp)
 target_link_libraries(pnpTest datagen p35p_solver_double p35p_solver_single p4p_solver_double p4p_solver_mixed gtest_main Eigen3::Eigen)
 add_test(NAME pnpTest COMMAND pnpTest)

# add_executable(simple_test simple_test.cpp)
# target_link_libraries(simple_test gtest_main)	
# add_test(NAME simple_test COMMAND simple_test)