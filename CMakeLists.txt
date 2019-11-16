cmake_minimum_required(VERSION 3.0)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
project(pnp_test)

include_directories(datagen)

include_directories(pnpf)
add_subdirectory(pnpf)


set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
include_directories(benchmark/include)
add_subdirectory(benchmark)

include_directories(googletest/googletest/include)
add_subdirectory(googletest)
enable_testing()

add_executable(bench bench.cpp)
target_link_libraries(bench pnpf Eigen3::Eigen benchmark_main gtest_main gtest)

add_executable(pnpTest test.cpp)
target_link_libraries(pnpTest pnpf gtest_main Eigen3::Eigen gtest)
add_test(NAME pnpTest COMMAND pnpTest)

# add_executable(simple_test simple_test.cpp)
# target_link_libraries(simple_test gtest_main)	
# add_test(NAME simple_test COMMAND simple_test)