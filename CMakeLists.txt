cmake_minimum_required(VERSION 3.0)
find_package(Eigen3 3.3 REQUIRED CONFIG)
project(pnpf)

add_compile_options(-march=native -Wall -Wextra -Werror)
add_compile_options(
    $<$<OR:$<CONFIG:RELEASE>,$<CONFIG:RELWITHDEBINFO>>:-O3>
)

if (NOT TARGET gtest_main)
    add_subdirectory(googletest)
    message(STATUS "Adding local google-test build")
endif()
if (NOT TARGET benchmark_main)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
    add_subdirectory(benchmark)
    message(STATUS "Adding local google-benchmark build")
endif()
enable_testing()

add_subdirectory(pnpf)
add_subdirectory(tests)
