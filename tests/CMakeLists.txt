add_executable(pnpf_benchmark pnpf_benchmark.cpp)
target_link_libraries(pnpf_benchmark pnpf benchmark_main)

add_executable(pnpf_test pnpf_test.cpp)
target_link_libraries(pnpf_test pnpf gtest_main)
add_test(NAME pnpf COMMAND pnpf_test)

