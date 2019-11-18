Pose+focal length solvers (P3.5P & P4Pf based on 3Q3)

# Building

```bash
git submodule update --init --recursive && mkdir -p build && cd build && cmake ../ -DCMAKE_BUILD_TYPE=Release && make -j
```

# Tests & benchmarks

```bash
 ./build/tests/pnpf_benchmark
 ./build/tests/pnpf_test
```
