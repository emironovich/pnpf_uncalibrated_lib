#include <Eigen/Dense>
#include <benchmark/benchmark.h>
#include <iostream>

#include "pnpf/solvers.hpp"

#include "data_generation.hpp"

using namespace pnpf;

template <class Solver> void BM_PNPf(benchmark::State &state) {
  using Type = typename Solver::Scalar;
  using Matrix34 = typename Solver::Matrix34;
  using Matrix24 = typename Solver::Matrix24;
  using Matrix33 = typename Solver::Matrix33;
  using Vector3 = typename Solver::Vector3;
  const int MaxSolutions = Solver::MaxSolutions;
  // allocate for generated data
  Matrix34 points_3d;
  Matrix24 points_2d;
  Type f_gen;
  Matrix33 R_gen;
  Vector3 C_gen;

  // allocate for estimated data
  int solution_num = 0;
  Type f_sol_data[MaxSolutions];
  Matrix33 R_sol_data[MaxSolutions];
  Vector3 T_sol_data[MaxSolutions];

  generateData(points_3d, points_2d, f_gen, R_gen, C_gen);
  Solver solver;
  for (auto _ : state) {
    solver.solve(points_3d, points_2d, &solution_num, f_sol_data, R_sol_data,
                 T_sol_data);
  }
  std::cout << "Solutions: " << solution_num << std::endl;
}
BENCHMARK_TEMPLATE(BM_PNPf, P35PSolver<double>);
BENCHMARK_TEMPLATE(BM_PNPf, P35PSolver<float>);
BENCHMARK_TEMPLATE(BM_PNPf, P4PSolver<double>);
BENCHMARK_TEMPLATE(BM_PNPf, P4PSolver<float>);
