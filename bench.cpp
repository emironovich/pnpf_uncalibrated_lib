#include "datagen.h"
#include "pnpf_solvers.h"
#include <Eigen/Dense>
#include <benchmark/benchmark.h>

template <class Type> void BM_P35P(benchmark::State &state) {
  Type eps = 1e-4;
  // allocate for generated data
  Matrix<Type, 3, 4> points_3d;
  Matrix<Type, 2, 4> points_2d;
  Type f_gen;
  Matrix<Type, 3, 3> R_gen;
  Matrix<Type, 3, 1> C_gen;

  // allocate for estimated data
  int solution_num = 0;
  Type f_sol_data[10];
  Matrix<Type, 3, 3> R_sol_data[10];
  Matrix<Type, 3, 1> T_sol_data[10];

  generateData(points_3d, points_2d, f_gen, R_gen, C_gen);

  for (auto _ : state) {
    p35p_solver(points_3d, points_2d, eps, &solution_num, f_sol_data,
                R_sol_data, T_sol_data);
  }
}
BENCHMARK_TEMPLATE(BM_P35P, double);
BENCHMARK_TEMPLATE1(BM_P35P, float);

template <class Type> void BM_P4P(benchmark::State &state) {
  Type eps = 1e-4;
  // allocate for generated data
  Matrix<Type, 3, 4> points_3d;
  Matrix<Type, 2, 4> points_2d;
  Type f_gen;
  Matrix<Type, 3, 3> R_gen;
  Matrix<Type, 3, 1> C_gen;

  // allocate for estimated data
  int solution_num = 0;
  Type f_sol_data[10];
  Matrix<Type, 3, 3> R_sol_data[10];
  Matrix<Type, 3, 1> T_sol_data[10];

  generateData(points_3d, points_2d, f_gen, R_gen, C_gen);

  for (auto _ : state) {
    p4p_solver(points_3d, points_2d, eps, &solution_num, f_sol_data, R_sol_data,
               T_sol_data);
  }
}
BENCHMARK_TEMPLATE(BM_P4P, double);
BENCHMARK_TEMPLATE1(BM_P4P, float);