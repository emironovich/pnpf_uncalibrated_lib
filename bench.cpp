#include "datagen.h"
#include "solvers.h"
#include <Eigen/Dense>
#include <benchmark/benchmark.h>

template <class Type> void BM_P35P(benchmark::State &state) {
  // allocate for generated data
  Eigen::Matrix<Type, 3, 4> points_3d;
  Eigen::Matrix<Type, 2, 4> points_2d;
  Type f_gen;
  Eigen::Matrix<Type, 3, 3> R_gen;
  Eigen::Matrix<Type, 3, 1> C_gen;

  // allocate for estimated data
  int solution_num = 0;
  Type f_sol_data[10];
  Eigen::Matrix<Type, 3, 3> R_sol_data[10];
  Eigen::Matrix<Type, 3, 1> T_sol_data[10];

  generateData(points_3d, points_2d, f_gen, R_gen, C_gen);
  P35PSolver<Type> p35pSolver;
  for (auto _ : state) {
    p35pSolver.solve(points_3d, points_2d, &solution_num, f_sol_data,
                     R_sol_data, T_sol_data);
  }
}
BENCHMARK_TEMPLATE(BM_P35P, double);
BENCHMARK_TEMPLATE1(BM_P35P, float);

template <class Type> void BM_P4P(benchmark::State &state) {
  Type eps = 1e-4;
  // allocate for generated data
  Eigen::Matrix<Type, 3, 4> points_3d;
  Eigen::Matrix<Type, 2, 4> points_2d;
  Type f_gen;
  Eigen::Matrix<Type, 3, 3> R_gen;
  Eigen::Matrix<Type, 3, 1> C_gen;

  // allocate for estimated data
  int solution_num = 0;
  Type f_sol_data[10];
  Eigen::Matrix<Type, 3, 3> R_sol_data[10];
  Eigen::Matrix<Type, 3, 1> T_sol_data[10];

  generateData(points_3d, points_2d, f_gen, R_gen, C_gen);
  P4PSolver<Type> p4pSolver;
  for (auto _ : state) {
    p4pSolver.solve(points_3d, points_2d, &solution_num, f_sol_data, R_sol_data,
                    T_sol_data);
  }
}
BENCHMARK_TEMPLATE(BM_P4P, double);
BENCHMARK_TEMPLATE1(BM_P4P, float);