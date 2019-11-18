#ifndef TEST_HELPER_HPP
#define TEST_HELPER_HPP

#include <Eigen/Dense>
#include <climits>
#include <cmath>
#include <iostream>

#include "pnpf/solvers.hpp"

#include "data_generation.hpp"

struct TestResult {
  int existSolutions;
  int belowThreshold;
};

template <typename Solver> TestResult runFunction(Solver &solver, int it_num) {
  using Scalar = typename Solver::Scalar;
  using Matrix34 = typename Solver::Matrix34;
  using Matrix24 = typename Solver::Matrix24;
  using Matrix33 = typename Solver::Matrix33;
  using Vector3 = typename Solver::Vector3;
  const int MaxSolutions = Solver::MaxSolutions;

  const double f_tolerance = 0.01;
  const double R_tolerance = 0.03;
  const double C_tolerance = 0.1;

  int succ_num = 0;
  int zero_solutions_num = 0;

  // allocate for generated data
  Matrix34 points_3d;
  Matrix24 points_2d;
  Scalar f_gen;
  Matrix33 R_gen;
  Vector3 C_gen;

  // allocate for estimated data
  int solution_num = 0;
  Scalar fs[MaxSolutions];
  Matrix33 Rs[MaxSolutions];
  Vector3 Cs[MaxSolutions];
  // p4p_solver_initialize();

  for (int curr_it = 0; curr_it < it_num; ++curr_it) {
    generateData(points_3d, points_2d, f_gen, R_gen, C_gen);
    Scalar diag =
        (points_2d.rowwise().maxCoeff() - points_2d.rowwise().minCoeff())
            .norm();
    solver.solve(points_3d, points_2d, &solution_num, fs, Rs, Cs, diag);
    // allocate for comparison
    Scalar min_diff, diff_C, diff_R;
    min_diff = diff_C = diff_R = std::numeric_limits<Scalar>::max();

    for (int i = 0; i < solution_num; ++i) {
      Scalar f = fs[i];
      Vector3 C;
      Matrix33 R;
      C = Cs[i];
      R = Rs[i];

      Scalar diff = abs(f - f_gen);

      Scalar diffR = (R - R_gen).norm() / 3.;
      Scalar diffC = (C - C_gen).norm();
      if (diff < min_diff) {
        min_diff = diff;
        diff_R = diffR;
        diff_C = diffC;
      } else if (diff == min_diff && diffR < diff_R) {
        diff_R = diffR;
        diff_C = diffC;
      }
    }
    if (min_diff / abs(f_gen) < f_tolerance && diff_R < R_tolerance &&
        diff_C / C_gen.norm() < C_tolerance && solution_num != 0)
      succ_num++;
    if (solution_num == 0)
      zero_solutions_num++;
  }
  return {(it_num - zero_solutions_num), succ_num};
}

#endif // PNP_TEST_TEST
