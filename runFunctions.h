/*
#ifndef PNP_TEST_TEST_H
#define PNP_TEST_TEST_H

#endif //PNP_TEST_TEST*/

#include "datagen.h"
#include "pnpf_solvers.h"
#include <Eigen/Dense>
#include <climits>
#include <cmath>
#include <iostream>

using namespace Eigen;

struct TestResult {
  int existSolutions;
  int belowThreshold;
};

template <class Type>
TestResult runFunction(const char *funcName, Type eps, int it_num = 1) {
  int succ_num = 0;
  int zero_solutions_num = 0;

  // allocate for generated data
  Matrix<Type, 3, 4> points_3d;
  Matrix<Type, 2, 4> points_2d;
  Type f_gen;
  Matrix<Type, 3, 3> R_gen;
  Matrix<Type, 3, 1> C_gen;

  // allocate for estimated data
  int solution_num = 0;
  Type fs[10];
  Matrix<Type, 3, 3> Rs[10];
  Matrix<Type, 3, 1> Cs[10];
  // p4p_solver_initialize();

  for (int curr_it = 0; curr_it < it_num; ++curr_it) {
    generateData(points_3d, points_2d, f_gen, R_gen, C_gen);
    if (strcmp(funcName, "p4p") == 0)
      p4p_solver(points_3d, points_2d, eps, &solution_num, fs, Rs, Cs);
    else if (strcmp(funcName, "p3.5p") == 0)
      p35p_solver(points_3d, points_2d, eps, &solution_num, fs, Rs, Cs);
    else {
      solution_num = 0; // todo: maybe add a warning
    }
    // allocate for comparison
    Type min_diff = std::numeric_limits<Type>::max();
    Type diff_C = std::numeric_limits<Type>::max();
    Type diff_R = std::numeric_limits<Type>::max();

    for (int i = 0; i < solution_num; ++i) {

      Type f = fs[i];
      Matrix<Type, 3, 1> C;
      Matrix<Type, 3, 3> R;
      C = Cs[i];
      R = Rs[i];
      Matrix<Type, 3, 3> K;
      K.setZero();
      K(0, 0) = K(1, 1) = f;
      K(2, 2) = 1;

      Type diff = abs(f - f_gen);

      Type diffR = (R - R_gen).norm() / 3;
      Type diffC = (C - C_gen).norm();
      if (diff < min_diff) {
        min_diff = diff;
        diff_R = diffR;
        diff_C = diffC;
      } else if (diff == min_diff && diffR < diff_R) {
        diff_R = diffR;
        diff_C = diffC;
      }
    }
    if (min_diff / abs(f_gen) < 0.01 && diff_R < 0.03 &&
        diff_C / C_gen.norm() < 0.1 && solution_num != 0)
      succ_num++;
    if (solution_num == 0)
      zero_solutions_num++;
  }
  /*    if(strcmp(funcName, "p4p"))
          p4p_solver_terminate();
      else if(strcmp(funcName, "p3.5p"))
          p35p_solver_terminate();*/

  return {(it_num - zero_solutions_num), succ_num};
}