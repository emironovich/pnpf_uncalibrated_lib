//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: quadruple_constraint.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "quadruple_constraint.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : double i
//                double j
//                double k
//                const double x[4]
//                const double y[4]
//                const double X[12]
//                const double R[54]
//                double F_row[18]
// Return Type  : void
//
void quadruple_constraint(double i, double j, double k, const double x[4], const
  double y[4], const double X[12], const double R[54], double F_row[18])
{
  int F_row_tmp_tmp_tmp;
  int b_i;
  int b_F_row_tmp_tmp_tmp;
  int F_row_tmp_tmp;
  double F_row_tmp[3];
  int c_F_row_tmp_tmp_tmp;
  int b_F_row_tmp_tmp;
  double b_F_row_tmp[3];
  int c_F_row_tmp_tmp;
  double a_tmp_tmp;
  double b[6];
  double a_tmp;
  double b_b[6];
  double a;

  // fc:
  F_row_tmp_tmp_tmp = static_cast<int>(i) - 1;
  b_i = 3 * F_row_tmp_tmp_tmp;
  b_F_row_tmp_tmp_tmp = static_cast<int>(j) - 1;
  F_row_tmp_tmp = 3 * b_F_row_tmp_tmp_tmp;
  F_row_tmp[0] = X[F_row_tmp_tmp] - X[b_i];
  c_F_row_tmp_tmp_tmp = static_cast<int>(k) - 1;
  b_F_row_tmp_tmp = 3 * c_F_row_tmp_tmp_tmp;
  b_F_row_tmp[0] = X[b_F_row_tmp_tmp] - X[b_i];
  c_F_row_tmp_tmp = b_i + 1;
  F_row_tmp[1] = X[F_row_tmp_tmp + 1] - X[c_F_row_tmp_tmp];
  b_F_row_tmp[1] = X[b_F_row_tmp_tmp + 1] - X[c_F_row_tmp_tmp];
  b_i += 2;
  F_row_tmp[2] = X[F_row_tmp_tmp + 2] - X[b_i];
  b_F_row_tmp[2] = X[b_F_row_tmp_tmp + 2] - X[b_i];
  a_tmp_tmp = y[F_row_tmp_tmp_tmp] - y[c_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * F_row_tmp[b_i];
    }
  }

  a_tmp = x[F_row_tmp_tmp_tmp] - x[b_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] *
        b_F_row_tmp[b_i];
    }
  }

  // fs:
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i] = a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] * F_row_tmp[b_i];
    }
  }

  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * b_F_row_tmp[b_i];
    }
  }

  // 1:
  a = -a_tmp_tmp * x[b_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i + 1] = -a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 2] * F_row_tmp[b_i];
    }
  }

  a_tmp_tmp = a_tmp * y[c_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 2] *
        b_F_row_tmp[b_i];
    }
  }

  for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
    F_row[3 * F_row_tmp_tmp + 2] = a * b[F_row_tmp_tmp] + a_tmp_tmp *
      b_b[F_row_tmp_tmp];
  }
}

//
// File trailer for quadruple_constraint.cpp
//
// [EOF]
//
