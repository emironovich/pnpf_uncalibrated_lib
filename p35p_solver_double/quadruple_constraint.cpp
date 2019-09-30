/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadruple_constraint.cpp
 *
 * Code generation for function 'quadruple_constraint'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "quadruple_constraint.h"

/* Function Definitions */
void quadruple_constraint(double i, double j, double k, const double x[4], const
  double y[4], const double X[12], const double R[54], double F_row[18])
{
  int F_row_tmp_tmp_tmp_tmp;
  int b_i;
  int F_row_tmp_tmp_tmp;
  int F_row_tmp_tmp;
  double F_row_tmp[3];
  int b_F_row_tmp_tmp_tmp;
  int b_F_row_tmp_tmp;
  double b_F_row_tmp[3];
  double c_F_row_tmp_tmp;
  double sum[6];
  double b_sum[6];
  double a_tmp;
  double a;

  /* fc: */
  F_row_tmp_tmp_tmp_tmp = static_cast<int>(i) - 1;
  b_i = 3 * F_row_tmp_tmp_tmp_tmp;
  F_row_tmp_tmp_tmp = static_cast<int>(j) - 1;
  F_row_tmp_tmp = 3 * F_row_tmp_tmp_tmp;
  F_row_tmp[0] = X[F_row_tmp_tmp] - X[b_i];
  b_F_row_tmp_tmp_tmp = static_cast<int>(k) - 1;
  b_F_row_tmp_tmp = 3 * b_F_row_tmp_tmp_tmp;
  b_F_row_tmp[0] = X[b_F_row_tmp_tmp] - X[b_i];
  c_F_row_tmp_tmp = X[1 + b_i];
  F_row_tmp[1] = X[1 + F_row_tmp_tmp] - c_F_row_tmp_tmp;
  b_F_row_tmp[1] = X[1 + b_F_row_tmp_tmp] - c_F_row_tmp_tmp;
  c_F_row_tmp_tmp = X[2 + b_i];
  F_row_tmp[2] = X[2 + F_row_tmp_tmp] - c_F_row_tmp_tmp;
  b_F_row_tmp[2] = X[2 + b_F_row_tmp_tmp] - c_F_row_tmp_tmp;
  for (b_i = 0; b_i < 6; b_i++) {
    sum[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      sum[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * F_row_tmp[b_i];
    }
  }

  c_F_row_tmp_tmp = y[F_row_tmp_tmp_tmp_tmp] - y[b_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    b_sum[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_sum[F_row_tmp_tmp] += R[1 + (3 * b_i + 9 * F_row_tmp_tmp)] *
        b_F_row_tmp[b_i];
    }
  }

  a_tmp = x[F_row_tmp_tmp_tmp_tmp] - x[F_row_tmp_tmp_tmp];

  /* fs: */
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i] = c_F_row_tmp_tmp * sum[b_i] - a_tmp * b_sum[b_i];
    sum[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      sum[F_row_tmp_tmp] += R[1 + (3 * b_i + 9 * F_row_tmp_tmp)] * F_row_tmp[b_i];
    }
  }

  for (b_i = 0; b_i < 6; b_i++) {
    b_sum[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_sum[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * b_F_row_tmp[b_i];
    }
  }

  /* 1: */
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[1 + 3 * b_i] = -c_F_row_tmp_tmp * sum[b_i] - a_tmp * b_sum[b_i];
    sum[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      sum[F_row_tmp_tmp] += R[2 + (3 * b_i + 9 * F_row_tmp_tmp)] * F_row_tmp[b_i];
    }
  }

  a = -(y[static_cast<int>(i) - 1] - y[static_cast<int>(k) - 1]) * x[static_cast<
    int>(j) - 1];
  for (b_i = 0; b_i < 6; b_i++) {
    b_sum[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_sum[F_row_tmp_tmp] += R[2 + (3 * b_i + 9 * F_row_tmp_tmp)] *
        b_F_row_tmp[b_i];
    }
  }

  c_F_row_tmp_tmp = a_tmp * y[static_cast<int>(k) - 1];
  for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
    F_row[2 + 3 * F_row_tmp_tmp] = a * sum[F_row_tmp_tmp] + c_F_row_tmp_tmp *
      b_sum[F_row_tmp_tmp];
  }
}

/* End of code generation (quadruple_constraint.cpp) */
