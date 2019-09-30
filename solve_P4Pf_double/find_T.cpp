/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find_T.cpp
 *
 * Code generation for function 'find_T'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "find_T.h"
#include "xgeqp3.h"

/* Function Definitions */
void find_T(const double X[12], const double u[4], const double v[4], const
            double R[9], double w, double T[3])
{
  double A[36];
  double b[12];
  int i;
  double A_tmp[3];
  int jpvt[3];
  int i11;
  int rankR;
  double tol;
  int j;
  int T_tmp;
  int b_A_tmp;
  double b_A[9];
  double c_A[3];
  memset(&A[0], 0, 36U * sizeof(double));
  memset(&b[0], 0, 12U * sizeof(double));
  for (i = 0; i < 4; i++) {
    i11 = 3 * (1 + i);
    A_tmp[0] = static_cast<double>(i11) + -2.0;
    A_tmp[1] = static_cast<double>(i11) + -1.0;
    A_tmp[2] = i11;
    A[(i11 + -2) - 1] = 0.0;
    A[(i11 + -2) + 11] = -1.0;
    A[(i11 + -2) + 23] = w * v[i];
    A[(i11 + -1) - 1] = 1.0;
    A[(i11 + -1) + 11] = 0.0;
    A[(i11 + -1) + 23] = -w * u[i];
    A[i11 - 1] = -v[i];
    A[i11 + 11] = u[i];
    A[i11 + 23] = 0.0;
    i11 -= 3;
    for (T_tmp = 0; T_tmp < 3; T_tmp++) {
      b_A_tmp = i11 + 12 * T_tmp;
      b_A[3 * T_tmp] = -A[b_A_tmp];
      b_A[1 + 3 * T_tmp] = -A[b_A_tmp + 1];
      b_A[2 + 3 * T_tmp] = -A[b_A_tmp + 2];
    }

    for (i11 = 0; i11 < 3; i11++) {
      c_A[i11] = 0.0;
      for (T_tmp = 0; T_tmp < 3; T_tmp++) {
        c_A[i11] += ((b_A[i11] * R[3 * T_tmp] + b_A[i11 + 3] * R[1 + 3 * T_tmp])
                     + b_A[i11 + 6] * R[2 + 3 * T_tmp]) * X[T_tmp + 3 * i];
      }

      b[static_cast<int>(A_tmp[i11]) - 1] = c_A[i11];
    }
  }

  xgeqp3(A, A_tmp, jpvt);
  rankR = 0;
  tol = 2.6645352591003757E-14 * std::abs(A[0]);
  while ((rankR < 3) && (!(std::abs(A[rankR + 12 * rankR]) <= tol))) {
    rankR++;
  }

  for (j = 0; j < 3; j++) {
    T[j] = 0.0;
    if (A_tmp[j] != 0.0) {
      tol = b[j];
      i11 = j + 2;
      for (i = i11; i < 13; i++) {
        tol += A[(i + 12 * j) - 1] * b[i - 1];
      }

      tol *= A_tmp[j];
      if (tol != 0.0) {
        b[j] -= tol;
        i11 = j + 2;
        for (i = i11; i < 13; i++) {
          b[i - 1] -= A[(i + 12 * j) - 1] * tol;
        }
      }
    }
  }

  for (i = 0; i < rankR; i++) {
    T[jpvt[i] - 1] = b[i];
  }

  for (j = rankR; j >= 1; j--) {
    T_tmp = jpvt[j - 1] - 1;
    b_A_tmp = 12 * (j - 1);
    T[T_tmp] /= A[(j + b_A_tmp) - 1];
    for (i = 0; i <= j - 2; i++) {
      T[jpvt[i] - 1] -= T[T_tmp] * A[i + b_A_tmp];
    }
  }
}

/* End of code generation (find_T.cpp) */
