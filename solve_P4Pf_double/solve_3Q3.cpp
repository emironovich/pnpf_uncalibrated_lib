/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * solve_3Q3.cpp
 *
 * Code generation for function 'solve_3Q3'
 *
 */

/* Include files */
#include <cmath>
#include <math.h>
#include "rt_nonfinite.h"
#include <string.h>
#include "solve_P4Pf.h"
#include "solve_3Q3.h"
#include "qr.h"
#include "roots.h"
#include "find_M.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d6;
  double d7;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d6 = std::abs(u0);
    d7 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d6 == 1.0) {
        y = 1.0;
      } else if (d6 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d7 == 0.0) {
      y = 1.0;
    } else if (d7 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void solve_3Q3(const double c[30], double e, double *n, double xs_data[], int
               xs_size[2], double ys_data[], int ys_size[2], double zs_data[],
               int zs_size[2])
{
  double A[9];
  int i;
  double a21;
  double P[27];
  double P_prime[27];
  double M[45];
  double b_A[9];
  int r1;
  int k;
  int r2;
  int r3;
  double maxval;
  int rtemp;
  int i5;
  double C[13];
  double y[20];
  double d1;
  double b_C[13];
  double d2;
  double d3;
  double d4;
  double d5;
  double c_C[13];
  double d_C[13];
  creal_T xs_complex_data[12];
  int xs_complex_size[1];
  double mons[5];
  double c_A[9];

  /* c -- 3x10 coefficients matrix */
  /* SOLVE_3Q3 Summary of this function goes here */
  /*    Detailed explanation goes here */
  A[0] = c[3];
  A[3] = c[6];
  A[6] = c[15];
  A[1] = c[4];
  A[4] = c[7];
  A[7] = c[16];
  A[2] = c[5];
  A[5] = c[8];
  A[8] = c[17];

  /* [x^2, x, 1] */
  for (i = 0; i < 3; i++) {
    P[i] = 0.0;
    P[i + 9] = -c[9 + i];
    P[i + 18] = -c[21 + i];

    /* y */
    P[3 + i] = 0.0;
    P[i + 12] = -c[12 + i];
    P[i + 21] = -c[24 + i];

    /* z */
    P[6 + i] = -c[i];
    P[i + 15] = -c[18 + i];
    P[i + 24] = -c[27 + i];

    /* 1 */
  }

  a21 = std::abs(c[4]);
  for (i = 0; i < 3; i++) {
    memcpy(&b_A[0], &A[0], 9U * sizeof(double));
    r1 = 0;
    r2 = 1;
    r3 = 2;
    maxval = std::abs(A[0]);
    if (a21 > maxval) {
      maxval = a21;
      r1 = 1;
      r2 = 0;
    }

    if (std::abs(A[2]) > maxval) {
      r1 = 2;
      r2 = 1;
      r3 = 0;
    }

    b_A[r2] = A[r2] / A[r1];
    b_A[r3] /= b_A[r1];
    b_A[3 + r2] -= b_A[r2] * b_A[3 + r1];
    b_A[3 + r3] -= b_A[r3] * b_A[3 + r1];
    b_A[6 + r2] -= b_A[r2] * b_A[6 + r1];
    b_A[6 + r3] -= b_A[r3] * b_A[6 + r1];
    if (std::abs(b_A[3 + r3]) > std::abs(b_A[3 + r2])) {
      rtemp = r2;
      r2 = r3;
      r3 = rtemp;
    }

    b_A[3 + r3] /= b_A[3 + r2];
    b_A[6 + r3] -= b_A[3 + r3] * b_A[6 + r2];
    i5 = r1 + 9 * i;
    rtemp = r2 + 9 * i;
    maxval = P[rtemp] - P[i5] * b_A[r2];
    k = r3 + 9 * i;
    d1 = ((P[k] - P[i5] * b_A[r3]) - maxval * b_A[3 + r3]) / b_A[6 + r3];
    P_prime[2 + 9 * i] = d1;
    maxval -= d1 * b_A[6 + r2];
    maxval /= b_A[3 + r2];
    P_prime[1 + 9 * i] = maxval;
    P_prime[9 * i] = ((P[i5] - d1 * b_A[6 + r1]) - maxval * b_A[3 + r1]) /
      b_A[r1];
    d2 = P[i5 + 3];
    maxval = P[rtemp + 3] - d2 * b_A[r2];
    d3 = b_A[3 + r3];
    d4 = b_A[6 + r3];
    d1 = ((P[k + 3] - d2 * b_A[r3]) - maxval * d3) / d4;
    P_prime[9 * i + 5] = d1;
    d5 = b_A[6 + r1];
    d2 -= d1 * d5;
    maxval -= d1 * b_A[6 + r2];
    maxval /= b_A[3 + r2];
    P_prime[9 * i + 4] = maxval;
    d2 -= maxval * b_A[3 + r1];
    d2 /= b_A[r1];
    P_prime[3 + 9 * i] = d2;
    d2 = P[i5 + 6];
    maxval = P[rtemp + 6] - d2 * b_A[r2];
    d1 = ((P[k + 6] - d2 * b_A[r3]) - maxval * d3) / d4;
    P_prime[9 * i + 8] = d1;
    d2 -= d1 * d5;
    maxval -= d1 * b_A[6 + r2];
    maxval /= b_A[3 + r2];
    P_prime[9 * i + 7] = maxval;
    d2 -= maxval * b_A[3 + r1];
    d2 /= b_A[r1];
    P_prime[6 + 9 * i] = d2;
  }

  find_M(P_prime, M);
  memset(&A[0], 0, 9U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 5; r1++) {
      rtemp = k + r1;
      A[rtemp] += M[40 + k] * M[20 + r1];
    }
  }

  memset(&b_A[0], 0, 9U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 5; r1++) {
      rtemp = k + r1;
      b_A[rtemp] += M[25 + k] * M[35 + r1];
    }
  }

  for (i5 = 0; i5 < 9; i5++) {
    A[i5] -= b_A[i5];
  }

  memset(&C[0], 0, 13U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 9; r1++) {
      rtemp = k + r1;
      C[rtemp] += M[k] * A[r1];
    }
  }

  rtemp = -1;
  for (r1 = 0; r1 < 10; r1++) {
    rtemp++;
    y[rtemp] = M[r1 % 5 + 5 * (1 + r1 / 5)];
  }

  for (r1 = 0; r1 < 10; r1++) {
    rtemp++;
    y[rtemp] = M[30 + (r1 % 5 + 5 * (1 + r1 / 5))];
  }

  memset(&A[0], 0, 9U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 5; r1++) {
      rtemp = k + r1;
      A[rtemp] += y[15 + k] * y[r1];
    }
  }

  memset(&b_A[0], 0, 9U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 5; r1++) {
      rtemp = k + r1;
      b_A[rtemp] += y[5 + k] * y[10 + r1];
    }
  }

  for (i5 = 0; i5 < 9; i5++) {
    A[i5] -= b_A[i5];
  }

  memset(&b_C[0], 0, 13U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 9; r1++) {
      rtemp = k + r1;
      b_C[rtemp] += M[15 + k] * A[r1];
    }
  }

  memset(&A[0], 0, 9U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 5; r1++) {
      rtemp = k + r1;
      A[rtemp] += M[25 + k] * M[5 + r1];
    }
  }

  memset(&b_A[0], 0, 9U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 5; r1++) {
      rtemp = k + r1;
      b_A[rtemp] += M[10 + k] * M[20 + r1];
    }
  }

  for (i5 = 0; i5 < 9; i5++) {
    A[i5] -= b_A[i5];
  }

  memset(&c_C[0], 0, 13U * sizeof(double));
  for (k = 0; k < 5; k++) {
    for (r1 = 0; r1 < 9; r1++) {
      rtemp = k + r1;
      c_C[rtemp] += M[30 + k] * A[r1];
    }
  }

  for (i5 = 0; i5 < 13; i5++) {
    d_C[i5] = (C[i5] - b_C[i5]) + c_C[i5];
  }

  roots(d_C, xs_complex_data, xs_complex_size);
  xs_size[0] = 1;
  if (0 <= xs_complex_size[0] - 1) {
    memset(&xs_data[0], 0, (unsigned int)(xs_complex_size[0] * static_cast<int>
            (sizeof(double))));
  }

  *n = 0.0;
  i5 = xs_complex_size[0];
  for (i = 0; i < i5; i++) {
    if (std::abs(xs_complex_data[i].im) < e) {
      (*n)++;
      xs_data[static_cast<int>(*n) - 1] = xs_complex_data[i].re;
    }
  }

  if (1.0 > *n) {
    xs_size[1] = 0;
  } else {
    xs_size[1] = static_cast<int>(*n);
  }

  i5 = static_cast<int>(*n);
  ys_size[0] = 1;
  rtemp = static_cast<signed char>(*n);
  ys_size[1] = rtemp;
  zs_size[0] = 1;
  zs_size[1] = rtemp;
  if (0 <= i5 - 1) {
    mons[4] = 1.0;
  }

  for (i = 0; i < i5; i++) {
    mons[0] = rt_powd_snf(xs_data[i], 4.0);
    mons[1] = rt_powd_snf(xs_data[i], 3.0);
    mons[2] = xs_data[i] * xs_data[i];
    mons[3] = xs_data[i];
    for (k = 0; k < 3; k++) {
      for (r1 = 0; r1 < 3; r1++) {
        maxval = 0.0;
        for (rtemp = 0; rtemp < 5; rtemp++) {
          maxval += mons[rtemp] * M[(rtemp + 5 * k) + 15 * r1];
        }

        c_A[r1 + 3 * k] = maxval;
      }
    }

    b_qr(c_A, A, b_A);
    ys_data[i] = A[6] / A[8];
    zs_data[i] = A[7] / A[8];
  }
}

/* End of code generation (solve_3Q3.cpp) */
