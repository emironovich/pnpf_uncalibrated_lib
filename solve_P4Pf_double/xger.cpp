/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xger.cpp
 *
 * Code generation for function 'xger'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xger.h"

/* Function Definitions */
void b_xger(int m, int n, double alpha1, int ix0, const double y[3], double A[9],
            int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i32;
  int i33;
  int ijA;
  if (!(alpha1 == 0.0)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i32 = jA + 1;
        i33 = m + jA;
        for (ijA = i32; ijA <= i33; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

void xger(int m, int n, double alpha1, int ix0, const double y[8], double A[64],
          int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i15;
  int i16;
  int ijA;
  if (!(alpha1 == 0.0)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i15 = jA + 1;
        i16 = m + jA;
        for (ijA = i15; ijA <= i16; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 8;
    }
  }
}

/* End of code generation (xger.cpp) */
