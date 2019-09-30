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
#include "p35p_solver.h"
#include "xger.h"

/* Function Definitions */
void xger(int m, int n, double alpha1, int ix0, const double y[4], double A[12],
          int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i35;
  int i36;
  int ijA;
  if (!(alpha1 == 0.0)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i35 = jA + 1;
        i36 = m + jA;
        for (ijA = i35; ijA <= i36; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

/* End of code generation (xger.cpp) */
