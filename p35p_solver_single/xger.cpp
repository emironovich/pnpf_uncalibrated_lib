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
void xger(int m, int n, float alpha1, int ix0, const float y[4], float A[12],
          int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i29;
  int i30;
  int ijA;
  if (!(alpha1 == 0.0F)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i29 = jA + 1;
        i30 = m + jA;
        for (ijA = i29; ijA <= i30; ijA++) {
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
