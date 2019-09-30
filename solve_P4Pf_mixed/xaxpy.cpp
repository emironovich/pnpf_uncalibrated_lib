/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xaxpy.cpp
 *
 * Code generation for function 'xaxpy'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xaxpy.h"

/* Function Definitions */
void b_xaxpy(int n, float a, const float x[8], int ix0, float y[2])
{
  if ((n < 1) || (a == 0.0F)) {
  } else {
    y[1] += a * x[ix0 - 1];
  }
}

void c_xaxpy(int n, float a, const float x[2], float y[8], int iy0)
{
  int iy;
  if ((n < 1) || (a == 0.0F)) {
  } else {
    iy = iy0 - 1;
    y[iy] += a * x[1];
  }
}

void xaxpy(int n, float a, int ix0, float y[8], int iy0)
{
  int ix;
  int iy;
  int i25;
  int k;
  if (!(a == 0.0F)) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    i25 = n - 1;
    for (k = 0; k <= i25; k++) {
      y[iy] += a * y[ix];
      ix++;
      iy++;
    }
  }
}

/* End of code generation (xaxpy.cpp) */
