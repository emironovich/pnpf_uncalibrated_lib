/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgerc.cpp
 *
 * Code generation for function 'xgerc'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xgerc.h"

/* Function Definitions */
void xgerc(int m, int n, double alpha1, int ix0, const double y[10], double A
           [100], int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i22;
  int i23;
  int ijA;
  if (!(alpha1 == 0.0)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i22 = jA + 1;
        i23 = m + jA;
        for (ijA = i22; ijA <= i23; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 10;
    }
  }
}

/* End of code generation (xgerc.cpp) */
