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
void xgerc(int m, int n, float alpha1, int ix0, const float y[10], float A[100],
           int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i17;
  int i18;
  int ijA;
  if (!(alpha1 == 0.0F)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i17 = jA + 1;
        i18 = m + jA;
        for (ijA = i17; ijA <= i18; ijA++) {
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
