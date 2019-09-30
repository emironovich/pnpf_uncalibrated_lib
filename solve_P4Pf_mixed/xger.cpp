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
void b_xger(int m, int n, float alpha1, int ix0, const float y[3], float A[9],
            int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i21;
  int i22;
  int ijA;
  if (!(alpha1 == 0.0F)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i21 = jA + 1;
        i22 = m + jA;
        for (ijA = i21; ijA <= i22; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

void xger(int m, int n, float alpha1, int ix0, const float y[8], float A[64],
          int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i9;
  int i10;
  int ijA;
  if (!(alpha1 == 0.0F)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i9 = jA + 1;
        i10 = m + jA;
        for (ijA = i9; ijA <= i10; ijA++) {
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
