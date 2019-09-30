/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgemv.cpp
 *
 * Code generation for function 'xgemv'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xgemv.h"

/* Function Definitions */
void b_xgemv(int m, int n, const float A[9], int ia0, const float x[9], int ix0,
             float y[3])
{
  int iy;
  int i19;
  int iac;
  int ix;
  float c;
  int i20;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(float))));
    }

    iy = 0;
    i19 = ia0 + 3 * (n - 1);
    for (iac = ia0; iac <= i19; iac += 3) {
      ix = ix0;
      c = 0.0F;
      i20 = (iac + m) - 1;
      for (ia = iac; ia <= i20; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

void xgemv(int m, int n, const float A[64], int ia0, const float x[64], int ix0,
           float y[8])
{
  int iy;
  int i7;
  int iac;
  int ix;
  float c;
  int i8;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(float))));
    }

    iy = 0;
    i7 = ia0 + ((n - 1) << 3);
    for (iac = ia0; iac <= i7; iac += 8) {
      ix = ix0;
      c = 0.0F;
      i8 = (iac + m) - 1;
      for (ia = iac; ia <= i8; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

/* End of code generation (xgemv.cpp) */
