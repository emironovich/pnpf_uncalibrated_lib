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
#include "p35p_solver.h"
#include "xgemv.h"

/* Function Definitions */
void b_xgemv(int m, int n, const float A[12], int ia0, const float x[12], int
             ix0, float y[4])
{
  int iy;
  int i27;
  int iac;
  int ix;
  float c;
  int i28;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(float))));
    }

    iy = 0;
    i27 = ia0 + 3 * (n - 1);
    for (iac = ia0; iac <= i27; iac += 3) {
      ix = ix0;
      c = 0.0F;
      i28 = (iac + m) - 1;
      for (ia = iac; ia <= i28; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

void xgemv(int m, int n, const float A[100], int ia0, const float x[100], int
           ix0, float y[10])
{
  int iy;
  int i15;
  int iac;
  int ix;
  float c;
  int i16;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(float))));
    }

    iy = 0;
    i15 = ia0 + 10 * (n - 1);
    for (iac = ia0; iac <= i15; iac += 10) {
      ix = ix0;
      c = 0.0F;
      i16 = (iac + m) - 1;
      for (ia = iac; ia <= i16; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

/* End of code generation (xgemv.cpp) */
