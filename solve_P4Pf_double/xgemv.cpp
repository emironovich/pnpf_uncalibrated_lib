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
void b_xgemv(int m, int n, const double A[9], int ia0, const double x[9], int
             ix0, double y[3])
{
  int iy;
  int i30;
  int iac;
  int ix;
  double c;
  int i31;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(double))));
    }

    iy = 0;
    i30 = ia0 + 3 * (n - 1);
    for (iac = ia0; iac <= i30; iac += 3) {
      ix = ix0;
      c = 0.0;
      i31 = (iac + m) - 1;
      for (ia = iac; ia <= i31; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

void xgemv(int m, int n, const double A[64], int ia0, const double x[64], int
           ix0, double y[8])
{
  int iy;
  int i13;
  int iac;
  int ix;
  double c;
  int i14;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(double))));
    }

    iy = 0;
    i13 = ia0 + ((n - 1) << 3);
    for (iac = ia0; iac <= i13; iac += 8) {
      ix = ix0;
      c = 0.0;
      i14 = (iac + m) - 1;
      for (ia = iac; ia <= i14; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

/* End of code generation (xgemv.cpp) */
