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
void b_xgemv(int m, int n, const double A[12], int ia0, const double x[12], int
             ix0, double y[4])
{
  int iy;
  int i33;
  int iac;
  int ix;
  double c;
  int i34;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(double))));
    }

    iy = 0;
    i33 = ia0 + 3 * (n - 1);
    for (iac = ia0; iac <= i33; iac += 3) {
      ix = ix0;
      c = 0.0;
      i34 = (iac + m) - 1;
      for (ia = iac; ia <= i34; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

void xgemv(int m, int n, const double A[100], int ia0, const double x[100], int
           ix0, double y[10])
{
  int iy;
  int i20;
  int iac;
  int ix;
  double c;
  int i21;
  int ia;
  if (n != 0) {
    if (0 <= n - 1) {
      memset(&y[0], 0, (unsigned int)(n * static_cast<int>(sizeof(double))));
    }

    iy = 0;
    i20 = ia0 + 10 * (n - 1);
    for (iac = ia0; iac <= i20; iac += 10) {
      ix = ix0;
      c = 0.0;
      i21 = (iac + m) - 1;
      for (ia = iac; ia <= i21; ia++) {
        c += A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[iy] += c;
      iy++;
    }
  }
}

/* End of code generation (xgemv.cpp) */
