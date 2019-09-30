/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xungorghr.cpp
 *
 * Code generation for function 'xungorghr'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xungorghr.h"
#include "xzlarf.h"

/* Function Definitions */
void xungorghr(double A[100], const double tau[9])
{
  int j;
  int ia;
  int i;
  int itau;
  double work[10];
  int i24;
  int iaii;
  int k;
  for (j = 8; j >= 0; j--) {
    ia = (j + 1) * 10;
    for (i = 0; i <= j; i++) {
      A[ia + i] = 0.0;
    }

    i24 = j + 3;
    for (i = i24; i < 11; i++) {
      k = ia + i;
      A[k - 1] = A[k - 11];
    }
  }

  memset(&A[0], 0, 10U * sizeof(double));
  A[0] = 1.0;
  itau = 8;
  memset(&work[0], 0, 10U * sizeof(double));
  for (i = 8; i >= 0; i--) {
    iaii = (i + i * 10) + 11;
    if (i + 1 < 9) {
      A[iaii] = 1.0;
      xzlarf(9 - i, 8 - i, iaii + 1, tau[itau], A, iaii + 11, work);
      ia = iaii + 2;
      i24 = (iaii - i) + 9;
      for (k = ia; k <= i24; k++) {
        A[k - 1] *= -tau[itau];
      }
    }

    A[iaii] = 1.0 - tau[itau];
    for (j = 0; j < i; j++) {
      A[(iaii - j) - 1] = 0.0;
    }

    itau--;
  }
}

/* End of code generation (xungorghr.cpp) */
