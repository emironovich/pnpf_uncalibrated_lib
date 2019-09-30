/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * schur.cpp
 *
 * Code generation for function 'schur'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "schur.h"
#include "triu.h"
#include "rsf2csf.h"
#include "xhseqr.h"
#include "xungorghr.h"
#include "xgehrd.h"
#include "anyNonFinite.h"

/* Function Definitions */
void schur(const double A[100], creal_T V[100], creal_T T[100])
{
  double b_A[100];
  int i7;
  double tau[9];
  double Vr[100];
  if (anyNonFinite(A)) {
    for (i7 = 0; i7 < 100; i7++) {
      V[i7].re = rtNaN;
      V[i7].im = 0.0;
    }

    triu(V);
    for (i7 = 0; i7 < 100; i7++) {
      T[i7].re = rtNaN;
      T[i7].im = 0.0;
    }
  } else {
    memcpy(&b_A[0], &A[0], 100U * sizeof(double));
    xgehrd(b_A, tau);
    memcpy(&Vr[0], &b_A[0], 100U * sizeof(double));
    xungorghr(Vr, tau);
    xhseqr(b_A, Vr);
    rsf2csf(Vr, b_A, V, T);
  }
}

/* End of code generation (schur.cpp) */
