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
#include "rsf2csf.h"
#include "xhseqr.h"
#include "xzlarf.h"
#include "xgehrd.h"

/* Function Definitions */
void schur(const float A[100], creal32_T V[100], creal32_T T[100])
{
  boolean_T p;
  int ia;
  float b_A[100];
  int i5;
  float tau[9];
  float Vr[100];
  int j;
  int i;
  int itau;
  float work[10];
  int Vr_tmp;
  int iaii;
  p = true;
  for (ia = 0; ia < 100; ia++) {
    if (p && ((!rtIsInfF(A[ia])) && (!rtIsNaNF(A[ia])))) {
      p = true;
    } else {
      p = false;
    }
  }

  if (!p) {
    for (i5 = 0; i5 < 100; i5++) {
      V[i5].re = rtNaNF;
      V[i5].im = 0.0F;
    }

    ia = 3;
    for (j = 0; j < 8; j++) {
      for (i = ia; i < 11; i++) {
        i5 = (i + 10 * j) - 1;
        V[i5].re = 0.0F;
        V[i5].im = 0.0F;
      }

      ia++;
    }

    for (i5 = 0; i5 < 100; i5++) {
      T[i5].re = rtNaNF;
      T[i5].im = 0.0F;
    }
  } else {
    memcpy(&b_A[0], &A[0], 100U * sizeof(float));
    xgehrd(b_A, tau);
    memcpy(&Vr[0], &b_A[0], 100U * sizeof(float));
    for (j = 8; j >= 0; j--) {
      ia = (j + 1) * 10;
      for (i = 0; i <= j; i++) {
        Vr[ia + i] = 0.0F;
      }

      i5 = j + 3;
      for (i = i5; i < 11; i++) {
        Vr_tmp = ia + i;
        Vr[Vr_tmp - 1] = Vr[Vr_tmp - 11];
      }
    }

    for (i = 0; i < 10; i++) {
      Vr[i] = 0.0F;
    }

    Vr[0] = 1.0F;
    itau = 8;
    for (i = 0; i < 10; i++) {
      work[i] = 0.0F;
    }

    for (i = 8; i >= 0; i--) {
      iaii = (i + i * 10) + 11;
      if (i + 1 < 9) {
        Vr[iaii] = 1.0F;
        xzlarf(9 - i, 8 - i, iaii + 1, tau[itau], Vr, iaii + 11, work);
        Vr_tmp = iaii + 2;
        i5 = (iaii - i) + 9;
        for (ia = Vr_tmp; ia <= i5; ia++) {
          Vr[ia - 1] *= -tau[itau];
        }
      }

      Vr[iaii] = 1.0F - tau[itau];
      for (j = 0; j < i; j++) {
        Vr[(iaii - j) - 1] = 0.0F;
      }

      itau--;
    }

    xhseqr(b_A, Vr);
    rsf2csf(Vr, b_A, V, T);
  }
}

/* End of code generation (schur.cpp) */
