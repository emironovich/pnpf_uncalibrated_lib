/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eig.cpp
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "eig.h"
#include "schur.h"
#include "xzgeev.h"

/* Function Definitions */
void eig(const float A[100], creal32_T V[100], creal32_T D[100])
{
  boolean_T p;
  int info;
  int i4;
  int j;
  boolean_T exitg2;
  creal32_T alpha1[10];
  creal32_T beta1[10];
  int exitg1;
  int D_re_tmp_tmp;
  float brm;
  float bim;
  float d;
  p = true;
  for (info = 0; info < 100; info++) {
    if (p && ((!rtIsInfF(A[info])) && (!rtIsNaNF(A[info])))) {
      p = true;
    } else {
      p = false;
    }
  }

  if (!p) {
    for (i4 = 0; i4 < 100; i4++) {
      V[i4].re = rtNaNF;
      V[i4].im = 0.0F;
      D[i4].re = 0.0F;
      D[i4].im = 0.0F;
    }

    for (info = 0; info < 10; info++) {
      i4 = info + 10 * info;
      D[i4].re = rtNaNF;
      D[i4].im = 0.0F;
    }
  } else {
    p = true;
    j = 0;
    exitg2 = false;
    while ((!exitg2) && (j < 10)) {
      info = 0;
      do {
        exitg1 = 0;
        if (info <= j) {
          if (!(A[info + 10 * j] == A[j + 10 * info])) {
            p = false;
            exitg1 = 1;
          } else {
            info++;
          }
        } else {
          j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (p) {
      schur(A, V, D);
      D[0].im = 0.0F;
      for (j = 0; j < 9; j++) {
        D_re_tmp_tmp = 10 * (j + 1);
        info = (j + D_re_tmp_tmp) + 1;
        D[info].im = 0.0F;
        i4 = (j + 10 * j) + 1;
        D[i4].re = 0.0F;
        D[i4].im = 0.0F;
        for (info = 0; info <= j; info++) {
          i4 = info + D_re_tmp_tmp;
          D[i4].re = 0.0F;
          D[i4].im = 0.0F;
        }
      }
    } else {
      xzgeev(A, &info, alpha1, beta1, V);
      for (i4 = 0; i4 < 100; i4++) {
        D[i4].re = 0.0F;
        D[i4].im = 0.0F;
      }

      for (info = 0; info < 10; info++) {
        if (beta1[info].im == 0.0F) {
          if (alpha1[info].im == 0.0F) {
            i4 = info + 10 * info;
            D[i4].re = alpha1[info].re / beta1[info].re;
            D[i4].im = 0.0F;
          } else if (alpha1[info].re == 0.0F) {
            D[info + 10 * info].re = 0.0F;
            D[info + 10 * info].im = alpha1[info].im / beta1[info].re;
          } else {
            D[info + 10 * info].re = alpha1[info].re / beta1[info].re;
            D[info + 10 * info].im = alpha1[info].im / beta1[info].re;
          }
        } else if (beta1[info].re == 0.0F) {
          if (alpha1[info].re == 0.0F) {
            D[info + 10 * info].re = alpha1[info].im / beta1[info].im;
            D[info + 10 * info].im = 0.0F;
          } else if (alpha1[info].im == 0.0F) {
            D[info + 10 * info].re = 0.0F;
            D[info + 10 * info].im = -(alpha1[info].re / beta1[info].im);
          } else {
            D[info + 10 * info].re = alpha1[info].im / beta1[info].im;
            D[info + 10 * info].im = -(alpha1[info].re / beta1[info].im);
          }
        } else {
          brm = std::abs(beta1[info].re);
          bim = std::abs(beta1[info].im);
          if (brm > bim) {
            bim = beta1[info].im / beta1[info].re;
            d = beta1[info].re + bim * beta1[info].im;
            D[info + 10 * info].re = (alpha1[info].re + bim * alpha1[info].im) /
              d;
            D[info + 10 * info].im = (alpha1[info].im - bim * alpha1[info].re) /
              d;
          } else if (bim == brm) {
            if (beta1[info].re > 0.0F) {
              bim = 0.5F;
            } else {
              bim = -0.5F;
            }

            if (beta1[info].im > 0.0F) {
              d = 0.5F;
            } else {
              d = -0.5F;
            }

            D[info + 10 * info].re = (alpha1[info].re * bim + alpha1[info].im *
              d) / brm;
            D[info + 10 * info].im = (alpha1[info].im * bim - alpha1[info].re *
              d) / brm;
          } else {
            bim = beta1[info].re / beta1[info].im;
            d = beta1[info].im + bim * beta1[info].re;
            D[info + 10 * info].re = (bim * alpha1[info].re + alpha1[info].im) /
              d;
            D[info + 10 * info].im = (bim * alpha1[info].im - alpha1[info].re) /
              d;
          }
        }
      }
    }
  }
}

/* End of code generation (eig.cpp) */
