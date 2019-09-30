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
#include "ishermitian.h"
#include "anyNonFinite.h"

/* Function Declarations */
static void diagDiagUpperHessNoImag(creal_T D[100]);
static void makeD(const creal_T alpha1[10], const creal_T beta1[10], creal_T D
                  [100]);

/* Function Definitions */
static void diagDiagUpperHessNoImag(creal_T D[100])
{
  int j;
  int D_re_tmp_tmp;
  int D_re_tmp;
  int i;
  D[0].im = 0.0;
  for (j = 0; j < 9; j++) {
    D_re_tmp_tmp = 10 * (j + 1);
    D_re_tmp = (j + D_re_tmp_tmp) + 1;
    D[D_re_tmp].im = 0.0;
    D_re_tmp = (j + 10 * j) + 1;
    D[D_re_tmp].re = 0.0;
    D[D_re_tmp].im = 0.0;
    for (i = 0; i <= j; i++) {
      D_re_tmp = i + D_re_tmp_tmp;
      D[D_re_tmp].re = 0.0;
      D[D_re_tmp].im = 0.0;
    }
  }
}

static void makeD(const creal_T alpha1[10], const creal_T beta1[10], creal_T D
                  [100])
{
  int i10;
  int k;
  double brm;
  double bim;
  double d;
  for (i10 = 0; i10 < 100; i10++) {
    D[i10].re = 0.0;
    D[i10].im = 0.0;
  }

  for (k = 0; k < 10; k++) {
    if (beta1[k].im == 0.0) {
      if (alpha1[k].im == 0.0) {
        i10 = k + 10 * k;
        D[i10].re = alpha1[k].re / beta1[k].re;
        D[i10].im = 0.0;
      } else if (alpha1[k].re == 0.0) {
        D[k + 10 * k].re = 0.0;
        D[k + 10 * k].im = alpha1[k].im / beta1[k].re;
      } else {
        D[k + 10 * k].re = alpha1[k].re / beta1[k].re;
        D[k + 10 * k].im = alpha1[k].im / beta1[k].re;
      }
    } else if (beta1[k].re == 0.0) {
      if (alpha1[k].re == 0.0) {
        D[k + 10 * k].re = alpha1[k].im / beta1[k].im;
        D[k + 10 * k].im = 0.0;
      } else if (alpha1[k].im == 0.0) {
        D[k + 10 * k].re = 0.0;
        D[k + 10 * k].im = -(alpha1[k].re / beta1[k].im);
      } else {
        D[k + 10 * k].re = alpha1[k].im / beta1[k].im;
        D[k + 10 * k].im = -(alpha1[k].re / beta1[k].im);
      }
    } else {
      brm = std::abs(beta1[k].re);
      bim = std::abs(beta1[k].im);
      if (brm > bim) {
        bim = beta1[k].im / beta1[k].re;
        d = beta1[k].re + bim * beta1[k].im;
        D[k + 10 * k].re = (alpha1[k].re + bim * alpha1[k].im) / d;
        D[k + 10 * k].im = (alpha1[k].im - bim * alpha1[k].re) / d;
      } else if (bim == brm) {
        if (beta1[k].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (beta1[k].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        D[k + 10 * k].re = (alpha1[k].re * bim + alpha1[k].im * d) / brm;
        D[k + 10 * k].im = (alpha1[k].im * bim - alpha1[k].re * d) / brm;
      } else {
        bim = beta1[k].re / beta1[k].im;
        d = beta1[k].im + bim * beta1[k].re;
        D[k + 10 * k].re = (bim * alpha1[k].re + alpha1[k].im) / d;
        D[k + 10 * k].im = (bim * alpha1[k].im - alpha1[k].re) / d;
      }
    }
  }
}

void eig(const double A[100], creal_T V[100], creal_T D[100])
{
  int info;
  creal_T alpha1[10];
  creal_T beta1[10];
  int k;
  if (anyNonFinite(A)) {
    for (info = 0; info < 100; info++) {
      V[info].re = rtNaN;
      V[info].im = 0.0;
      D[info].re = 0.0;
      D[info].im = 0.0;
    }

    for (k = 0; k < 10; k++) {
      info = k + 10 * k;
      D[info].re = rtNaN;
      D[info].im = 0.0;
    }
  } else if (ishermitian(A)) {
    schur(A, V, D);
    diagDiagUpperHessNoImag(D);
  } else {
    xzgeev(A, &info, alpha1, beta1, V);
    makeD(alpha1, beta1, D);
  }
}

/* End of code generation (eig.cpp) */
