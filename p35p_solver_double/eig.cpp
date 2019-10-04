//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: eig.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "eig.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include "schur.h"
#include "xzggev.h"
#include <cmath>
#include <cstring>

// Function Definitions

//
// Arguments    : const double A[100]
//                creal_T V[100]
//                creal_T D[100]
// Return Type  : void
//
void eig(const double A[100], creal_T V[100], creal_T D[100])
{
  boolean_T p;
  int k;
  int info;
  boolean_T exitg2;
  double b_V[100];
  double b_D[100];
  int exitg1;
  creal_T At[100];
  creal_T alpha1[10];
  creal_T beta1[10];
  int coltop;
  double colnorm;
  double scale;
  double t;
  double absxk;
  p = true;
  for (k = 0; k < 100; k++) {
    if ((!p) || (rtIsInf(A[k]) || rtIsNaN(A[k]))) {
      p = false;
    }
  }

  if (!p) {
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
  } else {
    p = true;
    k = 0;
    exitg2 = false;
    while ((!exitg2) && (k < 10)) {
      info = 0;
      do {
        exitg1 = 0;
        if (info <= k) {
          if (!(A[info + 10 * k] == A[k + 10 * info])) {
            p = false;
            exitg1 = 1;
          } else {
            info++;
          }
        } else {
          k++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (p) {
      schur(A, b_V, b_D);
      for (info = 0; info < 100; info++) {
        V[info].re = b_V[info];
        V[info].im = 0.0;
      }

      for (k = 0; k < 9; k++) {
        b_D[(k + 10 * k) + 1] = 0.0;
        for (info = 0; info <= k; info++) {
          b_D[info + 10 * (k + 1)] = 0.0;
        }
      }

      for (info = 0; info < 100; info++) {
        D[info].re = b_D[info];
        D[info].im = 0.0;
      }
    } else {
      for (info = 0; info < 100; info++) {
        At[info].re = A[info];
        At[info].im = 0.0;
      }

      xzggev(At, &info, alpha1, beta1, V);
      for (coltop = 0; coltop <= 90; coltop += 10) {
        colnorm = 0.0;
        scale = 3.3121686421112381E-170;
        info = coltop + 10;
        for (k = coltop + 1; k <= info; k++) {
          absxk = std::abs(V[k - 1].re);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }

          absxk = std::abs(V[k - 1].im);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }
        }

        colnorm = scale * std::sqrt(colnorm);
        info = coltop + 10;
        for (k = coltop + 1; k <= info; k++) {
          absxk = V[k - 1].re;
          scale = V[k - 1].im;
          if (scale == 0.0) {
            absxk /= colnorm;
            scale = 0.0;
          } else if (absxk == 0.0) {
            absxk = 0.0;
            scale /= colnorm;
          } else {
            absxk /= colnorm;
            scale /= colnorm;
          }

          V[k - 1].re = absxk;
          V[k - 1].im = scale;
        }
      }

      std::memset(&D[0], 0, 100U * sizeof(creal_T));
      for (k = 0; k < 10; k++) {
        if (beta1[k].im == 0.0) {
          if (alpha1[k].im == 0.0) {
            info = k + 10 * k;
            D[info].re = alpha1[k].re / beta1[k].re;
            D[info].im = 0.0;
          } else if (alpha1[k].re == 0.0) {
            info = k + 10 * k;
            D[info].re = 0.0;
            D[info].im = alpha1[k].im / beta1[k].re;
          } else {
            info = k + 10 * k;
            D[info].re = alpha1[k].re / beta1[k].re;
            D[info].im = alpha1[k].im / beta1[k].re;
          }
        } else if (beta1[k].re == 0.0) {
          if (alpha1[k].re == 0.0) {
            info = k + 10 * k;
            D[info].re = alpha1[k].im / beta1[k].im;
            D[info].im = 0.0;
          } else if (alpha1[k].im == 0.0) {
            info = k + 10 * k;
            D[info].re = 0.0;
            D[info].im = -(alpha1[k].re / beta1[k].im);
          } else {
            info = k + 10 * k;
            D[info].re = alpha1[k].im / beta1[k].im;
            D[info].im = -(alpha1[k].re / beta1[k].im);
          }
        } else {
          t = std::abs(beta1[k].re);
          scale = std::abs(beta1[k].im);
          if (t > scale) {
            scale = beta1[k].im / beta1[k].re;
            absxk = beta1[k].re + scale * beta1[k].im;
            info = k + 10 * k;
            D[info].re = (alpha1[k].re + scale * alpha1[k].im) / absxk;
            D[info].im = (alpha1[k].im - scale * alpha1[k].re) / absxk;
          } else if (scale == t) {
            if (beta1[k].re > 0.0) {
              scale = 0.5;
            } else {
              scale = -0.5;
            }

            if (beta1[k].im > 0.0) {
              absxk = 0.5;
            } else {
              absxk = -0.5;
            }

            info = k + 10 * k;
            D[info].re = (alpha1[k].re * scale + alpha1[k].im * absxk) / t;
            D[info].im = (alpha1[k].im * scale - alpha1[k].re * absxk) / t;
          } else {
            scale = beta1[k].re / beta1[k].im;
            absxk = beta1[k].im + scale * beta1[k].re;
            info = k + 10 * k;
            D[info].re = (scale * alpha1[k].re + alpha1[k].im) / absxk;
            D[info].im = (scale * alpha1[k].im - alpha1[k].re) / absxk;
          }
        }
      }
    }
  }
}

//
// File trailer for eig.cpp
//
// [EOF]
//
