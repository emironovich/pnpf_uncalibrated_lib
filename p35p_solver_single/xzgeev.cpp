/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeev.cpp
 *
 * Code generation for function 'xzgeev'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzgeev.h"
#include "xzggev.h"

/* Function Definitions */
void xzgeev(const float A[100], int *info, creal32_T alpha1[10], creal32_T
            beta1[10], creal32_T V[100])
{
  int kend;
  creal32_T At[100];
  int coltop;
  float colnorm;
  float scale;
  int k;
  float absxk;
  float t;
  for (kend = 0; kend < 100; kend++) {
    At[kend].re = A[kend];
    At[kend].im = 0.0F;
  }

  xzggev(At, info, alpha1, beta1, V);
  for (coltop = 0; coltop <= 90; coltop += 10) {
    colnorm = 0.0F;
    scale = 1.29246971E-26F;
    kend = coltop + 10;
    for (k = coltop + 1; k <= kend; k++) {
      absxk = std::abs(V[k - 1].re);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0F + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }

      absxk = std::abs(V[k - 1].im);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0F + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }
    }

    colnorm = scale * std::sqrt(colnorm);
    kend = coltop + 10;
    for (k = coltop + 1; k <= kend; k++) {
      scale = V[k - 1].re;
      absxk = V[k - 1].im;
      if (absxk == 0.0F) {
        V[k - 1].re = scale / colnorm;
        V[k - 1].im = 0.0F;
      } else if (scale == 0.0F) {
        V[k - 1].re = 0.0F;
        V[k - 1].im = absxk / colnorm;
      } else {
        V[k - 1].re = scale / colnorm;
        V[k - 1].im = absxk / colnorm;
      }
    }
  }
}

/* End of code generation (xzgeev.cpp) */
