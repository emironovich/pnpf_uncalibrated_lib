/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xnrm2.cpp
 *
 * Code generation for function 'xnrm2'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xnrm2.h"
#include "xgeqp3.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
float b_xnrm2(int n, const creal32_T x_data[], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = rt_hypotf_snf(x_data[ix0 - 1].re, x_data[ix0 - 1].im);
    } else {
      scale = 1.29246971E-26F;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x_data[k - 1].re);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0F + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }

        absxk = std::abs(x_data[k - 1].im);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0F + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

float c_xnrm2(int n, const float x[9], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 1.29246971E-26F;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0F + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

float d_xnrm2(int n, const float x[36], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  scale = 1.29246971E-26F;
  kend = (ix0 + n) - 1;
  for (k = ix0; k <= kend; k++) {
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}

float e_xnrm2(int n, const float x[4], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  scale = 1.29246971E-26F;
  kend = (ix0 + n) - 1;
  for (k = ix0; k <= kend; k++) {
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}

float xnrm2(int n, const float x[64], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 1.29246971E-26F;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0F + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

/* End of code generation (xnrm2.cpp) */
