/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sqrt.cpp
 *
 * Code generation for function 'sqrt'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "sqrt.h"
#include "xgeqp3.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
void b_sqrt(float *x)
{
  *x = std::sqrt(*x);
}

void c_sqrt(creal32_T *x)
{
  float xr;
  float xi;
  float yr;
  float absxr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0F) {
    if (xr < 0.0F) {
      yr = 0.0F;
      xr = std::sqrt(-xr);
    } else {
      yr = std::sqrt(xr);
      xr = 0.0F;
    }
  } else if (xr == 0.0F) {
    if (xi < 0.0F) {
      yr = std::sqrt(-xi / 2.0F);
      xr = -yr;
    } else {
      yr = std::sqrt(xi / 2.0F);
      xr = yr;
    }
  } else if (rtIsNaNF(xr)) {
    yr = xr;
  } else if (rtIsNaNF(xi)) {
    yr = xi;
    xr = xi;
  } else if (rtIsInfF(xi)) {
    yr = std::abs(xi);
    xr = xi;
  } else if (rtIsInfF(xr)) {
    if (xr < 0.0F) {
      yr = 0.0F;
      xr = xi * -xr;
    } else {
      yr = xr;
      xr = 0.0F;
    }
  } else {
    absxr = std::abs(xr);
    yr = std::abs(xi);
    if ((absxr > 8.50705867E+37F) || (yr > 8.50705867E+37F)) {
      absxr *= 0.5F;
      yr = rt_hypotf_snf(absxr, yr * 0.5F);
      if (yr > absxr) {
        yr = std::sqrt(yr) * std::sqrt(1.0F + absxr / yr);
      } else {
        yr = std::sqrt(yr) * 1.41421354F;
      }
    } else {
      yr = std::sqrt((rt_hypotf_snf(absxr, yr) + absxr) * 0.5F);
    }

    if (xr > 0.0F) {
      xr = 0.5F * (xi / yr);
    } else {
      if (xi < 0.0F) {
        xr = -yr;
      } else {
        xr = yr;
      }

      yr = 0.5F * (xi / xr);
    }
  }

  x->re = yr;
  x->im = xr;
}

/* End of code generation (sqrt.cpp) */
