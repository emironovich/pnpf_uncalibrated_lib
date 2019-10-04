//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xnrm2.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xnrm2.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions

//
// Arguments    : int n
//                const double x[3]
// Return Type  : double
//
double b_xnrm2(int n, const double x[3])
{
  double y;
  double scale;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[1]);
    } else {
      scale = 3.3121686421112381E-170;
      absxk = std::abs(x[1]);
      if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
      }

      absxk = std::abs(x[2]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                const double x[12]
//                int ix0
// Return Type  : double
//
double c_xnrm2(int n, const double x[12], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
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

//
// Arguments    : int n
//                const double x[100]
//                int ix0
// Return Type  : double
//
double xnrm2(int n, const double x[100], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
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

//
// File trailer for xnrm2.cpp
//
// [EOF]
//
