/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xdlapy3.cpp
 *
 * Code generation for function 'xdlapy3'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xdlapy3.h"

/* Function Definitions */
double xdlapy3(double x1, double x2, double x3)
{
  double y;
  double a;
  double b;
  double c;
  a = std::abs(x1);
  b = std::abs(x2);
  c = std::abs(x3);
  if ((a > b) || rtIsNaN(b)) {
    y = a;
  } else {
    y = b;
  }

  if (c > y) {
    y = c;
  }

  if ((y > 0.0) && (!rtIsInf(y))) {
    a /= y;
    b /= y;
    c /= y;
    y *= std::sqrt((a * a + c * c) + b * b);
  } else {
    y = (a + b) + c;
  }

  return y;
}

/* End of code generation (xdlapy3.cpp) */
