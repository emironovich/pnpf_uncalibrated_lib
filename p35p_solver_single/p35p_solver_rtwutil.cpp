/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * p35p_solver_rtwutil.cpp
 *
 * Code generation for function 'p35p_solver_rtwutil'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
float rt_hypotf_snf(float u0, float u1)
{
  float y;
  float a;
  float b;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0F);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0F);
  } else if (rtIsNaNF(b)) {
    y = b;
  } else {
    y = a * 1.41421354F;
  }

  return y;
}

/* End of code generation (p35p_solver_rtwutil.cpp) */
