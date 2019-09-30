/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * anyNonFinite.cpp
 *
 * Code generation for function 'anyNonFinite'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "anyNonFinite.h"

/* Function Definitions */
boolean_T anyNonFinite(const creal32_T x_data[], const int x_size[2])
{
  boolean_T p;
  int nx;
  int k;
  nx = x_size[0] * x_size[1];
  p = true;
  for (k = 0; k < nx; k++) {
    if (p && ((!rtIsInfF(x_data[k].re)) && (!rtIsInfF(x_data[k].im)) &&
              ((!rtIsNaNF(x_data[k].re)) && (!rtIsNaNF(x_data[k].im))))) {
      p = true;
    } else {
      p = false;
    }
  }

  return !p;
}

/* End of code generation (anyNonFinite.cpp) */
