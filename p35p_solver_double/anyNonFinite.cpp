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
#include "p35p_solver.h"
#include "anyNonFinite.h"

/* Function Definitions */
boolean_T anyNonFinite(const double x[100])
{
  boolean_T p;
  int k;
  p = true;
  for (k = 0; k < 100; k++) {
    if (p && ((!rtIsInf(x[k])) && (!rtIsNaN(x[k])))) {
      p = true;
    } else {
      p = false;
    }
  }

  return !p;
}

/* End of code generation (anyNonFinite.cpp) */
