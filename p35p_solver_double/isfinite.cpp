/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * isfinite.cpp
 *
 * Code generation for function 'isfinite'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "isfinite.h"

/* Function Definitions */
boolean_T b_isfinite(double x)
{
  return (!rtIsInf(x)) && (!rtIsNaN(x));
}

/* End of code generation (isfinite.cpp) */
