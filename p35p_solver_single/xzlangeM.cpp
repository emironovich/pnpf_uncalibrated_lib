/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlangeM.cpp
 *
 * Code generation for function 'xzlangeM'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzlangeM.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
float xzlangeM(const creal32_T x[100])
{
  float y;
  int k;
  boolean_T exitg1;
  float absxk;
  y = 0.0F;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 100)) {
    absxk = rt_hypotf_snf(x[k].re, x[k].im);
    if (rtIsNaNF(absxk)) {
      y = rtNaNF;
      exitg1 = true;
    } else {
      if (absxk > y) {
        y = absxk;
      }

      k++;
    }
  }

  return y;
}

/* End of code generation (xzlangeM.cpp) */
