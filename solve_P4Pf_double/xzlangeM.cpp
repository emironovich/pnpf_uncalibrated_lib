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
#include "solve_P4Pf.h"
#include "xzlangeM.h"
#include "xgeqp3.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
double xzlangeM(const creal_T x_data[], const int x_size[2])
{
  double y;
  int i8;
  int k;
  boolean_T exitg1;
  double absxk;
  y = 0.0;
  i8 = x_size[0] * x_size[1];
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= i8 - 1)) {
    absxk = rt_hypotd_snf(x_data[k].re, x_data[k].im);
    if (rtIsNaN(absxk)) {
      y = rtNaN;
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
