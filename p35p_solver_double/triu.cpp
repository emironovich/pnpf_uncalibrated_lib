/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * triu.cpp
 *
 * Code generation for function 'triu'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "triu.h"

/* Function Definitions */
void triu(creal_T x[100])
{
  int istart;
  int j;
  int i;
  int i17;
  istart = 3;
  for (j = 0; j < 8; j++) {
    for (i = istart; i < 11; i++) {
      i17 = (i + 10 * j) - 1;
      x[i17].re = 0.0;
      x[i17].im = 0.0;
    }

    istart++;
  }
}

/* End of code generation (triu.cpp) */
