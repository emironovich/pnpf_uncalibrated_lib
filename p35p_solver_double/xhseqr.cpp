/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xhseqr.cpp
 *
 * Code generation for function 'xhseqr'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xhseqr.h"
#include "xdhseqr.h"

/* Function Definitions */
int xhseqr(double h[100], double z[100])
{
  int info;
  int istart;
  int j;
  info = eml_dlahqr(h, z);
  istart = 4;
  for (j = 0; j < 7; j++) {
    if (istart <= 10) {
      memset(&h[(j * 10 + istart) + -1], 0, (unsigned int)((11 - istart) *
              static_cast<int>(sizeof(double))));
    }

    istart++;
  }

  return info;
}

/* End of code generation (xhseqr.cpp) */
