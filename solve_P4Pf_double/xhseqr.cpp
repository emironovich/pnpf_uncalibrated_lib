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
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xhseqr.h"
#include "xzhseqr.h"

/* Function Definitions */
int xhseqr(creal_T h_data[], int h_size[2])
{
  int info;
  int m;
  int istart;
  int jend;
  int j;
  int i;
  int i23;
  info = eml_zlahqr(h_data, h_size);
  m = h_size[0];
  if (3 < h_size[0]) {
    istart = 4;
    if (h_size[0] - 4 < h_size[1] - 1) {
      jend = h_size[0] - 3;
    } else {
      jend = h_size[1];
    }

    for (j = 0; j < jend; j++) {
      for (i = istart; i <= m; i++) {
        i23 = (i + h_size[0] * j) - 1;
        h_data[i23].re = 0.0;
        h_data[i23].im = 0.0;
      }

      istart++;
    }
  }

  return info;
}

/* End of code generation (xhseqr.cpp) */
