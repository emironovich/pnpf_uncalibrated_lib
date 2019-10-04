//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzggbak.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xzggbak.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : creal_T V[100]
//                int ilo
//                int ihi
//                const int rscale[10]
// Return Type  : void
//
void xzggbak(creal_T V[100], int ilo, int ihi, const int rscale[10])
{
  int i;
  int b_i;
  int k;
  int j;
  int tmp_re_tmp;
  double tmp_re;
  double tmp_im;
  int i1;
  if (ilo > 1) {
    for (i = ilo - 2; i + 1 >= 1; i--) {
      k = rscale[i] - 1;
      if (rscale[i] != i + 1) {
        for (j = 0; j < 10; j++) {
          tmp_re_tmp = i + 10 * j;
          tmp_re = V[tmp_re_tmp].re;
          tmp_im = V[tmp_re_tmp].im;
          b_i = k + 10 * j;
          V[tmp_re_tmp] = V[b_i];
          V[b_i].re = tmp_re;
          V[b_i].im = tmp_im;
        }
      }
    }
  }

  if (ihi < 10) {
    b_i = ihi + 1;
    for (i = b_i; i < 11; i++) {
      k = rscale[i - 1];
      if (k != i) {
        for (j = 0; j < 10; j++) {
          tmp_re_tmp = (i + 10 * j) - 1;
          tmp_re = V[tmp_re_tmp].re;
          tmp_im = V[tmp_re_tmp].im;
          i1 = (k + 10 * j) - 1;
          V[tmp_re_tmp] = V[i1];
          V[i1].re = tmp_re;
          V[i1].im = tmp_im;
        }
      }
    }
  }
}

//
// File trailer for xzggbak.cpp
//
// [EOF]
//
