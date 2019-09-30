/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzggbak.cpp
 *
 * Code generation for function 'xzggbak'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzggbak.h"

/* Function Definitions */
void xzggbak(creal_T V[100], int ilo, int ihi, const int rscale[10])
{
  int i;
  int i30;
  int k;
  int j;
  int tmp_re_tmp;
  double tmp_re;
  double tmp_im;
  int i31;
  if (ilo > 1) {
    for (i = ilo - 2; i + 1 >= 1; i--) {
      k = rscale[i] - 1;
      if (rscale[i] != i + 1) {
        for (j = 0; j < 10; j++) {
          tmp_re_tmp = i + 10 * j;
          tmp_re = V[tmp_re_tmp].re;
          tmp_im = V[i + 10 * j].im;
          i30 = k + 10 * j;
          V[tmp_re_tmp] = V[i30];
          V[i30].re = tmp_re;
          V[i30].im = tmp_im;
        }
      }
    }
  }

  if (ihi < 10) {
    i30 = ihi + 1;
    for (i = i30; i < 11; i++) {
      k = rscale[i - 1];
      if (k != i) {
        for (j = 0; j < 10; j++) {
          tmp_re_tmp = (i + 10 * j) - 1;
          tmp_re = V[tmp_re_tmp].re;
          tmp_im = V[tmp_re_tmp].im;
          i31 = (k + 10 * j) - 1;
          V[tmp_re_tmp] = V[i31];
          V[i31].re = tmp_re;
          V[i31].im = tmp_im;
        }
      }
    }
  }
}

/* End of code generation (xzggbak.cpp) */
