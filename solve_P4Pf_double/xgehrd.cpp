/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgehrd.cpp
 *
 * Code generation for function 'xgehrd'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xgehrd.h"
#include "xzlarf.h"
#include "xzlarfg.h"

/* Function Definitions */
void xgehrd(creal_T a_data[], int a_size[2])
{
  int n;
  int loop_ub;
  int i17;
  creal_T work_data[12];
  int im1n;
  int in;
  int alpha1_tmp;
  creal_T alpha1;
  int u0;
  creal_T dc0;
  creal_T dc1;
  n = a_size[0];
  loop_ub = static_cast<signed char>(a_size[0]);
  for (i17 = 0; i17 < loop_ub; i17++) {
    work_data[i17].re = 0.0;
    work_data[i17].im = 0.0;
  }

  i17 = a_size[0];
  for (loop_ub = 0; loop_ub <= i17 - 2; loop_ub++) {
    im1n = loop_ub * n + 2;
    in = (loop_ub + 1) * n;
    alpha1_tmp = (loop_ub + a_size[0] * loop_ub) + 1;
    alpha1 = a_data[alpha1_tmp];
    u0 = loop_ub + 3;
    if (u0 >= n) {
      u0 = n;
    }

    dc0 = xzlarfg((n - loop_ub) - 1, &alpha1, a_data, u0 + loop_ub * n);
    a_data[alpha1_tmp].re = 1.0;
    a_data[alpha1_tmp].im = 0.0;
    xzlarf(n, (n - loop_ub) - 1, loop_ub + im1n, dc0, a_data, in + 1, n,
           work_data);
    dc1.re = dc0.re;
    dc1.im = -dc0.im;
    b_xzlarf((n - loop_ub) - 1, (n - loop_ub) - 1, loop_ub + im1n, dc1, a_data,
             (loop_ub + in) + 2, n, work_data);
    a_data[alpha1_tmp] = alpha1;
  }
}

/* End of code generation (xgehrd.cpp) */
