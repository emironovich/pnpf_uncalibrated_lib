/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rsf2csf.cpp
 *
 * Code generation for function 'rsf2csf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "rsf2csf.h"
#include "xdlanv2.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
void rsf2csf(const double Ur[100], const double Tr[100], creal_T U[100], creal_T
             T[100])
{
  int b_tmp;
  int m;
  int mm1;
  int i8;
  double r;
  int b_tmp_tmp;
  double b;
  double s;
  double t1_re;
  double t1_im;
  double rt1r;
  double rt1i;
  double rt2r;
  double mu1_im;
  double mu1_re;
  double sn;
  int re_tmp;
  int b_re_tmp;
  for (b_tmp = 0; b_tmp < 100; b_tmp++) {
    T[b_tmp].re = Tr[b_tmp];
    T[b_tmp].im = 0.0;
    U[b_tmp].re = Ur[b_tmp];
    U[b_tmp].im = 0.0;
  }

  for (m = 8; m >= 0; m--) {
    mm1 = m + 1;
    b_tmp = m + 10 * m;
    i8 = b_tmp + 1;
    if (Tr[i8] != 0.0) {
      r = Tr[b_tmp];
      b_tmp_tmp = 10 * (m + 1);
      b_tmp = m + b_tmp_tmp;
      b = Tr[b_tmp];
      s = Tr[i8];
      t1_re = Tr[b_tmp + 1];
      t1_im = t1_re;
      xdlanv2(&r, &b, &s, &t1_im, &rt1r, &rt1i, &rt2r, &mu1_im, &mu1_re, &sn);
      mu1_re = rt1r - t1_re;
      r = rt_hypotd_snf(rt_hypotd_snf(mu1_re, rt1i), Tr[i8]);
      if (rt1i == 0.0) {
        mu1_re /= r;
        mu1_im = 0.0;
      } else if (mu1_re == 0.0) {
        mu1_re = 0.0;
        mu1_im = rt1i / r;
      } else {
        mu1_re /= r;
        mu1_im = rt1i / r;
      }

      s = Tr[i8] / r;
      for (b_tmp = mm1; b_tmp < 11; b_tmp++) {
        re_tmp = m + 10 * (b_tmp - 1);
        t1_im = T[m + 10 * (b_tmp - 1)].im;
        t1_re = T[re_tmp].re;
        b_re_tmp = re_tmp + 1;
        rt1r = T[b_re_tmp].re;
        rt2r = T[(m + 10 * (b_tmp - 1)) + 1].im;
        r = T[re_tmp].re;
        T[re_tmp].re = (mu1_re * T[re_tmp].re + mu1_im * t1_im) + s * T[b_re_tmp]
          .re;
        T[re_tmp].im = (mu1_re * t1_im - mu1_im * r) + s * rt2r;
        r = mu1_re * rt1r - mu1_im * rt2r;
        b = mu1_re * rt2r + mu1_im * rt1r;
        T[b_re_tmp].re = r - s * t1_re;
        T[b_re_tmp].im = b - s * t1_im;
      }

      for (b_tmp = 0; b_tmp <= m + 1; b_tmp++) {
        re_tmp = b_tmp + 10 * m;
        t1_im = T[b_tmp + 10 * m].im;
        t1_re = T[re_tmp].re;
        b_re_tmp = b_tmp + b_tmp_tmp;
        rt1r = T[b_re_tmp].re;
        rt2r = T[b_tmp + 10 * (m + 1)].im;
        r = mu1_re * T[re_tmp].re - mu1_im * t1_im;
        b = mu1_re * t1_im + mu1_im * T[re_tmp].re;
        T[re_tmp].re = r + s * T[b_re_tmp].re;
        T[re_tmp].im = b + s * rt2r;
        T[b_re_tmp].re = (mu1_re * rt1r + mu1_im * rt2r) - s * t1_re;
        T[b_re_tmp].im = (mu1_re * rt2r - mu1_im * rt1r) - s * t1_im;
      }

      for (b_tmp = 0; b_tmp < 10; b_tmp++) {
        re_tmp = b_tmp + 10 * m;
        t1_im = U[b_tmp + 10 * m].im;
        t1_re = U[re_tmp].re;
        b_re_tmp = b_tmp + b_tmp_tmp;
        rt1r = U[b_re_tmp].re;
        rt2r = U[b_tmp + 10 * (m + 1)].im;
        r = mu1_re * U[re_tmp].re - mu1_im * t1_im;
        b = mu1_re * t1_im + mu1_im * U[re_tmp].re;
        U[re_tmp].re = r + s * U[b_re_tmp].re;
        U[re_tmp].im = b + s * rt2r;
        U[b_re_tmp].re = (mu1_re * rt1r + mu1_im * rt2r) - s * t1_re;
        U[b_re_tmp].im = (mu1_re * rt2r - mu1_im * rt1r) - s * t1_im;
      }

      T[i8].re = 0.0;
      T[i8].im = 0.0;
    }
  }
}

/* End of code generation (rsf2csf.cpp) */
