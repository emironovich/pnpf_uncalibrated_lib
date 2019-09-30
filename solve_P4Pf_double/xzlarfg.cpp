/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlarfg.cpp
 *
 * Code generation for function 'xzlarfg'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzlarfg.h"
#include "recip.h"
#include "xdlapy3.h"
#include "xgeqp3.h"
#include "xnrm2.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
creal_T b_xzlarfg(creal_T *alpha1, creal_T *x)
{
  creal_T tau;
  double xnorm;
  double beta1;
  int knt;
  double ai;
  creal_T b_alpha1;
  double x_re;
  double x_im;
  int k;
  tau.re = 0.0;
  tau.im = 0.0;
  xnorm = rt_hypotd_snf(x->re, x->im);
  if ((xnorm != 0.0) || (alpha1->im != 0.0)) {
    beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
    if (alpha1->re >= 0.0) {
      beta1 = -beta1;
    }

    if (std::abs(beta1) < 1.0020841800044864E-292) {
      knt = -1;
      do {
        knt++;
        x->re *= 9.9792015476736E+291;
        x->im *= 9.9792015476736E+291;
        beta1 *= 9.9792015476736E+291;
        alpha1->re *= 9.9792015476736E+291;
        alpha1->im *= 9.9792015476736E+291;
      } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

      beta1 = xdlapy3(alpha1->re, alpha1->im, rt_hypotd_snf(x->re, x->im));
      if (alpha1->re >= 0.0) {
        beta1 = -beta1;
      }

      xnorm = beta1 - alpha1->re;
      ai = 0.0 - alpha1->im;
      if (ai == 0.0) {
        tau.re = xnorm / beta1;
        tau.im = 0.0;
      } else if (xnorm == 0.0) {
        tau.re = 0.0;
        tau.im = ai / beta1;
      } else {
        tau.re = xnorm / beta1;
        tau.im = ai / beta1;
      }

      b_alpha1.re = alpha1->re - beta1;
      b_alpha1.im = alpha1->im;
      *alpha1 = recip(b_alpha1);
      xnorm = alpha1->re;
      ai = alpha1->im;
      x_re = x->re;
      x_im = x->im;
      x->re = xnorm * x_re - ai * x_im;
      x->im = xnorm * x_im + ai * x_re;
      for (k = 0; k <= knt; k++) {
        beta1 *= 1.0020841800044864E-292;
      }

      alpha1->re = beta1;
      alpha1->im = 0.0;
    } else {
      xnorm = beta1 - alpha1->re;
      ai = 0.0 - alpha1->im;
      if (ai == 0.0) {
        tau.re = xnorm / beta1;
        tau.im = 0.0;
      } else if (xnorm == 0.0) {
        tau.re = 0.0;
        tau.im = ai / beta1;
      } else {
        tau.re = xnorm / beta1;
        tau.im = ai / beta1;
      }

      b_alpha1.re = alpha1->re - beta1;
      b_alpha1.im = alpha1->im;
      *alpha1 = recip(b_alpha1);
      xnorm = alpha1->re;
      ai = alpha1->im;
      x_re = x->re;
      x_im = x->im;
      x->re = xnorm * x_re - ai * x_im;
      x->im = xnorm * x_im + ai * x_re;
      alpha1->re = beta1;
      alpha1->im = 0.0;
    }
  }

  return tau;
}

creal_T xzlarfg(int n, creal_T *alpha1, creal_T x_data[], int ix0)
{
  creal_T tau;
  double xnorm;
  double beta1;
  int knt;
  double ai;
  int i18;
  int k;
  creal_T b_alpha1;
  double x_data_re;
  double x_data_im;
  tau.re = 0.0;
  tau.im = 0.0;
  if (n > 0) {
    xnorm = b_xnrm2(n - 1, x_data, ix0);
    if ((xnorm != 0.0) || (alpha1->im != 0.0)) {
      beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
      if (alpha1->re >= 0.0) {
        beta1 = -beta1;
      }

      if (std::abs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        i18 = (ix0 + n) - 2;
        do {
          knt++;
          for (k = ix0; k <= i18; k++) {
            x_data_re = x_data[k - 1].re;
            x_data_im = x_data[k - 1].im;
            x_data[k - 1].re = 9.9792015476736E+291 * x_data_re - 0.0 *
              x_data_im;
            x_data[k - 1].im = 9.9792015476736E+291 * x_data_im + 0.0 *
              x_data_re;
          }

          beta1 *= 9.9792015476736E+291;
          alpha1->re *= 9.9792015476736E+291;
          alpha1->im *= 9.9792015476736E+291;
        } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

        beta1 = xdlapy3(alpha1->re, alpha1->im, b_xnrm2(n - 1, x_data, ix0));
        if (alpha1->re >= 0.0) {
          beta1 = -beta1;
        }

        xnorm = beta1 - alpha1->re;
        ai = 0.0 - alpha1->im;
        if (ai == 0.0) {
          tau.re = xnorm / beta1;
          tau.im = 0.0;
        } else if (xnorm == 0.0) {
          tau.re = 0.0;
          tau.im = ai / beta1;
        } else {
          tau.re = xnorm / beta1;
          tau.im = ai / beta1;
        }

        b_alpha1.re = alpha1->re - beta1;
        b_alpha1.im = alpha1->im;
        *alpha1 = recip(b_alpha1);
        for (k = ix0; k <= i18; k++) {
          xnorm = alpha1->re;
          ai = alpha1->im;
          x_data_re = x_data[k - 1].re;
          x_data_im = x_data[k - 1].im;
          x_data[k - 1].re = xnorm * x_data_re - ai * x_data_im;
          x_data[k - 1].im = xnorm * x_data_im + ai * x_data_re;
        }

        for (k = 0; k <= knt; k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        alpha1->re = beta1;
        alpha1->im = 0.0;
      } else {
        xnorm = beta1 - alpha1->re;
        ai = 0.0 - alpha1->im;
        if (ai == 0.0) {
          tau.re = xnorm / beta1;
          tau.im = 0.0;
        } else if (xnorm == 0.0) {
          tau.re = 0.0;
          tau.im = ai / beta1;
        } else {
          tau.re = xnorm / beta1;
          tau.im = ai / beta1;
        }

        b_alpha1.re = alpha1->re - beta1;
        b_alpha1.im = alpha1->im;
        *alpha1 = recip(b_alpha1);
        i18 = (ix0 + n) - 2;
        for (k = ix0; k <= i18; k++) {
          xnorm = alpha1->re;
          ai = alpha1->im;
          x_data_re = x_data[k - 1].re;
          x_data_im = x_data[k - 1].im;
          x_data[k - 1].re = xnorm * x_data_re - ai * x_data_im;
          x_data[k - 1].im = xnorm * x_data_im + ai * x_data_re;
        }

        alpha1->re = beta1;
        alpha1->im = 0.0;
      }
    }
  }

  return tau;
}

/* End of code generation (xzlarfg.cpp) */
