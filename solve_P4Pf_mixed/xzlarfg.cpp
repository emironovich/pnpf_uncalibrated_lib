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
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
creal32_T xzlarfg(creal32_T *alpha1, creal32_T *x)
{
  creal32_T tau;
  float xnorm;
  float beta1;
  int knt;
  float ai;
  creal32_T b_alpha1;
  float x_re;
  float x_im;
  int k;
  tau.re = 0.0F;
  tau.im = 0.0F;
  xnorm = rt_hypotf_snf(x->re, x->im);
  if ((xnorm != 0.0F) || (alpha1->im != 0.0F)) {
    beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
    if (alpha1->re >= 0.0F) {
      beta1 = -beta1;
    }

    if (std::abs(beta1) < 9.86076132E-32F) {
      knt = -1;
      do {
        knt++;
        x->re *= 1.01412048E+31F;
        x->im *= 1.01412048E+31F;
        beta1 *= 1.01412048E+31F;
        alpha1->re *= 1.01412048E+31F;
        alpha1->im *= 1.01412048E+31F;
      } while (!(std::abs(beta1) >= 9.86076132E-32F));

      beta1 = xdlapy3(alpha1->re, alpha1->im, rt_hypotf_snf(x->re, x->im));
      if (alpha1->re >= 0.0F) {
        beta1 = -beta1;
      }

      xnorm = beta1 - alpha1->re;
      ai = 0.0F - alpha1->im;
      if (ai == 0.0F) {
        tau.re = xnorm / beta1;
        tau.im = 0.0F;
      } else if (xnorm == 0.0F) {
        tau.re = 0.0F;
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
        beta1 *= 9.86076132E-32F;
      }

      alpha1->re = beta1;
      alpha1->im = 0.0F;
    } else {
      xnorm = beta1 - alpha1->re;
      ai = 0.0F - alpha1->im;
      if (ai == 0.0F) {
        tau.re = xnorm / beta1;
        tau.im = 0.0F;
      } else if (xnorm == 0.0F) {
        tau.re = 0.0F;
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
      alpha1->im = 0.0F;
    }
  }

  return tau;
}

/* End of code generation (xzlarfg.cpp) */
