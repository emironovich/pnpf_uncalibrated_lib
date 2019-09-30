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
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xgehrd.h"
#include "xzlarf.h"
#include "xnrm2.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
void xgehrd(double a[100], double tau[9])
{
  double work[10];
  int i;
  int im1n;
  int in;
  int alpha1_tmp_tmp;
  double alpha1_tmp;
  int b_i;
  double xnorm;
  double beta1;
  int iv0;
  int ic0;
  int knt;
  int lastv;
  int i18;
  int lastc;
  int rowleft;
  boolean_T exitg2;
  int ix;
  int jA;
  int i19;
  int exitg1;
  memset(&work[0], 0, 10U * sizeof(double));
  for (i = 0; i < 9; i++) {
    im1n = i * 10 + 2;
    in = (i + 1) * 10;
    alpha1_tmp_tmp = (i + 10 * i) + 1;
    alpha1_tmp = a[alpha1_tmp_tmp];
    b_i = i + 3;
    if (b_i >= 10) {
      b_i = 10;
    }

    b_i += i * 10;
    tau[i] = 0.0;
    xnorm = xnrm2(8 - i, a, b_i);
    if (xnorm != 0.0) {
      beta1 = rt_hypotd_snf(alpha1_tmp, xnorm);
      if (a[alpha1_tmp_tmp] >= 0.0) {
        beta1 = -beta1;
      }

      if (std::abs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        i18 = (b_i - i) + 7;
        do {
          knt++;
          for (rowleft = b_i; rowleft <= i18; rowleft++) {
            a[rowleft - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          alpha1_tmp *= 9.9792015476736E+291;
        } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

        beta1 = rt_hypotd_snf(alpha1_tmp, xnrm2(8 - i, a, b_i));
        if (alpha1_tmp >= 0.0) {
          beta1 = -beta1;
        }

        tau[i] = (beta1 - alpha1_tmp) / beta1;
        xnorm = 1.0 / (alpha1_tmp - beta1);
        for (rowleft = b_i; rowleft <= i18; rowleft++) {
          a[rowleft - 1] *= xnorm;
        }

        for (rowleft = 0; rowleft <= knt; rowleft++) {
          beta1 *= 1.0020841800044864E-292;
        }

        alpha1_tmp = beta1;
      } else {
        tau[i] = (beta1 - a[(i + 10 * i) + 1]) / beta1;
        xnorm = 1.0 / (a[(i + 10 * i) + 1] - beta1);
        i18 = (b_i - i) + 7;
        for (rowleft = b_i; rowleft <= i18; rowleft++) {
          a[rowleft - 1] *= xnorm;
        }

        alpha1_tmp = beta1;
      }
    }

    a[alpha1_tmp_tmp] = 1.0;
    iv0 = i + im1n;
    ic0 = in + 1;
    if (tau[i] != 0.0) {
      lastv = 8 - i;
      b_i = (iv0 - i) + 7;
      while ((lastv + 1 > 0) && (a[b_i] == 0.0)) {
        lastv--;
        b_i--;
      }

      lastc = 10;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        rowleft = in + lastc;
        jA = rowleft;
        do {
          exitg1 = 0;
          if (jA <= rowleft + lastv * 10) {
            if (a[jA - 1] != 0.0) {
              exitg1 = 1;
            } else {
              jA += 10;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = -1;
      lastc = 0;
    }

    if (lastv + 1 > 0) {
      if (lastc != 0) {
        memset(&work[0], 0, (unsigned int)(lastc * static_cast<int>(sizeof
                 (double))));
        ix = iv0;
        i18 = (in + 10 * lastv) + 1;
        for (knt = ic0; knt <= i18; knt += 10) {
          b_i = 0;
          i19 = (knt + lastc) - 1;
          for (jA = knt; jA <= i19; jA++) {
            work[b_i] += a[jA - 1] * a[ix - 1];
            b_i++;
          }

          ix++;
        }
      }

      if (!(-tau[i] == 0.0)) {
        jA = in;
        knt = iv0 - 1;
        for (b_i = 0; b_i <= lastv; b_i++) {
          if (a[knt] != 0.0) {
            xnorm = a[knt] * -tau[i];
            ix = 0;
            i18 = jA + 1;
            i19 = lastc + jA;
            for (rowleft = i18; rowleft <= i19; rowleft++) {
              a[rowleft - 1] += work[ix] * xnorm;
              ix++;
            }
          }

          knt++;
          jA += 10;
        }
      }
    }

    xzlarf(9 - i, 9 - i, i + im1n, tau[i], a, (i + in) + 2, work);
    a[alpha1_tmp_tmp] = alpha1_tmp;
  }
}

/* End of code generation (xgehrd.cpp) */
