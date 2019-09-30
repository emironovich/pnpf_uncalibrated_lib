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
void xgehrd(float a[100], float tau[9])
{
  int i;
  float work[10];
  int im1n;
  int in;
  int alpha1_tmp_tmp;
  float alpha1_tmp;
  int b_i;
  float xnorm;
  float beta1;
  int iv0;
  int ic0;
  int knt;
  int lastv;
  int i13;
  int lastc;
  int rowleft;
  boolean_T exitg2;
  int ix;
  int jA;
  int i14;
  int exitg1;
  for (i = 0; i < 10; i++) {
    work[i] = 0.0F;
  }

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
    tau[i] = 0.0F;
    xnorm = xnrm2(8 - i, a, b_i);
    if (xnorm != 0.0F) {
      beta1 = rt_hypotf_snf(alpha1_tmp, xnorm);
      if (a[alpha1_tmp_tmp] >= 0.0F) {
        beta1 = -beta1;
      }

      if (std::abs(beta1) < 9.86076132E-32F) {
        knt = -1;
        i13 = (b_i - i) + 7;
        do {
          knt++;
          for (rowleft = b_i; rowleft <= i13; rowleft++) {
            a[rowleft - 1] *= 1.01412048E+31F;
          }

          beta1 *= 1.01412048E+31F;
          alpha1_tmp *= 1.01412048E+31F;
        } while (!(std::abs(beta1) >= 9.86076132E-32F));

        beta1 = rt_hypotf_snf(alpha1_tmp, xnrm2(8 - i, a, b_i));
        if (alpha1_tmp >= 0.0F) {
          beta1 = -beta1;
        }

        tau[i] = (beta1 - alpha1_tmp) / beta1;
        xnorm = 1.0F / (alpha1_tmp - beta1);
        for (rowleft = b_i; rowleft <= i13; rowleft++) {
          a[rowleft - 1] *= xnorm;
        }

        for (rowleft = 0; rowleft <= knt; rowleft++) {
          beta1 *= 9.86076132E-32F;
        }

        alpha1_tmp = beta1;
      } else {
        tau[i] = (beta1 - a[(i + 10 * i) + 1]) / beta1;
        xnorm = 1.0F / (a[(i + 10 * i) + 1] - beta1);
        i13 = (b_i - i) + 7;
        for (rowleft = b_i; rowleft <= i13; rowleft++) {
          a[rowleft - 1] *= xnorm;
        }

        alpha1_tmp = beta1;
      }
    }

    a[alpha1_tmp_tmp] = 1.0F;
    iv0 = i + im1n;
    ic0 = in + 1;
    if (tau[i] != 0.0F) {
      lastv = 8 - i;
      b_i = (iv0 - i) + 7;
      while ((lastv + 1 > 0) && (a[b_i] == 0.0F)) {
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
            if (a[jA - 1] != 0.0F) {
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
        memset(&work[0], 0, (unsigned int)(lastc * static_cast<int>(sizeof(float))));
        ix = iv0;
        i13 = (in + 10 * lastv) + 1;
        for (knt = ic0; knt <= i13; knt += 10) {
          b_i = 0;
          i14 = (knt + lastc) - 1;
          for (jA = knt; jA <= i14; jA++) {
            work[b_i] += a[jA - 1] * a[ix - 1];
            b_i++;
          }

          ix++;
        }
      }

      if (!(-tau[i] == 0.0F)) {
        jA = in;
        knt = iv0 - 1;
        for (b_i = 0; b_i <= lastv; b_i++) {
          if (a[knt] != 0.0F) {
            xnorm = a[knt] * -tau[i];
            ix = 0;
            i13 = jA + 1;
            i14 = lastc + jA;
            for (rowleft = i13; rowleft <= i14; rowleft++) {
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
