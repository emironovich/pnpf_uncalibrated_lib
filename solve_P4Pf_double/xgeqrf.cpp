/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgeqrf.cpp
 *
 * Code generation for function 'xgeqrf'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xgeqrf.h"
#include "xger.h"
#include "xgemv.h"
#include "xgeqp3.h"
#include "xnrm2.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
void xgeqrf(double A[64], double tau[8])
{
  double work[8];
  int i;
  int i_i;
  double atmp;
  int lastc;
  double xnorm;
  double beta1;
  int knt;
  int coltop;
  int lastv;
  int k;
  boolean_T exitg2;
  int exitg1;
  memset(&work[0], 0, sizeof(double) << 3);
  for (i = 0; i < 8; i++) {
    i_i = i + (i << 3);
    if (1 + i < 8) {
      atmp = A[i_i];
      lastc = i_i + 2;
      tau[i] = 0.0;
      xnorm = xnrm2(7 - i, A, i_i + 2);
      if (xnorm != 0.0) {
        beta1 = rt_hypotd_snf(A[i_i], xnorm);
        if (A[i_i] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          coltop = (i_i - i) + 8;
          do {
            knt++;
            for (k = lastc; k <= coltop; k++) {
              A[k - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd_snf(atmp, xnrm2(7 - i, A, i_i + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          xnorm = 1.0 / (atmp - beta1);
          for (k = lastc; k <= coltop; k++) {
            A[k - 1] *= xnorm;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - A[i_i]) / beta1;
          xnorm = 1.0 / (A[i_i] - beta1);
          coltop = (i_i - i) + 8;
          for (k = lastc; k <= coltop; k++) {
            A[k - 1] *= xnorm;
          }

          atmp = beta1;
        }
      }

      A[i_i] = atmp;
      atmp = A[i_i];
      A[i_i] = 1.0;
      knt = (i + ((1 + i) << 3)) + 1;
      if (tau[i] != 0.0) {
        lastv = 8 - i;
        lastc = (i_i - i) + 7;
        while ((lastv > 0) && (A[lastc] == 0.0)) {
          lastv--;
          lastc--;
        }

        lastc = 7 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          coltop = knt + ((lastc - 1) << 3);
          k = coltop;
          do {
            exitg1 = 0;
            if (k <= (coltop + lastv) - 1) {
              if (A[k - 1] != 0.0) {
                exitg1 = 1;
              } else {
                k++;
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
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        xgemv(lastv, lastc, A, knt, A, i_i + 1, work);
        xger(lastv, lastc, -tau[i], i_i + 1, work, A, knt);
      }

      A[i_i] = atmp;
    } else {
      tau[7] = 0.0;
    }
  }
}

/* End of code generation (xgeqrf.cpp) */
