/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * qr.cpp
 *
 * Code generation for function 'qr'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "qr.h"
#include "xger.h"
#include "xgemv.h"
#include "xnrm2.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
void qr(const double A[12], double Q[9], double R[12])
{
  double b_A[12];
  double work[4];
  int i;
  int ix0;
  int i_i;
  double tau[3];
  double atmp;
  int ia;
  int knt;
  double xnorm;
  double beta1;
  int lastv;
  int iaii;
  boolean_T exitg2;
  int coltop;
  int exitg1;
  memcpy(&b_A[0], &A[0], 12U * sizeof(double));
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 0; i < 3; i++) {
    i_i = i + i * 3;
    if (1 + i < 3) {
      atmp = b_A[i_i];
      ix0 = i_i + 2;
      tau[i] = 0.0;
      xnorm = d_xnrm2(2 - i, b_A, i_i + 2);
      if (xnorm != 0.0) {
        beta1 = rt_hypotd_snf(b_A[i_i], xnorm);
        if (b_A[i_i] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          ia = (i_i - i) + 3;
          do {
            knt++;
            for (coltop = ix0; coltop <= ia; coltop++) {
              b_A[coltop - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd_snf(atmp, d_xnrm2(2 - i, b_A, i_i + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          xnorm = 1.0 / (atmp - beta1);
          for (coltop = ix0; coltop <= ia; coltop++) {
            b_A[coltop - 1] *= xnorm;
          }

          for (coltop = 0; coltop <= knt; coltop++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - b_A[i_i]) / beta1;
          xnorm = 1.0 / (b_A[i_i] - beta1);
          ia = (i_i - i) + 3;
          for (coltop = ix0; coltop <= ia; coltop++) {
            b_A[coltop - 1] *= xnorm;
          }

          atmp = beta1;
        }
      }

      b_A[i_i] = atmp;
    } else {
      tau[2] = 0.0;
    }

    atmp = b_A[i_i];
    b_A[i_i] = 1.0;
    ix0 = (i + (1 + i) * 3) + 1;
    if (tau[i] != 0.0) {
      lastv = 3 - i;
      knt = (i_i - i) + 2;
      while ((lastv > 0) && (b_A[knt] == 0.0)) {
        lastv--;
        knt--;
      }

      knt = 3 - i;
      exitg2 = false;
      while ((!exitg2) && (knt > 0)) {
        coltop = ix0 + (knt - 1) * 3;
        ia = coltop;
        do {
          exitg1 = 0;
          if (ia <= (coltop + lastv) - 1) {
            if (b_A[ia - 1] != 0.0) {
              exitg1 = 1;
            } else {
              ia++;
            }
          } else {
            knt--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      knt = 0;
    }

    if (lastv > 0) {
      b_xgemv(lastv, knt, b_A, ix0, b_A, i_i + 1, work);
      xger(lastv, knt, -tau[i], i_i + 1, work, b_A, ix0);
    }

    b_A[i_i] = atmp;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    for (i = 0; i <= ix0; i++) {
      knt = i + 3 * ix0;
      R[knt] = b_A[knt];
    }

    ia = ix0 + 2;
    if (ia <= 3) {
      memset(&R[(ix0 * 3 + ia) + -1], 0, (unsigned int)((4 - ia) * static_cast<
              int>(sizeof(double))));
    }
  }

  R[9] = b_A[9];
  R[10] = b_A[10];
  R[11] = b_A[11];
  i_i = 2;
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 2; i >= 0; i--) {
    iaii = (i + i * 3) + 4;
    if (i + 1 < 3) {
      b_A[iaii - 4] = 1.0;
      if (tau[i_i] != 0.0) {
        lastv = 3 - i;
        knt = iaii - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0)) {
          lastv--;
          knt--;
        }

        knt = 2 - i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          coltop = iaii + (knt - 1) * 3;
          ia = coltop;
          do {
            exitg1 = 0;
            if (ia <= (coltop + lastv) - 1) {
              if (b_A[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        knt = 0;
      }

      if (lastv > 0) {
        b_xgemv(lastv, knt, b_A, iaii, b_A, iaii - 3, work);
        xger(lastv, knt, -tau[i_i], iaii - 3, work, b_A, iaii);
      }

      ix0 = iaii - 2;
      ia = (iaii - i) - 1;
      for (coltop = ix0; coltop <= ia; coltop++) {
        b_A[coltop - 1] *= -tau[i_i];
      }
    }

    b_A[iaii - 4] = 1.0 - tau[i_i];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(iaii - ix0) - 5] = 0.0;
    }

    i_i--;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    Q[3 * ix0] = b_A[3 * ix0];
    knt = 1 + 3 * ix0;
    Q[knt] = b_A[knt];
    knt = 2 + 3 * ix0;
    Q[knt] = b_A[knt];
  }
}

/* End of code generation (qr.cpp) */
