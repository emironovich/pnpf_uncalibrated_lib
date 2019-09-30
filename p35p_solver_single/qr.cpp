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
void qr(const float A[12], float Q[9], float R[12])
{
  int ia;
  float work[4];
  float b_A[12];
  int i;
  int ix0;
  int i_i;
  float tau[3];
  float atmp;
  int knt;
  float xnorm;
  float beta1;
  int lastv;
  int iaii;
  boolean_T exitg2;
  int coltop;
  int exitg1;
  for (ia = 0; ia < 12; ia++) {
    b_A[ia] = A[ia];
  }

  work[0] = 0.0F;
  work[1] = 0.0F;
  work[2] = 0.0F;
  work[3] = 0.0F;
  for (i = 0; i < 3; i++) {
    i_i = i + i * 3;
    if (1 + i < 3) {
      atmp = b_A[i_i];
      ix0 = i_i + 2;
      tau[i] = 0.0F;
      xnorm = c_xnrm2(2 - i, b_A, i_i + 2);
      if (xnorm != 0.0F) {
        beta1 = rt_hypotf_snf(b_A[i_i], xnorm);
        if (b_A[i_i] >= 0.0F) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 9.86076132E-32F) {
          knt = -1;
          ia = (i_i - i) + 3;
          do {
            knt++;
            for (coltop = ix0; coltop <= ia; coltop++) {
              b_A[coltop - 1] *= 1.01412048E+31F;
            }

            beta1 *= 1.01412048E+31F;
            atmp *= 1.01412048E+31F;
          } while (!(std::abs(beta1) >= 9.86076132E-32F));

          beta1 = rt_hypotf_snf(atmp, c_xnrm2(2 - i, b_A, i_i + 2));
          if (atmp >= 0.0F) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          xnorm = 1.0F / (atmp - beta1);
          for (coltop = ix0; coltop <= ia; coltop++) {
            b_A[coltop - 1] *= xnorm;
          }

          for (coltop = 0; coltop <= knt; coltop++) {
            beta1 *= 9.86076132E-32F;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - b_A[i_i]) / beta1;
          xnorm = 1.0F / (b_A[i_i] - beta1);
          ia = (i_i - i) + 3;
          for (coltop = ix0; coltop <= ia; coltop++) {
            b_A[coltop - 1] *= xnorm;
          }

          atmp = beta1;
        }
      }

      b_A[i_i] = atmp;
    } else {
      tau[2] = 0.0F;
    }

    atmp = b_A[i_i];
    b_A[i_i] = 1.0F;
    ix0 = (i + (1 + i) * 3) + 1;
    if (tau[i] != 0.0F) {
      lastv = 3 - i;
      knt = (i_i - i) + 2;
      while ((lastv > 0) && (b_A[knt] == 0.0F)) {
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
            if (b_A[ia - 1] != 0.0F) {
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
              int>(sizeof(float))));
    }
  }

  R[9] = b_A[9];
  R[10] = b_A[10];
  R[11] = b_A[11];
  i_i = 2;
  work[0] = 0.0F;
  work[1] = 0.0F;
  work[2] = 0.0F;
  work[3] = 0.0F;
  for (i = 2; i >= 0; i--) {
    iaii = (i + i * 3) + 4;
    if (i + 1 < 3) {
      b_A[iaii - 4] = 1.0F;
      if (tau[i_i] != 0.0F) {
        lastv = 3 - i;
        knt = iaii - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0F)) {
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
              if (b_A[ia - 1] != 0.0F) {
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

    b_A[iaii - 4] = 1.0F - tau[i_i];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(iaii - ix0) - 5] = 0.0F;
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
