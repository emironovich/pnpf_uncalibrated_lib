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
#include "solve_P4Pf.h"
#include "qr.h"
#include "xger.h"
#include "xgemv.h"
#include "xgeqp3.h"
#include "xnrm2.h"
#include "xgeqrf.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
void b_qr(const double A[9], double Q[9], double R[9])
{
  double b_A[9];
  double work[3];
  int i;
  int itau;
  int i_i;
  int ix0;
  double tau[3];
  double atmp;
  int ia;
  int knt;
  double xnorm;
  double beta1;
  int lastv;
  boolean_T exitg2;
  int coltop;
  int exitg1;
  memcpy(&b_A[0], &A[0], 9U * sizeof(double));
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  for (i = 0; i < 3; i++) {
    i_i = i + i * 3;
    if (1 + i < 3) {
      atmp = b_A[i_i];
      ix0 = i_i + 2;
      tau[i] = 0.0;
      xnorm = c_xnrm2(2 - i, b_A, i_i + 2);
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

          beta1 = rt_hypotd_snf(atmp, c_xnrm2(2 - i, b_A, i_i + 2));
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

      b_A[i_i] = 1.0;
      ix0 = (i + (1 + i) * 3) + 1;
      if (tau[i] != 0.0) {
        lastv = 3 - i;
        knt = (i_i - i) + 2;
        while ((lastv > 0) && (b_A[knt] == 0.0)) {
          lastv--;
          knt--;
        }

        knt = 2 - i;
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
        b_xger(lastv, knt, -tau[i], i_i + 1, work, b_A, ix0);
      }

      b_A[i_i] = atmp;
    } else {
      tau[2] = 0.0;
    }
  }

  itau = 2;
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

    work[ix0] = 0.0;
  }

  for (i = 2; i >= 0; i--) {
    i_i = (i + i * 3) + 4;
    if (i + 1 < 3) {
      b_A[i_i - 4] = 1.0;
      if (tau[itau] != 0.0) {
        lastv = 3 - i;
        knt = i_i - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0)) {
          lastv--;
          knt--;
        }

        knt = 2 - i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          coltop = i_i + (knt - 1) * 3;
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
        b_xgemv(lastv, knt, b_A, i_i, b_A, i_i - 3, work);
        b_xger(lastv, knt, -tau[itau], i_i - 3, work, b_A, i_i);
      }

      ix0 = i_i - 2;
      ia = (i_i - i) - 1;
      for (coltop = ix0; coltop <= ia; coltop++) {
        b_A[coltop - 1] *= -tau[itau];
      }
    }

    b_A[i_i - 4] = 1.0 - tau[itau];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(i_i - ix0) - 5] = 0.0;
    }

    itau--;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    Q[3 * ix0] = b_A[3 * ix0];
    knt = 1 + 3 * ix0;
    Q[knt] = b_A[knt];
    knt = 2 + 3 * ix0;
    Q[knt] = b_A[knt];
  }
}

void qr(const double A[32], double Q[64], double R[32])
{
  double tau[8];
  int lastc;
  double work[8];
  int lastv;
  boolean_T exitg2;
  int coltop;
  int ia;
  int exitg1;
  memcpy(&Q[0], &A[0], sizeof(double) << 3);
  memcpy(&Q[8], &A[8], sizeof(double) << 3);
  memcpy(&Q[16], &A[16], sizeof(double) << 3);
  memcpy(&Q[24], &A[24], sizeof(double) << 3);
  memset(&Q[32], 0, sizeof(double) << 3);
  memset(&Q[40], 0, sizeof(double) << 3);
  memset(&Q[48], 0, sizeof(double) << 3);
  memset(&Q[56], 0, sizeof(double) << 3);
  xgeqrf(Q, tau);
  for (lastc = 0; lastc < 1; lastc++) {
    R[0] = Q[0];
  }

  memset(&R[1], 0, (unsigned int)(7 * static_cast<int>(sizeof(double))));
  for (lastc = 0; lastc < 2; lastc++) {
    R[lastc + 8] = Q[lastc + 8];
  }

  memset(&R[10], 0, (unsigned int)(6 * static_cast<int>(sizeof(double))));
  for (lastc = 0; lastc < 3; lastc++) {
    R[lastc + 16] = Q[lastc + 16];
  }

  memset(&R[19], 0, (unsigned int)(5 * static_cast<int>(sizeof(double))));
  for (lastc = 0; lastc < 4; lastc++) {
    R[lastc + 24] = Q[lastc + 24];
  }

  memset(&R[28], 0, (unsigned int)(4 * static_cast<int>(sizeof(double))));
  memset(&Q[32], 0, sizeof(double) << 3);
  Q[36] = 1.0;
  memset(&Q[40], 0, sizeof(double) << 3);
  Q[45] = 1.0;
  memset(&Q[48], 0, sizeof(double) << 3);
  Q[54] = 1.0;
  memset(&Q[56], 0, sizeof(double) << 3);
  Q[63] = 1.0;
  memset(&work[0], 0, sizeof(double) << 3);
  Q[27] = 1.0;
  if (tau[3] != 0.0) {
    lastv = 5;
    lastc = 33;
    while ((lastv > 0) && (Q[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = 4;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = 36 + ((lastc - 1) << 3);
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
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
    xgemv(lastv, lastc, Q, 36, Q, 28, work);
    xger(lastv, lastc, -tau[3], 28, work, Q, 36);
  }

  for (lastc = 29; lastc < 33; lastc++) {
    Q[lastc - 1] *= -tau[3];
  }

  Q[27] = 1.0 - tau[3];
  for (lastc = 0; lastc < 3; lastc++) {
    Q[26 - lastc] = 0.0;
  }

  Q[18] = 1.0;
  if (tau[2] != 0.0) {
    lastv = 6;
    lastc = 25;
    while ((lastv > 0) && (Q[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = 5;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = 27 + ((lastc - 1) << 3);
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
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
    xgemv(lastv, lastc, Q, 27, Q, 19, work);
    xger(lastv, lastc, -tau[2], 19, work, Q, 27);
  }

  for (lastc = 20; lastc < 25; lastc++) {
    Q[lastc - 1] *= -tau[2];
  }

  Q[18] = 1.0 - tau[2];
  for (lastc = 0; lastc < 2; lastc++) {
    Q[17 - lastc] = 0.0;
  }

  Q[9] = 1.0;
  if (tau[1] != 0.0) {
    lastv = 7;
    lastc = 17;
    while ((lastv > 0) && (Q[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = 6;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = 18 + ((lastc - 1) << 3);
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
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
    xgemv(lastv, lastc, Q, 18, Q, 10, work);
    xger(lastv, lastc, -tau[1], 10, work, Q, 18);
  }

  for (lastc = 11; lastc < 17; lastc++) {
    Q[lastc - 1] *= -tau[1];
  }

  Q[9] = 1.0 - tau[1];
  for (lastc = 0; lastc < 1; lastc++) {
    Q[8] = 0.0;
  }

  Q[0] = 1.0;
  if (tau[0] != 0.0) {
    lastv = 8;
    lastc = 9;
    while ((lastv > 0) && (Q[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = 7;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = 9 + ((lastc - 1) << 3);
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
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
    xgemv(lastv, lastc, Q, 9, Q, 1, work);
    xger(lastv, lastc, -tau[0], 1, work, Q, 9);
  }

  for (lastc = 2; lastc < 9; lastc++) {
    Q[lastc - 1] *= -tau[0];
  }

  Q[0] = 1.0 - tau[0];
}

/* End of code generation (qr.cpp) */
