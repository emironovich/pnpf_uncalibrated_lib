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
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
void b_qr(const float A[9], float Q[9], float R[9])
{
  int ia;
  float work[3];
  float b_A[9];
  int i;
  int itau;
  int i_i;
  int ix0;
  float tau[3];
  float atmp;
  int knt;
  float xnorm;
  float beta1;
  int lastv;
  boolean_T exitg2;
  int coltop;
  int exitg1;
  for (ia = 0; ia < 9; ia++) {
    b_A[ia] = A[ia];
  }

  work[0] = 0.0F;
  work[1] = 0.0F;
  work[2] = 0.0F;
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

      b_A[i_i] = 1.0F;
      ix0 = (i + (1 + i) * 3) + 1;
      if (tau[i] != 0.0F) {
        lastv = 3 - i;
        knt = (i_i - i) + 2;
        while ((lastv > 0) && (b_A[knt] == 0.0F)) {
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
        b_xger(lastv, knt, -tau[i], i_i + 1, work, b_A, ix0);
      }

      b_A[i_i] = atmp;
    } else {
      tau[2] = 0.0F;
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
              int>(sizeof(float))));
    }

    work[ix0] = 0.0F;
  }

  for (i = 2; i >= 0; i--) {
    i_i = (i + i * 3) + 4;
    if (i + 1 < 3) {
      b_A[i_i - 4] = 1.0F;
      if (tau[itau] != 0.0F) {
        lastv = 3 - i;
        knt = i_i - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0F)) {
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
        b_xgemv(lastv, knt, b_A, i_i, b_A, i_i - 3, work);
        b_xger(lastv, knt, -tau[itau], i_i - 3, work, b_A, i_i);
      }

      ix0 = i_i - 2;
      ia = (i_i - i) - 1;
      for (coltop = ix0; coltop <= ia; coltop++) {
        b_A[coltop - 1] *= -tau[itau];
      }
    }

    b_A[i_i - 4] = 1.0F - tau[itau];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(i_i - ix0) - 5] = 0.0F;
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

void qr(const float A[32], float Q[64], float R[32])
{
  int lastc;
  int i;
  int b_i;
  float work[8];
  int i_i;
  float tau[8];
  float atmp;
  float xnorm;
  float beta1;
  int k;
  int knt;
  int ia;
  boolean_T exitg2;
  int exitg1;
  for (lastc = 0; lastc < 4; lastc++) {
    for (i = 0; i < 8; i++) {
      b_i = i + (lastc << 3);
      Q[b_i] = A[b_i];
    }
  }

  for (lastc = 0; lastc < 4; lastc++) {
    for (i = 0; i < 8; i++) {
      Q[i + ((lastc + 4) << 3)] = 0.0F;
    }
  }

  for (i = 0; i < 8; i++) {
    work[i] = 0.0F;
  }

  for (i = 0; i < 8; i++) {
    i_i = i + (i << 3);
    if (1 + i < 8) {
      atmp = Q[i_i];
      lastc = i_i + 2;
      tau[i] = 0.0F;
      xnorm = xnrm2(7 - i, Q, i_i + 2);
      if (xnorm != 0.0F) {
        beta1 = rt_hypotf_snf(Q[i_i], xnorm);
        if (Q[i_i] >= 0.0F) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 9.86076132E-32F) {
          knt = -1;
          b_i = (i_i - i) + 8;
          do {
            knt++;
            for (k = lastc; k <= b_i; k++) {
              Q[k - 1] *= 1.01412048E+31F;
            }

            beta1 *= 1.01412048E+31F;
            atmp *= 1.01412048E+31F;
          } while (!(std::abs(beta1) >= 9.86076132E-32F));

          beta1 = rt_hypotf_snf(atmp, xnrm2(7 - i, Q, i_i + 2));
          if (atmp >= 0.0F) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          xnorm = 1.0F / (atmp - beta1);
          for (k = lastc; k <= b_i; k++) {
            Q[k - 1] *= xnorm;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 9.86076132E-32F;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - Q[i_i]) / beta1;
          xnorm = 1.0F / (Q[i_i] - beta1);
          b_i = (i_i - i) + 8;
          for (k = lastc; k <= b_i; k++) {
            Q[k - 1] *= xnorm;
          }

          atmp = beta1;
        }
      }

      Q[i_i] = 1.0F;
      k = (i + ((1 + i) << 3)) + 1;
      if (tau[i] != 0.0F) {
        knt = 8 - i;
        b_i = (i_i - i) + 7;
        while ((knt > 0) && (Q[b_i] == 0.0F)) {
          knt--;
          b_i--;
        }

        lastc = 7 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          b_i = k + ((lastc - 1) << 3);
          ia = b_i;
          do {
            exitg1 = 0;
            if (ia <= (b_i + knt) - 1) {
              if (Q[ia - 1] != 0.0F) {
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
        knt = 0;
        lastc = 0;
      }

      if (knt > 0) {
        xgemv(knt, lastc, Q, k, Q, i_i + 1, work);
        xger(knt, lastc, -tau[i], i_i + 1, work, Q, k);
      }

      Q[i_i] = atmp;
    } else {
      tau[7] = 0.0F;
    }
  }

  for (i = 0; i < 1; i++) {
    R[0] = Q[0];
  }

  memset(&R[1], 0, (unsigned int)(7 * static_cast<int>(sizeof(float))));
  for (i = 0; i < 2; i++) {
    R[i + 8] = Q[i + 8];
  }

  memset(&R[10], 0, (unsigned int)(6 * static_cast<int>(sizeof(float))));
  for (i = 0; i < 3; i++) {
    R[i + 16] = Q[i + 16];
  }

  memset(&R[19], 0, (unsigned int)(5 * static_cast<int>(sizeof(float))));
  for (i = 0; i < 4; i++) {
    R[i + 24] = Q[i + 24];
  }

  memset(&R[28], 0, (unsigned int)(4 * static_cast<int>(sizeof(float))));
  for (lastc = 0; lastc < 4; lastc++) {
    ia = (lastc + 4) << 3;
    for (i = 0; i < 8; i++) {
      Q[ia + i] = 0.0F;
    }

    Q[(ia + lastc) + 4] = 1.0F;
  }

  for (i = 0; i < 8; i++) {
    work[i] = 0.0F;
  }

  Q[27] = 1.0F;
  if (tau[3] != 0.0F) {
    knt = 5;
    i = 33;
    while ((knt > 0) && (Q[i - 2] == 0.0F)) {
      knt--;
      i--;
    }

    lastc = 4;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      b_i = 36 + ((lastc - 1) << 3);
      ia = b_i;
      do {
        exitg1 = 0;
        if (ia <= (b_i + knt) - 1) {
          if (Q[ia - 1] != 0.0F) {
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
    knt = 0;
    lastc = 0;
  }

  if (knt > 0) {
    xgemv(knt, lastc, Q, 36, Q, 28, work);
    xger(knt, lastc, -tau[3], 28, work, Q, 36);
  }

  for (k = 29; k < 33; k++) {
    Q[k - 1] *= -tau[3];
  }

  Q[27] = 1.0F - tau[3];
  for (lastc = 0; lastc < 3; lastc++) {
    Q[26 - lastc] = 0.0F;
  }

  Q[18] = 1.0F;
  if (tau[2] != 0.0F) {
    knt = 6;
    i = 25;
    while ((knt > 0) && (Q[i - 2] == 0.0F)) {
      knt--;
      i--;
    }

    lastc = 5;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      b_i = 27 + ((lastc - 1) << 3);
      ia = b_i;
      do {
        exitg1 = 0;
        if (ia <= (b_i + knt) - 1) {
          if (Q[ia - 1] != 0.0F) {
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
    knt = 0;
    lastc = 0;
  }

  if (knt > 0) {
    xgemv(knt, lastc, Q, 27, Q, 19, work);
    xger(knt, lastc, -tau[2], 19, work, Q, 27);
  }

  for (k = 20; k < 25; k++) {
    Q[k - 1] *= -tau[2];
  }

  Q[18] = 1.0F - tau[2];
  for (lastc = 0; lastc < 2; lastc++) {
    Q[17 - lastc] = 0.0F;
  }

  Q[9] = 1.0F;
  if (tau[1] != 0.0F) {
    knt = 7;
    i = 17;
    while ((knt > 0) && (Q[i - 2] == 0.0F)) {
      knt--;
      i--;
    }

    lastc = 6;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      b_i = 18 + ((lastc - 1) << 3);
      ia = b_i;
      do {
        exitg1 = 0;
        if (ia <= (b_i + knt) - 1) {
          if (Q[ia - 1] != 0.0F) {
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
    knt = 0;
    lastc = 0;
  }

  if (knt > 0) {
    xgemv(knt, lastc, Q, 18, Q, 10, work);
    xger(knt, lastc, -tau[1], 10, work, Q, 18);
  }

  for (k = 11; k < 17; k++) {
    Q[k - 1] *= -tau[1];
  }

  Q[9] = 1.0F - tau[1];
  for (lastc = 0; lastc < 1; lastc++) {
    Q[8] = 0.0F;
  }

  Q[0] = 1.0F;
  if (tau[0] != 0.0F) {
    knt = 8;
    i = 9;
    while ((knt > 0) && (Q[i - 2] == 0.0F)) {
      knt--;
      i--;
    }

    lastc = 7;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      b_i = 9 + ((lastc - 1) << 3);
      ia = b_i;
      do {
        exitg1 = 0;
        if (ia <= (b_i + knt) - 1) {
          if (Q[ia - 1] != 0.0F) {
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
    knt = 0;
    lastc = 0;
  }

  if (knt > 0) {
    xgemv(knt, lastc, Q, 9, Q, 1, work);
    xger(knt, lastc, -tau[0], 1, work, Q, 9);
  }

  for (k = 2; k < 9; k++) {
    Q[k - 1] *= -tau[0];
  }

  Q[0] = 1.0F - tau[0];
}

/* End of code generation (qr.cpp) */
