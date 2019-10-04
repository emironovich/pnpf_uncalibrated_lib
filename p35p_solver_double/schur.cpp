//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: schur.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "schur.h"
#include "p35p_solver.h"
#include "p35p_solver_rtwutil.h"
#include "rt_nonfinite.h"
#include "xdhseqr.h"
#include "xnrm2.h"
#include "xzlarf.h"
#include <cmath>
#include <cstring>

// Function Definitions

//
// Arguments    : const double A[100]
//                double V[100]
//                double T[100]
// Return Type  : void
//
void schur(const double A[100], double V[100], double T[100])
{
  boolean_T p;
  int k;
  int i;
  double work[10];
  int b_i;
  int knt;
  int im1n_tmp;
  int in;
  int alpha1_tmp;
  int ix0;
  double alpha1;
  int c_i;
  double tau[9];
  double xnorm;
  double beta1;
  int jy;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int ix;
  int exitg1;
  int i1;
  p = true;
  for (k = 0; k < 100; k++) {
    if ((!p) || (rtIsInf(A[k]) || rtIsNaN(A[k]))) {
      p = false;
    }
  }

  if (!p) {
    for (i = 0; i < 100; i++) {
      V[i] = rtNaN;
    }

    knt = 2;
    for (k = 0; k < 9; k++) {
      if (knt <= 10) {
        std::memset(&V[(k * 10 + knt) + -1], 0, (11 - knt) * sizeof(double));
      }

      knt++;
    }

    for (i = 0; i < 100; i++) {
      T[i] = rtNaN;
    }
  } else {
    std::memcpy(&T[0], &A[0], 100U * sizeof(double));
    std::memset(&work[0], 0, 10U * sizeof(double));
    for (b_i = 0; b_i < 9; b_i++) {
      im1n_tmp = b_i * 10 + 2;
      in = (b_i + 1) * 10;
      alpha1_tmp = (b_i + 10 * b_i) + 1;
      alpha1 = T[alpha1_tmp];
      if (b_i + 3 < 10) {
        c_i = b_i + 1;
      } else {
        c_i = 8;
      }

      ix0 = c_i + im1n_tmp;
      tau[b_i] = 0.0;
      xnorm = xnrm2(8 - b_i, T, ix0);
      if (xnorm != 0.0) {
        beta1 = rt_hypotd_snf(T[alpha1_tmp], xnorm);
        if (T[alpha1_tmp] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i = (ix0 - b_i) + 7;
          do {
            knt++;
            for (k = ix0; k <= i; k++) {
              T[k - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1 *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd_snf(alpha1, xnrm2(8 - b_i, T, ix0));
          if (alpha1 >= 0.0) {
            beta1 = -beta1;
          }

          tau[b_i] = (beta1 - alpha1) / beta1;
          xnorm = 1.0 / (alpha1 - beta1);
          i = (ix0 - b_i) + 7;
          for (k = ix0; k <= i; k++) {
            T[k - 1] *= xnorm;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1 = beta1;
        } else {
          tau[b_i] = (beta1 - T[alpha1_tmp]) / beta1;
          xnorm = 1.0 / (T[alpha1_tmp] - beta1);
          i = (ix0 - b_i) + 7;
          for (k = ix0; k <= i; k++) {
            T[k - 1] *= xnorm;
          }

          alpha1 = beta1;
        }
      }

      T[alpha1_tmp] = 1.0;
      jy = (b_i + im1n_tmp) - 1;
      k = in + 1;
      if (tau[b_i] != 0.0) {
        lastv = 8 - b_i;
        c_i = (jy - b_i) + 8;
        while ((lastv + 1 > 0) && (T[c_i] == 0.0)) {
          lastv--;
          c_i--;
        }

        lastc = 10;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = in + lastc;
          ix0 = knt;
          do {
            exitg1 = 0;
            if (ix0 <= knt + lastv * 10) {
              if (T[ix0 - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ix0 += 10;
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
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(double));
          }

          ix = jy;
          i = (in + 10 * lastv) + 1;
          for (knt = k; knt <= i; knt += 10) {
            c_i = 0;
            i1 = (knt + lastc) - 1;
            for (ix0 = knt; ix0 <= i1; ix0++) {
              work[c_i] += T[ix0 - 1] * T[ix];
              c_i++;
            }

            ix++;
          }
        }

        if (!(-tau[b_i] == 0.0)) {
          knt = in;
          for (k = 0; k <= lastv; k++) {
            if (T[jy] != 0.0) {
              xnorm = T[jy] * -tau[b_i];
              ix = 0;
              i = knt + 1;
              i1 = lastc + knt;
              for (c_i = i; c_i <= i1; c_i++) {
                T[c_i - 1] += work[ix] * xnorm;
                ix++;
              }
            }

            jy++;
            knt += 10;
          }
        }
      }

      xzlarf(9 - b_i, 9 - b_i, b_i + im1n_tmp, tau[b_i], T, (b_i + in) + 2, work);
      T[alpha1_tmp] = alpha1;
    }

    std::memcpy(&V[0], &T[0], 100U * sizeof(double));
    for (k = 8; k >= 0; k--) {
      ix0 = (k + 1) * 10;
      for (b_i = 0; b_i <= k; b_i++) {
        V[ix0 + b_i] = 0.0;
      }

      i = k + 3;
      for (b_i = i; b_i < 11; b_i++) {
        knt = ix0 + b_i;
        V[knt - 1] = V[knt - 11];
      }
    }

    std::memset(&V[0], 0, 10U * sizeof(double));
    V[0] = 1.0;
    knt = 8;
    std::memset(&work[0], 0, 10U * sizeof(double));
    for (b_i = 8; b_i >= 0; b_i--) {
      c_i = (b_i + b_i * 10) + 11;
      if (b_i + 1 < 9) {
        V[c_i] = 1.0;
        xzlarf(9 - b_i, 8 - b_i, c_i + 1, tau[knt], V, c_i + 11, work);
        ix0 = c_i + 2;
        i = (c_i - b_i) + 9;
        for (k = ix0; k <= i; k++) {
          V[k - 1] *= -tau[knt];
        }
      }

      V[c_i] = 1.0 - tau[knt];
      for (k = 0; k < b_i; k++) {
        V[(c_i - k) - 1] = 0.0;
      }

      knt--;
    }

    eml_dlahqr(T, V);
    knt = 4;
    for (k = 0; k < 7; k++) {
      if (knt <= 10) {
        std::memset(&T[(k * 10 + knt) + -1], 0, (11 - knt) * sizeof(double));
      }

      knt++;
    }
  }
}

//
// File trailer for schur.cpp
//
// [EOF]
//
