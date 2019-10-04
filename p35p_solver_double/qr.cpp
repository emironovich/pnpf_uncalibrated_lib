//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: qr.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "qr.h"
#include "p35p_solver.h"
#include "p35p_solver_rtwutil.h"
#include "rt_nonfinite.h"
#include "xgerc.h"
#include "xnrm2.h"
#include <cmath>
#include <cstring>

// Function Definitions

//
// Arguments    : const double A[12]
//                double Q[9]
//                double R[12]
// Return Type  : void
//
void qr(const double A[12], double Q[9], double R[12])
{
  double b_A[12];
  double tau[3];
  double work[4];
  int i;
  int ix0;
  int ii;
  double atmp;
  int b_i;
  int knt;
  double c;
  int lastv;
  double beta1;
  int lastc;
  int iaii;
  boolean_T exitg2;
  int k;
  int ia;
  int exitg1;
  int ix;
  int i1;
  std::memcpy(&b_A[0], &A[0], 12U * sizeof(double));
  tau[0] = 0.0;
  tau[1] = 0.0;
  tau[2] = 0.0;
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 0; i < 3; i++) {
    ii = i * 3 + i;
    if (i + 1 < 3) {
      atmp = b_A[ii];
      ix0 = ii + 2;
      tau[i] = 0.0;
      c = c_xnrm2(2 - i, b_A, ii + 2);
      if (c != 0.0) {
        beta1 = rt_hypotd_snf(b_A[ii], c);
        if (b_A[ii] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          b_i = (ii - i) + 3;
          do {
            knt++;
            for (k = ix0; k <= b_i; k++) {
              b_A[k - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd_snf(atmp, c_xnrm2(2 - i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          c = 1.0 / (atmp - beta1);
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - b_A[ii]) / beta1;
          c = 1.0 / (b_A[ii] - beta1);
          b_i = (ii - i) + 3;
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = atmp;
    } else {
      tau[2] = 0.0;
    }

    atmp = b_A[ii];
    b_A[ii] = 1.0;
    if (tau[i] != 0.0) {
      lastv = 3 - i;
      knt = (ii - i) + 2;
      while ((lastv > 0) && (b_A[knt] == 0.0)) {
        lastv--;
        knt--;
      }

      lastc = 3 - i;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        knt = (ii + (lastc - 1) * 3) + 3;
        ia = knt;
        do {
          exitg1 = 0;
          if (ia + 1 <= knt + lastv) {
            if (b_A[ia] != 0.0) {
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
      knt = ii + 4;
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work[0], 0, lastc * sizeof(double));
        }

        ix0 = 0;
        b_i = (ii + 3 * (lastc - 1)) + 4;
        for (k = knt; k <= b_i; k += 3) {
          ix = ii;
          c = 0.0;
          i1 = (k + lastv) - 1;
          for (ia = k; ia <= i1; ia++) {
            c += b_A[ia - 1] * b_A[ix];
            ix++;
          }

          work[ix0] += c;
          ix0++;
        }
      }

      xgerc(lastv, lastc, -tau[i], ii + 1, work, b_A, ii + 4);
    }

    b_A[ii] = atmp;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    for (i = 0; i <= ix0; i++) {
      knt = i + 3 * ix0;
      R[knt] = b_A[knt];
    }

    b_i = ix0 + 2;
    if (b_i <= 3) {
      std::memset(&R[(ix0 * 3 + b_i) + -1], 0, (4 - b_i) * sizeof(double));
    }
  }

  R[9] = b_A[9];
  R[10] = b_A[10];
  R[11] = b_A[11];
  ii = 2;
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 2; i >= 0; i--) {
    iaii = (i + i * 3) + 4;
    if (i + 1 < 3) {
      b_A[iaii - 4] = 1.0;
      if (tau[ii] != 0.0) {
        lastv = 3 - i;
        knt = iaii - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0)) {
          lastv--;
          knt--;
        }

        lastc = 2 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = iaii + (lastc - 1) * 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= (knt + lastv) - 1) {
              if (b_A[ia - 1] != 0.0) {
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
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(double));
          }

          ix0 = 0;
          b_i = iaii + 3 * (lastc - 1);
          for (k = iaii; k <= b_i; k += 3) {
            ix = iaii;
            c = 0.0;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix - 4];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        xgerc(lastv, lastc, -tau[ii], iaii - 3, work, b_A, iaii);
      }

      ix0 = iaii - 2;
      b_i = (iaii - i) - 1;
      for (k = ix0; k <= b_i; k++) {
        b_A[k - 1] *= -tau[ii];
      }
    }

    b_A[iaii - 4] = 1.0 - tau[ii];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(iaii - ix0) - 5] = 0.0;
    }

    ii--;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    Q[3 * ix0] = b_A[3 * ix0];
    knt = 3 * ix0 + 1;
    Q[knt] = b_A[knt];
    knt = 3 * ix0 + 2;
    Q[knt] = b_A[knt];
  }
}

//
// File trailer for qr.cpp
//
// [EOF]
//
