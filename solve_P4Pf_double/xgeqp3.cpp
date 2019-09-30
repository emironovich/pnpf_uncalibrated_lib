/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgeqp3.cpp
 *
 * Code generation for function 'xgeqp3'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xgeqp3.h"
#include "xnrm2.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
void xgeqp3(double A[36], double tau[3], int jpvt[3])
{
  int k;
  int j;
  int i;
  double work[3];
  int ip1;
  double smax;
  int i_i;
  double scale;
  int kend;
  int pvt;
  int ix;
  double absxk;
  double vn1[3];
  double vn2[3];
  double t;
  int iy;
  int i34;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int exitg1;
  int i35;
  k = 1;
  for (j = 0; j < 3; j++) {
    jpvt[j] = 1 + j;
    work[j] = 0.0;
    smax = 0.0;
    scale = 3.3121686421112381E-170;
    kend = k + 11;
    for (pvt = k; pvt <= kend; pvt++) {
      absxk = std::abs(A[pvt - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        smax = 1.0 + smax * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        smax += t * t;
      }
    }

    smax = scale * std::sqrt(smax);
    vn1[j] = smax;
    vn2[j] = smax;
    k += 12;
  }

  for (i = 0; i < 3; i++) {
    ip1 = i + 2;
    i_i = i + i * 12;
    kend = 3 - i;
    pvt = 1;
    if (3 - i > 1) {
      ix = i;
      smax = std::abs(vn1[i]);
      for (k = 2; k <= kend; k++) {
        ix++;
        scale = std::abs(vn1[ix]);
        if (scale > smax) {
          pvt = k;
          smax = scale;
        }
      }
    }

    pvt = (i + pvt) - 1;
    if (pvt != i) {
      ix = 12 * pvt;
      iy = 12 * i;
      for (k = 0; k < 12; k++) {
        smax = A[ix];
        A[ix] = A[iy];
        A[iy] = smax;
        ix++;
        iy++;
      }

      kend = jpvt[pvt];
      jpvt[pvt] = jpvt[i];
      jpvt[i] = kend;
      vn1[pvt] = vn1[i];
      vn2[pvt] = vn2[i];
    }

    absxk = A[i_i];
    kend = i_i + 2;
    tau[i] = 0.0;
    smax = d_xnrm2(11 - i, A, i_i + 2);
    if (smax != 0.0) {
      scale = rt_hypotd_snf(A[i_i], smax);
      if (A[i_i] >= 0.0) {
        scale = -scale;
      }

      if (std::abs(scale) < 1.0020841800044864E-292) {
        pvt = -1;
        i34 = (i_i - i) + 12;
        do {
          pvt++;
          for (k = kend; k <= i34; k++) {
            A[k - 1] *= 9.9792015476736E+291;
          }

          scale *= 9.9792015476736E+291;
          absxk *= 9.9792015476736E+291;
        } while (!(std::abs(scale) >= 1.0020841800044864E-292));

        scale = rt_hypotd_snf(absxk, d_xnrm2(11 - i, A, i_i + 2));
        if (absxk >= 0.0) {
          scale = -scale;
        }

        tau[i] = (scale - absxk) / scale;
        smax = 1.0 / (absxk - scale);
        for (k = kend; k <= i34; k++) {
          A[k - 1] *= smax;
        }

        for (k = 0; k <= pvt; k++) {
          scale *= 1.0020841800044864E-292;
        }

        absxk = scale;
      } else {
        tau[i] = (scale - A[i_i]) / scale;
        smax = 1.0 / (A[i_i] - scale);
        i34 = (i_i - i) + 12;
        for (k = kend; k <= i34; k++) {
          A[k - 1] *= smax;
        }

        absxk = scale;
      }
    }

    A[i_i] = absxk;
    if (i + 1 < 3) {
      absxk = A[i_i];
      A[i_i] = 1.0;
      k = (i + (i + 1) * 12) + 1;
      if (tau[i] != 0.0) {
        lastv = 12 - i;
        kend = (i_i - i) + 11;
        while ((lastv > 0) && (A[kend] == 0.0)) {
          lastv--;
          kend--;
        }

        lastc = 1 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc + 1 > 0)) {
          kend = k + lastc * 12;
          pvt = kend;
          do {
            exitg1 = 0;
            if (pvt <= (kend + lastv) - 1) {
              if (A[pvt - 1] != 0.0) {
                exitg1 = 1;
              } else {
                pvt++;
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
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (0 <= lastc) {
            memset(&work[0], 0, (unsigned int)((lastc + 1) * static_cast<int>
                    (sizeof(double))));
          }

          iy = 0;
          i34 = k + 12 * lastc;
          for (kend = k; kend <= i34; kend += 12) {
            ix = i_i;
            smax = 0.0;
            i35 = (kend + lastv) - 1;
            for (pvt = kend; pvt <= i35; pvt++) {
              smax += A[pvt - 1] * A[ix];
              ix++;
            }

            work[iy] += smax;
            iy++;
          }
        }

        if (!(-tau[i] == 0.0)) {
          kend = k - 1;
          pvt = 0;
          for (j = 0; j <= lastc; j++) {
            if (work[pvt] != 0.0) {
              smax = work[pvt] * -tau[i];
              ix = i_i;
              i34 = kend + 1;
              i35 = lastv + kend;
              for (iy = i34; iy <= i35; iy++) {
                A[iy - 1] += A[ix] * smax;
                ix++;
              }
            }

            pvt++;
            kend += 12;
          }
        }
      }

      A[i_i] = absxk;
    }

    for (j = ip1; j < 4; j++) {
      smax = vn1[j - 1];
      if (smax != 0.0) {
        kend = i + 12 * (j - 1);
        scale = std::abs(A[kend]) / smax;
        scale = 1.0 - scale * scale;
        if (scale < 0.0) {
          scale = 0.0;
        }

        absxk = smax / vn2[j - 1];
        absxk = scale * (absxk * absxk);
        if (absxk <= 1.4901161193847656E-8) {
          smax = d_xnrm2(11 - i, A, kend + 2);
          vn1[j - 1] = smax;
          vn2[j - 1] = smax;
        } else {
          vn1[j - 1] = smax * std::sqrt(scale);
        }
      }
    }
  }
}

/* End of code generation (xgeqp3.cpp) */
