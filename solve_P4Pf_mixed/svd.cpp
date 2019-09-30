/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * svd.cpp
 *
 * Code generation for function 'svd'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "svd.h"
#include "xrotg.h"
#include "sqrt.h"
#include "xaxpy.h"
#include "xnrm2.h"

/* Function Definitions */
void svd(const float A[8], float U[2])
{
  int kend;
  float s[3];
  float b_A[8];
  float e[4];
  float work[2];
  int q;
  int qq;
  int qp1;
  int nmqp1;
  boolean_T apply_transform;
  int iy;
  float snorm;
  float nrm;
  float sm;
  float scale;
  float absxk;
  int qjj;
  int k;
  int ix;
  float t;
  int exitg1;
  boolean_T exitg2;
  float sqds;
  for (kend = 0; kend < 8; kend++) {
    b_A[kend] = A[kend];
  }

  s[0] = 0.0F;
  e[0] = 0.0F;
  e[1] = 0.0F;
  e[2] = 0.0F;
  e[3] = 0.0F;
  work[0] = 0.0F;
  for (q = 0; q < 2; q++) {
    qp1 = q + 2;
    qq = (q + (q << 1)) + 1;
    nmqp1 = 1 - q;
    apply_transform = false;
    if (q + 1 <= 1) {
      nrm = 0.0F;
      scale = 1.29246971E-26F;
      kend = qq + 1;
      for (k = qq; k <= kend; k++) {
        absxk = std::abs(b_A[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          nrm = 1.0F + nrm * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          nrm += t * t;
        }
      }

      nrm = scale * std::sqrt(nrm);
      if (nrm > 0.0F) {
        apply_transform = true;
        if (b_A[qq - 1] < 0.0F) {
          sm = -nrm;
        } else {
          sm = nrm;
        }

        if (std::abs(sm) >= 9.86076132E-32F) {
          nrm = 1.0F / sm;
          kend = qq + 1;
          for (k = qq; k <= kend; k++) {
            b_A[k - 1] *= nrm;
          }
        } else {
          kend = qq + 1;
          for (k = qq; k <= kend; k++) {
            b_A[k - 1] /= sm;
          }
        }

        b_A[qq - 1]++;
        s[0] = -sm;
      } else {
        s[0] = 0.0F;
      }
    }

    for (kend = qp1; kend < 5; kend++) {
      qjj = q + ((kend - 1) << 1);
      if (apply_transform) {
        ix = qq;
        iy = qjj + 1;
        t = 0.0F;
        for (k = 0; k <= nmqp1; k++) {
          t += b_A[ix - 1] * b_A[iy - 1];
          ix++;
          iy++;
        }

        xaxpy(2 - q, -(t / b_A[q + (q << 1)]), qq, b_A, qjj + 1);
      }

      e[kend - 1] = b_A[qjj];
    }

    nrm = e_xnrm2(3 - q, e, q + 2);
    if (nrm == 0.0F) {
      e[q] = 0.0F;
    } else {
      if (e[q + 1] < 0.0F) {
        e[q] = -nrm;
      } else {
        e[q] = nrm;
      }

      nrm = e[q];
      if (std::abs(e[q]) >= 9.86076132E-32F) {
        nrm = 1.0F / e[q];
        for (k = qp1; k < 5; k++) {
          e[k - 1] *= nrm;
        }
      } else {
        for (k = qp1; k < 5; k++) {
          e[k - 1] /= nrm;
        }
      }

      e[q + 1]++;
      e[q] = -e[q];
      if (q + 2 <= 2) {
        work[1] = 0.0F;
        b_xaxpy(1 - q, e[1], b_A, 4, work);
        b_xaxpy(1 - q, e[2], b_A, 6, work);
        b_xaxpy(1 - q, e[3], b_A, 8, work);
        c_xaxpy(1 - q, -e[1] / e[1], work, b_A, 4);
        c_xaxpy(1 - q, -e[2] / e[1], work, b_A, 6);
        c_xaxpy(1 - q, -e[3] / e[1], work, b_A, 8);
      }
    }
  }

  qq = 1;
  s[1] = b_A[3];
  s[2] = 0.0F;
  e[2] = 0.0F;
  iy = 0;
  snorm = 0.0F;
  sm = s[0];
  if (s[0] != 0.0F) {
    absxk = std::abs(s[0]);
    nrm = s[0] / absxk;
    sm = absxk;
    s[0] = absxk;
    e[0] /= nrm;
  }

  if (e[0] != 0.0F) {
    absxk = std::abs(e[0]);
    t = e[0];
    e[0] = absxk;
    s[1] = b_A[3] * (absxk / t);
  }

  nrm = std::abs(sm);
  absxk = std::abs(e[0]);
  if ((nrm > absxk) || rtIsNaNF(absxk)) {
    absxk = nrm;
  }

  if (!rtIsNaNF(absxk)) {
    snorm = absxk;
  }

  sm = s[1];
  if (s[1] != 0.0F) {
    absxk = std::abs(s[1]);
    nrm = s[1] / absxk;
    sm = absxk;
    s[1] = absxk;
    e[1] /= nrm;
  }

  if (e[1] != 0.0F) {
    nrm = std::abs(e[1]);
    t = e[1];
    e[1] = nrm;
    s[2] = 0.0F * (nrm / t);
  }

  nrm = std::abs(sm);
  absxk = std::abs(e[1]);
  if ((nrm > absxk) || rtIsNaNF(absxk)) {
    absxk = nrm;
  }

  if ((!(snorm > absxk)) && (!rtIsNaNF(absxk))) {
    snorm = absxk;
  }

  if (s[2] != 0.0F) {
    s[2] = rtNaNF;
  }

  if (!(snorm > 0.0F)) {
    snorm = 0.0F;
  }

  while ((qq + 2 > 0) && (iy < 75)) {
    kend = qq;
    do {
      exitg1 = 0;
      q = kend + 1;
      if (kend + 1 == 0) {
        exitg1 = 1;
      } else {
        nrm = std::abs(e[kend]);
        if ((nrm <= 1.1920929E-7F * (std::abs(s[kend]) + std::abs(s[kend + 1])))
            || (nrm <= 9.86076132E-32F) || ((iy > 20) && (nrm <= 1.1920929E-7F *
              snorm))) {
          e[kend] = 0.0F;
          exitg1 = 1;
        } else {
          kend--;
        }
      }
    } while (exitg1 == 0);

    if (kend + 1 == qq + 1) {
      kend = 4;
    } else {
      ix = qq + 2;
      qjj = qq + 2;
      exitg2 = false;
      while ((!exitg2) && (qjj >= kend + 1)) {
        ix = qjj;
        if (qjj == kend + 1) {
          exitg2 = true;
        } else {
          nrm = 0.0F;
          if (qjj < qq + 2) {
            nrm = std::abs(e[qjj - 1]);
          }

          if (qjj > kend + 2) {
            nrm += std::abs(e[qjj - 2]);
          }

          absxk = std::abs(s[qjj - 1]);
          if ((absxk <= 1.1920929E-7F * nrm) || (absxk <= 9.86076132E-32F)) {
            s[qjj - 1] = 0.0F;
            exitg2 = true;
          } else {
            qjj--;
          }
        }
      }

      if (ix == kend + 1) {
        kend = 3;
      } else if (ix == qq + 2) {
        kend = 1;
      } else {
        kend = 2;
        q = ix;
      }
    }

    switch (kend) {
     case 1:
      absxk = e[qq];
      e[qq] = 0.0F;
      kend = qq + 1;
      for (k = kend; k >= q + 1; k--) {
        xrotg(&s[k - 1], &absxk, &sm, &sqds);
        if (k > q + 1) {
          absxk = -sqds * e[0];
          e[0] *= sm;
        }
      }
      break;

     case 2:
      absxk = e[q - 1];
      e[q - 1] = 0.0F;
      for (k = q + 1; k <= qq + 2; k++) {
        xrotg(&s[k - 1], &absxk, &sm, &sqds);
        t = e[k - 1];
        absxk = -sqds * t;
        e[k - 1] = t * sm;
      }
      break;

     case 3:
      kend = qq + 1;
      scale = std::abs(s[qq + 1]);
      absxk = std::abs(s[qq]);
      if ((!(scale > absxk)) && (!rtIsNaNF(absxk))) {
        scale = absxk;
      }

      absxk = std::abs(e[qq]);
      if ((!(scale > absxk)) && (!rtIsNaNF(absxk))) {
        scale = absxk;
      }

      absxk = std::abs(s[q]);
      if ((!(scale > absxk)) && (!rtIsNaNF(absxk))) {
        scale = absxk;
      }

      absxk = std::abs(e[q]);
      if ((!(scale > absxk)) && (!rtIsNaNF(absxk))) {
        scale = absxk;
      }

      sm = s[qq + 1] / scale;
      nrm = s[qq] / scale;
      absxk = e[qq] / scale;
      sqds = s[q] / scale;
      t = ((nrm + sm) * (nrm - sm) + absxk * absxk) / 2.0F;
      nrm = sm * absxk;
      nrm *= nrm;
      if ((t != 0.0F) || (nrm != 0.0F)) {
        absxk = t * t + nrm;
        b_sqrt(&absxk);
        if (t < 0.0F) {
          absxk = -absxk;
        }

        absxk = nrm / (t + absxk);
      } else {
        absxk = 0.0F;
      }

      absxk += (sqds + sm) * (sqds - sm);
      nrm = sqds * (e[q] / scale);
      for (k = q + 1; k <= kend; k++) {
        xrotg(&absxk, &nrm, &sm, &sqds);
        if (k > q + 1) {
          e[0] = absxk;
        }

        t = e[k - 1];
        nrm = s[k - 1];
        e[k - 1] = sm * t - sqds * nrm;
        absxk = sqds * s[k];
        s[k] *= sm;
        s[k - 1] = sm * nrm + sqds * t;
        xrotg(&s[k - 1], &absxk, &sm, &sqds);
        t = e[k - 1];
        absxk = sm * t + sqds * s[k];
        s[k] = -sqds * t + sm * s[k];
        nrm = sqds * e[k];
        e[k] *= sm;
      }

      e[qq] = absxk;
      iy++;
      break;

     default:
      if (s[q] < 0.0F) {
        s[q] = -s[q];
      }

      qp1 = q + 1;
      while ((q + 1 < 3) && (s[q] < s[qp1])) {
        absxk = s[q];
        s[q] = s[qp1];
        s[qp1] = absxk;
        q = qp1;
        qp1++;
      }

      iy = 0;
      qq--;
      break;
    }
  }

  U[0] = s[0];
  U[1] = s[1];
}

/* End of code generation (svd.cpp) */
