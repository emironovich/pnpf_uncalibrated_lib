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
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "svd.h"
#include "xrotg.h"
#include "sqrt.h"
#include "xaxpy.h"
#include "xnrm2.h"

/* Function Definitions */
void svd(const double A[8], double U[2])
{
  double b_A[8];
  double s[3];
  double e[4];
  double work[2];
  int q;
  int m;
  int qp1;
  int qq;
  int nmqp1;
  int iy;
  boolean_T apply_transform;
  double snorm;
  double t;
  double nrm;
  int jj;
  double scale;
  int kend;
  double absxk;
  int k;
  int ix;
  int i12;
  int exitg1;
  boolean_T exitg2;
  double sm;
  double sqds;
  memcpy(&b_A[0], &A[0], sizeof(double) << 3);
  s[0] = 0.0;
  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 0.0;
  e[3] = 0.0;
  work[0] = 0.0;
  for (q = 0; q < 2; q++) {
    qp1 = q + 2;
    m = q + (q << 1);
    qq = m + 1;
    nmqp1 = 1 - q;
    apply_transform = false;
    if (q + 1 <= 1) {
      nrm = 0.0;
      scale = 3.3121686421112381E-170;
      kend = qq + 1;
      for (k = qq; k <= kend; k++) {
        absxk = std::abs(b_A[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          nrm = 1.0 + nrm * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          nrm += t * t;
        }
      }

      nrm = scale * std::sqrt(nrm);
      if (nrm > 0.0) {
        apply_transform = true;
        if (b_A[qq - 1] < 0.0) {
          t = -nrm;
        } else {
          t = nrm;
        }

        if (std::abs(t) >= 1.0020841800044864E-292) {
          nrm = 1.0 / t;
          i12 = qq + 1;
          for (k = qq; k <= i12; k++) {
            b_A[k - 1] *= nrm;
          }
        } else {
          i12 = qq + 1;
          for (k = qq; k <= i12; k++) {
            b_A[k - 1] /= t;
          }
        }

        b_A[qq - 1]++;
        s[0] = -t;
      } else {
        s[0] = 0.0;
      }
    }

    for (jj = qp1; jj < 5; jj++) {
      kend = q + ((jj - 1) << 1);
      if (apply_transform) {
        ix = qq;
        iy = kend;
        t = 0.0;
        for (k = 0; k <= nmqp1; k++) {
          t += b_A[ix - 1] * b_A[iy];
          ix++;
          iy++;
        }

        nrm = -(t / b_A[m]);
        if (!(nrm == 0.0)) {
          ix = qq - 1;
          iy = kend;
          i12 = 1 - q;
          for (k = 0; k <= i12; k++) {
            b_A[iy] += nrm * b_A[ix];
            ix++;
            iy++;
          }
        }
      }

      e[jj - 1] = b_A[kend];
    }

    nrm = e_xnrm2(3 - q, e, q + 2);
    if (nrm == 0.0) {
      e[q] = 0.0;
    } else {
      if (e[q + 1] < 0.0) {
        e[q] = -nrm;
      } else {
        e[q] = nrm;
      }

      nrm = e[q];
      if (std::abs(e[q]) >= 1.0020841800044864E-292) {
        nrm = 1.0 / e[q];
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
        work[1] = 0.0;
        xaxpy(1 - q, e[1], b_A, 4, work);
        xaxpy(1 - q, e[2], b_A, 6, work);
        xaxpy(1 - q, e[3], b_A, 8, work);
        b_xaxpy(1 - q, -e[1] / e[1], work, b_A, 4);
        b_xaxpy(1 - q, -e[2] / e[1], work, b_A, 6);
        b_xaxpy(1 - q, -e[3] / e[1], work, b_A, 8);
      }
    }
  }

  m = 1;
  s[1] = b_A[3];
  s[2] = 0.0;
  e[2] = 0.0;
  iy = 0;
  snorm = 0.0;
  t = s[0];
  if (s[0] != 0.0) {
    nrm = std::abs(s[0]);
    absxk = s[0] / nrm;
    t = nrm;
    s[0] = nrm;
    e[0] /= absxk;
  }

  if (e[0] != 0.0) {
    nrm = std::abs(e[0]);
    absxk = nrm / e[0];
    e[0] = nrm;
    s[1] = b_A[3] * absxk;
  }

  nrm = std::abs(t);
  absxk = std::abs(e[0]);
  if ((nrm > absxk) || rtIsNaN(absxk)) {
    absxk = nrm;
  }

  if (!rtIsNaN(absxk)) {
    snorm = absxk;
  }

  t = s[1];
  if (s[1] != 0.0) {
    nrm = std::abs(s[1]);
    absxk = s[1] / nrm;
    t = nrm;
    s[1] = nrm;
    e[1] /= absxk;
  }

  if (e[1] != 0.0) {
    nrm = std::abs(e[1]);
    absxk = nrm / e[1];
    e[1] = nrm;
    s[2] = 0.0 * absxk;
  }

  nrm = std::abs(t);
  absxk = std::abs(e[1]);
  if ((nrm > absxk) || rtIsNaN(absxk)) {
    absxk = nrm;
  }

  if ((!(snorm > absxk)) && (!rtIsNaN(absxk))) {
    snorm = absxk;
  }

  if (s[2] != 0.0) {
    s[2] = rtNaN;
  }

  if (!(snorm > 0.0)) {
    snorm = 0.0;
  }

  while ((m + 2 > 0) && (iy < 75)) {
    kend = m;
    do {
      exitg1 = 0;
      q = kend + 1;
      if (kend + 1 == 0) {
        exitg1 = 1;
      } else {
        nrm = std::abs(e[kend]);
        if ((nrm <= 2.2204460492503131E-16 * (std::abs(s[kend]) + std::abs
              (s[kend + 1]))) || (nrm <= 1.0020841800044864E-292) || ((iy > 20) &&
             (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[kend] = 0.0;
          exitg1 = 1;
        } else {
          kend--;
        }
      }
    } while (exitg1 == 0);

    if (kend + 1 == m + 1) {
      kend = 4;
    } else {
      ix = m + 2;
      jj = m + 2;
      exitg2 = false;
      while ((!exitg2) && (jj >= kend + 1)) {
        ix = jj;
        if (jj == kend + 1) {
          exitg2 = true;
        } else {
          nrm = 0.0;
          if (jj < m + 2) {
            nrm = std::abs(e[jj - 1]);
          }

          if (jj > kend + 2) {
            nrm += std::abs(e[jj - 2]);
          }

          absxk = std::abs(s[jj - 1]);
          if ((absxk <= 2.2204460492503131E-16 * nrm) || (absxk <=
               1.0020841800044864E-292)) {
            s[jj - 1] = 0.0;
            exitg2 = true;
          } else {
            jj--;
          }
        }
      }

      if (ix == kend + 1) {
        kend = 3;
      } else if (ix == m + 2) {
        kend = 1;
      } else {
        kend = 2;
        q = ix;
      }
    }

    switch (kend) {
     case 1:
      absxk = e[m];
      e[m] = 0.0;
      i12 = m + 1;
      for (k = i12; k >= q + 1; k--) {
        xrotg(&s[k - 1], &absxk, &sm, &sqds);
        if (k > q + 1) {
          absxk = -sqds * e[0];
          e[0] *= sm;
        }
      }
      break;

     case 2:
      absxk = e[q - 1];
      e[q - 1] = 0.0;
      for (k = q + 1; k <= m + 2; k++) {
        xrotg(&s[k - 1], &absxk, &sm, &sqds);
        t = e[k - 1];
        absxk = -sqds * t;
        e[k - 1] = t * sm;
      }
      break;

     case 3:
      kend = m + 1;
      scale = std::abs(s[m + 1]);
      absxk = std::abs(s[m]);
      if ((!(scale > absxk)) && (!rtIsNaN(absxk))) {
        scale = absxk;
      }

      absxk = std::abs(e[m]);
      if ((!(scale > absxk)) && (!rtIsNaN(absxk))) {
        scale = absxk;
      }

      absxk = std::abs(s[q]);
      if ((!(scale > absxk)) && (!rtIsNaN(absxk))) {
        scale = absxk;
      }

      absxk = std::abs(e[q]);
      if ((!(scale > absxk)) && (!rtIsNaN(absxk))) {
        scale = absxk;
      }

      sm = s[m + 1] / scale;
      nrm = s[m] / scale;
      absxk = e[m] / scale;
      sqds = s[q] / scale;
      t = ((nrm + sm) * (nrm - sm) + absxk * absxk) / 2.0;
      nrm = sm * absxk;
      nrm *= nrm;
      if ((t != 0.0) || (nrm != 0.0)) {
        absxk = t * t + nrm;
        b_sqrt(&absxk);
        if (t < 0.0) {
          absxk = -absxk;
        }

        absxk = nrm / (t + absxk);
      } else {
        absxk = 0.0;
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

      e[m] = absxk;
      iy++;
      break;

     default:
      if (s[q] < 0.0) {
        s[q] = -s[q];
      }

      qp1 = q + 1;
      while ((q + 1 < 3) && (s[q] < s[qp1])) {
        nrm = s[q];
        s[q] = s[qp1];
        s[qp1] = nrm;
        q = qp1;
        qp1++;
      }

      iy = 0;
      m--;
      break;
    }
  }

  U[0] = s[0];
  U[1] = s[1];
}

/* End of code generation (svd.cpp) */
