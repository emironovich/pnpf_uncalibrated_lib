/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeev.cpp
 *
 * Code generation for function 'xzgeev'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzgeev.h"
#include "xzlartg.h"
#include "xzhgeqz.h"
#include "xzggbal.h"
#include "xgeqp3.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
void xzgeev(const creal32_T A_data[], const int A_size[2], int *info, creal32_T
            alpha1_data[], int alpha1_size[1], creal32_T beta1_data[], int
            beta1_size[1])
{
  int At_size[2];
  int loop_ub_tmp;
  creal32_T At_data[144];
  float anrm;
  int k;
  boolean_T exitg1;
  float absxk;
  boolean_T ilascl;
  float anrmto;
  int ilo;
  int ihi;
  int rscale_data[12];
  int rscale_size[1];
  float ctoc;
  int n;
  boolean_T notdone;
  int jcol;
  float stemp_im;
  float cto1;
  int jcolp1;
  int jrow;
  float a;
  creal32_T s;
  int stemp_re_tmp;
  At_size[0] = A_size[0];
  At_size[1] = A_size[1];
  loop_ub_tmp = A_size[0] * A_size[1];
  if (0 <= loop_ub_tmp - 1) {
    memcpy(&At_data[0], &A_data[0], (unsigned int)(loop_ub_tmp * static_cast<int>
            (sizeof(creal32_T))));
  }

  *info = 0;
  anrm = 0.0F;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= loop_ub_tmp - 1)) {
    absxk = rt_hypotf_snf(A_data[k].re, A_data[k].im);
    if (rtIsNaNF(absxk)) {
      anrm = rtNaNF;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      k++;
    }
  }

  if (rtIsInfF(anrm) || rtIsNaNF(anrm)) {
    alpha1_size[0] = A_size[0];
    loop_ub_tmp = A_size[0];
    for (k = 0; k < loop_ub_tmp; k++) {
      alpha1_data[k].re = rtNaNF;
      alpha1_data[k].im = 0.0F;
    }

    beta1_size[0] = A_size[0];
    loop_ub_tmp = A_size[0];
    for (k = 0; k < loop_ub_tmp; k++) {
      beta1_data[k].re = rtNaNF;
      beta1_data[k].im = 0.0F;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0F) && (anrm < 9.09494702E-13F)) {
      anrmto = 9.09494702E-13F;
      ilascl = true;
    } else {
      if (anrm > 1.09951163E+12F) {
        anrmto = 1.09951163E+12F;
        ilascl = true;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        stemp_im = absxk * 1.97215226E-31F;
        cto1 = ctoc / 5.0706024E+30F;
        if ((stemp_im > ctoc) && (ctoc != 0.0F)) {
          a = 1.97215226E-31F;
          absxk = stemp_im;
        } else if (cto1 > absxk) {
          a = 5.0706024E+30F;
          ctoc = cto1;
        } else {
          a = ctoc / absxk;
          notdone = false;
        }

        loop_ub_tmp = At_size[0] * At_size[1] - 1;
        for (k = 0; k <= loop_ub_tmp; k++) {
          At_data[k].re *= a;
          At_data[k].im *= a;
        }
      }
    }

    xzggbal(At_data, At_size, &ilo, &ihi, rscale_data, rscale_size);
    n = At_size[0];
    if ((At_size[0] > 1) && (ihi >= ilo + 2)) {
      for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
        jcolp1 = jcol + 2;
        for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
          loop_ub_tmp = jrow + At_size[0] * jcol;
          xzlartg(At_data[loop_ub_tmp - 1], At_data[loop_ub_tmp], &absxk, &s,
                  &At_data[(jrow + At_size[0] * jcol) - 1]);
          At_data[loop_ub_tmp].re = 0.0F;
          At_data[loop_ub_tmp].im = 0.0F;
          for (loop_ub_tmp = jcolp1; loop_ub_tmp <= n; loop_ub_tmp++) {
            k = jrow + At_size[0] * (loop_ub_tmp - 1);
            stemp_re_tmp = k - 1;
            ctoc = absxk * At_data[stemp_re_tmp].re + (s.re * At_data[k].re -
              s.im * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)].im);
            stemp_im = absxk * At_data[(jrow + At_size[0] * (loop_ub_tmp - 1)) -
              1].im + (s.re * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)].im
                       + s.im * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)].
                       re);
            cto1 = At_data[(jrow + At_size[0] * (loop_ub_tmp - 1)) - 1].re;
            At_data[k].re = absxk * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)]
              .re - (s.re * At_data[(jrow + At_size[0] * (loop_ub_tmp - 1)) - 1]
                     .re + s.im * At_data[(jrow + At_size[0] * (loop_ub_tmp - 1))
                     - 1].im);
            At_data[k].im = absxk * At_data[k].im - (s.re * At_data[(jrow +
              At_size[0] * (loop_ub_tmp - 1)) - 1].im - s.im * cto1);
            At_data[stemp_re_tmp].re = ctoc;
            At_data[stemp_re_tmp].im = stemp_im;
          }

          s.re = -s.re;
          s.im = -s.im;
          for (loop_ub_tmp = 1; loop_ub_tmp <= ihi; loop_ub_tmp++) {
            k = (loop_ub_tmp + At_size[0] * (jrow - 1)) - 1;
            stemp_re_tmp = (loop_ub_tmp + At_size[0] * jrow) - 1;
            ctoc = absxk * At_data[stemp_re_tmp].re + (s.re * At_data[k].re -
              s.im * At_data[(loop_ub_tmp + At_size[0] * (jrow - 1)) - 1].im);
            stemp_im = absxk * At_data[(loop_ub_tmp + At_size[0] * jrow) - 1].im
              + (s.re * At_data[(loop_ub_tmp + At_size[0] * (jrow - 1)) - 1].im
                 + s.im * At_data[(loop_ub_tmp + At_size[0] * (jrow - 1)) - 1].
                 re);
            cto1 = At_data[stemp_re_tmp].re;
            At_data[k].re = absxk * At_data[k].re - (s.re * At_data[(loop_ub_tmp
              + At_size[0] * jrow) - 1].re + s.im * At_data[(loop_ub_tmp +
              At_size[0] * jrow) - 1].im);
            At_data[k].im = absxk * At_data[k].im - (s.re * At_data[(loop_ub_tmp
              + At_size[0] * jrow) - 1].im - s.im * cto1);
            At_data[stemp_re_tmp].re = ctoc;
            At_data[stemp_re_tmp].im = stemp_im;
          }
        }
      }
    }

    xzhgeqz(At_data, At_size, ilo, ihi, info, alpha1_data, alpha1_size,
            beta1_data, beta1_size);
    if ((*info == 0) && ilascl) {
      notdone = true;
      while (notdone) {
        stemp_im = anrmto * 1.97215226E-31F;
        cto1 = anrm / 5.0706024E+30F;
        if ((stemp_im > anrm) && (anrm != 0.0F)) {
          a = 1.97215226E-31F;
          anrmto = stemp_im;
        } else if (cto1 > anrmto) {
          a = 5.0706024E+30F;
          anrm = cto1;
        } else {
          a = anrm / anrmto;
          notdone = false;
        }

        loop_ub_tmp = alpha1_size[0];
        for (k = 0; k < loop_ub_tmp; k++) {
          alpha1_data[k].re *= a;
          alpha1_data[k].im *= a;
        }
      }
    }
  }
}

/* End of code generation (xzgeev.cpp) */
