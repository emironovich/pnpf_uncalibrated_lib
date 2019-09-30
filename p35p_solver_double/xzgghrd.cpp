/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgghrd.cpp
 *
 * Code generation for function 'xzgghrd'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzgghrd.h"
#include "xzlartg.h"

/* Function Definitions */
void xzgghrd(int ilo, int ihi, creal_T A[100], creal_T Z[100])
{
  signed char b_I[100];
  int k;
  int jcol;
  int jcolp1;
  int jrow;
  double c;
  creal_T s;
  int s_re_tmp;
  int stemp_re_tmp;
  double stemp_re;
  double stemp_im;
  double A_im;
  double A_re;
  memset(&b_I[0], 0, 100U * sizeof(signed char));
  for (k = 0; k < 10; k++) {
    b_I[k + 10 * k] = 1;
  }

  for (k = 0; k < 100; k++) {
    Z[k].re = b_I[k];
    Z[k].im = 0.0;
  }

  if (ihi >= ilo + 2) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      jcolp1 = jcol + 2;
      for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
        k = jrow + 10 * jcol;
        xzlartg(A[k - 1], A[k], &c, &s, &A[(jrow + 10 * jcol) - 1]);
        A[k].re = 0.0;
        A[k].im = 0.0;
        for (k = jcolp1; k < 11; k++) {
          s_re_tmp = jrow + 10 * (k - 1);
          stemp_re_tmp = s_re_tmp - 1;
          stemp_re = c * A[stemp_re_tmp].re + (s.re * A[s_re_tmp].re - s.im *
            A[jrow + 10 * (k - 1)].im);
          stemp_im = c * A[(jrow + 10 * (k - 1)) - 1].im + (s.re * A[jrow + 10 *
            (k - 1)].im + s.im * A[jrow + 10 * (k - 1)].re);
          A_im = A[(jrow + 10 * (k - 1)) - 1].im;
          A_re = A[(jrow + 10 * (k - 1)) - 1].re;
          A[s_re_tmp].re = c * A[jrow + 10 * (k - 1)].re - (s.re * A[(jrow + 10 *
            (k - 1)) - 1].re + s.im * A[(jrow + 10 * (k - 1)) - 1].im);
          A[s_re_tmp].im = c * A[jrow + 10 * (k - 1)].im - (s.re * A_im - s.im *
            A_re);
          A[stemp_re_tmp].re = stemp_re;
          A[stemp_re_tmp].im = stemp_im;
        }

        s.re = -s.re;
        s.im = -s.im;
        for (k = 1; k <= ihi; k++) {
          s_re_tmp = (k + 10 * (jrow - 1)) - 1;
          stemp_re_tmp = (k + 10 * jrow) - 1;
          stemp_re = c * A[stemp_re_tmp].re + (s.re * A[s_re_tmp].re - s.im * A
            [(k + 10 * (jrow - 1)) - 1].im);
          stemp_im = c * A[(k + 10 * jrow) - 1].im + (s.re * A[(k + 10 * (jrow -
            1)) - 1].im + s.im * A[(k + 10 * (jrow - 1)) - 1].re);
          A_im = A[(k + 10 * jrow) - 1].im;
          A_re = A[(k + 10 * jrow) - 1].re;
          A[s_re_tmp].re = c * A[(k + 10 * (jrow - 1)) - 1].re - (s.re * A[(k +
            10 * jrow) - 1].re + s.im * A[(k + 10 * jrow) - 1].im);
          A[s_re_tmp].im = c * A[(k + 10 * (jrow - 1)) - 1].im - (s.re * A_im -
            s.im * A_re);
          A[stemp_re_tmp].re = stemp_re;
          A[stemp_re_tmp].im = stemp_im;
        }

        for (k = 0; k < 10; k++) {
          s_re_tmp = k + 10 * (jrow - 1);
          stemp_re_tmp = k + 10 * jrow;
          stemp_re = c * Z[stemp_re_tmp].re + (s.re * Z[s_re_tmp].re - s.im *
            Z[k + 10 * (jrow - 1)].im);
          stemp_im = c * Z[k + 10 * jrow].im + (s.re * Z[k + 10 * (jrow - 1)].im
            + s.im * Z[k + 10 * (jrow - 1)].re);
          A_im = Z[k + 10 * jrow].re;
          Z[s_re_tmp].re = c * Z[k + 10 * (jrow - 1)].re - (s.re * Z[k + 10 *
            jrow].re + s.im * Z[k + 10 * jrow].im);
          Z[s_re_tmp].im = c * Z[s_re_tmp].im - (s.re * Z[k + 10 * jrow].im -
            s.im * A_im);
          Z[stemp_re_tmp].re = stemp_re;
          Z[stemp_re_tmp].im = stemp_im;
        }
      }
    }
  }
}

/* End of code generation (xzgghrd.cpp) */
