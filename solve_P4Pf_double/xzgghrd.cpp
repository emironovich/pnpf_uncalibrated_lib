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
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzgghrd.h"
#include "xzlartg.h"

/* Function Definitions */
void xzgghrd(int ilo, int ihi, creal_T A_data[], int A_size[2])
{
  int n;
  int jcol;
  int jcolp1;
  int jrow;
  int A_data_tmp;
  double c;
  creal_T s;
  int s_re_tmp;
  int stemp_re_tmp;
  double stemp_re;
  double stemp_im;
  double A_data_im;
  double A_data_re;
  n = A_size[0];
  if ((A_size[0] > 1) && (ihi >= ilo + 2)) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      jcolp1 = jcol + 2;
      for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
        A_data_tmp = jrow + A_size[0] * jcol;
        xzlartg(A_data[A_data_tmp - 1], A_data[A_data_tmp], &c, &s, &A_data
                [(jrow + A_size[0] * jcol) - 1]);
        A_data[A_data_tmp].re = 0.0;
        A_data[A_data_tmp].im = 0.0;
        for (A_data_tmp = jcolp1; A_data_tmp <= n; A_data_tmp++) {
          s_re_tmp = jrow + A_size[0] * (A_data_tmp - 1);
          stemp_re_tmp = s_re_tmp - 1;
          stemp_re = c * A_data[stemp_re_tmp].re + (s.re * A_data[s_re_tmp].re -
            s.im * A_data[jrow + A_size[0] * (A_data_tmp - 1)].im);
          stemp_im = c * A_data[(jrow + A_size[0] * (A_data_tmp - 1)) - 1].im +
            (s.re * A_data[jrow + A_size[0] * (A_data_tmp - 1)].im + s.im *
             A_data[jrow + A_size[0] * (A_data_tmp - 1)].re);
          A_data_im = A_data[(jrow + A_size[0] * (A_data_tmp - 1)) - 1].im;
          A_data_re = A_data[(jrow + A_size[0] * (A_data_tmp - 1)) - 1].re;
          A_data[s_re_tmp].re = c * A_data[jrow + A_size[0] * (A_data_tmp - 1)].
            re - (s.re * A_data[(jrow + A_size[0] * (A_data_tmp - 1)) - 1].re +
                  s.im * A_data[(jrow + A_size[0] * (A_data_tmp - 1)) - 1].im);
          A_data[s_re_tmp].im = c * A_data[jrow + A_size[0] * (A_data_tmp - 1)].
            im - (s.re * A_data_im - s.im * A_data_re);
          A_data[stemp_re_tmp].re = stemp_re;
          A_data[stemp_re_tmp].im = stemp_im;
        }

        s.re = -s.re;
        s.im = -s.im;
        for (A_data_tmp = 1; A_data_tmp <= ihi; A_data_tmp++) {
          s_re_tmp = (A_data_tmp + A_size[0] * (jrow - 1)) - 1;
          stemp_re_tmp = (A_data_tmp + A_size[0] * jrow) - 1;
          stemp_re = c * A_data[stemp_re_tmp].re + (s.re * A_data[s_re_tmp].re -
            s.im * A_data[(A_data_tmp + A_size[0] * (jrow - 1)) - 1].im);
          stemp_im = c * A_data[(A_data_tmp + A_size[0] * jrow) - 1].im + (s.re *
            A_data[(A_data_tmp + A_size[0] * (jrow - 1)) - 1].im + s.im *
            A_data[(A_data_tmp + A_size[0] * (jrow - 1)) - 1].re);
          A_data_im = A_data[(A_data_tmp + A_size[0] * jrow) - 1].im;
          A_data_re = A_data[(A_data_tmp + A_size[0] * jrow) - 1].re;
          A_data[s_re_tmp].re = c * A_data[(A_data_tmp + A_size[0] * (jrow - 1))
            - 1].re - (s.re * A_data[(A_data_tmp + A_size[0] * jrow) - 1].re +
                       s.im * A_data[(A_data_tmp + A_size[0] * jrow) - 1].im);
          A_data[s_re_tmp].im = c * A_data[(A_data_tmp + A_size[0] * (jrow - 1))
            - 1].im - (s.re * A_data_im - s.im * A_data_re);
          A_data[stemp_re_tmp].re = stemp_re;
          A_data[stemp_re_tmp].im = stemp_im;
        }
      }
    }
  }
}

/* End of code generation (xzgghrd.cpp) */
