//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzgghrd.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xzgghrd.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include "xzlartg.h"
#include <cstring>

// Function Definitions

//
// Arguments    : int ilo
//                int ihi
//                creal_T A[100]
//                creal_T Z[100]
// Return Type  : void
//
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
  double d;
  double d1;
  double d2;
  std::memset(&b_I[0], 0, 100U * sizeof(signed char));
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
            A[s_re_tmp].im);
          stemp_im = c * A[stemp_re_tmp].im + (s.re * A[s_re_tmp].im + s.im *
            A[s_re_tmp].re);
          d = A[stemp_re_tmp].im;
          d1 = A[stemp_re_tmp].re;
          A[s_re_tmp].re = c * A[s_re_tmp].re - (s.re * A[stemp_re_tmp].re +
            s.im * A[stemp_re_tmp].im);
          A[s_re_tmp].im = c * A[s_re_tmp].im - (s.re * d - s.im * d1);
          A[stemp_re_tmp].re = stemp_re;
          A[stemp_re_tmp].im = stemp_im;
        }

        s.re = -s.re;
        s.im = -s.im;
        for (k = 1; k <= ihi; k++) {
          s_re_tmp = (k + 10 * (jrow - 1)) - 1;
          stemp_re_tmp = (k + 10 * jrow) - 1;
          stemp_re = c * A[stemp_re_tmp].re + (s.re * A[s_re_tmp].re - s.im *
            A[s_re_tmp].im);
          stemp_im = c * A[stemp_re_tmp].im + (s.re * A[s_re_tmp].im + s.im *
            A[s_re_tmp].re);
          d = A[stemp_re_tmp].im;
          d1 = A[stemp_re_tmp].re;
          A[s_re_tmp].re = c * A[s_re_tmp].re - (s.re * A[stemp_re_tmp].re +
            s.im * A[stemp_re_tmp].im);
          A[s_re_tmp].im = c * A[s_re_tmp].im - (s.re * d - s.im * d1);
          A[stemp_re_tmp].re = stemp_re;
          A[stemp_re_tmp].im = stemp_im;
        }

        d = s.re;
        d1 = s.im;
        for (k = 0; k < 10; k++) {
          s_re_tmp = k + 10 * (jrow - 1);
          stemp_re_tmp = k + 10 * jrow;
          stemp_re = c * Z[stemp_re_tmp].re + (d * Z[s_re_tmp].re - d1 *
            Z[s_re_tmp].im);
          stemp_im = c * Z[stemp_re_tmp].im + (d * Z[s_re_tmp].im + d1 *
            Z[s_re_tmp].re);
          d2 = Z[stemp_re_tmp].re;
          Z[s_re_tmp].re = c * Z[s_re_tmp].re - (d * Z[stemp_re_tmp].re + d1 *
            Z[stemp_re_tmp].im);
          Z[s_re_tmp].im = c * Z[s_re_tmp].im - (d * Z[stemp_re_tmp].im - d1 *
            d2);
          Z[stemp_re_tmp].re = stemp_re;
          Z[stemp_re_tmp].im = stemp_im;
        }
      }
    }
  }
}

//
// File trailer for xzgghrd.cpp
//
// [EOF]
//
