/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xztgevc.cpp
 *
 * Code generation for function 'xztgevc'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xztgevc.h"

/* Function Definitions */
void xztgevc(const creal_T A[100], creal_T V[100])
{
  double rworka[10];
  double anorm;
  int j;
  double xmx;
  int i;
  double d_re;
  double ascale;
  int je;
  int x_tmp;
  double temp;
  double salpha_re;
  double salpha_im;
  double acoeff;
  boolean_T lscalea;
  double z;
  boolean_T lscaleb;
  double scale;
  int jr;
  creal_T work1[10];
  double dmin;
  int i28;
  creal_T work2[10];
  double d_im;
  int i29;
  memset(&rworka[0], 0, 10U * sizeof(double));
  anorm = std::abs(A[0].re) + std::abs(A[0].im);
  for (j = 0; j < 9; j++) {
    for (i = 0; i <= j; i++) {
      rworka[j + 1] += std::abs(A[i + 10 * (j + 1)].re) + std::abs(A[i + 10 * (j
        + 1)].im);
    }

    d_re = rworka[j + 1] + (std::abs(A[(j + 10 * (j + 1)) + 1].re) + std::abs(A
      [(j + 10 * (j + 1)) + 1].im));
    if (d_re > anorm) {
      anorm = d_re;
    }
  }

  xmx = anorm;
  if (2.2250738585072014E-308 > anorm) {
    xmx = 2.2250738585072014E-308;
  }

  ascale = 1.0 / xmx;
  for (je = 0; je < 10; je++) {
    x_tmp = 10 * (9 - je);
    xmx = (std::abs(A[(x_tmp - je) + 9].re) + std::abs(A[(10 * (9 - je) - je) +
            9].im)) * ascale;
    if (1.0 > xmx) {
      xmx = 1.0;
    }

    temp = 1.0 / xmx;
    salpha_re = ascale * (temp * A[(10 * (9 - je) - je) + 9].re);
    salpha_im = ascale * (temp * A[(10 * (9 - je) - je) + 9].im);
    acoeff = temp * ascale;
    if ((temp >= 2.2250738585072014E-308) && (acoeff < 1.0020841800044864E-291))
    {
      lscalea = true;
    } else {
      lscalea = false;
    }

    z = std::abs(salpha_re) + std::abs(salpha_im);
    if ((z >= 2.2250738585072014E-308) && (z < 1.0020841800044864E-291)) {
      lscaleb = true;
    } else {
      lscaleb = false;
    }

    scale = 1.0;
    if (lscalea) {
      xmx = anorm;
      if (9.9792015476736E+290 < anorm) {
        xmx = 9.9792015476736E+290;
      }

      scale = 1.0020841800044864E-291 / temp * xmx;
    }

    if (lscaleb) {
      d_re = 1.0020841800044864E-291 / z;
      if (d_re > scale) {
        scale = d_re;
      }
    }

    if (lscalea || lscaleb) {
      xmx = acoeff;
      if (1.0 > acoeff) {
        xmx = 1.0;
      }

      if (z > xmx) {
        xmx = z;
      }

      d_re = 1.0 / (2.2250738585072014E-308 * xmx);
      if (d_re < scale) {
        scale = d_re;
      }

      if (lscalea) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      salpha_re *= scale;
      salpha_im *= scale;
    }

    for (jr = 0; jr < 10; jr++) {
      work1[jr].re = 0.0;
      work1[jr].im = 0.0;
    }

    work1[9 - je].re = 1.0;
    work1[9 - je].im = 0.0;
    dmin = 2.2204460492503131E-16 * acoeff * anorm;
    d_re = 2.2204460492503131E-16 * (std::abs(salpha_re) + std::abs(salpha_im));
    if (d_re > dmin) {
      dmin = d_re;
    }

    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }

    i28 = 8 - je;
    for (jr = 0; jr <= i28; jr++) {
      work1[jr].re = acoeff * A[jr + x_tmp].re;
      work1[jr].im = acoeff * A[jr + 10 * (9 - je)].im;
    }

    work1[9 - je].re = 1.0;
    work1[9 - je].im = 0.0;
    i28 = static_cast<int>(((1.0 + (-1.0 - ((10.0 + -static_cast<double>(je)) -
      1.0))) / -1.0));
    for (j = 0; j < i28; j++) {
      i = 8 - (je + j);
      d_re = acoeff * A[i + 10 * i].re - salpha_re;
      d_im = acoeff * A[i + 10 * i].im - salpha_im;
      if (std::abs(d_re) + std::abs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0;
      }

      xmx = std::abs(d_re) + std::abs(d_im);
      if (xmx < 1.0) {
        z = std::abs(work1[i].re) + std::abs(work1[i].im);
        if (z >= 4.49423283715579E+306 * xmx) {
          temp = 1.0 / z;
          i29 = 9 - je;
          for (jr = 0; jr <= i29; jr++) {
            work1[jr].re *= temp;
            work1[jr].im *= temp;
          }
        }
      }

      temp = -work1[i].re;
      if (d_im == 0.0) {
        if (-work1[i].im == 0.0) {
          work1[i].re = -work1[i].re / d_re;
          work1[i].im = 0.0;
        } else if (-work1[i].re == 0.0) {
          work1[i].re = 0.0;
          work1[i].im = -work1[i].im / d_re;
        } else {
          work1[i].re = -work1[i].re / d_re;
          work1[i].im = -work1[i].im / d_re;
        }
      } else if (d_re == 0.0) {
        if (-work1[i].re == 0.0) {
          work1[i].re = -work1[i].im / d_im;
          work1[i].im = 0.0;
        } else if (-work1[i].im == 0.0) {
          work1[i].re = 0.0;
          work1[i].im = -(temp / d_im);
        } else {
          work1[i].re = -work1[i].im / d_im;
          work1[i].im = -(temp / d_im);
        }
      } else {
        scale = std::abs(d_re);
        xmx = std::abs(d_im);
        if (scale > xmx) {
          scale = d_im / d_re;
          z = d_re + scale * d_im;
          work1[i].re = (-work1[i].re + scale * -work1[i].im) / z;
          work1[i].im = (-work1[i].im - scale * temp) / z;
        } else if (xmx == scale) {
          if (d_re > 0.0) {
            z = 0.5;
          } else {
            z = -0.5;
          }

          if (d_im > 0.0) {
            xmx = 0.5;
          } else {
            xmx = -0.5;
          }

          work1[i].re = (-work1[i].re * z + -work1[i].im * xmx) / scale;
          work1[i].im = (-work1[i].im * z - temp * xmx) / scale;
        } else {
          scale = d_re / d_im;
          z = d_im + scale * d_re;
          xmx = scale * -work1[i].im - (-work1[i].re);
          work1[i].re = (scale * -work1[i].re + -work1[i].im) / z;
          work1[i].im = xmx / z;
        }
      }

      if (i + 1 > 1) {
        xmx = std::abs(work1[i].re) + std::abs(work1[i].im);
        if (xmx > 1.0) {
          temp = 1.0 / xmx;
          if (acoeff * rworka[i] >= 4.49423283715579E+306 * temp) {
            i29 = 9 - je;
            for (jr = 0; jr <= i29; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        d_re = acoeff * work1[i].re;
        d_im = acoeff * work1[i].im;
        for (jr = 0; jr < i; jr++) {
          work1[jr].re += d_re * A[jr + 10 * i].re - d_im * A[jr + 10 * i].im;
          work1[jr].im += d_re * A[jr + 10 * i].im + d_im * A[jr + 10 * i].re;
        }
      }
    }

    for (jr = 0; jr < 10; jr++) {
      work2[jr].re = 0.0;
      work2[jr].im = 0.0;
    }

    i28 = 9 - je;
    for (i = 0; i <= i28; i++) {
      for (jr = 0; jr < 10; jr++) {
        work2[jr].re += V[jr + 10 * i].re * work1[i].re - V[jr + 10 * i].im *
          work1[i].im;
        work2[jr].im += V[jr + 10 * i].re * work1[i].im + V[jr + 10 * i].im *
          work1[i].re;
      }
    }

    xmx = std::abs(work2[0].re) + std::abs(work2[0].im);
    for (jr = 0; jr < 9; jr++) {
      d_re = std::abs(work2[jr + 1].re) + std::abs(work2[jr + 1].im);
      if (d_re > xmx) {
        xmx = d_re;
      }
    }

    if (xmx > 2.2250738585072014E-308) {
      temp = 1.0 / xmx;
      for (jr = 0; jr < 10; jr++) {
        i28 = jr + x_tmp;
        V[i28].re = temp * work2[jr].re;
        V[i28].im = temp * work2[jr].im;
      }
    } else {
      for (jr = 0; jr < 10; jr++) {
        i28 = jr + x_tmp;
        V[i28].re = 0.0;
        V[i28].im = 0.0;
      }
    }
  }
}

/* End of code generation (xztgevc.cpp) */
