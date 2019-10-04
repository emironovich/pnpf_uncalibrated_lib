//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xztgevc.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xztgevc.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Function Definitions

//
// Arguments    : const creal_T A[100]
//                creal_T V[100]
// Return Type  : void
//
void xztgevc(const creal_T A[100], creal_T V[100])
{
  double rworka[10];
  double anorm;
  int j;
  double xmx;
  int i;
  int re_tmp;
  double ascale;
  double d_re;
  int je;
  int x_tmp_tmp;
  double temp;
  double salpha_re;
  double salpha_im;
  double acoeff;
  boolean_T lscalea;
  double z;
  boolean_T lscaleb;
  double scale;
  creal_T work1[10];
  double dmin;
  int b_i;
  int jr;
  creal_T work2[10];
  double d_im;
  double brm;
  std::memset(&rworka[0], 0, 10U * sizeof(double));
  anorm = std::abs(A[0].re) + std::abs(A[0].im);
  for (j = 0; j < 9; j++) {
    for (i = 0; i <= j; i++) {
      re_tmp = i + 10 * (j + 1);
      rworka[j + 1] += std::abs(A[re_tmp].re) + std::abs(A[re_tmp].im);
    }

    i = (j + 10 * (j + 1)) + 1;
    d_re = rworka[j + 1] + (std::abs(A[i].re) + std::abs(A[i].im));
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
    x_tmp_tmp = 10 * (9 - je);
    i = (x_tmp_tmp - je) + 9;
    xmx = (std::abs(A[i].re) + std::abs(A[i].im)) * ascale;
    if (1.0 > xmx) {
      xmx = 1.0;
    }

    temp = 1.0 / xmx;
    salpha_re = ascale * (temp * A[i].re);
    salpha_im = ascale * (temp * A[i].im);
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

    std::memset(&work1[0], 0, 10U * sizeof(creal_T));
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

    b_i = 8 - je;
    for (jr = 0; jr <= b_i; jr++) {
      i = jr + x_tmp_tmp;
      work1[jr].re = acoeff * A[i].re;
      work1[jr].im = acoeff * A[i].im;
    }

    work1[9 - je].re = 1.0;
    work1[9 - je].im = 0.0;
    b_i = static_cast<int>((((-1.0 - ((-static_cast<double>(je) + 10.0) - 1.0))
      + 1.0) / -1.0));
    for (j = 0; j < b_i; j++) {
      re_tmp = 8 - (je + j);
      i = re_tmp + 10 * re_tmp;
      d_re = acoeff * A[i].re - salpha_re;
      d_im = acoeff * A[i].im - salpha_im;
      if (std::abs(d_re) + std::abs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0;
      }

      brm = std::abs(d_re);
      scale = std::abs(d_im);
      xmx = brm + scale;
      if (xmx < 1.0) {
        z = std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im);
        if (z >= 4.49423283715579E+306 * xmx) {
          temp = 1.0 / z;
          i = 9 - je;
          for (jr = 0; jr <= i; jr++) {
            work1[jr].re *= temp;
            work1[jr].im *= temp;
          }
        }
      }

      if (d_im == 0.0) {
        if (-work1[re_tmp].im == 0.0) {
          scale = -work1[re_tmp].re / d_re;
          xmx = 0.0;
        } else if (-work1[re_tmp].re == 0.0) {
          scale = 0.0;
          xmx = -work1[re_tmp].im / d_re;
        } else {
          scale = -work1[re_tmp].re / d_re;
          xmx = -work1[re_tmp].im / d_re;
        }
      } else if (d_re == 0.0) {
        if (-work1[re_tmp].re == 0.0) {
          scale = -work1[re_tmp].im / d_im;
          xmx = 0.0;
        } else if (-work1[re_tmp].im == 0.0) {
          scale = 0.0;
          xmx = -(-work1[re_tmp].re / d_im);
        } else {
          scale = -work1[re_tmp].im / d_im;
          xmx = -(-work1[re_tmp].re / d_im);
        }
      } else if (brm > scale) {
        z = d_im / d_re;
        xmx = d_re + z * d_im;
        scale = (-work1[re_tmp].re + z * -work1[re_tmp].im) / xmx;
        xmx = (-work1[re_tmp].im - z * -work1[re_tmp].re) / xmx;
      } else if (scale == brm) {
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

        scale = (-work1[re_tmp].re * z + -work1[re_tmp].im * xmx) / brm;
        xmx = (-work1[re_tmp].im * z - -work1[re_tmp].re * xmx) / brm;
      } else {
        z = d_re / d_im;
        xmx = d_im + z * d_re;
        scale = (z * -work1[re_tmp].re + -work1[re_tmp].im) / xmx;
        xmx = (z * -work1[re_tmp].im - (-work1[re_tmp].re)) / xmx;
      }

      work1[re_tmp].re = scale;
      work1[re_tmp].im = xmx;
      if (re_tmp + 1 > 1) {
        if (std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im) > 1.0) {
          temp = 1.0 / (std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im));
          if (acoeff * rworka[re_tmp] >= 4.49423283715579E+306 * temp) {
            i = 9 - je;
            for (jr = 0; jr <= i; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        d_re = acoeff * work1[re_tmp].re;
        d_im = acoeff * work1[re_tmp].im;
        for (jr = 0; jr < re_tmp; jr++) {
          i = jr + 10 * re_tmp;
          work1[jr].re += d_re * A[i].re - d_im * A[i].im;
          work1[jr].im += d_re * A[i].im + d_im * A[i].re;
        }
      }
    }

    std::memset(&work2[0], 0, 10U * sizeof(creal_T));
    b_i = 9 - je;
    for (i = 0; i <= b_i; i++) {
      xmx = work1[i].re;
      z = work1[i].im;
      for (jr = 0; jr < 10; jr++) {
        re_tmp = jr + 10 * i;
        work2[jr].re += V[re_tmp].re * xmx - V[re_tmp].im * z;
        work2[jr].im += V[re_tmp].re * z + V[re_tmp].im * xmx;
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
        b_i = jr + x_tmp_tmp;
        V[b_i].re = temp * work2[jr].re;
        V[b_i].im = temp * work2[jr].im;
      }
    } else {
      std::memset(&V[x_tmp_tmp], 0, 10U * sizeof(creal_T));
    }
  }
}

//
// File trailer for xztgevc.cpp
//
// [EOF]
//
