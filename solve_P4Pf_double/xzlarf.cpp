/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlarf.cpp
 *
 * Code generation for function 'xzlarf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzlarf.h"

/* Function Definitions */
void b_xzlarf(int m, int n, int iv0, const creal_T tau, creal_T C_data[], int
              ic0, int ldc, creal_T work_data[])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  double c_re;
  double c_im;
  int i21;
  int jy;
  int ia;
  int ix;
  int exitg1;
  double temp_re;
  int i22;
  double temp_im;
  int ijA;
  double C_data_im;
  if ((tau.re != 0.0) || (tau.im != 0.0)) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && ((C_data[i - 2].re == 0.0) && (C_data[i - 2].im == 0.0)))
    {
      lastv--;
      i--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      i = ic0 + lastc * ldc;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if ((C_data[ia - 1].re != 0.0) || (C_data[ia - 1].im != 0.0)) {
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
    lastc = -1;
  }

  if (lastv > 0) {
    if (lastc + 1 != 0) {
      for (i = 0; i <= lastc; i++) {
        work_data[i].re = 0.0;
        work_data[i].im = 0.0;
      }

      i = 0;
      i21 = ic0 + ldc * lastc;
      for (jy = ic0; ldc < 0 ? jy >= i21 : jy <= i21; jy += ldc) {
        ix = iv0 - 1;
        c_re = 0.0;
        c_im = 0.0;
        i22 = (jy + lastv) - 1;
        for (ia = jy; ia <= i22; ia++) {
          c_re += C_data[ia - 1].re * C_data[ix].re + C_data[ia - 1].im *
            C_data[ix].im;
          c_im += C_data[ia - 1].re * C_data[ix].im - C_data[ia - 1].im *
            C_data[ix].re;
          ix++;
        }

        work_data[i].re += c_re - 0.0 * c_im;
        work_data[i].im += c_im + 0.0 * c_re;
        i++;
      }
    }

    c_re = -tau.re;
    c_im = -tau.im;
    if ((!(-tau.re == 0.0)) || (!(-tau.im == 0.0))) {
      i = ic0 - 1;
      jy = 0;
      for (ia = 0; ia <= lastc; ia++) {
        if ((work_data[jy].re != 0.0) || (work_data[jy].im != 0.0)) {
          temp_re = work_data[jy].re * c_re + work_data[jy].im * c_im;
          temp_im = work_data[jy].re * c_im - work_data[jy].im * c_re;
          ix = iv0;
          i21 = i + 1;
          i22 = lastv + i;
          for (ijA = i21; ijA <= i22; ijA++) {
            C_data_im = C_data[ix - 1].re * temp_im + C_data[ix - 1].im *
              temp_re;
            C_data[ijA - 1].re += C_data[ix - 1].re * temp_re - C_data[ix - 1].
              im * temp_im;
            C_data[ijA - 1].im += C_data_im;
            ix++;
          }
        }

        jy++;
        i += ldc;
      }
    }
  }
}

void xzlarf(int m, int n, int iv0, const creal_T tau, creal_T C_data[], int ic0,
            int ldc, creal_T work_data[])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  double c_re;
  double c_im;
  int ix;
  int i19;
  int rowright;
  int ia;
  int exitg1;
  double temp_re;
  int i20;
  double temp_im;
  int ijA;
  if ((tau.re != 0.0) || (tau.im != 0.0)) {
    lastv = n;
    i = iv0 + n;
    while ((lastv > 0) && ((C_data[i - 2].re == 0.0) && (C_data[i - 2].im == 0.0)))
    {
      lastv--;
      i--;
    }

    lastc = m;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = (ic0 + lastc) - 1;
      rowright = i + (lastv - 1) * ldc;
      do {
        exitg1 = 0;
        if (((ldc > 0) && (i <= rowright)) || ((ldc < 0) && (i >= rowright))) {
          if ((C_data[i - 1].re != 0.0) || (C_data[i - 1].im != 0.0)) {
            exitg1 = 1;
          } else {
            i += ldc;
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
      for (i = 0; i < lastc; i++) {
        work_data[i].re = 0.0;
        work_data[i].im = 0.0;
      }

      ix = iv0;
      i19 = ic0 + ldc * (lastv - 1);
      for (rowright = ic0; ldc < 0 ? rowright >= i19 : rowright <= i19; rowright
           += ldc) {
        c_re = C_data[ix - 1].re - 0.0 * C_data[ix - 1].im;
        c_im = C_data[ix - 1].im + 0.0 * C_data[ix - 1].re;
        i = 0;
        i20 = (rowright + lastc) - 1;
        for (ia = rowright; ia <= i20; ia++) {
          work_data[i].re += C_data[ia - 1].re * c_re - C_data[ia - 1].im * c_im;
          work_data[i].im += C_data[ia - 1].re * c_im + C_data[ia - 1].im * c_re;
          i++;
        }

        ix++;
      }
    }

    c_re = -tau.re;
    c_im = -tau.im;
    if ((!(-tau.re == 0.0)) || (!(-tau.im == 0.0))) {
      i = ic0 - 1;
      rowright = iv0 - 1;
      for (ia = 0; ia < lastv; ia++) {
        if ((C_data[rowright].re != 0.0) || (C_data[rowright].im != 0.0)) {
          temp_re = C_data[rowright].re * c_re + C_data[rowright].im * c_im;
          temp_im = C_data[rowright].re * c_im - C_data[rowright].im * c_re;
          ix = 0;
          i19 = i + 1;
          i20 = lastc + i;
          for (ijA = i19; ijA <= i20; ijA++) {
            C_data[ijA - 1].re += work_data[ix].re * temp_re - work_data[ix].im *
              temp_im;
            C_data[ijA - 1].im += work_data[ix].re * temp_im + work_data[ix].im *
              temp_re;
            ix++;
          }
        }

        rowright++;
        i += ldc;
      }
    }
  }
}

/* End of code generation (xzlarf.cpp) */
