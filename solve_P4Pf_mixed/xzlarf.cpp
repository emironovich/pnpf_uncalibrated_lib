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
void b_xzlarf(int m, int n, int iv0, const creal32_T tau, creal32_T C_data[],
              int ic0, int ldc, creal32_T work_data[])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  float c_re;
  float c_im;
  int i15;
  int jy;
  int ia;
  int ix;
  int exitg1;
  float temp_re;
  int i16;
  float temp_im;
  int ijA;
  float C_data_im;
  if ((tau.re != 0.0F) || (tau.im != 0.0F)) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && ((C_data[i - 2].re == 0.0F) && (C_data[i - 2].im ==
             0.0F))) {
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
          if ((C_data[ia - 1].re != 0.0F) || (C_data[ia - 1].im != 0.0F)) {
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
        work_data[i].re = 0.0F;
        work_data[i].im = 0.0F;
      }

      i = 0;
      i15 = ic0 + ldc * lastc;
      for (jy = ic0; ldc < 0 ? jy >= i15 : jy <= i15; jy += ldc) {
        ix = iv0 - 1;
        c_re = 0.0F;
        c_im = 0.0F;
        i16 = (jy + lastv) - 1;
        for (ia = jy; ia <= i16; ia++) {
          c_re += C_data[ia - 1].re * C_data[ix].re + C_data[ia - 1].im *
            C_data[ix].im;
          c_im += C_data[ia - 1].re * C_data[ix].im - C_data[ia - 1].im *
            C_data[ix].re;
          ix++;
        }

        work_data[i].re += c_re - 0.0F * c_im;
        work_data[i].im += c_im + 0.0F * c_re;
        i++;
      }
    }

    c_re = -tau.re;
    c_im = -tau.im;
    if ((!(-tau.re == 0.0F)) || (!(-tau.im == 0.0F))) {
      i = ic0 - 1;
      jy = 0;
      for (ia = 0; ia <= lastc; ia++) {
        if ((work_data[jy].re != 0.0F) || (work_data[jy].im != 0.0F)) {
          temp_re = work_data[jy].re * c_re + work_data[jy].im * c_im;
          temp_im = work_data[jy].re * c_im - work_data[jy].im * c_re;
          ix = iv0;
          i15 = i + 1;
          i16 = lastv + i;
          for (ijA = i15; ijA <= i16; ijA++) {
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

void xzlarf(int m, int n, int iv0, const creal32_T tau, creal32_T C_data[], int
            ic0, int ldc, creal32_T work_data[])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  float c_re;
  float c_im;
  int ix;
  int i13;
  int rowright;
  int ia;
  int exitg1;
  float temp_re;
  int i14;
  float temp_im;
  int ijA;
  if ((tau.re != 0.0F) || (tau.im != 0.0F)) {
    lastv = n;
    i = iv0 + n;
    while ((lastv > 0) && ((C_data[i - 2].re == 0.0F) && (C_data[i - 2].im ==
             0.0F))) {
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
          if ((C_data[i - 1].re != 0.0F) || (C_data[i - 1].im != 0.0F)) {
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
        work_data[i].re = 0.0F;
        work_data[i].im = 0.0F;
      }

      ix = iv0;
      i13 = ic0 + ldc * (lastv - 1);
      for (rowright = ic0; ldc < 0 ? rowright >= i13 : rowright <= i13; rowright
           += ldc) {
        c_re = C_data[ix - 1].re - 0.0F * C_data[ix - 1].im;
        c_im = C_data[ix - 1].im + 0.0F * C_data[ix - 1].re;
        i = 0;
        i14 = (rowright + lastc) - 1;
        for (ia = rowright; ia <= i14; ia++) {
          work_data[i].re += C_data[ia - 1].re * c_re - C_data[ia - 1].im * c_im;
          work_data[i].im += C_data[ia - 1].re * c_im + C_data[ia - 1].im * c_re;
          i++;
        }

        ix++;
      }
    }

    c_re = -tau.re;
    c_im = -tau.im;
    if ((!(-tau.re == 0.0F)) || (!(-tau.im == 0.0F))) {
      i = ic0 - 1;
      rowright = iv0 - 1;
      for (ia = 0; ia < lastv; ia++) {
        if ((C_data[rowright].re != 0.0F) || (C_data[rowright].im != 0.0F)) {
          temp_re = C_data[rowright].re * c_re + C_data[rowright].im * c_im;
          temp_im = C_data[rowright].re * c_im - C_data[rowright].im * c_re;
          ix = 0;
          i13 = i + 1;
          i14 = lastc + i;
          for (ijA = i13; ijA <= i14; ijA++) {
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
