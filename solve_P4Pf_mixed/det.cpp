/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * det.cpp
 *
 * Code generation for function 'det'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "det.h"

/* Function Definitions */
float det(const float x[9])
{
  float y;
  int i3;
  signed char ipiv[3];
  float b_x[9];
  int j;
  boolean_T isodd;
  int b_tmp;
  int jp1j;
  int n;
  int iy;
  int ix;
  float smax;
  int jA;
  float s;
  int i4;
  int ijA;
  for (i3 = 0; i3 < 9; i3++) {
    b_x[i3] = x[i3];
  }

  ipiv[0] = 1;
  ipiv[1] = 2;
  for (j = 0; j < 2; j++) {
    b_tmp = j << 2;
    jp1j = b_tmp + 2;
    n = 3 - j;
    iy = 0;
    ix = b_tmp;
    smax = std::abs(b_x[b_tmp]);
    for (jA = 2; jA <= n; jA++) {
      ix++;
      s = std::abs(b_x[ix]);
      if (s > smax) {
        iy = jA - 1;
        smax = s;
      }
    }

    if (b_x[b_tmp + iy] != 0.0F) {
      if (iy != 0) {
        iy += j;
        ipiv[j] = static_cast<signed char>((iy + 1));
        smax = b_x[j];
        b_x[j] = b_x[iy];
        b_x[iy] = smax;
        ix = j + 3;
        iy += 3;
        smax = b_x[ix];
        b_x[ix] = b_x[iy];
        b_x[iy] = smax;
        ix += 3;
        iy += 3;
        smax = b_x[ix];
        b_x[ix] = b_x[iy];
        b_x[iy] = smax;
      }

      i3 = b_tmp - j;
      for (iy = jp1j; iy <= i3 + 3; iy++) {
        b_x[iy - 1] /= b_x[b_tmp];
      }
    }

    n = 1 - j;
    iy = b_tmp + 3;
    jA = b_tmp + 4;
    for (jp1j = 0; jp1j <= n; jp1j++) {
      smax = b_x[iy];
      if (b_x[iy] != 0.0F) {
        ix = b_tmp + 1;
        i3 = jA + 1;
        i4 = (jA - j) + 2;
        for (ijA = i3; ijA <= i4; ijA++) {
          b_x[ijA - 1] += b_x[ix] * -smax;
          ix++;
        }
      }

      iy += 3;
      jA += 3;
    }
  }

  isodd = false;
  if (ipiv[0] > 1) {
    isodd = true;
  }

  y = b_x[0] * b_x[4] * b_x[8];
  if (ipiv[1] > 2) {
    isodd = !isodd;
  }

  if (isodd) {
    y = -y;
  }

  return y;
}

/* End of code generation (det.cpp) */
