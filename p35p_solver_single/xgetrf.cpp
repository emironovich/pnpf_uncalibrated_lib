/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgetrf.cpp
 *
 * Code generation for function 'xgetrf'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xgetrf.h"

/* Function Definitions */
void xgetrf(float A[400], int ipiv[20], int *info)
{
  int i12;
  int j;
  int b;
  int jj;
  int jp1j;
  int n;
  int iy;
  int ix;
  float smax;
  int jA;
  float s;
  int ijA;
  for (i12 = 0; i12 < 20; i12++) {
    ipiv[i12] = 1 + i12;
  }

  *info = 0;
  for (j = 0; j < 19; j++) {
    b = j * 21;
    jj = j * 21;
    jp1j = b + 2;
    n = 20 - j;
    iy = 0;
    ix = b;
    smax = std::abs(A[b]);
    for (jA = 2; jA <= n; jA++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        iy = jA - 1;
        smax = s;
      }
    }

    if (A[jj + iy] != 0.0F) {
      if (iy != 0) {
        iy += j;
        ipiv[j] = iy + 1;
        ix = j;
        for (jA = 0; jA < 20; jA++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 20;
          iy += 20;
        }
      }

      i12 = jj - j;
      for (iy = jp1j; iy <= i12 + 20; iy++) {
        A[iy - 1] /= A[jj];
      }
    } else {
      *info = j + 1;
    }

    n = 18 - j;
    iy = b + 20;
    jA = jj + 21;
    for (b = 0; b <= n; b++) {
      smax = A[iy];
      if (A[iy] != 0.0F) {
        ix = jj + 1;
        i12 = jA + 1;
        jp1j = (jA - j) + 19;
        for (ijA = i12; ijA <= jp1j; ijA++) {
          A[ijA - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 20;
      jA += 20;
    }
  }

  if ((*info == 0) && (!(A[399] != 0.0F))) {
    *info = 20;
  }
}

/* End of code generation (xgetrf.cpp) */
