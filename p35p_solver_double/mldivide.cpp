//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: mldivide.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "mldivide.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Function Definitions

//
// Arguments    : const double A[400]
//                double B[200]
// Return Type  : void
//
void mldivide(const double A[400], double B[200])
{
  double b_A[400];
  int i;
  int j;
  signed char ipiv[20];
  int mmj_tmp;
  int b;
  int jA;
  int jj;
  int k;
  int jp1j;
  int iy;
  int ix;
  double smax;
  int i1;
  double s;
  std::memcpy(&b_A[0], &A[0], 400U * sizeof(double));
  for (i = 0; i < 20; i++) {
    ipiv[i] = static_cast<signed char>((i + 1));
  }

  for (j = 0; j < 19; j++) {
    mmj_tmp = 18 - j;
    b = j * 21;
    jj = j * 21;
    jp1j = b + 2;
    iy = 20 - j;
    jA = 0;
    ix = b;
    smax = std::abs(b_A[jj]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = std::abs(b_A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (b_A[jj + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = static_cast<signed char>((iy + 1));
        ix = j;
        for (k = 0; k < 20; k++) {
          smax = b_A[ix];
          b_A[ix] = b_A[iy];
          b_A[iy] = smax;
          ix += 20;
          iy += 20;
        }
      }

      i = (jj - j) + 20;
      for (ix = jp1j; ix <= i; ix++) {
        b_A[ix - 1] /= b_A[jj];
      }
    }

    iy = b + 20;
    jA = jj;
    for (b = 0; b <= mmj_tmp; b++) {
      smax = b_A[iy];
      if (b_A[iy] != 0.0) {
        ix = jj + 1;
        i = jA + 22;
        i1 = (jA - j) + 40;
        for (jp1j = i; jp1j <= i1; jp1j++) {
          b_A[jp1j - 1] += b_A[ix] * -smax;
          ix++;
        }
      }

      iy += 20;
      jA += 20;
    }

    if (ipiv[j] != j + 1) {
      for (iy = 0; iy < 10; iy++) {
        jA = j + 20 * iy;
        smax = B[jA];
        i = (ipiv[j] + 20 * iy) - 1;
        B[jA] = B[i];
        B[i] = smax;
      }
    }
  }

  for (j = 0; j < 10; j++) {
    jA = 20 * j;
    for (k = 0; k < 20; k++) {
      iy = 20 * k;
      i = k + jA;
      if (B[i] != 0.0) {
        i1 = k + 2;
        for (ix = i1; ix < 21; ix++) {
          b = (ix + jA) - 1;
          B[b] -= B[i] * b_A[(ix + iy) - 1];
        }
      }
    }
  }

  for (j = 0; j < 10; j++) {
    jA = 20 * j;
    for (k = 19; k >= 0; k--) {
      iy = 20 * k;
      i = k + jA;
      if (B[i] != 0.0) {
        B[i] /= b_A[k + iy];
        for (ix = 0; ix < k; ix++) {
          i1 = ix + jA;
          B[i1] -= B[i] * b_A[ix + iy];
        }
      }
    }
  }
}

//
// File trailer for mldivide.cpp
//
// [EOF]
//
