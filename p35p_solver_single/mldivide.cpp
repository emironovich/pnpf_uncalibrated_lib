/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mldivide.cpp
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "mldivide.h"
#include "xgetrf.h"

/* Function Definitions */
void mldivide(const float A[400], float B[200])
{
  float b_A[400];
  int ipiv[20];
  int info;
  int jBcol;
  int temp_tmp;
  float temp;
  int kAcol;
  int i9;
  int i10;
  int i;
  int i11;
  memcpy(&b_A[0], &A[0], 400U * sizeof(float));
  xgetrf(b_A, ipiv, &info);
  for (info = 0; info < 19; info++) {
    if (ipiv[info] != info + 1) {
      for (jBcol = 0; jBcol < 10; jBcol++) {
        temp_tmp = info + 20 * jBcol;
        temp = B[temp_tmp];
        i9 = (ipiv[info] + 20 * jBcol) - 1;
        B[temp_tmp] = B[i9];
        B[i9] = temp;
      }
    }
  }

  for (info = 0; info < 10; info++) {
    jBcol = 20 * info;
    for (temp_tmp = 0; temp_tmp < 20; temp_tmp++) {
      kAcol = 20 * temp_tmp;
      i9 = temp_tmp + jBcol;
      if (B[i9] != 0.0F) {
        i10 = temp_tmp + 2;
        for (i = i10; i < 21; i++) {
          i11 = (i + jBcol) - 1;
          B[i11] -= B[i9] * b_A[(i + kAcol) - 1];
        }
      }
    }
  }

  for (info = 0; info < 10; info++) {
    jBcol = 20 * info;
    for (temp_tmp = 19; temp_tmp >= 0; temp_tmp--) {
      kAcol = 20 * temp_tmp;
      i9 = temp_tmp + jBcol;
      temp = B[i9];
      if (temp != 0.0F) {
        B[i9] = temp / b_A[temp_tmp + kAcol];
        for (i = 0; i < temp_tmp; i++) {
          i10 = i + jBcol;
          B[i10] -= B[i9] * b_A[i + kAcol];
        }
      }
    }
  }
}

/* End of code generation (mldivide.cpp) */
