/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * schur.cpp
 *
 * Code generation for function 'schur'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "schur.h"
#include "xzhseqr.h"
#include "xgehrd.h"
#include "anyNonFinite.h"

/* Function Definitions */
void schur(const creal32_T A_data[], const int A_size[2], creal32_T V_data[],
           int V_size[2])
{
  signed char unnamed_idx_0;
  signed char unnamed_idx_1;
  int m;
  int i6;
  int istart;
  int jend;
  int j;
  int i;
  if (anyNonFinite(A_data, A_size)) {
    unnamed_idx_0 = static_cast<signed char>(A_size[0]);
    unnamed_idx_1 = static_cast<signed char>(A_size[1]);
    V_size[0] = unnamed_idx_0;
    V_size[1] = unnamed_idx_1;
    m = unnamed_idx_0 * unnamed_idx_1;
    for (i6 = 0; i6 < m; i6++) {
      V_data[i6].re = rtNaNF;
      V_data[i6].im = 0.0F;
    }

    m = unnamed_idx_0;
    if (1 < unnamed_idx_0) {
      istart = 2;
      if (unnamed_idx_0 - 2 < unnamed_idx_1 - 1) {
        jend = unnamed_idx_0 - 1;
      } else {
        jend = unnamed_idx_1;
      }

      for (j = 0; j < jend; j++) {
        for (i = istart; i <= m; i++) {
          i6 = (i + unnamed_idx_0 * j) - 1;
          V_data[i6].re = 0.0F;
          V_data[i6].im = 0.0F;
        }

        istart++;
      }
    }
  } else {
    V_size[0] = A_size[0];
    V_size[1] = A_size[1];
    m = A_size[0] * A_size[1];
    if (0 <= m - 1) {
      memcpy(&V_data[0], &A_data[0], (unsigned int)(m * static_cast<int>(sizeof
               (creal32_T))));
    }

    xgehrd(V_data, V_size);
    eml_zlahqr(V_data, V_size);
    m = V_size[0];
    if (3 < V_size[0]) {
      istart = 4;
      if (V_size[0] - 4 < V_size[1] - 1) {
        jend = V_size[0] - 3;
      } else {
        jend = V_size[1];
      }

      for (j = 0; j < jend; j++) {
        for (i = istart; i <= m; i++) {
          i6 = (i + V_size[0] * j) - 1;
          V_data[i6].re = 0.0F;
          V_data[i6].im = 0.0F;
        }

        istart++;
      }
    }
  }
}

/* End of code generation (schur.cpp) */
