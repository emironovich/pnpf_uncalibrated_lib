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
#include "xhseqr.h"
#include "xgehrd.h"
#include "anyNonFinite.h"

/* Function Definitions */
void schur(creal_T A_data[], int A_size[2], creal_T V_data[], int V_size[2])
{
  signed char unnamed_idx_0;
  signed char unnamed_idx_1;
  int m;
  int i7;
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
    for (i7 = 0; i7 < m; i7++) {
      V_data[i7].re = rtNaN;
      V_data[i7].im = 0.0;
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
          i7 = (i + unnamed_idx_0 * j) - 1;
          V_data[i7].re = 0.0;
          V_data[i7].im = 0.0;
        }

        istart++;
      }
    }
  } else {
    xgehrd(A_data, A_size);
    V_size[0] = A_size[0];
    V_size[1] = A_size[1];
    m = A_size[0] * A_size[1];
    if (0 <= m - 1) {
      memcpy(&V_data[0], &A_data[0], (unsigned int)(m * static_cast<int>(sizeof
               (creal_T))));
    }

    xhseqr(V_data, V_size);
  }
}

/* End of code generation (schur.cpp) */
