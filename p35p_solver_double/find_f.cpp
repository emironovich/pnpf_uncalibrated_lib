/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find_f.cpp
 *
 * Code generation for function 'find_f'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "find_f.h"
#include "qr.h"

/* Function Definitions */
void find_f(const double F[72], double x, double y, double e, double fc_data[],
            int fc_size[2], double fs_data[], int fs_size[2], double *n)
{
  double F_eval[12];
  double mons[6];
  int i;
  double b_F_eval[12];
  double Q[9];
  int j;
  int i11;
  double Q_idx_0_tmp;
  int k;
  memset(&F_eval[0], 0, 12U * sizeof(double));
  mons[0] = x * x;
  mons[1] = x * y;
  mons[2] = y * y;
  mons[3] = x;
  mons[4] = y;
  mons[5] = 1.0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      i11 = i + (j << 2);
      Q_idx_0_tmp = F_eval[i11];
      for (k = 0; k < 6; k++) {
        Q_idx_0_tmp += F[i11 + 12 * k] * mons[k];
      }

      F_eval[i11] = Q_idx_0_tmp;
      b_F_eval[j + 3 * i] = Q_idx_0_tmp;
    }
  }

  qr(b_F_eval, Q, F_eval);
  if (std::abs(F_eval[4]) < e) {
    /* rankF = 1 */
    *n = 2.0;
    Q_idx_0_tmp = Q[3] / Q[5];
    fc_size[0] = 1;
    fc_size[1] = 2;
    fc_data[0] = Q_idx_0_tmp;
    fc_data[1] = Q_idx_0_tmp;
    fs_size[0] = 1;
    fs_size[1] = 2;
    fs_data[0] = Q[4] / Q[5];
    fs_data[1] = Q[7] / Q[8];
  } else {
    *n = 1.0;
    fc_size[0] = 1;
    fc_size[1] = 1;
    fc_data[0] = Q[6] / Q[8];
    fs_size[0] = 1;
    fs_size[1] = 1;
    fs_data[0] = Q[7] / Q[8];
  }
}

/* End of code generation (find_f.cpp) */
