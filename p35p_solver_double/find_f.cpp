//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: find_f.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "find_f.h"
#include "p35p_solver.h"
#include "qr.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Function Definitions

//
// Arguments    : const double F[72]
//                double x
//                double y
//                double e
//                double fc_data[]
//                int fc_size[2]
//                double fs_data[]
//                int fs_size[2]
//                double *n
// Return Type  : void
//
void find_f(const double F[72], double x, double y, double e, double fc_data[],
            int fc_size[2], double fs_data[], int fs_size[2], double *n)
{
  double F_eval[12];
  double mons[6];
  int i;
  double b_F_eval[12];
  double Q[9];
  int j;
  int b_i;
  double fc_tmp;
  int k;
  std::memset(&F_eval[0], 0, 12U * sizeof(double));
  mons[0] = x * x;
  mons[1] = x * y;
  mons[2] = y * y;
  mons[3] = x;
  mons[4] = y;
  mons[5] = 1.0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      b_i = i + (j << 2);
      fc_tmp = F_eval[b_i];
      for (k = 0; k < 6; k++) {
        fc_tmp += F[b_i + 12 * k] * mons[k];
      }

      F_eval[b_i] = fc_tmp;
      b_F_eval[j + 3 * i] = fc_tmp;
    }
  }

  qr(b_F_eval, Q, F_eval);
  if (std::abs(F_eval[4]) < e) {
    // rankF = 1
    *n = 2.0;
    fc_size[0] = 1;
    fc_size[1] = 2;
    fc_tmp = Q[3] / Q[5];
    fc_data[0] = fc_tmp;
    fc_data[1] = fc_tmp;
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

//
// File trailer for find_f.cpp
//
// [EOF]
//
