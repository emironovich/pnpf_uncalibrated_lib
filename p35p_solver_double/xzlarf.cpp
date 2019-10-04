//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzlarf.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xzlarf.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include <cstring>

// Function Definitions

//
// Arguments    : int m
//                int n
//                int iv0
//                double tau
//                double C[100]
//                int ic0
//                double work[10]
// Return Type  : void
//
void xzlarf(int m, int n, int iv0, double tau, double C[100], int ic0, double
            work[10])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int jy;
  int b_i;
  int j;
  int ia;
  int ix;
  int exitg1;
  double c;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      i = ic0 + lastc * 10;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0) {
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
      if (0 <= lastc) {
        std::memset(&work[0], 0, (lastc + 1) * sizeof(double));
      }

      i = 0;
      b_i = ic0 + 10 * lastc;
      for (jy = ic0; jy <= b_i; jy += 10) {
        ix = iv0;
        c = 0.0;
        j = (jy + lastv) - 1;
        for (ia = jy; ia <= j; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[i] += c;
        i++;
      }
    }

    if (!(-tau == 0.0)) {
      i = ic0;
      jy = 0;
      for (j = 0; j <= lastc; j++) {
        if (work[jy] != 0.0) {
          c = work[jy] * -tau;
          ix = iv0;
          b_i = lastv + i;
          for (ia = i; ia < b_i; ia++) {
            C[ia - 1] += C[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += 10;
      }
    }
  }
}

//
// File trailer for xzlarf.cpp
//
// [EOF]
//
