//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: mult_poly22.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "mult_poly22.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : const double a[6]
//                const double b[6]
//                double c[15]
// Return Type  : void
//
void mult_poly22(const double a[6], const double b[6], double c[15])
{
  c[0] = a[0] * b[0];

  // x^4
  c[1] = a[0] * b[1] + a[1] * b[0];

  // x^3*y
  c[2] = (a[0] * b[2] + a[1] * b[1]) + a[2] * b[0];

  // x^2*y^2
  c[3] = a[1] * b[2] + a[2] * b[1];

  // x*y^3
  c[4] = a[2] * b[2];

  // y^4
  c[5] = a[0] * b[3] + a[3] * b[0];

  // x^3
  c[6] = ((a[0] * b[4] + a[1] * b[3]) + a[3] * b[1]) + a[4] * b[0];

  // x^2*y
  c[7] = ((a[1] * b[4] + a[2] * b[3]) + a[3] * b[2]) + a[4] * b[1];

  // x*y^2
  c[8] = a[2] * b[4] + a[4] * b[2];

  // y^3
  c[9] = (a[0] * b[5] + a[5] * b[0]) + a[3] * b[3];

  // x^2
  c[10] = ((a[1] * b[5] + a[5] * b[1]) + a[3] * b[4]) + a[4] * b[3];

  // x*y
  c[11] = (a[2] * b[5] + a[5] * b[2]) + a[4] * b[4];

  // y^2
  c[12] = a[3] * b[5] + a[5] * b[3];

  // x
  c[13] = a[4] * b[5] + a[5] * b[4];

  // y
  c[14] = a[5] * b[5];

  // 1
}

//
// File trailer for mult_poly22.cpp
//
// [EOF]
//
