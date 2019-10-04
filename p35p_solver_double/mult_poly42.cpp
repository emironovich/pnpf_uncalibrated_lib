//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: mult_poly42.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "mult_poly42.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : const double d[15]
//                const double a[6]
//                double p[28]
// Return Type  : void
//
void mult_poly42(const double d[15], const double a[6], double p[28])
{
  p[0] = a[0] * d[0];

  // x^6
  p[1] = a[0] * d[1] + a[1] * d[0];

  // x^5*y
  p[2] = (a[0] * d[2] + a[1] * d[1]) + a[2] * d[0];

  // x^4*y^2
  p[3] = (a[0] * d[3] + a[1] * d[2]) + a[2] * d[1];

  // x^3*y^3
  p[4] = (a[0] * d[4] + a[1] * d[3]) + a[2] * d[2];

  // x^2*y^4
  p[5] = a[1] * d[4] + a[2] * d[3];

  // x*y^5
  p[6] = a[2] * d[4];

  // y^6
  p[7] = a[3] * d[0] + a[0] * d[5];

  // x^5
  p[8] = ((a[3] * d[1] + a[4] * d[0]) + a[0] * d[6]) + a[1] * d[5];

  // x^4*y
  p[9] = (((a[3] * d[2] + a[4] * d[1]) + a[0] * d[7]) + a[1] * d[6]) + a[2] * d
    [5];

  // x^3*y^2
  p[10] = (((a[3] * d[3] + a[4] * d[2]) + a[0] * d[8]) + a[1] * d[7]) + a[2] *
    d[6];

  // x^2*y^3
  p[11] = ((a[3] * d[4] + a[4] * d[3]) + a[1] * d[8]) + a[2] * d[7];

  // x*y^4
  p[12] = a[4] * d[4] + a[2] * d[8];

  // y^5
  p[13] = (a[5] * d[0] + a[3] * d[5]) + a[0] * d[9];

  // x^4
  p[14] = (((a[5] * d[1] + a[3] * d[6]) + a[4] * d[5]) + a[0] * d[10]) + a[1] *
    d[9];

  // x^3*y
  p[15] = ((((a[5] * d[2] + a[3] * d[7]) + a[4] * d[6]) + a[0] * d[11]) + a[1] *
           d[10]) + a[2] * d[9];

  // x^2*y^2
  p[16] = (((a[5] * d[3] + a[3] * d[8]) + a[4] * d[7]) + a[1] * d[11]) + a[2] *
    d[10];

  // x*y^3
  p[17] = (a[5] * d[4] + a[4] * d[8]) + a[2] * d[11];

  // y^4
  p[18] = (a[5] * d[5] + a[0] * d[12]) + a[3] * d[9];

  // x^3
  p[19] = (((a[5] * d[6] + a[0] * d[13]) + a[1] * d[12]) + a[3] * d[10]) + a[4] *
    d[9];

  // x^2*y
  p[20] = (((a[5] * d[7] + a[1] * d[13]) + a[2] * d[12]) + a[3] * d[11]) + a[4] *
    d[10];

  // x*y^2
  p[21] = (a[5] * d[8] + a[2] * d[13]) + a[4] * d[11];

  // y^3
  p[22] = (a[0] * d[14] + a[5] * d[9]) + a[3] * d[12];

  // x^2
  p[23] = ((a[1] * d[14] + a[5] * d[10]) + a[3] * d[13]) + a[4] * d[12];

  // x*y
  p[24] = (a[2] * d[14] + a[5] * d[11]) + a[4] * d[13];

  // y^2
  p[25] = a[3] * d[14] + a[5] * d[12];

  // x
  p[26] = a[4] * d[14] + a[5] * d[13];

  // y
  p[27] = a[5] * d[14];

  // 1
}

//
// File trailer for mult_poly42.cpp
//
// [EOF]
//
