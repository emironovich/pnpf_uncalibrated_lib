//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: equations_for_groebner.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "equations_for_groebner.h"
#include "mult_poly22.h"
#include "mult_poly42.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : const double F[72]
//                double G[112]
// Return Type  : void
//
void equations_for_groebner(const double F[72], double G[112])
{
  int i;
  double b_F[6];
  double c_F[6];
  double dv[15];
  int M_tmp;
  double M[24];
  double dv1[15];
  int b_M_tmp;
  double dv2[28];
  double dv3[28];
  double dv4[28];
  int i1;
  double b_M[54];
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    M[M_tmp] = F[12 * i + 2];
    M[M_tmp + 2] = F[12 * i + 10];
    M[M_tmp + 1] = F[12 * i + 3];
    b_M_tmp = 12 * i + 11;
    M[M_tmp + 3] = F[b_M_tmp];
    b_F[i] = F[12 * i + 6];
    c_F[i] = F[b_M_tmp];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 10];
    c_F[i] = F[12 * i + 7];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 1];
  }

  mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 5];
  }

  mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 2];
    c_F[i] = F[12 * i + 7];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 6];
    c_F[i] = F[12 * i + 3];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 9];
  }

  mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[i << 2] = (dv2[i] - dv3[i]) + dv4[i];
  }

  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * i;
      b_M_tmp = 3 * i1 + 9 * i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 2];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    M[M_tmp] = b_M[9 * i + 1];
    M[M_tmp + 2] = b_M[9 * i + 7];
    M[M_tmp + 1] = b_M[9 * i + 2];
    b_M_tmp = 9 * i + 8;
    M[M_tmp + 3] = b_M[b_M_tmp];
    b_F[i] = b_M[9 * i + 4];
    c_F[i] = b_M[b_M_tmp];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 7];
    c_F[i] = b_M[9 * i + 5];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i];
  }

  mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 3];
  }

  mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 1];
    c_F[i] = b_M[9 * i + 5];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 4];
    c_F[i] = b_M[9 * i + 2];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 6];
  }

  mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[(i << 2) + 1] = (dv2[i] - dv3[i]) + dv4[i];
  }

  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * i;
      b_M_tmp = 3 * i1 + 9 * i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 1];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    M[M_tmp] = b_M[9 * i + 1];
    M[M_tmp + 2] = b_M[9 * i + 7];
    M[M_tmp + 1] = b_M[9 * i + 2];
    b_M_tmp = 9 * i + 8;
    M[M_tmp + 3] = b_M[b_M_tmp];
    b_F[i] = b_M[9 * i + 4];
    c_F[i] = b_M[b_M_tmp];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 7];
    c_F[i] = b_M[9 * i + 5];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i];
  }

  mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 3];
  }

  mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 1];
    c_F[i] = b_M[9 * i + 5];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 4];
    c_F[i] = b_M[9 * i + 2];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 6];
  }

  mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[(i << 2) + 2] = (dv2[i] - dv3[i]) + dv4[i];
  }

  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    M[M_tmp] = F[12 * i + 1];
    M[M_tmp + 2] = F[12 * i + 9];
    M[M_tmp + 1] = F[12 * i + 2];
    b_M_tmp = 12 * i + 10;
    M[M_tmp + 3] = F[b_M_tmp];
    b_F[i] = F[12 * i + 5];
    c_F[i] = F[b_M_tmp];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 9];
    c_F[i] = F[12 * i + 6];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i];
  }

  mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 4];
  }

  mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 1];
    c_F[i] = F[12 * i + 6];
  }

  mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 5];
    c_F[i] = F[12 * i + 2];
  }

  mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 8];
  }

  mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[(i << 2) + 3] = (dv2[i] - dv3[i]) + dv4[i];
  }
}

//
// File trailer for equations_for_groebner.cpp
//
// [EOF]
//
