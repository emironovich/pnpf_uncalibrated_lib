/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * equations_for_groebner.cpp
 *
 * Code generation for function 'equations_for_groebner'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "equations_for_groebner.h"
#include "mult_poly42.h"
#include "mult_poly22.h"

/* Function Definitions */
void equations_for_groebner(const double F[72], double G[112])
{
  int i3;
  double b_F[6];
  double c_F[6];
  double dv6[15];
  int M_tmp;
  double M[24];
  double dv7[15];
  double b_M_tmp;
  double dv8[28];
  double dv9[28];
  double dv10[28];
  int i4;
  double b_M[54];
  int c_M_tmp;
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    M[M_tmp] = F[12 * i3 + 2];
    M[2 + M_tmp] = F[12 * i3 + 10];
    M[1 + M_tmp] = F[12 * i3 + 3];
    b_M_tmp = F[12 * i3 + 11];
    M[M_tmp + 3] = b_M_tmp;
    b_F[i3] = F[6 + 12 * i3];
    c_F[i3] = b_M_tmp;
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[10 + 12 * i3];
    c_F[i3] = F[7 + 12 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[1 + 12 * i3];
  }

  mult_poly42(dv6, b_F, dv8);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[M_tmp];
    c_F[i3] = M[3 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[2 + M_tmp];
    c_F[i3] = M[1 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[5 + 12 * i3];
  }

  mult_poly42(dv6, b_F, dv9);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[2 + 12 * i3];
    c_F[i3] = F[7 + 12 * i3];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[6 + 12 * i3];
    c_F[i3] = F[3 + 12 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[9 + 12 * i3];
  }

  mult_poly42(dv6, b_F, dv10);
  for (i3 = 0; i3 < 28; i3++) {
    G[i3 << 2] = (dv8[i3] - dv9[i3]) + dv10[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    for (i4 = 0; i4 < 3; i4++) {
      M_tmp = (i4 << 2) + 12 * i3;
      c_M_tmp = 3 * i4 + 9 * i3;
      b_M[c_M_tmp] = F[M_tmp];
      b_M[c_M_tmp + 1] = F[M_tmp + 2];
      b_M[c_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    M[M_tmp] = b_M[9 * i3 + 1];
    M[2 + M_tmp] = b_M[9 * i3 + 7];
    M[1 + M_tmp] = b_M[9 * i3 + 2];
    b_M_tmp = b_M[9 * i3 + 8];
    M[M_tmp + 3] = b_M_tmp;
    b_F[i3] = b_M[4 + 9 * i3];
    c_F[i3] = b_M_tmp;
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[7 + 9 * i3];
    c_F[i3] = b_M[5 + 9 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[9 * i3];
  }

  mult_poly42(dv6, b_F, dv8);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[M_tmp];
    c_F[i3] = M[3 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[2 + M_tmp];
    c_F[i3] = M[1 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[3 + 9 * i3];
  }

  mult_poly42(dv6, b_F, dv9);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[1 + 9 * i3];
    c_F[i3] = b_M[5 + 9 * i3];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[4 + 9 * i3];
    c_F[i3] = b_M[2 + 9 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[6 + 9 * i3];
  }

  mult_poly42(dv6, b_F, dv10);
  for (i3 = 0; i3 < 28; i3++) {
    G[1 + (i3 << 2)] = (dv8[i3] - dv9[i3]) + dv10[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    for (i4 = 0; i4 < 3; i4++) {
      M_tmp = (i4 << 2) + 12 * i3;
      c_M_tmp = 3 * i4 + 9 * i3;
      b_M[c_M_tmp] = F[M_tmp];
      b_M[c_M_tmp + 1] = F[M_tmp + 1];
      b_M[2 + c_M_tmp] = F[3 + M_tmp];
    }
  }

  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    M[M_tmp] = b_M[9 * i3 + 1];
    M[2 + M_tmp] = b_M[9 * i3 + 7];
    M[1 + M_tmp] = b_M[9 * i3 + 2];
    b_M_tmp = b_M[9 * i3 + 8];
    M[M_tmp + 3] = b_M_tmp;
    b_F[i3] = b_M[4 + 9 * i3];
    c_F[i3] = b_M_tmp;
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[7 + 9 * i3];
    c_F[i3] = b_M[5 + 9 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[9 * i3];
  }

  mult_poly42(dv6, b_F, dv8);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[M_tmp];
    c_F[i3] = M[3 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[2 + M_tmp];
    c_F[i3] = M[1 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[3 + 9 * i3];
  }

  mult_poly42(dv6, b_F, dv9);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[1 + 9 * i3];
    c_F[i3] = b_M[5 + 9 * i3];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[4 + 9 * i3];
    c_F[i3] = b_M[2 + 9 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = b_M[6 + 9 * i3];
  }

  mult_poly42(dv6, b_F, dv10);
  for (i3 = 0; i3 < 28; i3++) {
    G[2 + (i3 << 2)] = (dv8[i3] - dv9[i3]) + dv10[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    M[M_tmp] = F[12 * i3 + 1];
    M[2 + M_tmp] = F[12 * i3 + 9];
    M[1 + M_tmp] = F[12 * i3 + 2];
    b_M_tmp = F[12 * i3 + 10];
    M[M_tmp + 3] = b_M_tmp;
    b_F[i3] = F[5 + 12 * i3];
    c_F[i3] = b_M_tmp;
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[9 + 12 * i3];
    c_F[i3] = F[6 + 12 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[12 * i3];
  }

  mult_poly42(dv6, b_F, dv8);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[M_tmp];
    c_F[i3] = M[3 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    M_tmp = i3 << 2;
    b_F[i3] = M[2 + M_tmp];
    c_F[i3] = M[1 + M_tmp];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[4 + 12 * i3];
  }

  mult_poly42(dv6, b_F, dv9);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[1 + 12 * i3];
    c_F[i3] = F[6 + 12 * i3];
  }

  mult_poly22(b_F, c_F, dv6);
  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[5 + 12 * i3];
    c_F[i3] = F[2 + 12 * i3];
  }

  mult_poly22(b_F, c_F, dv7);
  for (i3 = 0; i3 < 15; i3++) {
    dv6[i3] -= dv7[i3];
  }

  for (i3 = 0; i3 < 6; i3++) {
    b_F[i3] = F[8 + 12 * i3];
  }

  mult_poly42(dv6, b_F, dv10);
  for (i3 = 0; i3 < 28; i3++) {
    G[3 + (i3 << 2)] = (dv8[i3] - dv9[i3]) + dv10[i3];
  }
}

/* End of code generation (equations_for_groebner.cpp) */
