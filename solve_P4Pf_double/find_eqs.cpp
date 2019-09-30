/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find_eqs.cpp
 *
 * Code generation for function 'find_eqs'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "find_eqs.h"

/* Function Definitions */
void find_eqs(const double ND[48], double eqs[40])
{
  double b_ND[10];
  double c_ND[10];
  double d_ND[10];
  int i4;
  double ND_tmp;
  double e_ND[10];
  double f_ND[10];
  double g_ND[10];

  /* ND = [N; D] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  b_ND[0] = ND[0] * ND[4];
  b_ND[1] = ND[12] * ND[16];
  b_ND[2] = ND[24] * ND[28];
  b_ND[3] = ND[0] * ND[16] + ND[12] * ND[4];
  b_ND[4] = ND[0] * ND[28] + ND[24] * ND[4];
  b_ND[5] = ND[12] * ND[28] + ND[24] * ND[16];
  b_ND[6] = ND[0] * ND[40] + ND[36] * ND[4];
  b_ND[7] = ND[12] * ND[40] + ND[36] * ND[16];
  b_ND[8] = ND[24] * ND[40] + ND[36] * ND[28];
  b_ND[9] = ND[36] * ND[40];
  c_ND[0] = ND[1] * ND[5];
  c_ND[1] = ND[13] * ND[17];
  c_ND[2] = ND[25] * ND[29];
  c_ND[3] = ND[1] * ND[17] + ND[13] * ND[5];
  c_ND[4] = ND[1] * ND[29] + ND[25] * ND[5];
  c_ND[5] = ND[13] * ND[29] + ND[25] * ND[17];
  c_ND[6] = ND[1] * ND[41] + ND[37] * ND[5];
  c_ND[7] = ND[13] * ND[41] + ND[37] * ND[17];
  c_ND[8] = ND[25] * ND[41] + ND[37] * ND[29];
  c_ND[9] = ND[37] * ND[41];
  d_ND[0] = ND[2] * ND[6];
  d_ND[1] = ND[14] * ND[18];
  d_ND[2] = ND[26] * ND[30];
  d_ND[3] = ND[2] * ND[18] + ND[14] * ND[6];
  d_ND[4] = ND[2] * ND[30] + ND[26] * ND[6];
  d_ND[5] = ND[14] * ND[30] + ND[26] * ND[18];
  d_ND[6] = ND[2] * ND[42] + ND[38] * ND[6];
  d_ND[7] = ND[14] * ND[42] + ND[38] * ND[18];
  d_ND[8] = ND[26] * ND[42] + ND[38] * ND[30];
  d_ND[9] = ND[38] * ND[42];
  for (i4 = 0; i4 < 10; i4++) {
    eqs[i4 << 2] = (b_ND[i4] + c_ND[i4]) + d_ND[i4];
  }

  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  b_ND[0] = ND[8] * ND[0];
  b_ND[1] = ND[20] * ND[12];
  b_ND[2] = ND[32] * ND[24];
  b_ND[3] = ND[8] * ND[12] + ND[20] * ND[0];
  b_ND[4] = ND[8] * ND[24] + ND[32] * ND[0];
  b_ND[5] = ND[20] * ND[24] + ND[32] * ND[12];
  b_ND[6] = ND[8] * ND[36] + ND[44] * ND[0];
  b_ND[7] = ND[20] * ND[36] + ND[44] * ND[12];
  b_ND[8] = ND[32] * ND[36] + ND[44] * ND[24];
  b_ND[9] = ND[44] * ND[36];
  c_ND[0] = ND[9] * ND[1];
  c_ND[1] = ND[21] * ND[13];
  c_ND[2] = ND[33] * ND[25];
  c_ND[3] = ND[9] * ND[13] + ND[21] * ND[1];
  c_ND[4] = ND[9] * ND[25] + ND[33] * ND[1];
  c_ND[5] = ND[21] * ND[25] + ND[33] * ND[13];
  c_ND[6] = ND[9] * ND[37] + ND[45] * ND[1];
  c_ND[7] = ND[21] * ND[37] + ND[45] * ND[13];
  c_ND[8] = ND[33] * ND[37] + ND[45] * ND[25];
  c_ND[9] = ND[45] * ND[37];
  d_ND[0] = ND[10] * ND[2];
  d_ND[1] = ND[22] * ND[14];
  d_ND[2] = ND[34] * ND[26];
  d_ND[3] = ND[10] * ND[14] + ND[22] * ND[2];
  d_ND[4] = ND[10] * ND[26] + ND[34] * ND[2];
  d_ND[5] = ND[22] * ND[26] + ND[34] * ND[14];
  d_ND[6] = ND[10] * ND[38] + ND[46] * ND[2];
  d_ND[7] = ND[22] * ND[38] + ND[46] * ND[14];
  d_ND[8] = ND[34] * ND[38] + ND[46] * ND[26];
  d_ND[9] = ND[46] * ND[38];
  for (i4 = 0; i4 < 10; i4++) {
    eqs[1 + (i4 << 2)] = (b_ND[i4] + c_ND[i4]) + d_ND[i4];
  }

  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  b_ND[0] = ND[8] * ND[4];
  b_ND[1] = ND[20] * ND[16];
  b_ND[2] = ND[32] * ND[28];
  b_ND[3] = ND[8] * ND[16] + ND[20] * ND[4];
  b_ND[4] = ND[8] * ND[28] + ND[32] * ND[4];
  b_ND[5] = ND[20] * ND[28] + ND[32] * ND[16];
  b_ND[6] = ND[8] * ND[40] + ND[44] * ND[4];
  b_ND[7] = ND[20] * ND[40] + ND[44] * ND[16];
  b_ND[8] = ND[32] * ND[40] + ND[44] * ND[28];
  b_ND[9] = ND[44] * ND[40];
  c_ND[0] = ND[9] * ND[5];
  c_ND[1] = ND[21] * ND[17];
  c_ND[2] = ND[33] * ND[29];
  c_ND[3] = ND[9] * ND[17] + ND[21] * ND[5];
  c_ND[4] = ND[9] * ND[29] + ND[33] * ND[5];
  c_ND[5] = ND[21] * ND[29] + ND[33] * ND[17];
  c_ND[6] = ND[9] * ND[41] + ND[45] * ND[5];
  c_ND[7] = ND[21] * ND[41] + ND[45] * ND[17];
  c_ND[8] = ND[33] * ND[41] + ND[45] * ND[29];
  c_ND[9] = ND[45] * ND[41];
  d_ND[0] = ND[10] * ND[6];
  d_ND[1] = ND[22] * ND[18];
  d_ND[2] = ND[34] * ND[30];
  d_ND[3] = ND[10] * ND[18] + ND[22] * ND[6];
  d_ND[4] = ND[10] * ND[30] + ND[34] * ND[6];
  d_ND[5] = ND[22] * ND[30] + ND[34] * ND[18];
  d_ND[6] = ND[10] * ND[42] + ND[46] * ND[6];
  d_ND[7] = ND[22] * ND[42] + ND[46] * ND[18];
  d_ND[8] = ND[34] * ND[42] + ND[46] * ND[30];
  d_ND[9] = ND[46] * ND[42];
  for (i4 = 0; i4 < 10; i4++) {
    eqs[2 + (i4 << 2)] = (b_ND[i4] + c_ND[i4]) + d_ND[i4];
  }

  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  /*  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)] */
  b_ND[0] = ND[0] * ND[0];
  b_ND[1] = ND[12] * ND[12];
  b_ND[2] = ND[24] * ND[24];
  ND_tmp = ND[0] * ND[12];
  b_ND[3] = ND_tmp + ND_tmp;
  ND_tmp = ND[0] * ND[24];
  b_ND[4] = ND_tmp + ND_tmp;
  ND_tmp = ND[12] * ND[24];
  b_ND[5] = ND_tmp + ND_tmp;
  ND_tmp = ND[0] * ND[36];
  b_ND[6] = ND_tmp + ND_tmp;
  ND_tmp = ND[12] * ND[36];
  b_ND[7] = ND_tmp + ND_tmp;
  ND_tmp = ND[24] * ND[36];
  b_ND[8] = ND_tmp + ND_tmp;
  b_ND[9] = ND[36] * ND[36];
  c_ND[0] = ND[1] * ND[1];
  c_ND[1] = ND[13] * ND[13];
  c_ND[2] = ND[25] * ND[25];
  ND_tmp = ND[1] * ND[13];
  c_ND[3] = ND_tmp + ND_tmp;
  ND_tmp = ND[1] * ND[25];
  c_ND[4] = ND_tmp + ND_tmp;
  ND_tmp = ND[13] * ND[25];
  c_ND[5] = ND_tmp + ND_tmp;
  ND_tmp = ND[1] * ND[37];
  c_ND[6] = ND_tmp + ND_tmp;
  ND_tmp = ND[13] * ND[37];
  c_ND[7] = ND_tmp + ND_tmp;
  ND_tmp = ND[25] * ND[37];
  c_ND[8] = ND_tmp + ND_tmp;
  c_ND[9] = ND[37] * ND[37];
  d_ND[0] = ND[2] * ND[2];
  d_ND[1] = ND[14] * ND[14];
  d_ND[2] = ND[26] * ND[26];
  ND_tmp = ND[2] * ND[14];
  d_ND[3] = ND_tmp + ND_tmp;
  ND_tmp = ND[2] * ND[26];
  d_ND[4] = ND_tmp + ND_tmp;
  ND_tmp = ND[14] * ND[26];
  d_ND[5] = ND_tmp + ND_tmp;
  ND_tmp = ND[2] * ND[38];
  d_ND[6] = ND_tmp + ND_tmp;
  ND_tmp = ND[14] * ND[38];
  d_ND[7] = ND_tmp + ND_tmp;
  ND_tmp = ND[26] * ND[38];
  d_ND[8] = ND_tmp + ND_tmp;
  d_ND[9] = ND[38] * ND[38];
  e_ND[0] = ND[4] * ND[4];
  e_ND[1] = ND[16] * ND[16];
  e_ND[2] = ND[28] * ND[28];
  ND_tmp = ND[4] * ND[16];
  e_ND[3] = ND_tmp + ND_tmp;
  ND_tmp = ND[4] * ND[28];
  e_ND[4] = ND_tmp + ND_tmp;
  ND_tmp = ND[16] * ND[28];
  e_ND[5] = ND_tmp + ND_tmp;
  ND_tmp = ND[4] * ND[40];
  e_ND[6] = ND_tmp + ND_tmp;
  ND_tmp = ND[16] * ND[40];
  e_ND[7] = ND_tmp + ND_tmp;
  ND_tmp = ND[28] * ND[40];
  e_ND[8] = ND_tmp + ND_tmp;
  e_ND[9] = ND[40] * ND[40];
  f_ND[0] = ND[5] * ND[5];
  f_ND[1] = ND[17] * ND[17];
  f_ND[2] = ND[29] * ND[29];
  ND_tmp = ND[5] * ND[17];
  f_ND[3] = ND_tmp + ND_tmp;
  ND_tmp = ND[5] * ND[29];
  f_ND[4] = ND_tmp + ND_tmp;
  ND_tmp = ND[17] * ND[29];
  f_ND[5] = ND_tmp + ND_tmp;
  ND_tmp = ND[5] * ND[41];
  f_ND[6] = ND_tmp + ND_tmp;
  ND_tmp = ND[17] * ND[41];
  f_ND[7] = ND_tmp + ND_tmp;
  ND_tmp = ND[29] * ND[41];
  f_ND[8] = ND_tmp + ND_tmp;
  f_ND[9] = ND[41] * ND[41];
  g_ND[0] = ND[6] * ND[6];
  g_ND[1] = ND[18] * ND[18];
  g_ND[2] = ND[30] * ND[30];
  ND_tmp = ND[6] * ND[18];
  g_ND[3] = ND_tmp + ND_tmp;
  ND_tmp = ND[6] * ND[30];
  g_ND[4] = ND_tmp + ND_tmp;
  ND_tmp = ND[18] * ND[30];
  g_ND[5] = ND_tmp + ND_tmp;
  ND_tmp = ND[6] * ND[42];
  g_ND[6] = ND_tmp + ND_tmp;
  ND_tmp = ND[18] * ND[42];
  g_ND[7] = ND_tmp + ND_tmp;
  ND_tmp = ND[30] * ND[42];
  g_ND[8] = ND_tmp + ND_tmp;
  g_ND[9] = ND[42] * ND[42];
  for (i4 = 0; i4 < 10; i4++) {
    eqs[3 + (i4 << 2)] = ((((b_ND[i4] + c_ND[i4]) + d_ND[i4]) - e_ND[i4]) -
                          f_ND[i4]) - g_ND[i4];
  }
}

/* End of code generation (find_eqs.cpp) */
