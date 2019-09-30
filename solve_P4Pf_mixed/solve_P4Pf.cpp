/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * solve_P4Pf.cpp
 *
 * Code generation for function 'solve_P4Pf'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "svd.h"
#include "xgeqp3.h"
#include "det.h"
#include "solve_3Q3.h"
#include "qr.h"

/* Function Declarations */
static void find_D(const float X[16], const float u[4], const float v[4], const
                   float ns[32], float e, float D[16]);

/* Function Definitions */
static void find_D(const float X[16], const float u[4], const float v[4], const
                   float ns[32], float e, float D[16])
{
  int i;
  int j;
  float smax;
  int X_tmp;
  int offset;
  int jj;
  int iy;
  int jp1j;
  int n;
  signed char ipiv[4];
  int ix;
  int jA;
  float B[16];
  float s;
  int i1;
  for (i = 0; i < 4; i++) {
    if (std::abs(u[i]) < e) {
      smax = v[i];
      offset = 4;
    } else {
      smax = u[i];
      offset = 0;
    }

    for (j = 0; j < 4; j++) {
      iy = i << 2;
      jA = i + (j << 2);
      B[jA] = smax * X[j + iy];
      X_tmp = offset + (j << 3);
      D[jA] = ((X[iy] * ns[X_tmp] + X[1 + iy] * ns[X_tmp + 1]) + X[2 + iy] *
               ns[X_tmp + 2]) + X[3 + iy] * ns[X_tmp + 3];
    }

    ipiv[i] = static_cast<signed char>((1 + i));
  }

  for (j = 0; j < 3; j++) {
    X_tmp = j * 5;
    jj = j * 5;
    jp1j = X_tmp + 2;
    n = 4 - j;
    offset = 0;
    ix = X_tmp;
    smax = std::abs(B[X_tmp]);
    for (jA = 2; jA <= n; jA++) {
      ix++;
      s = std::abs(B[ix]);
      if (s > smax) {
        offset = jA - 1;
        smax = s;
      }
    }

    if (B[jj + offset] != 0.0F) {
      if (offset != 0) {
        iy = j + offset;
        ipiv[j] = static_cast<signed char>((iy + 1));
        smax = B[j];
        B[j] = B[iy];
        B[iy] = smax;
        ix = j + 4;
        iy += 4;
        smax = B[ix];
        B[ix] = B[iy];
        B[iy] = smax;
        ix += 4;
        iy += 4;
        smax = B[ix];
        B[ix] = B[iy];
        B[iy] = smax;
        ix += 4;
        iy += 4;
        smax = B[ix];
        B[ix] = B[iy];
        B[iy] = smax;
      }

      i1 = jj - j;
      for (i = jp1j; i <= i1 + 4; i++) {
        B[i - 1] /= B[jj];
      }
    }

    n = 2 - j;
    offset = X_tmp + 4;
    jA = jj + 5;
    for (X_tmp = 0; X_tmp <= n; X_tmp++) {
      smax = B[offset];
      if (B[offset] != 0.0F) {
        ix = jj + 1;
        i1 = jA + 1;
        iy = (jA - j) + 3;
        for (i = i1; i <= iy; i++) {
          B[i - 1] += B[ix] * -smax;
          ix++;
        }
      }

      offset += 4;
      jA += 4;
    }

    if (ipiv[j] != j + 1) {
      smax = D[j];
      offset = ipiv[j] - 1;
      D[j] = D[offset];
      D[offset] = smax;
      smax = D[j + 4];
      offset = ipiv[j] + 3;
      D[j + 4] = D[offset];
      D[offset] = smax;
      smax = D[j + 8];
      offset = ipiv[j] + 7;
      D[j + 8] = D[offset];
      D[offset] = smax;
      smax = D[j + 12];
      offset = ipiv[j] + 11;
      D[j + 12] = D[offset];
      D[offset] = smax;
    }
  }

  for (j = 0; j < 4; j++) {
    iy = j << 2;
    if (D[iy] != 0.0F) {
      for (i = 2; i < 5; i++) {
        offset = (i + iy) - 1;
        D[offset] -= D[iy] * B[i - 1];
      }
    }

    if (D[1 + iy] != 0.0F) {
      for (i = 3; i < 5; i++) {
        D[(i + iy) - 1] -= D[1 + iy] * B[i + 3];
      }
    }

    if (D[2 + iy] != 0.0F) {
      for (i = 4; i < 5; i++) {
        D[iy + 3] -= D[2 + iy] * B[11];
      }
    }
  }

  for (j = 0; j < 4; j++) {
    iy = j << 2;
    smax = D[3 + iy];
    if (smax != 0.0F) {
      D[3 + iy] = smax / B[15];
      for (i = 0; i < 3; i++) {
        offset = i + iy;
        D[offset] -= D[3 + iy] * B[i + 12];
      }
    }

    smax = D[2 + iy];
    if (smax != 0.0F) {
      D[2 + iy] = smax / B[10];
      for (i = 0; i < 2; i++) {
        D[i + iy] -= D[2 + iy] * B[i + 8];
      }
    }

    smax = D[1 + iy];
    if (smax != 0.0F) {
      D[1 + iy] = smax / B[5];
      for (i = 0; i < 1; i++) {
        D[iy] -= D[1 + iy] * B[4];
      }
    }

    if (D[iy] != 0.0F) {
      D[iy] /= B[0];
    }
  }
}

void solve_P4Pf(const float X[12], const float u[4], const float v[4], float e,
                float *solution_num, float fs_data[], int fs_size[2], float
                Rs_data[], int Rs_size[3], float Ts_data[], int Ts_size[2])
{
  int i;
  float A[32];
  float Q[64];
  float b_A[32];
  int iy;
  float b_X[16];
  float D[16];
  int aoffset;
  int i0;
  int coffset;
  int c_data_tmp;
  float eqs[40];
  int boffset;
  float ND[48];
  float b_ND[10];
  float c_ND[10];
  float d_ND[10];
  float scale;
  float e_ND[10];
  float f_ND[10];
  float g_ND[10];
  float b_eqs[30];
  float xs_data[12];
  int xs_size[2];
  float ys_data[12];
  int ys_size[2];
  float zs_data[12];
  int zs_size[2];
  float absxk;
  float P1[3];
  float t;
  float xs_idx_2;
  float b_Q;
  float P3[3];
  float w;
  float y;
  float R2[3];
  float R[9];
  float c_A[36];
  int b_i;
  float d[9];
  float b[12];
  int j;
  float tmp_data[12];
  float a_data[36];
  float c_data[36];
  int jpvt[3];
  float b_R[12];
  float b_u[8];
  float fv0[2];
  float y_data[99];

  /* SOLVE_P4Pf Summary of this function goes here */
  /*        X = [p1, p2, p3, p4], pi = [4, 1]; X(:, i) <-> (u(i), v(i)) */
  /*        if f is a correct foal length, then [R, T] = [R, T] / sign(d)*abs(d)^(1/3); */
  /*        where d = det(R) */
  /* [p11, p12, p13, p14, p21, p22, p23, p24]' */
  for (i = 0; i < 4; i++) {
    iy = i << 2;
    b_X[iy] = X[3 * i];
    aoffset = 1 + iy;
    b_X[aoffset] = X[1 + 3 * i];
    coffset = 2 + iy;
    b_X[coffset] = X[2 + 3 * i];
    boffset = 3 + iy;
    b_X[boffset] = 1.0F;
    b_A[i] = -v[i] * b_X[iy];
    b_A[i + 16] = u[i] * b_X[iy];
    b_A[i + 4] = -v[i] * b_X[aoffset];
    b_A[i + 20] = u[i] * b_X[aoffset];
    b_A[i + 8] = -v[i] * b_X[coffset];
    b_A[i + 24] = u[i] * b_X[coffset];
    b_A[i + 12] = -v[i] * b_X[boffset];
    b_A[i + 28] = u[i] * b_X[boffset];
    for (i0 = 0; i0 < 8; i0++) {
      A[i0 + (i << 3)] = b_A[i + (i0 << 2)];
    }
  }

  qr(A, Q, b_A);

  /* nullspace */
  find_D(b_X, u, v, *(float (*)[32])&Q[32], e, D);
  for (i0 = 0; i0 < 4; i0++) {
    for (c_data_tmp = 0; c_data_tmp < 8; c_data_tmp++) {
      ND[c_data_tmp + 12 * i0] = Q[c_data_tmp + ((4 + i0) << 3)];
    }

    iy = i0 << 2;
    ND[12 * i0 + 8] = D[iy];
    ND[12 * i0 + 9] = D[1 + iy];
    ND[12 * i0 + 10] = D[2 + iy];
    ND[12 * i0 + 11] = D[3 + iy];
  }

  /* ND = [N; D] */
  memset(&eqs[0], 0, 40U * sizeof(float));

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
  for (i0 = 0; i0 < 10; i0++) {
    eqs[i0 << 2] = (b_ND[i0] + c_ND[i0]) + d_ND[i0];
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
  for (i0 = 0; i0 < 10; i0++) {
    eqs[1 + (i0 << 2)] = (b_ND[i0] + c_ND[i0]) + d_ND[i0];
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
  for (i0 = 0; i0 < 10; i0++) {
    eqs[2 + (i0 << 2)] = (b_ND[i0] + c_ND[i0]) + d_ND[i0];
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
  scale = ND[0] * ND[12];
  b_ND[3] = scale + scale;
  scale = ND[0] * ND[24];
  b_ND[4] = scale + scale;
  scale = ND[12] * ND[24];
  b_ND[5] = scale + scale;
  scale = ND[0] * ND[36];
  b_ND[6] = scale + scale;
  scale = ND[12] * ND[36];
  b_ND[7] = scale + scale;
  scale = ND[24] * ND[36];
  b_ND[8] = scale + scale;
  b_ND[9] = ND[36] * ND[36];
  c_ND[0] = ND[1] * ND[1];
  c_ND[1] = ND[13] * ND[13];
  c_ND[2] = ND[25] * ND[25];
  scale = ND[1] * ND[13];
  c_ND[3] = scale + scale;
  scale = ND[1] * ND[25];
  c_ND[4] = scale + scale;
  scale = ND[13] * ND[25];
  c_ND[5] = scale + scale;
  scale = ND[1] * ND[37];
  c_ND[6] = scale + scale;
  scale = ND[13] * ND[37];
  c_ND[7] = scale + scale;
  scale = ND[25] * ND[37];
  c_ND[8] = scale + scale;
  c_ND[9] = ND[37] * ND[37];
  d_ND[0] = ND[2] * ND[2];
  d_ND[1] = ND[14] * ND[14];
  d_ND[2] = ND[26] * ND[26];
  scale = ND[2] * ND[14];
  d_ND[3] = scale + scale;
  scale = ND[2] * ND[26];
  d_ND[4] = scale + scale;
  scale = ND[14] * ND[26];
  d_ND[5] = scale + scale;
  scale = ND[2] * ND[38];
  d_ND[6] = scale + scale;
  scale = ND[14] * ND[38];
  d_ND[7] = scale + scale;
  scale = ND[26] * ND[38];
  d_ND[8] = scale + scale;
  d_ND[9] = ND[38] * ND[38];
  e_ND[0] = ND[4] * ND[4];
  e_ND[1] = ND[16] * ND[16];
  e_ND[2] = ND[28] * ND[28];
  scale = ND[4] * ND[16];
  e_ND[3] = scale + scale;
  scale = ND[4] * ND[28];
  e_ND[4] = scale + scale;
  scale = ND[16] * ND[28];
  e_ND[5] = scale + scale;
  scale = ND[4] * ND[40];
  e_ND[6] = scale + scale;
  scale = ND[16] * ND[40];
  e_ND[7] = scale + scale;
  scale = ND[28] * ND[40];
  e_ND[8] = scale + scale;
  e_ND[9] = ND[40] * ND[40];
  f_ND[0] = ND[5] * ND[5];
  f_ND[1] = ND[17] * ND[17];
  f_ND[2] = ND[29] * ND[29];
  scale = ND[5] * ND[17];
  f_ND[3] = scale + scale;
  scale = ND[5] * ND[29];
  f_ND[4] = scale + scale;
  scale = ND[17] * ND[29];
  f_ND[5] = scale + scale;
  scale = ND[5] * ND[41];
  f_ND[6] = scale + scale;
  scale = ND[17] * ND[41];
  f_ND[7] = scale + scale;
  scale = ND[29] * ND[41];
  f_ND[8] = scale + scale;
  f_ND[9] = ND[41] * ND[41];
  g_ND[0] = ND[6] * ND[6];
  g_ND[1] = ND[18] * ND[18];
  g_ND[2] = ND[30] * ND[30];
  scale = ND[6] * ND[18];
  g_ND[3] = scale + scale;
  scale = ND[6] * ND[30];
  g_ND[4] = scale + scale;
  scale = ND[18] * ND[30];
  g_ND[5] = scale + scale;
  scale = ND[6] * ND[42];
  g_ND[6] = scale + scale;
  scale = ND[18] * ND[42];
  g_ND[7] = scale + scale;
  scale = ND[30] * ND[42];
  g_ND[8] = scale + scale;
  g_ND[9] = ND[42] * ND[42];
  for (i0 = 0; i0 < 10; i0++) {
    iy = i0 << 2;
    eqs[3 + iy] = ((((b_ND[i0] + c_ND[i0]) + d_ND[i0]) - e_ND[i0]) - f_ND[i0]) -
      g_ND[i0];
    b_eqs[3 * i0] = eqs[iy];
    b_eqs[1 + 3 * i0] = eqs[1 + iy];
    b_eqs[2 + 3 * i0] = eqs[2 + iy];
  }

  solve_3Q3(b_eqs, e, &scale, xs_data, xs_size, ys_data, ys_size, zs_data,
            zs_size);

  /* we may choes another 3 out of 4 eq-s */
  fs_size[0] = 1;
  fs_size[1] = 0;
  Rs_size[0] = 3;
  Rs_size[1] = 3;
  Rs_size[2] = 0;
  Ts_size[0] = 3;
  Ts_size[1] = 0;
  *solution_num = 0.0F;
  i0 = static_cast<int>(scale);
  for (i = 0; i < i0; i++) {
    iy = static_cast<int>((1.0F + static_cast<float>(i))) - 1;
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      P1[c_data_tmp] = ((Q[c_data_tmp + 32] * xs_data[iy] + Q[c_data_tmp + 40] *
                         ys_data[iy]) + Q[c_data_tmp + 48] * zs_data[iy]) +
        Q[c_data_tmp + 56];
    }

    absxk = xs_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    t = ys_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    xs_idx_2 = zs_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      P3[c_data_tmp] = ((D[c_data_tmp] * absxk + D[c_data_tmp + 4] * t) +
                        D[c_data_tmp + 8] * xs_idx_2) + D[c_data_tmp + 12];
    }

    b_Q = ((Q[36] * xs_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1]
            + Q[44] * ys_data[static_cast<int>((1.0F + static_cast<float>(i))) -
            1]) + Q[52] * zs_data[static_cast<int>((1.0F + static_cast<float>(i)))
           - 1]) + Q[60];
    w = std::sqrt(((P3[0] * P3[0] + P3[1] * P3[1]) + P3[2] * P3[2]) / ((P1[0] *
      P1[0] + P1[1] * P1[1]) + P1[2] * P1[2]));
    scale = 1.29246971E-26F;
    absxk = std::abs(P1[0]);
    if (absxk > 1.29246971E-26F) {
      y = 1.0F;
      scale = absxk;
    } else {
      t = absxk / 1.29246971E-26F;
      y = t * t;
    }

    absxk = std::abs(P1[1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }

    absxk = std::abs(P1[2]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }

    scale *= std::sqrt(y);
    absxk = xs_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    t = ys_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    xs_idx_2 = zs_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      P1[c_data_tmp] = (((Q[c_data_tmp + 32] * absxk + Q[c_data_tmp + 40] * t) +
                         Q[c_data_tmp + 48] * xs_idx_2) + Q[c_data_tmp + 56]) /
        scale;
    }

    y = w * scale;
    absxk = xs_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    t = ys_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    xs_idx_2 = zs_data[static_cast<int>((1.0F + static_cast<float>(i))) - 1];
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      P3[c_data_tmp] = (((D[c_data_tmp] * absxk + D[c_data_tmp + 4] * t) +
                         D[c_data_tmp + 8] * xs_idx_2) + D[c_data_tmp + 12]) / y;
    }

    R2[0] = P1[1] * P3[2] - P1[2] * P3[1];
    R2[1] = P1[2] * P3[0] - P1[0] * P3[2];
    R2[2] = P1[0] * P3[1] - P1[1] * P3[0];

    /* R2 = -R2/norm(R2); */
    scale = R2[0];
    if (R2[0] < 0.0F) {
      scale = -1.0F;
    } else if (R2[0] > 0.0F) {
      scale = 1.0F;
    } else {
      if (R2[0] == 0.0F) {
        scale = 0.0F;
      }
    }

    if (b_Q < 0.0F) {
      b_Q = -1.0F;
    } else if (b_Q > 0.0F) {
      b_Q = 1.0F;
    } else {
      if (b_Q == 0.0F) {
        b_Q = 0.0F;
      }
    }

    if (scale == b_Q) {
      R[0] = P1[0];
      R[1] = R2[0];
      R[2] = P3[0];
      R[3] = P1[1];
      R[4] = R2[1];
      R[5] = P3[1];
      R[6] = P1[2];
      R[7] = R2[2];
      R[8] = P3[2];
    } else {
      R[0] = P1[0];
      R[1] = -R2[0];
      R[2] = P3[0];
      R[3] = P1[1];
      R[4] = -R2[1];
      R[5] = P3[1];
      R[6] = P1[2];
      R[7] = -R2[2];
      R[8] = P3[2];
    }

    if (det(R) < 0.0F) {
      for (c_data_tmp = 0; c_data_tmp < 9; c_data_tmp++) {
        R[c_data_tmp] = -R[c_data_tmp];
      }
    }

    memset(&c_A[0], 0, 36U * sizeof(float));
    for (b_i = 0; b_i < 12; b_i++) {
      b[b_i] = 0.0F;
    }

    d[0] = 0.0F;
    d[3] = -1.0F;
    d[6] = w * v[0];
    d[1] = 1.0F;
    d[4] = 0.0F;
    d[7] = -w * u[0];
    d[2] = -v[0];
    d[5] = u[0];
    d[8] = 0.0F;
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        c_A[aoffset + 12 * c_data_tmp] = d[aoffset + 3 * c_data_tmp];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        a_data[aoffset + 3 * c_data_tmp] = -c_A[aoffset + 12 * c_data_tmp];
      }
    }

    for (j = 0; j < 3; j++) {
      coffset = j * 3;
      boffset = j * 3;
      memset(&c_data[coffset], 0, (unsigned int)(3 * static_cast<int>(sizeof
               (float))));
      for (iy = 0; iy < 3; iy++) {
        aoffset = iy * 3;
        scale = R[boffset + iy];
        for (b_i = 0; b_i < 3; b_i++) {
          c_data_tmp = coffset + b_i;
          c_data[c_data_tmp] += scale * a_data[aoffset + b_i];
        }
      }
    }

    memset(&tmp_data[0], 0, (unsigned int)(3 * static_cast<int>(sizeof(float))));
    for (iy = 0; iy < 3; iy++) {
      aoffset = iy * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        tmp_data[b_i] += b_X[iy] * c_data[aoffset + b_i];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      b[c_data_tmp] = tmp_data[c_data_tmp];
    }

    d[0] = 0.0F;
    d[3] = -1.0F;
    d[6] = w * v[1];
    d[1] = 1.0F;
    d[4] = 0.0F;
    d[7] = -w * u[1];
    d[2] = -v[1];
    d[5] = u[1];
    d[8] = 0.0F;
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        c_A[(aoffset + 12 * c_data_tmp) + 3] = d[aoffset + 3 * c_data_tmp];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        a_data[aoffset + 3 * c_data_tmp] = -c_A[(aoffset + 12 * c_data_tmp) + 3];
      }
    }

    for (j = 0; j < 3; j++) {
      coffset = j * 3;
      boffset = j * 3;
      memset(&c_data[coffset], 0, (unsigned int)(3 * static_cast<int>(sizeof
               (float))));
      for (iy = 0; iy < 3; iy++) {
        aoffset = iy * 3;
        scale = R[boffset + iy];
        for (b_i = 0; b_i < 3; b_i++) {
          c_data_tmp = coffset + b_i;
          c_data[c_data_tmp] += scale * a_data[aoffset + b_i];
        }
      }
    }

    memset(&tmp_data[0], 0, (unsigned int)(3 * static_cast<int>(sizeof(float))));
    for (iy = 0; iy < 3; iy++) {
      aoffset = iy * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        tmp_data[b_i] += b_X[iy + 4] * c_data[aoffset + b_i];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      b[c_data_tmp + 3] = tmp_data[c_data_tmp];
    }

    d[0] = 0.0F;
    d[3] = -1.0F;
    d[6] = w * v[2];
    d[1] = 1.0F;
    d[4] = 0.0F;
    d[7] = -w * u[2];
    d[2] = -v[2];
    d[5] = u[2];
    d[8] = 0.0F;
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        c_A[(aoffset + 12 * c_data_tmp) + 6] = d[aoffset + 3 * c_data_tmp];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        a_data[aoffset + 3 * c_data_tmp] = -c_A[(aoffset + 12 * c_data_tmp) + 6];
      }
    }

    for (j = 0; j < 3; j++) {
      coffset = j * 3;
      boffset = j * 3;
      memset(&c_data[coffset], 0, (unsigned int)(3 * static_cast<int>(sizeof
               (float))));
      for (iy = 0; iy < 3; iy++) {
        aoffset = iy * 3;
        scale = R[boffset + iy];
        for (b_i = 0; b_i < 3; b_i++) {
          c_data[coffset + b_i] += scale * a_data[aoffset + b_i];
        }
      }
    }

    memset(&tmp_data[0], 0, (unsigned int)(3 * static_cast<int>(sizeof(float))));
    for (iy = 0; iy < 3; iy++) {
      aoffset = iy * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        tmp_data[b_i] += b_X[iy + 8] * c_data[aoffset + b_i];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      b[c_data_tmp + 6] = tmp_data[c_data_tmp];
    }

    d[0] = 0.0F;
    d[3] = -1.0F;
    d[6] = w * v[3];
    d[1] = 1.0F;
    d[4] = 0.0F;
    d[7] = -w * u[3];
    d[2] = -v[3];
    d[5] = u[3];
    d[8] = 0.0F;
    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        c_A[(aoffset + 12 * c_data_tmp) + 9] = d[aoffset + 3 * c_data_tmp];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 3; aoffset++) {
        a_data[aoffset + 3 * c_data_tmp] = -c_A[(aoffset + 12 * c_data_tmp) + 9];
      }
    }

    for (j = 0; j < 3; j++) {
      coffset = j * 3;
      boffset = j * 3;
      memset(&c_data[coffset], 0, (unsigned int)(3 * static_cast<int>(sizeof
               (float))));
      for (iy = 0; iy < 3; iy++) {
        aoffset = iy * 3;
        scale = R[boffset + iy];
        for (b_i = 0; b_i < 3; b_i++) {
          c_data[coffset + b_i] += scale * a_data[aoffset + b_i];
        }
      }
    }

    memset(&tmp_data[0], 0, (unsigned int)(3 * static_cast<int>(sizeof(float))));
    for (iy = 0; iy < 3; iy++) {
      aoffset = iy * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        tmp_data[b_i] += b_X[iy + 12] * c_data[aoffset + b_i];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      b[c_data_tmp + 9] = tmp_data[c_data_tmp];
    }

    xgeqp3(c_A, P1, jpvt);
    aoffset = 0;
    scale = 1.43051147E-5F * std::abs(c_A[0]);
    while ((aoffset < 3) && (!(std::abs(c_A[aoffset + 12 * aoffset]) <= scale)))
    {
      aoffset++;
    }

    for (j = 0; j < 3; j++) {
      P3[j] = 0.0F;
      if (P1[j] != 0.0F) {
        scale = b[j];
        for (b_i = j + 2; b_i < 13; b_i++) {
          scale += c_A[(b_i + 12 * j) - 1] * b[b_i - 1];
        }

        scale *= P1[j];
        if (scale != 0.0F) {
          b[j] -= scale;
          for (b_i = j + 2; b_i < 13; b_i++) {
            b[b_i - 1] -= c_A[(b_i + 12 * j) - 1] * scale;
          }
        }
      }
    }

    for (b_i = 0; b_i < aoffset; b_i++) {
      P3[jpvt[b_i] - 1] = b[b_i];
    }

    for (j = aoffset; j >= 1; j--) {
      coffset = jpvt[j - 1] - 1;
      boffset = 12 * (j - 1);
      P3[coffset] /= c_A[(j + boffset) - 1];
      for (b_i = 0; b_i <= j - 2; b_i++) {
        P3[jpvt[b_i] - 1] -= P3[coffset] * c_A[b_i + boffset];
      }
    }

    /* 	 T = find_T((X(1:3, :)), (u),(v),(R),(w)); */
    P1[0] = 1.0F;
    P1[1] = 1.0F;
    P1[2] = w;
    for (c_data_tmp = 0; c_data_tmp < 9; c_data_tmp++) {
      d[c_data_tmp] = 0.0F;
    }

    for (j = 0; j < 3; j++) {
      d[j + 3 * j] = P1[j];
      b_R[3 * j] = R[3 * j];
      iy = 1 + 3 * j;
      b_R[iy] = R[iy];
      iy = 2 + 3 * j;
      b_R[iy] = R[iy];
      b_R[9 + j] = P3[j];
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 4; aoffset++) {
        tmp_data[c_data_tmp + 3 * aoffset] = (d[c_data_tmp] * b_R[3 * aoffset] +
          d[c_data_tmp + 3] * b_R[1 + 3 * aoffset]) + d[c_data_tmp + 6] * b_R[2
          + 3 * aoffset];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 3; c_data_tmp++) {
      for (aoffset = 0; aoffset < 4; aoffset++) {
        iy = aoffset << 2;
        b_R[c_data_tmp + 3 * aoffset] = ((tmp_data[c_data_tmp] * b_X[iy] +
          tmp_data[c_data_tmp + 3] * b_X[1 + iy]) + tmp_data[c_data_tmp + 6] *
          b_X[2 + iy]) + tmp_data[c_data_tmp + 9] * b_X[3 + iy];
      }
    }

    for (c_data_tmp = 0; c_data_tmp < 12; c_data_tmp++) {
      b[c_data_tmp] = b_R[c_data_tmp];
    }

    /* U_eval = bsxfun(@rdivide, U_eval, U_eval(3, :)); */
    /* reprojection error check */
    for (j = 0; j < 4; j++) {
      c_data_tmp = 2 + 3 * j;
      scale = b[c_data_tmp];
      iy = 1 + 3 * j;
      absxk = b[iy] / b[c_data_tmp];
      b[3 * j] /= b[c_data_tmp];
      b[iy] = absxk;
      b[c_data_tmp] = scale / scale;
      iy = j << 1;
      b_u[iy] = u[j];
      b_u[1 + iy] = v[j];
    }

    for (c_data_tmp = 0; c_data_tmp < 4; c_data_tmp++) {
      iy = c_data_tmp << 1;
      b_u[iy] = b[3 * c_data_tmp] - b_u[iy];
      iy++;
      b_u[iy] = b[1 + 3 * c_data_tmp] - b_u[iy];
    }

    y = 0.0F;
    for (j = 0; j < 4; j++) {
      iy = j << 1;
      scale = std::abs(b_u[iy]);
      if (rtIsNaNF(scale) || (scale > y)) {
        y = scale;
      }

      scale = std::abs(b_u[1 + iy]);
      if (rtIsNaNF(scale) || (scale > y)) {
        y = scale;
      }
    }

    if ((!rtIsInfF(y)) && (!rtIsNaNF(y))) {
      svd(b_u, fv0);
      y = fv0[0];
    }

    if (y < 0.01F / w) {
      (*solution_num)++;
      c_data_tmp = fs_size[1];
      fs_size[1]++;
      fs_data[c_data_tmp] = 1.0F / w;
      if (*solution_num == 1.0F) {
        Rs_size[0] = 3;
        Rs_size[1] = 3;
        Rs_size[2] = 1;
        for (c_data_tmp = 0; c_data_tmp < 9; c_data_tmp++) {
          Rs_data[c_data_tmp] = R[c_data_tmp];
        }
      } else {
        c_data_tmp = Rs_size[2];
        iy = -1;
        aoffset = 9 * Rs_size[2];
        for (j = 0; j < aoffset; j++) {
          iy++;
          y_data[iy] = Rs_data[j];
        }

        for (j = 0; j < 9; j++) {
          iy++;
          y_data[iy] = R[j];
        }

        Rs_size[0] = 3;
        Rs_size[1] = 3;
        Rs_size[2] = static_cast<signed char>((c_data_tmp + 1));
        aoffset = 9 * static_cast<signed char>((c_data_tmp + 1));
        if (0 <= aoffset - 1) {
          memcpy(&Rs_data[0], &y_data[0], (unsigned int)(aoffset * static_cast<
                  int>(sizeof(float))));
        }
      }

      aoffset = Ts_size[1];
      coffset = Ts_size[1] + 1;
      for (c_data_tmp = 0; c_data_tmp < aoffset; c_data_tmp++) {
        b_eqs[3 * c_data_tmp] = Ts_data[3 * c_data_tmp];
        boffset = 1 + 3 * c_data_tmp;
        b_eqs[boffset] = Ts_data[boffset];
        boffset = 2 + 3 * c_data_tmp;
        b_eqs[boffset] = Ts_data[boffset];
      }

      b_eqs[3 * Ts_size[1]] = P3[0];
      b_eqs[1 + 3 * Ts_size[1]] = P3[1];
      b_eqs[2 + 3 * Ts_size[1]] = P3[2];
      Ts_size[0] = 3;
      Ts_size[1] = coffset;
      aoffset = 3 * coffset;
      if (0 <= aoffset - 1) {
        memcpy(&Ts_data[0], &b_eqs[0], (unsigned int)(aoffset * static_cast<int>
                (sizeof(float))));
      }
    }
  }
}

/* End of code generation (solve_P4Pf.cpp) */
