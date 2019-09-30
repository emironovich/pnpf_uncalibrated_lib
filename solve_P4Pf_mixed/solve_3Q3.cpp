/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * solve_3Q3.cpp
 *
 * Code generation for function 'solve_3Q3'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include <string.h>
#include "solve_P4Pf.h"
#include "solve_3Q3.h"
#include "qr.h"
#include "roots.h"
#include "det.h"

/* Function Declarations */
static float rt_powf_snf(float u0, float u1);

/* Function Definitions */
static float rt_powf_snf(float u0, float u1)
{
  float y;
  float f0;
  float f1;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = rtNaNF;
  } else {
    f0 = std::abs(u0);
    f1 = std::abs(u1);
    if (rtIsInfF(u1)) {
      if (f0 == 1.0F) {
        y = 1.0F;
      } else if (f0 > 1.0F) {
        if (u1 > 0.0F) {
          y = rtInfF;
        } else {
          y = 0.0F;
        }
      } else if (u1 > 0.0F) {
        y = 0.0F;
      } else {
        y = rtInfF;
      }
    } else if (f1 == 0.0F) {
      y = 1.0F;
    } else if (f1 == 1.0F) {
      if (u1 > 0.0F) {
        y = u0;
      } else {
        y = 1.0F / u0;
      }
    } else if (u1 == 2.0F) {
      y = u0 * u0;
    } else if ((u1 == 0.5F) && (u0 >= 0.0F)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0F) && (u1 > std::floor(u1))) {
      y = rtNaNF;
    } else {
      y = std::pow(u0, u1);
    }
  }

  return y;
}

void solve_3Q3(const float c[30], float e, float *n, float xs_data[], int
               xs_size[2], float ys_data[], int ys_size[2], float zs_data[], int
               zs_size[2])
{
  float A[9];
  int i;
  float a21;
  float P[27];
  float M[45];
  int i2;
  float mons[5];
  int r1;
  float b_A[9];
  float maxval;
  float P_prime[27];
  int r2;
  int r3;
  float mons_tmp;
  int rtemp;
  int k;
  float b_mons_tmp;
  float c_mons_tmp;
  float d_mons_tmp;
  float e_mons_tmp;
  float f_mons_tmp;
  float g_mons_tmp;
  float h_mons_tmp;
  float i_mons_tmp;
  float j_mons_tmp;
  float pol[13];
  float y[20];
  float C[13];
  float b_C[13];
  boolean_T b_y;
  boolean_T exitg1;
  boolean_T b[13];
  boolean_T b_b[13];
  creal32_T xs_complex_data[12];
  int xs_complex_size[1];
  float b_xs_data[12];
  float c_A[9];

  /* c -- 3x10 coefficients matrix */
  /* SOLVE_3Q3 Summary of this function goes here */
  /*    Detailed explanation goes here */
  A[0] = c[3];
  A[3] = c[6];
  A[6] = c[15];
  A[1] = c[4];
  A[4] = c[7];
  A[7] = c[16];
  A[2] = c[5];
  A[5] = c[8];
  A[8] = c[17];
  if (std::abs(det(A)) < e) {
    *n = 0.0F;
    xs_size[0] = 0;
    xs_size[1] = 0;
    ys_size[0] = 0;
    ys_size[1] = 0;
    zs_size[0] = 0;
    zs_size[1] = 0;
  } else {
    /* [x^2, x, 1] */
    for (i = 0; i < 3; i++) {
      P[i] = 0.0F;
      P[i + 9] = -c[9 + i];
      P[i + 18] = -c[21 + i];

      /* y */
      P[3 + i] = 0.0F;
      P[i + 12] = -c[12 + i];
      P[i + 21] = -c[24 + i];

      /* z */
      P[6 + i] = -c[i];
      P[i + 15] = -c[18 + i];
      P[i + 24] = -c[27 + i];

      /* 1 */
    }

    a21 = std::abs(c[4]);
    for (i = 0; i < 3; i++) {
      for (i2 = 0; i2 < 9; i2++) {
        b_A[i2] = A[i2];
      }

      r1 = 0;
      r2 = 1;
      r3 = 2;
      maxval = std::abs(A[0]);
      if (a21 > maxval) {
        maxval = a21;
        r1 = 1;
        r2 = 0;
      }

      if (std::abs(A[2]) > maxval) {
        r1 = 2;
        r2 = 1;
        r3 = 0;
      }

      b_A[r2] = A[r2] / A[r1];
      b_A[r3] /= b_A[r1];
      b_A[3 + r2] -= b_A[r2] * b_A[3 + r1];
      b_A[3 + r3] -= b_A[r3] * b_A[3 + r1];
      b_A[6 + r2] -= b_A[r2] * b_A[6 + r1];
      b_A[6 + r3] -= b_A[r3] * b_A[6 + r1];
      if (std::abs(b_A[3 + r3]) > std::abs(b_A[3 + r2])) {
        rtemp = r2;
        r2 = r3;
        r3 = rtemp;
      }

      b_A[3 + r3] /= b_A[3 + r2];
      b_A[6 + r3] -= b_A[3 + r3] * b_A[6 + r2];
      i2 = r1 + 9 * i;
      rtemp = r2 + 9 * i;
      maxval = P[rtemp] - P[i2] * b_A[r2];
      k = r3 + 9 * i;
      mons_tmp = b_A[3 + r3];
      b_mons_tmp = b_A[6 + r3];
      c_mons_tmp = ((P[k] - P[i2] * b_A[r3]) - maxval * mons_tmp) / b_mons_tmp;
      P_prime[2 + 9 * i] = c_mons_tmp;
      d_mons_tmp = b_A[6 + r1];
      e_mons_tmp = b_A[6 + r2];
      maxval -= c_mons_tmp * e_mons_tmp;
      f_mons_tmp = b_A[3 + r2];
      maxval /= f_mons_tmp;
      P_prime[1 + 9 * i] = maxval;
      g_mons_tmp = b_A[3 + r1];
      P_prime[9 * i] = ((P[i2] - c_mons_tmp * d_mons_tmp) - maxval * g_mons_tmp)
        / b_A[r1];
      h_mons_tmp = P[i2 + 3];
      maxval = P[rtemp + 3] - h_mons_tmp * b_A[r2];
      c_mons_tmp = ((P[k + 3] - h_mons_tmp * b_A[r3]) - maxval * mons_tmp) /
        b_mons_tmp;
      P_prime[9 * i + 5] = c_mons_tmp;
      h_mons_tmp -= c_mons_tmp * d_mons_tmp;
      maxval -= c_mons_tmp * e_mons_tmp;
      maxval /= f_mons_tmp;
      P_prime[9 * i + 4] = maxval;
      h_mons_tmp -= maxval * g_mons_tmp;
      h_mons_tmp /= b_A[r1];
      P_prime[3 + 9 * i] = h_mons_tmp;
      h_mons_tmp = P[i2 + 6];
      maxval = P[rtemp + 6] - h_mons_tmp * b_A[r2];
      c_mons_tmp = ((P[k + 6] - h_mons_tmp * b_A[r3]) - maxval * mons_tmp) /
        b_mons_tmp;
      P_prime[9 * i + 8] = c_mons_tmp;
      h_mons_tmp -= c_mons_tmp * d_mons_tmp;
      maxval -= c_mons_tmp * e_mons_tmp;
      maxval /= f_mons_tmp;
      P_prime[9 * i + 7] = maxval;
      h_mons_tmp -= maxval * g_mons_tmp;
      h_mons_tmp /= b_A[r1];
      P_prime[6 + 9 * i] = h_mons_tmp;
    }

    memset(&M[0], 0, 45U * sizeof(float));
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    maxval = P_prime[9] - P_prime[14];
    mons[2] = ((P_prime[12] * P_prime[10] - P_prime[8]) - P_prime[9] * P_prime
               [11]) + P_prime[11] * maxval;
    mons_tmp = P_prime[18] - P_prime[23];
    mons[3] = (((((P_prime[12] * P_prime[19] - P_prime[17]) + P_prime[21] *
                  P_prime[10]) - P_prime[9] * P_prime[20]) - P_prime[18] *
                P_prime[11]) + P_prime[20] * maxval) + P_prime[11] * mons_tmp;
    mons[4] = ((P_prime[21] * P_prime[19] - P_prime[26]) - P_prime[18] *
               P_prime[20]) + P_prime[20] * mons_tmp;
    for (i2 = 0; i2 < 5; i2++) {
      M[i2] = mons[i2];
    }

    /*  degree = 2 */
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    mons[2] = ((P_prime[6] + P_prime[12] * P_prime[13]) - P_prime[12] * P_prime
               [11]) + P_prime[14] * maxval;
    mons[3] = (((((P_prime[15] + P_prime[12] * P_prime[22]) + P_prime[21] *
                  P_prime[13]) - P_prime[12] * P_prime[20]) - P_prime[21] *
                P_prime[11]) + P_prime[23] * maxval) + P_prime[14] * mons_tmp;
    mons[4] = ((P_prime[24] + P_prime[21] * P_prime[22]) - P_prime[21] *
               P_prime[20]) + P_prime[23] * mons_tmp;
    for (i2 = 0; i2 < 5; i2++) {
      M[15 + i2] = mons[i2];
    }

    /*  degree = 2 */
    mons[0] = 0.0F;
    maxval = P_prime[9] - P_prime[14];
    mons[1] = (P_prime[12] * P_prime[7] - P_prime[6] * P_prime[11]) + P_prime[8]
      * maxval;
    mons_tmp = P_prime[18] - P_prime[23];
    mons[2] = ((((P_prime[12] * P_prime[16] + P_prime[21] * P_prime[7]) -
                 P_prime[6] * P_prime[20]) - P_prime[15] * P_prime[11]) +
               P_prime[17] * maxval) + P_prime[8] * mons_tmp;
    mons[3] = ((((P_prime[12] * P_prime[25] + P_prime[21] * P_prime[16]) -
                 P_prime[15] * P_prime[20]) - P_prime[24] * P_prime[11]) +
               P_prime[26] * maxval) + P_prime[17] * mons_tmp;
    mons[4] = (P_prime[21] * P_prime[25] - P_prime[24] * P_prime[20]) + P_prime
      [26] * mons_tmp;
    for (i2 = 0; i2 < 5; i2++) {
      M[30 + i2] = mons[i2];
    }

    /*  degree = 3 */
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    maxval = P_prime[13] - P_prime[11];
    mons[2] = ((P_prime[10] * P_prime[14] - P_prime[9] * P_prime[10]) - P_prime
               [7]) - P_prime[11] * maxval;
    mons_tmp = P_prime[22] - P_prime[20];
    mons[3] = (((((P_prime[10] * P_prime[23] - P_prime[9] * P_prime[19]) -
                  P_prime[18] * P_prime[10]) - P_prime[16]) + P_prime[19] *
                P_prime[14]) - P_prime[20] * maxval) - P_prime[11] * mons_tmp;
    mons[4] = ((P_prime[19] * P_prime[23] - P_prime[18] * P_prime[19]) -
               P_prime[25]) - P_prime[20] * mons_tmp;
    for (i2 = 0; i2 < 5; i2++) {
      M[5 + i2] = mons[i2];
    }

    /*  degree = 2 */
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    mons[2] = ((P_prime[8] - P_prime[12] * P_prime[10]) + P_prime[13] * P_prime
               [14]) - P_prime[14] * (P_prime[13] - P_prime[11]);
    maxval = P_prime[22] - P_prime[20];
    mons[3] = (((((P_prime[17] - P_prime[12] * P_prime[19]) - P_prime[21] *
                  P_prime[10]) + P_prime[13] * P_prime[23]) + P_prime[22] *
                P_prime[14]) - P_prime[23] * (P_prime[13] - P_prime[11])) -
      P_prime[14] * maxval;
    mons[4] = ((P_prime[26] - P_prime[21] * P_prime[19]) + P_prime[22] *
               P_prime[23]) - P_prime[23] * maxval;
    for (i2 = 0; i2 < 5; i2++) {
      M[20 + i2] = mons[i2];
    }

    /*  degree = 2 */
    mons[0] = 0.0F;
    maxval = P_prime[13] - P_prime[11];
    mons_tmp = P_prime[6] * P_prime[10];
    mons[1] = (P_prime[7] * P_prime[14] - mons_tmp) - P_prime[8] * maxval;
    b_mons_tmp = P_prime[22] - P_prime[20];
    c_mons_tmp = P_prime[6] * P_prime[19];
    d_mons_tmp = P_prime[15] * P_prime[10];
    mons[2] = ((((P_prime[7] * P_prime[23] - d_mons_tmp) - c_mons_tmp) +
                P_prime[16] * P_prime[14]) - P_prime[17] * maxval) - P_prime[8] *
      b_mons_tmp;
    e_mons_tmp = P_prime[15] * P_prime[19];
    f_mons_tmp = P_prime[24] * P_prime[10];
    mons[3] = ((((P_prime[16] * P_prime[23] - f_mons_tmp) - e_mons_tmp) +
                P_prime[25] * P_prime[14]) - P_prime[26] * maxval) - P_prime[17]
      * b_mons_tmp;
    maxval = P_prime[24] * P_prime[19];
    mons[4] = (P_prime[25] * P_prime[23] - maxval) - P_prime[26] * b_mons_tmp;
    for (i2 = 0; i2 < 5; i2++) {
      M[35 + i2] = mons[i2];
    }

    /*  degree = 3 */
    mons[0] = 0.0F;
    b_mons_tmp = P_prime[9] * P_prime[10] - P_prime[11] * P_prime[11];
    g_mons_tmp = P_prime[12] * P_prime[13] - P_prime[14] * P_prime[14];
    h_mons_tmp = (P_prime[9] * P_prime[13] + P_prime[12] * P_prime[10]) - 2.0F *
      P_prime[11] * P_prime[14];
    mons[1] = ((((2.0F * P_prime[11] * P_prime[8] - mons_tmp) - P_prime[9] *
                 P_prime[7]) - P_prime[9] * b_mons_tmp) - P_prime[10] *
               g_mons_tmp) - P_prime[11] * h_mons_tmp;
    mons_tmp = ((((P_prime[9] * P_prime[22] + P_prime[18] * P_prime[13]) +
                  P_prime[12] * P_prime[19]) + P_prime[21] * P_prime[10]) - 2.0F
                * P_prime[11] * P_prime[23]) - 2.0F * P_prime[20] * P_prime[14];
    a21 = (P_prime[9] * P_prime[19] + P_prime[18] * P_prime[10]) - 2.0F *
      P_prime[11] * P_prime[20];
    i_mons_tmp = (P_prime[12] * P_prime[22] + P_prime[21] * P_prime[13]) - 2.0F *
      P_prime[14] * P_prime[23];
    mons[2] = ((((((((((2.0F * P_prime[11] * P_prime[17] - P_prime[9] * P_prime
                        [16]) - P_prime[18] * P_prime[7]) - c_mons_tmp) -
                     d_mons_tmp) - P_prime[11] * mons_tmp) + 2.0F * P_prime[20] *
                   P_prime[8]) - P_prime[18] * b_mons_tmp) - P_prime[19] *
                 g_mons_tmp) - P_prime[9] * a21) - P_prime[10] * i_mons_tmp) -
      P_prime[20] * h_mons_tmp;
    c_mons_tmp = P_prime[18] * P_prime[19] - P_prime[20] * P_prime[20];
    d_mons_tmp = P_prime[21] * P_prime[22] - P_prime[23] * P_prime[23];
    j_mons_tmp = (P_prime[18] * P_prime[22] + P_prime[21] * P_prime[19]) - 2.0F *
      P_prime[20] * P_prime[23];
    mons[3] = ((((((((((2.0F * P_prime[11] * P_prime[26] - P_prime[9] * P_prime
                        [25]) - P_prime[18] * P_prime[16]) - e_mons_tmp) -
                     f_mons_tmp) - P_prime[20] * mons_tmp) + 2.0F * P_prime[20] *
                   P_prime[17]) - P_prime[9] * c_mons_tmp) - P_prime[10] *
                 d_mons_tmp) - P_prime[18] * a21) - P_prime[19] * i_mons_tmp) -
      P_prime[11] * j_mons_tmp;
    mons[4] = ((((2.0F * P_prime[20] * P_prime[26] - maxval) - P_prime[18] *
                 P_prime[25]) - P_prime[18] * c_mons_tmp) - P_prime[19] *
               d_mons_tmp) - P_prime[20] * j_mons_tmp;
    for (i2 = 0; i2 < 5; i2++) {
      M[10 + i2] = mons[i2];
    }

    /*  degree = 3 */
    mons[0] = 0.0F;
    mons[1] = ((((2.0F * P_prime[14] * P_prime[8] - P_prime[6] * P_prime[13]) -
                 P_prime[12] * P_prime[7]) - P_prime[12] * b_mons_tmp) -
               P_prime[13] * g_mons_tmp) - P_prime[14] * h_mons_tmp;
    mons[2] = ((((((((((2.0F * P_prime[14] * P_prime[17] - P_prime[12] *
                        P_prime[16]) - P_prime[21] * P_prime[7]) - P_prime[6] *
                      P_prime[22]) - P_prime[15] * P_prime[13]) - P_prime[14] *
                    mons_tmp) + 2.0F * P_prime[23] * P_prime[8]) - P_prime[21] *
                  b_mons_tmp) - P_prime[22] * g_mons_tmp) - P_prime[12] * a21) -
               P_prime[13] * i_mons_tmp) - P_prime[23] * h_mons_tmp;
    mons[3] = ((((((((((2.0F * P_prime[14] * P_prime[26] - P_prime[12] *
                        P_prime[25]) - P_prime[21] * P_prime[16]) - P_prime[15] *
                      P_prime[22]) - P_prime[24] * P_prime[13]) - P_prime[23] *
                    mons_tmp) + 2.0F * P_prime[23] * P_prime[17]) - P_prime[12] *
                  c_mons_tmp) - P_prime[13] * d_mons_tmp) - P_prime[21] * a21) -
               P_prime[22] * i_mons_tmp) - P_prime[14] * j_mons_tmp;
    mons[4] = ((((2.0F * P_prime[23] * P_prime[26] - P_prime[24] * P_prime[22])
                 - P_prime[21] * P_prime[25]) - P_prime[21] * c_mons_tmp) -
               P_prime[22] * d_mons_tmp) - P_prime[23] * j_mons_tmp;
    for (i2 = 0; i2 < 5; i2++) {
      M[25 + i2] = mons[i2];
    }

    /*  degree = 3 */
    mons[0] = (((P_prime[8] * P_prime[8] - P_prime[6] * b_mons_tmp) - P_prime[7]
                * g_mons_tmp) - P_prime[8] * h_mons_tmp) - P_prime[6] * P_prime
      [7];
    mons[1] = (((((((2.0F * P_prime[8] * P_prime[17] - P_prime[6] * P_prime[16])
                    - P_prime[15] * P_prime[7]) - P_prime[8] * mons_tmp) -
                  P_prime[15] * b_mons_tmp) - P_prime[16] * g_mons_tmp) -
                P_prime[6] * a21) - P_prime[7] * i_mons_tmp) - P_prime[17] *
      h_mons_tmp;
    mons[2] = ((((((((((((2.0F * P_prime[8] * P_prime[26] - P_prime[6] *
                          P_prime[25]) - P_prime[15] * P_prime[16]) - P_prime[24]
                        * P_prime[7]) - P_prime[17] * mons_tmp) - P_prime[24] *
                      b_mons_tmp) - P_prime[6] * c_mons_tmp) - P_prime[25] *
                    g_mons_tmp) - P_prime[7] * d_mons_tmp) - P_prime[15] * a21)
                 - P_prime[16] * i_mons_tmp) - P_prime[26] * h_mons_tmp) -
               P_prime[8] * j_mons_tmp) + P_prime[17] * P_prime[17];
    mons[3] = (((((((2.0F * P_prime[17] * P_prime[26] - P_prime[15] * P_prime[25])
                    - P_prime[24] * P_prime[16]) - P_prime[26] * mons_tmp) -
                  P_prime[15] * c_mons_tmp) - P_prime[16] * d_mons_tmp) -
                P_prime[24] * a21) - P_prime[25] * i_mons_tmp) - P_prime[17] *
      j_mons_tmp;
    mons[4] = (((P_prime[26] * P_prime[26] - P_prime[24] * c_mons_tmp) -
                P_prime[25] * d_mons_tmp) - P_prime[26] * j_mons_tmp) - P_prime
      [24] * P_prime[25];
    for (i2 = 0; i2 < 5; i2++) {
      M[40 + i2] = mons[i2];
    }

    /*  degree = 4 */
    for (i2 = 0; i2 < 9; i2++) {
      A[i2] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 5; r1++) {
        rtemp = k + r1;
        A[rtemp] += M[40 + k] * M[20 + r1];
      }
    }

    for (i2 = 0; i2 < 9; i2++) {
      b_A[i2] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 5; r1++) {
        rtemp = k + r1;
        b_A[rtemp] += M[25 + k] * M[35 + r1];
      }
    }

    for (i2 = 0; i2 < 9; i2++) {
      A[i2] -= b_A[i2];
    }

    for (i = 0; i < 13; i++) {
      pol[i] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 9; r1++) {
        rtemp = k + r1;
        pol[rtemp] += M[k] * A[r1];
      }
    }

    rtemp = -1;
    for (r1 = 0; r1 < 10; r1++) {
      rtemp++;
      y[rtemp] = M[r1 % 5 + 5 * (1 + r1 / 5)];
    }

    for (r1 = 0; r1 < 10; r1++) {
      rtemp++;
      y[rtemp] = M[30 + (r1 % 5 + 5 * (1 + r1 / 5))];
    }

    for (i2 = 0; i2 < 9; i2++) {
      A[i2] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 5; r1++) {
        rtemp = k + r1;
        A[rtemp] += y[15 + k] * y[r1];
      }
    }

    for (i2 = 0; i2 < 9; i2++) {
      b_A[i2] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 5; r1++) {
        rtemp = k + r1;
        b_A[rtemp] += y[5 + k] * y[10 + r1];
      }
    }

    for (i2 = 0; i2 < 9; i2++) {
      A[i2] -= b_A[i2];
    }

    for (i = 0; i < 13; i++) {
      C[i] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 9; r1++) {
        rtemp = k + r1;
        C[rtemp] += M[15 + k] * A[r1];
      }
    }

    for (i2 = 0; i2 < 9; i2++) {
      A[i2] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 5; r1++) {
        rtemp = k + r1;
        A[rtemp] += M[25 + k] * M[5 + r1];
      }
    }

    for (i2 = 0; i2 < 9; i2++) {
      b_A[i2] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 5; r1++) {
        rtemp = k + r1;
        b_A[rtemp] += M[10 + k] * M[20 + r1];
      }
    }

    for (i2 = 0; i2 < 9; i2++) {
      A[i2] -= b_A[i2];
    }

    for (i = 0; i < 13; i++) {
      b_C[i] = 0.0F;
    }

    for (k = 0; k < 5; k++) {
      for (r1 = 0; r1 < 9; r1++) {
        rtemp = k + r1;
        b_C[rtemp] += M[30 + k] * A[r1];
      }
    }

    for (i = 0; i < 13; i++) {
      maxval = (pol[i] - C[i]) + b_C[i];
      pol[i] = maxval;
      b[i] = rtIsInfF(maxval);
      b_b[i] = rtIsNaNF(maxval);
    }

    b_y = true;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 13)) {
      if ((!b[k]) && (!b_b[k])) {
        b_y = false;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (b_y) {
      *n = 0.0F;
      xs_size[0] = 0;
      xs_size[1] = 0;
      ys_size[0] = 0;
      ys_size[1] = 0;
      zs_size[0] = 0;
      zs_size[1] = 0;
    } else {
      roots(pol, xs_complex_data, xs_complex_size);
      if (0 <= xs_complex_size[0] - 1) {
        memset(&b_xs_data[0], 0, (unsigned int)(xs_complex_size[0] * static_cast<
                int>(sizeof(float))));
      }

      *n = 0.0F;
      i2 = xs_complex_size[0];
      for (i = 0; i < i2; i++) {
        if (std::abs(xs_complex_data[i].im) < e) {
          (*n)++;
          b_xs_data[static_cast<int>(*n) - 1] = xs_complex_data[i].re;
        }
      }

      if (1.0F > *n) {
        rtemp = 0;
      } else {
        rtemp = static_cast<int>(*n);
      }

      xs_size[0] = 1;
      xs_size[1] = rtemp;
      if (0 <= rtemp - 1) {
        memcpy(&xs_data[0], &b_xs_data[0], (unsigned int)(rtemp * static_cast<
                int>(sizeof(float))));
      }

      ys_size[0] = 1;
      rtemp = static_cast<int>(*n);
      ys_size[1] = rtemp;
      if (0 <= rtemp - 1) {
        memset(&ys_data[0], 0, (unsigned int)(rtemp * static_cast<int>(sizeof
                 (float))));
      }

      zs_size[0] = 1;
      zs_size[1] = rtemp;
      if (0 <= rtemp - 1) {
        memset(&zs_data[0], 0, (unsigned int)(rtemp * static_cast<int>(sizeof
                 (float))));
      }

      if (0 <= static_cast<int>(*n) - 1) {
        mons[4] = 1.0F;
      }

      for (i = 0; i < rtemp; i++) {
        mons[0] = rt_powf_snf(xs_data[i], 4.0F);
        mons[1] = rt_powf_snf(xs_data[i], 3.0F);
        mons[2] = xs_data[i] * xs_data[i];
        mons[3] = xs_data[i];
        for (k = 0; k < 3; k++) {
          for (r1 = 0; r1 < 3; r1++) {
            maxval = 0.0F;
            for (i2 = 0; i2 < 5; i2++) {
              maxval += mons[i2] * M[(i2 + 5 * k) + 15 * r1];
            }

            c_A[r1 + 3 * k] = maxval;
          }
        }

        b_qr(c_A, A, b_A);
        ys_data[i] = A[6] / A[8];
        zs_data[i] = A[7] / A[8];
      }
    }
  }
}

/* End of code generation (solve_3Q3.cpp) */
