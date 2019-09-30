/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find_M.cpp
 *
 * Code generation for function 'find_M'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "find_M.h"

/* Function Definitions */
void find_M(const double P[27], double M[45])
{
  double b_P[5];
  double P_tmp;
  double b_P_tmp;
  double c_P_tmp;
  double d_P_tmp;
  double e_P_tmp;
  double f_P_tmp;
  int i6;
  double g_P_tmp;
  double h_P_tmp;
  double i_P_tmp;
  double j_P_tmp;
  double k_P_tmp;
  double l_P_tmp;
  double m_P_tmp;
  double n_P_tmp;
  double o_P_tmp;
  double p_P_tmp;
  double q_P_tmp;
  double r_P_tmp;
  double s_P_tmp;
  double t_P_tmp;
  double u_P_tmp;
  double v_P_tmp;
  double w_P_tmp;
  double x_P_tmp;
  double y_P_tmp;
  memset(&M[0], 0, 45U * sizeof(double));
  b_P[0] = 0.0;
  b_P[1] = 0.0;
  P_tmp = P[9] - P[14];
  b_P_tmp = P[12] * P[10];
  b_P[2] = ((b_P_tmp - P[8]) - P[9] * P[11]) + P[11] * P_tmp;
  c_P_tmp = P[18] - P[23];
  d_P_tmp = P[12] * P[19];
  e_P_tmp = P[21] * P[10];
  b_P[3] = (((((d_P_tmp - P[17]) + e_P_tmp) - P[9] * P[20]) - P[18] * P[11]) +
            P[20] * P_tmp) + P[11] * c_P_tmp;
  f_P_tmp = P[21] * P[19];
  b_P[4] = ((f_P_tmp - P[26]) - P[18] * P[20]) + P[20] * c_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[i6] = b_P[i6];
  }

  /*  degree = 2 */
  b_P[0] = 0.0;
  b_P[1] = 0.0;
  g_P_tmp = P[12] * P[13];
  b_P[2] = ((P[6] + g_P_tmp) - P[12] * P[11]) + P[14] * P_tmp;
  h_P_tmp = P[12] * P[22];
  i_P_tmp = P[21] * P[13];
  b_P[3] = (((((P[15] + h_P_tmp) + i_P_tmp) - P[12] * P[20]) - P[21] * P[11]) +
            P[23] * P_tmp) + P[14] * c_P_tmp;
  j_P_tmp = P[21] * P[22];
  b_P[4] = ((P[24] + j_P_tmp) - P[21] * P[20]) + P[23] * c_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[15 + i6] = b_P[i6];
  }

  /*  degree = 2 */
  b_P[0] = 0.0;
  k_P_tmp = P[12] * P[7];
  b_P[1] = (k_P_tmp - P[6] * P[11]) + P[8] * P_tmp;
  l_P_tmp = P[12] * P[16];
  m_P_tmp = P[21] * P[7];
  b_P[2] = ((((l_P_tmp + m_P_tmp) - P[6] * P[20]) - P[15] * P[11]) + P[17] *
            P_tmp) + P[8] * c_P_tmp;
  n_P_tmp = P[12] * P[25];
  o_P_tmp = P[21] * P[16];
  b_P[3] = ((((n_P_tmp + o_P_tmp) - P[15] * P[20]) - P[24] * P[11]) + P[26] *
            P_tmp) + P[17] * c_P_tmp;
  P_tmp = P[21] * P[25];
  b_P[4] = (P_tmp - P[24] * P[20]) + P[26] * c_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[30 + i6] = b_P[i6];
  }

  /*  degree = 3 */
  b_P[0] = 0.0;
  b_P[1] = 0.0;
  c_P_tmp = P[13] - P[11];
  p_P_tmp = P[9] * P[10];
  b_P[2] = ((P[10] * P[14] - p_P_tmp) - P[7]) - P[11] * c_P_tmp;
  q_P_tmp = P[22] - P[20];
  r_P_tmp = P[9] * P[19];
  s_P_tmp = P[18] * P[10];
  b_P[3] = (((((P[10] * P[23] - r_P_tmp) - s_P_tmp) - P[16]) + P[19] * P[14]) -
            P[20] * c_P_tmp) - P[11] * q_P_tmp;
  t_P_tmp = P[18] * P[19];
  b_P[4] = ((P[19] * P[23] - t_P_tmp) - P[25]) - P[20] * q_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[5 + i6] = b_P[i6];
  }

  /*  degree = 2 */
  b_P[0] = 0.0;
  b_P[1] = 0.0;
  b_P[2] = ((P[8] - b_P_tmp) + P[13] * P[14]) - P[14] * c_P_tmp;
  b_P[3] = (((((P[17] - d_P_tmp) - e_P_tmp) + P[13] * P[23]) + P[22] * P[14]) -
            P[23] * c_P_tmp) - P[14] * q_P_tmp;
  b_P[4] = ((P[26] - f_P_tmp) + P[22] * P[23]) - P[23] * q_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[20 + i6] = b_P[i6];
  }

  /*  degree = 2 */
  b_P[0] = 0.0;
  u_P_tmp = P[6] * P[10];
  b_P[1] = (P[7] * P[14] - u_P_tmp) - P[8] * c_P_tmp;
  v_P_tmp = P[6] * P[19];
  w_P_tmp = P[15] * P[10];
  b_P[2] = ((((P[7] * P[23] - w_P_tmp) - v_P_tmp) + P[16] * P[14]) - P[17] *
            c_P_tmp) - P[8] * q_P_tmp;
  x_P_tmp = P[15] * P[19];
  y_P_tmp = P[24] * P[10];
  b_P[3] = ((((P[16] * P[23] - y_P_tmp) - x_P_tmp) + P[25] * P[14]) - P[26] *
            c_P_tmp) - P[17] * q_P_tmp;
  c_P_tmp = P[24] * P[19];
  b_P[4] = (P[25] * P[23] - c_P_tmp) - P[26] * q_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[35 + i6] = b_P[i6];
  }

  /*  degree = 3 */
  b_P[0] = 0.0;
  p_P_tmp -= P[11] * P[11];
  g_P_tmp -= P[14] * P[14];
  b_P_tmp = (P[9] * P[13] + b_P_tmp) - 2.0 * P[11] * P[14];
  b_P[1] = ((((2.0 * P[11] * P[8] - u_P_tmp) - P[9] * P[7]) - P[9] * p_P_tmp) -
            P[10] * g_P_tmp) - P[11] * b_P_tmp;
  d_P_tmp = ((((P[9] * P[22] + P[18] * P[13]) + d_P_tmp) + e_P_tmp) - 2.0 * P[11]
             * P[23]) - 2.0 * P[20] * P[14];
  e_P_tmp = (r_P_tmp + s_P_tmp) - 2.0 * P[11] * P[20];
  h_P_tmp = (h_P_tmp + i_P_tmp) - 2.0 * P[14] * P[23];
  b_P[2] = ((((((((((2.0 * P[11] * P[17] - P[9] * P[16]) - P[18] * P[7]) -
                   v_P_tmp) - w_P_tmp) - P[11] * d_P_tmp) + 2.0 * P[20] * P[8])
               - P[18] * p_P_tmp) - P[19] * g_P_tmp) - P[9] * e_P_tmp) - P[10] *
            h_P_tmp) - P[20] * b_P_tmp;
  i_P_tmp = t_P_tmp - P[20] * P[20];
  j_P_tmp -= P[23] * P[23];
  f_P_tmp = (P[18] * P[22] + f_P_tmp) - 2.0 * P[20] * P[23];
  b_P[3] = ((((((((((2.0 * P[11] * P[26] - P[9] * P[25]) - P[18] * P[16]) -
                   x_P_tmp) - y_P_tmp) - P[20] * d_P_tmp) + 2.0 * P[20] * P[17])
               - P[9] * i_P_tmp) - P[10] * j_P_tmp) - P[18] * e_P_tmp) - P[19] *
            h_P_tmp) - P[11] * f_P_tmp;
  b_P[4] = ((((2.0 * P[20] * P[26] - c_P_tmp) - P[18] * P[25]) - P[18] * i_P_tmp)
            - P[19] * j_P_tmp) - P[20] * f_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[10 + i6] = b_P[i6];
  }

  /*  degree = 3 */
  b_P[0] = 0.0;
  b_P[1] = ((((2.0 * P[14] * P[8] - P[6] * P[13]) - k_P_tmp) - P[12] * p_P_tmp)
            - P[13] * g_P_tmp) - P[14] * b_P_tmp;
  b_P[2] = ((((((((((2.0 * P[14] * P[17] - l_P_tmp) - m_P_tmp) - P[6] * P[22]) -
                  P[15] * P[13]) - P[14] * d_P_tmp) + 2.0 * P[23] * P[8]) - P[21]
               * p_P_tmp) - P[22] * g_P_tmp) - P[12] * e_P_tmp) - P[13] *
            h_P_tmp) - P[23] * b_P_tmp;
  b_P[3] = ((((((((((2.0 * P[14] * P[26] - n_P_tmp) - o_P_tmp) - P[15] * P[22])
                  - P[24] * P[13]) - P[23] * d_P_tmp) + 2.0 * P[23] * P[17]) -
               P[12] * i_P_tmp) - P[13] * j_P_tmp) - P[21] * e_P_tmp) - P[22] *
            h_P_tmp) - P[14] * f_P_tmp;
  b_P[4] = ((((2.0 * P[23] * P[26] - P[24] * P[22]) - P_tmp) - P[21] * i_P_tmp)
            - P[22] * j_P_tmp) - P[23] * f_P_tmp;
  for (i6 = 0; i6 < 5; i6++) {
    M[25 + i6] = b_P[i6];
  }

  /*  degree = 3 */
  b_P[0] = (((P[8] * P[8] - P[6] * p_P_tmp) - P[7] * g_P_tmp) - P[8] * b_P_tmp)
    - P[6] * P[7];
  b_P[1] = (((((((2.0 * P[8] * P[17] - P[6] * P[16]) - P[15] * P[7]) - P[8] *
                d_P_tmp) - P[15] * p_P_tmp) - P[16] * g_P_tmp) - P[6] * e_P_tmp)
            - P[7] * h_P_tmp) - P[17] * b_P_tmp;
  b_P[2] = ((((((((((((2.0 * P[8] * P[26] - P[6] * P[25]) - P[15] * P[16]) - P
                     [24] * P[7]) - P[17] * d_P_tmp) - P[24] * p_P_tmp) - P[6] *
                  i_P_tmp) - P[25] * g_P_tmp) - P[7] * j_P_tmp) - P[15] *
               e_P_tmp) - P[16] * h_P_tmp) - P[26] * b_P_tmp) - P[8] * f_P_tmp)
    + P[17] * P[17];
  b_P[3] = (((((((2.0 * P[17] * P[26] - P[15] * P[25]) - P[24] * P[16]) - P[26] *
                d_P_tmp) - P[15] * i_P_tmp) - P[16] * j_P_tmp) - P[24] * e_P_tmp)
            - P[25] * h_P_tmp) - P[17] * f_P_tmp;
  b_P[4] = (((P[26] * P[26] - P[24] * i_P_tmp) - P[25] * j_P_tmp) - P[26] *
            f_P_tmp) - P[24] * P[25];
  for (i6 = 0; i6 < 5; i6++) {
    M[40 + i6] = b_P[i6];
  }

  /*  degree = 4 */
}

/* End of code generation (find_M.cpp) */
