/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * p35p_solver.cpp
 *
 * Code generation for function 'p35p_solver'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "cat.h"
#include "find_f.h"
#include "eig.h"
#include "make_mult_matrix.h"
#include "mldivide.h"
#include "mult_for_groebner.h"
#include "equations_for_groebner.h"
#include "init_F.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
void p35p_solver(const double X[12], const double x[4], const double y[4],
                 double e, double *solution_num, double f_sol_data[], int
                 f_sol_size[2], double R_sol_data[], int R_sol_size[3], double
                 T_sol_data[], int T_sol_size[2])
{
  static const double dv0[54] = { 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0,
    0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, -2.0,
    0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  double F[72];
  double dv1[112];
  double G20[900];
  int i0;
  double b_G20[400];
  double dv2[200];
  double dv3[100];
  double dv4[100];
  creal_T W[100];
  creal_T D[100];
  int i1;
  int i;
  double ar;
  double ai;
  double br;
  double bi;
  double qx;
  double qy_re;
  double qy_im;
  double s;
  double sgnbr;
  double fc_set_data[2];
  int fc_set_size[2];
  double fs_set_data[2];
  int fs_set_size[2];
  int j;
  double R_xy[9];
  double b_x;
  double T[3];
  double fc_set[9];
  int fc_set_tmp;
  double b_fc_set[12];
  double R_curr[9];
  double p4[3];
  double b_X[4];
  double tmp_data[99];
  int tmp_size[3];
  int T_sol_size_idx_1;
  double b_T_sol_data[30];
  int T_sol_data_tmp;
  init_F(x, y, X, dv0, F);

  /* 4x3x6 */
  /* 4x28 */
  equations_for_groebner(F, dv1);
  mult_for_groebner(dv1, G20);

  /* 20x45 */
  for (i0 = 0; i0 < 10; i0++) {
    memcpy(&dv2[i0 * 20], &G20[i0 * 20 + 700], 20U * sizeof(double));
  }

  memcpy(&b_G20[0], &G20[0], 20U * sizeof(double));
  memcpy(&b_G20[20], &G20[20], 20U * sizeof(double));
  memcpy(&b_G20[40], &G20[180], 20U * sizeof(double));
  memcpy(&b_G20[60], &G20[200], 20U * sizeof(double));
  memcpy(&b_G20[80], &G20[220], 20U * sizeof(double));
  memcpy(&b_G20[100], &G20[240], 20U * sizeof(double));
  for (i0 = 0; i0 < 5; i0++) {
    memcpy(&b_G20[i0 * 20 + 120], &G20[i0 * 20 + 340], 20U * sizeof(double));
  }

  memcpy(&b_G20[220], &G20[480], 20U * sizeof(double));
  memcpy(&b_G20[240], &G20[500], 20U * sizeof(double));
  memcpy(&b_G20[260], &G20[520], 20U * sizeof(double));
  memcpy(&b_G20[280], &G20[540], 20U * sizeof(double));
  for (i0 = 0; i0 < 5; i0++) {
    memcpy(&b_G20[i0 * 20 + 300], &G20[i0 * 20 + 600], 20U * sizeof(double));
  }

  mldivide(b_G20, dv2);
  make_mult_matrix(dv2, dv3);
  for (i0 = 0; i0 < 10; i0++) {
    for (i1 = 0; i1 < 10; i1++) {
      dv4[i1 + 10 * i0] = dv3[i0 + 10 * i1];
    }
  }

  eig(dv4, W, D);
  *solution_num = 0.0;
  f_sol_size[0] = 1;
  f_sol_size[1] = 0;
  R_sol_size[0] = 3;
  R_sol_size[1] = 3;
  R_sol_size[2] = 0;
  T_sol_size[0] = 3;
  T_sol_size[1] = 0;
  for (i = 0; i < 10; i++) {
    ar = W[8 + 10 * i].re;
    ai = W[8 + 10 * i].im;
    br = W[9 + 10 * i].re;
    bi = W[9 + 10 * i].im;
    if (bi == 0.0) {
      if (ai == 0.0) {
        qy_re = ar / br;
        qy_im = 0.0;
      } else if (ar == 0.0) {
        qy_re = 0.0;
        qy_im = ai / br;
      } else {
        qy_re = ar / br;
        qy_im = ai / br;
      }
    } else if (br == 0.0) {
      if (ar == 0.0) {
        qy_re = ai / bi;
        qy_im = 0.0;
      } else if (ai == 0.0) {
        qy_re = 0.0;
        qy_im = -(ar / bi);
      } else {
        qy_re = ai / bi;
        qy_im = -(ar / bi);
      }
    } else {
      qx = std::abs(br);
      qy_im = std::abs(bi);
      if (qx > qy_im) {
        s = bi / br;
        qy_im = br + s * bi;
        qy_re = (ar + s * ai) / qy_im;
        qy_im = (ai - s * ar) / qy_im;
      } else if (qy_im == qx) {
        if (br > 0.0) {
          sgnbr = 0.5;
        } else {
          sgnbr = -0.5;
        }

        if (bi > 0.0) {
          qy_im = 0.5;
        } else {
          qy_im = -0.5;
        }

        qy_re = (ar * sgnbr + ai * qy_im) / qx;
        qy_im = (ai * sgnbr - ar * qy_im) / qx;
      } else {
        s = br / bi;
        qy_im = bi + s * br;
        qy_re = (s * ar + ai) / qy_im;
        qy_im = (s * ai - ar) / qy_im;
      }
    }

    if ((std::abs(D[i + 10 * i].im) > e) || (std::abs(qy_im) > e)) {
    } else {
      qx = D[i + 10 * i].re;
      find_f(F, qx, qy_re, e, fc_set_data, fc_set_size, fs_set_data, fs_set_size,
             &qy_im);
      i0 = static_cast<int>(qy_im);
      if (0 <= i0 - 1) {
        qy_im = qy_re * qy_re;
        sgnbr = qx * qx;
        s = 1.0 / ((1.0 + sgnbr) + qy_im);
        R_xy[0] = 1.0 - 2.0 * s * qy_im;
        br = 2.0 * s * qx;
        bi = br * qy_re;
        R_xy[3] = bi;
        R_xy[6] = 2.0 * s * qy_re;
        R_xy[1] = bi;
        R_xy[4] = 1.0 - 2.0 * s * sgnbr;
        R_xy[7] = -2.0 * s * qx;
        R_xy[2] = -2.0 * s * qy_re;
        R_xy[5] = br;
        R_xy[8] = 1.0 - 2.0 * s * (sgnbr + qy_im);
        b_x = x[1];
      }

      for (j = 0; j < i0; j++) {
        /* li = find_lamda(R, fc, fs, Xi, Xj, xi, xj) -- signature */
        qy_im = ((-((fc_set_data[j] * R_xy[0] - fs_set_data[j] * R_xy[1]) - b_x *
                    R_xy[2]) * (X[3] - X[0]) + -((fc_set_data[j] * R_xy[3] -
                    fs_set_data[j] * R_xy[4]) - b_x * R_xy[5]) * (X[4] - X[1]))
                 + -((fc_set_data[j] * R_xy[6] - fs_set_data[j] * R_xy[7]) - b_x
                     * R_xy[8]) * (X[5] - X[2])) / (x[0] - x[1]);

        /* T = find_translation(R, fc, fs, li, Xi, xi, yi) */
        T[0] = qy_im * x[0] - (((fc_set_data[j] * R_xy[0] - fs_set_data[j] *
          R_xy[1]) * X[0] + (fc_set_data[j] * R_xy[3] - fs_set_data[j] * R_xy[4])
          * X[1]) + (fc_set_data[j] * R_xy[6] - fs_set_data[j] * R_xy[7]) * X[2]);
        T[1] = qy_im * y[0] - (((fs_set_data[j] * R_xy[0] + fc_set_data[j] *
          R_xy[1]) * X[0] + (fs_set_data[j] * R_xy[3] + fc_set_data[j] * R_xy[4])
          * X[1]) + (fs_set_data[j] * R_xy[6] + fc_set_data[j] * R_xy[7]) * X[2]);
        T[2] = qy_im - ((R_xy[2] * X[0] + R_xy[5] * X[1]) + R_xy[8] * X[2]);
        qy_im = rt_hypotd_snf(fs_set_data[j], fc_set_data[j]);
        fc_set[0] = fc_set_data[j];
        fc_set[3] = -fs_set_data[j];
        fc_set[6] = 0.0;
        fc_set[1] = fs_set_data[j];
        fc_set[4] = fc_set_data[j];
        fc_set[7] = 0.0;
        fc_set[2] = 0.0;
        fc_set[5] = 0.0;
        fc_set[8] = 1.0;
        for (i1 = 0; i1 < 3; i1++) {
          for (fc_set_tmp = 0; fc_set_tmp < 3; fc_set_tmp++) {
            R_curr[i1 + 3 * fc_set_tmp] = (fc_set[i1] * R_xy[3 * fc_set_tmp] +
              fc_set[i1 + 3] * R_xy[1 + 3 * fc_set_tmp]) + fc_set[i1 + 6] *
              R_xy[2 + 3 * fc_set_tmp];
          }
        }

        for (i1 = 0; i1 < 3; i1++) {
          b_fc_set[3 * i1] = R_curr[3 * i1];
          fc_set_tmp = 1 + 3 * i1;
          b_fc_set[fc_set_tmp] = R_curr[fc_set_tmp];
          fc_set_tmp = 2 + 3 * i1;
          b_fc_set[fc_set_tmp] = R_curr[fc_set_tmp];
          b_fc_set[9 + i1] = T[i1];
          b_X[i1] = X[9 + i1];
        }

        for (i1 = 0; i1 < 3; i1++) {
          p4[i1] = ((b_fc_set[i1] * b_X[0] + b_fc_set[i1 + 3] * b_X[1]) +
                    b_fc_set[i1 + 6] * b_X[2]) + b_fc_set[i1 + 9];
        }

        if (std::abs(p4[1] / p4[2] - y[3]) < 0.01 * qy_im) {
          (*solution_num)++;
          sgnbr = fc_set_data[j] / qy_im;
          fc_set[0] = sgnbr;
          fc_set[3] = -fs_set_data[j] / qy_im;
          fc_set[6] = 0.0;
          fc_set[1] = fs_set_data[j] / qy_im;
          fc_set[4] = sgnbr;
          fc_set[7] = 0.0;
          fc_set[2] = 0.0;
          fc_set[5] = 0.0;
          fc_set[8] = 1.0;
          for (i1 = 0; i1 < 3; i1++) {
            for (fc_set_tmp = 0; fc_set_tmp < 3; fc_set_tmp++) {
              R_curr[i1 + 3 * fc_set_tmp] = (fc_set[i1] * R_xy[3 * fc_set_tmp] +
                fc_set[i1 + 3] * R_xy[1 + 3 * fc_set_tmp]) + fc_set[i1 + 6] *
                R_xy[2 + 3 * fc_set_tmp];
            }
          }

          i1 = f_sol_size[1];
          f_sol_size[1]++;
          f_sol_data[i1] = qy_im;
          if (*solution_num == 1.0) {
            R_sol_size[0] = 3;
            R_sol_size[1] = 3;
            R_sol_size[2] = 1;
            memcpy(&R_sol_data[0], &R_curr[0], 9U * sizeof(double));
          } else {
            cat(R_sol_data, R_sol_size, R_curr, tmp_data, tmp_size);
            R_sol_size[0] = 3;
            R_sol_size[1] = 3;
            R_sol_size[2] = tmp_size[2];
            fc_set_tmp = tmp_size[0] * tmp_size[1] * tmp_size[2];
            if (0 <= fc_set_tmp - 1) {
              memcpy(&R_sol_data[0], &tmp_data[0], (unsigned int)(fc_set_tmp *
                      static_cast<int>(sizeof(double))));
            }
          }

          fc_set_tmp = T_sol_size[1];
          T_sol_size_idx_1 = T_sol_size[1] + 1;
          for (i1 = 0; i1 < fc_set_tmp; i1++) {
            b_T_sol_data[3 * i1] = T_sol_data[3 * i1];
            T_sol_data_tmp = 1 + 3 * i1;
            b_T_sol_data[T_sol_data_tmp] = T_sol_data[T_sol_data_tmp];
            T_sol_data_tmp = 2 + 3 * i1;
            b_T_sol_data[T_sol_data_tmp] = T_sol_data[T_sol_data_tmp];
          }

          b_T_sol_data[3 * T_sol_size[1]] = T[0];
          b_T_sol_data[1 + 3 * T_sol_size[1]] = T[1];
          b_T_sol_data[2 + 3 * T_sol_size[1]] = T[2];
          T_sol_size[0] = 3;
          T_sol_size[1] = T_sol_size_idx_1;
          fc_set_tmp = 3 * T_sol_size_idx_1;
          if (0 <= fc_set_tmp - 1) {
            memcpy(&T_sol_data[0], &b_T_sol_data[0], (unsigned int)(fc_set_tmp *
                    static_cast<int>(sizeof(double))));
          }
        }
      }
    }
  }
}

/* End of code generation (p35p_solver.cpp) */
