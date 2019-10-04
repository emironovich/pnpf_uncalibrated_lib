//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: p35p_solver.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "p35p_solver.h"
#include "cat.h"
#include "eig.h"
#include "equations_for_groebner.h"
#include "find_f.h"
#include "mldivide.h"
#include "mult_for_groebner.h"
#include "p35p_solver_data.h"
#include "p35p_solver_initialize.h"
#include "p35p_solver_rtwutil.h"
#include "quadruple_constraint.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Function Definitions

//
// Arguments    : const double X[12]
//                const double x[4]
//                const double y[4]
//                double e
//                double *solution_num
//                double f_sol_data[]
//                int f_sol_size[2]
//                double R_sol_data[]
//                int R_sol_size[3]
//                double T_sol_data[]
//                int T_sol_size[2]
// Return Type  : void
//
void p35p_solver(const double X[12], const double x[4], const double y[4],
                 double e, double *solution_num, double f_sol_data[], int
                 f_sol_size[2], double R_sol_data[], int R_sol_size[3], double
                 T_sol_data[], int T_sol_size[2])
{
  double F[72];
  static const double R[54] = { 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0,
    0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, -2.0,
    0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  double dv[18];
  int i;
  double dv1[112];
  double G20[900];
  double C[200];
  double b_G20[400];
  double M[100];
  int b_i;
  double b_M[100];
  creal_T W[100];
  creal_T D[100];
  int i1;
  int ar_tmp;
  int br_tmp;
  double brm;
  double qy_re;
  double qy_im;
  double s;
  double fc_set_data[2];
  int fc_set_size[2];
  double fs_set_data[2];
  int fs_set_size[2];
  double f;
  int j;
  double R_xy[9];
  double R_xy_tmp;
  double b_x;
  double T[3];
  double fc_set[9];
  double b_fc_set[12];
  double R_curr[9];
  double p4[3];
  double b_X[4];
  double tmp_data[99];
  int tmp_size[3];
  double b_T_sol_data[30];
  int T_sol_data_tmp;
  if (isInitialized_p35p_solver == false) {
    p35p_solver_initialize();
  }

  std::memset(&F[0], 0, 72U * sizeof(double));
  quadruple_constraint(1.0, 2.0, 3.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i] = dv[3 * i];
    F[12 * i + 4] = dv[3 * i + 1];
    F[12 * i + 8] = dv[3 * i + 2];
  }

  quadruple_constraint(1.0, 3.0, 2.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 1] = dv[3 * i];
    F[12 * i + 5] = dv[3 * i + 1];
    F[12 * i + 9] = dv[3 * i + 2];
  }

  quadruple_constraint(2.0, 4.0, 3.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 2] = dv[3 * i];
    F[12 * i + 6] = dv[3 * i + 1];
    F[12 * i + 10] = dv[3 * i + 2];
  }

  quadruple_constraint(3.0, 4.0, 2.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 3] = dv[3 * i];
    F[12 * i + 7] = dv[3 * i + 1];
    F[12 * i + 11] = dv[3 * i + 2];
  }

  // 4x3x6
  // 4x28
  equations_for_groebner(F, dv1);
  mult_for_groebner(dv1, G20);

  // 20x45
  for (i = 0; i < 10; i++) {
    std::memcpy(&C[i * 20], &G20[i * 20 + 700], 20U * sizeof(double));
  }

  for (i = 0; i < 2; i++) {
    std::memcpy(&b_G20[i * 20], &G20[i * 20], 20U * sizeof(double));
  }

  for (i = 0; i < 4; i++) {
    std::memcpy(&b_G20[i * 20 + 40], &G20[i * 20 + 180], 20U * sizeof(double));
  }

  for (i = 0; i < 5; i++) {
    std::memcpy(&b_G20[i * 20 + 120], &G20[i * 20 + 340], 20U * sizeof(double));
  }

  for (i = 0; i < 4; i++) {
    std::memcpy(&b_G20[i * 20 + 220], &G20[i * 20 + 480], 20U * sizeof(double));
  }

  for (i = 0; i < 5; i++) {
    std::memcpy(&b_G20[i * 20 + 300], &G20[i * 20 + 600], 20U * sizeof(double));
  }

  mldivide(b_G20, C);

  // this function creates a matrix for multiplication by x in the monomial
  // basis B
  // monomial basis B = {x^3, ...., 1} -- monomials up to the 3d degree, #B = 10 
  // x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1
  std::memset(&M[0], 0, 100U * sizeof(double));
  for (b_i = 0; b_i < 4; b_i++) {
    for (i = 0; i < 10; i++) {
      M[i + 10 * b_i] = -C[(b_i + 20 * i) + 15];
    }
  }

  M[40] = 1.0;
  M[51] = 1.0;
  M[62] = 1.0;
  M[74] = 1.0;
  M[85] = 1.0;
  M[97] = 1.0;
  for (i = 0; i < 10; i++) {
    for (i1 = 0; i1 < 10; i1++) {
      b_M[i1 + 10 * i] = M[i + 10 * i1];
    }
  }

  eig(b_M, W, D);
  *solution_num = 0.0;
  f_sol_size[0] = 1;
  f_sol_size[1] = 0;
  R_sol_size[0] = 3;
  R_sol_size[1] = 3;
  R_sol_size[2] = 0;
  T_sol_size[0] = 3;
  T_sol_size[1] = 0;
  for (b_i = 0; b_i < 10; b_i++) {
    ar_tmp = 10 * b_i + 8;
    br_tmp = 10 * b_i + 9;
    if (W[br_tmp].im == 0.0) {
      if (W[ar_tmp].im == 0.0) {
        qy_re = W[ar_tmp].re / W[br_tmp].re;
        qy_im = 0.0;
      } else if (W[ar_tmp].re == 0.0) {
        qy_re = 0.0;
        qy_im = W[ar_tmp].im / W[br_tmp].re;
      } else {
        qy_re = W[ar_tmp].re / W[br_tmp].re;
        qy_im = W[ar_tmp].im / W[br_tmp].re;
      }
    } else if (W[br_tmp].re == 0.0) {
      if (W[ar_tmp].re == 0.0) {
        qy_re = W[ar_tmp].im / W[br_tmp].im;
        qy_im = 0.0;
      } else if (W[ar_tmp].im == 0.0) {
        qy_re = 0.0;
        qy_im = -(W[ar_tmp].re / W[br_tmp].im);
      } else {
        qy_re = W[ar_tmp].im / W[br_tmp].im;
        qy_im = -(W[ar_tmp].re / W[br_tmp].im);
      }
    } else {
      brm = std::abs(W[br_tmp].re);
      qy_im = std::abs(W[br_tmp].im);
      if (brm > qy_im) {
        s = W[br_tmp].im / W[br_tmp].re;
        qy_im = W[br_tmp].re + s * W[br_tmp].im;
        qy_re = (W[ar_tmp].re + s * W[ar_tmp].im) / qy_im;
        qy_im = (W[ar_tmp].im - s * W[ar_tmp].re) / qy_im;
      } else if (qy_im == brm) {
        if (W[br_tmp].re > 0.0) {
          qy_im = 0.5;
        } else {
          qy_im = -0.5;
        }

        if (W[br_tmp].im > 0.0) {
          f = 0.5;
        } else {
          f = -0.5;
        }

        qy_re = (W[ar_tmp].re * qy_im + W[ar_tmp].im * f) / brm;
        qy_im = (W[ar_tmp].im * qy_im - W[ar_tmp].re * f) / brm;
      } else {
        s = W[br_tmp].re / W[br_tmp].im;
        qy_im = W[br_tmp].im + s * W[br_tmp].re;
        qy_re = (s * W[ar_tmp].re + W[ar_tmp].im) / qy_im;
        qy_im = (s * W[ar_tmp].im - W[ar_tmp].re) / qy_im;
      }
    }

    i = b_i + 10 * b_i;
    if ((!(std::abs(D[i].im) > e)) && (!(std::abs(qy_im) > e))) {
      find_f(F, D[i].re, qy_re, e, fc_set_data, fc_set_size, fs_set_data,
             fs_set_size, &qy_im);
      i1 = static_cast<int>(qy_im);
      if (0 <= i1 - 1) {
        qy_im = qy_re * qy_re;
        f = D[i].re * D[i].re;
        s = 1.0 / ((f + 1.0) + qy_im);
        R_xy[0] = 1.0 - 2.0 * s * qy_im;
        brm = 2.0 * s * D[i].re;
        R_xy_tmp = brm * qy_re;
        R_xy[3] = R_xy_tmp;
        R_xy[6] = 2.0 * s * qy_re;
        R_xy[1] = R_xy_tmp;
        R_xy[4] = 1.0 - 2.0 * s * f;
        R_xy[7] = -2.0 * s * D[i].re;
        R_xy[2] = -2.0 * s * qy_re;
        R_xy[5] = brm;
        R_xy[8] = 1.0 - 2.0 * s * (f + qy_im);
        b_x = x[1];
      }

      for (j = 0; j < i1; j++) {
        // li = find_lamda(R, fc, fs, Xi, Xj, xi, xj) -- signature
        R_xy_tmp = fc_set_data[j] * R_xy[0] - fs_set_data[j] * R_xy[1];
        qy_im = fc_set_data[j] * R_xy[3] - fs_set_data[j] * R_xy[4];
        f = fc_set_data[j] * R_xy[6] - fs_set_data[j] * R_xy[7];
        brm = ((-(R_xy_tmp - b_x * R_xy[2]) * (X[3] - X[0]) + -(qy_im - b_x *
                 R_xy[5]) * (X[4] - X[1])) + -(f - b_x * R_xy[8]) * (X[5] - X[2]))
          / (x[0] - x[1]);

        // T = find_translation(R, fc, fs, li, Xi, xi, yi)
        T[0] = brm * x[0] - ((R_xy_tmp * X[0] + qy_im * X[1]) + f * X[2]);
        T[1] = brm * y[0] - (((fs_set_data[j] * R_xy[0] + fc_set_data[j] * R_xy
          [1]) * X[0] + (fs_set_data[j] * R_xy[3] + fc_set_data[j] * R_xy[4]) *
                              X[1]) + (fs_set_data[j] * R_xy[6] + fc_set_data[j]
          * R_xy[7]) * X[2]);
        T[2] = brm - ((R_xy[2] * X[0] + R_xy[5] * X[1]) + R_xy[8] * X[2]);
        f = rt_hypotd_snf(fs_set_data[j], fc_set_data[j]);
        fc_set[0] = fc_set_data[j];
        fc_set[3] = -fs_set_data[j];
        fc_set[6] = 0.0;
        fc_set[1] = fs_set_data[j];
        fc_set[4] = fc_set_data[j];
        fc_set[7] = 0.0;
        fc_set[2] = 0.0;
        fc_set[5] = 0.0;
        fc_set[8] = 1.0;
        for (i = 0; i < 3; i++) {
          qy_im = fc_set[i + 3];
          ar_tmp = static_cast<int>(fc_set[i + 6]);
          for (br_tmp = 0; br_tmp < 3; br_tmp++) {
            R_curr[i + 3 * br_tmp] = (fc_set[i] * R_xy[3 * br_tmp] + qy_im *
              R_xy[3 * br_tmp + 1]) + static_cast<double>(ar_tmp) * R_xy[3 *
              br_tmp + 2];
          }
        }

        for (i = 0; i < 3; i++) {
          b_fc_set[3 * i] = R_curr[3 * i];
          ar_tmp = 3 * i + 1;
          b_fc_set[ar_tmp] = R_curr[ar_tmp];
          ar_tmp = 3 * i + 2;
          b_fc_set[ar_tmp] = R_curr[ar_tmp];
          b_fc_set[i + 9] = T[i];
          b_X[i] = X[i + 9];
        }

        for (i = 0; i < 3; i++) {
          p4[i] = ((b_fc_set[i] * b_X[0] + b_fc_set[i + 3] * b_X[1]) +
                   b_fc_set[i + 6] * b_X[2]) + b_fc_set[i + 9];
        }

        if (std::abs(p4[1] / p4[2] - y[3]) < 0.01 * f) {
          (*solution_num)++;
          R_xy_tmp = fc_set_data[j] / f;
          fc_set[0] = R_xy_tmp;
          fc_set[3] = -fs_set_data[j] / f;
          fc_set[6] = 0.0;
          fc_set[1] = fs_set_data[j] / f;
          fc_set[4] = R_xy_tmp;
          fc_set[7] = 0.0;
          fc_set[2] = 0.0;
          fc_set[5] = 0.0;
          fc_set[8] = 1.0;
          for (i = 0; i < 3; i++) {
            qy_im = fc_set[i + 3];
            ar_tmp = static_cast<int>(fc_set[i + 6]);
            for (br_tmp = 0; br_tmp < 3; br_tmp++) {
              R_curr[i + 3 * br_tmp] = (fc_set[i] * R_xy[3 * br_tmp] + qy_im *
                R_xy[3 * br_tmp + 1]) + static_cast<double>(ar_tmp) * R_xy[3 *
                br_tmp + 2];
            }
          }

          i = f_sol_size[1];
          f_sol_size[1]++;
          f_sol_data[i] = f;
          if (*solution_num == 1.0) {
            R_sol_size[0] = 3;
            R_sol_size[1] = 3;
            R_sol_size[2] = 1;
            std::memcpy(&R_sol_data[0], &R_curr[0], 9U * sizeof(double));
          } else {
            cat(R_sol_data, R_sol_size, R_curr, tmp_data, tmp_size);
            R_sol_size[0] = 3;
            R_sol_size[1] = 3;
            R_sol_size[2] = tmp_size[2];
            ar_tmp = tmp_size[0] * tmp_size[1] * tmp_size[2];
            if (0 <= ar_tmp - 1) {
              std::memcpy(&R_sol_data[0], &tmp_data[0], ar_tmp * sizeof(double));
            }
          }

          ar_tmp = T_sol_size[1];
          br_tmp = T_sol_size[1] + 1;
          for (i = 0; i < ar_tmp; i++) {
            b_T_sol_data[3 * i] = T_sol_data[3 * i];
            T_sol_data_tmp = 3 * i + 1;
            b_T_sol_data[T_sol_data_tmp] = T_sol_data[T_sol_data_tmp];
            T_sol_data_tmp = 3 * i + 2;
            b_T_sol_data[T_sol_data_tmp] = T_sol_data[T_sol_data_tmp];
          }

          b_T_sol_data[3 * T_sol_size[1]] = T[0];
          b_T_sol_data[3 * T_sol_size[1] + 1] = T[1];
          b_T_sol_data[3 * T_sol_size[1] + 2] = T[2];
          T_sol_size[0] = 3;
          T_sol_size[1] = br_tmp;
          ar_tmp = 3 * br_tmp;
          if (0 <= ar_tmp - 1) {
            std::memcpy(&T_sol_data[0], &b_T_sol_data[0], ar_tmp * sizeof(double));
          }
        }
      }
    }
  }
}

//
// File trailer for p35p_solver.cpp
//
// [EOF]
//
