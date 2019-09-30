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
#include "find_T.h"
#include "det.h"
#include "solve_3Q3.h"
#include "find_eqs.h"
#include "qr.h"

/* Function Declarations */
static void find_A(const double X[16], const double u[4], const double v[4],
                   double A[32]);
static void find_D(const double X[16], const double u[4], const double v[4],
                   const double ns[32], double e, double D[16]);

/* Function Definitions */
static void find_A(const double X[16], const double u[4], const double v[4],
                   double A[32])
{
  int i;
  int i2;
  double d0;

  /* [p11, p12, p13, p14, p21, p22, p23, p24]' */
  for (i = 0; i < 4; i++) {
    i2 = i << 2;
    A[i] = -v[i] * X[i2];
    A[i + 16] = u[i] * X[i2];
    d0 = X[1 + i2];
    A[i + 4] = -v[i] * d0;
    A[i + 20] = u[i] * d0;
    d0 = X[2 + i2];
    A[i + 8] = -v[i] * d0;
    A[i + 24] = u[i] * d0;
    d0 = X[3 + i2];
    A[i + 12] = -v[i] * d0;
    A[i + 28] = u[i] * d0;
  }
}

static void find_D(const double X[16], const double u[4], const double v[4],
                   const double ns[32], double e, double D[16])
{
  int i;
  int j;
  double smax;
  int X_tmp;
  int offset;
  int jj;
  int iy;
  int jp1j;
  int n;
  signed char ipiv[4];
  int ix;
  int jA;
  double B[16];
  double s;
  int i3;
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

    if (B[jj + offset] != 0.0) {
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

      i3 = jj - j;
      for (i = jp1j; i <= i3 + 4; i++) {
        B[i - 1] /= B[jj];
      }
    }

    n = 2 - j;
    offset = X_tmp + 4;
    jA = jj + 5;
    for (X_tmp = 0; X_tmp <= n; X_tmp++) {
      smax = B[offset];
      if (B[offset] != 0.0) {
        ix = jj + 1;
        i3 = jA + 1;
        iy = (jA - j) + 3;
        for (i = i3; i <= iy; i++) {
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
    if (D[iy] != 0.0) {
      for (i = 2; i < 5; i++) {
        offset = (i + iy) - 1;
        D[offset] -= D[iy] * B[i - 1];
      }
    }

    if (D[1 + iy] != 0.0) {
      for (i = 3; i < 5; i++) {
        D[(i + iy) - 1] -= D[1 + iy] * B[i + 3];
      }
    }

    if (D[2 + iy] != 0.0) {
      for (i = 4; i < 5; i++) {
        D[iy + 3] -= D[2 + iy] * B[11];
      }
    }
  }

  for (j = 0; j < 4; j++) {
    iy = j << 2;
    smax = D[3 + iy];
    if (smax != 0.0) {
      D[3 + iy] = smax / B[15];
      for (i = 0; i < 3; i++) {
        offset = i + iy;
        D[offset] -= D[3 + iy] * B[i + 12];
      }
    }

    smax = D[2 + iy];
    if (smax != 0.0) {
      D[2 + iy] = smax / B[10];
      for (i = 0; i < 2; i++) {
        D[i + iy] -= D[2 + iy] * B[i + 8];
      }
    }

    smax = D[1 + iy];
    if (smax != 0.0) {
      D[1 + iy] = smax / B[5];
      for (i = 0; i < 1; i++) {
        D[iy] -= D[1 + iy] * B[4];
      }
    }

    if (D[iy] != 0.0) {
      D[iy] /= B[0];
    }
  }
}

void solve_P4Pf(const double X[12], const double u[4], const double v[4], double
                e, double *solution_num, double fs_data[], int fs_size[2],
                double Rs_data[], int Rs_size[3], double Ts_data[], int Ts_size
                [2])
{
  int i0;
  double b_X[16];
  double unusedU0[32];
  int iy;
  double b_unusedU0[32];
  double Q[64];
  int i1;
  double D[16];
  double b_Q[48];
  double dv0[40];
  double b_Ts_data[30];
  double scale;
  double xs_data[12];
  int xs_size[2];
  double ys_data[12];
  int ys_size[2];
  double zs_data[12];
  int zs_size[2];
  int i;
  double P1[3];
  double c_Q;
  double P3[3];
  double w;
  double absxk;
  double t;
  double alpha;
  double R2[3];
  double R[9];
  double U_eval[12];
  int U_eval_tmp;
  double d[9];
  int j;
  int u_tmp;
  double b_d[12];
  double b_u[8];
  double dv1[2];
  double b_Rs_data[99];

  /* SOLVE_P4Pf Summary of this function goes here */
  /*        X = [p1, p2, p3, p4], pi = [4, 1]; X(:, i) <-> (u(i), v(i)) */
  /*        if f is a correct foal length, then [R, T] = [R, T] / sign(d)*abs(d)^(1/3); */
  /*        where d = det(R) */
  for (i0 = 0; i0 < 4; i0++) {
    iy = i0 << 2;
    b_X[iy] = X[3 * i0];
    b_X[1 + iy] = X[1 + 3 * i0];
    b_X[2 + iy] = X[2 + 3 * i0];
    b_X[3 + iy] = 1.0;
  }

  find_A(b_X, u, v, unusedU0);
  for (i0 = 0; i0 < 4; i0++) {
    for (i1 = 0; i1 < 8; i1++) {
      b_unusedU0[i1 + (i0 << 3)] = unusedU0[i0 + (i1 << 2)];
    }
  }

  qr(b_unusedU0, Q, unusedU0);

  /* nullspace */
  find_D(b_X, u, v, *(double (*)[32])&Q[32], e, D);
  for (i0 = 0; i0 < 4; i0++) {
    memcpy(&b_Q[i0 * 12], &Q[(i0 << 3) + 32], sizeof(double) << 3);
    iy = i0 << 2;
    b_Q[12 * i0 + 8] = D[iy];
    b_Q[12 * i0 + 9] = D[1 + iy];
    b_Q[12 * i0 + 10] = D[2 + iy];
    b_Q[12 * i0 + 11] = D[3 + iy];
  }

  find_eqs(b_Q, dv0);
  for (i0 = 0; i0 < 10; i0++) {
    i1 = i0 << 2;
    b_Ts_data[3 * i0] = dv0[i1];
    b_Ts_data[1 + 3 * i0] = dv0[1 + i1];
    b_Ts_data[2 + 3 * i0] = dv0[2 + i1];
  }

  solve_3Q3(b_Ts_data, e, &scale, xs_data, xs_size, ys_data, ys_size, zs_data,
            zs_size);

  /* we may choes another 3 out of 4 eq-s */
  fs_size[0] = 1;
  fs_size[1] = 0;
  Rs_size[0] = 3;
  Rs_size[1] = 3;
  Rs_size[2] = 0;
  Ts_size[0] = 3;
  Ts_size[1] = 0;
  *solution_num = 0.0;
  i0 = static_cast<int>(scale);
  for (i = 0; i < i0; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      P1[i1] = ((Q[i1 + 32] * xs_data[i] + Q[i1 + 40] * ys_data[i]) + Q[i1 + 48]
                * zs_data[i]) + Q[i1 + 56];
    }

    for (i1 = 0; i1 < 3; i1++) {
      P3[i1] = ((D[i1] * xs_data[i] + D[i1 + 4] * ys_data[i]) + D[i1 + 8] *
                zs_data[i]) + D[i1 + 12];
    }

    c_Q = ((Q[36] * xs_data[i] + Q[44] * ys_data[i]) + Q[52] * zs_data[i]) + Q
      [60];
    w = std::sqrt(((P3[0] * P3[0] + P3[1] * P3[1]) + P3[2] * P3[2]) / ((P1[0] *
      P1[0] + P1[1] * P1[1]) + P1[2] * P1[2]));
    scale = 3.3121686421112381E-170;
    absxk = std::abs(P1[0]);
    if (absxk > 3.3121686421112381E-170) {
      alpha = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      alpha = t * t;
    }

    absxk = std::abs(P1[1]);
    if (absxk > scale) {
      t = scale / absxk;
      alpha = 1.0 + alpha * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      alpha += t * t;
    }

    absxk = std::abs(P1[2]);
    if (absxk > scale) {
      t = scale / absxk;
      alpha = 1.0 + alpha * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      alpha += t * t;
    }

    alpha = scale * std::sqrt(alpha);
    for (i1 = 0; i1 < 3; i1++) {
      P1[i1] = (((Q[i1 + 32] * xs_data[i] + Q[i1 + 40] * ys_data[i]) + Q[i1 + 48]
                 * zs_data[i]) + Q[i1 + 56]) / alpha;
    }

    scale = w * alpha;
    for (i1 = 0; i1 < 3; i1++) {
      P3[i1] = (((D[i1] * xs_data[i] + D[i1 + 4] * ys_data[i]) + D[i1 + 8] *
                 zs_data[i]) + D[i1 + 12]) / scale;
    }

    R2[0] = P1[1] * P3[2] - P1[2] * P3[1];
    R2[1] = P1[2] * P3[0] - P1[0] * P3[2];
    R2[2] = P1[0] * P3[1] - P1[1] * P3[0];

    /* R2 = -R2/norm(R2); */
    scale = R2[0];
    if (R2[0] < 0.0) {
      scale = -1.0;
    } else if (R2[0] > 0.0) {
      scale = 1.0;
    } else {
      if (R2[0] == 0.0) {
        scale = 0.0;
      }
    }

    if (c_Q < 0.0) {
      c_Q = -1.0;
    } else if (c_Q > 0.0) {
      c_Q = 1.0;
    } else {
      if (c_Q == 0.0) {
        c_Q = 0.0;
      }
    }

    if (scale == c_Q) {
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

    if (det(R) < 0.0) {
      for (i1 = 0; i1 < 9; i1++) {
        R[i1] = -R[i1];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      U_eval_tmp = i1 << 2;
      U_eval[3 * i1] = b_X[U_eval_tmp];
      U_eval[1 + 3 * i1] = b_X[1 + U_eval_tmp];
      U_eval[2 + 3 * i1] = b_X[2 + U_eval_tmp];
    }

    find_T(U_eval, u, v, R, w, R2);
    P1[0] = R2[0];
    P1[1] = R2[1];
    P1[2] = R2[2];
    P3[0] = 1.0;
    P3[1] = 1.0;
    P3[2] = w;
    memset(&d[0], 0, 9U * sizeof(double));
    for (j = 0; j < 3; j++) {
      d[j + 3 * j] = P3[j];
      U_eval[3 * j] = R[3 * j];
      U_eval_tmp = 1 + 3 * j;
      U_eval[U_eval_tmp] = R[U_eval_tmp];
      U_eval_tmp = 2 + 3 * j;
      U_eval[U_eval_tmp] = R[U_eval_tmp];
      U_eval[9 + j] = P1[j];
    }

    for (i1 = 0; i1 < 3; i1++) {
      for (u_tmp = 0; u_tmp < 4; u_tmp++) {
        b_d[i1 + 3 * u_tmp] = (d[i1] * U_eval[3 * u_tmp] + d[i1 + 3] * U_eval[1
          + 3 * u_tmp]) + d[i1 + 6] * U_eval[2 + 3 * u_tmp];
      }
    }

    for (i1 = 0; i1 < 3; i1++) {
      for (u_tmp = 0; u_tmp < 4; u_tmp++) {
        iy = u_tmp << 2;
        U_eval[i1 + 3 * u_tmp] = ((b_d[i1] * b_X[iy] + b_d[i1 + 3] * b_X[1 + iy])
          + b_d[i1 + 6] * b_X[2 + iy]) + b_d[i1 + 9] * b_X[3 + iy];
      }
    }

    /* U_eval = bsxfun(@rdivide, U_eval, U_eval(3, :)); */
    /* reprojection error check */
    for (j = 0; j < 4; j++) {
      i1 = 2 + 3 * j;
      scale = U_eval[i1];
      U_eval[3 * j] /= U_eval[i1];
      U_eval_tmp = 1 + 3 * j;
      U_eval[U_eval_tmp] /= scale;
      U_eval[i1] /= scale;
      iy = j << 1;
      b_u[iy] = u[j];
      u_tmp = 1 + iy;
      b_u[u_tmp] = v[j];
      b_u[iy] = U_eval[3 * j] - b_u[iy];
      b_u[u_tmp] = U_eval[U_eval_tmp] - b_u[u_tmp];
    }

    scale = 0.0;
    for (j = 0; j < 4; j++) {
      iy = j << 1;
      absxk = std::abs(b_u[iy]);
      if (rtIsNaN(absxk) || (absxk > scale)) {
        scale = absxk;
      }

      absxk = std::abs(b_u[1 + iy]);
      if (rtIsNaN(absxk) || (absxk > scale)) {
        scale = absxk;
      }
    }

    if ((!rtIsInf(scale)) && (!rtIsNaN(scale))) {
      svd(b_u, dv1);
      scale = dv1[0];
    }

    if (scale < 0.01 / w) {
      (*solution_num)++;
      i1 = fs_size[1];
      fs_size[1]++;
      fs_data[i1] = 1.0 / w;
      if (*solution_num == 1.0) {
        Rs_size[0] = 3;
        Rs_size[1] = 3;
        Rs_size[2] = 1;
        memcpy(&Rs_data[0], &R[0], 9U * sizeof(double));
      } else {
        i1 = Rs_size[2];
        iy = -1;
        u_tmp = 9 * Rs_size[2];
        for (j = 0; j < u_tmp; j++) {
          iy++;
          b_Rs_data[iy] = Rs_data[j];
        }

        for (j = 0; j < 9; j++) {
          iy++;
          b_Rs_data[iy] = R[j];
        }

        Rs_size[0] = 3;
        Rs_size[1] = 3;
        Rs_size[2] = static_cast<signed char>((i1 + 1));
        u_tmp = 9 * static_cast<signed char>((i1 + 1));
        if (0 <= u_tmp - 1) {
          memcpy(&Rs_data[0], &b_Rs_data[0], (unsigned int)(u_tmp * static_cast<
                  int>(sizeof(double))));
        }
      }

      u_tmp = Ts_size[1];
      iy = Ts_size[1] + 1;
      for (i1 = 0; i1 < u_tmp; i1++) {
        b_Ts_data[3 * i1] = Ts_data[3 * i1];
        U_eval_tmp = 1 + 3 * i1;
        b_Ts_data[U_eval_tmp] = Ts_data[U_eval_tmp];
        U_eval_tmp = 2 + 3 * i1;
        b_Ts_data[U_eval_tmp] = Ts_data[U_eval_tmp];
      }

      b_Ts_data[3 * Ts_size[1]] = P1[0];
      b_Ts_data[1 + 3 * Ts_size[1]] = P1[1];
      b_Ts_data[2 + 3 * Ts_size[1]] = P1[2];
      Ts_size[0] = 3;
      Ts_size[1] = iy;
      u_tmp = 3 * iy;
      if (0 <= u_tmp - 1) {
        memcpy(&Ts_data[0], &b_Ts_data[0], (unsigned int)(u_tmp * static_cast<
                int>(sizeof(double))));
      }
    }
  }
}

/* End of code generation (solve_P4Pf.cpp) */
