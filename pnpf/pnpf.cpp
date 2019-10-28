//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: pnpf.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Oct-2019 19:14:02
//

// Include Files
#include "pnpf.h"
#include <cmath>
#include <cstring>
#include <math.h>

// Function Declarations
static void b_cat(const float varargin_1_data[], const int varargin_1_size[3],
                  const float varargin_2[9], float y_data[], int y_size[3]);
static void b_eig(const float A[100], creal32_T V[100], creal32_T D[100]);
static int b_eml_dlahqr(float h[100], float z[100]);
static int b_eml_zlahqr(creal32_T h_data[], const int h_size[2]);
static void b_equations_for_groebner(const float F[72], float G[112]);
static void b_find_f(const float F[72], float x, float y, float e, float
                     fc_data[], int fc_size[2], float fs_data[], int fs_size[2],
                     double *n);
static void b_mldivide(const double A[400], double B[200]);
static void b_mult_for_groebner(const float G4[112], float G20[900]);
static void b_mult_poly22(const float a[6], const float b[6], double c[15]);
static void b_mult_poly42(const double d[15], const float a[6], double p[28]);
static double b_norm(const double x[8]);
static void b_qr(const float A[12], float Q[9], float R[12]);
static void b_quadruple_constraint(double i, double j, double k, const float x[4],
  const float y[4], const float X[12], const float R[54], float F_row[18]);
static float b_rcond(const float A[400]);
static creal32_T b_recip(const creal32_T y);
static void b_schur(float A[100], float V[100]);
static void b_solve_3Q3(const float c[30], double *n, float xs_data[], int
  xs_size[2], float ys_data[], int ys_size[2], float zs_data[], int zs_size[2]);
static void b_sqrt(creal_T *x);
static void b_xaxpy(int n, float a, const float x[2], float y[8], int iy0);
static void b_xdlanv2(float *a, float *b, float *c, float *d, float *rt1r, float
                      *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);
static float b_xdlapy3(float x1, float x2, float x3);
static void b_xgehrd(creal32_T a_data[], const int a_size[2]);
static void b_xgerc(int m, int n, float alpha1, int ix0, const float y[4], float
                    A[12], int ia0);
static double b_xnrm2(int n, const double x[3]);
static void b_xrot(double x[100], int ix0, int iy0, double c, double s);
static void b_xrotg(float *a, float *b, float *c, float *s);
static void b_xzgeev(const creal32_T A_data[], const int A_size[2], int *info,
                     creal32_T alpha1_data[], int alpha1_size[1], creal32_T
                     beta1_data[], int beta1_size[1]);
static void b_xzgetrf(float A[400], int ipiv[20], int *info);
static void b_xzggev(creal32_T A[100], int *info, creal32_T alpha1[10],
                     creal32_T beta1[10], creal32_T V[100]);
static void b_xzhgeqz(const creal32_T A_data[], const int A_size[2], int ilo,
                      int ihi, int *info, creal32_T alpha1_data[], int
                      alpha1_size[1], creal32_T beta1_data[], int beta1_size[1]);
static void b_xzlarf(int m, int n, int iv0, float tau, float C[100], int ic0,
                     float work[10]);
static void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn);
static void b_xztgevc(const creal32_T A[100], creal32_T V[100]);
static void c_mldivide(const float A[400], float B[200]);
static float c_norm(const float x[8]);
static void c_qr(const double A[32], double Q[64], double R[32]);
static double c_rcond(const double A[9]);
static void c_sqrt(creal32_T *x);
static void c_xgerc(int m, int n, double alpha1, int ix0, const double y[8],
                    double A[64], int ia0);
static double c_xnrm2(int n, const double x[12], int ix0);
static void c_xrot(int n, float x[100], int ix0, int iy0, float c, float s);
static void c_xzgetrf(double A[9], int ipiv[3], int *info);
static void c_xzhgeqz(creal_T A[100], int ilo, int ihi, creal_T Z[100], int
                      *info, creal_T alpha1[10], creal_T beta1[10]);
static void c_xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
                      *sn, creal32_T *r);
static void cat(const double varargin_1_data[], const int varargin_1_size[3],
                const double varargin_2[9], double y_data[], int y_size[3]);
static void d_mldivide(const double A[16], double B[16]);
static void d_qr(const double A[9], double Q[9], double R[9]);
static float d_rcond(const float A[9]);
static void d_xgerc(int m, int n, double alpha1, int ix0, const double y[3],
                    double A[9], int ia0);
static float d_xnrm2(int n, const float x[100], int ix0);
static void d_xrot(float x[100], int ix0, int iy0, float c, float s);
static void d_xzgetrf(float A[9], int ipiv[3], int *info);
static void d_xzhgeqz(creal32_T A[100], int ilo, int ihi, creal32_T Z[100], int *
                      info, creal32_T alpha1[10], creal32_T beta1[10]);
static void d_xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
                      *sn);
static void e_qr(const float A[32], float Q[64], float R[32]);
static void e_xgerc(int m, int n, float alpha1, int ix0, const float y[8], float
                    A[64], int ia0);
static float e_xnrm2(int n, const float x[3]);
static void eig(const double A[100], creal_T V[100], creal_T D[100]);
static int eml_dlahqr(double h[100], double z[100]);
static int eml_zlahqr(creal_T h_data[], const int h_size[2]);
static void equations_for_groebner(const double F[72], double G[112]);
static void f_qr(const float A[9], float Q[9], float R[9]);
static void f_xgerc(int m, int n, float alpha1, int ix0, const float y[3], float
                    A[9], int ia0);
static float f_xnrm2(int n, const float x[12], int ix0);
static void find_f(const double F[72], double x, double y, double e, double
                   fc_data[], int fc_size[2], double fs_data[], int fs_size[2],
                   double *n);
static double g_xnrm2(int n, const double x[64], int ix0);
static double h_xnrm2(int n, const creal_T x_data[], int ix0);
static double i_xnrm2(int n, const double x[9], int ix0);
static double j_xnrm2(int n, const double x[36], int ix0);
static double k_xnrm2(int n, const double x[8], int ix0);
static double l_xnrm2(int n, const double x[4], int ix0);
static float m_xnrm2(int n, const float x[64], int ix0);
static void mldivide(const double A[36], const double B[12], double Y[3]);
static void mult_for_groebner(const double G4[112], double G20[900]);
static void mult_poly22(const double a[6], const double b[6], double c[15]);
static void mult_poly42(const double d[15], const double a[6], double p[28]);
static float n_xnrm2(int n, const creal32_T x_data[], int ix0);
static float o_xnrm2(int n, const float x[9], int ix0);
static float p_xnrm2(int n, const float x[8], int ix0);
static float q_xnrm2(int n, const float x[4], int ix0);
static void qr(const double A[12], double Q[9], double R[12]);
static void qrpf(double A[36], double tau[3], int jpvt[3]);
static void quadruple_constraint(double i, double j, double k, const double x[4],
  const double y[4], const double X[12], const double R[54], double F_row[18]);
static double rcond(const double A[400]);
static creal_T recip(const creal_T y);
static void roots(const double c[9], creal_T r_data[], int r_size[1]);
static double rt_hypotd(double u0, double u1);
static float rt_hypotf(float u0, float u1);
static void schur(const double A[100], double V[100], double T[100]);
static void solve_3Q3(const double c[30], double *n, double xs_data[], int
                      xs_size[2], double ys_data[], int ys_size[2], double
                      zs_data[], int zs_size[2]);
static void xaxpy(int n, double a, const double x[2], double y[8], int iy0);
static void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn);
static double xdlapy3(double x1, double x2, double x3);
static void xgehrd(creal_T a_data[], const int a_size[2]);
static void xgerc(int m, int n, double alpha1, int ix0, const double y[4],
                  double A[12], int ia0);
static double xnrm2(int n, const double x[100], int ix0);
static void xrot(int n, double x[100], int ix0, int iy0, double c, double s);
static void xrotg(double *a, double *b, double *c, double *s);
static void xzgeev(const creal_T A_data[], const int A_size[2], int *info,
                   creal_T alpha1_data[], int alpha1_size[1], creal_T
                   beta1_data[], int beta1_size[1]);
static void xzgetrf(double A[400], int ipiv[20], int *info);
static void xzggbak(creal_T V[100], int ilo, int ihi, const int rscale[10]);
static void xzggbal(creal_T A[100], int *ilo, int *ihi, int rscale[10]);
static void xzggev(creal_T A[100], int *info, creal_T alpha1[10], creal_T beta1
                   [10], creal_T V[100]);
static void xzgghrd(int ilo, int ihi, creal_T A[100], creal_T Z[100]);
static void xzhgeqz(const creal_T A_data[], const int A_size[2], int ilo, int
                    ihi, int *info, creal_T alpha1_data[], int alpha1_size[1],
                    creal_T beta1_data[], int beta1_size[1]);
static void xzlarf(int m, int n, int iv0, double tau, double C[100], int ic0,
                   double work[10]);
static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r);
static void xztgevc(const creal_T A[100], creal_T V[100]);

// Function Definitions

//
// Arguments    : const float varargin_1_data[]
//                const int varargin_1_size[3]
//                const float varargin_2[9]
//                float y_data[]
//                int y_size[3]
// Return Type  : void
//
static void b_cat(const float varargin_1_data[], const int varargin_1_size[3],
                  const float varargin_2[9], float y_data[], int y_size[3])
{
  int iy;
  int i;
  int j;
  y_size[0] = 3;
  y_size[1] = 3;
  y_size[2] = static_cast<signed char>((varargin_1_size[2] + 1));
  iy = -1;
  i = 9 * varargin_1_size[2];
  for (j = 0; j < i; j++) {
    iy++;
    y_data[iy] = varargin_1_data[j];
  }

  for (j = 0; j < 9; j++) {
    iy++;
    y_data[iy] = varargin_2[j];
  }
}

//
// Arguments    : const float A[100]
//                creal32_T V[100]
//                creal32_T D[100]
// Return Type  : void
//
static void b_eig(const float A[100], creal32_T V[100], creal32_T D[100])
{
  boolean_T p;
  int j;
  boolean_T exitg2;
  int info;
  float b_D[100];
  float b_V[100];
  int exitg1;
  creal32_T At[100];
  creal32_T alpha1[10];
  creal32_T beta1[10];
  int coltop;
  float colnorm;
  float scale;
  float t;
  float absxk;
  p = true;
  j = 0;
  exitg2 = false;
  while ((!exitg2) && (j < 10)) {
    info = 0;
    do {
      exitg1 = 0;
      if (info <= j) {
        if (A[info + 10 * j] != A[j + 10 * info]) {
          p = false;
          exitg1 = 1;
        } else {
          info++;
        }
      } else {
        j++;
        exitg1 = 2;
      }
    } while (exitg1 == 0);

    if (exitg1 == 1) {
      exitg2 = true;
    }
  }

  if (p) {
    std::memcpy(&b_D[0], &A[0], 100U * sizeof(float));
    b_schur(b_D, b_V);
    for (info = 0; info < 100; info++) {
      V[info].re = b_V[info];
      V[info].im = 0.0F;
    }

    for (j = 0; j < 9; j++) {
      b_D[(j + 10 * j) + 1] = 0.0F;
      for (info = 0; info <= j; info++) {
        b_D[info + 10 * (j + 1)] = 0.0F;
      }
    }

    for (info = 0; info < 100; info++) {
      D[info].re = b_D[info];
      D[info].im = 0.0F;
    }
  } else {
    for (info = 0; info < 100; info++) {
      At[info].re = A[info];
      At[info].im = 0.0F;
    }

    b_xzggev(At, &info, alpha1, beta1, V);
    for (coltop = 0; coltop <= 90; coltop += 10) {
      colnorm = 0.0F;
      scale = 1.29246971E-26F;
      info = coltop + 10;
      for (j = coltop + 1; j <= info; j++) {
        absxk = std::abs(V[j - 1].re);
        if (absxk > scale) {
          t = scale / absxk;
          colnorm = colnorm * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          colnorm += t * t;
        }

        absxk = std::abs(V[j - 1].im);
        if (absxk > scale) {
          t = scale / absxk;
          colnorm = colnorm * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          colnorm += t * t;
        }
      }

      colnorm = scale * std::sqrt(colnorm);
      info = coltop + 10;
      for (j = coltop + 1; j <= info; j++) {
        absxk = V[j - 1].re;
        scale = V[j - 1].im;
        if (scale == 0.0F) {
          absxk /= colnorm;
          scale = 0.0F;
        } else if (absxk == 0.0F) {
          absxk = 0.0F;
          scale /= colnorm;
        } else {
          absxk /= colnorm;
          scale /= colnorm;
        }

        V[j - 1].re = absxk;
        V[j - 1].im = scale;
      }
    }

    std::memset(&D[0], 0, 100U * sizeof(creal32_T));
    for (j = 0; j < 10; j++) {
      if (beta1[j].im == 0.0F) {
        if (alpha1[j].im == 0.0F) {
          info = j + 10 * j;
          D[info].re = alpha1[j].re / beta1[j].re;
          D[info].im = 0.0F;
        } else if (alpha1[j].re == 0.0F) {
          info = j + 10 * j;
          D[info].re = 0.0F;
          D[info].im = alpha1[j].im / beta1[j].re;
        } else {
          info = j + 10 * j;
          D[info].re = alpha1[j].re / beta1[j].re;
          D[info].im = alpha1[j].im / beta1[j].re;
        }
      } else if (beta1[j].re == 0.0F) {
        if (alpha1[j].re == 0.0F) {
          info = j + 10 * j;
          D[info].re = alpha1[j].im / beta1[j].im;
          D[info].im = 0.0F;
        } else if (alpha1[j].im == 0.0F) {
          info = j + 10 * j;
          D[info].re = 0.0F;
          D[info].im = -(alpha1[j].re / beta1[j].im);
        } else {
          info = j + 10 * j;
          D[info].re = alpha1[j].im / beta1[j].im;
          D[info].im = -(alpha1[j].re / beta1[j].im);
        }
      } else {
        t = std::abs(beta1[j].re);
        scale = std::abs(beta1[j].im);
        if (t > scale) {
          scale = beta1[j].im / beta1[j].re;
          absxk = beta1[j].re + scale * beta1[j].im;
          info = j + 10 * j;
          D[info].re = (alpha1[j].re + scale * alpha1[j].im) / absxk;
          D[info].im = (alpha1[j].im - scale * alpha1[j].re) / absxk;
        } else if (scale == t) {
          if (beta1[j].re > 0.0F) {
            scale = 0.5F;
          } else {
            scale = -0.5F;
          }

          if (beta1[j].im > 0.0F) {
            absxk = 0.5F;
          } else {
            absxk = -0.5F;
          }

          info = j + 10 * j;
          D[info].re = (alpha1[j].re * scale + alpha1[j].im * absxk) / t;
          D[info].im = (alpha1[j].im * scale - alpha1[j].re * absxk) / t;
        } else {
          scale = beta1[j].re / beta1[j].im;
          absxk = beta1[j].im + scale * beta1[j].re;
          info = j + 10 * j;
          D[info].re = (scale * alpha1[j].re + alpha1[j].im) / absxk;
          D[info].im = (scale * alpha1[j].im - alpha1[j].re) / absxk;
        }
      }
    }
  }
}

//
// Arguments    : float h[100]
//                float z[100]
// Return Type  : int
//
static int b_eml_dlahqr(float h[100], float z[100])
{
  int info;
  float v[3];
  int j;
  int i;
  int b_i;
  boolean_T exitg1;
  int L;
  boolean_T goto150;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  int ix;
  float s;
  int nr;
  float ba;
  int hoffset;
  float f;
  float tst;
  float bb;
  float aa;
  float ab;
  float h22;
  float rt1r;
  int knt;
  int m;
  int b_k;
  info = 0;
  v[0] = 0.0F;
  v[1] = 0.0F;
  v[2] = 0.0F;
  for (j = 0; j < 7; j++) {
    i = j + 10 * j;
    h[i + 2] = 0.0F;
    h[i + 3] = 0.0F;
  }

  h[79] = 0.0F;
  b_i = 9;
  exitg1 = false;
  while ((!exitg1) && (b_i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      k = b_i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        i = k + 10 * (k - 1);
        ba = std::abs(h[i]);
        if (ba <= 9.86076132E-31F) {
          exitg3 = true;
        } else {
          ix = k + 10 * k;
          bb = std::abs(h[ix]);
          hoffset = i - 1;
          tst = std::abs(h[hoffset]) + bb;
          if (tst == 0.0F) {
            if (k - 1 >= 1) {
              tst = std::abs(h[(k + 10 * (k - 2)) - 1]);
            }

            if (k + 2 <= 10) {
              tst += std::abs(h[ix + 1]);
            }
          }

          if (ba <= 1.1920929E-7F * tst) {
            tst = std::abs(h[ix - 1]);
            if (ba > tst) {
              ab = ba;
              ba = tst;
            } else {
              ab = tst;
            }

            tst = std::abs(h[hoffset] - h[ix]);
            if (bb > tst) {
              aa = bb;
              bb = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            tst = 1.1920929E-7F * (bb * (aa / s));
            if (9.86076132E-31F > tst) {
              tst = 9.86076132E-31F;
            }

            if (ba * (ab / s) <= tst) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + 10 * (k - 1)] = 0.0F;
      }

      if (k + 1 >= b_i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          hoffset = k + 10 * k;
          s = std::abs(h[hoffset + 1]) + std::abs(h[(k + 10 * (k + 1)) + 2]);
          tst = 0.75F * s + h[hoffset];
          aa = -0.4375F * s;
          ab = s;
          h22 = tst;
        } else if (its == 20) {
          s = std::abs(h[b_i + 10 * (b_i - 1)]) + std::abs(h[(b_i + 10 * (b_i -
            2)) - 1]);
          tst = 0.75F * s + h[b_i + 10 * b_i];
          aa = -0.4375F * s;
          ab = s;
          h22 = tst;
        } else {
          ix = b_i + 10 * (b_i - 1);
          tst = h[ix - 1];
          ab = h[ix];
          aa = h[(b_i + 10 * b_i) - 1];
          h22 = h[b_i + 10 * b_i];
        }

        s = ((std::abs(tst) + std::abs(aa)) + std::abs(ab)) + std::abs(h22);
        if (s == 0.0F) {
          rt1r = 0.0F;
          tst = 0.0F;
          ba = 0.0F;
          aa = 0.0F;
        } else {
          tst /= s;
          ab /= s;
          aa /= s;
          h22 /= s;
          bb = (tst + h22) / 2.0F;
          tst = (tst - bb) * (h22 - bb) - aa * ab;
          aa = std::sqrt(std::abs(tst));
          if (tst >= 0.0F) {
            rt1r = bb * s;
            ba = rt1r;
            tst = aa * s;
            aa = -tst;
          } else {
            rt1r = bb + aa;
            ba = bb - aa;
            if (std::abs(rt1r - h22) <= std::abs(ba - h22)) {
              rt1r *= s;
              ba = rt1r;
            } else {
              ba *= s;
              rt1r = ba;
            }

            tst = 0.0F;
            aa = 0.0F;
          }
        }

        m = b_i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          hoffset = m + 10 * (m - 1);
          ix = hoffset - 1;
          ab = h[ix] - ba;
          s = (std::abs(ab) + std::abs(aa)) + std::abs(h[hoffset]);
          bb = h[hoffset] / s;
          hoffset = m + 10 * m;
          v[0] = (bb * h[hoffset - 1] + (h[ix] - rt1r) * (ab / s)) - tst * (aa /
            s);
          v[1] = bb * (((h[ix] + h[hoffset]) - rt1r) - ba);
          v[2] = bb * h[hoffset + 1];
          s = (std::abs(v[0]) + std::abs(v[1])) + std::abs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if (m == k + 1) {
            exitg3 = true;
          } else {
            i = m + 10 * (m - 2);
            if (std::abs(h[i - 1]) * (std::abs(v[1]) + std::abs(v[2])) <=
                1.1920929E-7F * std::abs(v[0]) * ((std::abs(h[i - 2]) + std::abs
                  (h[ix])) + std::abs(h[hoffset]))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (b_k = m; b_k <= b_i; b_k++) {
          nr = (b_i - b_k) + 2;
          if (3 < nr) {
            nr = 3;
          }

          if (b_k > m) {
            hoffset = (b_k + 10 * (b_k - 2)) - 1;
            for (j = 0; j < nr; j++) {
              v[j] = h[j + hoffset];
            }
          }

          aa = v[0];
          bb = 0.0F;
          if (nr > 0) {
            tst = e_xnrm2(nr - 1, v);
            if (tst != 0.0F) {
              ab = rt_hypotf(v[0], tst);
              if (v[0] >= 0.0F) {
                ab = -ab;
              }

              if (std::abs(ab) < 9.86076132E-32F) {
                knt = -1;
                do {
                  knt++;
                  for (ix = 2; ix <= nr; ix++) {
                    v[ix - 1] *= 1.01412048E+31F;
                  }

                  ab *= 1.01412048E+31F;
                  aa *= 1.01412048E+31F;
                } while (!(std::abs(ab) >= 9.86076132E-32F));

                ab = rt_hypotf(aa, e_xnrm2(nr - 1, v));
                if (aa >= 0.0F) {
                  ab = -ab;
                }

                bb = (ab - aa) / ab;
                tst = 1.0F / (aa - ab);
                for (ix = 2; ix <= nr; ix++) {
                  v[ix - 1] *= tst;
                }

                for (ix = 0; ix <= knt; ix++) {
                  ab *= 9.86076132E-32F;
                }

                aa = ab;
              } else {
                bb = (ab - v[0]) / ab;
                tst = 1.0F / (v[0] - ab);
                for (ix = 2; ix <= nr; ix++) {
                  v[ix - 1] *= tst;
                }

                aa = ab;
              }
            }
          }

          v[0] = aa;
          if (b_k > m) {
            h[(b_k + 10 * (b_k - 2)) - 1] = aa;
            i = b_k + 10 * (b_k - 2);
            h[i] = 0.0F;
            if (b_k < b_i) {
              h[i + 1] = 0.0F;
            }
          } else {
            if (m > k + 1) {
              h[(b_k + 10 * (b_k - 2)) - 1] *= 1.0F - bb;
            }
          }

          s = v[1];
          tst = bb * v[1];
          if (nr == 3) {
            f = v[2];
            ab = bb * v[2];
            for (j = b_k; j < 11; j++) {
              ix = b_k + 10 * (j - 1);
              hoffset = ix - 1;
              knt = ix + 1;
              aa = (h[hoffset] + s * h[ix]) + f * h[knt];
              h[hoffset] -= aa * bb;
              h[ix] -= aa * tst;
              h[knt] -= aa * ab;
            }

            if (b_k + 3 < b_i + 1) {
              i = b_k + 2;
            } else {
              i = b_i;
            }

            for (j = 0; j <= i; j++) {
              ix = j + 10 * (b_k - 1);
              hoffset = j + 10 * b_k;
              knt = j + 10 * (b_k + 1);
              aa = (h[ix] + s * h[hoffset]) + f * h[knt];
              h[ix] -= aa * bb;
              h[hoffset] -= aa * tst;
              h[knt] -= aa * ab;
            }

            for (j = 0; j < 10; j++) {
              ix = j + 10 * (b_k - 1);
              hoffset = j + 10 * b_k;
              knt = j + 10 * (b_k + 1);
              aa = (z[ix] + s * z[hoffset]) + f * z[knt];
              z[ix] -= aa * bb;
              z[hoffset] -= aa * tst;
              z[knt] -= aa * ab;
            }
          } else {
            if (nr == 2) {
              for (j = b_k; j < 11; j++) {
                ix = b_k + 10 * (j - 1);
                hoffset = ix - 1;
                aa = h[hoffset] + s * h[ix];
                h[hoffset] -= aa * bb;
                h[ix] -= aa * tst;
              }

              for (j = 0; j <= b_i; j++) {
                ix = j + 10 * (b_k - 1);
                hoffset = j + 10 * b_k;
                aa = h[ix] + s * h[hoffset];
                h[ix] -= aa * bb;
                h[hoffset] -= aa * tst;
              }

              for (j = 0; j < 10; j++) {
                ix = j + 10 * (b_k - 1);
                hoffset = j + 10 * b_k;
                aa = z[ix] + s * z[hoffset];
                z[ix] -= aa * bb;
                z[hoffset] -= aa * tst;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = b_i + 1;
      exitg1 = true;
    } else {
      if ((L != b_i + 1) && (L == b_i)) {
        i = b_i + 10 * b_i;
        ix = i - 1;
        s = h[ix];
        nr = 10 * (b_i - 1);
        hoffset = b_i + nr;
        f = h[hoffset];
        tst = h[i];
        b_xdlanv2(&h[(b_i + 10 * (b_i - 1)) - 1], &s, &f, &tst, &aa, &ab, &bb,
                  &ba, &h22, &rt1r);
        h[ix] = s;
        h[hoffset] = f;
        h[i] = tst;
        if (10 > b_i + 1) {
          hoffset = 8 - b_i;
          knt = b_i + (b_i + 1) * 10;
          ix = knt - 1;
          for (k = 0; k <= hoffset; k++) {
            tst = h22 * h[ix] + rt1r * h[knt];
            h[knt] = h22 * h[knt] - rt1r * h[ix];
            h[ix] = tst;
            knt += 10;
            ix += 10;
          }
        }

        ix = nr + 1;
        hoffset = b_i * 10 + 1;
        c_xrot(b_i - 1, h, ix, hoffset, h22, rt1r);
        d_xrot(z, ix, hoffset, h22, rt1r);
      }

      b_i = L - 2;
    }
  }

  return info;
}

//
// Arguments    : creal32_T h_data[]
//                const int h_size[2]
// Return Type  : int
//
static int b_eml_zlahqr(creal32_T h_data[], const int h_size[2])
{
  int info;
  int n;
  int ldh;
  int i;
  int j;
  int knt;
  int h_tmp;
  int b_i;
  float SMLNUM;
  boolean_T exitg1;
  float tst;
  float aa;
  float br;
  int L;
  boolean_T goto140;
  creal32_T sc;
  int its;
  boolean_T exitg2;
  int k;
  float re;
  boolean_T exitg3;
  float im;
  float ba;
  int ix0;
  float bb;
  float t_re;
  float ab;
  creal32_T v;
  boolean_T goto70;
  int m;
  float u_re;
  float u_im;
  float s;
  int b_k;
  creal32_T b_v[2];
  n = h_size[0];
  ldh = h_size[0];
  info = 0;
  if (1 != h_size[0]) {
    i = h_size[0];
    for (j = 0; j <= i - 4; j++) {
      knt = j + h_size[0] * j;
      h_tmp = knt + 2;
      h_data[h_tmp].re = 0.0F;
      h_data[h_tmp].im = 0.0F;
      knt += 3;
      h_data[knt].re = 0.0F;
      h_data[knt].im = 0.0F;
    }

    if (1 <= n - 2) {
      i = (n + h_size[0] * (n - 3)) - 1;
      h_data[i].re = 0.0F;
      h_data[i].im = 0.0F;
    }

    for (b_i = 2; b_i <= n; b_i++) {
      i = (b_i + h_size[0] * (b_i - 2)) - 1;
      if (h_data[i].im != 0.0F) {
        tst = h_data[i].re;
        aa = h_data[i].im;
        br = std::abs(h_data[i].re) + std::abs(h_data[i].im);
        if (aa == 0.0F) {
          sc.re = tst / br;
          sc.im = 0.0F;
        } else if (tst == 0.0F) {
          sc.re = 0.0F;
          sc.im = aa / br;
        } else {
          sc.re = tst / br;
          sc.im = aa / br;
        }

        br = rt_hypotf(sc.re, sc.im);
        if (-sc.im == 0.0F) {
          re = sc.re / br;
          im = 0.0F;
        } else if (sc.re == 0.0F) {
          re = 0.0F;
          im = -sc.im / br;
        } else {
          re = sc.re / br;
          im = -sc.im / br;
        }

        h_data[i].re = rt_hypotf(h_data[i].re, h_data[i].im);
        h_data[i].im = 0.0F;
        h_tmp = (b_i - 1) * ldh;
        ix0 = b_i + h_tmp;
        i = ix0 + ldh * (n - b_i);
        for (k = ix0; ldh < 0 ? k >= i : k <= i; k += ldh) {
          aa = re * h_data[k - 1].im + im * h_data[k - 1].re;
          h_data[k - 1].re = re * h_data[k - 1].re - im * h_data[k - 1].im;
          h_data[k - 1].im = aa;
        }

        ix0 = h_tmp + 1;
        knt = b_i + 1;
        if (n < knt) {
          knt = n;
        }

        i = h_tmp + knt;
        for (k = ix0; k <= i; k++) {
          aa = re * h_data[k - 1].im + -im * h_data[k - 1].re;
          h_data[k - 1].re = re * h_data[k - 1].re - -im * h_data[k - 1].im;
          h_data[k - 1].im = aa;
        }
      }
    }

    SMLNUM = 1.17549435E-38F * (static_cast<float>(n) / 1.1920929E-7F);
    b_i = n - 1;
    exitg1 = false;
    while ((!exitg1) && (b_i + 1 >= 1)) {
      L = -1;
      goto140 = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its < 301)) {
        k = b_i;
        exitg3 = false;
        while ((!exitg3) && (k + 1 > L + 2)) {
          i = k + h_size[0] * (k - 1);
          aa = std::abs(h_data[i].re);
          ba = aa + std::abs(h_data[i].im);
          if (ba <= SMLNUM) {
            exitg3 = true;
          } else {
            ix0 = k + h_size[0] * k;
            bb = std::abs(h_data[ix0].re) + std::abs(h_data[ix0].im);
            knt = i - 1;
            tst = (std::abs(h_data[knt].re) + std::abs(h_data[knt].im)) + bb;
            if (tst == 0.0F) {
              if (k - 1 >= 1) {
                tst = std::abs(h_data[(k + h_size[0] * (k - 2)) - 1].re);
              }

              if (k + 2 <= n) {
                tst += std::abs(h_data[ix0 + 1].re);
              }
            }

            if (aa <= 1.1920929E-7F * tst) {
              h_tmp = ix0 - 1;
              tst = std::abs(h_data[h_tmp].re) + std::abs(h_data[h_tmp].im);
              if (ba > tst) {
                ab = ba;
                ba = tst;
              } else {
                ab = tst;
              }

              tst = std::abs(h_data[knt].re - h_data[ix0].re) + std::abs
                (h_data[knt].im - h_data[ix0].im);
              if (bb > tst) {
                aa = bb;
                bb = tst;
              } else {
                aa = tst;
              }

              s = aa + ab;
              aa = 1.1920929E-7F * (bb * (aa / s));
              if (SMLNUM > aa) {
                aa = SMLNUM;
              }

              if (ba * (ab / s) <= aa) {
                exitg3 = true;
              } else {
                k--;
              }
            } else {
              k--;
            }
          }
        }

        L = k - 1;
        if (k + 1 > 1) {
          i = k + h_size[0] * (k - 1);
          h_data[i].re = 0.0F;
          h_data[i].im = 0.0F;
        }

        if (k + 1 >= b_i + 1) {
          goto140 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            h_tmp = k + h_size[0] * k;
            t_re = 0.75F * std::abs(h_data[(k + h_size[0] * k) + 1].re) +
              h_data[h_tmp].re;
            ab = h_data[h_tmp].im;
          } else if (its == 20) {
            h_tmp = b_i + h_size[0] * b_i;
            t_re = 0.75F * std::abs(h_data[b_i + h_size[0] * (b_i - 1)].re) +
              h_data[h_tmp].re;
            ab = h_data[h_tmp].im;
          } else {
            h_tmp = b_i + h_size[0] * b_i;
            t_re = h_data[h_tmp].re;
            ab = h_data[h_tmp].im;
            v = h_data[h_tmp - 1];
            c_sqrt(&v);
            ix0 = b_i + h_size[0] * (b_i - 1);
            sc = h_data[ix0];
            c_sqrt(&sc);
            u_re = v.re * sc.re - v.im * sc.im;
            u_im = v.re * sc.im + v.im * sc.re;
            s = std::abs(u_re) + std::abs(u_im);
            if (s != 0.0F) {
              knt = ix0 - 1;
              ba = 0.5F * (h_data[knt].re - h_data[h_tmp].re);
              im = 0.5F * (h_data[knt].im - h_data[h_tmp].im);
              bb = std::abs(ba) + std::abs(im);
              if (s <= bb) {
                s = bb;
              }

              if (im == 0.0F) {
                t_re = ba / s;
                ab = 0.0F;
              } else if (ba == 0.0F) {
                t_re = 0.0F;
                ab = im / s;
              } else {
                t_re = ba / s;
                ab = im / s;
              }

              re = t_re * t_re - ab * ab;
              tst = t_re * ab;
              if (u_im == 0.0F) {
                sc.re = u_re / s;
                sc.im = 0.0F;
              } else if (u_re == 0.0F) {
                sc.re = 0.0F;
                sc.im = u_im / s;
              } else {
                sc.re = u_re / s;
                sc.im = u_im / s;
              }

              aa = sc.re * sc.re - sc.im * sc.im;
              ab = sc.re * sc.im;
              v.re = re + aa;
              v.im = (tst + tst) + (ab + ab);
              c_sqrt(&v);
              sc.re = s * v.re;
              sc.im = s * v.im;
              if (bb > 0.0F) {
                if (im == 0.0F) {
                  t_re = ba / bb;
                  ab = 0.0F;
                } else if (ba == 0.0F) {
                  t_re = 0.0F;
                  ab = im / bb;
                } else {
                  t_re = ba / bb;
                  ab = im / bb;
                }

                if (t_re * sc.re + ab * sc.im < 0.0F) {
                  sc.re = -sc.re;
                  sc.im = -sc.im;
                }
              }

              br = ba + sc.re;
              ab = im + sc.im;
              if (ab == 0.0F) {
                if (u_im == 0.0F) {
                  ba = u_re / br;
                  tst = 0.0F;
                } else if (u_re == 0.0F) {
                  ba = 0.0F;
                  tst = u_im / br;
                } else {
                  ba = u_re / br;
                  tst = u_im / br;
                }
              } else if (br == 0.0F) {
                if (u_re == 0.0F) {
                  ba = u_im / ab;
                  tst = 0.0F;
                } else if (u_im == 0.0F) {
                  ba = 0.0F;
                  tst = -(u_re / ab);
                } else {
                  ba = u_im / ab;
                  tst = -(u_re / ab);
                }
              } else {
                bb = std::abs(br);
                tst = std::abs(ab);
                if (bb > tst) {
                  s = ab / br;
                  tst = br + s * ab;
                  ba = (u_re + s * u_im) / tst;
                  tst = (u_im - s * u_re) / tst;
                } else if (tst == bb) {
                  if (br > 0.0F) {
                    aa = 0.5F;
                  } else {
                    aa = -0.5F;
                  }

                  if (ab > 0.0F) {
                    tst = 0.5F;
                  } else {
                    tst = -0.5F;
                  }

                  ba = (u_re * aa + u_im * tst) / bb;
                  tst = (u_im * aa - u_re * tst) / bb;
                } else {
                  s = br / ab;
                  tst = ab + s * br;
                  ba = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h_data[h_tmp].re - (u_re * ba - u_im * tst);
              ab = h_data[h_tmp].im - (u_re * tst + u_im * ba);
            }
          }

          goto70 = false;
          m = b_i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            h_tmp = m + h_size[0] * (m - 1);
            ix0 = h_tmp - 1;
            sc.re = h_data[ix0].re - t_re;
            sc.im = h_data[ix0].im - ab;
            tst = h_data[h_tmp].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            if (sc.im == 0.0F) {
              re = sc.re / s;
              im = 0.0F;
            } else if (sc.re == 0.0F) {
              re = 0.0F;
              im = sc.im / s;
            } else {
              re = sc.re / s;
              im = sc.im / s;
            }

            sc.re = re;
            sc.im = im;
            tst /= s;
            b_v[0] = sc;
            b_v[1].re = tst;
            b_v[1].im = 0.0F;
            i = m + h_size[0] * m;
            if (std::abs(h_data[(m + h_size[0] * (m - 2)) - 1].re) * std::abs
                (tst) <= 1.1920929E-7F * ((std::abs(re) + std::abs(im)) * ((std::
                   abs(h_data[ix0].re) + std::abs(h_data[ix0].im)) + (std::abs
                   (h_data[i].re) + std::abs(h_data[i].im))))) {
              goto70 = true;
              exitg3 = true;
            } else {
              m--;
            }
          }

          if (!goto70) {
            ix0 = k + h_size[0] * k;
            sc.re = h_data[ix0].re - t_re;
            sc.im = h_data[ix0].im - ab;
            tst = h_data[(k + h_size[0] * k) + 1].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            if (sc.im == 0.0F) {
              b_v[0].re = sc.re / s;
              b_v[0].im = 0.0F;
            } else if (sc.re == 0.0F) {
              b_v[0].re = 0.0F;
              b_v[0].im = sc.im / s;
            } else {
              b_v[0].re = sc.re / s;
              b_v[0].im = sc.im / s;
            }

            tst /= s;
            b_v[1].re = tst;
            b_v[1].im = 0.0F;
          }

          for (b_k = m; b_k <= b_i; b_k++) {
            if (b_k > m) {
              knt = b_k + h_size[0] * (b_k - 2);
              b_v[0] = h_data[knt - 1];
              b_v[1] = h_data[knt];
            }

            ba = b_v[1].re;
            im = b_v[1].im;
            sc = b_v[0];
            t_re = 0.0F;
            ab = 0.0F;
            tst = rt_hypotf(b_v[1].re, b_v[1].im);
            if ((tst != 0.0F) || (b_v[0].im != 0.0F)) {
              aa = b_xdlapy3(b_v[0].re, b_v[0].im, tst);
              if (b_v[0].re >= 0.0F) {
                aa = -aa;
              }

              if (std::abs(aa) < 9.86076132E-32F) {
                knt = -1;
                do {
                  knt++;
                  ba *= 1.01412048E+31F;
                  im *= 1.01412048E+31F;
                  aa *= 1.01412048E+31F;
                  sc.re *= 1.01412048E+31F;
                  sc.im *= 1.01412048E+31F;
                } while (!(std::abs(aa) >= 9.86076132E-32F));

                aa = b_xdlapy3(sc.re, sc.im, rt_hypotf(ba, im));
                if (sc.re >= 0.0F) {
                  aa = -aa;
                }

                tst = aa - sc.re;
                if (0.0F - sc.im == 0.0F) {
                  t_re = tst / aa;
                  ab = 0.0F;
                } else if (tst == 0.0F) {
                  t_re = 0.0F;
                  ab = (0.0F - sc.im) / aa;
                } else {
                  t_re = tst / aa;
                  ab = (0.0F - sc.im) / aa;
                }

                v.re = sc.re - aa;
                v.im = sc.im;
                sc = b_recip(v);
                re = sc.re * ba - sc.im * im;
                im = sc.re * im + sc.im * ba;
                ba = re;
                for (h_tmp = 0; h_tmp <= knt; h_tmp++) {
                  aa *= 9.86076132E-32F;
                }

                sc.re = aa;
                sc.im = 0.0F;
              } else {
                tst = aa - b_v[0].re;
                if (0.0F - b_v[0].im == 0.0F) {
                  t_re = tst / aa;
                  ab = 0.0F;
                } else if (tst == 0.0F) {
                  t_re = 0.0F;
                  ab = (0.0F - b_v[0].im) / aa;
                } else {
                  t_re = tst / aa;
                  ab = (0.0F - b_v[0].im) / aa;
                }

                v.re = b_v[0].re - aa;
                v.im = b_v[0].im;
                v = b_recip(v);
                ba = v.re * b_v[1].re - v.im * b_v[1].im;
                im = v.re * b_v[1].im + v.im * b_v[1].re;
                sc.re = aa;
                sc.im = 0.0F;
              }
            }

            b_v[0] = sc;
            b_v[1].re = ba;
            b_v[1].im = im;
            if (b_k > m) {
              h_data[(b_k + h_size[0] * (b_k - 2)) - 1] = sc;
              i = b_k + h_size[0] * (b_k - 2);
              h_data[i].re = 0.0F;
              h_data[i].im = 0.0F;
            }

            tst = t_re * ba - ab * im;
            for (j = b_k; j <= n; j++) {
              h_tmp = b_k + h_size[0] * (j - 1);
              ix0 = h_tmp - 1;
              sc.re = (t_re * h_data[ix0].re - -ab * h_data[ix0].im) + tst *
                h_data[h_tmp].re;
              sc.im = (t_re * h_data[ix0].im + -ab * h_data[ix0].re) + tst *
                h_data[h_tmp].im;
              h_data[ix0].re -= sc.re;
              h_data[ix0].im -= sc.im;
              h_data[h_tmp].re -= sc.re * ba - sc.im * im;
              h_data[h_tmp].im -= sc.re * im + sc.im * ba;
            }

            if (b_k + 2 < b_i + 1) {
              i = b_k + 1;
            } else {
              i = b_i;
            }

            for (j = 0; j <= i; j++) {
              ix0 = j + h_size[0] * (b_k - 1);
              h_tmp = j + h_size[0] * b_k;
              sc.re = (t_re * h_data[ix0].re - ab * h_data[ix0].im) + tst *
                h_data[h_tmp].re;
              sc.im = (t_re * h_data[ix0].im + ab * h_data[ix0].re) + tst *
                h_data[h_tmp].im;
              h_data[ix0].re -= sc.re;
              h_data[ix0].im -= sc.im;
              h_data[h_tmp].re -= sc.re * ba - sc.im * -im;
              h_data[h_tmp].im -= sc.re * -im + sc.im * ba;
            }

            if ((b_k == m) && (m > k + 1)) {
              t_re = 1.0F - t_re;
              ab = 0.0F - ab;
              br = rt_hypotf(t_re, ab);
              if (ab == 0.0F) {
                re = t_re / br;
                im = 0.0F;
              } else if (t_re == 0.0F) {
                re = 0.0F;
                im = ab / br;
              } else {
                re = t_re / br;
                im = ab / br;
              }

              knt = m + h_size[0] * (m - 1);
              aa = h_data[knt].re * -im + h_data[knt].im * re;
              h_data[knt].re = h_data[knt].re * re - h_data[knt].im * -im;
              h_data[knt].im = aa;
              if (m + 2 <= b_i + 1) {
                knt = (m + h_size[0] * m) + 1;
                aa = h_data[knt].re * im + h_data[knt].im * re;
                h_data[knt].re = h_data[knt].re * re - h_data[knt].im * im;
                h_data[knt].im = aa;
              }

              for (j = m; j <= b_i + 1; j++) {
                if (j != m + 1) {
                  if (n > j) {
                    ix0 = j + j * ldh;
                    i = ix0 + ldh * ((n - j) - 1);
                    for (h_tmp = ix0; ldh < 0 ? h_tmp >= i : h_tmp <= i; h_tmp +=
                         ldh) {
                      aa = re * h_data[h_tmp - 1].im + im * h_data[h_tmp - 1].re;
                      h_data[h_tmp - 1].re = re * h_data[h_tmp - 1].re - im *
                        h_data[h_tmp - 1].im;
                      h_data[h_tmp - 1].im = aa;
                    }
                  }

                  h_tmp = (j - 1) * ldh;
                  ix0 = h_tmp + 1;
                  i = (h_tmp + j) - 1;
                  for (h_tmp = ix0; h_tmp <= i; h_tmp++) {
                    aa = re * h_data[h_tmp - 1].im + -im * h_data[h_tmp - 1].re;
                    h_data[h_tmp - 1].re = re * h_data[h_tmp - 1].re - -im *
                      h_data[h_tmp - 1].im;
                    h_data[h_tmp - 1].im = aa;
                  }
                }
              }
            }
          }

          h_tmp = b_i + h_size[0] * (b_i - 1);
          t_re = h_data[h_tmp].re;
          ab = h_data[h_tmp].im;
          if (h_data[h_tmp].im != 0.0F) {
            aa = rt_hypotf(h_data[h_tmp].re, h_data[h_tmp].im);
            h_data[h_tmp].re = aa;
            h_data[h_tmp].im = 0.0F;
            if (ab == 0.0F) {
              re = t_re / aa;
              im = 0.0F;
            } else if (t_re == 0.0F) {
              re = 0.0F;
              im = ab / aa;
            } else {
              re = t_re / aa;
              im = ab / aa;
            }

            if (n > b_i + 1) {
              ix0 = (b_i + (b_i + 1) * ldh) + 1;
              i = ix0 + ldh * ((n - b_i) - 2);
              for (k = ix0; ldh < 0 ? k >= i : k <= i; k += ldh) {
                aa = re * h_data[k - 1].im + -im * h_data[k - 1].re;
                h_data[k - 1].re = re * h_data[k - 1].re - -im * h_data[k - 1].
                  im;
                h_data[k - 1].im = aa;
              }
            }

            h_tmp = b_i * ldh;
            ix0 = h_tmp + 1;
            i = h_tmp + b_i;
            for (k = ix0; k <= i; k++) {
              aa = re * h_data[k - 1].im + im * h_data[k - 1].re;
              h_data[k - 1].re = re * h_data[k - 1].re - im * h_data[k - 1].im;
              h_data[k - 1].im = aa;
            }
          }

          its++;
        }
      }

      if (!goto140) {
        info = b_i + 1;
        exitg1 = true;
      } else {
        b_i = L;
      }
    }
  }

  return info;
}

//
// function G = equations_for_groebner(F)
// Arguments    : const float F[72]
//                float G[112]
// Return Type  : void
//
static void b_equations_for_groebner(const float F[72], float G[112])
{
  int i;
  float b_F[6];
  float c_F[6];
  double dv[15];
  int M_tmp;
  float M[24];
  double dv1[15];
  int b_M_tmp;
  double dv2[28];
  double dv3[28];
  double dv4[28];
  int i1;
  float b_M[54];

  // 'equations_for_groebner:2' G = zeros(4, 1, 28, 'like', F);
  // 'equations_for_groebner:3' G(1, 1, :) = find_det3(F(2:4, :, :));
  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 10];
    c_F[i] = F[12 * i + 7];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 1];
  }

  b_mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 5];
  }

  b_mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 2];
    c_F[i] = F[12 * i + 7];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 6];
    c_F[i] = F[12 * i + 3];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 9];
  }

  b_mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[i << 2] = static_cast<float>(((dv2[i] - dv3[i]) + dv4[i]));
  }

  // 'equations_for_groebner:4' G(2, 1, :) = find_det3([F(1, :, :); F(3:4, :, :)]); 
  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * i;
      b_M_tmp = 3 * i1 + 9 * i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 2];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 7];
    c_F[i] = b_M[9 * i + 5];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i];
  }

  b_mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 3];
  }

  b_mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 1];
    c_F[i] = b_M[9 * i + 5];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 4];
    c_F[i] = b_M[9 * i + 2];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 6];
  }

  b_mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[(i << 2) + 1] = static_cast<float>(((dv2[i] - dv3[i]) + dv4[i]));
  }

  // 'equations_for_groebner:5' G(3, 1, :) = find_det3([F(1:2, :, :); F(4, :, :)]); 
  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * i;
      b_M_tmp = 3 * i1 + 9 * i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 1];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 7];
    c_F[i] = b_M[9 * i + 5];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i];
  }

  b_mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 3];
  }

  b_mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 1];
    c_F[i] = b_M[9 * i + 5];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 4];
    c_F[i] = b_M[9 * i + 2];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = b_M[9 * i + 6];
  }

  b_mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[(i << 2) + 2] = static_cast<float>(((dv2[i] - dv3[i]) + dv4[i]));
  }

  // 'equations_for_groebner:6' G(4, 1, :) = find_det3(F(1:3, :, :));
  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 9];
    c_F[i] = F[12 * i + 6];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i];
  }

  b_mult_poly42(dv, b_F, dv2);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp];
    c_F[i] = M[M_tmp + 3];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    M_tmp = i << 2;
    b_F[i] = M[M_tmp + 2];
    c_F[i] = M[M_tmp + 1];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 4];
  }

  b_mult_poly42(dv, b_F, dv3);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 1];
    c_F[i] = F[12 * i + 6];
  }

  b_mult_poly22(b_F, c_F, dv);
  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 5];
    c_F[i] = F[12 * i + 2];
  }

  b_mult_poly22(b_F, c_F, dv1);
  for (i = 0; i < 15; i++) {
    dv[i] -= dv1[i];
  }

  for (i = 0; i < 6; i++) {
    b_F[i] = F[12 * i + 8];
  }

  b_mult_poly42(dv, b_F, dv4);
  for (i = 0; i < 28; i++) {
    G[(i << 2) + 3] = static_cast<float>(((dv2[i] - dv3[i]) + dv4[i]));
  }

  // 'equations_for_groebner:7' G = squeeze(G);
}

//
// function [fc, fs, n] = find_f(F, x, y, e)
// Arguments    : const float F[72]
//                float x
//                float y
//                float e
//                float fc_data[]
//                int fc_size[2]
//                float fs_data[]
//                int fs_size[2]
//                double *n
// Return Type  : void
//
static void b_find_f(const float F[72], float x, float y, float e, float
                     fc_data[], int fc_size[2], float fs_data[], int fs_size[2],
                     double *n)
{
  int i;
  float F_eval[12];
  float mons[6];
  int b_i;
  float b_F_eval[12];
  float Q[9];
  int j;
  float fc_tmp;
  int k;

  // 'find_f:2' F_eval = zeros(4, 3, 'like', F);
  for (i = 0; i < 12; i++) {
    F_eval[i] = 0.0F;
  }

  // 'find_f:3' mons =  [x^2, x*y, y^2, x, y, 1];
  mons[0] = x * x;
  mons[1] = x * y;
  mons[2] = y * y;
  mons[3] = x;
  mons[4] = y;
  mons[5] = 1.0F;

  // 'find_f:4' for i = 1 : 4
  // 'find_f:12' [Q, R] = qr(F_eval');
  for (b_i = 0; b_i < 4; b_i++) {
    // 'find_f:5' for j = 1 : 3
    for (j = 0; j < 3; j++) {
      // 'find_f:6' for k = 1 : 6
      i = b_i + (j << 2);
      fc_tmp = F_eval[i];
      for (k = 0; k < 6; k++) {
        // 'find_f:7' F_eval(i, j) = F_eval(i, j) + F(i, j, k)*mons(k);
        fc_tmp += F[i + 12 * k] * mons[k];
      }

      F_eval[i] = fc_tmp;
      b_F_eval[j + 3 * b_i] = fc_tmp;
    }
  }

  b_qr(b_F_eval, Q, F_eval);

  // 'find_f:13' if abs(R(2, 2) - 0) < e
  if (std::abs(F_eval[4]) < e) {
    // rankF = 1
    // 'find_f:14' n = 2;
    *n = 2.0;

    // 'find_f:15' fc = [Q(1, 2) / Q(3, 2), Q(1, 2) / Q(3, 2)];
    fc_tmp = Q[3] / Q[5];
    fc_size[0] = 1;
    fc_size[1] = 2;
    fc_data[0] = fc_tmp;
    fc_data[1] = fc_tmp;

    // 'find_f:16' fs = [Q(2, 2) / Q(3, 2), Q(2, 3) / Q(3, 3)];
    fs_size[0] = 1;
    fs_size[1] = 2;
    fs_data[0] = Q[4] / Q[5];
    fs_data[1] = Q[7] / Q[8];
  } else {
    // 'find_f:17' else
    // 'find_f:18' n = 1;
    *n = 1.0;

    // 'find_f:19' fc = Q(1, 3) / Q(3, 3);
    fc_size[0] = 1;
    fc_size[1] = 1;
    fc_data[0] = Q[6] / Q[8];

    // 'find_f:20' fs = Q(2, 3) / Q(3, 3);
    fs_size[0] = 1;
    fs_size[1] = 1;
    fs_data[0] = Q[7] / Q[8];
  }
}

//
// Arguments    : const double A[400]
//                double B[200]
// Return Type  : void
//
static void b_mldivide(const double A[400], double B[200])
{
  double b_A[400];
  int ipiv[20];
  int info;
  int jBcol;
  int temp_tmp;
  double temp;
  int kAcol;
  int i;
  int i1;
  int b_i;
  int i2;
  std::memcpy(&b_A[0], &A[0], 400U * sizeof(double));
  xzgetrf(b_A, ipiv, &info);
  for (info = 0; info < 19; info++) {
    if (ipiv[info] != info + 1) {
      for (jBcol = 0; jBcol < 10; jBcol++) {
        temp_tmp = info + 20 * jBcol;
        temp = B[temp_tmp];
        i = (ipiv[info] + 20 * jBcol) - 1;
        B[temp_tmp] = B[i];
        B[i] = temp;
      }
    }
  }

  for (info = 0; info < 10; info++) {
    jBcol = 20 * info;
    for (temp_tmp = 0; temp_tmp < 20; temp_tmp++) {
      kAcol = 20 * temp_tmp;
      i = temp_tmp + jBcol;
      if (B[i] != 0.0) {
        i1 = temp_tmp + 2;
        for (b_i = i1; b_i < 21; b_i++) {
          i2 = (b_i + jBcol) - 1;
          B[i2] -= B[i] * b_A[(b_i + kAcol) - 1];
        }
      }
    }
  }

  for (info = 0; info < 10; info++) {
    jBcol = 20 * info;
    for (temp_tmp = 19; temp_tmp >= 0; temp_tmp--) {
      kAcol = 20 * temp_tmp;
      i = temp_tmp + jBcol;
      if (B[i] != 0.0) {
        B[i] /= b_A[temp_tmp + kAcol];
        for (b_i = 0; b_i < temp_tmp; b_i++) {
          i1 = b_i + jBcol;
          B[i1] -= B[i] * b_A[b_i + kAcol];
        }
      }
    }
  }
}

//
// function G20 = mult_for_groebner(G4)
// Arguments    : const float G4[112]
//                float G20[900]
// Return Type  : void
//
static void b_mult_for_groebner(const float G4[112], float G20[900])
{
  int i;
  float G20_tmp;
  float b_G20_tmp;
  float c_G20_tmp;
  float d_G20_tmp;
  float e_G20_tmp;
  float f_G20_tmp;
  float g_G20_tmp;
  float h_G20_tmp;
  float i_G20_tmp;
  float j_G20_tmp;
  float k_G20_tmp;
  float l_G20_tmp;
  float m_G20_tmp;
  float n_G20_tmp;
  float o_G20_tmp;
  float p_G20_tmp;
  float q_G20_tmp;
  float r_G20_tmp;
  float s_G20_tmp;
  float t_G20_tmp;
  float u_G20_tmp;
  float v_G20_tmp;
  float w_G20_tmp;
  float x_G20_tmp;
  float y_G20_tmp;
  float ab_G20_tmp;
  float bb_G20_tmp;
  int b_i;

  // 'mult_for_groebner:2' G20 = zeros(20, 45, 'like', G4);
  std::memset(&G20[0], 0, 900U * sizeof(float));

  // multiply by qx^2
  // 'mult_for_groebner:5' for i = 1 : 4
  // multiply by qxqy
  // 'mult_for_groebner:37' for i = 1 : 4
  // multiply by qx
  // 'mult_for_groebner:69' for i = 1 : 4
  // multiply by qy
  // 'mult_for_groebner:101' for i = 1 : 4
  // multiply by 1
  // 'mult_for_groebner:133' for i = 1 : 4
  for (i = 0; i < 4; i++) {
    // 'mult_for_groebner:6' G20(i, 1) = G4(i, 1);
    G20[i] = G4[i];

    // 'mult_for_groebner:7' G20(i, 2) = G4(i, 2);
    G20_tmp = G4[i + 4];
    G20[i + 20] = G20_tmp;

    // 'mult_for_groebner:8' G20(i, 3) = G4(i, 3);
    b_G20_tmp = G4[i + 8];
    G20[i + 40] = b_G20_tmp;

    // 'mult_for_groebner:9' G20(i, 4) = G4(i, 4);
    c_G20_tmp = G4[i + 12];
    G20[i + 60] = c_G20_tmp;

    // 'mult_for_groebner:10' G20(i, 5) = G4(i, 5);
    d_G20_tmp = G4[i + 16];
    G20[i + 80] = d_G20_tmp;

    // 'mult_for_groebner:11' G20(i, 6) = G4(i, 6);
    e_G20_tmp = G4[i + 20];
    G20[i + 100] = e_G20_tmp;

    // 'mult_for_groebner:12' G20(i, 7) = G4(i, 7);
    f_G20_tmp = G4[i + 24];
    G20[i + 120] = f_G20_tmp;

    // 'mult_for_groebner:13' G20(i, 10) = G4(i, 8);
    g_G20_tmp = G4[i + 28];
    G20[i + 180] = g_G20_tmp;

    // 'mult_for_groebner:14' G20(i, 11) = G4(i, 9);
    h_G20_tmp = G4[i + 32];
    G20[i + 200] = h_G20_tmp;

    // 'mult_for_groebner:15' G20(i, 12) = G4(i, 10);
    i_G20_tmp = G4[i + 36];
    G20[i + 220] = i_G20_tmp;

    // 'mult_for_groebner:16' G20(i, 13) = G4(i, 11);
    j_G20_tmp = G4[i + 40];
    G20[i + 240] = j_G20_tmp;

    // 'mult_for_groebner:17' G20(i, 14) = G4(i, 12);
    k_G20_tmp = G4[i + 44];
    G20[i + 260] = k_G20_tmp;

    // 'mult_for_groebner:18' G20(i, 15) = G4(i, 13);
    l_G20_tmp = G4[i + 48];
    G20[i + 280] = l_G20_tmp;

    // 'mult_for_groebner:19' G20(i, 18) = G4(i, 14);
    m_G20_tmp = G4[i + 52];
    G20[i + 340] = m_G20_tmp;

    // 'mult_for_groebner:20' G20(i, 19) = G4(i, 15);
    n_G20_tmp = G4[i + 56];
    G20[i + 360] = n_G20_tmp;

    // 'mult_for_groebner:21' G20(i, 20) = G4(i, 16);
    o_G20_tmp = G4[i + 60];
    G20[i + 380] = o_G20_tmp;

    // 'mult_for_groebner:22' G20(i, 21) = G4(i, 17);
    p_G20_tmp = G4[i + 64];
    G20[i + 400] = p_G20_tmp;

    // 'mult_for_groebner:23' G20(i, 22) = G4(i, 18);
    q_G20_tmp = G4[i + 68];
    G20[i + 420] = q_G20_tmp;

    // 'mult_for_groebner:24' G20(i, 25) = G4(i, 19);
    r_G20_tmp = G4[i + 72];
    G20[i + 480] = r_G20_tmp;

    // 'mult_for_groebner:25' G20(i, 26) = G4(i, 20);
    s_G20_tmp = G4[i + 76];
    G20[i + 500] = s_G20_tmp;

    // 'mult_for_groebner:26' G20(i, 27) = G4(i, 21);
    t_G20_tmp = G4[i + 80];
    G20[i + 520] = t_G20_tmp;

    // 'mult_for_groebner:27' G20(i, 28) = G4(i, 22);
    u_G20_tmp = G4[i + 84];
    G20[i + 540] = u_G20_tmp;

    // 'mult_for_groebner:28' G20(i, 31) = G4(i, 23);
    v_G20_tmp = G4[i + 88];
    G20[i + 600] = v_G20_tmp;

    // 'mult_for_groebner:29' G20(i, 32) = G4(i, 24);
    w_G20_tmp = G4[i + 92];
    G20[i + 620] = w_G20_tmp;

    // 'mult_for_groebner:30' G20(i, 33) = G4(i, 25);
    x_G20_tmp = G4[i + 96];
    G20[i + 640] = x_G20_tmp;

    // 'mult_for_groebner:31' G20(i, 36) = G4(i, 26);
    y_G20_tmp = G4[i + 100];
    G20[i + 700] = y_G20_tmp;

    // 'mult_for_groebner:32' G20(i, 37) = G4(i, 27);
    ab_G20_tmp = G4[i + 104];
    G20[i + 720] = ab_G20_tmp;

    // 'mult_for_groebner:33' G20(i, 40) = G4(i, 28);
    bb_G20_tmp = G4[i + 108];
    G20[i + 780] = bb_G20_tmp;

    // 'mult_for_groebner:38' G20(i + 4, 2) = G4(i, 1);
    G20[i + 24] = G4[i];

    // 'mult_for_groebner:39' G20(i + 4, 3) = G4(i, 2);
    G20[i + 44] = G20_tmp;

    // 'mult_for_groebner:40' G20(i + 4, 4) = G4(i, 3);
    G20[i + 64] = b_G20_tmp;

    // 'mult_for_groebner:41' G20(i + 4, 5) = G4(i, 4);
    G20[i + 84] = c_G20_tmp;

    // 'mult_for_groebner:42' G20(i + 4, 6) = G4(i, 5);
    G20[i + 104] = d_G20_tmp;

    // 'mult_for_groebner:43' G20(i + 4, 7) = G4(i, 6);
    G20[i + 124] = e_G20_tmp;

    // 'mult_for_groebner:44' G20(i + 4, 8) = G4(i, 7);
    G20[i + 144] = f_G20_tmp;

    // 'mult_for_groebner:45' G20(i + 4, 11) = G4(i, 8);
    G20[i + 204] = g_G20_tmp;

    // 'mult_for_groebner:46' G20(i + 4, 12) = G4(i, 9);
    G20[i + 224] = h_G20_tmp;

    // 'mult_for_groebner:47' G20(i + 4, 13) = G4(i, 10);
    G20[i + 244] = i_G20_tmp;

    // 'mult_for_groebner:48' G20(i + 4, 14) = G4(i, 11);
    G20[i + 264] = j_G20_tmp;

    // 'mult_for_groebner:49' G20(i + 4, 15) = G4(i, 12);
    G20[i + 284] = k_G20_tmp;

    // 'mult_for_groebner:50' G20(i + 4, 16) = G4(i, 13);
    G20[i + 304] = l_G20_tmp;

    // 'mult_for_groebner:51' G20(i + 4, 19) = G4(i, 14);
    G20[i + 364] = m_G20_tmp;

    // 'mult_for_groebner:52' G20(i + 4, 20) = G4(i, 15);
    G20[i + 384] = n_G20_tmp;

    // 'mult_for_groebner:53' G20(i + 4, 21) = G4(i, 16);
    G20[i + 404] = o_G20_tmp;

    // 'mult_for_groebner:54' G20(i + 4, 22) = G4(i, 17);
    G20[i + 424] = p_G20_tmp;

    // 'mult_for_groebner:55' G20(i + 4, 23) = G4(i, 18);
    G20[i + 444] = q_G20_tmp;

    // 'mult_for_groebner:56' G20(i + 4, 26) = G4(i, 19);
    G20[i + 504] = r_G20_tmp;

    // 'mult_for_groebner:57' G20(i + 4, 27) = G4(i, 20);
    G20[i + 524] = s_G20_tmp;

    // 'mult_for_groebner:58' G20(i + 4, 28) = G4(i, 21);
    G20[i + 544] = t_G20_tmp;

    // 'mult_for_groebner:59' G20(i + 4, 29) = G4(i, 22);
    G20[i + 564] = u_G20_tmp;

    // 'mult_for_groebner:60' G20(i + 4, 32) = G4(i, 23);
    G20[i + 624] = v_G20_tmp;

    // 'mult_for_groebner:61' G20(i + 4, 33) = G4(i, 24);
    G20[i + 644] = w_G20_tmp;

    // 'mult_for_groebner:62' G20(i + 4, 34) = G4(i, 25);
    G20[i + 664] = x_G20_tmp;

    // 'mult_for_groebner:63' G20(i + 4, 37) = G4(i, 26);
    G20[i + 724] = y_G20_tmp;

    // 'mult_for_groebner:64' G20(i + 4, 38) = G4(i, 27);
    G20[i + 744] = ab_G20_tmp;

    // 'mult_for_groebner:65' G20(i + 4, 41) = G4(i, 28);
    G20[i + 804] = bb_G20_tmp;

    // 'mult_for_groebner:70' G20(i + 8, 10) = G4(i, 1);
    G20[i + 188] = G4[i];

    // 'mult_for_groebner:71' G20(i + 8, 11) = G4(i, 2);
    G20[i + 208] = G20_tmp;

    // 'mult_for_groebner:72' G20(i + 8, 12) = G4(i, 3);
    G20[i + 228] = b_G20_tmp;

    // 'mult_for_groebner:73' G20(i + 8, 13) = G4(i, 4);
    G20[i + 248] = c_G20_tmp;

    // 'mult_for_groebner:74' G20(i + 8, 14) = G4(i, 5);
    G20[i + 268] = d_G20_tmp;

    // 'mult_for_groebner:75' G20(i + 8, 15) = G4(i, 6);
    G20[i + 288] = e_G20_tmp;

    // 'mult_for_groebner:76' G20(i + 8, 16) = G4(i, 7);
    G20[i + 308] = f_G20_tmp;

    // 'mult_for_groebner:77' G20(i + 8, 18) = G4(i, 8);
    G20[i + 348] = g_G20_tmp;

    // 'mult_for_groebner:78' G20(i + 8, 19) = G4(i, 9);
    G20[i + 368] = h_G20_tmp;

    // 'mult_for_groebner:79' G20(i + 8, 20) = G4(i, 10);
    G20[i + 388] = i_G20_tmp;

    // 'mult_for_groebner:80' G20(i + 8, 21) = G4(i, 11);
    G20[i + 408] = j_G20_tmp;

    // 'mult_for_groebner:81' G20(i + 8, 22) = G4(i, 12);
    G20[i + 428] = k_G20_tmp;

    // 'mult_for_groebner:82' G20(i + 8, 23) = G4(i, 13);
    G20[i + 448] = l_G20_tmp;

    // 'mult_for_groebner:83' G20(i + 8, 25) = G4(i, 14);
    G20[i + 488] = m_G20_tmp;

    // 'mult_for_groebner:84' G20(i + 8, 26) = G4(i, 15);
    G20[i + 508] = n_G20_tmp;

    // 'mult_for_groebner:85' G20(i + 8, 27) = G4(i, 16);
    G20[i + 528] = o_G20_tmp;

    // 'mult_for_groebner:86' G20(i + 8, 28) = G4(i, 17);
    G20[i + 548] = p_G20_tmp;

    // 'mult_for_groebner:87' G20(i + 8, 29) = G4(i, 18);
    G20[i + 568] = q_G20_tmp;

    // 'mult_for_groebner:88' G20(i + 8, 31) = G4(i, 19);
    G20[i + 608] = r_G20_tmp;

    // 'mult_for_groebner:89' G20(i + 8, 32) = G4(i, 20);
    G20[i + 628] = s_G20_tmp;

    // 'mult_for_groebner:90' G20(i + 8, 33) = G4(i, 21);
    G20[i + 648] = t_G20_tmp;

    // 'mult_for_groebner:91' G20(i + 8, 34) = G4(i, 22);
    G20[i + 668] = u_G20_tmp;

    // 'mult_for_groebner:92' G20(i + 8, 36) = G4(i, 23);
    G20[i + 708] = v_G20_tmp;

    // 'mult_for_groebner:93' G20(i + 8, 37) = G4(i, 24);
    G20[i + 728] = w_G20_tmp;

    // 'mult_for_groebner:94' G20(i + 8, 38) = G4(i, 25);
    G20[i + 748] = x_G20_tmp;

    // 'mult_for_groebner:95' G20(i + 8, 40) = G4(i, 26);
    G20[i + 788] = y_G20_tmp;

    // 'mult_for_groebner:96' G20(i + 8, 41) = G4(i, 27);
    G20[i + 808] = ab_G20_tmp;

    // 'mult_for_groebner:97' G20(i + 8, 43) = G4(i, 28);
    G20[i + 848] = bb_G20_tmp;

    // 'mult_for_groebner:102' G20(i + 12, 11) = G4(i, 1);
    G20[i + 212] = G4[i];

    // 'mult_for_groebner:103' G20(i + 12, 12) = G4(i, 2);
    G20[i + 232] = G20_tmp;

    // 'mult_for_groebner:104' G20(i + 12, 13) = G4(i, 3);
    G20[i + 252] = b_G20_tmp;

    // 'mult_for_groebner:105' G20(i + 12, 14) = G4(i, 4);
    G20[i + 272] = c_G20_tmp;

    // 'mult_for_groebner:106' G20(i + 12, 15) = G4(i, 5);
    G20[i + 292] = d_G20_tmp;

    // 'mult_for_groebner:107' G20(i + 12, 16) = G4(i, 6);
    G20[i + 312] = e_G20_tmp;

    // 'mult_for_groebner:108' G20(i + 12, 17) = G4(i, 7);
    G20[i + 332] = f_G20_tmp;

    // 'mult_for_groebner:109' G20(i + 12, 19) = G4(i, 8);
    G20[i + 372] = g_G20_tmp;

    // 'mult_for_groebner:110' G20(i + 12, 20) = G4(i, 9);
    G20[i + 392] = h_G20_tmp;

    // 'mult_for_groebner:111' G20(i + 12, 21) = G4(i, 10);
    G20[i + 412] = i_G20_tmp;

    // 'mult_for_groebner:112' G20(i + 12, 22) = G4(i, 11);
    G20[i + 432] = j_G20_tmp;

    // 'mult_for_groebner:113' G20(i + 12, 23) = G4(i, 12);
    G20[i + 452] = k_G20_tmp;

    // 'mult_for_groebner:114' G20(i + 12, 24) = G4(i, 13);
    G20[i + 472] = l_G20_tmp;

    // 'mult_for_groebner:115' G20(i + 12, 26) = G4(i, 14);
    G20[i + 512] = m_G20_tmp;

    // 'mult_for_groebner:116' G20(i + 12, 27) = G4(i, 15);
    G20[i + 532] = n_G20_tmp;

    // 'mult_for_groebner:117' G20(i + 12, 28) = G4(i, 16);
    G20[i + 552] = o_G20_tmp;

    // 'mult_for_groebner:118' G20(i + 12, 29) = G4(i, 17);
    G20[i + 572] = p_G20_tmp;

    // 'mult_for_groebner:119' G20(i + 12, 30) = G4(i, 18);
    G20[i + 592] = q_G20_tmp;

    // 'mult_for_groebner:120' G20(i + 12, 32) = G4(i, 19);
    G20[i + 632] = r_G20_tmp;

    // 'mult_for_groebner:121' G20(i + 12, 33) = G4(i, 20);
    G20[i + 652] = s_G20_tmp;

    // 'mult_for_groebner:122' G20(i + 12, 34) = G4(i, 21);
    G20[i + 672] = t_G20_tmp;

    // 'mult_for_groebner:123' G20(i + 12, 35) = G4(i, 22);
    G20[i + 692] = u_G20_tmp;

    // 'mult_for_groebner:124' G20(i + 12, 37) = G4(i, 23);
    G20[i + 732] = v_G20_tmp;

    // 'mult_for_groebner:125' G20(i + 12, 38) = G4(i, 24);
    G20[i + 752] = w_G20_tmp;

    // 'mult_for_groebner:126' G20(i + 12, 39) = G4(i, 25);
    G20[i + 772] = x_G20_tmp;

    // 'mult_for_groebner:127' G20(i + 12, 41) = G4(i, 26);
    G20[i + 812] = y_G20_tmp;

    // 'mult_for_groebner:128' G20(i + 12, 42) = G4(i, 27);
    G20[i + 832] = ab_G20_tmp;

    // 'mult_for_groebner:129' G20(i + 12, 44) = G4(i, 28);
    G20[i + 872] = bb_G20_tmp;

    // 'mult_for_groebner:134' G20(i + 16, 18:end) = G4(i, :);
    for (b_i = 0; b_i < 28; b_i++) {
      G20[(i + 20 * (b_i + 17)) + 16] = G4[i + (b_i << 2)];
    }
  }
}

//
// function c = mult_poly22(a, b)
// Arguments    : const float a[6]
//                const float b[6]
//                double c[15]
// Return Type  : void
//
static void b_mult_poly22(const float a[6], const float b[6], double c[15])
{
  // 'mult_poly22:2' c = zeros(1, 1, 15);
  // 'mult_poly22:3' c(1) = a(:, :, 1)*b(:, :, 1);
  c[0] = a[0] * b[0];

  // x^4
  // 'mult_poly22:4' c(2) = a(:, :, 1)*b(:, :, 2) + a(:, :, 2)*b(:, :, 1);
  c[1] = a[0] * b[1] + a[1] * b[0];

  // x^3*y
  // 'mult_poly22:5' c(3) = a(:, :, 1)*b(:, :, 3) + a(:, :, 2)*b(:, :, 2) + a(:, :, 3)*b(:, :, 1); 
  c[2] = (a[0] * b[2] + a[1] * b[1]) + a[2] * b[0];

  // x^2*y^2
  // 'mult_poly22:6' c(4) = a(:, :, 2)*b(:, :, 3) + a(:, :, 3)*b(:, :, 2);
  c[3] = a[1] * b[2] + a[2] * b[1];

  // x*y^3
  // 'mult_poly22:7' c(5) = a(:, :, 3)*b(:, :, 3);
  c[4] = a[2] * b[2];

  // y^4
  // 'mult_poly22:8' c(6) = a(:, :, 1)*b(:, :, 4) + a(:, :, 4)*b(:, :, 1);
  c[5] = a[0] * b[3] + a[3] * b[0];

  // x^3
  // 'mult_poly22:9' c(7) = a(:, :, 1)*b(:, :, 5) + a(:, :, 2)*b(:, :, 4) + a(:, :, 4)*b(:, :, 2) + a(:, :, 5)*b(:, :, 1); 
  c[6] = ((a[0] * b[4] + a[1] * b[3]) + a[3] * b[1]) + a[4] * b[0];

  // x^2*y
  // 'mult_poly22:10' c(8) = a(:, :, 2)*b(:, :, 5) + a(:, :, 3)*b(:, :, 4) + a(:, :, 4)*b(:, :, 3) + a(:, :, 5)*b(:, :, 2); 
  c[7] = ((a[1] * b[4] + a[2] * b[3]) + a[3] * b[2]) + a[4] * b[1];

  // x*y^2
  // 'mult_poly22:11' c(9) = a(:, :, 3)*b(:, :, 5) + a(:, :, 5)*b(:, :, 3);
  c[8] = a[2] * b[4] + a[4] * b[2];

  // y^3
  // 'mult_poly22:12' c(10) = a(:, :, 1)*b(:, :, 6) + a(:, :, 6)*b(:, :, 1) + a(:, :, 4)*b(:, :, 4); 
  c[9] = (a[0] * b[5] + a[5] * b[0]) + a[3] * b[3];

  // x^2
  // 'mult_poly22:13' c(11) = a(:, :, 2)*b(:, :, 6) + a(:, :, 6)*b(:, :, 2) + a(:, :, 4)*b(:, :, 5) + a(:, :, 5)*b(:, :, 4); 
  c[10] = ((a[1] * b[5] + a[5] * b[1]) + a[3] * b[4]) + a[4] * b[3];

  // x*y
  // 'mult_poly22:14' c(12) = a(:, :, 3)*b(:, :, 6) + a(:, :, 6)*b(:, :, 3) + a(:, :, 5)*b(:, :, 5); 
  c[11] = (a[2] * b[5] + a[5] * b[2]) + a[4] * b[4];

  // y^2
  // 'mult_poly22:15' c(13) = a(:, :, 4)*b(:, :, 6) + a(:, :, 6)*b(:, :, 4);
  c[12] = a[3] * b[5] + a[5] * b[3];

  // x
  // 'mult_poly22:16' c(14) = a(:, :, 5)*b(:, :, 6) + a(:, :, 6)*b(:, :, 5);
  c[13] = a[4] * b[5] + a[5] * b[4];

  // y
  // 'mult_poly22:17' c(15) = a(:, :, 6)*b(:, :, 6);
  c[14] = a[5] * b[5];

  // 1
}

//
// function p = mult_poly42(d, a)
// Arguments    : const double d[15]
//                const float a[6]
//                double p[28]
// Return Type  : void
//
static void b_mult_poly42(const double d[15], const float a[6], double p[28])
{
  // 'mult_poly42:2' p = zeros(1, 1, 28);
  // 'mult_poly42:3' p(1) = a(:, :, 1)*d(:, :, 1);
  p[0] = a[0] * static_cast<float>(d[0]);

  // x^6
  // 'mult_poly42:4' p(2) = a(:, :, 1)*d(:, :, 2) + a(:, :, 2)*d(:, :, 1);
  p[1] = a[0] * static_cast<float>(d[1]) + a[1] * static_cast<float>(d[0]);

  // x^5*y
  // 'mult_poly42:5' p(3) = a(:, :, 1)*d(:, :, 3) + a(:, :, 2)*d(:, :, 2) + a(:, :, 3)*d(:, :, 1); 
  p[2] = (a[0] * static_cast<float>(d[2]) + a[1] * static_cast<float>(d[1])) +
    a[2] * static_cast<float>(d[0]);

  // x^4*y^2
  // 'mult_poly42:6' p(4) = a(:, :, 1)*d(:, :, 4) + a(:, :, 2)*d(:, :, 3) + a(:, :, 3)*d(:, :, 2); 
  p[3] = (a[0] * static_cast<float>(d[3]) + a[1] * static_cast<float>(d[2])) +
    a[2] * static_cast<float>(d[1]);

  // x^3*y^3
  // 'mult_poly42:7' p(5) = a(:, :, 1)*d(:, :, 5) + a(:, :, 2)*d(:, :, 4) + a(:, :, 3)*d(:, :, 3); 
  p[4] = (a[0] * static_cast<float>(d[4]) + a[1] * static_cast<float>(d[3])) +
    a[2] * static_cast<float>(d[2]);

  // x^2*y^4
  // 'mult_poly42:8' p(6) = a(:, :, 2)*d(:, :, 5) + a(:, :, 3)*d(:, :, 4);
  p[5] = a[1] * static_cast<float>(d[4]) + a[2] * static_cast<float>(d[3]);

  // x*y^5
  // 'mult_poly42:9' p(7) = a(:, :, 3)*d(:, :, 5);
  p[6] = a[2] * static_cast<float>(d[4]);

  // y^6
  // 'mult_poly42:10' p(8) = a(:, :, 4)*d(:, :, 1) + a(:, :, 1)*d(:, :, 6);
  p[7] = a[3] * static_cast<float>(d[0]) + a[0] * static_cast<float>(d[5]);

  // x^5
  // 'mult_poly42:11' p(9) = a(:, :, 4)*d(:, :, 2) + a(:, :, 5)*d(:, :, 1) + a(:, :, 1)*d(:, :, 7) + a(:, :, 2)*d(:, :, 6); 
  p[8] = ((a[3] * static_cast<float>(d[1]) + a[4] * static_cast<float>(d[0])) +
          a[0] * static_cast<float>(d[6])) + a[1] * static_cast<float>(d[5]);

  // x^4*y
  // 'mult_poly42:12' p(10) = a(:, :, 4)*d(:, :, 3) + a(:, :, 5)*d(:, :, 2) + a(:, :, 1)*d(:, :, 8) + a(:, :, 2)*d(:, :, 7) + a(:, :, 3)*d(:, :, 6); 
  p[9] = (((a[3] * static_cast<float>(d[2]) + a[4] * static_cast<float>(d[1])) +
           a[0] * static_cast<float>(d[7])) + a[1] * static_cast<float>(d[6])) +
    a[2] * static_cast<float>(d[5]);

  // x^3*y^2
  // 'mult_poly42:13' p(11) = a(:, :, 4)*d(:, :, 4) + a(:, :, 5)*d(:, :, 3) + a(:, :, 1)*d(:, :, 9) + a(:, :, 2)*d(:, :, 8) + a(:, :, 3)*d(:, :, 7); 
  p[10] = (((a[3] * static_cast<float>(d[3]) + a[4] * static_cast<float>(d[2]))
            + a[0] * static_cast<float>(d[8])) + a[1] * static_cast<float>(d[7]))
    + a[2] * static_cast<float>(d[6]);

  // x^2*y^3
  // 'mult_poly42:14' p(12) = a(:, :, 4)*d(:, :, 5) + a(:, :, 5)*d(:, :, 4) + a(:, :, 2)*d(:, :, 9) + a(:, :, 3)*d(:, :, 8); 
  p[11] = ((a[3] * static_cast<float>(d[4]) + a[4] * static_cast<float>(d[3])) +
           a[1] * static_cast<float>(d[8])) + a[2] * static_cast<float>(d[7]);

  // x*y^4
  // 'mult_poly42:15' p(13) = a(:, :, 5)*d(:, :, 5) + a(:, :, 3)*d(:, :, 9);
  p[12] = a[4] * static_cast<float>(d[4]) + a[2] * static_cast<float>(d[8]);

  // y^5
  // 'mult_poly42:16' p(14) = a(:, :, 6)*d(:, :, 1) + a(:, :, 4)*d(:, :, 6) + a(:, :, 1)*d(:, :, 10); 
  p[13] = (a[5] * static_cast<float>(d[0]) + a[3] * static_cast<float>(d[5])) +
    a[0] * static_cast<float>(d[9]);

  // x^4
  // 'mult_poly42:17' p(15) = a(:, :, 6)*d(:, :, 2) + a(:, :, 4)*d(:, :, 7) + a(:, :, 5)*d(:, :, 6) + a(:, :, 1)*d(:, :, 11) + a(:, :, 2)*d(:, :, 10); 
  p[14] = (((a[5] * static_cast<float>(d[1]) + a[3] * static_cast<float>(d[6]))
            + a[4] * static_cast<float>(d[5])) + a[0] * static_cast<float>(d[10]))
    + a[1] * static_cast<float>(d[9]);

  // x^3*y
  // 'mult_poly42:18' p(16) = a(:, :, 6)*d(:, :, 3) + a(:, :, 4)*d(:, :, 8) + a(:, :, 5)*d(:, :, 7) + a(:, :, 1)*d(:, :, 12) + a(:, :, 2)*d(:, :, 11) + a(:, :, 3)*d(:, :, 10); 
  p[15] = ((((a[5] * static_cast<float>(d[2]) + a[3] * static_cast<float>(d[7]))
             + a[4] * static_cast<float>(d[6])) + a[0] * static_cast<float>(d[11]))
           + a[1] * static_cast<float>(d[10])) + a[2] * static_cast<float>(d[9]);

  // x^2*y^2
  // 'mult_poly42:19' p(17) = a(:, :, 6)*d(:, :, 4) + a(:, :, 4)*d(:, :, 9) + a(:, :, 5)*d(:, :, 8) + a(:, :, 2)*d(:, :, 12) + a(:, :, 3)*d(:, :, 11); 
  p[16] = (((a[5] * static_cast<float>(d[3]) + a[3] * static_cast<float>(d[8]))
            + a[4] * static_cast<float>(d[7])) + a[1] * static_cast<float>(d[11]))
    + a[2] * static_cast<float>(d[10]);

  // x*y^3
  // 'mult_poly42:20' p(18) = a(:, :, 6)*d(:, :, 5) + a(:, :, 5)*d(:, :, 9) + a(:, :, 3)*d(:, :, 12); 
  p[17] = (a[5] * static_cast<float>(d[4]) + a[4] * static_cast<float>(d[8])) +
    a[2] * static_cast<float>(d[11]);

  // y^4
  // 'mult_poly42:21' p(19) = a(:, :, 6)*d(:, :, 6) + a(:, :, 1)*d(:, :, 13) + a(:, :, 4)*d(:, :, 10); 
  p[18] = (a[5] * static_cast<float>(d[5]) + a[0] * static_cast<float>(d[12])) +
    a[3] * static_cast<float>(d[9]);

  // x^3
  // 'mult_poly42:22' p(20) = a(:, :, 6)*d(:, :, 7) + a(:, :, 1)*d(:, :, 14) + a(:, :, 2)*d(:, :, 13) + a(:, :, 4)*d(:, :, 11) + a(:, :, 5)*d(:, :, 10); 
  p[19] = (((a[5] * static_cast<float>(d[6]) + a[0] * static_cast<float>(d[13]))
            + a[1] * static_cast<float>(d[12])) + a[3] * static_cast<float>(d[10]))
    + a[4] * static_cast<float>(d[9]);

  // x^2*y
  // 'mult_poly42:23' p(21) = a(:, :, 6)*d(:, :, 8) + a(:, :, 2)*d(:, :, 14) + a(:, :, 3)*d(:, :, 13) + a(:, :, 4)*d(:, :, 12) + a(:, :, 5)*d(:, :, 11); 
  p[20] = (((a[5] * static_cast<float>(d[7]) + a[1] * static_cast<float>(d[13]))
            + a[2] * static_cast<float>(d[12])) + a[3] * static_cast<float>(d[11]))
    + a[4] * static_cast<float>(d[10]);

  // x*y^2
  // 'mult_poly42:24' p(22) = a(:, :, 6)*d(:, :, 9) + a(:, :, 3)*d(:, :, 14) + a(:, :, 5)*d(:, :, 12); 
  p[21] = (a[5] * static_cast<float>(d[8]) + a[2] * static_cast<float>(d[13])) +
    a[4] * static_cast<float>(d[11]);

  // y^3
  // 'mult_poly42:25' p(23) = a(:, :, 1)*d(:, :, 15) + a(:, :, 6)*d(:, :, 10) + a(:, :, 4)*d(:, :, 13); 
  p[22] = (a[0] * static_cast<float>(d[14]) + a[5] * static_cast<float>(d[9])) +
    a[3] * static_cast<float>(d[12]);

  // x^2
  // 'mult_poly42:26' p(24) = a(:, :, 2)*d(:, :, 15) + a(:, :, 6)*d(:, :, 11) + a(:, :, 4)*d(:, :, 14) + a(:, :, 5)*d(:, :, 13); 
  p[23] = ((a[1] * static_cast<float>(d[14]) + a[5] * static_cast<float>(d[10]))
           + a[3] * static_cast<float>(d[13])) + a[4] * static_cast<float>(d[12]);

  // x*y
  // 'mult_poly42:27' p(25) = a(:, :, 3)*d(:, :, 15) + a(:, :, 6)*d(:, :, 12) + a(:, :, 5)*d(:, :, 14); 
  p[24] = (a[2] * static_cast<float>(d[14]) + a[5] * static_cast<float>(d[11]))
    + a[4] * static_cast<float>(d[13]);

  // y^2
  // 'mult_poly42:28' p(26) = a(:, :, 4)*d(:, :, 15) + a(:, :, 6)*d(:, :, 13);
  p[25] = a[3] * static_cast<float>(d[14]) + a[5] * static_cast<float>(d[12]);

  // x
  // 'mult_poly42:29' p(27) = a(:, :, 5)*d(:, :, 15) + a(:, :, 6)*d(:, :, 14);
  p[26] = a[4] * static_cast<float>(d[14]) + a[5] * static_cast<float>(d[13]);

  // y
  // 'mult_poly42:30' p(28) = a(:, :, 6)*d(:, :, 15);
  p[27] = a[5] * static_cast<float>(d[14]);

  // 1
}

//
// Arguments    : const double x[8]
// Return Type  : double
//
static double b_norm(const double x[8])
{
  double A[8];
  double s[3];
  double e[4];
  double work[2];
  int q;
  int qjj;
  int qp1;
  int qq_tmp;
  int qq;
  boolean_T apply_transform;
  int iter;
  double b;
  double nrm;
  double r;
  int kase;
  double snorm;
  int ix;
  int iy;
  int k;
  int exitg1;
  boolean_T exitg2;
  double scale;
  double sm;
  double sqds;
  std::memcpy(&A[0], &x[0], 8U * sizeof(double));
  s[0] = 0.0;
  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 0.0;
  e[3] = 0.0;
  work[0] = 0.0;
  for (q = 0; q < 2; q++) {
    qp1 = q + 2;
    qq_tmp = q + (q << 1);
    qq = qq_tmp + 1;
    apply_transform = false;
    if (q + 1 <= 1) {
      nrm = k_xnrm2(2, A, qq);
      if (nrm > 0.0) {
        apply_transform = true;
        if (A[qq - 1] < 0.0) {
          b = -nrm;
        } else {
          b = nrm;
        }

        if (std::abs(b) >= 1.0020841800044864E-292) {
          nrm = 1.0 / b;
          kase = qq + 1;
          for (k = qq; k <= kase; k++) {
            A[k - 1] *= nrm;
          }
        } else {
          kase = qq + 1;
          for (k = qq; k <= kase; k++) {
            A[k - 1] /= b;
          }
        }

        A[qq - 1]++;
        s[0] = -b;
      } else {
        s[0] = 0.0;
      }
    }

    for (iter = qp1; iter < 5; iter++) {
      qjj = q + ((iter - 1) << 1);
      if (apply_transform) {
        kase = 1 - q;
        ix = qq;
        iy = qjj;
        nrm = 0.0;
        for (k = 0; k <= kase; k++) {
          nrm += A[ix - 1] * A[iy];
          ix++;
          iy++;
        }

        nrm = -(nrm / A[qq_tmp]);
        if (nrm != 0.0) {
          ix = qq - 1;
          iy = qjj;
          kase = 1 - q;
          for (k = 0; k <= kase; k++) {
            A[iy] += nrm * A[ix];
            ix++;
            iy++;
          }
        }
      }

      e[iter - 1] = A[qjj];
    }

    nrm = l_xnrm2(3 - q, e, q + 2);
    if (nrm == 0.0) {
      e[q] = 0.0;
    } else {
      if (e[q + 1] < 0.0) {
        e[q] = -nrm;
      } else {
        e[q] = nrm;
      }

      nrm = e[q];
      if (std::abs(e[q]) >= 1.0020841800044864E-292) {
        nrm = 1.0 / e[q];
        for (k = qp1; k < 5; k++) {
          e[k - 1] *= nrm;
        }
      } else {
        for (k = qp1; k < 5; k++) {
          e[k - 1] /= nrm;
        }
      }

      e[q + 1]++;
      e[q] = -e[q];
      if (q + 2 <= 2) {
        b = 0.0;
        work[1] = 0.0;
        if ((1 - q >= 1) && (e[1] != 0.0)) {
          b = e[1] * A[3];
          work[1] = b;
        }

        if ((1 - q >= 1) && (e[2] != 0.0)) {
          b += e[2] * A[5];
          work[1] = b;
        }

        if ((1 - q >= 1) && (e[3] != 0.0)) {
          b += e[3] * A[7];
          work[1] = b;
        }

        xaxpy(1 - q, -e[1] / e[1], work, A, 4);
        xaxpy(1 - q, -e[2] / e[1], work, A, 6);
        xaxpy(1 - q, -e[3] / e[1], work, A, 8);
      }
    }
  }

  qjj = 1;
  s[1] = A[3];
  s[2] = 0.0;
  e[2] = 0.0;
  iter = 0;
  b = s[0];
  if (s[0] != 0.0) {
    nrm = std::abs(s[0]);
    r = s[0] / nrm;
    b = nrm;
    s[0] = nrm;
    e[0] /= r;
  }

  if (e[0] != 0.0) {
    nrm = std::abs(e[0]);
    r = nrm / e[0];
    e[0] = nrm;
    s[1] = A[3] * r;
  }

  nrm = std::abs(b);
  snorm = e[0];
  if (nrm > snorm) {
    snorm = nrm;
  }

  b = s[1];
  if (s[1] != 0.0) {
    nrm = std::abs(s[1]);
    r = s[1] / nrm;
    b = nrm;
    s[1] = nrm;
    e[1] /= r;
  }

  if (e[1] != 0.0) {
    e[1] = std::abs(e[1]);
    s[2] = 0.0;
  }

  nrm = std::abs(b);
  r = e[1];
  if (nrm > r) {
    r = nrm;
  }

  if (snorm <= r) {
    snorm = r;
  }

  if (snorm <= 0.0) {
    snorm = 0.0;
  }

  while ((qjj + 2 > 0) && (iter < 75)) {
    kase = qjj;
    do {
      exitg1 = 0;
      q = kase + 1;
      if (kase + 1 == 0) {
        exitg1 = 1;
      } else {
        nrm = std::abs(e[kase]);
        if ((nrm <= 2.2204460492503131E-16 * (std::abs(s[kase]) + std::abs
              (s[kase + 1]))) || (nrm <= 1.0020841800044864E-292) || ((iter > 20)
             && (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[kase] = 0.0;
          exitg1 = 1;
        } else {
          kase--;
        }
      }
    } while (exitg1 == 0);

    if (kase + 1 == qjj + 1) {
      kase = 4;
    } else {
      qq = qjj + 2;
      qq_tmp = qjj + 2;
      exitg2 = false;
      while ((!exitg2) && (qq_tmp >= kase + 1)) {
        qq = qq_tmp;
        if (qq_tmp == kase + 1) {
          exitg2 = true;
        } else {
          nrm = 0.0;
          if (qq_tmp < qjj + 2) {
            nrm = std::abs(e[qq_tmp - 1]);
          }

          if (qq_tmp > kase + 2) {
            nrm += std::abs(e[qq_tmp - 2]);
          }

          r = std::abs(s[qq_tmp - 1]);
          if ((r <= 2.2204460492503131E-16 * nrm) || (r <=
               1.0020841800044864E-292)) {
            s[qq_tmp - 1] = 0.0;
            exitg2 = true;
          } else {
            qq_tmp--;
          }
        }
      }

      if (qq == kase + 1) {
        kase = 3;
      } else if (qq == qjj + 2) {
        kase = 1;
      } else {
        kase = 2;
        q = qq;
      }
    }

    switch (kase) {
     case 1:
      r = e[qjj];
      e[qjj] = 0.0;
      kase = qjj + 1;
      for (k = kase; k >= q + 1; k--) {
        xrotg(&s[k - 1], &r, &sm, &sqds);
        if (k > q + 1) {
          r = -sqds * e[0];
          e[0] *= sm;
        }
      }
      break;

     case 2:
      r = e[q - 1];
      e[q - 1] = 0.0;
      for (k = q + 1; k <= qjj + 2; k++) {
        xrotg(&s[k - 1], &r, &sm, &sqds);
        b = e[k - 1];
        r = -sqds * b;
        e[k - 1] = b * sm;
      }
      break;

     case 3:
      kase = qjj + 1;
      nrm = s[qjj + 1];
      scale = std::abs(nrm);
      r = std::abs(s[qjj]);
      if (scale <= r) {
        scale = r;
      }

      r = std::abs(e[qjj]);
      if (scale <= r) {
        scale = r;
      }

      r = std::abs(s[q]);
      if (scale <= r) {
        scale = r;
      }

      r = std::abs(e[q]);
      if (scale <= r) {
        scale = r;
      }

      sm = nrm / scale;
      nrm = s[qjj] / scale;
      r = e[qjj] / scale;
      sqds = s[q] / scale;
      b = ((nrm + sm) * (nrm - sm) + r * r) / 2.0;
      nrm = sm * r;
      nrm *= nrm;
      if ((b != 0.0) || (nrm != 0.0)) {
        r = std::sqrt(b * b + nrm);
        if (b < 0.0) {
          r = -r;
        }

        r = nrm / (b + r);
      } else {
        r = 0.0;
      }

      r += (sqds + sm) * (sqds - sm);
      nrm = sqds * (e[q] / scale);
      for (k = q + 1; k <= kase; k++) {
        xrotg(&r, &nrm, &sm, &sqds);
        if (k > q + 1) {
          e[0] = r;
        }

        nrm = e[k - 1];
        b = s[k - 1];
        e[k - 1] = sm * nrm - sqds * b;
        r = sqds * s[k];
        s[k] *= sm;
        s[k - 1] = sm * b + sqds * nrm;
        xrotg(&s[k - 1], &r, &sm, &sqds);
        r = sm * e[k - 1] + sqds * s[k];
        s[k] = -sqds * e[k - 1] + sm * s[k];
        nrm = sqds * e[k];
        e[k] *= sm;
      }

      e[qjj] = r;
      iter++;
      break;

     default:
      if (s[q] < 0.0) {
        s[q] = -s[q];
      }

      qp1 = q + 1;
      while ((q + 1 < 3) && (s[q] < s[qp1])) {
        nrm = s[q];
        s[q] = s[qp1];
        s[qp1] = nrm;
        q = qp1;
        qp1++;
      }

      iter = 0;
      qjj--;
      break;
    }
  }

  return s[0];
}

//
// Arguments    : const float A[12]
//                float Q[9]
//                float R[12]
// Return Type  : void
//
static void b_qr(const float A[12], float Q[9], float R[12])
{
  int i;
  float tau[3];
  float b_A[12];
  float work[4];
  int b_i;
  int ix0;
  int ii;
  float atmp;
  int knt;
  float c;
  int lastv;
  float beta1;
  int lastc;
  int iaii;
  boolean_T exitg2;
  int k;
  int ia;
  int exitg1;
  int ix;
  int i1;
  for (i = 0; i < 12; i++) {
    b_A[i] = A[i];
  }

  tau[0] = 0.0F;
  tau[1] = 0.0F;
  tau[2] = 0.0F;
  work[0] = 0.0F;
  work[1] = 0.0F;
  work[2] = 0.0F;
  work[3] = 0.0F;
  for (b_i = 0; b_i < 3; b_i++) {
    ii = b_i * 3 + b_i;
    if (b_i + 1 < 3) {
      atmp = b_A[ii];
      ix0 = ii + 2;
      tau[b_i] = 0.0F;
      c = f_xnrm2(2 - b_i, b_A, ii + 2);
      if (c != 0.0F) {
        beta1 = rt_hypotf(b_A[ii], c);
        if (b_A[ii] >= 0.0F) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 9.86076132E-32F) {
          knt = -1;
          i = (ii - b_i) + 3;
          do {
            knt++;
            for (k = ix0; k <= i; k++) {
              b_A[k - 1] *= 1.01412048E+31F;
            }

            beta1 *= 1.01412048E+31F;
            atmp *= 1.01412048E+31F;
          } while (!(std::abs(beta1) >= 9.86076132E-32F));

          beta1 = rt_hypotf(atmp, f_xnrm2(2 - b_i, b_A, ii + 2));
          if (atmp >= 0.0F) {
            beta1 = -beta1;
          }

          tau[b_i] = (beta1 - atmp) / beta1;
          c = 1.0F / (atmp - beta1);
          for (k = ix0; k <= i; k++) {
            b_A[k - 1] *= c;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 9.86076132E-32F;
          }

          atmp = beta1;
        } else {
          tau[b_i] = (beta1 - b_A[ii]) / beta1;
          c = 1.0F / (b_A[ii] - beta1);
          i = (ii - b_i) + 3;
          for (k = ix0; k <= i; k++) {
            b_A[k - 1] *= c;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = atmp;
    } else {
      tau[2] = 0.0F;
    }

    atmp = b_A[ii];
    b_A[ii] = 1.0F;
    if (tau[b_i] != 0.0F) {
      lastv = 3 - b_i;
      knt = (ii - b_i) + 2;
      while ((lastv > 0) && (b_A[knt] == 0.0F)) {
        lastv--;
        knt--;
      }

      lastc = 3 - b_i;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        knt = (ii + (lastc - 1) * 3) + 3;
        ia = knt;
        do {
          exitg1 = 0;
          if (ia + 1 <= knt + lastv) {
            if (b_A[ia] != 0.0F) {
              exitg1 = 1;
            } else {
              ia++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      knt = ii + 4;
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work[0], 0, lastc * sizeof(float));
        }

        ix0 = 0;
        i = (ii + 3 * (lastc - 1)) + 4;
        for (k = knt; k <= i; k += 3) {
          ix = ii;
          c = 0.0F;
          i1 = (k + lastv) - 1;
          for (ia = k; ia <= i1; ia++) {
            c += b_A[ia - 1] * b_A[ix];
            ix++;
          }

          work[ix0] += c;
          ix0++;
        }
      }

      b_xgerc(lastv, lastc, -tau[b_i], ii + 1, work, b_A, ii + 4);
    }

    b_A[ii] = atmp;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    for (b_i = 0; b_i <= ix0; b_i++) {
      knt = b_i + 3 * ix0;
      R[knt] = b_A[knt];
    }

    i = ix0 + 2;
    if (i <= 3) {
      std::memset(&R[(ix0 * 3 + i) + -1], 0, (4 - i) * sizeof(float));
    }
  }

  R[9] = b_A[9];
  R[10] = b_A[10];
  R[11] = b_A[11];
  ii = 2;
  work[0] = 0.0F;
  work[1] = 0.0F;
  work[2] = 0.0F;
  work[3] = 0.0F;
  for (b_i = 2; b_i >= 0; b_i--) {
    iaii = (b_i + b_i * 3) + 4;
    if (b_i + 1 < 3) {
      b_A[iaii - 4] = 1.0F;
      if (tau[ii] != 0.0F) {
        lastv = 3 - b_i;
        knt = iaii - b_i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0F)) {
          lastv--;
          knt--;
        }

        lastc = 2 - b_i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = iaii + (lastc - 1) * 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= (knt + lastv) - 1) {
              if (b_A[ia - 1] != 0.0F) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(float));
          }

          ix0 = 0;
          i = iaii + 3 * (lastc - 1);
          for (k = iaii; k <= i; k += 3) {
            ix = iaii;
            c = 0.0F;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix - 4];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        b_xgerc(lastv, lastc, -tau[ii], iaii - 3, work, b_A, iaii);
      }

      ix0 = iaii - 2;
      i = (iaii - b_i) - 1;
      for (k = ix0; k <= i; k++) {
        b_A[k - 1] *= -tau[ii];
      }
    }

    b_A[iaii - 4] = 1.0F - tau[ii];
    for (ix0 = 0; ix0 < b_i; ix0++) {
      b_A[(iaii - ix0) - 5] = 0.0F;
    }

    ii--;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    Q[3 * ix0] = b_A[3 * ix0];
    knt = 3 * ix0 + 1;
    Q[knt] = b_A[knt];
    knt = 3 * ix0 + 2;
    Q[knt] = b_A[knt];
  }
}

//
// function F_row = quadruple_constraint(i, j, k, x, y, X, R)
// Arguments    : double i
//                double j
//                double k
//                const float x[4]
//                const float y[4]
//                const float X[12]
//                const float R[54]
//                float F_row[18]
// Return Type  : void
//
static void b_quadruple_constraint(double i, double j, double k, const float x[4],
  const float y[4], const float X[12], const float R[54], float F_row[18])
{
  int F_row_tmp_tmp_tmp;
  int b_i;
  int b_F_row_tmp_tmp_tmp;
  int F_row_tmp_tmp;
  float F_row_tmp[3];
  int c_F_row_tmp_tmp_tmp;
  int b_F_row_tmp_tmp;
  float b_F_row_tmp[3];
  int c_F_row_tmp_tmp;
  float a_tmp_tmp;
  float b[6];
  float a_tmp;
  float b_b[6];
  float a;

  // 'quadruple_constraint:2' F_row = zeros(1, 3, 6, 'like', X);
  // fc:
  // 'quadruple_constraint:4' F_row(1, 1, :) = (y(i) - y(k)) * mult_R(R, 1, X(:, j) - X(:, i)) ... 
  // 'quadruple_constraint:5'                         -(x(i) - x(j)) * mult_R(R, 2, X(:, k) - X(:, i)); 
  F_row_tmp_tmp_tmp = static_cast<int>(i) - 1;
  b_i = 3 * F_row_tmp_tmp_tmp;
  b_F_row_tmp_tmp_tmp = static_cast<int>(j) - 1;
  F_row_tmp_tmp = 3 * b_F_row_tmp_tmp_tmp;
  F_row_tmp[0] = X[F_row_tmp_tmp] - X[b_i];
  c_F_row_tmp_tmp_tmp = static_cast<int>(k) - 1;
  b_F_row_tmp_tmp = 3 * c_F_row_tmp_tmp_tmp;
  b_F_row_tmp[0] = X[b_F_row_tmp_tmp] - X[b_i];
  c_F_row_tmp_tmp = b_i + 1;
  F_row_tmp[1] = X[F_row_tmp_tmp + 1] - X[c_F_row_tmp_tmp];
  b_F_row_tmp[1] = X[b_F_row_tmp_tmp + 1] - X[c_F_row_tmp_tmp];
  b_i += 2;
  F_row_tmp[2] = X[F_row_tmp_tmp + 2] - X[b_i];
  b_F_row_tmp[2] = X[b_F_row_tmp_tmp + 2] - X[b_i];
  a_tmp_tmp = y[F_row_tmp_tmp_tmp] - y[c_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b[b_i] = 0.0F;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * F_row_tmp[b_i];
    }
  }

  a_tmp = x[F_row_tmp_tmp_tmp] - x[b_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0F;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] *
        b_F_row_tmp[b_i];
    }
  }

  // fs:
  // 'quadruple_constraint:7' F_row(1, 2, :) = -(y(i) - y(k)) * mult_R(R, 2, X(:, j) - X(:, i)) ... 
  // 'quadruple_constraint:8'                          -(x(i) - x(j)) * mult_R(R, 1, X(:, k) - X(:, i)); 
  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i] = a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0F;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] * F_row_tmp[b_i];
    }
  }

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0F;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * b_F_row_tmp[b_i];
    }
  }

  // 1:
  // 'quadruple_constraint:10' F_row(1, 3, :) = -(y(i) - y(k)) * x(j) * mult_R(R, 3, X(:, j) - X(:, i)) ... 
  // 'quadruple_constraint:11'                          +(x(i) - x(j)) * y(k) * mult_R(R, 3, X(:, k) - X(:, i)); 
  a = -a_tmp_tmp * x[b_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i + 1] = -a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0F;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 2] * F_row_tmp[b_i];
    }
  }

  a_tmp_tmp = a_tmp * y[c_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0F;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 2] *
        b_F_row_tmp[b_i];
    }
  }

  for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
    F_row[3 * F_row_tmp_tmp + 2] = a * b[F_row_tmp_tmp] + a_tmp_tmp *
      b_b[F_row_tmp_tmp];
  }
}

//
// Arguments    : const float A[400]
// Return Type  : float
//
static float b_rcond(const float A[400])
{
  float result;
  float normA;
  int j;
  float s;
  float b_A[400];
  int i;
  int ipiv[20];
  int jA;
  int jjA;
  int exitg1;
  float ainvnm;
  int iter;
  int kase;
  int jump;
  float x[20];
  int b_j;
  int b_i;
  float absrexk;
  result = 0.0F;
  normA = 0.0F;
  for (j = 0; j < 20; j++) {
    s = 0.0F;
    for (i = 0; i < 20; i++) {
      s += std::abs(A[i + 20 * j]);
    }

    if (s > normA) {
      normA = s;
    }
  }

  if (normA != 0.0F) {
    std::memcpy(&b_A[0], &A[0], 400U * sizeof(float));
    b_xzgetrf(b_A, ipiv, &jA);
    jjA = 19;
    do {
      exitg1 = 0;
      if (jjA + 1 > 0) {
        if (b_A[jjA + 20 * jjA] == 0.0F) {
          exitg1 = 1;
        } else {
          jjA--;
        }
      } else {
        ainvnm = 0.0F;
        iter = 2;
        kase = 1;
        jump = 1;
        j = 0;
        for (i = 0; i < 20; i++) {
          x[i] = 0.05F;
        }

        while (kase != 0) {
          if (kase == 1) {
            for (b_j = 0; b_j < 20; b_j++) {
              jjA = b_j + b_j * 20;
              b_i = 18 - b_j;
              for (i = 0; i <= b_i; i++) {
                jA = (b_j + i) + 1;
                x[jA] -= x[b_j] * b_A[(jjA + i) + 1];
              }
            }

            for (b_j = 19; b_j >= 0; b_j--) {
              jjA = b_j + b_j * 20;
              x[b_j] /= b_A[jjA];
              for (i = 0; i < b_j; i++) {
                jA = (b_j - i) - 1;
                x[jA] -= x[b_j] * b_A[(jjA - i) - 1];
              }
            }
          } else {
            for (b_j = 0; b_j < 20; b_j++) {
              jA = b_j * 20;
              s = x[b_j];
              for (i = 0; i < b_j; i++) {
                s -= b_A[jA + i] * x[i];
              }

              x[b_j] = s / b_A[jA + b_j];
            }

            for (b_j = 19; b_j >= 0; b_j--) {
              jA = b_j * 20;
              s = x[b_j];
              b_i = b_j + 2;
              for (i = 20; i >= b_i; i--) {
                s -= b_A[(jA + i) - 1] * x[i - 1];
              }

              x[b_j] = s;
            }
          }

          if (jump == 1) {
            ainvnm = 0.0F;
            for (jjA = 0; jjA < 20; jjA++) {
              s = std::abs(x[jjA]);
              ainvnm += s;
              if (s > 1.17549435E-38F) {
                x[jjA] /= s;
              } else {
                x[jjA] = 1.0F;
              }
            }

            kase = 2;
            jump = 2;
          } else if (jump == 2) {
            j = 0;
            s = std::abs(x[0]);
            for (jjA = 0; jjA < 19; jjA++) {
              absrexk = std::abs(x[jjA + 1]);
              if (absrexk > s) {
                j = jjA + 1;
                s = absrexk;
              }
            }

            iter = 2;
            std::memset(&x[0], 0, 20U * sizeof(float));
            x[j] = 1.0F;
            kase = 1;
            jump = 3;
          } else if (jump == 3) {
            ainvnm = 0.0F;
            for (jjA = 0; jjA < 20; jjA++) {
              ainvnm += std::abs(x[jjA]);
            }

            if (ainvnm <= x[0]) {
              s = 1.0F;
              for (jjA = 0; jjA < 20; jjA++) {
                x[jjA] = s * (((static_cast<float>(jjA) + 1.0F) - 1.0F) / 19.0F
                              + 1.0F);
                s = -s;
              }

              kase = 1;
              jump = 5;
            } else {
              for (jjA = 0; jjA < 20; jjA++) {
                s = std::abs(x[jjA]);
                if (s > 1.17549435E-38F) {
                  x[jjA] /= s;
                } else {
                  x[jjA] = 1.0F;
                }
              }

              kase = 2;
              jump = 4;
            }
          } else if (jump == 4) {
            jA = j;
            j = 0;
            s = std::abs(x[0]);
            for (jjA = 0; jjA < 19; jjA++) {
              absrexk = std::abs(x[jjA + 1]);
              if (absrexk > s) {
                j = jjA + 1;
                s = absrexk;
              }
            }

            if ((std::abs(x[jA]) != std::abs(x[j])) && (iter <= 5)) {
              iter++;
              std::memset(&x[0], 0, 20U * sizeof(float));
              x[j] = 1.0F;
              kase = 1;
              jump = 3;
            } else {
              s = 1.0F;
              for (jjA = 0; jjA < 20; jjA++) {
                x[jjA] = s * (((static_cast<float>(jjA) + 1.0F) - 1.0F) / 19.0F
                              + 1.0F);
                s = -s;
              }

              kase = 1;
              jump = 5;
            }
          } else {
            if (jump == 5) {
              s = 0.0F;
              for (jjA = 0; jjA < 20; jjA++) {
                s += std::abs(x[jjA]);
              }

              s = 2.0F * s / 3.0F / 20.0F;
              if (s > ainvnm) {
                ainvnm = s;
              }

              kase = 0;
            }
          }
        }

        if (ainvnm != 0.0F) {
          result = 1.0F / ainvnm / normA;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return result;
}

//
// Arguments    : const creal32_T y
// Return Type  : creal32_T
//
static creal32_T b_recip(const creal32_T y)
{
  creal32_T z;
  float brm;
  float bim;
  float d;
  brm = std::abs(y.re);
  bim = std::abs(y.im);
  if (y.im == 0.0F) {
    z.re = 1.0F / y.re;
    z.im = 0.0F;
  } else if (y.re == 0.0F) {
    z.re = 0.0F;
    z.im = -1.0F / y.im;
  } else if (brm > bim) {
    bim = y.im / y.re;
    d = y.re + bim * y.im;
    z.re = 1.0F / d;
    z.im = -bim / d;
  } else if (brm == bim) {
    bim = 0.5F;
    if (y.re < 0.0F) {
      bim = -0.5F;
    }

    d = 0.5F;
    if (y.im < 0.0F) {
      d = -0.5F;
    }

    z.re = bim / brm;
    z.im = -d / brm;
  } else {
    bim = y.re / y.im;
    d = y.im + bim * y.re;
    z.re = bim / d;
    z.im = -1.0F / d;
  }

  return z;
}

//
// Arguments    : float A[100]
//                float V[100]
// Return Type  : void
//
static void b_schur(float A[100], float V[100])
{
  int i;
  float work[10];
  int im1n_tmp;
  int ix0;
  int in;
  int alpha1_tmp;
  int k;
  float alpha1;
  int b_i;
  int knt;
  int c_i;
  float tau[9];
  float xnorm;
  float beta1;
  int jy;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int ix;
  int exitg1;
  int i1;
  for (i = 0; i < 10; i++) {
    work[i] = 0.0F;
  }

  for (i = 0; i < 9; i++) {
    im1n_tmp = i * 10 + 2;
    in = (i + 1) * 10;
    alpha1_tmp = (i + 10 * i) + 1;
    alpha1 = A[alpha1_tmp];
    if (i + 3 < 10) {
      b_i = i + 1;
    } else {
      b_i = 8;
    }

    ix0 = b_i + im1n_tmp;
    tau[i] = 0.0F;
    xnorm = d_xnrm2(8 - i, A, ix0);
    if (xnorm != 0.0F) {
      beta1 = rt_hypotf(A[alpha1_tmp], xnorm);
      if (A[alpha1_tmp] >= 0.0F) {
        beta1 = -beta1;
      }

      if (std::abs(beta1) < 9.86076132E-32F) {
        knt = -1;
        c_i = (ix0 - i) + 7;
        do {
          knt++;
          for (k = ix0; k <= c_i; k++) {
            A[k - 1] *= 1.01412048E+31F;
          }

          beta1 *= 1.01412048E+31F;
          alpha1 *= 1.01412048E+31F;
        } while (!(std::abs(beta1) >= 9.86076132E-32F));

        beta1 = rt_hypotf(alpha1, d_xnrm2(8 - i, A, ix0));
        if (alpha1 >= 0.0F) {
          beta1 = -beta1;
        }

        tau[i] = (beta1 - alpha1) / beta1;
        xnorm = 1.0F / (alpha1 - beta1);
        c_i = (ix0 - i) + 7;
        for (k = ix0; k <= c_i; k++) {
          A[k - 1] *= xnorm;
        }

        for (k = 0; k <= knt; k++) {
          beta1 *= 9.86076132E-32F;
        }

        alpha1 = beta1;
      } else {
        tau[i] = (beta1 - A[alpha1_tmp]) / beta1;
        xnorm = 1.0F / (A[alpha1_tmp] - beta1);
        c_i = (ix0 - i) + 7;
        for (k = ix0; k <= c_i; k++) {
          A[k - 1] *= xnorm;
        }

        alpha1 = beta1;
      }
    }

    A[alpha1_tmp] = 1.0F;
    jy = (i + im1n_tmp) - 1;
    ix0 = in + 1;
    if (tau[i] != 0.0F) {
      lastv = 8 - i;
      b_i = (jy - i) + 8;
      while ((lastv + 1 > 0) && (A[b_i] == 0.0F)) {
        lastv--;
        b_i--;
      }

      lastc = 10;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        knt = in + lastc;
        k = knt;
        do {
          exitg1 = 0;
          if (k <= knt + lastv * 10) {
            if (A[k - 1] != 0.0F) {
              exitg1 = 1;
            } else {
              k += 10;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = -1;
      lastc = 0;
    }

    if (lastv + 1 > 0) {
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work[0], 0, lastc * sizeof(float));
        }

        ix = jy;
        c_i = (in + 10 * lastv) + 1;
        for (knt = ix0; knt <= c_i; knt += 10) {
          b_i = 0;
          i1 = (knt + lastc) - 1;
          for (k = knt; k <= i1; k++) {
            work[b_i] += A[k - 1] * A[ix];
            b_i++;
          }

          ix++;
        }
      }

      if (-tau[i] != 0.0F) {
        knt = in;
        for (ix0 = 0; ix0 <= lastv; ix0++) {
          if (A[jy] != 0.0F) {
            xnorm = A[jy] * -tau[i];
            ix = 0;
            c_i = knt + 1;
            i1 = lastc + knt;
            for (b_i = c_i; b_i <= i1; b_i++) {
              A[b_i - 1] += work[ix] * xnorm;
              ix++;
            }
          }

          jy++;
          knt += 10;
        }
      }
    }

    b_xzlarf(9 - i, 9 - i, i + im1n_tmp, tau[i], A, (i + in) + 2, work);
    A[alpha1_tmp] = alpha1;
  }

  std::memcpy(&V[0], &A[0], 100U * sizeof(float));
  for (ix0 = 8; ix0 >= 0; ix0--) {
    k = (ix0 + 1) * 10;
    for (i = 0; i <= ix0; i++) {
      V[k + i] = 0.0F;
    }

    c_i = ix0 + 3;
    for (i = c_i; i < 11; i++) {
      knt = k + i;
      V[knt - 1] = V[knt - 11];
    }
  }

  for (i = 0; i < 10; i++) {
    V[i] = 0.0F;
  }

  V[0] = 1.0F;
  knt = 8;
  for (i = 0; i < 10; i++) {
    work[i] = 0.0F;
  }

  for (i = 8; i >= 0; i--) {
    b_i = (i + i * 10) + 11;
    if (i + 1 < 9) {
      V[b_i] = 1.0F;
      b_xzlarf(9 - i, 8 - i, b_i + 1, tau[knt], V, b_i + 11, work);
      ix0 = b_i + 2;
      c_i = (b_i - i) + 9;
      for (k = ix0; k <= c_i; k++) {
        V[k - 1] *= -tau[knt];
      }
    }

    V[b_i] = 1.0F - tau[knt];
    for (ix0 = 0; ix0 < i; ix0++) {
      V[(b_i - ix0) - 1] = 0.0F;
    }

    knt--;
  }

  b_eml_dlahqr(A, V);
  knt = 4;
  for (ix0 = 0; ix0 < 7; ix0++) {
    if (knt <= 10) {
      std::memset(&A[(ix0 * 10 + knt) + -1], 0, (11 - knt) * sizeof(float));
    }

    knt++;
  }
}

//
// function [n, xs, ys, zs] = solve_3Q3(c, e)
// c -- 3x10 coefficients matrix
// SOLVE_3Q3 Summary of this function goes here
//    Detailed explanation goes here
// Arguments    : const float c[30]
//                double *n
//                float xs_data[]
//                int xs_size[2]
//                float ys_data[]
//                int ys_size[2]
//                float zs_data[]
//                int zs_size[2]
// Return Type  : void
//
static void b_solve_3Q3(const float c[30], double *n, float xs_data[], int
  xs_size[2], float ys_data[], int ys_size[2], float zs_data[], int zs_size[2])
{
  float A[9];
  int i;
  float P[27];
  float a21;
  int b_i;
  float M[45];
  int r1;
  float ctmp[9];
  float mons[5];
  int r2;
  int r3;
  float mons_tmp;
  float P_prime[27];
  float maxval;
  float b_mons_tmp;
  float c_mons_tmp;
  float d_mons_tmp;
  float e_mons_tmp;
  float f_mons_tmp;
  int rtemp;
  float g_mons_tmp;
  float h_mons_tmp;
  int nTrailingZeros;
  float M_tmp;
  float i_mons_tmp;
  float b_M_tmp;
  float j_mons_tmp;
  float k_mons_tmp;
  float l_mons_tmp;
  float m_mons_tmp;
  float n_mons_tmp;
  float o_mons_tmp;
  float p_mons_tmp;
  float q_mons_tmp;
  float r_mons_tmp;
  float s_mons_tmp;
  float t_mons_tmp;
  float u_mons_tmp;
  float v_mons_tmp;
  float w_mons_tmp;
  float C[13];
  float y[20];
  float b_C[13];
  float c_C[13];
  creal32_T xs_complex_data[8];
  int k2;
  int companDim;
  boolean_T exitg1;
  int a_size[2];
  creal32_T a_data[64];
  float b_xs_data[8];
  boolean_T p;
  creal32_T eiga_data[8];
  int eiga_size[1];
  creal32_T beta1_data[8];
  int beta1_size[1];
  int exitg2;
  float b_A[9];

  // 'solve_3Q3:4' xs = zeros(1, 0, 'like', c);
  xs_size[0] = 1;
  xs_size[1] = 0;

  // 'solve_3Q3:5' coder.varsize('xs', [1 13], [0 1]);
  // 'solve_3Q3:6' ys = zeros(1, 0, 'like', c);
  ys_size[0] = 1;
  ys_size[1] = 0;

  // 'solve_3Q3:7' coder.varsize('ys', [1 13], [0 1]);
  // 'solve_3Q3:8' zs = zeros(1, 0, 'like', c);
  zs_size[0] = 1;
  zs_size[1] = 0;

  // 'solve_3Q3:9' coder.varsize('zs', [1 13], [0 1]);
  // 'solve_3Q3:11' A = find_A(c);
  // 'solve_3Q3:47' A = [c(1, 2), c(1, 3), c(1, 6);
  // 'solve_3Q3:48'          c(2, 2), c(2, 3), c(2, 6);
  // 'solve_3Q3:49'          c(3, 2), c(3, 3), c(3, 6)];
  A[0] = c[3];
  A[3] = c[6];
  A[6] = c[15];
  A[1] = c[4];
  A[4] = c[7];
  A[7] = c[16];
  A[2] = c[5];
  A[5] = c[8];
  A[8] = c[17];

  // 'solve_3Q3:12' if rcond(A) < eps
  if (d_rcond(A) < 2.22044605E-16F) { //todo: change back
  //if(false){
    // 'solve_3Q3:13' n = 0;
    *n = 0.0;
  } else {
    // 'solve_3Q3:16' P = find_P(c);
    // 'solve_3Q3:53' P = zeros(3, 3, 3, 'like', c);
    // [x^2, x, 1]
    // 'solve_3Q3:54' for i = 1 : 3
    for (i = 0; i < 3; i++) {
      // 'solve_3Q3:55' P(i, 1, :) = [0, -c(i, 4), -c(i, 8)];
      P[i] = 0.0F;
      P[i + 9] = -c[i + 9];
      P[i + 18] = -c[i + 21];

      // y
      // 'solve_3Q3:56' P(i, 2, :) = [0, -c(i, 5), -c(i, 9)];
      P[i + 3] = 0.0F;
      P[i + 12] = -c[i + 12];
      P[i + 21] = -c[i + 24];

      // z
      // 'solve_3Q3:57' P(i, 3, :) = [-c(i, 1), -c(i, 7), -c(i, 10)];
      P[i + 6] = -c[i];
      P[i + 15] = -c[i + 18];
      P[i + 24] = -c[i + 27];

      // 1
    }

    // 'solve_3Q3:17' P_prime = zeros(3, 3, 3, 'like', c);
    // 'solve_3Q3:18' for i = 1 : 3
    a21 = std::abs(c[4]);
    for (i = 0; i < 3; i++) {
      // 'solve_3Q3:19' P_prime(:, :, i) = A\P(:, :, i);
      for (b_i = 0; b_i < 9; b_i++) {
        ctmp[b_i] = A[b_i];
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

      ctmp[r2] = A[r2] / A[r1];
      ctmp[r3] /= ctmp[r1];
      ctmp[r2 + 3] -= ctmp[r2] * ctmp[r1 + 3];
      ctmp[r3 + 3] -= ctmp[r3] * ctmp[r1 + 3];
      ctmp[r2 + 6] -= ctmp[r2] * ctmp[r1 + 6];
      ctmp[r3 + 6] -= ctmp[r3] * ctmp[r1 + 6];
      if (std::abs(ctmp[r3 + 3]) > std::abs(ctmp[r2 + 3])) {
        rtemp = r2;
        r2 = r3;
        r3 = rtemp;
      }

      ctmp[r3 + 3] /= ctmp[r2 + 3];
      ctmp[r3 + 6] -= ctmp[r3 + 3] * ctmp[r2 + 6];
      b_i = r1 + 9 * i;
      rtemp = r2 + 9 * i;
      maxval = P[rtemp] - P[b_i] * ctmp[r2];
      nTrailingZeros = r3 + 9 * i;
      M_tmp = ctmp[r3 + 3];
      b_M_tmp = ctmp[r3 + 6];
      mons_tmp = ((P[nTrailingZeros] - P[b_i] * ctmp[r3]) - maxval * M_tmp) /
        b_M_tmp;
      P_prime[9 * i + 2] = mons_tmp;
      b_mons_tmp = ctmp[r1 + 6];
      c_mons_tmp = ctmp[r2 + 6];
      maxval -= mons_tmp * c_mons_tmp;
      d_mons_tmp = ctmp[r2 + 3];
      maxval /= d_mons_tmp;
      P_prime[9 * i + 1] = maxval;
      e_mons_tmp = ctmp[r1 + 3];
      P_prime[9 * i] = ((P[b_i] - mons_tmp * b_mons_tmp) - maxval * e_mons_tmp) /
        ctmp[r1];
      f_mons_tmp = P[b_i + 3];
      maxval = P[rtemp + 3] - f_mons_tmp * ctmp[r2];
      mons_tmp = ((P[nTrailingZeros + 3] - f_mons_tmp * ctmp[r3]) - maxval *
                  M_tmp) / b_M_tmp;
      P_prime[9 * i + 5] = mons_tmp;
      f_mons_tmp -= mons_tmp * b_mons_tmp;
      maxval -= mons_tmp * c_mons_tmp;
      maxval /= d_mons_tmp;
      P_prime[9 * i + 4] = maxval;
      f_mons_tmp -= maxval * e_mons_tmp;
      f_mons_tmp /= ctmp[r1];
      P_prime[9 * i + 3] = f_mons_tmp;
      f_mons_tmp = P[b_i + 6];
      maxval = P[rtemp + 6] - f_mons_tmp * ctmp[r2];
      mons_tmp = ((P[nTrailingZeros + 6] - f_mons_tmp * ctmp[r3]) - maxval *
                  M_tmp) / b_M_tmp;
      P_prime[9 * i + 8] = mons_tmp;
      f_mons_tmp -= mons_tmp * b_mons_tmp;
      maxval -= mons_tmp * c_mons_tmp;
      maxval /= d_mons_tmp;
      P_prime[9 * i + 7] = maxval;
      f_mons_tmp -= maxval * e_mons_tmp;
      f_mons_tmp /= ctmp[r1];
      P_prime[9 * i + 6] = f_mons_tmp;
    }

    // 'solve_3Q3:21' M = find_M(P_prime);
    // 'find_M:2' M = zeros(5, 3, 3, 'like', P);
    std::memset(&M[0], 0, 45U * sizeof(float));

    // 'find_M:3' M(:, 1, 1) = [0, 0, P(1, 2, 2)*P(2, 1, 2) - P(3, 3, 1) - P(1, 1, 2)*P(3, 1, 2) + P(3, 1, 2)*(P(1, 1, 2) - P(3, 2, 2)), P(1, 2, 2)*P(2, 1, 3) - P(3, 3, 2) + P(1, 2, 3)*P(2, 1, 2) - P(1, 1, 2)*P(3, 1, 3) - P(1, 1, 3)*P(3, 1, 2) + P(3, 1, 3)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 1, 2)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 2, 3)*P(2, 1, 3) - P(3, 3, 3) - P(1, 1, 3)*P(3, 1, 3) + P(3, 1, 3)*(P(1, 1, 3) - P(3, 2, 3))]; 
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    mons_tmp = P_prime[9] - P_prime[14];
    b_mons_tmp = P_prime[12] * P_prime[10];
    mons[2] = ((b_mons_tmp - P_prime[8]) - P_prime[9] * P_prime[11]) + P_prime
      [11] * mons_tmp;
    c_mons_tmp = P_prime[18] - P_prime[23];
    d_mons_tmp = P_prime[12] * P_prime[19];
    e_mons_tmp = P_prime[21] * P_prime[10];
    mons[3] = (((((d_mons_tmp - P_prime[17]) + e_mons_tmp) - P_prime[9] *
                 P_prime[20]) - P_prime[18] * P_prime[11]) + P_prime[20] *
               mons_tmp) + P_prime[11] * c_mons_tmp;
    f_mons_tmp = P_prime[21] * P_prime[19];
    mons[4] = ((f_mons_tmp - P_prime[26]) - P_prime[18] * P_prime[20]) +
      P_prime[20] * c_mons_tmp;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i] = mons[b_i];
    }

    //  degree = 2
    // 'find_M:5' M(:, 1, 2) = [0, 0, P(1, 3, 1) + P(1, 2, 2)*P(2, 2, 2) - P(1, 2, 2)*P(3, 1, 2) + P(3, 2, 2)*(P(1, 1, 2) - P(3, 2, 2)), P(1, 3, 2) + P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - P(1, 2, 2)*P(3, 1, 3) - P(1, 2, 3)*P(3, 1, 2) + P(3, 2, 3)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 2, 2)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 3, 3) + P(1, 2, 3)*P(2, 2, 3) - P(1, 2, 3)*P(3, 1, 3) + P(3, 2, 3)*(P(1, 1, 3) - P(3, 2, 3))]; 
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    a21 = P_prime[12] * P_prime[13];
    mons[2] = ((P_prime[6] + a21) - P_prime[12] * P_prime[11]) + P_prime[14] *
      mons_tmp;
    g_mons_tmp = P_prime[12] * P_prime[22];
    h_mons_tmp = P_prime[21] * P_prime[13];
    mons[3] = (((((P_prime[15] + g_mons_tmp) + h_mons_tmp) - P_prime[12] *
                 P_prime[20]) - P_prime[21] * P_prime[11]) + P_prime[23] *
               mons_tmp) + P_prime[14] * c_mons_tmp;
    i_mons_tmp = P_prime[21] * P_prime[22];
    mons[4] = ((P_prime[24] + i_mons_tmp) - P_prime[21] * P_prime[20]) +
      P_prime[23] * c_mons_tmp;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 15] = mons[b_i];
    }

    //  degree = 2
    // 'find_M:7' M(:, 1, 3) = [0, P(1, 2, 2)*P(2, 3, 1) - P(1, 3, 1)*P(3, 1, 2) + P(3, 3, 1)*(P(1, 1, 2) - P(3, 2, 2)), P(1, 2, 2)*P(2, 3, 2) + P(1, 2, 3)*P(2, 3, 1) - P(1, 3, 1)*P(3, 1, 3) - P(1, 3, 2)*P(3, 1, 2) + P(3, 3, 2)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 3, 1)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 2, 2)*P(2, 3, 3) + P(1, 2, 3)*P(2, 3, 2) - P(1, 3, 2)*P(3, 1, 3) - P(1, 3, 3)*P(3, 1, 2) + P(3, 3, 3)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 3, 2)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 2, 3)*P(2, 3, 3) - P(1, 3, 3)*P(3, 1, 3) + P(3, 3, 3)*(P(1, 1, 3) - P(3, 2, 3))]; 
    mons[0] = 0.0F;
    j_mons_tmp = P_prime[12] * P_prime[7];
    mons[1] = (j_mons_tmp - P_prime[6] * P_prime[11]) + P_prime[8] * mons_tmp;
    k_mons_tmp = P_prime[12] * P_prime[16];
    l_mons_tmp = P_prime[21] * P_prime[7];
    mons[2] = ((((k_mons_tmp + l_mons_tmp) - P_prime[6] * P_prime[20]) -
                P_prime[15] * P_prime[11]) + P_prime[17] * mons_tmp) + P_prime[8]
      * c_mons_tmp;
    m_mons_tmp = P_prime[12] * P_prime[25];
    n_mons_tmp = P_prime[21] * P_prime[16];
    mons[3] = ((((m_mons_tmp + n_mons_tmp) - P_prime[15] * P_prime[20]) -
                P_prime[24] * P_prime[11]) + P_prime[26] * mons_tmp) + P_prime
      [17] * c_mons_tmp;
    mons_tmp = P_prime[21] * P_prime[25];
    mons[4] = (mons_tmp - P_prime[24] * P_prime[20]) + P_prime[26] * c_mons_tmp;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 30] = mons[b_i];
    }

    //  degree = 3
    // 'find_M:9' M(:, 2, 1) = [0, 0, P(2, 1, 2)*P(3, 2, 2) - P(1, 1, 2)*P(2, 1, 2) - P(2, 3, 1) - P(3, 1, 2)*(P(2, 2, 2) - P(3, 1, 2)), P(2, 1, 2)*P(3, 2, 3) - P(1, 1, 2)*P(2, 1, 3) - P(1, 1, 3)*P(2, 1, 2) - P(2, 3, 2) + P(2, 1, 3)*P(3, 2, 2) - P(3, 1, 3)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 1, 2)*(P(2, 2, 3) - P(3, 1, 3)), P(2, 1, 3)*P(3, 2, 3) - P(1, 1, 3)*P(2, 1, 3) - P(2, 3, 3) - P(3, 1, 3)*(P(2, 2, 3) - P(3, 1, 3))]; 
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    c_mons_tmp = P_prime[13] - P_prime[11];
    o_mons_tmp = P_prime[9] * P_prime[10];
    mons[2] = ((P_prime[10] * P_prime[14] - o_mons_tmp) - P_prime[7]) - P_prime
      [11] * c_mons_tmp;
    maxval = P_prime[22] - P_prime[20];
    p_mons_tmp = P_prime[9] * P_prime[19];
    q_mons_tmp = P_prime[18] * P_prime[10];
    mons[3] = (((((P_prime[10] * P_prime[23] - p_mons_tmp) - q_mons_tmp) -
                 P_prime[16]) + P_prime[19] * P_prime[14]) - P_prime[20] *
               c_mons_tmp) - P_prime[11] * maxval;
    r_mons_tmp = P_prime[18] * P_prime[19];
    mons[4] = ((P_prime[19] * P_prime[23] - r_mons_tmp) - P_prime[25]) -
      P_prime[20] * maxval;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 5] = mons[b_i];
    }

    //  degree = 2
    // 'find_M:11' M(:, 2, 2) = [0, 0, P(3, 3, 1) - P(1, 2, 2)*P(2, 1, 2) + P(2, 2, 2)*P(3, 2, 2) - P(3, 2, 2)*(P(2, 2, 2) - P(3, 1, 2)), P(3, 3, 2) - P(1, 2, 2)*P(2, 1, 3) - P(1, 2, 3)*P(2, 1, 2) + P(2, 2, 2)*P(3, 2, 3) + P(2, 2, 3)*P(3, 2, 2) - P(3, 2, 3)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 2, 2)*(P(2, 2, 3) - P(3, 1, 3)), P(3, 3, 3) - P(1, 2, 3)*P(2, 1, 3) + P(2, 2, 3)*P(3, 2, 3) - P(3, 2, 3)*(P(2, 2, 3) - P(3, 1, 3))]; 
    mons[0] = 0.0F;
    mons[1] = 0.0F;
    mons[2] = ((P_prime[8] - b_mons_tmp) + P_prime[13] * P_prime[14]) - P_prime
      [14] * c_mons_tmp;
    mons[3] = (((((P_prime[17] - d_mons_tmp) - e_mons_tmp) + P_prime[13] *
                 P_prime[23]) + P_prime[22] * P_prime[14]) - P_prime[23] *
               c_mons_tmp) - P_prime[14] * maxval;
    mons[4] = ((P_prime[26] - f_mons_tmp) + P_prime[22] * P_prime[23]) -
      P_prime[23] * maxval;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 20] = mons[b_i];
    }

    //  degree = 2
    // 'find_M:13' M(:, 2, 3) = [0, P(2, 3, 1)*P(3, 2, 2) - P(1, 3, 1)*P(2, 1, 2) - P(3, 3, 1)*(P(2, 2, 2) - P(3, 1, 2)), P(2, 3, 1)*P(3, 2, 3) - P(1, 3, 2)*P(2, 1, 2) - P(1, 3, 1)*P(2, 1, 3) + P(2, 3, 2)*P(3, 2, 2) - P(3, 3, 2)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 3, 1)*(P(2, 2, 3) - P(3, 1, 3)), P(2, 3, 2)*P(3, 2, 3) - P(1, 3, 3)*P(2, 1, 2) - P(1, 3, 2)*P(2, 1, 3) + P(2, 3, 3)*P(3, 2, 2) - P(3, 3, 3)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 3, 2)*(P(2, 2, 3) - P(3, 1, 3)), P(2, 3, 3)*P(3, 2, 3) - P(1, 3, 3)*P(2, 1, 3) - P(3, 3, 3)*(P(2, 2, 3) - P(3, 1, 3))]; 
    mons[0] = 0.0F;
    s_mons_tmp = P_prime[6] * P_prime[10];
    mons[1] = (P_prime[7] * P_prime[14] - s_mons_tmp) - P_prime[8] * c_mons_tmp;
    t_mons_tmp = P_prime[6] * P_prime[19];
    u_mons_tmp = P_prime[15] * P_prime[10];
    mons[2] = ((((P_prime[7] * P_prime[23] - u_mons_tmp) - t_mons_tmp) +
                P_prime[16] * P_prime[14]) - P_prime[17] * c_mons_tmp) -
      P_prime[8] * maxval;
    v_mons_tmp = P_prime[15] * P_prime[19];
    w_mons_tmp = P_prime[24] * P_prime[10];
    mons[3] = ((((P_prime[16] * P_prime[23] - w_mons_tmp) - v_mons_tmp) +
                P_prime[25] * P_prime[14]) - P_prime[26] * c_mons_tmp) -
      P_prime[17] * maxval;
    c_mons_tmp = P_prime[24] * P_prime[19];
    mons[4] = (P_prime[25] * P_prime[23] - c_mons_tmp) - P_prime[26] * maxval;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 35] = mons[b_i];
    }

    //  degree = 3
    // 'find_M:15' M(:, 3, 1) = [0, 2*P(3, 1, 2)*P(3, 3, 1) - P(1, 3, 1)*P(2, 1, 2) - P(1, 1, 2)*P(2, 3, 1) - P(1, 1, 2)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 1, 2)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(3, 1, 2)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 1, 2)*P(3, 3, 2) - P(1, 1, 2)*P(2, 3, 2) - P(1, 1, 3)*P(2, 3, 1) - P(1, 3, 1)*P(2, 1, 3) - P(1, 3, 2)*P(2, 1, 2) - P(3, 1, 2)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 1, 3)*P(3, 3, 1) - P(1, 1, 3)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 1, 3)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(1, 1, 2)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 1, 2)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 1, 3)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 1, 2)*P(3, 3, 3) - P(1, 1, 2)*P(2, 3, 3) - P(1, 1, 3)*P(2, 3, 2) - P(1, 3, 2)*P(2, 1, 3) - P(1, 3, 3)*P(2, 1, 2) - P(3, 1, 3)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 1, 3)*P(3, 3, 2) - P(1, 1, 2)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 1, 2)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 1, 3)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 1, 3)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 1, 2)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)), 2*P(3, 1, 3)*P(3, 3, 3) - P(1, 3, 3)*P(2, 1, 3) - P(1, 1, 3)*P(2, 3, 3) - P(1, 1, 3)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 1, 3)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(3, 1, 3)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3))]; 
    maxval = 2.0F * P_prime[11];
    M_tmp = 2.0F * P_prime[20];
    b_M_tmp = 2.0F * P_prime[14];
    mons[0] = 0.0F;
    o_mons_tmp -= P_prime[11] * P_prime[11];
    a21 -= P_prime[14] * P_prime[14];
    mons[1] = ((((maxval * P_prime[8] - s_mons_tmp) - P_prime[9] * P_prime[7]) -
                P_prime[9] * o_mons_tmp) - P_prime[10] * a21) - P_prime[11] *
      ((P_prime[9] * P_prime[13] + b_mons_tmp) - maxval * P_prime[14]);
    b_mons_tmp = (P_prime[9] * P_prime[13] + P_prime[12] * P_prime[10]) - 2.0F *
      P_prime[11] * P_prime[14];
    mons[2] = ((((((((((maxval * P_prime[17] - P_prime[9] * P_prime[16]) -
                       P_prime[18] * P_prime[7]) - t_mons_tmp) - u_mons_tmp) -
                    P_prime[11] * (((((P_prime[9] * P_prime[22] + P_prime[18] *
      P_prime[13]) + d_mons_tmp) + e_mons_tmp) - maxval * P_prime[23]) - M_tmp *
      P_prime[14])) + M_tmp * P_prime[8]) - P_prime[18] * o_mons_tmp) - P_prime
                 [19] * a21) - P_prime[9] * ((p_mons_tmp + q_mons_tmp) - maxval *
      P_prime[20])) - P_prime[10] * ((g_mons_tmp + h_mons_tmp) - b_M_tmp *
                P_prime[23])) - P_prime[20] * b_mons_tmp;
    d_mons_tmp = r_mons_tmp - P_prime[20] * P_prime[20];
    e_mons_tmp = i_mons_tmp - P_prime[23] * P_prime[23];
    g_mons_tmp = ((((P_prime[9] * P_prime[22] + P_prime[18] * P_prime[13]) +
                    P_prime[12] * P_prime[19]) + P_prime[21] * P_prime[10]) -
                  2.0F * P_prime[11] * P_prime[23]) - 2.0F * P_prime[20] *
      P_prime[14];
    h_mons_tmp = (P_prime[9] * P_prime[19] + P_prime[18] * P_prime[10]) - 2.0F *
      P_prime[11] * P_prime[20];
    i_mons_tmp = (P_prime[12] * P_prime[22] + P_prime[21] * P_prime[13]) - 2.0F *
      P_prime[14] * P_prime[23];
    mons[3] = ((((((((((maxval * P_prime[26] - P_prime[9] * P_prime[25]) -
                       P_prime[18] * P_prime[16]) - v_mons_tmp) - w_mons_tmp) -
                    P_prime[20] * g_mons_tmp) + M_tmp * P_prime[17]) - P_prime[9]
                  * d_mons_tmp) - P_prime[10] * e_mons_tmp) - P_prime[18] *
                h_mons_tmp) - P_prime[19] * i_mons_tmp) - P_prime[11] *
      ((P_prime[18] * P_prime[22] + f_mons_tmp) - M_tmp * P_prime[23]);
    f_mons_tmp = (P_prime[18] * P_prime[22] + P_prime[21] * P_prime[19]) - 2.0F *
      P_prime[20] * P_prime[23];
    mons[4] = ((((M_tmp * P_prime[26] - c_mons_tmp) - P_prime[18] * P_prime[25])
                - P_prime[18] * d_mons_tmp) - P_prime[19] * e_mons_tmp) -
      P_prime[20] * f_mons_tmp;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 10] = mons[b_i];
    }

    //  degree = 3
    // 'find_M:17' M(:, 3, 2) = [0, 2*P(3, 2, 2)*P(3, 3, 1) - P(1, 3, 1)*P(2, 2, 2) - P(1, 2, 2)*P(2, 3, 1) - P(1, 2, 2)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 2, 2)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(3, 2, 2)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 2, 2)*P(3, 3, 2) - P(1, 2, 2)*P(2, 3, 2) - P(1, 2, 3)*P(2, 3, 1) - P(1, 3, 1)*P(2, 2, 3) - P(1, 3, 2)*P(2, 2, 2) - P(3, 2, 2)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 2, 3)*P(3, 3, 1) - P(1, 2, 3)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 2, 3)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(1, 2, 2)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 2, 2)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 2, 3)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 2, 2)*P(3, 3, 3) - P(1, 2, 2)*P(2, 3, 3) - P(1, 2, 3)*P(2, 3, 2) - P(1, 3, 2)*P(2, 2, 3) - P(1, 3, 3)*P(2, 2, 2) - P(3, 2, 3)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 2, 3)*P(3, 3, 2) - P(1, 2, 2)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 2, 2)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 2, 3)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 2, 3)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 2, 2)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)), 2*P(3, 2, 3)*P(3, 3, 3) - P(1, 3, 3)*P(2, 2, 3) - P(1, 2, 3)*P(2, 3, 3) - P(1, 2, 3)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 2, 3)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(3, 2, 3)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3))]; 
    maxval = 2.0F * P_prime[23];
    mons[0] = 0.0F;
    mons[1] = ((((b_M_tmp * P_prime[8] - P_prime[6] * P_prime[13]) - j_mons_tmp)
                - P_prime[12] * o_mons_tmp) - P_prime[13] * a21) - P_prime[14] *
      b_mons_tmp;
    mons[2] = ((((((((((b_M_tmp * P_prime[17] - k_mons_tmp) - l_mons_tmp) -
                      P_prime[6] * P_prime[22]) - P_prime[15] * P_prime[13]) -
                    P_prime[14] * g_mons_tmp) + maxval * P_prime[8]) - P_prime
                  [21] * o_mons_tmp) - P_prime[22] * a21) - P_prime[12] *
                h_mons_tmp) - P_prime[13] * i_mons_tmp) - P_prime[23] *
      b_mons_tmp;
    mons[3] = ((((((((((b_M_tmp * P_prime[26] - m_mons_tmp) - n_mons_tmp) -
                      P_prime[15] * P_prime[22]) - P_prime[24] * P_prime[13]) -
                    P_prime[23] * g_mons_tmp) + maxval * P_prime[17]) - P_prime
                  [12] * d_mons_tmp) - P_prime[13] * e_mons_tmp) - P_prime[21] *
                h_mons_tmp) - P_prime[22] * i_mons_tmp) - P_prime[14] *
      f_mons_tmp;
    mons[4] = ((((maxval * P_prime[26] - P_prime[24] * P_prime[22]) - mons_tmp)
                - P_prime[21] * d_mons_tmp) - P_prime[22] * e_mons_tmp) -
      P_prime[23] * f_mons_tmp;
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 25] = mons[b_i];
    }

    //  degree = 3
    // 'find_M:19' M(:, 3, 3) = [P(3, 3, 1)^2 - P(1, 3, 1)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 3, 1)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(3, 3, 1)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)) - P(1, 3, 1)*P(2, 3, 1), 2*P(3, 3, 1)*P(3, 3, 2) - P(1, 3, 1)*P(2, 3, 2) - P(1, 3, 2)*P(2, 3, 1) - P(3, 3, 1)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) - P(1, 3, 2)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 3, 2)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(1, 3, 1)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 3, 1)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 3, 2)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 3, 1)*P(3, 3, 3) - P(1, 3, 1)*P(2, 3, 3) - P(1, 3, 2)*P(2, 3, 2) - P(1, 3, 3)*P(2, 3, 1) - P(3, 3, 2)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) - P(1, 3, 3)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(1, 3, 1)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 3, 3)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(2, 3, 1)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 3, 2)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 3, 2)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 3, 3)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)) - P(3, 3, 1)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)) + P(3, 3, 2)^2, 2*P(3, 3, 2)*P(3, 3, 3) - P(1, 3, 2)*P(2, 3, 3) - P(1, 3, 3)*P(2, 3, 2) - P(3, 3, 3)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) - P(1, 3, 2)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 3, 2)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 3, 3)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 3, 3)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 3, 2)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)), P(3, 3, 3)^2 - P(1, 3, 3)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 3, 3)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(3, 3, 3)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)) - P(1, 3, 3)*P(2, 3, 3)]; 
    maxval = 2.0F * P_prime[8];
    mons[0] = (((P_prime[8] * P_prime[8] - P_prime[6] * o_mons_tmp) - P_prime[7]
                * a21) - P_prime[8] * b_mons_tmp) - P_prime[6] * P_prime[7];
    mons[1] = (((((((maxval * P_prime[17] - P_prime[6] * P_prime[16]) - P_prime
                    [15] * P_prime[7]) - P_prime[8] * g_mons_tmp) - P_prime[15] *
                  o_mons_tmp) - P_prime[16] * a21) - P_prime[6] * h_mons_tmp) -
               P_prime[7] * i_mons_tmp) - P_prime[17] * b_mons_tmp;
    mons[2] = ((((((((((((maxval * P_prime[26] - P_prime[6] * P_prime[25]) -
                         P_prime[15] * P_prime[16]) - P_prime[24] * P_prime[7])
                       - P_prime[17] * g_mons_tmp) - P_prime[24] * o_mons_tmp) -
                     P_prime[6] * d_mons_tmp) - P_prime[25] * a21) - P_prime[7] *
                   e_mons_tmp) - P_prime[15] * h_mons_tmp) - P_prime[16] *
                 i_mons_tmp) - P_prime[26] * b_mons_tmp) - P_prime[8] *
               f_mons_tmp) + P_prime[17] * P_prime[17];
    mons[3] = (((((((2.0F * P_prime[17] * P_prime[26] - P_prime[15] * P_prime[25])
                    - P_prime[24] * P_prime[16]) - P_prime[26] * g_mons_tmp) -
                  P_prime[15] * d_mons_tmp) - P_prime[16] * e_mons_tmp) -
                P_prime[24] * h_mons_tmp) - P_prime[25] * i_mons_tmp) - P_prime
      [17] * f_mons_tmp;
    mons[4] = (((P_prime[26] * P_prime[26] - P_prime[24] * d_mons_tmp) -
                P_prime[25] * e_mons_tmp) - P_prime[26] * f_mons_tmp) - P_prime
      [24] * P_prime[25];
    for (b_i = 0; b_i < 5; b_i++) {
      M[b_i + 40] = mons[b_i];
    }

    //  degree = 4
    // 'solve_3Q3:22' pol = find_det_M(M);
    // 'solve_3Q3:62' d_ = conv(M(:, 1, 1), find_det2(M(:, 2:3, 2:3))) - ...
    // 'solve_3Q3:63'          conv(M(:, 1, 2), find_det2(cat(3, M(:, 2:3, 1), M(:, 2:3, 3)))) + ... 
    // 'solve_3Q3:64'          conv(M(:, 1, 3), find_det2(M(:, 2:3, 1:2)));
    // 'solve_3Q3:70' d = conv(M(:, 1, 1), M(:, 2, 2)) - conv(M(:, 1, 2), M(:, 2, 1)); 
    for (b_i = 0; b_i < 9; b_i++) {
      A[b_i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 5; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        A[rtemp] += M[r1 + 40] * M[nTrailingZeros + 20];
      }
    }

    for (b_i = 0; b_i < 9; b_i++) {
      ctmp[b_i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 5; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        ctmp[rtemp] += M[r1 + 25] * M[nTrailingZeros + 35];
      }
    }

    for (b_i = 0; b_i < 9; b_i++) {
      A[b_i] -= ctmp[b_i];
    }

    for (i = 0; i < 13; i++) {
      C[i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 9; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        C[rtemp] += M[r1] * A[nTrailingZeros];
      }
    }

    rtemp = -1;
    for (r3 = 0; r3 < 10; r3++) {
      rtemp++;
      y[rtemp] = M[r3 % 5 + 5 * (r3 / 5 + 1)];
    }

    for (r3 = 0; r3 < 10; r3++) {
      rtemp++;
      y[rtemp] = M[(r3 % 5 + 5 * (r3 / 5 + 1)) + 30];
    }

    // 'solve_3Q3:70' d = conv(M(:, 1, 1), M(:, 2, 2)) - conv(M(:, 1, 2), M(:, 2, 1)); 
    for (b_i = 0; b_i < 9; b_i++) {
      A[b_i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 5; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        A[rtemp] += y[r1 + 15] * y[nTrailingZeros];
      }
    }

    for (b_i = 0; b_i < 9; b_i++) {
      ctmp[b_i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 5; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        ctmp[rtemp] += y[r1 + 5] * y[nTrailingZeros + 10];
      }
    }

    for (b_i = 0; b_i < 9; b_i++) {
      A[b_i] -= ctmp[b_i];
    }

    for (i = 0; i < 13; i++) {
      b_C[i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 9; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        b_C[rtemp] += M[r1 + 15] * A[nTrailingZeros];
      }
    }

    // 'solve_3Q3:70' d = conv(M(:, 1, 1), M(:, 2, 2)) - conv(M(:, 1, 2), M(:, 2, 1)); 
    for (b_i = 0; b_i < 9; b_i++) {
      A[b_i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 5; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        A[rtemp] += M[r1 + 25] * M[nTrailingZeros + 5];
      }
    }

    for (b_i = 0; b_i < 9; b_i++) {
      ctmp[b_i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 5; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        ctmp[rtemp] += M[r1 + 10] * M[nTrailingZeros + 20];
      }
    }

    for (b_i = 0; b_i < 9; b_i++) {
      A[b_i] -= ctmp[b_i];
    }

    for (i = 0; i < 13; i++) {
      c_C[i] = 0.0F;
    }

    for (r1 = 0; r1 < 5; r1++) {
      for (nTrailingZeros = 0; nTrailingZeros < 9; nTrailingZeros++) {
        rtemp = r1 + nTrailingZeros;
        c_C[rtemp] += M[r1 + 30] * A[nTrailingZeros];
      }
    }

    //  Due to the structure of M, d is 8-degree polynomial
    // 'solve_3Q3:66' d = d_(end-8:end);
    // 'solve_3Q3:23' assert(numel(pol)==9);
    // 'solve_3Q3:24' if ~isfinite(pol)
    // 'solve_3Q3:29' xs_complex = roots(pol');
    for (b_i = 0; b_i < 13; b_i++) {
      C[b_i] = (C[b_i] - b_C[b_i]) + c_C[b_i];
    }

    for (b_i = 0; b_i < 9; b_i++) {
      A[b_i] = C[b_i + 4];
    }

    std::memset(&xs_complex_data[0], 0, 8U * sizeof(creal32_T));
    rtemp = 1;
    while ((rtemp <= 9) && (A[rtemp - 1] == 0.0F)) {
      rtemp++;
    }

    k2 = 9;
    while ((k2 >= rtemp) && (A[k2 - 1] == 0.0F)) {
      k2--;
    }

    nTrailingZeros = 8 - k2;
    if (rtemp < k2) {
      companDim = k2 - rtemp;
      exitg1 = false;
      while ((!exitg1) && (companDim > 0)) {
        for (r3 = 1; r3 <= companDim; r3++) {
          ctmp[r3 - 1] = A[(rtemp + r3) - 1] / A[rtemp - 1];
        }

        if (r3 > companDim) {
          exitg1 = true;
        } else {
          rtemp++;
          companDim--;
        }
      }

      if (companDim < 1) {
        if (1 > 9 - k2) {
          rtemp = 0;
        } else {
          rtemp = 9 - k2;
        }
      } else {
        a_size[0] = companDim;
        a_size[1] = companDim;
        r1 = companDim * companDim;
        if (0 <= r1 - 1) {
          std::memset(&a_data[0], 0, r1 * sizeof(creal32_T));
        }

        for (r1 = 0; r1 <= companDim - 2; r1++) {
          b_i = companDim * r1;
          a_data[b_i].re = -ctmp[r1];
          a_data[b_i].im = 0.0F;
          b_i = (r1 + b_i) + 1;
          a_data[b_i].re = 1.0F;
          a_data[b_i].im = 0.0F;
        }

        b_i = companDim * (companDim - 1);
        a_data[b_i].re = -ctmp[companDim - 1];
        a_data[b_i].im = 0.0F;
        if (0 <= nTrailingZeros) {
          std::memset(&xs_complex_data[0], 0, (nTrailingZeros + 1) * sizeof
                      (creal32_T));
        }

        if (companDim == 1) {
          eiga_data[0] = a_data[0];
        } else {
          p = true;
          r3 = 0;
          exitg1 = false;
          while ((!exitg1) && (r3 <= companDim - 1)) {
            i = 0;
            do {
              exitg2 = 0;
              if (i <= r3) {
                if (a_data[i + companDim * r3].re != a_data[r3 + companDim * i].
                    re) {
                  p = false;
                  exitg2 = 1;
                } else {
                  i++;
                }
              } else {
                r3++;
                exitg2 = 2;
              }
            } while (exitg2 == 0);

            if (exitg2 == 1) {
              exitg1 = true;
            }
          }

          if (p) {
            b_xgehrd(a_data, a_size);
            b_eml_zlahqr(a_data, a_size);
            rtemp = a_size[0];
            if (3 < a_size[0]) {
              r1 = 4;
              if (a_size[0] - 4 < a_size[1] - 1) {
                r2 = a_size[0] - 3;
              } else {
                r2 = a_size[1];
              }

              for (r3 = 0; r3 < r2; r3++) {
                for (i = r1; i <= rtemp; i++) {
                  b_i = (i + a_size[0] * r3) - 1;
                  a_data[b_i].re = 0.0F;
                  a_data[b_i].im = 0.0F;
                }

                r1++;
              }
            }

            rtemp = a_size[0];
            for (r1 = 0; r1 < rtemp; r1++) {
              eiga_data[r1].re = a_data[r1 + a_size[0] * r1].re;
              eiga_data[r1].im = 0.0F;
            }
          } else {
            b_xzgeev(a_data, a_size, &rtemp, eiga_data, eiga_size, beta1_data,
                     beta1_size);
            r1 = eiga_size[0];
            for (b_i = 0; b_i < r1; b_i++) {
              if (beta1_data[b_i].im == 0.0F) {
                if (eiga_data[b_i].im == 0.0F) {
                  mons_tmp = eiga_data[b_i].re / beta1_data[b_i].re;
                  maxval = 0.0F;
                } else if (eiga_data[b_i].re == 0.0F) {
                  mons_tmp = 0.0F;
                  maxval = eiga_data[b_i].im / beta1_data[b_i].re;
                } else {
                  mons_tmp = eiga_data[b_i].re / beta1_data[b_i].re;
                  maxval = eiga_data[b_i].im / beta1_data[b_i].re;
                }
              } else if (beta1_data[b_i].re == 0.0F) {
                if (eiga_data[b_i].re == 0.0F) {
                  mons_tmp = eiga_data[b_i].im / beta1_data[b_i].im;
                  maxval = 0.0F;
                } else if (eiga_data[b_i].im == 0.0F) {
                  mons_tmp = 0.0F;
                  maxval = -(eiga_data[b_i].re / beta1_data[b_i].im);
                } else {
                  mons_tmp = eiga_data[b_i].im / beta1_data[b_i].im;
                  maxval = -(eiga_data[b_i].re / beta1_data[b_i].im);
                }
              } else {
                b_M_tmp = std::abs(beta1_data[b_i].re);
                maxval = std::abs(beta1_data[b_i].im);
                if (b_M_tmp > maxval) {
                  maxval = beta1_data[b_i].im / beta1_data[b_i].re;
                  M_tmp = beta1_data[b_i].re + maxval * beta1_data[b_i].im;
                  mons_tmp = (eiga_data[b_i].re + maxval * eiga_data[b_i].im) /
                    M_tmp;
                  maxval = (eiga_data[b_i].im - maxval * eiga_data[b_i].re) /
                    M_tmp;
                } else if (maxval == b_M_tmp) {
                  if (beta1_data[b_i].re > 0.0F) {
                    maxval = 0.5F;
                  } else {
                    maxval = -0.5F;
                  }

                  if (beta1_data[b_i].im > 0.0F) {
                    M_tmp = 0.5F;
                  } else {
                    M_tmp = -0.5F;
                  }

                  mons_tmp = (eiga_data[b_i].re * maxval + eiga_data[b_i].im *
                              M_tmp) / b_M_tmp;
                  maxval = (eiga_data[b_i].im * maxval - eiga_data[b_i].re *
                            M_tmp) / b_M_tmp;
                } else {
                  maxval = beta1_data[b_i].re / beta1_data[b_i].im;
                  M_tmp = beta1_data[b_i].im + maxval * beta1_data[b_i].re;
                  mons_tmp = (maxval * eiga_data[b_i].re + eiga_data[b_i].im) /
                    M_tmp;
                  maxval = (maxval * eiga_data[b_i].im - eiga_data[b_i].re) /
                    M_tmp;
                }
              }

              eiga_data[b_i].re = mons_tmp;
              eiga_data[b_i].im = maxval;
            }
          }
        }

        for (r1 = 0; r1 < companDim; r1++) {
          xs_complex_data[(r1 - k2) + 9] = eiga_data[r1];
        }

        rtemp = (companDim - k2) + 9;
      }
    } else if (1 > 9 - k2) {
      rtemp = 0;
    } else {
      rtemp = 9 - k2;
    }

    // 'solve_3Q3:30' xs = zeros(1, length(xs_complex), 'like', c);
    if (0 <= rtemp - 1) {
      std::memset(&xs_data[0], 0, rtemp * sizeof(float));
    }

    // 'solve_3Q3:31' n = 0;
    *n = 0.0;

    // 'solve_3Q3:33' for i = 1 : length(xs_complex)
    for (i = 0; i < rtemp; i++) {
      // 'solve_3Q3:34' n = n + 1;
      (*n)++;

      // 'solve_3Q3:35' xs(n) = real(xs_complex(i));
      xs_data[static_cast<int>(*n) - 1] = xs_complex_data[i].re;
    }

    // 'solve_3Q3:38' xs = xs(1:n);
    if (1.0 > *n) {
      r1 = 0;
    } else {
      r1 = static_cast<int>(*n);
    }

    if (0 <= r1 - 1) {
      std::memcpy(&b_xs_data[0], &xs_data[0], r1 * sizeof(float));
      std::memcpy(&C[0], &xs_data[0], r1 * sizeof(float));
    }

    xs_size[0] = 1;
    xs_size[1] = r1;
    if (0 <= r1 - 1) {
      std::memcpy(&xs_data[0], &C[0], r1 * sizeof(float));
    }

    // 'solve_3Q3:39' ys = zeros(1, n, 'like', c);
    ys_size[0] = 1;
    r1 = static_cast<int>(*n);
    ys_size[1] = r1;

    // 'solve_3Q3:40' zs = zeros(1, n, 'like', c);
    zs_size[0] = 1;
    zs_size[1] = r1;
    if (0 <= r1 - 1) {
      std::memset(&ys_data[0], 0, r1 * sizeof(float));
      std::memset(&zs_data[0], 0, r1 * sizeof(float));
    }

    // 'solve_3Q3:41' for i = 1 : n
    if (0 <= r1 - 1) {
      mons[4] = 1.0F;
    }

    for (i = 0; i < r1; i++) {
      // 'solve_3Q3:42' [ys(i), zs(i)] = find_yz(M, xs(i));
      // 'find_yz:2' M_subs = zeros(3, 'like', x);
      // 'find_yz:3' mons = [x^4, x^3, x^2, x, 1];
      mons[0] = std::pow(b_xs_data[i], 4.0F);
      mons[1] = std::pow(b_xs_data[i], 3.0F);
      mons[2] = b_xs_data[i] * b_xs_data[i];
      mons[3] = xs_data[i];

      // 'find_yz:4' for i = 1 : 3
      // 'find_yz:9' [Q, ~] = qr(M_subs');
      for (rtemp = 0; rtemp < 3; rtemp++) {
        // 'find_yz:5' for j = 1 : 3
        for (r3 = 0; r3 < 3; r3++) {
          // 'find_yz:6' M_subs(i, j) = mons * M(:, i, j);
          maxval = 0.0F;
          for (b_i = 0; b_i < 5; b_i++) {
            maxval += mons[b_i] * M[(b_i + 5 * rtemp) + 15 * r3];
          }

          b_A[r3 + 3 * rtemp] = maxval;
        }
      }

      f_qr(b_A, A, ctmp);

      // 'find_yz:9' ~
      // 'find_yz:10' y = Q(1, 3) / Q(3, 3);
      ys_data[i] = A[6] / A[8];

      // 'find_yz:11' z = Q(2, 3) / Q(3, 3);
      zs_data[i] = A[7] / A[8];
    }
  }
}

//
// Arguments    : creal_T *x
// Return Type  : void
//
static void b_sqrt(creal_T *x)
{
  double xr;
  double xi;
  double absxr;
  double absxi;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxr = 0.0;
      absxi = std::sqrt(-xr);
    } else {
      absxr = std::sqrt(xr);
      absxi = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxr = std::sqrt(-xi / 2.0);
      absxi = -absxr;
    } else {
      absxr = std::sqrt(xi / 2.0);
      absxi = absxr;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi = rt_hypotd(absxr, absxi * 0.5);
      if (absxi > absxr) {
        absxr = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0);
      } else {
        absxr = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxr = std::sqrt((rt_hypotd(absxr, absxi) + absxr) * 0.5);
    }

    if (xr > 0.0) {
      absxi = 0.5 * (xi / absxr);
    } else {
      if (xi < 0.0) {
        absxi = -absxr;
      } else {
        absxi = absxr;
      }

      absxr = 0.5 * (xi / absxi);
    }
  }

  x->re = absxr;
  x->im = absxi;
}

//
// Arguments    : int n
//                float a
//                const float x[2]
//                float y[8]
//                int iy0
// Return Type  : void
//
static void b_xaxpy(int n, float a, const float x[2], float y[8], int iy0)
{
  int iy;
  if ((n >= 1) && (a != 0.0F)) {
    iy = iy0 - 1;
    y[iy] += a * x[1];
  }
}

//
// Arguments    : float *a
//                float *b
//                float *c
//                float *d
//                float *rt1r
//                float *rt1i
//                float *rt2r
//                float *rt2i
//                float *cs
//                float *sn
// Return Type  : void
//
static void b_xdlanv2(float *a, float *b, float *c, float *d, float *rt1r, float
                      *rt1i, float *rt2r, float *rt2i, float *cs, float *sn)
{
  float tau;
  float p;
  float bcmax;
  float bcmis;
  float scale;
  int b_b;
  int b_c;
  float z;
  if (*c == 0.0F) {
    *cs = 1.0F;
    *sn = 0.0F;
  } else if (*b == 0.0F) {
    *cs = 0.0F;
    *sn = 1.0F;
    bcmax = *d;
    *d = *a;
    *a = bcmax;
    *b = -*c;
    *c = 0.0F;
  } else {
    tau = *a - *d;
    if ((tau == 0.0F) && ((*b < 0.0F) != (*c < 0.0F))) {
      *cs = 1.0F;
      *sn = 0.0F;
    } else {
      p = 0.5F * tau;
      bcmis = std::abs(*b);
      scale = std::abs(*c);
      if (bcmis > scale) {
        bcmax = bcmis;
      } else {
        bcmax = scale;
      }

      if (bcmis < scale) {
        scale = bcmis;
      }

      if (*b >= 0.0F) {
        b_b = 1;
      } else {
        b_b = -1;
      }

      if (*c >= 0.0F) {
        b_c = 1;
      } else {
        b_c = -1;
      }

      bcmis = scale * static_cast<float>(b_b) * static_cast<float>(b_c);
      scale = std::abs(p);
      if (scale <= bcmax) {
        scale = bcmax;
      }

      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817842E-16F) {
        *a = std::sqrt(scale) * std::sqrt(z);
        if (p < 0.0F) {
          *a = -*a;
        }

        z = p + *a;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        tau = rt_hypotf(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0F;
      } else {
        bcmis = *b + *c;
        tau = rt_hypotf(bcmis, tau);
        *cs = std::sqrt(0.5F * (std::abs(bcmis) / tau + 1.0F));
        if (bcmis >= 0.0F) {
          b_b = 1;
        } else {
          b_b = -1;
        }

        *sn = -(p / (tau * *cs)) * static_cast<float>(b_b);
        bcmax = *a * *cs + *b * *sn;
        scale = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        bcmis = -*c * *sn + *d * *cs;
        *b = scale * *cs + bcmis * *sn;
        *c = -bcmax * *sn + z * *cs;
        bcmax = 0.5F * ((bcmax * *cs + z * *sn) + (-scale * *sn + bcmis * *cs));
        *a = bcmax;
        *d = bcmax;
        if (*c != 0.0F) {
          if (*b != 0.0F) {
            if ((*b < 0.0F) == (*c < 0.0F)) {
              bcmis = std::sqrt(std::abs(*b));
              z = std::sqrt(std::abs(*c));
              *a = bcmis * z;
              if (*c >= 0.0F) {
                p = *a;
              } else {
                p = -*a;
              }

              tau = 1.0F / std::sqrt(std::abs(*b + *c));
              *a = bcmax + p;
              *d = bcmax - p;
              *b -= *c;
              *c = 0.0F;
              scale = bcmis * tau;
              bcmis = z * tau;
              bcmax = *cs * scale - *sn * bcmis;
              *sn = *cs * bcmis + *sn * scale;
              *cs = bcmax;
            }
          } else {
            *b = -*c;
            *c = 0.0F;
            bcmax = *cs;
            *cs = -*sn;
            *sn = bcmax;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0F) {
    *rt1i = 0.0F;
    *rt2i = 0.0F;
  } else {
    *rt1i = std::sqrt(std::abs(*b)) * std::sqrt(std::abs(*c));
    *rt2i = -*rt1i;
  }
}

//
// Arguments    : float x1
//                float x2
//                float x3
// Return Type  : float
//
static float b_xdlapy3(float x1, float x2, float x3)
{
  float y;
  float a;
  float b;
  float c;
  a = std::abs(x1);
  b = std::abs(x2);
  c = std::abs(x3);
  if (a > b) {
    y = a;
  } else {
    y = b;
  }

  if (c > y) {
    y = c;
  }

  if (y > 0.0F) {
    a /= y;
    b /= y;
    c /= y;
    y *= std::sqrt((a * a + c * c) + b * b);
  } else {
    y = (a + b) + c;
  }

  return y;
}

//
// Arguments    : creal32_T a_data[]
//                const int a_size[2]
// Return Type  : void
//
static void b_xgehrd(creal32_T a_data[], const int a_size[2])
{
  int n;
  int i;
  creal32_T work_data[8];
  int b_i;
  int c_i;
  int in;
  int alpha1_tmp;
  creal32_T alpha1;
  int n_tmp_tmp;
  int n_tmp;
  creal32_T tau_data[7];
  float c_re;
  float beta1;
  int iv0_tmp;
  int jA;
  int lastv;
  int knt;
  int lastc;
  float c_im;
  int i1;
  int rowleft;
  boolean_T exitg2;
  int ix;
  int ia;
  int exitg1;
  int i2;
  float temp_im;
  creal32_T b_alpha1;
  float a_im;
  n = a_size[0];
  i = static_cast<signed char>(a_size[0]);
  if (0 <= i - 1) {
    std::memset(&work_data[0], 0, i * sizeof(creal32_T));
  }

  b_i = a_size[0];
  for (c_i = 0; c_i <= b_i - 2; c_i++) {
    in = (c_i + 1) * n;
    alpha1_tmp = (c_i + a_size[0] * c_i) + 1;
    alpha1 = a_data[alpha1_tmp];
    i = c_i + 3;
    if (i >= n) {
      i = n;
    }

    i += c_i * n;
    n_tmp_tmp = n - c_i;
    n_tmp = n_tmp_tmp - 3;
    tau_data[c_i].re = 0.0F;
    tau_data[c_i].im = 0.0F;
    if (n_tmp + 2 > 0) {
      c_re = n_xnrm2(n_tmp + 1, a_data, i);
      if ((c_re != 0.0F) || (a_data[alpha1_tmp].im != 0.0F)) {
        beta1 = b_xdlapy3(a_data[alpha1_tmp].re, a_data[alpha1_tmp].im, c_re);
        if (a_data[alpha1_tmp].re >= 0.0F) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 9.86076132E-32F) {
          knt = -1;
          i1 = i + n_tmp;
          do {
            knt++;
            for (rowleft = i; rowleft <= i1; rowleft++) {
              c_re = 1.01412048E+31F * a_data[rowleft - 1].im;
              a_data[rowleft - 1].re *= 1.01412048E+31F;
              a_data[rowleft - 1].im = c_re;
            }

            beta1 *= 1.01412048E+31F;
            alpha1.re *= 1.01412048E+31F;
            alpha1.im *= 1.01412048E+31F;
          } while (!(std::abs(beta1) >= 9.86076132E-32F));

          beta1 = b_xdlapy3(alpha1.re, alpha1.im, n_xnrm2(n_tmp + 1, a_data, i));
          if (alpha1.re >= 0.0F) {
            beta1 = -beta1;
          }

          c_re = beta1 - alpha1.re;
          if (0.0F - alpha1.im == 0.0F) {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = 0.0F;
          } else if (c_re == 0.0F) {
            tau_data[c_i].re = 0.0F;
            tau_data[c_i].im = (0.0F - alpha1.im) / beta1;
          } else {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = (0.0F - alpha1.im) / beta1;
          }

          b_alpha1.re = alpha1.re - beta1;
          b_alpha1.im = alpha1.im;
          alpha1 = b_recip(b_alpha1);
          for (rowleft = i; rowleft <= i1; rowleft++) {
            c_re = alpha1.re * a_data[rowleft - 1].im + alpha1.im *
              a_data[rowleft - 1].re;
            a_data[rowleft - 1].re = alpha1.re * a_data[rowleft - 1].re -
              alpha1.im * a_data[rowleft - 1].im;
            a_data[rowleft - 1].im = c_re;
          }

          for (rowleft = 0; rowleft <= knt; rowleft++) {
            beta1 *= 9.86076132E-32F;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0F;
        } else {
          c_re = beta1 - a_data[alpha1_tmp].re;
          c_im = 0.0F - a_data[alpha1_tmp].im;
          if (c_im == 0.0F) {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = 0.0F;
          } else if (c_re == 0.0F) {
            tau_data[c_i].re = 0.0F;
            tau_data[c_i].im = c_im / beta1;
          } else {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = c_im / beta1;
          }

          alpha1.re = a_data[alpha1_tmp].re - beta1;
          alpha1.im = a_data[alpha1_tmp].im;
          alpha1 = b_recip(alpha1);
          i1 = i + n_tmp;
          for (rowleft = i; rowleft <= i1; rowleft++) {
            c_re = alpha1.re * a_data[rowleft - 1].im + alpha1.im *
              a_data[rowleft - 1].re;
            a_data[rowleft - 1].re = alpha1.re * a_data[rowleft - 1].re -
              alpha1.im * a_data[rowleft - 1].im;
            a_data[rowleft - 1].im = c_re;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0F;
        }
      }
    }

    a_data[alpha1_tmp].re = 1.0F;
    a_data[alpha1_tmp].im = 0.0F;
    iv0_tmp = (c_i + c_i * n) + 1;
    jA = in + 1;
    if ((tau_data[c_i].re != 0.0F) || (tau_data[c_i].im != 0.0F)) {
      lastv = n_tmp + 1;
      i = iv0_tmp + n_tmp;
      while ((lastv + 1 > 0) && ((a_data[i + 1].re == 0.0F) && (a_data[i + 1].im
               == 0.0F))) {
        lastv--;
        i--;
      }

      lastc = n;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        rowleft = in + lastc;
        ia = rowleft;
        do {
          exitg1 = 0;
          if ((n > 0) && (ia <= rowleft + lastv * n)) {
            if ((a_data[ia - 1].re != 0.0F) || (a_data[ia - 1].im != 0.0F)) {
              exitg1 = 1;
            } else {
              ia += n;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = -1;
      lastc = 0;
    }

    if (lastv + 1 > 0) {
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work_data[0], 0, lastc * sizeof(creal32_T));
        }

        ix = iv0_tmp;
        i1 = (in + n * lastv) + 1;
        for (knt = jA; n < 0 ? knt >= i1 : knt <= i1; knt += n) {
          c_re = a_data[ix].re;
          c_im = a_data[ix].im;
          i = 0;
          i2 = (knt + lastc) - 1;
          for (ia = knt; ia <= i2; ia++) {
            work_data[i].re += a_data[ia - 1].re * c_re - a_data[ia - 1].im *
              c_im;
            work_data[i].im += a_data[ia - 1].re * c_im + a_data[ia - 1].im *
              c_re;
            i++;
          }

          ix++;
        }
      }

      c_re = tau_data[c_i].re;
      c_im = tau_data[c_i].im;
      if ((-tau_data[c_i].re != 0.0F) || (-tau_data[c_i].im != 0.0F)) {
        jA = in;
        knt = iv0_tmp;
        for (i = 0; i <= lastv; i++) {
          if ((a_data[knt].re != 0.0F) || (a_data[knt].im != 0.0F)) {
            beta1 = a_data[knt].re * -c_re + a_data[knt].im * -c_im;
            temp_im = a_data[knt].re * -c_im - a_data[knt].im * -c_re;
            ix = 0;
            i1 = jA + 1;
            i2 = lastc + jA;
            for (rowleft = i1; rowleft <= i2; rowleft++) {
              a_data[rowleft - 1].re += work_data[ix].re * beta1 - work_data[ix]
                .im * temp_im;
              a_data[rowleft - 1].im += work_data[ix].re * temp_im +
                work_data[ix].im * beta1;
              ix++;
            }
          }

          knt++;
          jA += n;
        }
      }
    }

    jA = (c_i + in) + 2;
    if ((tau_data[c_i].re != 0.0F) || (-tau_data[c_i].im != 0.0F)) {
      lastv = n_tmp + 2;
      i = iv0_tmp + n_tmp;
      while ((lastv > 0) && ((a_data[i + 1].re == 0.0F) && (a_data[i + 1].im ==
               0.0F))) {
        lastv--;
        i--;
      }

      lastc = n_tmp_tmp - 2;
      exitg2 = false;
      while ((!exitg2) && (lastc + 1 > 0)) {
        i = jA + lastc * n;
        ia = i;
        do {
          exitg1 = 0;
          if (ia <= (i + lastv) - 1) {
            if ((a_data[ia - 1].re != 0.0F) || (a_data[ia - 1].im != 0.0F)) {
              exitg1 = 1;
            } else {
              ia++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = -1;
    }

    if (lastv > 0) {
      if (lastc + 1 != 0) {
        if (0 <= lastc) {
          std::memset(&work_data[0], 0, (lastc + 1) * sizeof(creal32_T));
        }

        i = 0;
        i1 = jA + n * lastc;
        for (knt = jA; n < 0 ? knt >= i1 : knt <= i1; knt += n) {
          ix = iv0_tmp;
          c_re = 0.0F;
          c_im = 0.0F;
          i2 = (knt + lastv) - 1;
          for (ia = knt; ia <= i2; ia++) {
            c_re += a_data[ia - 1].re * a_data[ix].re + a_data[ia - 1].im *
              a_data[ix].im;
            c_im += a_data[ia - 1].re * a_data[ix].im - a_data[ia - 1].im *
              a_data[ix].re;
            ix++;
          }

          work_data[i].re += c_re;
          work_data[i].im += c_im;
          i++;
        }
      }

      c_re = tau_data[c_i].re;
      c_im = tau_data[c_i].im;
      if ((-tau_data[c_i].re != 0.0F) || (tau_data[c_i].im != 0.0F)) {
        knt = 0;
        for (i = 0; i <= lastc; i++) {
          if ((work_data[knt].re != 0.0F) || (work_data[knt].im != 0.0F)) {
            beta1 = work_data[knt].re * -c_re + work_data[knt].im * c_im;
            temp_im = work_data[knt].re * c_im - work_data[knt].im * -c_re;
            ix = iv0_tmp;
            i1 = lastv + jA;
            for (rowleft = jA; rowleft < i1; rowleft++) {
              a_im = a_data[ix].re * temp_im + a_data[ix].im * beta1;
              a_data[rowleft - 1].re += a_data[ix].re * beta1 - a_data[ix].im *
                temp_im;
              a_data[rowleft - 1].im += a_im;
              ix++;
            }
          }

          knt++;
          jA += n;
        }
      }
    }

    a_data[alpha1_tmp] = alpha1;
  }
}

//
// Arguments    : int m
//                int n
//                float alpha1
//                int ix0
//                const float y[4]
//                float A[12]
//                int ia0
// Return Type  : void
//
static void b_xgerc(int m, int n, float alpha1, int ix0, const float y[4], float
                    A[12], int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i;
  int ijA;
  if (alpha1 != 0.0F) {
    jA = ia0;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

//
// Arguments    : int n
//                const double x[3]
// Return Type  : double
//
static double b_xnrm2(int n, const double x[3])
{
  double y;
  double scale;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[1]);
    } else {
      scale = 3.3121686421112381E-170;
      absxk = std::abs(x[1]);
      if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
      }

      absxk = std::abs(x[2]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : double x[100]
//                int ix0
//                int iy0
//                double c
//                double s
// Return Type  : void
//
static void b_xrot(double x[100], int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 0; k < 10; k++) {
    temp = c * x[ix] + s * x[iy];
    x[iy] = c * x[iy] - s * x[ix];
    x[ix] = temp;
    iy++;
    ix++;
  }
}

//
// Arguments    : float *a
//                float *b
//                float *c
//                float *s
// Return Type  : void
//
static void b_xrotg(float *a, float *b, float *c, float *s)
{
  float roe;
  float absa;
  float absb;
  float scale;
  float ads;
  float bds;
  roe = *b;
  absa = std::abs(*a);
  absb = std::abs(*b);
  if (absa > absb) {
    roe = *a;
  }

  scale = absa + absb;
  if (scale == 0.0F) {
    *s = 0.0F;
    *c = 1.0F;
    *a = 0.0F;
    *b = 0.0F;
  } else {
    ads = absa / scale;
    bds = absb / scale;
    scale *= std::sqrt(ads * ads + bds * bds);
    if (roe < 0.0F) {
      scale = -scale;
    }

    *c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (*c != 0.0F) {
      *b = 1.0F / *c;
    } else {
      *b = 1.0F;
    }

    *a = scale;
  }
}

//
// Arguments    : const creal32_T A_data[]
//                const int A_size[2]
//                int *info
//                creal32_T alpha1_data[]
//                int alpha1_size[1]
//                creal32_T beta1_data[]
//                int beta1_size[1]
// Return Type  : void
//
static void b_xzgeev(const creal32_T A_data[], const int A_size[2], int *info,
                     creal32_T alpha1_data[], int alpha1_size[1], creal32_T
                     beta1_data[], int beta1_size[1])
{
  int At_size[2];
  int ii;
  creal32_T At_data[64];
  float anrm;
  int jcolp1;
  boolean_T ilascl;
  float absxk;
  float anrmto;
  boolean_T guard1 = false;
  int ilo;
  float ctoc;
  int ihi;
  boolean_T notdone;
  int exitg2;
  float stemp_im;
  int i;
  int n;
  float cto1;
  int j;
  int jcol;
  float a;
  boolean_T exitg3;
  int b_i;
  int jrow;
  int nzcount;
  creal32_T atmp;
  boolean_T exitg4;
  int exitg1;
  At_size[0] = A_size[0];
  At_size[1] = A_size[1];
  ii = A_size[0] * A_size[1];
  if (0 <= ii - 1) {
    std::memcpy(&At_data[0], &A_data[0], ii * sizeof(creal32_T));
  }

  anrm = 0.0F;
  ii = A_size[0] * A_size[1];
  for (jcolp1 = 0; jcolp1 < ii; jcolp1++) {
    absxk = rt_hypotf(At_data[jcolp1].re, At_data[jcolp1].im);
    if (absxk > anrm) {
      anrm = absxk;
    }
  }

  ilascl = false;
  anrmto = anrm;
  guard1 = false;
  if ((anrm > 0.0F) && (anrm < 9.09494702E-13F)) {
    anrmto = 9.09494702E-13F;
    ilascl = true;
    guard1 = true;
  } else {
    if (anrm > 1.09951163E+12F) {
      anrmto = 1.09951163E+12F;
      ilascl = true;
      guard1 = true;
    }
  }

  if (guard1) {
    absxk = anrm;
    ctoc = anrmto;
    notdone = true;
    while (notdone) {
      stemp_im = absxk * 1.97215226E-31F;
      cto1 = ctoc / 5.0706024E+30F;
      if ((stemp_im > ctoc) && (ctoc != 0.0F)) {
        a = 1.97215226E-31F;
        absxk = stemp_im;
      } else if (cto1 > absxk) {
        a = 5.0706024E+30F;
        ctoc = cto1;
      } else {
        a = ctoc / absxk;
        notdone = false;
      }

      for (b_i = 0; b_i < ii; b_i++) {
        At_data[b_i].re *= a;
        At_data[b_i].im *= a;
      }
    }
  }

  ilo = 1;
  ihi = A_size[0];
  if (A_size[0] <= 1) {
    ihi = 1;
  } else {
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      notdone = false;
      ii = ihi;
      exitg3 = false;
      while ((!exitg3) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        jcol = 0;
        exitg4 = false;
        while ((!exitg4) && (jcol <= ihi - 1)) {
          b_i = (ii + At_size[0] * jcol) - 1;
          if ((At_data[b_i].re != 0.0F) || (At_data[b_i].im != 0.0F) || (ii ==
               jcol + 1)) {
            if (nzcount == 0) {
              j = jcol + 1;
              nzcount = 1;
              jcol++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            jcol++;
          }
        }

        if (nzcount < 2) {
          notdone = true;
          exitg3 = true;
        } else {
          ii--;
        }
      }

      if (!notdone) {
        exitg2 = 2;
      } else {
        n = At_size[0];
        if (i != ihi) {
          for (jcolp1 = 1; jcolp1 <= n; jcolp1++) {
            ii = At_size[0] * (jcolp1 - 1);
            nzcount = (i + ii) - 1;
            atmp = At_data[nzcount];
            jcol = (ihi + ii) - 1;
            At_data[nzcount] = At_data[jcol];
            At_data[jcol] = atmp;
          }
        }

        if (j != ihi) {
          for (jcolp1 = 0; jcolp1 < ihi; jcolp1++) {
            ii = jcolp1 + At_size[0] * (j - 1);
            atmp = At_data[ii];
            jcol = jcolp1 + At_size[0] * (ihi - 1);
            At_data[ii] = At_data[jcol];
            At_data[jcol] = atmp;
          }
        }

        ihi--;
        if (ihi == 1) {
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0);

    if (exitg2 != 1) {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        notdone = false;
        jcol = ilo;
        exitg3 = false;
        while ((!exitg3) && (jcol <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = jcol;
          ii = ilo;
          exitg4 = false;
          while ((!exitg4) && (ii <= ihi)) {
            b_i = (ii + At_size[0] * (jcol - 1)) - 1;
            if ((At_data[b_i].re != 0.0F) || (At_data[b_i].im != 0.0F) || (ii ==
                 jcol)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                ii++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              ii++;
            }
          }

          if (nzcount < 2) {
            notdone = true;
            exitg3 = true;
          } else {
            jcol++;
          }
        }

        if (!notdone) {
          exitg1 = 1;
        } else {
          n = At_size[0];
          if (i != ilo) {
            for (jcolp1 = ilo; jcolp1 <= n; jcolp1++) {
              ii = At_size[0] * (jcolp1 - 1);
              nzcount = (i + ii) - 1;
              atmp = At_data[nzcount];
              jcol = (ilo + ii) - 1;
              At_data[nzcount] = At_data[jcol];
              At_data[jcol] = atmp;
            }
          }

          if (j != ilo) {
            for (jcolp1 = 0; jcolp1 < ihi; jcolp1++) {
              ii = jcolp1 + At_size[0] * (j - 1);
              atmp = At_data[ii];
              jcol = jcolp1 + At_size[0] * (ilo - 1);
              At_data[ii] = At_data[jcol];
              At_data[jcol] = atmp;
            }
          }

          ilo++;
          if (ilo == ihi) {
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }
  }

  n = A_size[0];
  if ((A_size[0] > 1) && (ihi >= ilo + 2)) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      jcolp1 = jcol + 2;
      for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
        c_xzlartg(At_data[(jrow + At_size[0] * jcol) - 1], At_data[jrow +
                  At_size[0] * jcol], &absxk, &atmp, &At_data[(jrow + At_size[0]
                   * jcol) - 1]);
        b_i = jrow + At_size[0] * jcol;
        At_data[b_i].re = 0.0F;
        At_data[b_i].im = 0.0F;
        for (j = jcolp1; j <= n; j++) {
          ii = jrow + At_size[0] * (j - 1);
          nzcount = ii - 1;
          ctoc = absxk * At_data[nzcount].re + (atmp.re * At_data[ii].re -
            atmp.im * At_data[ii].im);
          stemp_im = absxk * At_data[nzcount].im + (atmp.re * At_data[ii].im +
            atmp.im * At_data[ii].re);
          cto1 = At_data[nzcount].re;
          At_data[ii].re = absxk * At_data[ii].re - (atmp.re * At_data[nzcount].
            re + atmp.im * At_data[nzcount].im);
          At_data[ii].im = absxk * At_data[ii].im - (atmp.re * At_data[nzcount].
            im - atmp.im * cto1);
          At_data[nzcount].re = ctoc;
          At_data[nzcount].im = stemp_im;
        }

        atmp.re = -atmp.re;
        atmp.im = -atmp.im;
        for (i = 1; i <= ihi; i++) {
          ii = (i + At_size[0] * (jrow - 1)) - 1;
          nzcount = (i + At_size[0] * jrow) - 1;
          ctoc = absxk * At_data[nzcount].re + (atmp.re * At_data[ii].re -
            atmp.im * At_data[ii].im);
          stemp_im = absxk * At_data[nzcount].im + (atmp.re * At_data[ii].im +
            atmp.im * At_data[ii].re);
          cto1 = At_data[nzcount].re;
          At_data[ii].re = absxk * At_data[ii].re - (atmp.re * At_data[nzcount].
            re + atmp.im * At_data[nzcount].im);
          At_data[ii].im = absxk * At_data[ii].im - (atmp.re * At_data[nzcount].
            im - atmp.im * cto1);
          At_data[nzcount].re = ctoc;
          At_data[nzcount].im = stemp_im;
        }
      }
    }
  }

  b_xzhgeqz(At_data, At_size, ilo, ihi, info, alpha1_data, alpha1_size,
            beta1_data, beta1_size);
  if ((*info == 0) && ilascl) {
    notdone = true;
    while (notdone) {
      stemp_im = anrmto * 1.97215226E-31F;
      cto1 = anrm / 5.0706024E+30F;
      if ((stemp_im > anrm) && (anrm != 0.0F)) {
        a = 1.97215226E-31F;
        anrmto = stemp_im;
      } else if (cto1 > anrmto) {
        a = 5.0706024E+30F;
        anrm = cto1;
      } else {
        a = anrm / anrmto;
        notdone = false;
      }

      ii = alpha1_size[0];
      for (b_i = 0; b_i < ii; b_i++) {
        alpha1_data[b_i].re *= a;
        alpha1_data[b_i].im *= a;
      }
    }
  }
}

//
// Arguments    : float A[400]
//                int ipiv[20]
//                int *info
// Return Type  : void
//
static void b_xzgetrf(float A[400], int ipiv[20], int *info)
{
  int i;
  int j;
  int mmj_tmp;
  int b;
  int jj;
  int jp1j;
  int iy;
  int jA;
  int ix;
  float smax;
  int k;
  float s;
  for (i = 0; i < 20; i++) {
    ipiv[i] = i + 1;
  }

  *info = 0;
  for (j = 0; j < 19; j++) {
    mmj_tmp = 18 - j;
    b = j * 21;
    jj = j * 21;
    jp1j = b + 2;
    iy = 20 - j;
    jA = 0;
    ix = b;
    smax = std::abs(A[jj]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (A[jj + jA] != 0.0F) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = iy + 1;
        ix = j;
        for (k = 0; k < 20; k++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 20;
          iy += 20;
        }
      }

      i = (jj - j) + 20;
      for (iy = jp1j; iy <= i; iy++) {
        A[iy - 1] /= A[jj];
      }
    } else {
      *info = j + 1;
    }

    iy = b + 20;
    jA = jj;
    for (k = 0; k <= mmj_tmp; k++) {
      smax = A[iy];
      if (A[iy] != 0.0F) {
        ix = jj + 1;
        i = jA + 22;
        b = (jA - j) + 40;
        for (jp1j = i; jp1j <= b; jp1j++) {
          A[jp1j - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 20;
      jA += 20;
    }
  }

  if ((*info == 0) && (A[399] == 0.0F)) {
    *info = 20;
  }
}

//
// Arguments    : creal32_T A[100]
//                int *info
//                creal32_T alpha1[10]
//                creal32_T beta1[10]
//                creal32_T V[100]
// Return Type  : void
//
static void b_xzggev(creal32_T A[100], int *info, creal32_T alpha1[10],
                     creal32_T beta1[10], creal32_T V[100])
{
  float anrm;
  int k;
  boolean_T ilascl;
  float absxk;
  float anrmto;
  boolean_T guard1 = false;
  int i;
  float ctoc;
  boolean_T notdone;
  int ilo;
  int rscale[10];
  int ihi;
  float stemp_im;
  int exitg2;
  float cto1;
  int j;
  float a;
  int ii;
  boolean_T exitg3;
  int jcolp1;
  int nzcount;
  int jcol;
  boolean_T exitg4;
  creal32_T atmp;
  signed char b_I[100];
  int exitg1;
  int jrow;
  float f;
  anrm = 0.0F;
  for (k = 0; k < 100; k++) {
    absxk = rt_hypotf(A[k].re, A[k].im);
    if (absxk > anrm) {
      anrm = absxk;
    }
  }

  ilascl = false;
  anrmto = anrm;
  guard1 = false;
  if ((anrm > 0.0F) && (anrm < 9.09494702E-13F)) {
    anrmto = 9.09494702E-13F;
    ilascl = true;
    guard1 = true;
  } else {
    if (anrm > 1.09951163E+12F) {
      anrmto = 1.09951163E+12F;
      ilascl = true;
      guard1 = true;
    }
  }

  if (guard1) {
    absxk = anrm;
    ctoc = anrmto;
    notdone = true;
    while (notdone) {
      stemp_im = absxk * 1.97215226E-31F;
      cto1 = ctoc / 5.0706024E+30F;
      if ((stemp_im > ctoc) && (ctoc != 0.0F)) {
        a = 1.97215226E-31F;
        absxk = stemp_im;
      } else if (cto1 > absxk) {
        a = 5.0706024E+30F;
        ctoc = cto1;
      } else {
        a = ctoc / absxk;
        notdone = false;
      }

      for (jcolp1 = 0; jcolp1 < 100; jcolp1++) {
        A[jcolp1].re *= a;
        A[jcolp1].im *= a;
      }
    }
  }

  for (i = 0; i < 10; i++) {
    rscale[i] = 1;
  }

  ilo = 1;
  ihi = 10;
  do {
    exitg2 = 0;
    i = 0;
    j = 0;
    notdone = false;
    ii = ihi;
    exitg3 = false;
    while ((!exitg3) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = ihi;
      k = 0;
      exitg4 = false;
      while ((!exitg4) && (k <= ihi - 1)) {
        jcol = (ii + 10 * k) - 1;
        if ((A[jcol].re != 0.0F) || (A[jcol].im != 0.0F) || (ii == k + 1)) {
          if (nzcount == 0) {
            j = k + 1;
            nzcount = 1;
            k++;
          } else {
            nzcount = 2;
            exitg4 = true;
          }
        } else {
          k++;
        }
      }

      if (nzcount < 2) {
        notdone = true;
        exitg3 = true;
      } else {
        ii--;
      }
    }

    if (!notdone) {
      exitg2 = 2;
    } else {
      if (i != ihi) {
        for (k = 0; k < 10; k++) {
          jcol = (i + 10 * k) - 1;
          atmp = A[jcol];
          jcolp1 = (ihi + 10 * k) - 1;
          A[jcol] = A[jcolp1];
          A[jcolp1] = atmp;
        }
      }

      if (j != ihi) {
        for (k = 0; k < ihi; k++) {
          jcol = k + 10 * (j - 1);
          atmp = A[jcol];
          jcolp1 = k + 10 * (ihi - 1);
          A[jcol] = A[jcolp1];
          A[jcolp1] = atmp;
        }
      }

      rscale[ihi - 1] = j;
      ihi--;
      if (ihi == 1) {
        rscale[0] = 1;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 != 1) {
    do {
      exitg1 = 0;
      i = 0;
      j = 0;
      notdone = false;
      k = ilo;
      exitg3 = false;
      while ((!exitg3) && (k <= ihi)) {
        nzcount = 0;
        i = ihi;
        j = k;
        ii = ilo;
        exitg4 = false;
        while ((!exitg4) && (ii <= ihi)) {
          jcol = (ii + 10 * (k - 1)) - 1;
          if ((A[jcol].re != 0.0F) || (A[jcol].im != 0.0F) || (ii == k)) {
            if (nzcount == 0) {
              i = ii;
              nzcount = 1;
              ii++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            ii++;
          }
        }

        if (nzcount < 2) {
          notdone = true;
          exitg3 = true;
        } else {
          k++;
        }
      }

      if (!notdone) {
        exitg1 = 1;
      } else {
        if (i != ilo) {
          for (k = ilo; k < 11; k++) {
            jcol = 10 * (k - 1);
            ii = (i + jcol) - 1;
            atmp = A[ii];
            jcolp1 = (ilo + jcol) - 1;
            A[ii] = A[jcolp1];
            A[jcolp1] = atmp;
          }
        }

        if (j != ilo) {
          for (k = 0; k < ihi; k++) {
            jcol = k + 10 * (j - 1);
            atmp = A[jcol];
            jcolp1 = k + 10 * (ilo - 1);
            A[jcol] = A[jcolp1];
            A[jcolp1] = atmp;
          }
        }

        rscale[ilo - 1] = j;
        ilo++;
        if (ilo == ihi) {
          rscale[ilo - 1] = ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }

  std::memset(&b_I[0], 0, 100U * sizeof(signed char));
  for (k = 0; k < 10; k++) {
    b_I[k + 10 * k] = 1;
  }

  for (jcolp1 = 0; jcolp1 < 100; jcolp1++) {
    V[jcolp1].re = b_I[jcolp1];
    V[jcolp1].im = 0.0F;
  }

  if (ihi >= ilo + 2) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      jcolp1 = jcol + 2;
      for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
        k = jrow + 10 * jcol;
        c_xzlartg(A[k - 1], A[k], &absxk, &atmp, &A[(jrow + 10 * jcol) - 1]);
        A[k].re = 0.0F;
        A[k].im = 0.0F;
        for (j = jcolp1; j < 11; j++) {
          ii = jrow + 10 * (j - 1);
          nzcount = ii - 1;
          ctoc = absxk * A[nzcount].re + (atmp.re * A[ii].re - atmp.im * A[ii].
            im);
          stemp_im = absxk * A[nzcount].im + (atmp.re * A[ii].im + atmp.im *
            A[ii].re);
          cto1 = A[nzcount].im;
          a = A[nzcount].re;
          A[ii].re = absxk * A[ii].re - (atmp.re * A[nzcount].re + atmp.im *
            A[nzcount].im);
          A[ii].im = absxk * A[ii].im - (atmp.re * cto1 - atmp.im * a);
          A[nzcount].re = ctoc;
          A[nzcount].im = stemp_im;
        }

        atmp.re = -atmp.re;
        atmp.im = -atmp.im;
        for (i = 1; i <= ihi; i++) {
          ii = (i + 10 * (jrow - 1)) - 1;
          nzcount = (i + 10 * jrow) - 1;
          ctoc = absxk * A[nzcount].re + (atmp.re * A[ii].re - atmp.im * A[ii].
            im);
          stemp_im = absxk * A[nzcount].im + (atmp.re * A[ii].im + atmp.im *
            A[ii].re);
          cto1 = A[nzcount].im;
          a = A[nzcount].re;
          A[ii].re = absxk * A[ii].re - (atmp.re * A[nzcount].re + atmp.im *
            A[nzcount].im);
          A[ii].im = absxk * A[ii].im - (atmp.re * cto1 - atmp.im * a);
          A[nzcount].re = ctoc;
          A[nzcount].im = stemp_im;
        }

        cto1 = atmp.re;
        a = atmp.im;
        for (i = 0; i < 10; i++) {
          ii = i + 10 * (jrow - 1);
          nzcount = i + 10 * jrow;
          ctoc = absxk * V[nzcount].re + (cto1 * V[ii].re - a * V[ii].im);
          stemp_im = absxk * V[nzcount].im + (cto1 * V[ii].im + a * V[ii].re);
          f = V[nzcount].re;
          V[ii].re = absxk * V[ii].re - (cto1 * V[nzcount].re + a * V[nzcount].
            im);
          V[ii].im = absxk * V[ii].im - (cto1 * V[nzcount].im - a * f);
          V[nzcount].re = ctoc;
          V[nzcount].im = stemp_im;
        }
      }
    }
  }

  d_xzhgeqz(A, ilo, ihi, V, info, alpha1, beta1);
  if (*info == 0) {
    b_xztgevc(A, V);
    if (ilo > 1) {
      for (i = ilo - 2; i + 1 >= 1; i--) {
        k = rscale[i] - 1;
        if (rscale[i] != i + 1) {
          for (j = 0; j < 10; j++) {
            jcol = i + 10 * j;
            atmp = V[jcol];
            nzcount = k + 10 * j;
            V[jcol] = V[nzcount];
            V[nzcount] = atmp;
          }
        }
      }
    }

    if (ihi < 10) {
      jcolp1 = ihi + 1;
      for (i = jcolp1; i < 11; i++) {
        ii = rscale[i - 1];
        if (ii != i) {
          for (j = 0; j < 10; j++) {
            jcol = (i + 10 * j) - 1;
            atmp = V[jcol];
            nzcount = (ii + 10 * j) - 1;
            V[jcol] = V[nzcount];
            V[nzcount] = atmp;
          }
        }
      }
    }

    for (ii = 0; ii < 10; ii++) {
      absxk = std::abs(V[10 * ii].re) + std::abs(V[10 * ii].im);
      for (nzcount = 0; nzcount < 9; nzcount++) {
        k = (nzcount + 10 * ii) + 1;
        ctoc = std::abs(V[k].re) + std::abs(V[k].im);
        if (ctoc > absxk) {
          absxk = ctoc;
        }
      }

      if (absxk >= 9.09494702E-13F) {
        absxk = 1.0F / absxk;
        for (nzcount = 0; nzcount < 10; nzcount++) {
          jcolp1 = nzcount + 10 * ii;
          V[jcolp1].re *= absxk;
          V[jcolp1].im *= absxk;
        }
      }
    }

    if (ilascl) {
      notdone = true;
      while (notdone) {
        stemp_im = anrmto * 1.97215226E-31F;
        cto1 = anrm / 5.0706024E+30F;
        if ((stemp_im > anrm) && (anrm != 0.0F)) {
          a = 1.97215226E-31F;
          anrmto = stemp_im;
        } else if (cto1 > anrmto) {
          a = 5.0706024E+30F;
          anrm = cto1;
        } else {
          a = anrm / anrmto;
          notdone = false;
        }

        for (jcolp1 = 0; jcolp1 < 10; jcolp1++) {
          alpha1[jcolp1].re *= a;
          alpha1[jcolp1].im *= a;
        }
      }
    }
  }
}

//
// Arguments    : const creal32_T A_data[]
//                const int A_size[2]
//                int ilo
//                int ihi
//                int *info
//                creal32_T alpha1_data[]
//                int alpha1_size[1]
//                creal32_T beta1_data[]
//                int beta1_size[1]
// Return Type  : void
//
static void b_xzhgeqz(const creal32_T A_data[], const int A_size[2], int ilo,
                      int ihi, int *info, creal32_T alpha1_data[], int
                      alpha1_size[1], creal32_T beta1_data[], int beta1_size[1])
{
  int A_size_idx_0;
  int jm1;
  creal32_T b_A_data[64];
  int n;
  int ctemp_tmp;
  float eshift_re;
  float eshift_im;
  creal32_T ctemp;
  float anorm;
  float scale;
  float b_atol;
  boolean_T firstNonZero;
  int j;
  float ascale;
  int i;
  float bscale;
  float temp1;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  float temp2;
  int ifirst;
  int istart;
  int ilast;
  int ilastm1;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T b_guard1 = false;
  boolean_T guard3 = false;
  boolean_T exitg2;
  float ascale_re;
  creal32_T shift;
  float ascale_im;
  float ad22_re;
  float ad22_im;
  float t1_re;
  float t1_im;
  float t1_im_tmp;
  creal32_T b_ascale;
  A_size_idx_0 = A_size[0];
  jm1 = A_size[0] * A_size[1];
  if (0 <= jm1 - 1) {
    std::memcpy(&b_A_data[0], &A_data[0], jm1 * sizeof(creal32_T));
  }

  *info = 0;
  if ((A_size[0] == 1) && (A_size[1] == 1)) {
    ihi = 1;
  }

  n = A_size[0];
  alpha1_size[0] = A_size[0];
  jm1 = A_size[0];
  if (0 <= jm1 - 1) {
    std::memset(&alpha1_data[0], 0, jm1 * sizeof(creal32_T));
  }

  beta1_size[0] = A_size[0];
  jm1 = A_size[0];
  for (ctemp_tmp = 0; ctemp_tmp < jm1; ctemp_tmp++) {
    beta1_data[ctemp_tmp].re = 1.0F;
    beta1_data[ctemp_tmp].im = 0.0F;
  }

  eshift_re = 0.0F;
  eshift_im = 0.0F;
  ctemp.re = 0.0F;
  ctemp.im = 0.0F;
  anorm = 0.0F;
  if (ilo <= ihi) {
    scale = 0.0F;
    anorm = 0.0F;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ctemp_tmp = j + 1;
      if (ihi < j + 1) {
        ctemp_tmp = ihi;
      }

      for (i = ilo; i <= ctemp_tmp; i++) {
        jm1 = (i + A_size[0] * (j - 1)) - 1;
        if (A_data[jm1].re != 0.0F) {
          temp1 = std::abs(A_data[jm1].re);
          if (firstNonZero) {
            anorm = 1.0F;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0F;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }

        if (A_data[jm1].im != 0.0F) {
          temp1 = std::abs(A_data[jm1].im);
          if (firstNonZero) {
            anorm = 1.0F;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0F;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * std::sqrt(anorm);
  }

  scale = 1.1920929E-7F * anorm;
  b_atol = 1.17549435E-38F;
  if (scale > 1.17549435E-38F) {
    b_atol = scale;
  }

  scale = 1.17549435E-38F;
  if (anorm > 1.17549435E-38F) {
    scale = anorm;
  }

  ascale = 1.0F / scale;
  bscale = 1.0F / std::sqrt(static_cast<float>(A_size[0]));
  firstNonZero = true;
  ctemp_tmp = ihi + 1;
  for (j = ctemp_tmp; j <= n; j++) {
    alpha1_data[j - 1] = A_data[(j + A_size[0] * (j - 1)) - 1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          ctemp_tmp = ilast + A_size_idx_0 * ilastm1;
          if (std::abs(b_A_data[ctemp_tmp].re) + std::abs(b_A_data[ctemp_tmp].im)
              <= b_atol) {
            b_A_data[ctemp_tmp].re = 0.0F;
            b_A_data[ctemp_tmp].im = 0.0F;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                ctemp_tmp = j + A_size_idx_0 * (j - 1);
                if (std::abs(b_A_data[ctemp_tmp].re) + std::abs
                    (b_A_data[ctemp_tmp].im) <= b_atol) {
                  b_A_data[ctemp_tmp].re = 0.0F;
                  b_A_data[ctemp_tmp].im = 0.0F;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }

            if (goto70) {
              b_guard1 = true;
            } else {
              jm1 = alpha1_size[0];
              if (0 <= jm1 - 1) {
                std::memset(&alpha1_data[0], 0, jm1 * sizeof(creal32_T));
              }

              jm1 = beta1_size[0];
              if (0 <= jm1 - 1) {
                std::memset(&beta1_data[0], 0, jm1 * sizeof(creal32_T));
              }

              *info = 1;
              exitg1 = 1;
            }
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1_data[ilast] = b_A_data[ilast + A_size_idx_0 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0F;
              eshift_im = 0.0F;
              ilastm = ilast + 1;
              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                jm1 = ilastm1 + A_size_idx_0 * ilastm1;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0F) {
                  shift.re = anorm / bscale;
                  shift.im = 0.0F;
                } else if (anorm == 0.0F) {
                  shift.re = 0.0F;
                  shift.im = scale / bscale;
                } else {
                  shift.re = anorm / bscale;
                  shift.im = scale / bscale;
                }

                jm1 = ilast + A_size_idx_0 * ilast;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0F) {
                  ad22_re = anorm / bscale;
                  ad22_im = 0.0F;
                } else if (anorm == 0.0F) {
                  ad22_re = 0.0F;
                  ad22_im = scale / bscale;
                } else {
                  ad22_re = anorm / bscale;
                  ad22_im = scale / bscale;
                }

                t1_re = 0.5F * (shift.re + ad22_re);
                t1_im = 0.5F * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                jm1 = ilastm1 + A_size_idx_0 * ilast;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0F) {
                  ascale_re = anorm / bscale;
                  ascale_im = 0.0F;
                } else if (anorm == 0.0F) {
                  ascale_re = 0.0F;
                  ascale_im = scale / bscale;
                } else {
                  ascale_re = anorm / bscale;
                  ascale_im = scale / bscale;
                }

                jm1 = ilast + A_size_idx_0 * ilastm1;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0F) {
                  temp2 = anorm / bscale;
                  anorm = 0.0F;
                } else if (anorm == 0.0F) {
                  temp2 = 0.0F;
                  anorm = scale / bscale;
                } else {
                  temp2 = anorm / bscale;
                  anorm = scale / bscale;
                }

                scale = shift.re * ad22_re - shift.im * ad22_im;
                temp1 = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (ascale_re * temp2
                  - ascale_im * anorm)) - scale;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (ascale_re * anorm +
                  ascale_im * temp2)) - temp1;
                c_sqrt(&shift);
                if ((t1_re - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                    0.0F) {
                  shift.re += t1_re;
                  shift.im += t1_im;
                } else {
                  shift.re = t1_re - shift.re;
                  shift.im = t1_im - shift.im;
                }
              } else {
                jm1 = ilast + A_size_idx_0 * ilastm1;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0F) {
                  ascale_re = anorm / bscale;
                  ascale_im = 0.0F;
                } else if (anorm == 0.0F) {
                  ascale_re = 0.0F;
                  ascale_im = scale / bscale;
                } else {
                  ascale_re = anorm / bscale;
                  ascale_im = scale / bscale;
                }

                eshift_re += ascale_re;
                eshift_im += ascale_im;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              n = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                jm1 = A_size_idx_0 * j;
                ctemp_tmp = j + jm1;
                ctemp.re = ascale * b_A_data[ctemp_tmp].re - shift.re * bscale;
                ctemp.im = ascale * b_A_data[ctemp_tmp].im - shift.im * bscale;
                anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
                jm1 += n;
                temp2 = ascale * (std::abs(b_A_data[jm1].re) + std::abs
                                  (b_A_data[jm1].im));
                scale = anorm;
                if (temp2 > anorm) {
                  scale = temp2;
                }

                if ((scale < 1.0F) && (scale != 0.0F)) {
                  anorm /= scale;
                  temp2 /= scale;
                }

                ctemp_tmp = j + A_size_idx_0 * (j - 1);
                if ((std::abs(b_A_data[ctemp_tmp].re) + std::abs
                     (b_A_data[ctemp_tmp].im)) * temp2 <= anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  n = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp_tmp = (ifirst + A_size_idx_0 * (ifirst - 1)) - 1;
                ctemp.re = ascale * b_A_data[ctemp_tmp].re - shift.re * bscale;
                ctemp.im = ascale * b_A_data[ctemp_tmp].im - shift.im * bscale;
              }

              goto90 = false;
              jm1 = istart + A_size_idx_0 * (istart - 1);
              b_ascale.re = ascale * b_A_data[jm1].re;
              b_ascale.im = ascale * b_A_data[jm1].im;
              d_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              jm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  c_xzlartg(b_A_data[(j + A_size_idx_0 * jm1) - 1], b_A_data[j +
                            A_size_idx_0 * jm1], &anorm, &shift, &b_A_data[(j +
                             A_size_idx_0 * jm1) - 1]);
                  ctemp_tmp = j + A_size_idx_0 * jm1;
                  b_A_data[ctemp_tmp].re = 0.0F;
                  b_A_data[ctemp_tmp].im = 0.0F;
                }

                for (n = j; n <= ilastm; n++) {
                  jm1 = j + A_size_idx_0 * (n - 1);
                  ctemp_tmp = jm1 - 1;
                  ad22_re = anorm * b_A_data[ctemp_tmp].re + (shift.re *
                    b_A_data[jm1].re - shift.im * b_A_data[jm1].im);
                  ad22_im = anorm * b_A_data[ctemp_tmp].im + (shift.re *
                    b_A_data[jm1].im + shift.im * b_A_data[jm1].re);
                  scale = b_A_data[ctemp_tmp].re;
                  b_A_data[jm1].re = anorm * b_A_data[jm1].re - (shift.re *
                    b_A_data[ctemp_tmp].re + shift.im * b_A_data[ctemp_tmp].im);
                  b_A_data[jm1].im = anorm * b_A_data[jm1].im - (shift.re *
                    b_A_data[ctemp_tmp].im - shift.im * scale);
                  b_A_data[ctemp_tmp].re = ad22_re;
                  b_A_data[ctemp_tmp].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                n = j;
                if (ilast + 1 < j + 2) {
                  n = ilast - 1;
                }

                for (i = ifirst; i <= n + 2; i++) {
                  jm1 = (i + A_size_idx_0 * (j - 1)) - 1;
                  ctemp_tmp = (i + A_size_idx_0 * j) - 1;
                  ad22_re = anorm * b_A_data[ctemp_tmp].re + (shift.re *
                    b_A_data[jm1].re - shift.im * b_A_data[jm1].im);
                  ad22_im = anorm * b_A_data[ctemp_tmp].im + (shift.re *
                    b_A_data[jm1].im + shift.im * b_A_data[jm1].re);
                  scale = b_A_data[ctemp_tmp].re;
                  b_A_data[jm1].re = anorm * b_A_data[jm1].re - (shift.re *
                    b_A_data[ctemp_tmp].re + shift.im * b_A_data[ctemp_tmp].im);
                  b_A_data[jm1].im = anorm * b_A_data[jm1].im - (shift.re *
                    b_A_data[ctemp_tmp].im - shift.im * scale);
                  b_A_data[ctemp_tmp].re = ad22_re;
                  b_A_data[ctemp_tmp].im = ad22_im;
                }

                jm1 = j - 1;
                j++;
              }
            }

            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (firstNonZero) {
      *info = ilast + 1;
      if (0 <= ilast) {
        std::memset(&alpha1_data[0], 0, (ilast + 1) * sizeof(creal32_T));
        std::memset(&beta1_data[0], 0, (ilast + 1) * sizeof(creal32_T));
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1_data[j] = b_A_data[j + A_size_idx_0 * j];
    }
  }
}

//
// Arguments    : int m
//                int n
//                int iv0
//                float tau
//                float C[100]
//                int ic0
//                float work[10]
// Return Type  : void
//
static void b_xzlarf(int m, int n, int iv0, float tau, float C[100], int ic0,
                     float work[10])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int jy;
  int b_i;
  int j;
  int ia;
  int ix;
  int exitg1;
  float c;
  if (tau != 0.0F) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C[i - 2] == 0.0F)) {
      lastv--;
      i--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      i = ic0 + lastc * 10;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0F) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    if (lastc + 1 != 0) {
      if (0 <= lastc) {
        std::memset(&work[0], 0, (lastc + 1) * sizeof(float));
      }

      i = 0;
      b_i = ic0 + 10 * lastc;
      for (jy = ic0; jy <= b_i; jy += 10) {
        ix = iv0;
        c = 0.0F;
        j = (jy + lastv) - 1;
        for (ia = jy; ia <= j; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[i] += c;
        i++;
      }
    }

    if (-tau != 0.0F) {
      i = ic0;
      jy = 0;
      for (j = 0; j <= lastc; j++) {
        if (work[jy] != 0.0F) {
          c = work[jy] * -tau;
          ix = iv0;
          b_i = lastv + i;
          for (ia = i; ia < b_i; ia++) {
            C[ia - 1] += C[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += 10;
      }
    }
  }
}

//
// Arguments    : const creal_T f
//                const creal_T g
//                double *cs
//                creal_T *sn
// Return Type  : void
//
static void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn)
{
  double scale_tmp;
  double f2;
  double scale;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  boolean_T guard1 = false;
  double g2;
  double g2s;
  scale_tmp = std::abs(f.re);
  f2 = std::abs(f.im);
  if (f2 > scale_tmp) {
    scale_tmp = f2;
  }

  f2 = std::abs(g.re);
  scale = std::abs(g.im);
  if (scale > f2) {
    f2 = scale;
  }

  scale = scale_tmp;
  if (f2 > scale_tmp) {
    scale = f2;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0 > g2) {
      scale = 1.0;
    }

    if (f2 <= scale * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        g2 = rt_hypotd(gs_re, gs_im);
        sn->re = gs_re / g2;
        sn->im = -gs_im / g2;
      } else {
        g2s = std::sqrt(g2);
        *cs = rt_hypotd(fs_re, fs_im) / g2s;
        if (scale_tmp > 1.0) {
          g2 = rt_hypotd(f.re, f.im);
          fs_re = f.re / g2;
          fs_im = f.im / g2;
        } else {
          f2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          g2 = rt_hypotd(f2, scale);
          fs_re = f2 / g2;
          fs_im = scale / g2;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      scale = std::sqrt(g2 / f2 + 1.0);
      *cs = 1.0 / scale;
      g2 += f2;
      fs_re = scale * fs_re / g2;
      fs_im = scale * fs_im / g2;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

//
// Arguments    : const creal32_T A[100]
//                creal32_T V[100]
// Return Type  : void
//
static void b_xztgevc(const creal32_T A[100], creal32_T V[100])
{
  int i;
  float anorm;
  float rworka[10];
  int j;
  float xmx;
  int re_tmp;
  float ascale;
  float d_re;
  int je;
  int x_tmp_tmp;
  float temp;
  float salpha_re;
  float salpha_im;
  float acoeff;
  boolean_T lscalea;
  float z;
  boolean_T lscaleb;
  float scale;
  creal32_T work1[10];
  float dmin;
  int b_i;
  int jr;
  creal32_T work2[10];
  float d_im;
  float brm;
  for (i = 0; i < 10; i++) {
    rworka[i] = 0.0F;
  }

  anorm = std::abs(A[0].re) + std::abs(A[0].im);
  for (j = 0; j < 9; j++) {
    for (i = 0; i <= j; i++) {
      re_tmp = i + 10 * (j + 1);
      rworka[j + 1] += std::abs(A[re_tmp].re) + std::abs(A[re_tmp].im);
    }

    i = (j + 10 * (j + 1)) + 1;
    d_re = rworka[j + 1] + (std::abs(A[i].re) + std::abs(A[i].im));
    if (d_re > anorm) {
      anorm = d_re;
    }
  }

  xmx = anorm;
  if (1.17549435E-38F > anorm) {
    xmx = 1.17549435E-38F;
  }

  ascale = 1.0F / xmx;
  for (je = 0; je < 10; je++) {
    x_tmp_tmp = 10 * (9 - je);
    i = (x_tmp_tmp - je) + 9;
    xmx = (std::abs(A[i].re) + std::abs(A[i].im)) * ascale;
    if (1.0F > xmx) {
      xmx = 1.0F;
    }

    temp = 1.0F / xmx;
    salpha_re = ascale * (temp * A[i].re);
    salpha_im = ascale * (temp * A[i].im);
    acoeff = temp * ascale;
    if ((temp >= 1.17549435E-38F) && (acoeff < 9.86076132E-31F)) {
      lscalea = true;
    } else {
      lscalea = false;
    }

    z = std::abs(salpha_re) + std::abs(salpha_im);
    if ((z >= 1.17549435E-38F) && (z < 9.86076132E-31F)) {
      lscaleb = true;
    } else {
      lscaleb = false;
    }

    scale = 1.0F;
    if (lscalea) {
      xmx = anorm;
      if (1.0141205E+30F < anorm) {
        xmx = 1.0141205E+30F;
      }

      scale = 9.86076132E-31F / temp * xmx;
    }

    if (lscaleb) {
      d_re = 9.86076132E-31F / z;
      if (d_re > scale) {
        scale = d_re;
      }
    }

    if (lscalea || lscaleb) {
      xmx = acoeff;
      if (1.0F > acoeff) {
        xmx = 1.0F;
      }

      if (z > xmx) {
        xmx = z;
      }

      d_re = 1.0F / (1.17549435E-38F * xmx);
      if (d_re < scale) {
        scale = d_re;
      }

      if (lscalea) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      salpha_re *= scale;
      salpha_im *= scale;
    }

    std::memset(&work1[0], 0, 10U * sizeof(creal32_T));
    work1[9 - je].re = 1.0F;
    work1[9 - je].im = 0.0F;
    d_re = 1.1920929E-7F * (std::abs(salpha_re) + std::abs(salpha_im));
    dmin = 1.1920929E-7F * acoeff * anorm;
    if (d_re > dmin) {
      dmin = d_re;
    }

    if (1.17549435E-38F > dmin) {
      dmin = 1.17549435E-38F;
    }

    b_i = 8 - je;
    for (jr = 0; jr <= b_i; jr++) {
      i = jr + x_tmp_tmp;
      work1[jr].re = acoeff * A[i].re;
      work1[jr].im = acoeff * A[i].im;
    }

    work1[9 - je].re = 1.0F;
    work1[9 - je].im = 0.0F;
    b_i = static_cast<int>((((-1.0 - ((-static_cast<double>(je) + 10.0) - 1.0))
      + 1.0) / -1.0));
    for (j = 0; j < b_i; j++) {
      re_tmp = 8 - (je + j);
      i = re_tmp + 10 * re_tmp;
      d_re = acoeff * A[i].re - salpha_re;
      d_im = acoeff * A[i].im - salpha_im;
      if (std::abs(d_re) + std::abs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0F;
      }

      brm = std::abs(d_re);
      scale = std::abs(d_im);
      xmx = brm + scale;
      if (xmx < 1.0F) {
        z = std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im);
        if (z >= 8.5070593E+36F * xmx) {
          temp = 1.0F / z;
          i = 9 - je;
          for (jr = 0; jr <= i; jr++) {
            work1[jr].re *= temp;
            work1[jr].im *= temp;
          }
        }
      }

      if (d_im == 0.0F) {
        if (-work1[re_tmp].im == 0.0F) {
          scale = -work1[re_tmp].re / d_re;
          xmx = 0.0F;
        } else if (-work1[re_tmp].re == 0.0F) {
          scale = 0.0F;
          xmx = -work1[re_tmp].im / d_re;
        } else {
          scale = -work1[re_tmp].re / d_re;
          xmx = -work1[re_tmp].im / d_re;
        }
      } else if (d_re == 0.0F) {
        if (-work1[re_tmp].re == 0.0F) {
          scale = -work1[re_tmp].im / d_im;
          xmx = 0.0F;
        } else if (-work1[re_tmp].im == 0.0F) {
          scale = 0.0F;
          xmx = -(-work1[re_tmp].re / d_im);
        } else {
          scale = -work1[re_tmp].im / d_im;
          xmx = -(-work1[re_tmp].re / d_im);
        }
      } else if (brm > scale) {
        z = d_im / d_re;
        xmx = d_re + z * d_im;
        scale = (-work1[re_tmp].re + z * -work1[re_tmp].im) / xmx;
        xmx = (-work1[re_tmp].im - z * -work1[re_tmp].re) / xmx;
      } else if (scale == brm) {
        if (d_re > 0.0F) {
          z = 0.5F;
        } else {
          z = -0.5F;
        }

        if (d_im > 0.0F) {
          xmx = 0.5F;
        } else {
          xmx = -0.5F;
        }

        scale = (-work1[re_tmp].re * z + -work1[re_tmp].im * xmx) / brm;
        xmx = (-work1[re_tmp].im * z - -work1[re_tmp].re * xmx) / brm;
      } else {
        z = d_re / d_im;
        xmx = d_im + z * d_re;
        scale = (z * -work1[re_tmp].re + -work1[re_tmp].im) / xmx;
        xmx = (z * -work1[re_tmp].im - (-work1[re_tmp].re)) / xmx;
      }

      work1[re_tmp].re = scale;
      work1[re_tmp].im = xmx;
      if (re_tmp + 1 > 1) {
        if (std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im) > 1.0F) {
          temp = 1.0F / (std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im));
          if (acoeff * rworka[re_tmp] >= 8.5070593E+36F * temp) {
            i = 9 - je;
            for (jr = 0; jr <= i; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        d_re = acoeff * work1[re_tmp].re;
        d_im = acoeff * work1[re_tmp].im;
        for (jr = 0; jr < re_tmp; jr++) {
          i = jr + 10 * re_tmp;
          work1[jr].re += d_re * A[i].re - d_im * A[i].im;
          work1[jr].im += d_re * A[i].im + d_im * A[i].re;
        }
      }
    }

    std::memset(&work2[0], 0, 10U * sizeof(creal32_T));
    b_i = 9 - je;
    for (i = 0; i <= b_i; i++) {
      xmx = work1[i].re;
      z = work1[i].im;
      for (jr = 0; jr < 10; jr++) {
        re_tmp = jr + 10 * i;
        work2[jr].re += V[re_tmp].re * xmx - V[re_tmp].im * z;
        work2[jr].im += V[re_tmp].re * z + V[re_tmp].im * xmx;
      }
    }

    xmx = std::abs(work2[0].re) + std::abs(work2[0].im);
    for (jr = 0; jr < 9; jr++) {
      d_re = std::abs(work2[jr + 1].re) + std::abs(work2[jr + 1].im);
      if (d_re > xmx) {
        xmx = d_re;
      }
    }

    if (xmx > 1.17549435E-38F) {
      temp = 1.0F / xmx;
      for (jr = 0; jr < 10; jr++) {
        b_i = jr + x_tmp_tmp;
        V[b_i].re = temp * work2[jr].re;
        V[b_i].im = temp * work2[jr].im;
      }
    } else {
      std::memset(&V[x_tmp_tmp], 0, 10U * sizeof(creal32_T));
    }
  }
}

//
// Arguments    : const float A[400]
//                float B[200]
// Return Type  : void
//
static void c_mldivide(const float A[400], float B[200])
{
  float b_A[400];
  int ipiv[20];
  int info;
  int jBcol;
  int temp_tmp;
  float temp;
  int kAcol;
  int i;
  int i1;
  int b_i;
  int i2;
  std::memcpy(&b_A[0], &A[0], 400U * sizeof(float));
  b_xzgetrf(b_A, ipiv, &info);
  for (info = 0; info < 19; info++) {
    if (ipiv[info] != info + 1) {
      for (jBcol = 0; jBcol < 10; jBcol++) {
        temp_tmp = info + 20 * jBcol;
        temp = B[temp_tmp];
        i = (ipiv[info] + 20 * jBcol) - 1;
        B[temp_tmp] = B[i];
        B[i] = temp;
      }
    }
  }

  for (info = 0; info < 10; info++) {
    jBcol = 20 * info;
    for (temp_tmp = 0; temp_tmp < 20; temp_tmp++) {
      kAcol = 20 * temp_tmp;
      i = temp_tmp + jBcol;
      if (B[i] != 0.0F) {
        i1 = temp_tmp + 2;
        for (b_i = i1; b_i < 21; b_i++) {
          i2 = (b_i + jBcol) - 1;
          B[i2] -= B[i] * b_A[(b_i + kAcol) - 1];
        }
      }
    }
  }

  for (info = 0; info < 10; info++) {
    jBcol = 20 * info;
    for (temp_tmp = 19; temp_tmp >= 0; temp_tmp--) {
      kAcol = 20 * temp_tmp;
      i = temp_tmp + jBcol;
      if (B[i] != 0.0F) {
        B[i] /= b_A[temp_tmp + kAcol];
        for (b_i = 0; b_i < temp_tmp; b_i++) {
          i1 = b_i + jBcol;
          B[i1] -= B[i] * b_A[b_i + kAcol];
        }
      }
    }
  }
}

//
// Arguments    : const float x[8]
// Return Type  : float
//
static float c_norm(const float x[8])
{
  float y;
  int kase;
  float s[3];
  float A[8];
  float e[4];
  float work[2];
  int q;
  int qjj;
  int qp1;
  int qq_tmp;
  int qq;
  boolean_T apply_transform;
  int iter;
  float b;
  float nrm;
  float rt;
  float snorm;
  int ix;
  int iy;
  int k;
  int exitg1;
  boolean_T exitg2;
  float scale;
  float sm;
  float sqds;
  for (kase = 0; kase < 8; kase++) {
    A[kase] = x[kase];
  }

  s[0] = 0.0F;
  e[0] = 0.0F;
  e[1] = 0.0F;
  e[2] = 0.0F;
  e[3] = 0.0F;
  work[0] = 0.0F;
  for (q = 0; q < 2; q++) {
    qp1 = q + 2;
    qq_tmp = q + (q << 1);
    qq = qq_tmp + 1;
    apply_transform = false;
    if (q + 1 <= 1) {
      nrm = p_xnrm2(2, A, qq);
      if (nrm > 0.0F) {
        apply_transform = true;
        if (A[qq - 1] < 0.0F) {
          b = -nrm;
        } else {
          b = nrm;
        }

        if (std::abs(b) >= 9.86076132E-32F) {
          nrm = 1.0F / b;
          kase = qq + 1;
          for (k = qq; k <= kase; k++) {
            A[k - 1] *= nrm;
          }
        } else {
          kase = qq + 1;
          for (k = qq; k <= kase; k++) {
            A[k - 1] /= b;
          }
        }

        A[qq - 1]++;
        s[0] = -b;
      } else {
        s[0] = 0.0F;
      }
    }

    for (iter = qp1; iter < 5; iter++) {
      qjj = q + ((iter - 1) << 1);
      if (apply_transform) {
        kase = 1 - q;
        ix = qq;
        iy = qjj;
        nrm = 0.0F;
        for (k = 0; k <= kase; k++) {
          nrm += A[ix - 1] * A[iy];
          ix++;
          iy++;
        }

        nrm = -(nrm / A[qq_tmp]);
        if (nrm != 0.0F) {
          ix = qq - 1;
          iy = qjj;
          kase = 1 - q;
          for (k = 0; k <= kase; k++) {
            A[iy] += nrm * A[ix];
            ix++;
            iy++;
          }
        }
      }

      e[iter - 1] = A[qjj];
    }

    nrm = q_xnrm2(3 - q, e, q + 2);
    if (nrm == 0.0F) {
      e[q] = 0.0F;
    } else {
      if (e[q + 1] < 0.0F) {
        e[q] = -nrm;
      } else {
        e[q] = nrm;
      }

      nrm = e[q];
      if (std::abs(e[q]) >= 9.86076132E-32F) {
        nrm = 1.0F / e[q];
        for (k = qp1; k < 5; k++) {
          e[k - 1] *= nrm;
        }
      } else {
        for (k = qp1; k < 5; k++) {
          e[k - 1] /= nrm;
        }
      }

      e[q + 1]++;
      e[q] = -e[q];
      if (q + 2 <= 2) {
        b = 0.0F;
        work[1] = 0.0F;
        if ((1 - q >= 1) && (e[1] != 0.0F)) {
          b = e[1] * A[3];
          work[1] = b;
        }

        if ((1 - q >= 1) && (e[2] != 0.0F)) {
          b += e[2] * A[5];
          work[1] = b;
        }

        if ((1 - q >= 1) && (e[3] != 0.0F)) {
          b += e[3] * A[7];
          work[1] = b;
        }

        b_xaxpy(1 - q, -e[1] / e[1], work, A, 4);
        b_xaxpy(1 - q, -e[2] / e[1], work, A, 6);
        b_xaxpy(1 - q, -e[3] / e[1], work, A, 8);
      }
    }
  }

  qjj = 1;
  s[1] = A[3];
  s[2] = 0.0F;
  e[2] = 0.0F;
  iter = 0;
  b = s[0];
  if (s[0] != 0.0F) {
    rt = std::abs(s[0]);
    nrm = s[0] / rt;
    b = rt;
    s[0] = rt;
    e[0] /= nrm;
  }

  if (e[0] != 0.0F) {
    rt = std::abs(e[0]);
    y = e[0];
    e[0] = rt;
    s[1] = A[3] * (rt / y);
  }

  nrm = std::abs(b);
  snorm = e[0];
  if (nrm > snorm) {
    snorm = nrm;
  }

  b = s[1];
  if (s[1] != 0.0F) {
    rt = std::abs(s[1]);
    nrm = s[1] / rt;
    b = rt;
    s[1] = rt;
    e[1] /= nrm;
  }

  if (e[1] != 0.0F) {
    e[1] = std::abs(e[1]);
    s[2] = 0.0F;
  }

  nrm = std::abs(b);
  rt = e[1];
  if (nrm > rt) {
    rt = nrm;
  }

  if (snorm <= rt) {
    snorm = rt;
  }

  if (snorm <= 0.0F) {
    snorm = 0.0F;
  }

  while ((qjj + 2 > 0) && (iter < 75)) {
    kase = qjj;
    do {
      exitg1 = 0;
      q = kase + 1;
      if (kase + 1 == 0) {
        exitg1 = 1;
      } else {
        nrm = std::abs(e[kase]);
        if ((nrm <= 1.1920929E-7F * (std::abs(s[kase]) + std::abs(s[kase + 1])))
            || (nrm <= 9.86076132E-32F) || ((iter > 20) && (nrm <= 1.1920929E-7F
              * snorm))) {
          e[kase] = 0.0F;
          exitg1 = 1;
        } else {
          kase--;
        }
      }
    } while (exitg1 == 0);

    if (kase + 1 == qjj + 1) {
      kase = 4;
    } else {
      qq = qjj + 2;
      qq_tmp = qjj + 2;
      exitg2 = false;
      while ((!exitg2) && (qq_tmp >= kase + 1)) {
        qq = qq_tmp;
        if (qq_tmp == kase + 1) {
          exitg2 = true;
        } else {
          nrm = 0.0F;
          if (qq_tmp < qjj + 2) {
            nrm = std::abs(e[qq_tmp - 1]);
          }

          if (qq_tmp > kase + 2) {
            nrm += std::abs(e[qq_tmp - 2]);
          }

          rt = std::abs(s[qq_tmp - 1]);
          if ((rt <= 1.1920929E-7F * nrm) || (rt <= 9.86076132E-32F)) {
            s[qq_tmp - 1] = 0.0F;
            exitg2 = true;
          } else {
            qq_tmp--;
          }
        }
      }

      if (qq == kase + 1) {
        kase = 3;
      } else if (qq == qjj + 2) {
        kase = 1;
      } else {
        kase = 2;
        q = qq;
      }
    }

    switch (kase) {
     case 1:
      rt = e[qjj];
      e[qjj] = 0.0F;
      kase = qjj + 1;
      for (k = kase; k >= q + 1; k--) {
        b_xrotg(&s[k - 1], &rt, &sm, &sqds);
        if (k > q + 1) {
          rt = -sqds * e[0];
          e[0] *= sm;
        }
      }
      break;

     case 2:
      rt = e[q - 1];
      e[q - 1] = 0.0F;
      for (k = q + 1; k <= qjj + 2; k++) {
        b_xrotg(&s[k - 1], &rt, &sm, &sqds);
        b = e[k - 1];
        rt = -sqds * b;
        e[k - 1] = b * sm;
      }
      break;

     case 3:
      kase = qjj + 1;
      nrm = s[qjj + 1];
      scale = std::abs(nrm);
      rt = std::abs(s[qjj]);
      if (scale <= rt) {
        scale = rt;
      }

      rt = std::abs(e[qjj]);
      if (scale <= rt) {
        scale = rt;
      }

      rt = std::abs(s[q]);
      if (scale <= rt) {
        scale = rt;
      }

      rt = std::abs(e[q]);
      if (scale <= rt) {
        scale = rt;
      }

      sm = nrm / scale;
      nrm = s[qjj] / scale;
      rt = e[qjj] / scale;
      sqds = s[q] / scale;
      b = ((nrm + sm) * (nrm - sm) + rt * rt) / 2.0F;
      nrm = sm * rt;
      nrm *= nrm;
      if ((b != 0.0F) || (nrm != 0.0F)) {
        rt = std::sqrt(b * b + nrm);
        if (b < 0.0F) {
          rt = -rt;
        }

        rt = nrm / (b + rt);
      } else {
        rt = 0.0F;
      }

      rt += (sqds + sm) * (sqds - sm);
      nrm = sqds * (e[q] / scale);
      for (k = q + 1; k <= kase; k++) {
        b_xrotg(&rt, &nrm, &sm, &sqds);
        if (k > q + 1) {
          e[0] = rt;
        }

        nrm = e[k - 1];
        b = s[k - 1];
        e[k - 1] = sm * nrm - sqds * b;
        rt = sqds * s[k];
        s[k] *= sm;
        s[k - 1] = sm * b + sqds * nrm;
        b_xrotg(&s[k - 1], &rt, &sm, &sqds);
        rt = sm * e[k - 1] + sqds * s[k];
        s[k] = -sqds * e[k - 1] + sm * s[k];
        nrm = sqds * e[k];
        e[k] *= sm;
      }

      e[qjj] = rt;
      iter++;
      break;

     default:
      if (s[q] < 0.0F) {
        s[q] = -s[q];
      }

      qp1 = q + 1;
      while ((q + 1 < 3) && (s[q] < s[qp1])) {
        rt = s[q];
        s[q] = s[qp1];
        s[qp1] = rt;
        q = qp1;
        qp1++;
      }

      iter = 0;
      qjj--;
      break;
    }
  }

  return s[0];
}

//
// Arguments    : const double A[32]
//                double Q[64]
//                double R[32]
// Return Type  : void
//
static void c_qr(const double A[32], double Q[64], double R[32])
{
  int coltop;
  double tau[8];
  int i;
  double work[8];
  int b_i;
  int ii;
  double atmp;
  double c;
  double beta1;
  int lastv;
  int lastc;
  int knt;
  int c_i;
  boolean_T exitg2;
  int ia;
  int exitg1;
  int ix;
  int i1;
  for (coltop = 0; coltop < 4; coltop++) {
    for (i = 0; i < 8; i++) {
      b_i = i + (coltop << 3);
      Q[b_i] = A[b_i];
      Q[i + ((coltop + 4) << 3)] = 0.0;
    }
  }

  std::memset(&tau[0], 0, 8U * sizeof(double));
  std::memset(&work[0], 0, 8U * sizeof(double));
  for (i = 0; i < 8; i++) {
    ii = (i << 3) + i;
    if (i + 1 < 8) {
      atmp = Q[ii];
      coltop = ii + 2;
      tau[i] = 0.0;
      c = g_xnrm2(7 - i, Q, ii + 2);
      if (c != 0.0) {
        beta1 = rt_hypotd(Q[ii], c);
        if (Q[ii] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          c_i = (ii - i) + 8;
          do {
            knt++;
            for (b_i = coltop; b_i <= c_i; b_i++) {
              Q[b_i - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd(atmp, g_xnrm2(7 - i, Q, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          c = 1.0 / (atmp - beta1);
          for (b_i = coltop; b_i <= c_i; b_i++) {
            Q[b_i - 1] *= c;
          }

          for (b_i = 0; b_i <= knt; b_i++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - Q[ii]) / beta1;
          c = 1.0 / (Q[ii] - beta1);
          c_i = (ii - i) + 8;
          for (b_i = coltop; b_i <= c_i; b_i++) {
            Q[b_i - 1] *= c;
          }

          atmp = beta1;
        }
      }

      Q[ii] = 1.0;
      if (tau[i] != 0.0) {
        lastv = 8 - i;
        b_i = (ii - i) + 7;
        while ((lastv > 0) && (Q[b_i] == 0.0)) {
          lastv--;
          b_i--;
        }

        lastc = 7 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          coltop = (ii + ((lastc - 1) << 3)) + 8;
          ia = coltop;
          do {
            exitg1 = 0;
            if (ia + 1 <= coltop + lastv) {
              if (Q[ia] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        b_i = ii + 9;
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(double));
          }

          knt = 0;
          c_i = (ii + ((lastc - 1) << 3)) + 9;
          for (coltop = b_i; coltop <= c_i; coltop += 8) {
            ix = ii;
            c = 0.0;
            i1 = (coltop + lastv) - 1;
            for (ia = coltop; ia <= i1; ia++) {
              c += Q[ia - 1] * Q[ix];
              ix++;
            }

            work[knt] += c;
            knt++;
          }
        }

        c_xgerc(lastv, lastc, -tau[i], ii + 1, work, Q, ii + 9);
      }

      Q[ii] = atmp;
    } else {
      tau[7] = 0.0;
    }
  }

  R[0] = Q[0];
  for (i = 2; i < 9; i++) {
    R[i - 1] = 0.0;
  }

  for (i = 0; i < 2; i++) {
    R[i + 8] = Q[i + 8];
  }

  for (i = 3; i < 9; i++) {
    R[i + 7] = 0.0;
  }

  for (i = 0; i < 3; i++) {
    R[i + 16] = Q[i + 16];
  }

  for (i = 4; i < 9; i++) {
    R[i + 15] = 0.0;
  }

  for (i = 0; i < 4; i++) {
    R[i + 24] = Q[i + 24];
  }

  for (i = 5; i < 9; i++) {
    R[i + 23] = 0.0;
  }

  for (coltop = 0; coltop < 4; coltop++) {
    ia = (coltop + 4) << 3;
    std::memset(&Q[ia], 0, 8U * sizeof(double));
    Q[(ia + coltop) + 4] = 1.0;
  }

  std::memset(&work[0], 0, 8U * sizeof(double));
  Q[27] = 1.0;
  if (tau[3] != 0.0) {
    lastv = 5;
    i = 33;
    while ((lastv > 0) && (Q[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = 4;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 36;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(double));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 36;
      for (coltop = 36; coltop <= c_i; coltop += 8) {
        ix = 36;
        c = 0.0;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    c_xgerc(lastv, lastc, -tau[3], 28, work, Q, 36);
  }

  for (b_i = 29; b_i < 33; b_i++) {
    Q[b_i - 1] *= -tau[3];
  }

  Q[27] = 1.0 - tau[3];
  for (coltop = 0; coltop < 3; coltop++) {
    Q[26 - coltop] = 0.0;
  }

  Q[18] = 1.0;
  if (tau[2] != 0.0) {
    lastv = 6;
    i = 25;
    while ((lastv > 0) && (Q[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = 5;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 27;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(double));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 27;
      for (coltop = 27; coltop <= c_i; coltop += 8) {
        ix = 27;
        c = 0.0;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    c_xgerc(lastv, lastc, -tau[2], 19, work, Q, 27);
  }

  for (b_i = 20; b_i < 25; b_i++) {
    Q[b_i - 1] *= -tau[2];
  }

  Q[18] = 1.0 - tau[2];
  for (coltop = 0; coltop < 2; coltop++) {
    Q[17 - coltop] = 0.0;
  }

  Q[9] = 1.0;
  if (tau[1] != 0.0) {
    lastv = 7;
    i = 17;
    while ((lastv > 0) && (Q[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = 6;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 18;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(double));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 18;
      for (coltop = 18; coltop <= c_i; coltop += 8) {
        ix = 18;
        c = 0.0;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    c_xgerc(lastv, lastc, -tau[1], 10, work, Q, 18);
  }

  for (b_i = 11; b_i < 17; b_i++) {
    Q[b_i - 1] *= -tau[1];
  }

  Q[9] = 1.0 - tau[1];
  Q[8] = 0.0;
  Q[0] = 1.0;
  if (tau[0] != 0.0) {
    lastv = 8;
    i = 9;
    while ((lastv > 0) && (Q[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = 7;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 9;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(double));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 9;
      for (coltop = 9; coltop <= c_i; coltop += 8) {
        ix = 9;
        c = 0.0;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    c_xgerc(lastv, lastc, -tau[0], 1, work, Q, 9);
  }

  for (b_i = 2; b_i < 9; b_i++) {
    Q[b_i - 1] *= -tau[0];
  }

  Q[0] = 1.0 - tau[0];
}

//
// Arguments    : const double A[9]
// Return Type  : double
//
static double c_rcond(const double A[9])
{
  double result;
  double normA;
  int j;
  double s;
  double b_A[9];
  int ipiv[3];
  int jA;
  int exitg1;
  double ainvnm;
  int iter;
  int kase;
  int jump;
  double x[3];
  int exitg2;
  int exitg3;
  int exitg4;
  int exitg5;
  int b_j;
  int i;
  int b_i;
  int ix;
  double absrexk;
  double ainvnm_tmp;
  result = 0.0;
  normA = 0.0;
  for (j = 0; j < 3; j++) {
    s = (std::abs(A[3 * j]) + std::abs(A[3 * j + 1])) + std::abs(A[3 * j + 2]);
    if (s > normA) {
      normA = s;
    }
  }

  if (normA != 0.0) {
    std::memcpy(&b_A[0], &A[0], 9U * sizeof(double));
    c_xzgetrf(b_A, ipiv, &jA);
    jA = 2;
    do {
      exitg1 = 0;
      if (jA + 1 > 0) {
        if (b_A[jA + 3 * jA] == 0.0) {
          exitg1 = 1;
        } else {
          jA--;
        }
      } else {
        ainvnm = 0.0;
        iter = 2;
        kase = 1;
        jump = 1;
        j = 0;
        x[0] = 0.33333333333333331;
        x[1] = 0.33333333333333331;
        x[2] = 0.33333333333333331;
        do {
          exitg2 = 0;
          do {
            exitg3 = 0;
            do {
              exitg4 = 0;
              do {
                exitg5 = 0;
                if (kase == 1) {
                  for (b_j = 0; b_j < 3; b_j++) {
                    jA = b_j + b_j * 3;
                    i = 1 - b_j;
                    for (b_i = 0; b_i <= i; b_i++) {
                      ix = (b_j + b_i) + 1;
                      x[ix] -= x[b_j] * b_A[(jA + b_i) + 1];
                    }
                  }

                  for (b_j = 2; b_j >= 0; b_j--) {
                    jA = b_j + b_j * 3;
                    x[b_j] /= b_A[jA];
                    for (b_i = 0; b_i < b_j; b_i++) {
                      ix = (b_j - b_i) - 1;
                      x[ix] -= x[b_j] * b_A[(jA - b_i) - 1];
                    }
                  }
                } else {
                  for (b_j = 0; b_j < 3; b_j++) {
                    jA = b_j * 3;
                    s = x[b_j];
                    for (b_i = 0; b_i < b_j; b_i++) {
                      s -= b_A[jA + b_i] * x[b_i];
                    }

                    x[b_j] = s / b_A[jA + b_j];
                  }

                  for (b_j = 2; b_j >= 0; b_j--) {
                    jA = b_j * 3;
                    s = x[b_j];
                    i = b_j + 2;
                    for (b_i = 3; b_i >= i; b_i--) {
                      s -= b_A[(jA + b_i) - 1] * x[b_i - 1];
                    }

                    x[b_j] = s;
                  }
                }

                if (jump == 1) {
                  s = std::abs(x[0]);
                  if (s > 2.2250738585072014E-308) {
                    x[0] /= s;
                  } else {
                    x[0] = 1.0;
                  }

                  ainvnm = s + std::abs(x[1]);
                  s = std::abs(x[1]);
                  if (s > 2.2250738585072014E-308) {
                    x[1] /= s;
                  } else {
                    x[1] = 1.0;
                  }

                  ainvnm += std::abs(x[2]);
                  s = std::abs(x[2]);
                  if (s > 2.2250738585072014E-308) {
                    x[2] /= s;
                  } else {
                    x[2] = 1.0;
                  }

                  kase = 2;
                  jump = 2;
                } else {
                  exitg5 = 1;
                }
              } while (exitg5 == 0);

              if (jump == 2) {
                j = 0;
                s = std::abs(x[0]);
                absrexk = std::abs(x[1]);
                if (absrexk > s) {
                  j = 1;
                  s = absrexk;
                }

                if (std::abs(x[2]) > s) {
                  j = 2;
                }

                iter = 2;
                x[0] = 0.0;
                x[1] = 0.0;
                x[2] = 0.0;
                x[j] = 1.0;
                kase = 1;
                jump = 3;
              } else {
                exitg4 = 1;
              }
            } while (exitg4 == 0);

            if (jump == 3) {
              s = std::abs(x[0]);
              absrexk = std::abs(x[1]);
              ainvnm_tmp = std::abs(x[2]);
              ainvnm = (s + absrexk) + ainvnm_tmp;
              if (ainvnm <= x[0]) {
                x[0] = 1.0;
                x[1] = -1.5;
                x[2] = 2.0;
                kase = 1;
                jump = 5;
              } else {
                if (s > 2.2250738585072014E-308) {
                  x[0] /= s;
                } else {
                  x[0] = 1.0;
                }

                if (absrexk > 2.2250738585072014E-308) {
                  x[1] /= absrexk;
                } else {
                  x[1] = 1.0;
                }

                if (ainvnm_tmp > 2.2250738585072014E-308) {
                  x[2] /= ainvnm_tmp;
                } else {
                  x[2] = 1.0;
                }

                kase = 2;
                jump = 4;
              }
            } else {
              exitg3 = 1;
            }
          } while (exitg3 == 0);

          if (jump == 4) {
            jA = j;
            j = 0;
            s = std::abs(x[0]);
            absrexk = std::abs(x[1]);
            if (absrexk > s) {
              j = 1;
              s = absrexk;
            }

            if (std::abs(x[2]) > s) {
              j = 2;
            }

            if ((std::abs(x[jA]) != std::abs(x[j])) && (iter <= 5)) {
              iter++;
              x[0] = 0.0;
              x[1] = 0.0;
              x[2] = 0.0;
              x[j] = 1.0;
              kase = 1;
              jump = 3;
            } else {
              x[0] = 1.0;
              x[1] = -1.5;
              x[2] = 2.0;
              kase = 1;
              jump = 5;
            }
          } else {
            exitg2 = 1;
          }
        } while ((exitg2 == 0) || (!(jump == 5)));

        s = 2.0 * ((std::abs(x[0]) + std::abs(x[1])) + std::abs(x[2])) / 3.0 /
          3.0;
        if (s > ainvnm) {
          ainvnm = s;
        }

        if (ainvnm != 0.0) {
          result = 1.0 / ainvnm / normA;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return result;
}

//
// Arguments    : creal32_T *x
// Return Type  : void
//
static void c_sqrt(creal32_T *x)
{
  float xr;
  float xi;
  float absxr;
  float absxi;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0F) {
    if (xr < 0.0F) {
      absxr = 0.0F;
      absxi = std::sqrt(-xr);
    } else {
      absxr = std::sqrt(xr);
      absxi = 0.0F;
    }
  } else if (xr == 0.0F) {
    if (xi < 0.0F) {
      absxr = std::sqrt(-xi / 2.0F);
      absxi = -absxr;
    } else {
      absxr = std::sqrt(xi / 2.0F);
      absxi = absxr;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 8.50705867E+37F) || (absxi > 8.50705867E+37F)) {
      absxr *= 0.5F;
      absxi = rt_hypotf(absxr, absxi * 0.5F);
      if (absxi > absxr) {
        absxr = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0F);
      } else {
        absxr = std::sqrt(absxi) * 1.41421354F;
      }
    } else {
      absxr = std::sqrt((rt_hypotf(absxr, absxi) + absxr) * 0.5F);
    }

    if (xr > 0.0F) {
      absxi = 0.5F * (xi / absxr);
    } else {
      if (xi < 0.0F) {
        absxi = -absxr;
      } else {
        absxi = absxr;
      }

      absxr = 0.5F * (xi / absxi);
    }
  }

  x->re = absxr;
  x->im = absxi;
}

//
// Arguments    : int m
//                int n
//                double alpha1
//                int ix0
//                const double y[8]
//                double A[64]
//                int ia0
// Return Type  : void
//
static void c_xgerc(int m, int n, double alpha1, int ix0, const double y[8],
                    double A[64], int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i;
  int ijA;
  if (alpha1 != 0.0) {
    jA = ia0;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 8;
    }
  }
}

//
// Arguments    : int n
//                const double x[12]
//                int ix0
// Return Type  : double
//
static double c_xnrm2(int n, const double x[12], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                float x[100]
//                int ix0
//                int iy0
//                float c
//                float s
// Return Type  : void
//
static void c_xrot(int n, float x[100], int ix0, int iy0, float c, float s)
{
  int ix;
  int iy;
  int k;
  float temp;
  if (n >= 1) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      temp = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = temp;
      iy++;
      ix++;
    }
  }
}

//
// Arguments    : double A[9]
//                int ipiv[3]
//                int *info
// Return Type  : void
//
static void c_xzgetrf(double A[9], int ipiv[3], int *info)
{
  int j;
  int mmj_tmp;
  int b_tmp;
  int jp1j;
  int iy;
  int jA;
  int ix;
  double smax;
  int k;
  double s;
  int i;
  int ijA;
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  *info = 0;
  for (j = 0; j < 2; j++) {
    mmj_tmp = 1 - j;
    b_tmp = j << 2;
    jp1j = b_tmp + 2;
    iy = 3 - j;
    jA = 0;
    ix = b_tmp;
    smax = std::abs(A[b_tmp]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (A[b_tmp + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = iy + 1;
        smax = A[j];
        A[j] = A[iy];
        A[iy] = smax;
        ix = j + 3;
        iy += 3;
        smax = A[ix];
        A[ix] = A[iy];
        A[iy] = smax;
        ix += 3;
        iy += 3;
        smax = A[ix];
        A[ix] = A[iy];
        A[iy] = smax;
      }

      i = (b_tmp - j) + 3;
      for (iy = jp1j; iy <= i; iy++) {
        A[iy - 1] /= A[b_tmp];
      }
    } else {
      *info = j + 1;
    }

    iy = b_tmp + 3;
    jA = b_tmp;
    for (k = 0; k <= mmj_tmp; k++) {
      smax = A[iy];
      if (A[iy] != 0.0) {
        ix = b_tmp + 1;
        i = jA + 5;
        jp1j = (jA - j) + 6;
        for (ijA = i; ijA <= jp1j; ijA++) {
          A[ijA - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 3;
      jA += 3;
    }
  }

  if ((*info == 0) && (A[8] == 0.0)) {
    *info = 3;
  }
}

//
// Arguments    : creal_T A[100]
//                int ilo
//                int ihi
//                creal_T Z[100]
//                int *info
//                creal_T alpha1[10]
//                creal_T beta1[10]
// Return Type  : void
//
static void c_xzhgeqz(creal_T A[100], int ilo, int ihi, creal_T Z[100], int
                      *info, creal_T alpha1[10], creal_T beta1[10])
{
  int i;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double reAij;
  double sumsq;
  double b_atol;
  boolean_T firstNonZero;
  int j;
  int ctemp_tmp;
  double ascale;
  int jp1;
  double imAij;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
  double temp2;
  int ilast;
  int ilastm1;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T b_guard1 = false;
  boolean_T guard3 = false;
  boolean_T exitg2;
  creal_T shift;
  double ad22_re;
  double ad22_im;
  double t1_re;
  double t1_im;
  double t1_im_tmp;
  creal_T b_ascale;
  int ad22_re_tmp;
  *info = 0;
  for (i = 0; i < 10; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (ilo <= ihi) {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ctemp_tmp = j + 1;
      if (ihi < j + 1) {
        ctemp_tmp = ihi;
      }

      for (i = ilo; i <= ctemp_tmp; i++) {
        jp1 = (i + 10 * (j - 1)) - 1;
        reAij = A[jp1].re;
        imAij = A[jp1].im;
        if (reAij != 0.0) {
          anorm = std::abs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = sumsq * temp2 * temp2 + 1.0;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          anorm = std::abs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = sumsq * temp2 * temp2 + 1.0;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * std::sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  reAij = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    reAij = anorm;
  }

  ascale = 1.0 / reAij;
  firstNonZero = true;
  ctemp_tmp = ihi + 1;
  for (j = ctemp_tmp; j < 11; j++) {
    alpha1[j - 1] = A[(j + 10 * (j - 1)) - 1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          ctemp_tmp = ilast + 10 * ilastm1;
          if (std::abs(A[ctemp_tmp].re) + std::abs(A[ctemp_tmp].im) <= b_atol) {
            A[ctemp_tmp].re = 0.0;
            A[ctemp_tmp].im = 0.0;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                ctemp_tmp = j + 10 * (j - 1);
                if (std::abs(A[ctemp_tmp].re) + std::abs(A[ctemp_tmp].im) <=
                    b_atol) {
                  A[ctemp_tmp].re = 0.0;
                  A[ctemp_tmp].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }

            if (goto70) {
              b_guard1 = true;
            } else {
              std::memset(&alpha1[0], 0, 10U * sizeof(creal_T));
              std::memset(&beta1[0], 0, 10U * sizeof(creal_T));
              std::memset(&Z[0], 0, 100U * sizeof(creal_T));
              *info = 1;
              exitg1 = 1;
            }
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1[ilast] = A[ilast + 10 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                jp1 = ilastm1 + 10 * ilastm1;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0) {
                  shift.re = anorm / 0.31622776601683794;
                  shift.im = 0.0;
                } else if (anorm == 0.0) {
                  shift.re = 0.0;
                  shift.im = reAij / 0.31622776601683794;
                } else {
                  shift.re = anorm / 0.31622776601683794;
                  shift.im = reAij / 0.31622776601683794;
                }

                jp1 = ilast + 10 * ilast;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0) {
                  ad22_re = anorm / 0.31622776601683794;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = reAij / 0.31622776601683794;
                } else {
                  ad22_re = anorm / 0.31622776601683794;
                  ad22_im = reAij / 0.31622776601683794;
                }

                t1_re = 0.5 * (shift.re + ad22_re);
                t1_im = 0.5 * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                jp1 = ilastm1 + 10 * ilast;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0) {
                  imAij = anorm / 0.31622776601683794;
                  temp2 = 0.0;
                } else if (anorm == 0.0) {
                  imAij = 0.0;
                  temp2 = reAij / 0.31622776601683794;
                } else {
                  imAij = anorm / 0.31622776601683794;
                  temp2 = reAij / 0.31622776601683794;
                }

                jp1 = ilast + 10 * ilastm1;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0) {
                  sumsq = anorm / 0.31622776601683794;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  sumsq = 0.0;
                  anorm = reAij / 0.31622776601683794;
                } else {
                  sumsq = anorm / 0.31622776601683794;
                  anorm = reAij / 0.31622776601683794;
                }

                reAij = shift.re * ad22_re - shift.im * ad22_im;
                scale = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (imAij * sumsq -
                  temp2 * anorm)) - reAij;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (imAij * anorm + temp2 *
                  sumsq)) - scale;
                b_sqrt(&shift);
                if ((t1_re - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                    0.0) {
                  shift.re += t1_re;
                  shift.im += t1_im;
                } else {
                  shift.re = t1_re - shift.re;
                  shift.im = t1_im - shift.im;
                }
              } else {
                jp1 = ilast + 10 * ilastm1;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0) {
                  imAij = anorm / 0.31622776601683794;
                  temp2 = 0.0;
                } else if (anorm == 0.0) {
                  imAij = 0.0;
                  temp2 = reAij / 0.31622776601683794;
                } else {
                  imAij = anorm / 0.31622776601683794;
                  temp2 = reAij / 0.31622776601683794;
                }

                eshift_re += imAij;
                eshift_im += temp2;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp_tmp = j + 10 * j;
                ctemp.re = ascale * A[ctemp_tmp].re - shift.re *
                  0.31622776601683794;
                ctemp.im = ascale * A[ctemp_tmp].im - shift.im *
                  0.31622776601683794;
                anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
                jp1 += 10 * j;
                temp2 = ascale * (std::abs(A[jp1].re) + std::abs(A[jp1].im));
                reAij = anorm;
                if (temp2 > anorm) {
                  reAij = temp2;
                }

                if ((reAij < 1.0) && (reAij != 0.0)) {
                  anorm /= reAij;
                  temp2 /= reAij;
                }

                ctemp_tmp = j + 10 * (j - 1);
                if ((std::abs(A[ctemp_tmp].re) + std::abs(A[ctemp_tmp].im)) *
                    temp2 <= anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp_tmp = (ifirst + 10 * (ifirst - 1)) - 1;
                ctemp.re = ascale * A[ctemp_tmp].re - shift.re *
                  0.31622776601683794;
                ctemp.im = ascale * A[ctemp_tmp].im - shift.im *
                  0.31622776601683794;
              }

              goto90 = false;
              jp1 = istart + 10 * (istart - 1);
              b_ascale.re = ascale * A[jp1].re;
              b_ascale.im = ascale * A[jp1].im;
              b_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              jp1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  ctemp_tmp = j + 10 * jp1;
                  xzlartg(A[ctemp_tmp - 1], A[ctemp_tmp], &anorm, &shift, &A[(j
                           + 10 * jp1) - 1]);
                  A[ctemp_tmp].re = 0.0;
                  A[ctemp_tmp].im = 0.0;
                }

                for (jp1 = j; jp1 < 11; jp1++) {
                  ctemp_tmp = j + 10 * (jp1 - 1);
                  ad22_re_tmp = ctemp_tmp - 1;
                  ad22_re = anorm * A[ad22_re_tmp].re + (shift.re * A[ctemp_tmp]
                    .re - shift.im * A[ctemp_tmp].im);
                  ad22_im = anorm * A[ad22_re_tmp].im + (shift.re * A[ctemp_tmp]
                    .im + shift.im * A[ctemp_tmp].re);
                  reAij = A[ad22_re_tmp].im;
                  scale = A[ad22_re_tmp].re;
                  A[ctemp_tmp].re = anorm * A[ctemp_tmp].re - (shift.re *
                    A[ad22_re_tmp].re + shift.im * A[ad22_re_tmp].im);
                  A[ctemp_tmp].im = anorm * A[ctemp_tmp].im - (shift.re * reAij
                    - shift.im * scale);
                  A[ad22_re_tmp].re = ad22_re;
                  A[ad22_re_tmp].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                jp1 = j;
                if (ilast + 1 < j + 2) {
                  jp1 = ilast - 1;
                }

                for (i = 1; i <= jp1 + 2; i++) {
                  ctemp_tmp = (i + 10 * (j - 1)) - 1;
                  ad22_re_tmp = (i + 10 * j) - 1;
                  ad22_re = anorm * A[ad22_re_tmp].re + (shift.re * A[ctemp_tmp]
                    .re - shift.im * A[ctemp_tmp].im);
                  ad22_im = anorm * A[ad22_re_tmp].im + (shift.re * A[ctemp_tmp]
                    .im + shift.im * A[ctemp_tmp].re);
                  reAij = A[ad22_re_tmp].im;
                  scale = A[ad22_re_tmp].re;
                  A[ctemp_tmp].re = anorm * A[ctemp_tmp].re - (shift.re *
                    A[ad22_re_tmp].re + shift.im * A[ad22_re_tmp].im);
                  A[ctemp_tmp].im = anorm * A[ctemp_tmp].im - (shift.re * reAij
                    - shift.im * scale);
                  A[ad22_re_tmp].re = ad22_re;
                  A[ad22_re_tmp].im = ad22_im;
                }

                reAij = shift.re;
                scale = shift.im;
                for (i = 0; i < 10; i++) {
                  ctemp_tmp = i + 10 * (j - 1);
                  ad22_re_tmp = i + 10 * j;
                  ad22_re = anorm * Z[ad22_re_tmp].re + (reAij * Z[ctemp_tmp].re
                    - scale * Z[ctemp_tmp].im);
                  ad22_im = anorm * Z[ad22_re_tmp].im + (reAij * Z[ctemp_tmp].im
                    + scale * Z[ctemp_tmp].re);
                  sumsq = Z[ad22_re_tmp].im;
                  imAij = Z[ad22_re_tmp].re;
                  Z[ctemp_tmp].re = anorm * Z[ctemp_tmp].re - (reAij *
                    Z[ad22_re_tmp].re + scale * Z[ad22_re_tmp].im);
                  Z[ctemp_tmp].im = anorm * Z[ctemp_tmp].im - (reAij * sumsq -
                    scale * imAij);
                  Z[ad22_re_tmp].re = ad22_re;
                  Z[ad22_re_tmp].im = ad22_im;
                }

                jp1 = j - 1;
                j++;
              }
            }

            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (firstNonZero) {
      *info = ilast + 1;
      if (0 <= ilast) {
        std::memset(&alpha1[0], 0, (ilast + 1) * sizeof(creal_T));
        std::memset(&beta1[0], 0, (ilast + 1) * sizeof(creal_T));
      }

      std::memset(&Z[0], 0, 100U * sizeof(creal_T));
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1[j] = A[j + 10 * j];
    }
  }
}

//
// Arguments    : const creal32_T f
//                const creal32_T g
//                float *cs
//                creal32_T *sn
//                creal32_T *r
// Return Type  : void
//
static void c_xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
                      *sn, creal32_T *r)
{
  float scale_tmp;
  float f2s;
  float scale;
  float fs_re;
  float fs_im;
  float gs_re;
  float gs_im;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  float f2;
  float g2;
  scale_tmp = std::abs(f.re);
  f2s = std::abs(f.im);
  if (f2s > scale_tmp) {
    scale_tmp = f2s;
  }

  f2s = std::abs(g.re);
  scale = std::abs(g.im);
  if (scale > f2s) {
    f2s = scale;
  }

  scale = scale_tmp;
  if (f2s > scale_tmp) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = -1;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 5.49755814E+11F) {
    do {
      count++;
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0F > g2) {
      scale = 1.0F;
    }

    if (f2 <= scale * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        r->re = rt_hypotf(g.re, g.im);
        r->im = 0.0F;
        f2 = rt_hypotf(gs_re, gs_im);
        sn->re = gs_re / f2;
        sn->im = -gs_im / f2;
      } else {
        g2 = std::sqrt(g2);
        *cs = rt_hypotf(fs_re, fs_im) / g2;
        if (scale_tmp > 1.0F) {
          f2 = rt_hypotf(f.re, f.im);
          fs_re = f.re / f2;
          fs_im = f.im / f2;
        } else {
          scale = 5.49755814E+11F * f.re;
          f2s = 5.49755814E+11F * f.im;
          f2 = rt_hypotf(scale, f2s);
          fs_re = scale / f2;
          fs_im = f2s / f2;
        }

        gs_re /= g2;
        gs_im = -gs_im / g2;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = std::sqrt(g2 / f2 + 1.0F);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0F / f2s;
      f2 += g2;
      f2s = r->re / f2;
      scale = r->im / f2;
      sn->re = f2s * gs_re - scale * -gs_im;
      sn->im = f2s * -gs_im + scale * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 5.49755814E+11F;
          r->im *= 5.49755814E+11F;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 0; rescaledir <= count; rescaledir++) {
            r->re *= 1.8189894E-12F;
            r->im *= 1.8189894E-12F;
          }
        }
      }
    }
  }
}

//
// Arguments    : const double varargin_1_data[]
//                const int varargin_1_size[3]
//                const double varargin_2[9]
//                double y_data[]
//                int y_size[3]
// Return Type  : void
//
static void cat(const double varargin_1_data[], const int varargin_1_size[3],
                const double varargin_2[9], double y_data[], int y_size[3])
{
  int iy;
  int i;
  int j;
  y_size[0] = 3;
  y_size[1] = 3;
  y_size[2] = static_cast<signed char>((varargin_1_size[2] + 1));
  iy = -1;
  i = 9 * varargin_1_size[2];
  for (j = 0; j < i; j++) {
    iy++;
    y_data[iy] = varargin_1_data[j];
  }

  for (j = 0; j < 9; j++) {
    iy++;
    y_data[iy] = varargin_2[j];
  }
}

//
// Arguments    : const double A[16]
//                double B[16]
// Return Type  : void
//
static void d_mldivide(const double A[16], double B[16])
{
  double b_A[16];
  signed char ipiv[4];
  int j;
  int mmj_tmp;
  int b;
  int iy;
  int jj;
  int jp1j;
  int jA;
  double smax;
  int i;
  int ix;
  int k;
  double s;
  std::memcpy(&b_A[0], &A[0], 16U * sizeof(double));
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  ipiv[3] = 4;
  for (j = 0; j < 3; j++) {
    mmj_tmp = 2 - j;
    b = j * 5;
    jj = j * 5;
    jp1j = b + 2;
    iy = 4 - j;
    jA = 0;
    ix = b;
    smax = std::abs(b_A[jj]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = std::abs(b_A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (b_A[jj + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = static_cast<signed char>((iy + 1));
        smax = b_A[j];
        b_A[j] = b_A[iy];
        b_A[iy] = smax;
        ix = j + 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
      }

      i = (jj - j) + 4;
      for (jA = jp1j; jA <= i; jA++) {
        b_A[jA - 1] /= b_A[jj];
      }
    }

    iy = b + 4;
    jA = jj;
    for (k = 0; k <= mmj_tmp; k++) {
      smax = b_A[iy];
      if (b_A[iy] != 0.0) {
        ix = jj + 1;
        i = jA + 6;
        b = (jA - j) + 8;
        for (jp1j = i; jp1j <= b; jp1j++) {
          b_A[jp1j - 1] += b_A[ix] * -smax;
          ix++;
        }
      }

      iy += 4;
      jA += 4;
    }

    if (ipiv[j] != j + 1) {
      smax = B[j];
      i = ipiv[j] - 1;
      B[j] = B[i];
      B[i] = smax;
      smax = B[j + 4];
      i = ipiv[j] + 3;
      B[j + 4] = B[i];
      B[i] = smax;
      smax = B[j + 8];
      i = ipiv[j] + 7;
      B[j + 8] = B[i];
      B[i] = smax;
      smax = B[j + 12];
      i = ipiv[j] + 11;
      B[j + 12] = B[i];
      B[i] = smax;
    }
  }

  for (j = 0; j < 4; j++) {
    iy = j << 2;
    if (B[iy] != 0.0) {
      for (jA = 2; jA < 5; jA++) {
        i = (jA + iy) - 1;
        B[i] -= B[iy] * b_A[jA - 1];
      }
    }

    if (B[iy + 1] != 0.0) {
      for (jA = 3; jA < 5; jA++) {
        i = (jA + iy) - 1;
        B[i] -= B[iy + 1] * b_A[jA + 3];
      }
    }

    if (B[iy + 2] != 0.0) {
      i = iy + 3;
      for (jA = 4; jA < 5; jA++) {
        B[i] -= B[iy + 2] * b_A[11];
      }
    }
  }

  for (j = 0; j < 4; j++) {
    iy = j << 2;
    smax = B[iy + 3];
    if (smax != 0.0) {
      B[iy + 3] = smax / b_A[15];
      for (jA = 0; jA < 3; jA++) {
        i = jA + iy;
        B[i] -= B[iy + 3] * b_A[jA + 12];
      }
    }

    smax = B[iy + 2];
    if (smax != 0.0) {
      B[iy + 2] = smax / b_A[10];
      for (jA = 0; jA < 2; jA++) {
        i = jA + iy;
        B[i] -= B[iy + 2] * b_A[jA + 8];
      }
    }

    smax = B[iy + 1];
    if (smax != 0.0) {
      B[iy + 1] = smax / b_A[5];
      for (jA = 0; jA < 1; jA++) {
        B[iy] -= B[iy + 1] * b_A[4];
      }
    }

    if (B[iy] != 0.0) {
      B[iy] /= b_A[0];
    }
  }
}

//
// Arguments    : const double A[9]
//                double Q[9]
//                double R[9]
// Return Type  : void
//
static void d_qr(const double A[9], double Q[9], double R[9])
{
  double b_A[9];
  double tau[3];
  double work[3];
  int i;
  int itau;
  int ii;
  int ix0;
  double atmp;
  int b_i;
  int knt;
  double c;
  double beta1;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int k;
  int ia;
  int ix;
  int exitg1;
  int i1;
  std::memcpy(&b_A[0], &A[0], 9U * sizeof(double));
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  tau[2] = 0.0;
  work[2] = 0.0;
  for (i = 0; i < 3; i++) {
    ii = i * 3 + i;
    if (i + 1 < 3) {
      atmp = b_A[ii];
      ix0 = ii + 2;
      tau[i] = 0.0;
      c = i_xnrm2(2 - i, b_A, ii + 2);
      if (c != 0.0) {
        beta1 = rt_hypotd(b_A[ii], c);
        if (b_A[ii] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          b_i = (ii - i) + 3;
          do {
            knt++;
            for (k = ix0; k <= b_i; k++) {
              b_A[k - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd(atmp, i_xnrm2(2 - i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          c = 1.0 / (atmp - beta1);
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - b_A[ii]) / beta1;
          c = 1.0 / (b_A[ii] - beta1);
          b_i = (ii - i) + 3;
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = 1.0;
      if (tau[i] != 0.0) {
        lastv = 3 - i;
        knt = (ii - i) + 2;
        while ((lastv > 0) && (b_A[knt] == 0.0)) {
          lastv--;
          knt--;
        }

        lastc = 2 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = (ii + (lastc - 1) * 3) + 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia + 1 <= knt + lastv) {
              if (b_A[ia] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        knt = ii + 4;
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(double));
          }

          ix0 = 0;
          b_i = (ii + 3 * (lastc - 1)) + 4;
          for (k = knt; k <= b_i; k += 3) {
            ix = ii;
            c = 0.0;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        d_xgerc(lastv, lastc, -tau[i], ii + 1, work, b_A, ii + 4);
      }

      b_A[ii] = atmp;
    } else {
      tau[2] = 0.0;
    }
  }

  itau = 2;
  for (ix0 = 0; ix0 < 3; ix0++) {
    for (i = 0; i <= ix0; i++) {
      knt = i + 3 * ix0;
      R[knt] = b_A[knt];
    }

    b_i = ix0 + 2;
    if (b_i <= 3) {
      std::memset(&R[(ix0 * 3 + b_i) + -1], 0, (4 - b_i) * sizeof(double));
    }

    work[ix0] = 0.0;
  }

  for (i = 2; i >= 0; i--) {
    ii = (i + i * 3) + 4;
    if (i + 1 < 3) {
      b_A[ii - 4] = 1.0;
      if (tau[itau] != 0.0) {
        lastv = 3 - i;
        knt = ii - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0)) {
          lastv--;
          knt--;
        }

        lastc = 2 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = ii + (lastc - 1) * 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= (knt + lastv) - 1) {
              if (b_A[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(double));
          }

          ix0 = 0;
          b_i = ii + 3 * (lastc - 1);
          for (k = ii; k <= b_i; k += 3) {
            ix = ii;
            c = 0.0;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix - 4];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        d_xgerc(lastv, lastc, -tau[itau], ii - 3, work, b_A, ii);
      }

      ix0 = ii - 2;
      b_i = (ii - i) - 1;
      for (k = ix0; k <= b_i; k++) {
        b_A[k - 1] *= -tau[itau];
      }
    }

    b_A[ii - 4] = 1.0 - tau[itau];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(ii - ix0) - 5] = 0.0;
    }

    itau--;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    Q[3 * ix0] = b_A[3 * ix0];
    knt = 3 * ix0 + 1;
    Q[knt] = b_A[knt];
    knt = 3 * ix0 + 2;
    Q[knt] = b_A[knt];
  }
}

//
// Arguments    : const float A[9]
// Return Type  : float
//
static float d_rcond(const float A[9])
{
  float result;
  float normA;
  int j;
  float s;
  int i;
  float b_A[9];
  int ipiv[3];
  int jA;
  int exitg1;
  float ainvnm;
  int iter;
  int kase;
  int jump;
  float x[3];
  int exitg2;
  int exitg3;
  int exitg4;
  int exitg5;
  int b_j;
  int b_i;
  int ix;
  float absrexk;
  float ainvnm_tmp;
  result = 0.0F;
  normA = 0.0F;
  for (j = 0; j < 3; j++) {
    s = (std::abs(A[3 * j]) + std::abs(A[3 * j + 1])) + std::abs(A[3 * j + 2]);
    if (s > normA) {
      normA = s;
    }
  }

  if (normA != 0.0F) {
    for (i = 0; i < 9; i++) {
      b_A[i] = A[i];
    }

    d_xzgetrf(b_A, ipiv, &jA);
    jA = 2;
    do {
      exitg1 = 0;
      if (jA + 1 > 0) {
        if (b_A[jA + 3 * jA] == 0.0F) {
          exitg1 = 1;
        } else {
          jA--;
        }
      } else {
        ainvnm = 0.0F;
        iter = 2;
        kase = 1;
        jump = 1;
        j = 0;
        x[0] = 0.333333343F;
        x[1] = 0.333333343F;
        x[2] = 0.333333343F;
        do {
          exitg2 = 0;
          do {
            exitg3 = 0;
            do {
              exitg4 = 0;
              do {
                exitg5 = 0;
                if (kase == 1) {
                  for (b_j = 0; b_j < 3; b_j++) {
                    jA = b_j + b_j * 3;
                    i = 1 - b_j;
                    for (b_i = 0; b_i <= i; b_i++) {
                      ix = (b_j + b_i) + 1;
                      x[ix] -= x[b_j] * b_A[(jA + b_i) + 1];
                    }
                  }

                  for (b_j = 2; b_j >= 0; b_j--) {
                    jA = b_j + b_j * 3;
                    x[b_j] /= b_A[jA];
                    for (b_i = 0; b_i < b_j; b_i++) {
                      ix = (b_j - b_i) - 1;
                      x[ix] -= x[b_j] * b_A[(jA - b_i) - 1];
                    }
                  }
                } else {
                  for (b_j = 0; b_j < 3; b_j++) {
                    jA = b_j * 3;
                    s = x[b_j];
                    for (b_i = 0; b_i < b_j; b_i++) {
                      s -= b_A[jA + b_i] * x[b_i];
                    }

                    x[b_j] = s / b_A[jA + b_j];
                  }

                  for (b_j = 2; b_j >= 0; b_j--) {
                    jA = b_j * 3;
                    s = x[b_j];
                    i = b_j + 2;
                    for (b_i = 3; b_i >= i; b_i--) {
                      s -= b_A[(jA + b_i) - 1] * x[b_i - 1];
                    }

                    x[b_j] = s;
                  }
                }

                if (jump == 1) {
                  s = std::abs(x[0]);
                  if (s > 1.17549435E-38F) {
                    x[0] /= s;
                  } else {
                    x[0] = 1.0F;
                  }

                  ainvnm = s + std::abs(x[1]);
                  s = std::abs(x[1]);
                  if (s > 1.17549435E-38F) {
                    x[1] /= s;
                  } else {
                    x[1] = 1.0F;
                  }

                  ainvnm += std::abs(x[2]);
                  s = std::abs(x[2]);
                  if (s > 1.17549435E-38F) {
                    x[2] /= s;
                  } else {
                    x[2] = 1.0F;
                  }

                  kase = 2;
                  jump = 2;
                } else {
                  exitg5 = 1;
                }
              } while (exitg5 == 0);

              if (jump == 2) {
                j = 0;
                s = std::abs(x[0]);
                absrexk = std::abs(x[1]);
                if (absrexk > s) {
                  j = 1;
                  s = absrexk;
                }

                if (std::abs(x[2]) > s) {
                  j = 2;
                }

                iter = 2;
                x[0] = 0.0F;
                x[1] = 0.0F;
                x[2] = 0.0F;
                x[j] = 1.0F;
                kase = 1;
                jump = 3;
              } else {
                exitg4 = 1;
              }
            } while (exitg4 == 0);

            if (jump == 3) {
              s = std::abs(x[0]);
              absrexk = std::abs(x[1]);
              ainvnm_tmp = std::abs(x[2]);
              ainvnm = (s + absrexk) + ainvnm_tmp;
              if (ainvnm <= x[0]) {
                x[0] = 1.0F;
                x[1] = -1.5F;
                x[2] = 2.0F;
                kase = 1;
                jump = 5;
              } else {
                if (s > 1.17549435E-38F) {
                  x[0] /= s;
                } else {
                  x[0] = 1.0F;
                }

                if (absrexk > 1.17549435E-38F) {
                  x[1] /= absrexk;
                } else {
                  x[1] = 1.0F;
                }

                if (ainvnm_tmp > 1.17549435E-38F) {
                  x[2] /= ainvnm_tmp;
                } else {
                  x[2] = 1.0F;
                }

                kase = 2;
                jump = 4;
              }
            } else {
              exitg3 = 1;
            }
          } while (exitg3 == 0);

          if (jump == 4) {
            jA = j;
            j = 0;
            s = std::abs(x[0]);
            absrexk = std::abs(x[1]);
            if (absrexk > s) {
              j = 1;
              s = absrexk;
            }

            if (std::abs(x[2]) > s) {
              j = 2;
            }

            if ((std::abs(x[jA]) != std::abs(x[j])) && (iter <= 5)) {
              iter++;
              x[0] = 0.0F;
              x[1] = 0.0F;
              x[2] = 0.0F;
              x[j] = 1.0F;
              kase = 1;
              jump = 3;
            } else {
              x[0] = 1.0F;
              x[1] = -1.5F;
              x[2] = 2.0F;
              kase = 1;
              jump = 5;
            }
          } else {
            exitg2 = 1;
          }
        } while ((exitg2 == 0) || (!(jump == 5)));

        s = 2.0F * ((std::abs(x[0]) + std::abs(x[1])) + std::abs(x[2])) / 3.0F /
          3.0F;
        if (s > ainvnm) {
          ainvnm = s;
        }

        if (ainvnm != 0.0F) {
          result = 1.0F / ainvnm / normA;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return result;
}

//
// Arguments    : int m
//                int n
//                double alpha1
//                int ix0
//                const double y[3]
//                double A[9]
//                int ia0
// Return Type  : void
//
static void d_xgerc(int m, int n, double alpha1, int ix0, const double y[3],
                    double A[9], int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i;
  int ijA;
  if (alpha1 != 0.0) {
    jA = ia0;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

//
// Arguments    : int n
//                const float x[100]
//                int ix0
// Return Type  : float
//
static float d_xnrm2(int n, const float x[100], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 1.29246971E-26F;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : float x[100]
//                int ix0
//                int iy0
//                float c
//                float s
// Return Type  : void
//
static void d_xrot(float x[100], int ix0, int iy0, float c, float s)
{
  int ix;
  int iy;
  int k;
  float temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 0; k < 10; k++) {
    temp = c * x[ix] + s * x[iy];
    x[iy] = c * x[iy] - s * x[ix];
    x[ix] = temp;
    iy++;
    ix++;
  }
}

//
// Arguments    : float A[9]
//                int ipiv[3]
//                int *info
// Return Type  : void
//
static void d_xzgetrf(float A[9], int ipiv[3], int *info)
{
  int j;
  int mmj_tmp;
  int b_tmp;
  int jp1j;
  int iy;
  int jA;
  int ix;
  float smax;
  int k;
  float s;
  int i;
  int ijA;
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  *info = 0;
  for (j = 0; j < 2; j++) {
    mmj_tmp = 1 - j;
    b_tmp = j << 2;
    jp1j = b_tmp + 2;
    iy = 3 - j;
    jA = 0;
    ix = b_tmp;
    smax = std::abs(A[b_tmp]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (A[b_tmp + jA] != 0.0F) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = iy + 1;
        smax = A[j];
        A[j] = A[iy];
        A[iy] = smax;
        ix = j + 3;
        iy += 3;
        smax = A[ix];
        A[ix] = A[iy];
        A[iy] = smax;
        ix += 3;
        iy += 3;
        smax = A[ix];
        A[ix] = A[iy];
        A[iy] = smax;
      }

      i = (b_tmp - j) + 3;
      for (iy = jp1j; iy <= i; iy++) {
        A[iy - 1] /= A[b_tmp];
      }
    } else {
      *info = j + 1;
    }

    iy = b_tmp + 3;
    jA = b_tmp;
    for (k = 0; k <= mmj_tmp; k++) {
      smax = A[iy];
      if (A[iy] != 0.0F) {
        ix = b_tmp + 1;
        i = jA + 5;
        jp1j = (jA - j) + 6;
        for (ijA = i; ijA <= jp1j; ijA++) {
          A[ijA - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 3;
      jA += 3;
    }
  }

  if ((*info == 0) && (A[8] == 0.0F)) {
    *info = 3;
  }
}

//
// Arguments    : creal32_T A[100]
//                int ilo
//                int ihi
//                creal32_T Z[100]
//                int *info
//                creal32_T alpha1[10]
//                creal32_T beta1[10]
// Return Type  : void
//
static void d_xzhgeqz(creal32_T A[100], int ilo, int ihi, creal32_T Z[100], int *
                      info, creal32_T alpha1[10], creal32_T beta1[10])
{
  int i;
  float eshift_re;
  float eshift_im;
  creal32_T ctemp;
  float anorm;
  float scale;
  float reAij;
  float sumsq;
  float b_atol;
  boolean_T firstNonZero;
  int j;
  int ctemp_tmp;
  float ascale;
  int jp1;
  float imAij;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
  float temp2;
  int ilast;
  int ilastm1;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T b_guard1 = false;
  boolean_T guard3 = false;
  boolean_T exitg2;
  creal32_T shift;
  float ad22_re;
  float ad22_im;
  float t1_re;
  float t1_im;
  float t1_im_tmp;
  creal32_T b_ascale;
  int ad22_re_tmp;
  *info = 0;
  for (i = 0; i < 10; i++) {
    alpha1[i].re = 0.0F;
    alpha1[i].im = 0.0F;
    beta1[i].re = 1.0F;
    beta1[i].im = 0.0F;
  }

  eshift_re = 0.0F;
  eshift_im = 0.0F;
  ctemp.re = 0.0F;
  ctemp.im = 0.0F;
  anorm = 0.0F;
  if (ilo <= ihi) {
    scale = 0.0F;
    sumsq = 0.0F;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ctemp_tmp = j + 1;
      if (ihi < j + 1) {
        ctemp_tmp = ihi;
      }

      for (i = ilo; i <= ctemp_tmp; i++) {
        jp1 = (i + 10 * (j - 1)) - 1;
        reAij = A[jp1].re;
        imAij = A[jp1].im;
        if (reAij != 0.0F) {
          anorm = std::abs(reAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = sumsq * temp2 * temp2 + 1.0F;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0F) {
          anorm = std::abs(imAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = sumsq * temp2 * temp2 + 1.0F;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * std::sqrt(sumsq);
  }

  reAij = 1.1920929E-7F * anorm;
  b_atol = 1.17549435E-38F;
  if (reAij > 1.17549435E-38F) {
    b_atol = reAij;
  }

  reAij = 1.17549435E-38F;
  if (anorm > 1.17549435E-38F) {
    reAij = anorm;
  }

  ascale = 1.0F / reAij;
  firstNonZero = true;
  ctemp_tmp = ihi + 1;
  for (j = ctemp_tmp; j < 11; j++) {
    alpha1[j - 1] = A[(j + 10 * (j - 1)) - 1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          ctemp_tmp = ilast + 10 * ilastm1;
          if (std::abs(A[ctemp_tmp].re) + std::abs(A[ctemp_tmp].im) <= b_atol) {
            A[ctemp_tmp].re = 0.0F;
            A[ctemp_tmp].im = 0.0F;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                ctemp_tmp = j + 10 * (j - 1);
                if (std::abs(A[ctemp_tmp].re) + std::abs(A[ctemp_tmp].im) <=
                    b_atol) {
                  A[ctemp_tmp].re = 0.0F;
                  A[ctemp_tmp].im = 0.0F;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }

            if (goto70) {
              b_guard1 = true;
            } else {
              std::memset(&alpha1[0], 0, 10U * sizeof(creal32_T));
              std::memset(&beta1[0], 0, 10U * sizeof(creal32_T));
              std::memset(&Z[0], 0, 100U * sizeof(creal32_T));
              *info = 1;
              exitg1 = 1;
            }
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1[ilast] = A[ilast + 10 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0F;
              eshift_im = 0.0F;
              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                jp1 = ilastm1 + 10 * ilastm1;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0F) {
                  shift.re = anorm / 0.316227764F;
                  shift.im = 0.0F;
                } else if (anorm == 0.0F) {
                  shift.re = 0.0F;
                  shift.im = reAij / 0.316227764F;
                } else {
                  shift.re = anorm / 0.316227764F;
                  shift.im = reAij / 0.316227764F;
                }

                jp1 = ilast + 10 * ilast;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0F) {
                  ad22_re = anorm / 0.316227764F;
                  ad22_im = 0.0F;
                } else if (anorm == 0.0F) {
                  ad22_re = 0.0F;
                  ad22_im = reAij / 0.316227764F;
                } else {
                  ad22_re = anorm / 0.316227764F;
                  ad22_im = reAij / 0.316227764F;
                }

                t1_re = 0.5F * (shift.re + ad22_re);
                t1_im = 0.5F * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                jp1 = ilastm1 + 10 * ilast;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0F) {
                  imAij = anorm / 0.316227764F;
                  temp2 = 0.0F;
                } else if (anorm == 0.0F) {
                  imAij = 0.0F;
                  temp2 = reAij / 0.316227764F;
                } else {
                  imAij = anorm / 0.316227764F;
                  temp2 = reAij / 0.316227764F;
                }

                jp1 = ilast + 10 * ilastm1;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0F) {
                  sumsq = anorm / 0.316227764F;
                  anorm = 0.0F;
                } else if (anorm == 0.0F) {
                  sumsq = 0.0F;
                  anorm = reAij / 0.316227764F;
                } else {
                  sumsq = anorm / 0.316227764F;
                  anorm = reAij / 0.316227764F;
                }

                reAij = shift.re * ad22_re - shift.im * ad22_im;
                scale = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (imAij * sumsq -
                  temp2 * anorm)) - reAij;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (imAij * anorm + temp2 *
                  sumsq)) - scale;
                c_sqrt(&shift);
                if ((t1_re - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                    0.0F) {
                  shift.re += t1_re;
                  shift.im += t1_im;
                } else {
                  shift.re = t1_re - shift.re;
                  shift.im = t1_im - shift.im;
                }
              } else {
                jp1 = ilast + 10 * ilastm1;
                anorm = ascale * A[jp1].re;
                reAij = ascale * A[jp1].im;
                if (reAij == 0.0F) {
                  imAij = anorm / 0.316227764F;
                  temp2 = 0.0F;
                } else if (anorm == 0.0F) {
                  imAij = 0.0F;
                  temp2 = reAij / 0.316227764F;
                } else {
                  imAij = anorm / 0.316227764F;
                  temp2 = reAij / 0.316227764F;
                }

                eshift_re += imAij;
                eshift_im += temp2;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp_tmp = j + 10 * j;
                ctemp.re = ascale * A[ctemp_tmp].re - shift.re * 0.316227764F;
                ctemp.im = ascale * A[ctemp_tmp].im - shift.im * 0.316227764F;
                anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
                jp1 += 10 * j;
                temp2 = ascale * (std::abs(A[jp1].re) + std::abs(A[jp1].im));
                reAij = anorm;
                if (temp2 > anorm) {
                  reAij = temp2;
                }

                if ((reAij < 1.0F) && (reAij != 0.0F)) {
                  anorm /= reAij;
                  temp2 /= reAij;
                }

                ctemp_tmp = j + 10 * (j - 1);
                if ((std::abs(A[ctemp_tmp].re) + std::abs(A[ctemp_tmp].im)) *
                    temp2 <= anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp_tmp = (ifirst + 10 * (ifirst - 1)) - 1;
                ctemp.re = ascale * A[ctemp_tmp].re - shift.re * 0.316227764F;
                ctemp.im = ascale * A[ctemp_tmp].im - shift.im * 0.316227764F;
              }

              goto90 = false;
              jp1 = istart + 10 * (istart - 1);
              b_ascale.re = ascale * A[jp1].re;
              b_ascale.im = ascale * A[jp1].im;
              d_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              jp1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  ctemp_tmp = j + 10 * jp1;
                  c_xzlartg(A[ctemp_tmp - 1], A[ctemp_tmp], &anorm, &shift, &A
                            [(j + 10 * jp1) - 1]);
                  A[ctemp_tmp].re = 0.0F;
                  A[ctemp_tmp].im = 0.0F;
                }

                for (jp1 = j; jp1 < 11; jp1++) {
                  ctemp_tmp = j + 10 * (jp1 - 1);
                  ad22_re_tmp = ctemp_tmp - 1;
                  ad22_re = anorm * A[ad22_re_tmp].re + (shift.re * A[ctemp_tmp]
                    .re - shift.im * A[ctemp_tmp].im);
                  ad22_im = anorm * A[ad22_re_tmp].im + (shift.re * A[ctemp_tmp]
                    .im + shift.im * A[ctemp_tmp].re);
                  reAij = A[ad22_re_tmp].im;
                  scale = A[ad22_re_tmp].re;
                  A[ctemp_tmp].re = anorm * A[ctemp_tmp].re - (shift.re *
                    A[ad22_re_tmp].re + shift.im * A[ad22_re_tmp].im);
                  A[ctemp_tmp].im = anorm * A[ctemp_tmp].im - (shift.re * reAij
                    - shift.im * scale);
                  A[ad22_re_tmp].re = ad22_re;
                  A[ad22_re_tmp].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                jp1 = j;
                if (ilast + 1 < j + 2) {
                  jp1 = ilast - 1;
                }

                for (i = 1; i <= jp1 + 2; i++) {
                  ctemp_tmp = (i + 10 * (j - 1)) - 1;
                  ad22_re_tmp = (i + 10 * j) - 1;
                  ad22_re = anorm * A[ad22_re_tmp].re + (shift.re * A[ctemp_tmp]
                    .re - shift.im * A[ctemp_tmp].im);
                  ad22_im = anorm * A[ad22_re_tmp].im + (shift.re * A[ctemp_tmp]
                    .im + shift.im * A[ctemp_tmp].re);
                  reAij = A[ad22_re_tmp].im;
                  scale = A[ad22_re_tmp].re;
                  A[ctemp_tmp].re = anorm * A[ctemp_tmp].re - (shift.re *
                    A[ad22_re_tmp].re + shift.im * A[ad22_re_tmp].im);
                  A[ctemp_tmp].im = anorm * A[ctemp_tmp].im - (shift.re * reAij
                    - shift.im * scale);
                  A[ad22_re_tmp].re = ad22_re;
                  A[ad22_re_tmp].im = ad22_im;
                }

                reAij = shift.re;
                scale = shift.im;
                for (i = 0; i < 10; i++) {
                  ctemp_tmp = i + 10 * (j - 1);
                  ad22_re_tmp = i + 10 * j;
                  ad22_re = anorm * Z[ad22_re_tmp].re + (reAij * Z[ctemp_tmp].re
                    - scale * Z[ctemp_tmp].im);
                  ad22_im = anorm * Z[ad22_re_tmp].im + (reAij * Z[ctemp_tmp].im
                    + scale * Z[ctemp_tmp].re);
                  sumsq = Z[ad22_re_tmp].im;
                  imAij = Z[ad22_re_tmp].re;
                  Z[ctemp_tmp].re = anorm * Z[ctemp_tmp].re - (reAij *
                    Z[ad22_re_tmp].re + scale * Z[ad22_re_tmp].im);
                  Z[ctemp_tmp].im = anorm * Z[ctemp_tmp].im - (reAij * sumsq -
                    scale * imAij);
                  Z[ad22_re_tmp].re = ad22_re;
                  Z[ad22_re_tmp].im = ad22_im;
                }

                jp1 = j - 1;
                j++;
              }
            }

            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (firstNonZero) {
      *info = ilast + 1;
      if (0 <= ilast) {
        std::memset(&alpha1[0], 0, (ilast + 1) * sizeof(creal32_T));
        std::memset(&beta1[0], 0, (ilast + 1) * sizeof(creal32_T));
      }

      std::memset(&Z[0], 0, 100U * sizeof(creal32_T));
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1[j] = A[j + 10 * j];
    }
  }
}

//
// Arguments    : const creal32_T f
//                const creal32_T g
//                float *cs
//                creal32_T *sn
// Return Type  : void
//
static void d_xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
                      *sn)
{
  float scale_tmp;
  float f2;
  float scale;
  float fs_re;
  float fs_im;
  float gs_re;
  float gs_im;
  boolean_T guard1 = false;
  float g2;
  float g2s;
  scale_tmp = std::abs(f.re);
  f2 = std::abs(f.im);
  if (f2 > scale_tmp) {
    scale_tmp = f2;
  }

  f2 = std::abs(g.re);
  scale = std::abs(g.im);
  if (scale > f2) {
    f2 = scale;
  }

  scale = scale_tmp;
  if (f2 > scale_tmp) {
    scale = f2;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 5.49755814E+11F) {
    do {
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    guard1 = true;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
    } else {
      do {
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0F > g2) {
      scale = 1.0F;
    }

    if (f2 <= scale * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        g2 = rt_hypotf(gs_re, gs_im);
        sn->re = gs_re / g2;
        sn->im = -gs_im / g2;
      } else {
        g2s = std::sqrt(g2);
        *cs = rt_hypotf(fs_re, fs_im) / g2s;
        if (scale_tmp > 1.0F) {
          g2 = rt_hypotf(f.re, f.im);
          fs_re = f.re / g2;
          fs_im = f.im / g2;
        } else {
          f2 = 5.49755814E+11F * f.re;
          scale = 5.49755814E+11F * f.im;
          g2 = rt_hypotf(f2, scale);
          fs_re = f2 / g2;
          fs_im = scale / g2;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      scale = std::sqrt(g2 / f2 + 1.0F);
      *cs = 1.0F / scale;
      g2 += f2;
      fs_re = scale * fs_re / g2;
      fs_im = scale * fs_im / g2;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

//
// Arguments    : const float A[32]
//                float Q[64]
//                float R[32]
// Return Type  : void
//
static void e_qr(const float A[32], float Q[64], float R[32])
{
  int coltop;
  int i;
  float tau[8];
  int b_i;
  float work[8];
  int ii;
  float atmp;
  float c;
  float beta1;
  int lastv;
  int lastc;
  int knt;
  int c_i;
  boolean_T exitg2;
  int ia;
  int exitg1;
  int ix;
  int i1;
  for (coltop = 0; coltop < 4; coltop++) {
    for (i = 0; i < 8; i++) {
      b_i = i + (coltop << 3);
      Q[b_i] = A[b_i];
      Q[i + ((coltop + 4) << 3)] = 0.0F;
    }
  }

  for (i = 0; i < 8; i++) {
    tau[i] = 0.0F;
    work[i] = 0.0F;
  }

  for (i = 0; i < 8; i++) {
    ii = (i << 3) + i;
    if (i + 1 < 8) {
      atmp = Q[ii];
      coltop = ii + 2;
      tau[i] = 0.0F;
      c = m_xnrm2(7 - i, Q, ii + 2);
      if (c != 0.0F) {
        beta1 = rt_hypotf(Q[ii], c);
        if (Q[ii] >= 0.0F) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 9.86076132E-32F) {
          knt = -1;
          c_i = (ii - i) + 8;
          do {
            knt++;
            for (b_i = coltop; b_i <= c_i; b_i++) {
              Q[b_i - 1] *= 1.01412048E+31F;
            }

            beta1 *= 1.01412048E+31F;
            atmp *= 1.01412048E+31F;
          } while (!(std::abs(beta1) >= 9.86076132E-32F));

          beta1 = rt_hypotf(atmp, m_xnrm2(7 - i, Q, ii + 2));
          if (atmp >= 0.0F) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          c = 1.0F / (atmp - beta1);
          for (b_i = coltop; b_i <= c_i; b_i++) {
            Q[b_i - 1] *= c;
          }

          for (b_i = 0; b_i <= knt; b_i++) {
            beta1 *= 9.86076132E-32F;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - Q[ii]) / beta1;
          c = 1.0F / (Q[ii] - beta1);
          c_i = (ii - i) + 8;
          for (b_i = coltop; b_i <= c_i; b_i++) {
            Q[b_i - 1] *= c;
          }

          atmp = beta1;
        }
      }

      Q[ii] = 1.0F;
      if (tau[i] != 0.0F) {
        lastv = 8 - i;
        b_i = (ii - i) + 7;
        while ((lastv > 0) && (Q[b_i] == 0.0F)) {
          lastv--;
          b_i--;
        }

        lastc = 7 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          coltop = (ii + ((lastc - 1) << 3)) + 8;
          ia = coltop;
          do {
            exitg1 = 0;
            if (ia + 1 <= coltop + lastv) {
              if (Q[ia] != 0.0F) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        b_i = ii + 9;
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(float));
          }

          knt = 0;
          c_i = (ii + ((lastc - 1) << 3)) + 9;
          for (coltop = b_i; coltop <= c_i; coltop += 8) {
            ix = ii;
            c = 0.0F;
            i1 = (coltop + lastv) - 1;
            for (ia = coltop; ia <= i1; ia++) {
              c += Q[ia - 1] * Q[ix];
              ix++;
            }

            work[knt] += c;
            knt++;
          }
        }

        e_xgerc(lastv, lastc, -tau[i], ii + 1, work, Q, ii + 9);
      }

      Q[ii] = atmp;
    } else {
      tau[7] = 0.0F;
    }
  }

  R[0] = Q[0];
  for (i = 2; i < 9; i++) {
    R[i - 1] = 0.0F;
  }

  for (i = 0; i < 2; i++) {
    R[i + 8] = Q[i + 8];
  }

  for (i = 3; i < 9; i++) {
    R[i + 7] = 0.0F;
  }

  for (i = 0; i < 3; i++) {
    R[i + 16] = Q[i + 16];
  }

  for (i = 4; i < 9; i++) {
    R[i + 15] = 0.0F;
  }

  for (i = 0; i < 4; i++) {
    R[i + 24] = Q[i + 24];
  }

  for (i = 5; i < 9; i++) {
    R[i + 23] = 0.0F;
  }

  for (coltop = 0; coltop < 4; coltop++) {
    ia = (coltop + 4) << 3;
    for (i = 0; i < 8; i++) {
      Q[ia + i] = 0.0F;
    }

    Q[(ia + coltop) + 4] = 1.0F;
  }

  for (i = 0; i < 8; i++) {
    work[i] = 0.0F;
  }

  Q[27] = 1.0F;
  if (tau[3] != 0.0F) {
    lastv = 5;
    i = 33;
    while ((lastv > 0) && (Q[i - 2] == 0.0F)) {
      lastv--;
      i--;
    }

    lastc = 4;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 36;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0F) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(float));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 36;
      for (coltop = 36; coltop <= c_i; coltop += 8) {
        ix = 36;
        c = 0.0F;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    e_xgerc(lastv, lastc, -tau[3], 28, work, Q, 36);
  }

  for (b_i = 29; b_i < 33; b_i++) {
    Q[b_i - 1] *= -tau[3];
  }

  Q[27] = 1.0F - tau[3];
  for (coltop = 0; coltop < 3; coltop++) {
    Q[26 - coltop] = 0.0F;
  }

  Q[18] = 1.0F;
  if (tau[2] != 0.0F) {
    lastv = 6;
    i = 25;
    while ((lastv > 0) && (Q[i - 2] == 0.0F)) {
      lastv--;
      i--;
    }

    lastc = 5;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 27;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0F) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(float));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 27;
      for (coltop = 27; coltop <= c_i; coltop += 8) {
        ix = 27;
        c = 0.0F;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    e_xgerc(lastv, lastc, -tau[2], 19, work, Q, 27);
  }

  for (b_i = 20; b_i < 25; b_i++) {
    Q[b_i - 1] *= -tau[2];
  }

  Q[18] = 1.0F - tau[2];
  for (coltop = 0; coltop < 2; coltop++) {
    Q[17 - coltop] = 0.0F;
  }

  Q[9] = 1.0F;
  if (tau[1] != 0.0F) {
    lastv = 7;
    i = 17;
    while ((lastv > 0) && (Q[i - 2] == 0.0F)) {
      lastv--;
      i--;
    }

    lastc = 6;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 18;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0F) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(float));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 18;
      for (coltop = 18; coltop <= c_i; coltop += 8) {
        ix = 18;
        c = 0.0F;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    e_xgerc(lastv, lastc, -tau[1], 10, work, Q, 18);
  }

  for (b_i = 11; b_i < 17; b_i++) {
    Q[b_i - 1] *= -tau[1];
  }

  Q[9] = 1.0F - tau[1];
  Q[8] = 0.0F;
  Q[0] = 1.0F;
  if (tau[0] != 0.0F) {
    lastv = 8;
    i = 9;
    while ((lastv > 0) && (Q[i - 2] == 0.0F)) {
      lastv--;
      i--;
    }

    lastc = 7;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      coltop = ((lastc - 1) << 3) + 9;
      ia = coltop;
      do {
        exitg1 = 0;
        if (ia <= (coltop + lastv) - 1) {
          if (Q[ia - 1] != 0.0F) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      if (0 <= lastc - 1) {
        std::memset(&work[0], 0, lastc * sizeof(float));
      }

      knt = 0;
      c_i = ((lastc - 1) << 3) + 9;
      for (coltop = 9; coltop <= c_i; coltop += 8) {
        ix = 9;
        c = 0.0F;
        i1 = (coltop + lastv) - 1;
        for (ia = coltop; ia <= i1; ia++) {
          c += Q[ia - 1] * Q[ix - 9];
          ix++;
        }

        work[knt] += c;
        knt++;
      }
    }

    e_xgerc(lastv, lastc, -tau[0], 1, work, Q, 9);
  }

  for (b_i = 2; b_i < 9; b_i++) {
    Q[b_i - 1] *= -tau[0];
  }

  Q[0] = 1.0F - tau[0];
}

//
// Arguments    : int m
//                int n
//                float alpha1
//                int ix0
//                const float y[8]
//                float A[64]
//                int ia0
// Return Type  : void
//
static void e_xgerc(int m, int n, float alpha1, int ix0, const float y[8], float
                    A[64], int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i;
  int ijA;
  if (alpha1 != 0.0F) {
    jA = ia0;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 8;
    }
  }
}

//
// Arguments    : int n
//                const float x[3]
// Return Type  : float
//
static float e_xnrm2(int n, const float x[3])
{
  float y;
  float scale;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[1]);
    } else {
      scale = 1.29246971E-26F;
      absxk = std::abs(x[1]);
      if (absxk > 1.29246971E-26F) {
        y = 1.0F;
        scale = absxk;
      } else {
        t = absxk / 1.29246971E-26F;
        y = t * t;
      }

      absxk = std::abs(x[2]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0F;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : const double A[100]
//                creal_T V[100]
//                creal_T D[100]
// Return Type  : void
//
static void eig(const double A[100], creal_T V[100], creal_T D[100])
{
  boolean_T p;
  int j;
  boolean_T exitg2;
  int info;
  double b_V[100];
  double b_D[100];
  int exitg1;
  creal_T At[100];
  creal_T alpha1[10];
  creal_T beta1[10];
  int coltop;
  double colnorm;
  double scale;
  double t;
  double absxk;
  p = true;
  j = 0;
  exitg2 = false;
  while ((!exitg2) && (j < 10)) {
    info = 0;
    do {
      exitg1 = 0;
      if (info <= j) {
        if (A[info + 10 * j] != A[j + 10 * info]) {
          p = false;
          exitg1 = 1;
        } else {
          info++;
        }
      } else {
        j++;
        exitg1 = 2;
      }
    } while (exitg1 == 0);

    if (exitg1 == 1) {
      exitg2 = true;
    }
  }

  if (p) {
    schur(A, b_V, b_D);
    for (info = 0; info < 100; info++) {
      V[info].re = b_V[info];
      V[info].im = 0.0;
    }

    for (j = 0; j < 9; j++) {
      b_D[(j + 10 * j) + 1] = 0.0;
      for (info = 0; info <= j; info++) {
        b_D[info + 10 * (j + 1)] = 0.0;
      }
    }

    for (info = 0; info < 100; info++) {
      D[info].re = b_D[info];
      D[info].im = 0.0;
    }
  } else {
    for (info = 0; info < 100; info++) {
      At[info].re = A[info];
      At[info].im = 0.0;
    }

    xzggev(At, &info, alpha1, beta1, V);
    for (coltop = 0; coltop <= 90; coltop += 10) {
      colnorm = 0.0;
      scale = 3.3121686421112381E-170;
      info = coltop + 10;
      for (j = coltop + 1; j <= info; j++) {
        absxk = std::abs(V[j - 1].re);
        if (absxk > scale) {
          t = scale / absxk;
          colnorm = colnorm * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          colnorm += t * t;
        }

        absxk = std::abs(V[j - 1].im);
        if (absxk > scale) {
          t = scale / absxk;
          colnorm = colnorm * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          colnorm += t * t;
        }
      }

      colnorm = scale * std::sqrt(colnorm);
      info = coltop + 10;
      for (j = coltop + 1; j <= info; j++) {
        absxk = V[j - 1].re;
        scale = V[j - 1].im;
        if (scale == 0.0) {
          absxk /= colnorm;
          scale = 0.0;
        } else if (absxk == 0.0) {
          absxk = 0.0;
          scale /= colnorm;
        } else {
          absxk /= colnorm;
          scale /= colnorm;
        }

        V[j - 1].re = absxk;
        V[j - 1].im = scale;
      }
    }

    std::memset(&D[0], 0, 100U * sizeof(creal_T));
    for (j = 0; j < 10; j++) {
      if (beta1[j].im == 0.0) {
        if (alpha1[j].im == 0.0) {
          info = j + 10 * j;
          D[info].re = alpha1[j].re / beta1[j].re;
          D[info].im = 0.0;
        } else if (alpha1[j].re == 0.0) {
          info = j + 10 * j;
          D[info].re = 0.0;
          D[info].im = alpha1[j].im / beta1[j].re;
        } else {
          info = j + 10 * j;
          D[info].re = alpha1[j].re / beta1[j].re;
          D[info].im = alpha1[j].im / beta1[j].re;
        }
      } else if (beta1[j].re == 0.0) {
        if (alpha1[j].re == 0.0) {
          info = j + 10 * j;
          D[info].re = alpha1[j].im / beta1[j].im;
          D[info].im = 0.0;
        } else if (alpha1[j].im == 0.0) {
          info = j + 10 * j;
          D[info].re = 0.0;
          D[info].im = -(alpha1[j].re / beta1[j].im);
        } else {
          info = j + 10 * j;
          D[info].re = alpha1[j].im / beta1[j].im;
          D[info].im = -(alpha1[j].re / beta1[j].im);
        }
      } else {
        t = std::abs(beta1[j].re);
        scale = std::abs(beta1[j].im);
        if (t > scale) {
          scale = beta1[j].im / beta1[j].re;
          absxk = beta1[j].re + scale * beta1[j].im;
          info = j + 10 * j;
          D[info].re = (alpha1[j].re + scale * alpha1[j].im) / absxk;
          D[info].im = (alpha1[j].im - scale * alpha1[j].re) / absxk;
        } else if (scale == t) {
          if (beta1[j].re > 0.0) {
            scale = 0.5;
          } else {
            scale = -0.5;
          }

          if (beta1[j].im > 0.0) {
            absxk = 0.5;
          } else {
            absxk = -0.5;
          }

          info = j + 10 * j;
          D[info].re = (alpha1[j].re * scale + alpha1[j].im * absxk) / t;
          D[info].im = (alpha1[j].im * scale - alpha1[j].re * absxk) / t;
        } else {
          scale = beta1[j].re / beta1[j].im;
          absxk = beta1[j].im + scale * beta1[j].re;
          info = j + 10 * j;
          D[info].re = (scale * alpha1[j].re + alpha1[j].im) / absxk;
          D[info].im = (scale * alpha1[j].im - alpha1[j].re) / absxk;
        }
      }
    }
  }
}

//
// Arguments    : double h[100]
//                double z[100]
// Return Type  : int
//
static int eml_dlahqr(double h[100], double z[100])
{
  int info;
  double v[3];
  int j;
  int i;
  int b_i;
  boolean_T exitg1;
  int L;
  boolean_T goto150;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  int ix;
  double s;
  int nr;
  double ba;
  int hoffset;
  double d;
  double tst;
  double bb;
  double aa;
  double ab;
  double h22;
  double rt1r;
  int knt;
  int m;
  int b_k;
  info = 0;
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
  for (j = 0; j < 7; j++) {
    i = j + 10 * j;
    h[i + 2] = 0.0;
    h[i + 3] = 0.0;
  }

  h[79] = 0.0;
  b_i = 9;
  exitg1 = false;
  while ((!exitg1) && (b_i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      k = b_i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        i = k + 10 * (k - 1);
        ba = std::abs(h[i]);
        if (ba <= 1.0020841800044864E-291) {
          exitg3 = true;
        } else {
          ix = k + 10 * k;
          bb = std::abs(h[ix]);
          hoffset = i - 1;
          tst = std::abs(h[hoffset]) + bb;
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = std::abs(h[(k + 10 * (k - 2)) - 1]);
            }

            if (k + 2 <= 10) {
              tst += std::abs(h[ix + 1]);
            }
          }

          if (ba <= 2.2204460492503131E-16 * tst) {
            tst = std::abs(h[ix - 1]);
            if (ba > tst) {
              ab = ba;
              ba = tst;
            } else {
              ab = tst;
            }

            tst = std::abs(h[hoffset] - h[ix]);
            if (bb > tst) {
              aa = bb;
              bb = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            tst = 2.2204460492503131E-16 * (bb * (aa / s));
            if (1.0020841800044864E-291 > tst) {
              tst = 1.0020841800044864E-291;
            }

            if (ba * (ab / s) <= tst) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + 10 * (k - 1)] = 0.0;
      }

      if (k + 1 >= b_i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          hoffset = k + 10 * k;
          s = std::abs(h[hoffset + 1]) + std::abs(h[(k + 10 * (k + 1)) + 2]);
          tst = 0.75 * s + h[hoffset];
          aa = -0.4375 * s;
          ab = s;
          h22 = tst;
        } else if (its == 20) {
          s = std::abs(h[b_i + 10 * (b_i - 1)]) + std::abs(h[(b_i + 10 * (b_i -
            2)) - 1]);
          tst = 0.75 * s + h[b_i + 10 * b_i];
          aa = -0.4375 * s;
          ab = s;
          h22 = tst;
        } else {
          ix = b_i + 10 * (b_i - 1);
          tst = h[ix - 1];
          ab = h[ix];
          aa = h[(b_i + 10 * b_i) - 1];
          h22 = h[b_i + 10 * b_i];
        }

        s = ((std::abs(tst) + std::abs(aa)) + std::abs(ab)) + std::abs(h22);
        if (s == 0.0) {
          rt1r = 0.0;
          tst = 0.0;
          ba = 0.0;
          aa = 0.0;
        } else {
          tst /= s;
          ab /= s;
          aa /= s;
          h22 /= s;
          bb = (tst + h22) / 2.0;
          tst = (tst - bb) * (h22 - bb) - aa * ab;
          aa = std::sqrt(std::abs(tst));
          if (tst >= 0.0) {
            rt1r = bb * s;
            ba = rt1r;
            tst = aa * s;
            aa = -tst;
          } else {
            rt1r = bb + aa;
            ba = bb - aa;
            if (std::abs(rt1r - h22) <= std::abs(ba - h22)) {
              rt1r *= s;
              ba = rt1r;
            } else {
              ba *= s;
              rt1r = ba;
            }

            tst = 0.0;
            aa = 0.0;
          }
        }

        m = b_i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          hoffset = m + 10 * (m - 1);
          ix = hoffset - 1;
          ab = h[ix] - ba;
          s = (std::abs(ab) + std::abs(aa)) + std::abs(h[hoffset]);
          bb = h[hoffset] / s;
          hoffset = m + 10 * m;
          v[0] = (bb * h[hoffset - 1] + (h[ix] - rt1r) * (ab / s)) - tst * (aa /
            s);
          v[1] = bb * (((h[ix] + h[hoffset]) - rt1r) - ba);
          v[2] = bb * h[hoffset + 1];
          s = (std::abs(v[0]) + std::abs(v[1])) + std::abs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if (m == k + 1) {
            exitg3 = true;
          } else {
            i = m + 10 * (m - 2);
            if (std::abs(h[i - 1]) * (std::abs(v[1]) + std::abs(v[2])) <=
                2.2204460492503131E-16 * std::abs(v[0]) * ((std::abs(h[i - 2]) +
                  std::abs(h[ix])) + std::abs(h[hoffset]))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (b_k = m; b_k <= b_i; b_k++) {
          nr = (b_i - b_k) + 2;
          if (3 < nr) {
            nr = 3;
          }

          if (b_k > m) {
            hoffset = (b_k + 10 * (b_k - 2)) - 1;
            for (j = 0; j < nr; j++) {
              v[j] = h[j + hoffset];
            }
          }

          aa = v[0];
          bb = 0.0;
          if (nr > 0) {
            tst = b_xnrm2(nr - 1, v);
            if (tst != 0.0) {
              ab = rt_hypotd(v[0], tst);
              if (v[0] >= 0.0) {
                ab = -ab;
              }

              if (std::abs(ab) < 1.0020841800044864E-292) {
                knt = -1;
                do {
                  knt++;
                  for (ix = 2; ix <= nr; ix++) {
                    v[ix - 1] *= 9.9792015476736E+291;
                  }

                  ab *= 9.9792015476736E+291;
                  aa *= 9.9792015476736E+291;
                } while (!(std::abs(ab) >= 1.0020841800044864E-292));

                ab = rt_hypotd(aa, b_xnrm2(nr - 1, v));
                if (aa >= 0.0) {
                  ab = -ab;
                }

                bb = (ab - aa) / ab;
                tst = 1.0 / (aa - ab);
                for (ix = 2; ix <= nr; ix++) {
                  v[ix - 1] *= tst;
                }

                for (ix = 0; ix <= knt; ix++) {
                  ab *= 1.0020841800044864E-292;
                }

                aa = ab;
              } else {
                bb = (ab - v[0]) / ab;
                tst = 1.0 / (v[0] - ab);
                for (ix = 2; ix <= nr; ix++) {
                  v[ix - 1] *= tst;
                }

                aa = ab;
              }
            }
          }

          v[0] = aa;
          if (b_k > m) {
            h[(b_k + 10 * (b_k - 2)) - 1] = aa;
            i = b_k + 10 * (b_k - 2);
            h[i] = 0.0;
            if (b_k < b_i) {
              h[i + 1] = 0.0;
            }
          } else {
            if (m > k + 1) {
              h[(b_k + 10 * (b_k - 2)) - 1] *= 1.0 - bb;
            }
          }

          s = v[1];
          tst = bb * v[1];
          if (nr == 3) {
            d = v[2];
            ab = bb * v[2];
            for (j = b_k; j < 11; j++) {
              ix = b_k + 10 * (j - 1);
              hoffset = ix - 1;
              knt = ix + 1;
              aa = (h[hoffset] + s * h[ix]) + d * h[knt];
              h[hoffset] -= aa * bb;
              h[ix] -= aa * tst;
              h[knt] -= aa * ab;
            }

            if (b_k + 3 < b_i + 1) {
              i = b_k + 2;
            } else {
              i = b_i;
            }

            for (j = 0; j <= i; j++) {
              ix = j + 10 * (b_k - 1);
              hoffset = j + 10 * b_k;
              knt = j + 10 * (b_k + 1);
              aa = (h[ix] + s * h[hoffset]) + d * h[knt];
              h[ix] -= aa * bb;
              h[hoffset] -= aa * tst;
              h[knt] -= aa * ab;
            }

            for (j = 0; j < 10; j++) {
              ix = j + 10 * (b_k - 1);
              hoffset = j + 10 * b_k;
              knt = j + 10 * (b_k + 1);
              aa = (z[ix] + s * z[hoffset]) + d * z[knt];
              z[ix] -= aa * bb;
              z[hoffset] -= aa * tst;
              z[knt] -= aa * ab;
            }
          } else {
            if (nr == 2) {
              for (j = b_k; j < 11; j++) {
                ix = b_k + 10 * (j - 1);
                hoffset = ix - 1;
                aa = h[hoffset] + s * h[ix];
                h[hoffset] -= aa * bb;
                h[ix] -= aa * tst;
              }

              for (j = 0; j <= b_i; j++) {
                ix = j + 10 * (b_k - 1);
                hoffset = j + 10 * b_k;
                aa = h[ix] + s * h[hoffset];
                h[ix] -= aa * bb;
                h[hoffset] -= aa * tst;
              }

              for (j = 0; j < 10; j++) {
                ix = j + 10 * (b_k - 1);
                hoffset = j + 10 * b_k;
                aa = z[ix] + s * z[hoffset];
                z[ix] -= aa * bb;
                z[hoffset] -= aa * tst;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = b_i + 1;
      exitg1 = true;
    } else {
      if ((L != b_i + 1) && (L == b_i)) {
        i = b_i + 10 * b_i;
        ix = i - 1;
        s = h[ix];
        nr = 10 * (b_i - 1);
        hoffset = b_i + nr;
        d = h[hoffset];
        tst = h[i];
        xdlanv2(&h[(b_i + 10 * (b_i - 1)) - 1], &s, &d, &tst, &aa, &ab, &bb, &ba,
                &h22, &rt1r);
        h[ix] = s;
        h[hoffset] = d;
        h[i] = tst;
        if (10 > b_i + 1) {
          hoffset = 8 - b_i;
          knt = b_i + (b_i + 1) * 10;
          ix = knt - 1;
          for (k = 0; k <= hoffset; k++) {
            tst = h22 * h[ix] + rt1r * h[knt];
            h[knt] = h22 * h[knt] - rt1r * h[ix];
            h[ix] = tst;
            knt += 10;
            ix += 10;
          }
        }

        ix = nr + 1;
        hoffset = b_i * 10 + 1;
        xrot(b_i - 1, h, ix, hoffset, h22, rt1r);
        b_xrot(z, ix, hoffset, h22, rt1r);
      }

      b_i = L - 2;
    }
  }

  return info;
}

//
// Arguments    : creal_T h_data[]
//                const int h_size[2]
// Return Type  : int
//
static int eml_zlahqr(creal_T h_data[], const int h_size[2])
{
  int info;
  int n;
  int ldh;
  int i;
  int j;
  int knt;
  int h_tmp;
  int b_i;
  double SMLNUM;
  boolean_T exitg1;
  double tst;
  double aa;
  double br;
  int L;
  boolean_T goto140;
  creal_T sc;
  int its;
  boolean_T exitg2;
  int k;
  double re;
  boolean_T exitg3;
  double im;
  double ba;
  int ix0;
  double bb;
  double t_re;
  double ab;
  creal_T v;
  boolean_T goto70;
  int m;
  double u_re;
  double u_im;
  double s;
  int b_k;
  creal_T b_v[2];
  n = h_size[0];
  ldh = h_size[0];
  info = 0;
  if (1 != h_size[0]) {
    i = h_size[0];
    for (j = 0; j <= i - 4; j++) {
      knt = j + h_size[0] * j;
      h_tmp = knt + 2;
      h_data[h_tmp].re = 0.0;
      h_data[h_tmp].im = 0.0;
      knt += 3;
      h_data[knt].re = 0.0;
      h_data[knt].im = 0.0;
    }

    if (1 <= n - 2) {
      i = (n + h_size[0] * (n - 3)) - 1;
      h_data[i].re = 0.0;
      h_data[i].im = 0.0;
    }

    for (b_i = 2; b_i <= n; b_i++) {
      i = (b_i + h_size[0] * (b_i - 2)) - 1;
      if (h_data[i].im != 0.0) {
        tst = h_data[i].re;
        aa = h_data[i].im;
        br = std::abs(h_data[i].re) + std::abs(h_data[i].im);
        if (aa == 0.0) {
          sc.re = tst / br;
          sc.im = 0.0;
        } else if (tst == 0.0) {
          sc.re = 0.0;
          sc.im = aa / br;
        } else {
          sc.re = tst / br;
          sc.im = aa / br;
        }

        br = rt_hypotd(sc.re, sc.im);
        if (-sc.im == 0.0) {
          re = sc.re / br;
          im = 0.0;
        } else if (sc.re == 0.0) {
          re = 0.0;
          im = -sc.im / br;
        } else {
          re = sc.re / br;
          im = -sc.im / br;
        }

        h_data[i].re = rt_hypotd(h_data[i].re, h_data[i].im);
        h_data[i].im = 0.0;
        h_tmp = (b_i - 1) * ldh;
        ix0 = b_i + h_tmp;
        i = ix0 + ldh * (n - b_i);
        for (k = ix0; ldh < 0 ? k >= i : k <= i; k += ldh) {
          aa = re * h_data[k - 1].im + im * h_data[k - 1].re;
          h_data[k - 1].re = re * h_data[k - 1].re - im * h_data[k - 1].im;
          h_data[k - 1].im = aa;
        }

        ix0 = h_tmp + 1;
        knt = b_i + 1;
        if (n < knt) {
          knt = n;
        }

        i = h_tmp + knt;
        for (k = ix0; k <= i; k++) {
          aa = re * h_data[k - 1].im + -im * h_data[k - 1].re;
          h_data[k - 1].re = re * h_data[k - 1].re - -im * h_data[k - 1].im;
          h_data[k - 1].im = aa;
        }
      }
    }

    SMLNUM = 2.2250738585072014E-308 * (static_cast<double>(n) /
      2.2204460492503131E-16);
    b_i = n - 1;
    exitg1 = false;
    while ((!exitg1) && (b_i + 1 >= 1)) {
      L = -1;
      goto140 = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its < 301)) {
        k = b_i;
        exitg3 = false;
        while ((!exitg3) && (k + 1 > L + 2)) {
          i = k + h_size[0] * (k - 1);
          aa = std::abs(h_data[i].re);
          ba = aa + std::abs(h_data[i].im);
          if (ba <= SMLNUM) {
            exitg3 = true;
          } else {
            ix0 = k + h_size[0] * k;
            bb = std::abs(h_data[ix0].re) + std::abs(h_data[ix0].im);
            knt = i - 1;
            tst = (std::abs(h_data[knt].re) + std::abs(h_data[knt].im)) + bb;
            if (tst == 0.0) {
              if (k - 1 >= 1) {
                tst = std::abs(h_data[(k + h_size[0] * (k - 2)) - 1].re);
              }

              if (k + 2 <= n) {
                tst += std::abs(h_data[ix0 + 1].re);
              }
            }

            if (aa <= 2.2204460492503131E-16 * tst) {
              h_tmp = ix0 - 1;
              tst = std::abs(h_data[h_tmp].re) + std::abs(h_data[h_tmp].im);
              if (ba > tst) {
                ab = ba;
                ba = tst;
              } else {
                ab = tst;
              }

              tst = std::abs(h_data[knt].re - h_data[ix0].re) + std::abs
                (h_data[knt].im - h_data[ix0].im);
              if (bb > tst) {
                aa = bb;
                bb = tst;
              } else {
                aa = tst;
              }

              s = aa + ab;
              aa = 2.2204460492503131E-16 * (bb * (aa / s));
              if (SMLNUM > aa) {
                aa = SMLNUM;
              }

              if (ba * (ab / s) <= aa) {
                exitg3 = true;
              } else {
                k--;
              }
            } else {
              k--;
            }
          }
        }

        L = k - 1;
        if (k + 1 > 1) {
          i = k + h_size[0] * (k - 1);
          h_data[i].re = 0.0;
          h_data[i].im = 0.0;
        }

        if (k + 1 >= b_i + 1) {
          goto140 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            h_tmp = k + h_size[0] * k;
            t_re = 0.75 * std::abs(h_data[(k + h_size[0] * k) + 1].re) +
              h_data[h_tmp].re;
            ab = h_data[h_tmp].im;
          } else if (its == 20) {
            h_tmp = b_i + h_size[0] * b_i;
            t_re = 0.75 * std::abs(h_data[b_i + h_size[0] * (b_i - 1)].re) +
              h_data[h_tmp].re;
            ab = h_data[h_tmp].im;
          } else {
            h_tmp = b_i + h_size[0] * b_i;
            t_re = h_data[h_tmp].re;
            ab = h_data[h_tmp].im;
            v = h_data[h_tmp - 1];
            b_sqrt(&v);
            ix0 = b_i + h_size[0] * (b_i - 1);
            sc = h_data[ix0];
            b_sqrt(&sc);
            u_re = v.re * sc.re - v.im * sc.im;
            u_im = v.re * sc.im + v.im * sc.re;
            s = std::abs(u_re) + std::abs(u_im);
            if (s != 0.0) {
              knt = ix0 - 1;
              ba = 0.5 * (h_data[knt].re - h_data[h_tmp].re);
              im = 0.5 * (h_data[knt].im - h_data[h_tmp].im);
              bb = std::abs(ba) + std::abs(im);
              if (s <= bb) {
                s = bb;
              }

              if (im == 0.0) {
                t_re = ba / s;
                ab = 0.0;
              } else if (ba == 0.0) {
                t_re = 0.0;
                ab = im / s;
              } else {
                t_re = ba / s;
                ab = im / s;
              }

              re = t_re * t_re - ab * ab;
              tst = t_re * ab;
              if (u_im == 0.0) {
                sc.re = u_re / s;
                sc.im = 0.0;
              } else if (u_re == 0.0) {
                sc.re = 0.0;
                sc.im = u_im / s;
              } else {
                sc.re = u_re / s;
                sc.im = u_im / s;
              }

              aa = sc.re * sc.re - sc.im * sc.im;
              ab = sc.re * sc.im;
              v.re = re + aa;
              v.im = (tst + tst) + (ab + ab);
              b_sqrt(&v);
              sc.re = s * v.re;
              sc.im = s * v.im;
              if (bb > 0.0) {
                if (im == 0.0) {
                  t_re = ba / bb;
                  ab = 0.0;
                } else if (ba == 0.0) {
                  t_re = 0.0;
                  ab = im / bb;
                } else {
                  t_re = ba / bb;
                  ab = im / bb;
                }

                if (t_re * sc.re + ab * sc.im < 0.0) {
                  sc.re = -sc.re;
                  sc.im = -sc.im;
                }
              }

              br = ba + sc.re;
              ab = im + sc.im;
              if (ab == 0.0) {
                if (u_im == 0.0) {
                  ba = u_re / br;
                  tst = 0.0;
                } else if (u_re == 0.0) {
                  ba = 0.0;
                  tst = u_im / br;
                } else {
                  ba = u_re / br;
                  tst = u_im / br;
                }
              } else if (br == 0.0) {
                if (u_re == 0.0) {
                  ba = u_im / ab;
                  tst = 0.0;
                } else if (u_im == 0.0) {
                  ba = 0.0;
                  tst = -(u_re / ab);
                } else {
                  ba = u_im / ab;
                  tst = -(u_re / ab);
                }
              } else {
                bb = std::abs(br);
                tst = std::abs(ab);
                if (bb > tst) {
                  s = ab / br;
                  tst = br + s * ab;
                  ba = (u_re + s * u_im) / tst;
                  tst = (u_im - s * u_re) / tst;
                } else if (tst == bb) {
                  if (br > 0.0) {
                    aa = 0.5;
                  } else {
                    aa = -0.5;
                  }

                  if (ab > 0.0) {
                    tst = 0.5;
                  } else {
                    tst = -0.5;
                  }

                  ba = (u_re * aa + u_im * tst) / bb;
                  tst = (u_im * aa - u_re * tst) / bb;
                } else {
                  s = br / ab;
                  tst = ab + s * br;
                  ba = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h_data[h_tmp].re - (u_re * ba - u_im * tst);
              ab = h_data[h_tmp].im - (u_re * tst + u_im * ba);
            }
          }

          goto70 = false;
          m = b_i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            h_tmp = m + h_size[0] * (m - 1);
            ix0 = h_tmp - 1;
            sc.re = h_data[ix0].re - t_re;
            sc.im = h_data[ix0].im - ab;
            tst = h_data[h_tmp].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            if (sc.im == 0.0) {
              re = sc.re / s;
              im = 0.0;
            } else if (sc.re == 0.0) {
              re = 0.0;
              im = sc.im / s;
            } else {
              re = sc.re / s;
              im = sc.im / s;
            }

            sc.re = re;
            sc.im = im;
            tst /= s;
            b_v[0] = sc;
            b_v[1].re = tst;
            b_v[1].im = 0.0;
            i = m + h_size[0] * m;
            if (std::abs(h_data[(m + h_size[0] * (m - 2)) - 1].re) * std::abs
                (tst) <= 2.2204460492503131E-16 * ((std::abs(re) + std::abs(im))
                 * ((std::abs(h_data[ix0].re) + std::abs(h_data[ix0].im)) + (std::
                   abs(h_data[i].re) + std::abs(h_data[i].im))))) {
              goto70 = true;
              exitg3 = true;
            } else {
              m--;
            }
          }

          if (!goto70) {
            ix0 = k + h_size[0] * k;
            sc.re = h_data[ix0].re - t_re;
            sc.im = h_data[ix0].im - ab;
            tst = h_data[(k + h_size[0] * k) + 1].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            if (sc.im == 0.0) {
              b_v[0].re = sc.re / s;
              b_v[0].im = 0.0;
            } else if (sc.re == 0.0) {
              b_v[0].re = 0.0;
              b_v[0].im = sc.im / s;
            } else {
              b_v[0].re = sc.re / s;
              b_v[0].im = sc.im / s;
            }

            tst /= s;
            b_v[1].re = tst;
            b_v[1].im = 0.0;
          }

          for (b_k = m; b_k <= b_i; b_k++) {
            if (b_k > m) {
              knt = b_k + h_size[0] * (b_k - 2);
              b_v[0] = h_data[knt - 1];
              b_v[1] = h_data[knt];
            }

            ba = b_v[1].re;
            im = b_v[1].im;
            sc = b_v[0];
            t_re = 0.0;
            ab = 0.0;
            tst = rt_hypotd(b_v[1].re, b_v[1].im);
            if ((tst != 0.0) || (b_v[0].im != 0.0)) {
              aa = xdlapy3(b_v[0].re, b_v[0].im, tst);
              if (b_v[0].re >= 0.0) {
                aa = -aa;
              }

              if (std::abs(aa) < 1.0020841800044864E-292) {
                knt = -1;
                do {
                  knt++;
                  ba *= 9.9792015476736E+291;
                  im *= 9.9792015476736E+291;
                  aa *= 9.9792015476736E+291;
                  sc.re *= 9.9792015476736E+291;
                  sc.im *= 9.9792015476736E+291;
                } while (!(std::abs(aa) >= 1.0020841800044864E-292));

                aa = xdlapy3(sc.re, sc.im, rt_hypotd(ba, im));
                if (sc.re >= 0.0) {
                  aa = -aa;
                }

                tst = aa - sc.re;
                if (0.0 - sc.im == 0.0) {
                  t_re = tst / aa;
                  ab = 0.0;
                } else if (tst == 0.0) {
                  t_re = 0.0;
                  ab = (0.0 - sc.im) / aa;
                } else {
                  t_re = tst / aa;
                  ab = (0.0 - sc.im) / aa;
                }

                v.re = sc.re - aa;
                v.im = sc.im;
                sc = recip(v);
                re = sc.re * ba - sc.im * im;
                im = sc.re * im + sc.im * ba;
                ba = re;
                for (h_tmp = 0; h_tmp <= knt; h_tmp++) {
                  aa *= 1.0020841800044864E-292;
                }

                sc.re = aa;
                sc.im = 0.0;
              } else {
                tst = aa - b_v[0].re;
                if (0.0 - b_v[0].im == 0.0) {
                  t_re = tst / aa;
                  ab = 0.0;
                } else if (tst == 0.0) {
                  t_re = 0.0;
                  ab = (0.0 - b_v[0].im) / aa;
                } else {
                  t_re = tst / aa;
                  ab = (0.0 - b_v[0].im) / aa;
                }

                v.re = b_v[0].re - aa;
                v.im = b_v[0].im;
                v = recip(v);
                ba = v.re * b_v[1].re - v.im * b_v[1].im;
                im = v.re * b_v[1].im + v.im * b_v[1].re;
                sc.re = aa;
                sc.im = 0.0;
              }
            }

            b_v[0] = sc;
            b_v[1].re = ba;
            b_v[1].im = im;
            if (b_k > m) {
              h_data[(b_k + h_size[0] * (b_k - 2)) - 1] = sc;
              i = b_k + h_size[0] * (b_k - 2);
              h_data[i].re = 0.0;
              h_data[i].im = 0.0;
            }

            tst = t_re * ba - ab * im;
            for (j = b_k; j <= n; j++) {
              h_tmp = b_k + h_size[0] * (j - 1);
              ix0 = h_tmp - 1;
              sc.re = (t_re * h_data[ix0].re - -ab * h_data[ix0].im) + tst *
                h_data[h_tmp].re;
              sc.im = (t_re * h_data[ix0].im + -ab * h_data[ix0].re) + tst *
                h_data[h_tmp].im;
              h_data[ix0].re -= sc.re;
              h_data[ix0].im -= sc.im;
              h_data[h_tmp].re -= sc.re * ba - sc.im * im;
              h_data[h_tmp].im -= sc.re * im + sc.im * ba;
            }

            if (b_k + 2 < b_i + 1) {
              i = b_k + 1;
            } else {
              i = b_i;
            }

            for (j = 0; j <= i; j++) {
              ix0 = j + h_size[0] * (b_k - 1);
              h_tmp = j + h_size[0] * b_k;
              sc.re = (t_re * h_data[ix0].re - ab * h_data[ix0].im) + tst *
                h_data[h_tmp].re;
              sc.im = (t_re * h_data[ix0].im + ab * h_data[ix0].re) + tst *
                h_data[h_tmp].im;
              h_data[ix0].re -= sc.re;
              h_data[ix0].im -= sc.im;
              h_data[h_tmp].re -= sc.re * ba - sc.im * -im;
              h_data[h_tmp].im -= sc.re * -im + sc.im * ba;
            }

            if ((b_k == m) && (m > k + 1)) {
              t_re = 1.0 - t_re;
              ab = 0.0 - ab;
              br = rt_hypotd(t_re, ab);
              if (ab == 0.0) {
                re = t_re / br;
                im = 0.0;
              } else if (t_re == 0.0) {
                re = 0.0;
                im = ab / br;
              } else {
                re = t_re / br;
                im = ab / br;
              }

              knt = m + h_size[0] * (m - 1);
              aa = h_data[knt].re * -im + h_data[knt].im * re;
              h_data[knt].re = h_data[knt].re * re - h_data[knt].im * -im;
              h_data[knt].im = aa;
              if (m + 2 <= b_i + 1) {
                knt = (m + h_size[0] * m) + 1;
                aa = h_data[knt].re * im + h_data[knt].im * re;
                h_data[knt].re = h_data[knt].re * re - h_data[knt].im * im;
                h_data[knt].im = aa;
              }

              for (j = m; j <= b_i + 1; j++) {
                if (j != m + 1) {
                  if (n > j) {
                    ix0 = j + j * ldh;
                    i = ix0 + ldh * ((n - j) - 1);
                    for (h_tmp = ix0; ldh < 0 ? h_tmp >= i : h_tmp <= i; h_tmp +=
                         ldh) {
                      aa = re * h_data[h_tmp - 1].im + im * h_data[h_tmp - 1].re;
                      h_data[h_tmp - 1].re = re * h_data[h_tmp - 1].re - im *
                        h_data[h_tmp - 1].im;
                      h_data[h_tmp - 1].im = aa;
                    }
                  }

                  h_tmp = (j - 1) * ldh;
                  ix0 = h_tmp + 1;
                  i = (h_tmp + j) - 1;
                  for (h_tmp = ix0; h_tmp <= i; h_tmp++) {
                    aa = re * h_data[h_tmp - 1].im + -im * h_data[h_tmp - 1].re;
                    h_data[h_tmp - 1].re = re * h_data[h_tmp - 1].re - -im *
                      h_data[h_tmp - 1].im;
                    h_data[h_tmp - 1].im = aa;
                  }
                }
              }
            }
          }

          h_tmp = b_i + h_size[0] * (b_i - 1);
          t_re = h_data[h_tmp].re;
          ab = h_data[h_tmp].im;
          if (h_data[h_tmp].im != 0.0) {
            aa = rt_hypotd(h_data[h_tmp].re, h_data[h_tmp].im);
            h_data[h_tmp].re = aa;
            h_data[h_tmp].im = 0.0;
            if (ab == 0.0) {
              re = t_re / aa;
              im = 0.0;
            } else if (t_re == 0.0) {
              re = 0.0;
              im = ab / aa;
            } else {
              re = t_re / aa;
              im = ab / aa;
            }

            if (n > b_i + 1) {
              ix0 = (b_i + (b_i + 1) * ldh) + 1;
              i = ix0 + ldh * ((n - b_i) - 2);
              for (k = ix0; ldh < 0 ? k >= i : k <= i; k += ldh) {
                aa = re * h_data[k - 1].im + -im * h_data[k - 1].re;
                h_data[k - 1].re = re * h_data[k - 1].re - -im * h_data[k - 1].
                  im;
                h_data[k - 1].im = aa;
              }
            }

            h_tmp = b_i * ldh;
            ix0 = h_tmp + 1;
            i = h_tmp + b_i;
            for (k = ix0; k <= i; k++) {
              aa = re * h_data[k - 1].im + im * h_data[k - 1].re;
              h_data[k - 1].re = re * h_data[k - 1].re - im * h_data[k - 1].im;
              h_data[k - 1].im = aa;
            }
          }

          its++;
        }
      }

      if (!goto140) {
        info = b_i + 1;
        exitg1 = true;
      } else {
        b_i = L;
      }
    }
  }

  return info;
}

//
// function G = equations_for_groebner(F)
// Arguments    : const double F[72]
//                double G[112]
// Return Type  : void
//
static void equations_for_groebner(const double F[72], double G[112])
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

  // 'equations_for_groebner:2' G = zeros(4, 1, 28, 'like', F);
  // 'equations_for_groebner:3' G(1, 1, :) = find_det3(F(2:4, :, :));
  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  // 'equations_for_groebner:4' G(2, 1, :) = find_det3([F(1, :, :); F(3:4, :, :)]); 
  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * i;
      b_M_tmp = 3 * i1 + 9 * i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 2];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  // 'equations_for_groebner:5' G(3, 1, :) = find_det3([F(1:2, :, :); F(4, :, :)]); 
  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * i;
      b_M_tmp = 3 * i1 + 9 * i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 1];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  // 'equations_for_groebner:6' G(4, 1, :) = find_det3(F(1:3, :, :));
  // 'equations_for_groebner:11' d = mult_poly42(find_det2(M(2:3, 2:3, :)), M(1, 1, :)) - ... 
  // 'equations_for_groebner:12'         mult_poly42(find_det2([M(2:3,1,:) M(2:3,3,:)]), M(1, 2, :)) + ... 
  // 'equations_for_groebner:13'         mult_poly42(find_det2(M(2:3, 1:2, :)), M(1, 3, :)); 
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
  // 'equations_for_groebner:17' d = mult_poly22(M(1, 1, :), M(2, 2, :)) - ...
  // 'equations_for_groebner:18'         mult_poly22(M(1, 2, :), M(2, 1, :));
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

  // 'equations_for_groebner:7' G = squeeze(G);
}

//
// Arguments    : const float A[9]
//                float Q[9]
//                float R[9]
// Return Type  : void
//
static void f_qr(const float A[9], float Q[9], float R[9])
{
  int i;
  float tau[3];
  float b_A[9];
  float work[3];
  int b_i;
  int itau;
  int ii;
  int ix0;
  float atmp;
  int knt;
  float c;
  float beta1;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int k;
  int ia;
  int ix;
  int exitg1;
  int i1;
  for (i = 0; i < 9; i++) {
    b_A[i] = A[i];
  }

  tau[0] = 0.0F;
  work[0] = 0.0F;
  tau[1] = 0.0F;
  work[1] = 0.0F;
  tau[2] = 0.0F;
  work[2] = 0.0F;
  for (b_i = 0; b_i < 3; b_i++) {
    ii = b_i * 3 + b_i;
    if (b_i + 1 < 3) {
      atmp = b_A[ii];
      ix0 = ii + 2;
      tau[b_i] = 0.0F;
      c = o_xnrm2(2 - b_i, b_A, ii + 2);
      if (c != 0.0F) {
        beta1 = rt_hypotf(b_A[ii], c);
        if (b_A[ii] >= 0.0F) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 9.86076132E-32F) {
          knt = -1;
          i = (ii - b_i) + 3;
          do {
            knt++;
            for (k = ix0; k <= i; k++) {
              b_A[k - 1] *= 1.01412048E+31F;
            }

            beta1 *= 1.01412048E+31F;
            atmp *= 1.01412048E+31F;
          } while (!(std::abs(beta1) >= 9.86076132E-32F));

          beta1 = rt_hypotf(atmp, o_xnrm2(2 - b_i, b_A, ii + 2));
          if (atmp >= 0.0F) {
            beta1 = -beta1;
          }

          tau[b_i] = (beta1 - atmp) / beta1;
          c = 1.0F / (atmp - beta1);
          for (k = ix0; k <= i; k++) {
            b_A[k - 1] *= c;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 9.86076132E-32F;
          }

          atmp = beta1;
        } else {
          tau[b_i] = (beta1 - b_A[ii]) / beta1;
          c = 1.0F / (b_A[ii] - beta1);
          i = (ii - b_i) + 3;
          for (k = ix0; k <= i; k++) {
            b_A[k - 1] *= c;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = 1.0F;
      if (tau[b_i] != 0.0F) {
        lastv = 3 - b_i;
        knt = (ii - b_i) + 2;
        while ((lastv > 0) && (b_A[knt] == 0.0F)) {
          lastv--;
          knt--;
        }

        lastc = 2 - b_i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = (ii + (lastc - 1) * 3) + 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia + 1 <= knt + lastv) {
              if (b_A[ia] != 0.0F) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        knt = ii + 4;
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(float));
          }

          ix0 = 0;
          i = (ii + 3 * (lastc - 1)) + 4;
          for (k = knt; k <= i; k += 3) {
            ix = ii;
            c = 0.0F;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        f_xgerc(lastv, lastc, -tau[b_i], ii + 1, work, b_A, ii + 4);
      }

      b_A[ii] = atmp;
    } else {
      tau[2] = 0.0F;
    }
  }

  itau = 2;
  for (ix0 = 0; ix0 < 3; ix0++) {
    for (b_i = 0; b_i <= ix0; b_i++) {
      knt = b_i + 3 * ix0;
      R[knt] = b_A[knt];
    }

    i = ix0 + 2;
    if (i <= 3) {
      std::memset(&R[(ix0 * 3 + i) + -1], 0, (4 - i) * sizeof(float));
    }

    work[ix0] = 0.0F;
  }

  for (b_i = 2; b_i >= 0; b_i--) {
    ii = (b_i + b_i * 3) + 4;
    if (b_i + 1 < 3) {
      b_A[ii - 4] = 1.0F;
      if (tau[itau] != 0.0F) {
        lastv = 3 - b_i;
        knt = ii - b_i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0F)) {
          lastv--;
          knt--;
        }

        lastc = 2 - b_i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = ii + (lastc - 1) * 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= (knt + lastv) - 1) {
              if (b_A[ia - 1] != 0.0F) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(float));
          }

          ix0 = 0;
          i = ii + 3 * (lastc - 1);
          for (k = ii; k <= i; k += 3) {
            ix = ii;
            c = 0.0F;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix - 4];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        f_xgerc(lastv, lastc, -tau[itau], ii - 3, work, b_A, ii);
      }

      ix0 = ii - 2;
      i = (ii - b_i) - 1;
      for (k = ix0; k <= i; k++) {
        b_A[k - 1] *= -tau[itau];
      }
    }

    b_A[ii - 4] = 1.0F - tau[itau];
    for (ix0 = 0; ix0 < b_i; ix0++) {
      b_A[(ii - ix0) - 5] = 0.0F;
    }

    itau--;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    Q[3 * ix0] = b_A[3 * ix0];
    knt = 3 * ix0 + 1;
    Q[knt] = b_A[knt];
    knt = 3 * ix0 + 2;
    Q[knt] = b_A[knt];
  }
}

//
// Arguments    : int m
//                int n
//                float alpha1
//                int ix0
//                const float y[3]
//                float A[9]
//                int ia0
// Return Type  : void
//
static void f_xgerc(int m, int n, float alpha1, int ix0, const float y[3], float
                    A[9], int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i;
  int ijA;
  if (alpha1 != 0.0F) {
    jA = ia0;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

//
// Arguments    : int n
//                const float x[12]
//                int ix0
// Return Type  : float
//
static float f_xnrm2(int n, const float x[12], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 1.29246971E-26F;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// function [fc, fs, n] = find_f(F, x, y, e)
// Arguments    : const double F[72]
//                double x
//                double y
//                double e
//                double fc_data[]
//                int fc_size[2]
//                double fs_data[]
//                int fs_size[2]
//                double *n
// Return Type  : void
//
static void find_f(const double F[72], double x, double y, double e, double
                   fc_data[], int fc_size[2], double fs_data[], int fs_size[2],
                   double *n)
{
  double F_eval[12];
  double mons[6];
  int i;
  double b_F_eval[12];
  double Q[9];
  int j;
  int b_i;
  double fc_tmp;
  int k;

  // 'find_f:2' F_eval = zeros(4, 3, 'like', F);
  std::memset(&F_eval[0], 0, 12U * sizeof(double));

  // 'find_f:3' mons =  [x^2, x*y, y^2, x, y, 1];
  mons[0] = x * x;
  mons[1] = x * y;
  mons[2] = y * y;
  mons[3] = x;
  mons[4] = y;
  mons[5] = 1.0;

  // 'find_f:4' for i = 1 : 4
  // 'find_f:12' [Q, R] = qr(F_eval');
  for (i = 0; i < 4; i++) {
    // 'find_f:5' for j = 1 : 3
    for (j = 0; j < 3; j++) {
      // 'find_f:6' for k = 1 : 6
      b_i = i + (j << 2);
      fc_tmp = F_eval[b_i];
      for (k = 0; k < 6; k++) {
        // 'find_f:7' F_eval(i, j) = F_eval(i, j) + F(i, j, k)*mons(k);
        fc_tmp += F[b_i + 12 * k] * mons[k];
      }

      F_eval[b_i] = fc_tmp;
      b_F_eval[j + 3 * i] = fc_tmp;
    }
  }

  qr(b_F_eval, Q, F_eval);

  // 'find_f:13' if abs(R(2, 2) - 0) < e
  if (std::abs(F_eval[4]) < e) {
    // rankF = 1
    // 'find_f:14' n = 2;
    *n = 2.0;

    // 'find_f:15' fc = [Q(1, 2) / Q(3, 2), Q(1, 2) / Q(3, 2)];
    fc_size[0] = 1;
    fc_size[1] = 2;
    fc_tmp = Q[3] / Q[5];
    fc_data[0] = fc_tmp;
    fc_data[1] = fc_tmp;

    // 'find_f:16' fs = [Q(2, 2) / Q(3, 2), Q(2, 3) / Q(3, 3)];
    fs_size[0] = 1;
    fs_size[1] = 2;
    fs_data[0] = Q[4] / Q[5];
    fs_data[1] = Q[7] / Q[8];
  } else {
    // 'find_f:17' else
    // 'find_f:18' n = 1;
    *n = 1.0;

    // 'find_f:19' fc = Q(1, 3) / Q(3, 3);
    fc_size[0] = 1;
    fc_size[1] = 1;
    fc_data[0] = Q[6] / Q[8];

    // 'find_f:20' fs = Q(2, 3) / Q(3, 3);
    fs_size[0] = 1;
    fs_size[1] = 1;
    fs_data[0] = Q[7] / Q[8];
  }
}

//
// Arguments    : int n
//                const double x[64]
//                int ix0
// Return Type  : double
//
static double g_xnrm2(int n, const double x[64], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                const creal_T x_data[]
//                int ix0
// Return Type  : double
//
static double h_xnrm2(int n, const creal_T x_data[], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = rt_hypotd(x_data[ix0 - 1].re, x_data[ix0 - 1].im);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x_data[k - 1].re);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }

        absxk = std::abs(x_data[k - 1].im);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                const double x[9]
//                int ix0
// Return Type  : double
//
static double i_xnrm2(int n, const double x[9], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                const double x[36]
//                int ix0
// Return Type  : double
//
static double j_xnrm2(int n, const double x[36], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = (ix0 + n) - 1;
  for (k = ix0; k <= kend; k++) {
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}

//
// Arguments    : int n
//                const double x[8]
//                int ix0
// Return Type  : double
//
static double k_xnrm2(int n, const double x[8], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n == 1) {
    y = std::abs(x[ix0 - 1]);
  } else {
    scale = 3.3121686421112381E-170;
    kend = ix0 + 1;
    for (k = ix0; k <= kend; k++) {
      absxk = std::abs(x[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * std::sqrt(y);
  }

  return y;
}

//
// Arguments    : int n
//                const double x[4]
//                int ix0
// Return Type  : double
//
static double l_xnrm2(int n, const double x[4], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = (ix0 + n) - 1;
  for (k = ix0; k <= kend; k++) {
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}

//
// Arguments    : int n
//                const float x[64]
//                int ix0
// Return Type  : float
//
static float m_xnrm2(int n, const float x[64], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 1.29246971E-26F;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : const double A[36]
//                const double B[12]
//                double Y[3]
// Return Type  : void
//
static void mldivide(const double A[36], const double B[12], double Y[3])
{
  int jpvt[3];
  double b_A[36];
  double tau[3];
  int rankR;
  double tol;
  double b_B[12];
  int j;
  int i;
  int Y_tmp;
  int b_Y_tmp;
  jpvt[0] = 1;
  jpvt[1] = 2;
  jpvt[2] = 3;
  std::memcpy(&b_A[0], &A[0], 36U * sizeof(double));
  tau[0] = 0.0;
  tau[1] = 0.0;
  tau[2] = 0.0;
  qrpf(b_A, tau, jpvt);
  rankR = 0;
  tol = 2.6645352591003757E-14 * std::abs(b_A[0]);
  while ((rankR < 3) && (std::abs(b_A[rankR + 12 * rankR]) > tol)) {
    rankR++;
  }

  std::memcpy(&b_B[0], &B[0], 12U * sizeof(double));
  for (j = 0; j < 3; j++) {
    Y[j] = 0.0;
    if (tau[j] != 0.0) {
      tol = b_B[j];
      Y_tmp = j + 2;
      for (i = Y_tmp; i < 13; i++) {
        tol += b_A[(i + 12 * j) - 1] * b_B[i - 1];
      }

      tol *= tau[j];
      if (tol != 0.0) {
        b_B[j] -= tol;
        Y_tmp = j + 2;
        for (i = Y_tmp; i < 13; i++) {
          b_B[i - 1] -= b_A[(i + 12 * j) - 1] * tol;
        }
      }
    }
  }

  for (i = 0; i < rankR; i++) {
    Y[jpvt[i] - 1] = b_B[i];
  }

  for (j = rankR; j >= 1; j--) {
    Y_tmp = jpvt[j - 1] - 1;
    b_Y_tmp = 12 * (j - 1);
    Y[Y_tmp] /= b_A[(j + b_Y_tmp) - 1];
    for (i = 0; i <= j - 2; i++) {
      Y[jpvt[i] - 1] -= Y[Y_tmp] * b_A[i + b_Y_tmp];
    }
  }
}

//
// function G20 = mult_for_groebner(G4)
// Arguments    : const double G4[112]
//                double G20[900]
// Return Type  : void
//
static void mult_for_groebner(const double G4[112], double G20[900])
{
  int i;
  double G20_tmp;
  double b_G20_tmp;
  double c_G20_tmp;
  double d_G20_tmp;
  double e_G20_tmp;
  double f_G20_tmp;
  double g_G20_tmp;
  double h_G20_tmp;
  double i_G20_tmp;
  double j_G20_tmp;
  double k_G20_tmp;
  double l_G20_tmp;
  double m_G20_tmp;
  double n_G20_tmp;
  double o_G20_tmp;
  double p_G20_tmp;
  double q_G20_tmp;
  double r_G20_tmp;
  double s_G20_tmp;
  double t_G20_tmp;
  double u_G20_tmp;
  double v_G20_tmp;
  double w_G20_tmp;
  double x_G20_tmp;
  double y_G20_tmp;
  double ab_G20_tmp;
  double bb_G20_tmp;
  int b_i;

  // 'mult_for_groebner:2' G20 = zeros(20, 45, 'like', G4);
  std::memset(&G20[0], 0, 900U * sizeof(double));

  // multiply by qx^2
  // 'mult_for_groebner:5' for i = 1 : 4
  // multiply by qxqy
  // 'mult_for_groebner:37' for i = 1 : 4
  // multiply by qx
  // 'mult_for_groebner:69' for i = 1 : 4
  // multiply by qy
  // 'mult_for_groebner:101' for i = 1 : 4
  // multiply by 1
  // 'mult_for_groebner:133' for i = 1 : 4
  for (i = 0; i < 4; i++) {
    // 'mult_for_groebner:6' G20(i, 1) = G4(i, 1);
    G20[i] = G4[i];

    // 'mult_for_groebner:7' G20(i, 2) = G4(i, 2);
    G20_tmp = G4[i + 4];
    G20[i + 20] = G20_tmp;

    // 'mult_for_groebner:8' G20(i, 3) = G4(i, 3);
    b_G20_tmp = G4[i + 8];
    G20[i + 40] = b_G20_tmp;

    // 'mult_for_groebner:9' G20(i, 4) = G4(i, 4);
    c_G20_tmp = G4[i + 12];
    G20[i + 60] = c_G20_tmp;

    // 'mult_for_groebner:10' G20(i, 5) = G4(i, 5);
    d_G20_tmp = G4[i + 16];
    G20[i + 80] = d_G20_tmp;

    // 'mult_for_groebner:11' G20(i, 6) = G4(i, 6);
    e_G20_tmp = G4[i + 20];
    G20[i + 100] = e_G20_tmp;

    // 'mult_for_groebner:12' G20(i, 7) = G4(i, 7);
    f_G20_tmp = G4[i + 24];
    G20[i + 120] = f_G20_tmp;

    // 'mult_for_groebner:13' G20(i, 10) = G4(i, 8);
    g_G20_tmp = G4[i + 28];
    G20[i + 180] = g_G20_tmp;

    // 'mult_for_groebner:14' G20(i, 11) = G4(i, 9);
    h_G20_tmp = G4[i + 32];
    G20[i + 200] = h_G20_tmp;

    // 'mult_for_groebner:15' G20(i, 12) = G4(i, 10);
    i_G20_tmp = G4[i + 36];
    G20[i + 220] = i_G20_tmp;

    // 'mult_for_groebner:16' G20(i, 13) = G4(i, 11);
    j_G20_tmp = G4[i + 40];
    G20[i + 240] = j_G20_tmp;

    // 'mult_for_groebner:17' G20(i, 14) = G4(i, 12);
    k_G20_tmp = G4[i + 44];
    G20[i + 260] = k_G20_tmp;

    // 'mult_for_groebner:18' G20(i, 15) = G4(i, 13);
    l_G20_tmp = G4[i + 48];
    G20[i + 280] = l_G20_tmp;

    // 'mult_for_groebner:19' G20(i, 18) = G4(i, 14);
    m_G20_tmp = G4[i + 52];
    G20[i + 340] = m_G20_tmp;

    // 'mult_for_groebner:20' G20(i, 19) = G4(i, 15);
    n_G20_tmp = G4[i + 56];
    G20[i + 360] = n_G20_tmp;

    // 'mult_for_groebner:21' G20(i, 20) = G4(i, 16);
    o_G20_tmp = G4[i + 60];
    G20[i + 380] = o_G20_tmp;

    // 'mult_for_groebner:22' G20(i, 21) = G4(i, 17);
    p_G20_tmp = G4[i + 64];
    G20[i + 400] = p_G20_tmp;

    // 'mult_for_groebner:23' G20(i, 22) = G4(i, 18);
    q_G20_tmp = G4[i + 68];
    G20[i + 420] = q_G20_tmp;

    // 'mult_for_groebner:24' G20(i, 25) = G4(i, 19);
    r_G20_tmp = G4[i + 72];
    G20[i + 480] = r_G20_tmp;

    // 'mult_for_groebner:25' G20(i, 26) = G4(i, 20);
    s_G20_tmp = G4[i + 76];
    G20[i + 500] = s_G20_tmp;

    // 'mult_for_groebner:26' G20(i, 27) = G4(i, 21);
    t_G20_tmp = G4[i + 80];
    G20[i + 520] = t_G20_tmp;

    // 'mult_for_groebner:27' G20(i, 28) = G4(i, 22);
    u_G20_tmp = G4[i + 84];
    G20[i + 540] = u_G20_tmp;

    // 'mult_for_groebner:28' G20(i, 31) = G4(i, 23);
    v_G20_tmp = G4[i + 88];
    G20[i + 600] = v_G20_tmp;

    // 'mult_for_groebner:29' G20(i, 32) = G4(i, 24);
    w_G20_tmp = G4[i + 92];
    G20[i + 620] = w_G20_tmp;

    // 'mult_for_groebner:30' G20(i, 33) = G4(i, 25);
    x_G20_tmp = G4[i + 96];
    G20[i + 640] = x_G20_tmp;

    // 'mult_for_groebner:31' G20(i, 36) = G4(i, 26);
    y_G20_tmp = G4[i + 100];
    G20[i + 700] = y_G20_tmp;

    // 'mult_for_groebner:32' G20(i, 37) = G4(i, 27);
    ab_G20_tmp = G4[i + 104];
    G20[i + 720] = ab_G20_tmp;

    // 'mult_for_groebner:33' G20(i, 40) = G4(i, 28);
    bb_G20_tmp = G4[i + 108];
    G20[i + 780] = bb_G20_tmp;

    // 'mult_for_groebner:38' G20(i + 4, 2) = G4(i, 1);
    G20[i + 24] = G4[i];

    // 'mult_for_groebner:39' G20(i + 4, 3) = G4(i, 2);
    G20[i + 44] = G20_tmp;

    // 'mult_for_groebner:40' G20(i + 4, 4) = G4(i, 3);
    G20[i + 64] = b_G20_tmp;

    // 'mult_for_groebner:41' G20(i + 4, 5) = G4(i, 4);
    G20[i + 84] = c_G20_tmp;

    // 'mult_for_groebner:42' G20(i + 4, 6) = G4(i, 5);
    G20[i + 104] = d_G20_tmp;

    // 'mult_for_groebner:43' G20(i + 4, 7) = G4(i, 6);
    G20[i + 124] = e_G20_tmp;

    // 'mult_for_groebner:44' G20(i + 4, 8) = G4(i, 7);
    G20[i + 144] = f_G20_tmp;

    // 'mult_for_groebner:45' G20(i + 4, 11) = G4(i, 8);
    G20[i + 204] = g_G20_tmp;

    // 'mult_for_groebner:46' G20(i + 4, 12) = G4(i, 9);
    G20[i + 224] = h_G20_tmp;

    // 'mult_for_groebner:47' G20(i + 4, 13) = G4(i, 10);
    G20[i + 244] = i_G20_tmp;

    // 'mult_for_groebner:48' G20(i + 4, 14) = G4(i, 11);
    G20[i + 264] = j_G20_tmp;

    // 'mult_for_groebner:49' G20(i + 4, 15) = G4(i, 12);
    G20[i + 284] = k_G20_tmp;

    // 'mult_for_groebner:50' G20(i + 4, 16) = G4(i, 13);
    G20[i + 304] = l_G20_tmp;

    // 'mult_for_groebner:51' G20(i + 4, 19) = G4(i, 14);
    G20[i + 364] = m_G20_tmp;

    // 'mult_for_groebner:52' G20(i + 4, 20) = G4(i, 15);
    G20[i + 384] = n_G20_tmp;

    // 'mult_for_groebner:53' G20(i + 4, 21) = G4(i, 16);
    G20[i + 404] = o_G20_tmp;

    // 'mult_for_groebner:54' G20(i + 4, 22) = G4(i, 17);
    G20[i + 424] = p_G20_tmp;

    // 'mult_for_groebner:55' G20(i + 4, 23) = G4(i, 18);
    G20[i + 444] = q_G20_tmp;

    // 'mult_for_groebner:56' G20(i + 4, 26) = G4(i, 19);
    G20[i + 504] = r_G20_tmp;

    // 'mult_for_groebner:57' G20(i + 4, 27) = G4(i, 20);
    G20[i + 524] = s_G20_tmp;

    // 'mult_for_groebner:58' G20(i + 4, 28) = G4(i, 21);
    G20[i + 544] = t_G20_tmp;

    // 'mult_for_groebner:59' G20(i + 4, 29) = G4(i, 22);
    G20[i + 564] = u_G20_tmp;

    // 'mult_for_groebner:60' G20(i + 4, 32) = G4(i, 23);
    G20[i + 624] = v_G20_tmp;

    // 'mult_for_groebner:61' G20(i + 4, 33) = G4(i, 24);
    G20[i + 644] = w_G20_tmp;

    // 'mult_for_groebner:62' G20(i + 4, 34) = G4(i, 25);
    G20[i + 664] = x_G20_tmp;

    // 'mult_for_groebner:63' G20(i + 4, 37) = G4(i, 26);
    G20[i + 724] = y_G20_tmp;

    // 'mult_for_groebner:64' G20(i + 4, 38) = G4(i, 27);
    G20[i + 744] = ab_G20_tmp;

    // 'mult_for_groebner:65' G20(i + 4, 41) = G4(i, 28);
    G20[i + 804] = bb_G20_tmp;

    // 'mult_for_groebner:70' G20(i + 8, 10) = G4(i, 1);
    G20[i + 188] = G4[i];

    // 'mult_for_groebner:71' G20(i + 8, 11) = G4(i, 2);
    G20[i + 208] = G20_tmp;

    // 'mult_for_groebner:72' G20(i + 8, 12) = G4(i, 3);
    G20[i + 228] = b_G20_tmp;

    // 'mult_for_groebner:73' G20(i + 8, 13) = G4(i, 4);
    G20[i + 248] = c_G20_tmp;

    // 'mult_for_groebner:74' G20(i + 8, 14) = G4(i, 5);
    G20[i + 268] = d_G20_tmp;

    // 'mult_for_groebner:75' G20(i + 8, 15) = G4(i, 6);
    G20[i + 288] = e_G20_tmp;

    // 'mult_for_groebner:76' G20(i + 8, 16) = G4(i, 7);
    G20[i + 308] = f_G20_tmp;

    // 'mult_for_groebner:77' G20(i + 8, 18) = G4(i, 8);
    G20[i + 348] = g_G20_tmp;

    // 'mult_for_groebner:78' G20(i + 8, 19) = G4(i, 9);
    G20[i + 368] = h_G20_tmp;

    // 'mult_for_groebner:79' G20(i + 8, 20) = G4(i, 10);
    G20[i + 388] = i_G20_tmp;

    // 'mult_for_groebner:80' G20(i + 8, 21) = G4(i, 11);
    G20[i + 408] = j_G20_tmp;

    // 'mult_for_groebner:81' G20(i + 8, 22) = G4(i, 12);
    G20[i + 428] = k_G20_tmp;

    // 'mult_for_groebner:82' G20(i + 8, 23) = G4(i, 13);
    G20[i + 448] = l_G20_tmp;

    // 'mult_for_groebner:83' G20(i + 8, 25) = G4(i, 14);
    G20[i + 488] = m_G20_tmp;

    // 'mult_for_groebner:84' G20(i + 8, 26) = G4(i, 15);
    G20[i + 508] = n_G20_tmp;

    // 'mult_for_groebner:85' G20(i + 8, 27) = G4(i, 16);
    G20[i + 528] = o_G20_tmp;

    // 'mult_for_groebner:86' G20(i + 8, 28) = G4(i, 17);
    G20[i + 548] = p_G20_tmp;

    // 'mult_for_groebner:87' G20(i + 8, 29) = G4(i, 18);
    G20[i + 568] = q_G20_tmp;

    // 'mult_for_groebner:88' G20(i + 8, 31) = G4(i, 19);
    G20[i + 608] = r_G20_tmp;

    // 'mult_for_groebner:89' G20(i + 8, 32) = G4(i, 20);
    G20[i + 628] = s_G20_tmp;

    // 'mult_for_groebner:90' G20(i + 8, 33) = G4(i, 21);
    G20[i + 648] = t_G20_tmp;

    // 'mult_for_groebner:91' G20(i + 8, 34) = G4(i, 22);
    G20[i + 668] = u_G20_tmp;

    // 'mult_for_groebner:92' G20(i + 8, 36) = G4(i, 23);
    G20[i + 708] = v_G20_tmp;

    // 'mult_for_groebner:93' G20(i + 8, 37) = G4(i, 24);
    G20[i + 728] = w_G20_tmp;

    // 'mult_for_groebner:94' G20(i + 8, 38) = G4(i, 25);
    G20[i + 748] = x_G20_tmp;

    // 'mult_for_groebner:95' G20(i + 8, 40) = G4(i, 26);
    G20[i + 788] = y_G20_tmp;

    // 'mult_for_groebner:96' G20(i + 8, 41) = G4(i, 27);
    G20[i + 808] = ab_G20_tmp;

    // 'mult_for_groebner:97' G20(i + 8, 43) = G4(i, 28);
    G20[i + 848] = bb_G20_tmp;

    // 'mult_for_groebner:102' G20(i + 12, 11) = G4(i, 1);
    G20[i + 212] = G4[i];

    // 'mult_for_groebner:103' G20(i + 12, 12) = G4(i, 2);
    G20[i + 232] = G20_tmp;

    // 'mult_for_groebner:104' G20(i + 12, 13) = G4(i, 3);
    G20[i + 252] = b_G20_tmp;

    // 'mult_for_groebner:105' G20(i + 12, 14) = G4(i, 4);
    G20[i + 272] = c_G20_tmp;

    // 'mult_for_groebner:106' G20(i + 12, 15) = G4(i, 5);
    G20[i + 292] = d_G20_tmp;

    // 'mult_for_groebner:107' G20(i + 12, 16) = G4(i, 6);
    G20[i + 312] = e_G20_tmp;

    // 'mult_for_groebner:108' G20(i + 12, 17) = G4(i, 7);
    G20[i + 332] = f_G20_tmp;

    // 'mult_for_groebner:109' G20(i + 12, 19) = G4(i, 8);
    G20[i + 372] = g_G20_tmp;

    // 'mult_for_groebner:110' G20(i + 12, 20) = G4(i, 9);
    G20[i + 392] = h_G20_tmp;

    // 'mult_for_groebner:111' G20(i + 12, 21) = G4(i, 10);
    G20[i + 412] = i_G20_tmp;

    // 'mult_for_groebner:112' G20(i + 12, 22) = G4(i, 11);
    G20[i + 432] = j_G20_tmp;

    // 'mult_for_groebner:113' G20(i + 12, 23) = G4(i, 12);
    G20[i + 452] = k_G20_tmp;

    // 'mult_for_groebner:114' G20(i + 12, 24) = G4(i, 13);
    G20[i + 472] = l_G20_tmp;

    // 'mult_for_groebner:115' G20(i + 12, 26) = G4(i, 14);
    G20[i + 512] = m_G20_tmp;

    // 'mult_for_groebner:116' G20(i + 12, 27) = G4(i, 15);
    G20[i + 532] = n_G20_tmp;

    // 'mult_for_groebner:117' G20(i + 12, 28) = G4(i, 16);
    G20[i + 552] = o_G20_tmp;

    // 'mult_for_groebner:118' G20(i + 12, 29) = G4(i, 17);
    G20[i + 572] = p_G20_tmp;

    // 'mult_for_groebner:119' G20(i + 12, 30) = G4(i, 18);
    G20[i + 592] = q_G20_tmp;

    // 'mult_for_groebner:120' G20(i + 12, 32) = G4(i, 19);
    G20[i + 632] = r_G20_tmp;

    // 'mult_for_groebner:121' G20(i + 12, 33) = G4(i, 20);
    G20[i + 652] = s_G20_tmp;

    // 'mult_for_groebner:122' G20(i + 12, 34) = G4(i, 21);
    G20[i + 672] = t_G20_tmp;

    // 'mult_for_groebner:123' G20(i + 12, 35) = G4(i, 22);
    G20[i + 692] = u_G20_tmp;

    // 'mult_for_groebner:124' G20(i + 12, 37) = G4(i, 23);
    G20[i + 732] = v_G20_tmp;

    // 'mult_for_groebner:125' G20(i + 12, 38) = G4(i, 24);
    G20[i + 752] = w_G20_tmp;

    // 'mult_for_groebner:126' G20(i + 12, 39) = G4(i, 25);
    G20[i + 772] = x_G20_tmp;

    // 'mult_for_groebner:127' G20(i + 12, 41) = G4(i, 26);
    G20[i + 812] = y_G20_tmp;

    // 'mult_for_groebner:128' G20(i + 12, 42) = G4(i, 27);
    G20[i + 832] = ab_G20_tmp;

    // 'mult_for_groebner:129' G20(i + 12, 44) = G4(i, 28);
    G20[i + 872] = bb_G20_tmp;

    // 'mult_for_groebner:134' G20(i + 16, 18:end) = G4(i, :);
    for (b_i = 0; b_i < 28; b_i++) {
      G20[(i + 20 * (b_i + 17)) + 16] = G4[i + (b_i << 2)];
    }
  }
}

//
// function c = mult_poly22(a, b)
// Arguments    : const double a[6]
//                const double b[6]
//                double c[15]
// Return Type  : void
//
static void mult_poly22(const double a[6], const double b[6], double c[15])
{
  // 'mult_poly22:2' c = zeros(1, 1, 15);
  // 'mult_poly22:3' c(1) = a(:, :, 1)*b(:, :, 1);
  c[0] = a[0] * b[0];

  // x^4
  // 'mult_poly22:4' c(2) = a(:, :, 1)*b(:, :, 2) + a(:, :, 2)*b(:, :, 1);
  c[1] = a[0] * b[1] + a[1] * b[0];

  // x^3*y
  // 'mult_poly22:5' c(3) = a(:, :, 1)*b(:, :, 3) + a(:, :, 2)*b(:, :, 2) + a(:, :, 3)*b(:, :, 1); 
  c[2] = (a[0] * b[2] + a[1] * b[1]) + a[2] * b[0];

  // x^2*y^2
  // 'mult_poly22:6' c(4) = a(:, :, 2)*b(:, :, 3) + a(:, :, 3)*b(:, :, 2);
  c[3] = a[1] * b[2] + a[2] * b[1];

  // x*y^3
  // 'mult_poly22:7' c(5) = a(:, :, 3)*b(:, :, 3);
  c[4] = a[2] * b[2];

  // y^4
  // 'mult_poly22:8' c(6) = a(:, :, 1)*b(:, :, 4) + a(:, :, 4)*b(:, :, 1);
  c[5] = a[0] * b[3] + a[3] * b[0];

  // x^3
  // 'mult_poly22:9' c(7) = a(:, :, 1)*b(:, :, 5) + a(:, :, 2)*b(:, :, 4) + a(:, :, 4)*b(:, :, 2) + a(:, :, 5)*b(:, :, 1); 
  c[6] = ((a[0] * b[4] + a[1] * b[3]) + a[3] * b[1]) + a[4] * b[0];

  // x^2*y
  // 'mult_poly22:10' c(8) = a(:, :, 2)*b(:, :, 5) + a(:, :, 3)*b(:, :, 4) + a(:, :, 4)*b(:, :, 3) + a(:, :, 5)*b(:, :, 2); 
  c[7] = ((a[1] * b[4] + a[2] * b[3]) + a[3] * b[2]) + a[4] * b[1];

  // x*y^2
  // 'mult_poly22:11' c(9) = a(:, :, 3)*b(:, :, 5) + a(:, :, 5)*b(:, :, 3);
  c[8] = a[2] * b[4] + a[4] * b[2];

  // y^3
  // 'mult_poly22:12' c(10) = a(:, :, 1)*b(:, :, 6) + a(:, :, 6)*b(:, :, 1) + a(:, :, 4)*b(:, :, 4); 
  c[9] = (a[0] * b[5] + a[5] * b[0]) + a[3] * b[3];

  // x^2
  // 'mult_poly22:13' c(11) = a(:, :, 2)*b(:, :, 6) + a(:, :, 6)*b(:, :, 2) + a(:, :, 4)*b(:, :, 5) + a(:, :, 5)*b(:, :, 4); 
  c[10] = ((a[1] * b[5] + a[5] * b[1]) + a[3] * b[4]) + a[4] * b[3];

  // x*y
  // 'mult_poly22:14' c(12) = a(:, :, 3)*b(:, :, 6) + a(:, :, 6)*b(:, :, 3) + a(:, :, 5)*b(:, :, 5); 
  c[11] = (a[2] * b[5] + a[5] * b[2]) + a[4] * b[4];

  // y^2
  // 'mult_poly22:15' c(13) = a(:, :, 4)*b(:, :, 6) + a(:, :, 6)*b(:, :, 4);
  c[12] = a[3] * b[5] + a[5] * b[3];

  // x
  // 'mult_poly22:16' c(14) = a(:, :, 5)*b(:, :, 6) + a(:, :, 6)*b(:, :, 5);
  c[13] = a[4] * b[5] + a[5] * b[4];

  // y
  // 'mult_poly22:17' c(15) = a(:, :, 6)*b(:, :, 6);
  c[14] = a[5] * b[5];

  // 1
}

//
// function p = mult_poly42(d, a)
// Arguments    : const double d[15]
//                const double a[6]
//                double p[28]
// Return Type  : void
//
static void mult_poly42(const double d[15], const double a[6], double p[28])
{
  // 'mult_poly42:2' p = zeros(1, 1, 28);
  // 'mult_poly42:3' p(1) = a(:, :, 1)*d(:, :, 1);
  p[0] = a[0] * d[0];

  // x^6
  // 'mult_poly42:4' p(2) = a(:, :, 1)*d(:, :, 2) + a(:, :, 2)*d(:, :, 1);
  p[1] = a[0] * d[1] + a[1] * d[0];

  // x^5*y
  // 'mult_poly42:5' p(3) = a(:, :, 1)*d(:, :, 3) + a(:, :, 2)*d(:, :, 2) + a(:, :, 3)*d(:, :, 1); 
  p[2] = (a[0] * d[2] + a[1] * d[1]) + a[2] * d[0];

  // x^4*y^2
  // 'mult_poly42:6' p(4) = a(:, :, 1)*d(:, :, 4) + a(:, :, 2)*d(:, :, 3) + a(:, :, 3)*d(:, :, 2); 
  p[3] = (a[0] * d[3] + a[1] * d[2]) + a[2] * d[1];

  // x^3*y^3
  // 'mult_poly42:7' p(5) = a(:, :, 1)*d(:, :, 5) + a(:, :, 2)*d(:, :, 4) + a(:, :, 3)*d(:, :, 3); 
  p[4] = (a[0] * d[4] + a[1] * d[3]) + a[2] * d[2];

  // x^2*y^4
  // 'mult_poly42:8' p(6) = a(:, :, 2)*d(:, :, 5) + a(:, :, 3)*d(:, :, 4);
  p[5] = a[1] * d[4] + a[2] * d[3];

  // x*y^5
  // 'mult_poly42:9' p(7) = a(:, :, 3)*d(:, :, 5);
  p[6] = a[2] * d[4];

  // y^6
  // 'mult_poly42:10' p(8) = a(:, :, 4)*d(:, :, 1) + a(:, :, 1)*d(:, :, 6);
  p[7] = a[3] * d[0] + a[0] * d[5];

  // x^5
  // 'mult_poly42:11' p(9) = a(:, :, 4)*d(:, :, 2) + a(:, :, 5)*d(:, :, 1) + a(:, :, 1)*d(:, :, 7) + a(:, :, 2)*d(:, :, 6); 
  p[8] = ((a[3] * d[1] + a[4] * d[0]) + a[0] * d[6]) + a[1] * d[5];

  // x^4*y
  // 'mult_poly42:12' p(10) = a(:, :, 4)*d(:, :, 3) + a(:, :, 5)*d(:, :, 2) + a(:, :, 1)*d(:, :, 8) + a(:, :, 2)*d(:, :, 7) + a(:, :, 3)*d(:, :, 6); 
  p[9] = (((a[3] * d[2] + a[4] * d[1]) + a[0] * d[7]) + a[1] * d[6]) + a[2] * d
    [5];

  // x^3*y^2
  // 'mult_poly42:13' p(11) = a(:, :, 4)*d(:, :, 4) + a(:, :, 5)*d(:, :, 3) + a(:, :, 1)*d(:, :, 9) + a(:, :, 2)*d(:, :, 8) + a(:, :, 3)*d(:, :, 7); 
  p[10] = (((a[3] * d[3] + a[4] * d[2]) + a[0] * d[8]) + a[1] * d[7]) + a[2] *
    d[6];

  // x^2*y^3
  // 'mult_poly42:14' p(12) = a(:, :, 4)*d(:, :, 5) + a(:, :, 5)*d(:, :, 4) + a(:, :, 2)*d(:, :, 9) + a(:, :, 3)*d(:, :, 8); 
  p[11] = ((a[3] * d[4] + a[4] * d[3]) + a[1] * d[8]) + a[2] * d[7];

  // x*y^4
  // 'mult_poly42:15' p(13) = a(:, :, 5)*d(:, :, 5) + a(:, :, 3)*d(:, :, 9);
  p[12] = a[4] * d[4] + a[2] * d[8];

  // y^5
  // 'mult_poly42:16' p(14) = a(:, :, 6)*d(:, :, 1) + a(:, :, 4)*d(:, :, 6) + a(:, :, 1)*d(:, :, 10); 
  p[13] = (a[5] * d[0] + a[3] * d[5]) + a[0] * d[9];

  // x^4
  // 'mult_poly42:17' p(15) = a(:, :, 6)*d(:, :, 2) + a(:, :, 4)*d(:, :, 7) + a(:, :, 5)*d(:, :, 6) + a(:, :, 1)*d(:, :, 11) + a(:, :, 2)*d(:, :, 10); 
  p[14] = (((a[5] * d[1] + a[3] * d[6]) + a[4] * d[5]) + a[0] * d[10]) + a[1] *
    d[9];

  // x^3*y
  // 'mult_poly42:18' p(16) = a(:, :, 6)*d(:, :, 3) + a(:, :, 4)*d(:, :, 8) + a(:, :, 5)*d(:, :, 7) + a(:, :, 1)*d(:, :, 12) + a(:, :, 2)*d(:, :, 11) + a(:, :, 3)*d(:, :, 10); 
  p[15] = ((((a[5] * d[2] + a[3] * d[7]) + a[4] * d[6]) + a[0] * d[11]) + a[1] *
           d[10]) + a[2] * d[9];

  // x^2*y^2
  // 'mult_poly42:19' p(17) = a(:, :, 6)*d(:, :, 4) + a(:, :, 4)*d(:, :, 9) + a(:, :, 5)*d(:, :, 8) + a(:, :, 2)*d(:, :, 12) + a(:, :, 3)*d(:, :, 11); 
  p[16] = (((a[5] * d[3] + a[3] * d[8]) + a[4] * d[7]) + a[1] * d[11]) + a[2] *
    d[10];

  // x*y^3
  // 'mult_poly42:20' p(18) = a(:, :, 6)*d(:, :, 5) + a(:, :, 5)*d(:, :, 9) + a(:, :, 3)*d(:, :, 12); 
  p[17] = (a[5] * d[4] + a[4] * d[8]) + a[2] * d[11];

  // y^4
  // 'mult_poly42:21' p(19) = a(:, :, 6)*d(:, :, 6) + a(:, :, 1)*d(:, :, 13) + a(:, :, 4)*d(:, :, 10); 
  p[18] = (a[5] * d[5] + a[0] * d[12]) + a[3] * d[9];

  // x^3
  // 'mult_poly42:22' p(20) = a(:, :, 6)*d(:, :, 7) + a(:, :, 1)*d(:, :, 14) + a(:, :, 2)*d(:, :, 13) + a(:, :, 4)*d(:, :, 11) + a(:, :, 5)*d(:, :, 10); 
  p[19] = (((a[5] * d[6] + a[0] * d[13]) + a[1] * d[12]) + a[3] * d[10]) + a[4] *
    d[9];

  // x^2*y
  // 'mult_poly42:23' p(21) = a(:, :, 6)*d(:, :, 8) + a(:, :, 2)*d(:, :, 14) + a(:, :, 3)*d(:, :, 13) + a(:, :, 4)*d(:, :, 12) + a(:, :, 5)*d(:, :, 11); 
  p[20] = (((a[5] * d[7] + a[1] * d[13]) + a[2] * d[12]) + a[3] * d[11]) + a[4] *
    d[10];

  // x*y^2
  // 'mult_poly42:24' p(22) = a(:, :, 6)*d(:, :, 9) + a(:, :, 3)*d(:, :, 14) + a(:, :, 5)*d(:, :, 12); 
  p[21] = (a[5] * d[8] + a[2] * d[13]) + a[4] * d[11];

  // y^3
  // 'mult_poly42:25' p(23) = a(:, :, 1)*d(:, :, 15) + a(:, :, 6)*d(:, :, 10) + a(:, :, 4)*d(:, :, 13); 
  p[22] = (a[0] * d[14] + a[5] * d[9]) + a[3] * d[12];

  // x^2
  // 'mult_poly42:26' p(24) = a(:, :, 2)*d(:, :, 15) + a(:, :, 6)*d(:, :, 11) + a(:, :, 4)*d(:, :, 14) + a(:, :, 5)*d(:, :, 13); 
  p[23] = ((a[1] * d[14] + a[5] * d[10]) + a[3] * d[13]) + a[4] * d[12];

  // x*y
  // 'mult_poly42:27' p(25) = a(:, :, 3)*d(:, :, 15) + a(:, :, 6)*d(:, :, 12) + a(:, :, 5)*d(:, :, 14); 
  p[24] = (a[2] * d[14] + a[5] * d[11]) + a[4] * d[13];

  // y^2
  // 'mult_poly42:28' p(26) = a(:, :, 4)*d(:, :, 15) + a(:, :, 6)*d(:, :, 13);
  p[25] = a[3] * d[14] + a[5] * d[12];

  // x
  // 'mult_poly42:29' p(27) = a(:, :, 5)*d(:, :, 15) + a(:, :, 6)*d(:, :, 14);
  p[26] = a[4] * d[14] + a[5] * d[13];

  // y
  // 'mult_poly42:30' p(28) = a(:, :, 6)*d(:, :, 15);
  p[27] = a[5] * d[14];

  // 1
}

//
// Arguments    : int n
//                const creal32_T x_data[]
//                int ix0
// Return Type  : float
//
static float n_xnrm2(int n, const creal32_T x_data[], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = rt_hypotf(x_data[ix0 - 1].re, x_data[ix0 - 1].im);
    } else {
      scale = 1.29246971E-26F;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x_data[k - 1].re);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }

        absxk = std::abs(x_data[k - 1].im);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                const float x[9]
//                int ix0
// Return Type  : float
//
static float o_xnrm2(int n, const float x[9], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 1.29246971E-26F;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0F;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                const float x[8]
//                int ix0
// Return Type  : float
//
static float p_xnrm2(int n, const float x[8], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  if (n == 1) {
    y = std::abs(x[ix0 - 1]);
  } else {
    scale = 1.29246971E-26F;
    kend = ix0 + 1;
    for (k = ix0; k <= kend; k++) {
      absxk = std::abs(x[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0F;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * std::sqrt(y);
  }

  return y;
}

//
// Arguments    : int n
//                const float x[4]
//                int ix0
// Return Type  : float
//
static float q_xnrm2(int n, const float x[4], int ix0)
{
  float y;
  float scale;
  int kend;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  scale = 1.29246971E-26F;
  kend = (ix0 + n) - 1;
  for (k = ix0; k <= kend; k++) {
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0F;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}

//
// Arguments    : const double A[12]
//                double Q[9]
//                double R[12]
// Return Type  : void
//
static void qr(const double A[12], double Q[9], double R[12])
{
  double b_A[12];
  double tau[3];
  double work[4];
  int i;
  int ix0;
  int ii;
  double atmp;
  int b_i;
  int knt;
  double c;
  int lastv;
  double beta1;
  int lastc;
  int iaii;
  boolean_T exitg2;
  int k;
  int ia;
  int exitg1;
  int ix;
  int i1;
  std::memcpy(&b_A[0], &A[0], 12U * sizeof(double));
  tau[0] = 0.0;
  tau[1] = 0.0;
  tau[2] = 0.0;
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 0; i < 3; i++) {
    ii = i * 3 + i;
    if (i + 1 < 3) {
      atmp = b_A[ii];
      ix0 = ii + 2;
      tau[i] = 0.0;
      c = c_xnrm2(2 - i, b_A, ii + 2);
      if (c != 0.0) {
        beta1 = rt_hypotd(b_A[ii], c);
        if (b_A[ii] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          b_i = (ii - i) + 3;
          do {
            knt++;
            for (k = ix0; k <= b_i; k++) {
              b_A[k - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd(atmp, c_xnrm2(2 - i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          c = 1.0 / (atmp - beta1);
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - b_A[ii]) / beta1;
          c = 1.0 / (b_A[ii] - beta1);
          b_i = (ii - i) + 3;
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = atmp;
    } else {
      tau[2] = 0.0;
    }

    atmp = b_A[ii];
    b_A[ii] = 1.0;
    if (tau[i] != 0.0) {
      lastv = 3 - i;
      knt = (ii - i) + 2;
      while ((lastv > 0) && (b_A[knt] == 0.0)) {
        lastv--;
        knt--;
      }

      lastc = 3 - i;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        knt = (ii + (lastc - 1) * 3) + 3;
        ia = knt;
        do {
          exitg1 = 0;
          if (ia + 1 <= knt + lastv) {
            if (b_A[ia] != 0.0) {
              exitg1 = 1;
            } else {
              ia++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      knt = ii + 4;
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work[0], 0, lastc * sizeof(double));
        }

        ix0 = 0;
        b_i = (ii + 3 * (lastc - 1)) + 4;
        for (k = knt; k <= b_i; k += 3) {
          ix = ii;
          c = 0.0;
          i1 = (k + lastv) - 1;
          for (ia = k; ia <= i1; ia++) {
            c += b_A[ia - 1] * b_A[ix];
            ix++;
          }

          work[ix0] += c;
          ix0++;
        }
      }

      xgerc(lastv, lastc, -tau[i], ii + 1, work, b_A, ii + 4);
    }

    b_A[ii] = atmp;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    for (i = 0; i <= ix0; i++) {
      knt = i + 3 * ix0;
      R[knt] = b_A[knt];
    }

    b_i = ix0 + 2;
    if (b_i <= 3) {
      std::memset(&R[(ix0 * 3 + b_i) + -1], 0, (4 - b_i) * sizeof(double));
    }
  }

  R[9] = b_A[9];
  R[10] = b_A[10];
  R[11] = b_A[11];
  ii = 2;
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 2; i >= 0; i--) {
    iaii = (i + i * 3) + 4;
    if (i + 1 < 3) {
      b_A[iaii - 4] = 1.0;
      if (tau[ii] != 0.0) {
        lastv = 3 - i;
        knt = iaii - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0)) {
          lastv--;
          knt--;
        }

        lastc = 2 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = iaii + (lastc - 1) * 3;
          ia = knt;
          do {
            exitg1 = 0;
            if (ia <= (knt + lastv) - 1) {
              if (b_A[ia - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            std::memset(&work[0], 0, lastc * sizeof(double));
          }

          ix0 = 0;
          b_i = iaii + 3 * (lastc - 1);
          for (k = iaii; k <= b_i; k += 3) {
            ix = iaii;
            c = 0.0;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix - 4];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        xgerc(lastv, lastc, -tau[ii], iaii - 3, work, b_A, iaii);
      }

      ix0 = iaii - 2;
      b_i = (iaii - i) - 1;
      for (k = ix0; k <= b_i; k++) {
        b_A[k - 1] *= -tau[ii];
      }
    }

    b_A[iaii - 4] = 1.0 - tau[ii];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(iaii - ix0) - 5] = 0.0;
    }

    ii--;
  }

  for (ix0 = 0; ix0 < 3; ix0++) {
    Q[3 * ix0] = b_A[3 * ix0];
    knt = 3 * ix0 + 1;
    Q[knt] = b_A[knt];
    knt = 3 * ix0 + 2;
    Q[knt] = b_A[knt];
  }
}

//
// Arguments    : double A[36]
//                double tau[3]
//                int jpvt[3]
// Return Type  : void
//
static void qrpf(double A[36], double tau[3], int jpvt[3])
{
  int j;
  int i;
  double work[3];
  int pvt;
  int ip1;
  double smax;
  int ii;
  double scale;
  int kend;
  int k;
  int ix;
  double absxk;
  double vn1[3];
  double vn2[3];
  double t;
  int b_i;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int exitg1;
  int i1;
  for (j = 0; j < 3; j++) {
    work[j] = 0.0;
    pvt = j * 12;
    smax = 0.0;
    scale = 3.3121686421112381E-170;
    kend = pvt + 12;
    for (k = pvt + 1; k <= kend; k++) {
      absxk = std::abs(A[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        smax = smax * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        smax += t * t;
      }
    }

    smax = scale * std::sqrt(smax);
    vn1[j] = smax;
    vn2[j] = smax;
  }

  for (i = 0; i < 3; i++) {
    ip1 = i + 2;
    ii = i * 12 + i;
    kend = 3 - i;
    pvt = 0;
    if (3 - i > 1) {
      ix = i;
      smax = std::abs(vn1[i]);
      for (k = 2; k <= kend; k++) {
        ix++;
        scale = std::abs(vn1[ix]);
        if (scale > smax) {
          pvt = k - 1;
          smax = scale;
        }
      }
    }

    pvt += i;
    if (pvt != i) {
      ix = pvt * 12;
      j = i * 12;
      for (k = 0; k < 12; k++) {
        smax = A[ix];
        A[ix] = A[j];
        A[j] = smax;
        ix++;
        j++;
      }

      kend = jpvt[pvt];
      jpvt[pvt] = jpvt[i];
      jpvt[i] = kend;
      vn1[pvt] = vn1[i];
      vn2[pvt] = vn2[i];
    }

    absxk = A[ii];
    pvt = ii + 2;
    tau[i] = 0.0;
    smax = j_xnrm2(11 - i, A, ii + 2);
    if (smax != 0.0) {
      scale = rt_hypotd(A[ii], smax);
      if (A[ii] >= 0.0) {
        scale = -scale;
      }

      if (std::abs(scale) < 1.0020841800044864E-292) {
        kend = -1;
        b_i = (ii - i) + 12;
        do {
          kend++;
          for (k = pvt; k <= b_i; k++) {
            A[k - 1] *= 9.9792015476736E+291;
          }

          scale *= 9.9792015476736E+291;
          absxk *= 9.9792015476736E+291;
        } while (!(std::abs(scale) >= 1.0020841800044864E-292));

        scale = rt_hypotd(absxk, j_xnrm2(11 - i, A, ii + 2));
        if (absxk >= 0.0) {
          scale = -scale;
        }

        tau[i] = (scale - absxk) / scale;
        smax = 1.0 / (absxk - scale);
        for (k = pvt; k <= b_i; k++) {
          A[k - 1] *= smax;
        }

        for (k = 0; k <= kend; k++) {
          scale *= 1.0020841800044864E-292;
        }

        absxk = scale;
      } else {
        tau[i] = (scale - A[ii]) / scale;
        smax = 1.0 / (A[ii] - scale);
        b_i = (ii - i) + 12;
        for (k = pvt; k <= b_i; k++) {
          A[k - 1] *= smax;
        }

        absxk = scale;
      }
    }

    A[ii] = absxk;
    if (i + 1 < 3) {
      absxk = A[ii];
      A[ii] = 1.0;
      k = ii + 13;
      if (tau[i] != 0.0) {
        lastv = 12 - i;
        kend = (ii - i) + 11;
        while ((lastv > 0) && (A[kend] == 0.0)) {
          lastv--;
          kend--;
        }

        lastc = 1 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc + 1 > 0)) {
          kend = (ii + lastc * 12) + 12;
          pvt = kend;
          do {
            exitg1 = 0;
            if (pvt + 1 <= kend + lastv) {
              if (A[pvt] != 0.0) {
                exitg1 = 1;
              } else {
                pvt++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (0 <= lastc) {
            std::memset(&work[0], 0, (lastc + 1) * sizeof(double));
          }

          j = 0;
          b_i = (ii + 12 * lastc) + 13;
          for (kend = k; kend <= b_i; kend += 12) {
            ix = ii;
            smax = 0.0;
            i1 = (kend + lastv) - 1;
            for (pvt = kend; pvt <= i1; pvt++) {
              smax += A[pvt - 1] * A[ix];
              ix++;
            }

            work[j] += smax;
            j++;
          }
        }

        if (-tau[i] != 0.0) {
          kend = ii;
          pvt = 0;
          for (j = 0; j <= lastc; j++) {
            if (work[pvt] != 0.0) {
              smax = work[pvt] * -tau[i];
              ix = ii;
              b_i = kend + 13;
              i1 = lastv + kend;
              for (k = b_i; k <= i1 + 12; k++) {
                A[k - 1] += A[ix] * smax;
                ix++;
              }
            }

            pvt++;
            kend += 12;
          }
        }
      }

      A[ii] = absxk;
    }

    for (j = ip1; j < 4; j++) {
      kend = i + (j - 1) * 12;
      smax = vn1[j - 1];
      if (smax != 0.0) {
        scale = std::abs(A[kend]) / smax;
        scale = 1.0 - scale * scale;
        if (scale < 0.0) {
          scale = 0.0;
        }

        absxk = smax / vn2[j - 1];
        absxk = scale * (absxk * absxk);
        if (absxk <= 1.4901161193847656E-8) {
          smax = j_xnrm2(11 - i, A, kend + 2);
          vn1[j - 1] = smax;
          vn2[j - 1] = smax;
        } else {
          vn1[j - 1] = smax * std::sqrt(scale);
        }
      }
    }
  }
}

//
// function F_row = quadruple_constraint(i, j, k, x, y, X, R)
// Arguments    : double i
//                double j
//                double k
//                const double x[4]
//                const double y[4]
//                const double X[12]
//                const double R[54]
//                double F_row[18]
// Return Type  : void
//
static void quadruple_constraint(double i, double j, double k, const double x[4],
  const double y[4], const double X[12], const double R[54], double F_row[18])
{
  int F_row_tmp_tmp_tmp;
  int b_i;
  int b_F_row_tmp_tmp_tmp;
  int F_row_tmp_tmp;
  double F_row_tmp[3];
  int c_F_row_tmp_tmp_tmp;
  int b_F_row_tmp_tmp;
  double b_F_row_tmp[3];
  int c_F_row_tmp_tmp;
  double a_tmp_tmp;
  double b[6];
  double a_tmp;
  double b_b[6];
  double a;

  // 'quadruple_constraint:2' F_row = zeros(1, 3, 6, 'like', X);
  // fc:
  // 'quadruple_constraint:4' F_row(1, 1, :) = (y(i) - y(k)) * mult_R(R, 1, X(:, j) - X(:, i)) ... 
  // 'quadruple_constraint:5'                         -(x(i) - x(j)) * mult_R(R, 2, X(:, k) - X(:, i)); 
  F_row_tmp_tmp_tmp = static_cast<int>(i) - 1;
  b_i = 3 * F_row_tmp_tmp_tmp;
  b_F_row_tmp_tmp_tmp = static_cast<int>(j) - 1;
  F_row_tmp_tmp = 3 * b_F_row_tmp_tmp_tmp;
  F_row_tmp[0] = X[F_row_tmp_tmp] - X[b_i];
  c_F_row_tmp_tmp_tmp = static_cast<int>(k) - 1;
  b_F_row_tmp_tmp = 3 * c_F_row_tmp_tmp_tmp;
  b_F_row_tmp[0] = X[b_F_row_tmp_tmp] - X[b_i];
  c_F_row_tmp_tmp = b_i + 1;
  F_row_tmp[1] = X[F_row_tmp_tmp + 1] - X[c_F_row_tmp_tmp];
  b_F_row_tmp[1] = X[b_F_row_tmp_tmp + 1] - X[c_F_row_tmp_tmp];
  b_i += 2;
  F_row_tmp[2] = X[F_row_tmp_tmp + 2] - X[b_i];
  b_F_row_tmp[2] = X[b_F_row_tmp_tmp + 2] - X[b_i];
  a_tmp_tmp = y[F_row_tmp_tmp_tmp] - y[c_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b[b_i] = 0.0;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * F_row_tmp[b_i];
    }
  }

  a_tmp = x[F_row_tmp_tmp_tmp] - x[b_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] *
        b_F_row_tmp[b_i];
    }
  }

  // fs:
  // 'quadruple_constraint:7' F_row(1, 2, :) = -(y(i) - y(k)) * mult_R(R, 2, X(:, j) - X(:, i)) ... 
  // 'quadruple_constraint:8'                          -(x(i) - x(j)) * mult_R(R, 1, X(:, k) - X(:, i)); 
  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i] = a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] * F_row_tmp[b_i];
    }
  }

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * b_F_row_tmp[b_i];
    }
  }

  // 1:
  // 'quadruple_constraint:10' F_row(1, 3, :) = -(y(i) - y(k)) * x(j) * mult_R(R, 3, X(:, j) - X(:, i)) ... 
  // 'quadruple_constraint:11'                          +(x(i) - x(j)) * y(k) * mult_R(R, 3, X(:, k) - X(:, i)); 
  a = -a_tmp_tmp * x[b_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i + 1] = -a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 2] * F_row_tmp[b_i];
    }
  }

  a_tmp_tmp = a_tmp * y[c_F_row_tmp_tmp_tmp];

  // 'quadruple_constraint:15' sum = zeros(6,1, 'like', X);
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  // 'quadruple_constraint:16' for j = 1 : 3
  for (b_i = 0; b_i < 3; b_i++) {
    // 'quadruple_constraint:17' sum = sum + squeeze(R(i, j, :)) * X(j);
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 2] *
        b_F_row_tmp[b_i];
    }
  }

  for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
    F_row[3 * F_row_tmp_tmp + 2] = a * b[F_row_tmp_tmp] + a_tmp_tmp *
      b_b[F_row_tmp_tmp];
  }
}

//
// Arguments    : const double A[400]
// Return Type  : double
//
static double rcond(const double A[400])
{
  double result;
  double normA;
  int j;
  double s;
  double b_A[400];
  int i;
  int ipiv[20];
  int jA;
  int jjA;
  int exitg1;
  double ainvnm;
  int iter;
  int kase;
  int jump;
  double x[20];
  int b_j;
  int b_i;
  double absrexk;
  result = 0.0;
  normA = 0.0;
  for (j = 0; j < 20; j++) {
    s = 0.0;
    for (i = 0; i < 20; i++) {
      s += std::abs(A[i + 20 * j]);
    }

    if (s > normA) {
      normA = s;
    }
  }

  if (normA != 0.0) {
    std::memcpy(&b_A[0], &A[0], 400U * sizeof(double));
    xzgetrf(b_A, ipiv, &jA);
    jjA = 19;
    do {
      exitg1 = 0;
      if (jjA + 1 > 0) {
        if (b_A[jjA + 20 * jjA] == 0.0) {
          exitg1 = 1;
        } else {
          jjA--;
        }
      } else {
        ainvnm = 0.0;
        iter = 2;
        kase = 1;
        jump = 1;
        j = 0;
        for (i = 0; i < 20; i++) {
          x[i] = 0.05;
        }

        while (kase != 0) {
          if (kase == 1) {
            for (b_j = 0; b_j < 20; b_j++) {
              jjA = b_j + b_j * 20;
              b_i = 18 - b_j;
              for (i = 0; i <= b_i; i++) {
                jA = (b_j + i) + 1;
                x[jA] -= x[b_j] * b_A[(jjA + i) + 1];
              }
            }

            for (b_j = 19; b_j >= 0; b_j--) {
              jjA = b_j + b_j * 20;
              x[b_j] /= b_A[jjA];
              for (i = 0; i < b_j; i++) {
                jA = (b_j - i) - 1;
                x[jA] -= x[b_j] * b_A[(jjA - i) - 1];
              }
            }
          } else {
            for (b_j = 0; b_j < 20; b_j++) {
              jA = b_j * 20;
              s = x[b_j];
              for (i = 0; i < b_j; i++) {
                s -= b_A[jA + i] * x[i];
              }

              x[b_j] = s / b_A[jA + b_j];
            }

            for (b_j = 19; b_j >= 0; b_j--) {
              jA = b_j * 20;
              s = x[b_j];
              b_i = b_j + 2;
              for (i = 20; i >= b_i; i--) {
                s -= b_A[(jA + i) - 1] * x[i - 1];
              }

              x[b_j] = s;
            }
          }

          if (jump == 1) {
            ainvnm = 0.0;
            for (jjA = 0; jjA < 20; jjA++) {
              s = std::abs(x[jjA]);
              ainvnm += s;
              if (s > 2.2250738585072014E-308) {
                x[jjA] /= s;
              } else {
                x[jjA] = 1.0;
              }
            }

            kase = 2;
            jump = 2;
          } else if (jump == 2) {
            j = 0;
            s = std::abs(x[0]);
            for (jjA = 0; jjA < 19; jjA++) {
              absrexk = std::abs(x[jjA + 1]);
              if (absrexk > s) {
                j = jjA + 1;
                s = absrexk;
              }
            }

            iter = 2;
            std::memset(&x[0], 0, 20U * sizeof(double));
            x[j] = 1.0;
            kase = 1;
            jump = 3;
          } else if (jump == 3) {
            ainvnm = 0.0;
            for (jjA = 0; jjA < 20; jjA++) {
              ainvnm += std::abs(x[jjA]);
            }

            if (ainvnm <= x[0]) {
              s = 1.0;
              for (jjA = 0; jjA < 20; jjA++) {
                x[jjA] = s * (((static_cast<double>(jjA) + 1.0) - 1.0) / 19.0 +
                              1.0);
                s = -s;
              }

              kase = 1;
              jump = 5;
            } else {
              for (jjA = 0; jjA < 20; jjA++) {
                s = std::abs(x[jjA]);
                if (s > 2.2250738585072014E-308) {
                  x[jjA] /= s;
                } else {
                  x[jjA] = 1.0;
                }
              }

              kase = 2;
              jump = 4;
            }
          } else if (jump == 4) {
            jA = j;
            j = 0;
            s = std::abs(x[0]);
            for (jjA = 0; jjA < 19; jjA++) {
              absrexk = std::abs(x[jjA + 1]);
              if (absrexk > s) {
                j = jjA + 1;
                s = absrexk;
              }
            }

            if ((std::abs(x[jA]) != std::abs(x[j])) && (iter <= 5)) {
              iter++;
              std::memset(&x[0], 0, 20U * sizeof(double));
              x[j] = 1.0;
              kase = 1;
              jump = 3;
            } else {
              s = 1.0;
              for (jjA = 0; jjA < 20; jjA++) {
                x[jjA] = s * (((static_cast<double>(jjA) + 1.0) - 1.0) / 19.0 +
                              1.0);
                s = -s;
              }

              kase = 1;
              jump = 5;
            }
          } else {
            if (jump == 5) {
              s = 0.0;
              for (jjA = 0; jjA < 20; jjA++) {
                s += std::abs(x[jjA]);
              }

              s = 2.0 * s / 3.0 / 20.0;
              if (s > ainvnm) {
                ainvnm = s;
              }

              kase = 0;
            }
          }
        }

        if (ainvnm != 0.0) {
          result = 1.0 / ainvnm / normA;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return result;
}

//
// Arguments    : const creal_T y
// Return Type  : creal_T
//
static creal_T recip(const creal_T y)
{
  creal_T z;
  double brm;
  double bim;
  double d;
  brm = std::abs(y.re);
  bim = std::abs(y.im);
  if (y.im == 0.0) {
    z.re = 1.0 / y.re;
    z.im = 0.0;
  } else if (y.re == 0.0) {
    z.re = 0.0;
    z.im = -1.0 / y.im;
  } else if (brm > bim) {
    bim = y.im / y.re;
    d = y.re + bim * y.im;
    z.re = 1.0 / d;
    z.im = -bim / d;
  } else if (brm == bim) {
    bim = 0.5;
    if (y.re < 0.0) {
      bim = -0.5;
    }

    d = 0.5;
    if (y.im < 0.0) {
      d = -0.5;
    }

    z.re = bim / brm;
    z.im = -d / brm;
  } else {
    bim = y.re / y.im;
    d = y.im + bim * y.re;
    z.re = bim / d;
    z.im = -1.0 / d;
  }

  return z;
}

//
// Arguments    : const double c[9]
//                creal_T r_data[]
//                int r_size[1]
// Return Type  : void
//
static void roots(const double c[9], creal_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int nTrailingZeros;
  int companDim;
  boolean_T exitg1;
  int j;
  int a_size[2];
  double ctmp[9];
  creal_T a_data[64];
  int m;
  int i;
  boolean_T p;
  creal_T eiga_data[8];
  int eiga_size[1];
  creal_T beta1_data[8];
  int beta1_size[1];
  int exitg2;
  double brm;
  double re;
  int jend;
  double bim;
  double d;
  std::memset(&r_data[0], 0, 8U * sizeof(creal_T));
  k1 = 1;
  while ((k1 <= 9) && (c[k1 - 1] == 0.0)) {
    k1++;
  }

  k2 = 9;
  while ((k2 >= k1) && (c[k2 - 1] == 0.0)) {
    k2--;
  }

  nTrailingZeros = 8 - k2;
  if (k1 < k2) {
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      for (j = 1; j <= companDim; j++) {
        ctmp[j - 1] = c[(k1 + j) - 1] / c[k1 - 1];
      }

      if (j > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }

    if (companDim < 1) {
      if (1 > 9 - k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = 9 - k2;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      k1 = companDim * companDim;
      if (0 <= k1 - 1) {
        std::memset(&a_data[0], 0, k1 * sizeof(creal_T));
      }

      for (m = 0; m <= companDim - 2; m++) {
        i = companDim * m;
        a_data[i].re = -ctmp[m];
        a_data[i].im = 0.0;
        i = (m + i) + 1;
        a_data[i].re = 1.0;
        a_data[i].im = 0.0;
      }

      i = companDim * (companDim - 1);
      a_data[i].re = -ctmp[companDim - 1];
      a_data[i].im = 0.0;
      if (0 <= nTrailingZeros) {
        std::memset(&r_data[0], 0, (nTrailingZeros + 1) * sizeof(creal_T));
      }

      if (companDim == 1) {
        eiga_data[0] = a_data[0];
      } else {
        p = true;
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j <= companDim - 1)) {
          k1 = 0;
          do {
            exitg2 = 0;
            if (k1 <= j) {
              if (a_data[k1 + companDim * j].re != a_data[j + companDim * k1].re)
              {
                p = false;
                exitg2 = 1;
              } else {
                k1++;
              }
            } else {
              j++;
              exitg2 = 2;
            }
          } while (exitg2 == 0);

          if (exitg2 == 1) {
            exitg1 = true;
          }
        }

        if (p) {
          xgehrd(a_data, a_size);
          eml_zlahqr(a_data, a_size);
          m = a_size[0];
          if (3 < a_size[0]) {
            nTrailingZeros = 4;
            if (a_size[0] - 4 < a_size[1] - 1) {
              jend = a_size[0] - 3;
            } else {
              jend = a_size[1];
            }

            for (j = 0; j < jend; j++) {
              for (k1 = nTrailingZeros; k1 <= m; k1++) {
                i = (k1 + a_size[0] * j) - 1;
                a_data[i].re = 0.0;
                a_data[i].im = 0.0;
              }

              nTrailingZeros++;
            }
          }

          k1 = a_size[0];
          for (m = 0; m < k1; m++) {
            eiga_data[m].re = a_data[m + a_size[0] * m].re;
            eiga_data[m].im = 0.0;
          }
        } else {
          xzgeev(a_data, a_size, &k1, eiga_data, eiga_size, beta1_data,
                 beta1_size);
          k1 = eiga_size[0];
          for (i = 0; i < k1; i++) {
            if (beta1_data[i].im == 0.0) {
              if (eiga_data[i].im == 0.0) {
                re = eiga_data[i].re / beta1_data[i].re;
                bim = 0.0;
              } else if (eiga_data[i].re == 0.0) {
                re = 0.0;
                bim = eiga_data[i].im / beta1_data[i].re;
              } else {
                re = eiga_data[i].re / beta1_data[i].re;
                bim = eiga_data[i].im / beta1_data[i].re;
              }
            } else if (beta1_data[i].re == 0.0) {
              if (eiga_data[i].re == 0.0) {
                re = eiga_data[i].im / beta1_data[i].im;
                bim = 0.0;
              } else if (eiga_data[i].im == 0.0) {
                re = 0.0;
                bim = -(eiga_data[i].re / beta1_data[i].im);
              } else {
                re = eiga_data[i].im / beta1_data[i].im;
                bim = -(eiga_data[i].re / beta1_data[i].im);
              }
            } else {
              brm = std::abs(beta1_data[i].re);
              bim = std::abs(beta1_data[i].im);
              if (brm > bim) {
                bim = beta1_data[i].im / beta1_data[i].re;
                d = beta1_data[i].re + bim * beta1_data[i].im;
                re = (eiga_data[i].re + bim * eiga_data[i].im) / d;
                bim = (eiga_data[i].im - bim * eiga_data[i].re) / d;
              } else if (bim == brm) {
                if (beta1_data[i].re > 0.0) {
                  bim = 0.5;
                } else {
                  bim = -0.5;
                }

                if (beta1_data[i].im > 0.0) {
                  d = 0.5;
                } else {
                  d = -0.5;
                }

                re = (eiga_data[i].re * bim + eiga_data[i].im * d) / brm;
                bim = (eiga_data[i].im * bim - eiga_data[i].re * d) / brm;
              } else {
                bim = beta1_data[i].re / beta1_data[i].im;
                d = beta1_data[i].im + bim * beta1_data[i].re;
                re = (bim * eiga_data[i].re + eiga_data[i].im) / d;
                bim = (bim * eiga_data[i].im - eiga_data[i].re) / d;
              }
            }

            eiga_data[i].re = re;
            eiga_data[i].im = bim;
          }
        }
      }

      for (m = 0; m < companDim; m++) {
        r_data[(m - k2) + 9] = eiga_data[m];
      }

      r_size[0] = (companDim - k2) + 9;
    }
  } else if (1 > 9 - k2) {
    r_size[0] = 0;
  } else {
    r_size[0] = 9 - k2;
  }
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_hypotd(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

//
// Arguments    : float u0
//                float u1
// Return Type  : float
//
static float rt_hypotf(float u0, float u1)
{
  float y;
  float a;
  float b;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0F);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0F);
  } else {
    y = a * 1.41421354F;
  }

  return y;
}

//
// Arguments    : const double A[100]
//                double V[100]
//                double T[100]
// Return Type  : void
//
static void schur(const double A[100], double V[100], double T[100])
{
  double work[10];
  int i;
  int im1n_tmp;
  int ix0;
  int in;
  int alpha1_tmp;
  int k;
  double alpha1;
  int knt;
  int b_i;
  int c_i;
  double tau[9];
  double xnorm;
  double beta1;
  int jy;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int ix;
  int exitg1;
  int i1;
  std::memcpy(&T[0], &A[0], 100U * sizeof(double));
  std::memset(&work[0], 0, 10U * sizeof(double));
  for (i = 0; i < 9; i++) {
    im1n_tmp = i * 10 + 2;
    in = (i + 1) * 10;
    alpha1_tmp = (i + 10 * i) + 1;
    alpha1 = T[alpha1_tmp];
    if (i + 3 < 10) {
      b_i = i + 1;
    } else {
      b_i = 8;
    }

    ix0 = b_i + im1n_tmp;
    tau[i] = 0.0;
    xnorm = xnrm2(8 - i, T, ix0);
    if (xnorm != 0.0) {
      beta1 = rt_hypotd(T[alpha1_tmp], xnorm);
      if (T[alpha1_tmp] >= 0.0) {
        beta1 = -beta1;
      }

      if (std::abs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        c_i = (ix0 - i) + 7;
        do {
          knt++;
          for (k = ix0; k <= c_i; k++) {
            T[k - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

        beta1 = rt_hypotd(alpha1, xnrm2(8 - i, T, ix0));
        if (alpha1 >= 0.0) {
          beta1 = -beta1;
        }

        tau[i] = (beta1 - alpha1) / beta1;
        xnorm = 1.0 / (alpha1 - beta1);
        c_i = (ix0 - i) + 7;
        for (k = ix0; k <= c_i; k++) {
          T[k - 1] *= xnorm;
        }

        for (k = 0; k <= knt; k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        alpha1 = beta1;
      } else {
        tau[i] = (beta1 - T[alpha1_tmp]) / beta1;
        xnorm = 1.0 / (T[alpha1_tmp] - beta1);
        c_i = (ix0 - i) + 7;
        for (k = ix0; k <= c_i; k++) {
          T[k - 1] *= xnorm;
        }

        alpha1 = beta1;
      }
    }

    T[alpha1_tmp] = 1.0;
    jy = (i + im1n_tmp) - 1;
    ix0 = in + 1;
    if (tau[i] != 0.0) {
      lastv = 8 - i;
      b_i = (jy - i) + 8;
      while ((lastv + 1 > 0) && (T[b_i] == 0.0)) {
        lastv--;
        b_i--;
      }

      lastc = 10;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        knt = in + lastc;
        k = knt;
        do {
          exitg1 = 0;
          if (k <= knt + lastv * 10) {
            if (T[k - 1] != 0.0) {
              exitg1 = 1;
            } else {
              k += 10;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = -1;
      lastc = 0;
    }

    if (lastv + 1 > 0) {
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work[0], 0, lastc * sizeof(double));
        }

        ix = jy;
        c_i = (in + 10 * lastv) + 1;
        for (knt = ix0; knt <= c_i; knt += 10) {
          b_i = 0;
          i1 = (knt + lastc) - 1;
          for (k = knt; k <= i1; k++) {
            work[b_i] += T[k - 1] * T[ix];
            b_i++;
          }

          ix++;
        }
      }

      if (-tau[i] != 0.0) {
        knt = in;
        for (ix0 = 0; ix0 <= lastv; ix0++) {
          if (T[jy] != 0.0) {
            xnorm = T[jy] * -tau[i];
            ix = 0;
            c_i = knt + 1;
            i1 = lastc + knt;
            for (b_i = c_i; b_i <= i1; b_i++) {
              T[b_i - 1] += work[ix] * xnorm;
              ix++;
            }
          }

          jy++;
          knt += 10;
        }
      }
    }

    xzlarf(9 - i, 9 - i, i + im1n_tmp, tau[i], T, (i + in) + 2, work);
    T[alpha1_tmp] = alpha1;
  }

  std::memcpy(&V[0], &T[0], 100U * sizeof(double));
  for (ix0 = 8; ix0 >= 0; ix0--) {
    k = (ix0 + 1) * 10;
    for (i = 0; i <= ix0; i++) {
      V[k + i] = 0.0;
    }

    c_i = ix0 + 3;
    for (i = c_i; i < 11; i++) {
      knt = k + i;
      V[knt - 1] = V[knt - 11];
    }
  }

  std::memset(&V[0], 0, 10U * sizeof(double));
  V[0] = 1.0;
  knt = 8;
  std::memset(&work[0], 0, 10U * sizeof(double));
  for (i = 8; i >= 0; i--) {
    b_i = (i + i * 10) + 11;
    if (i + 1 < 9) {
      V[b_i] = 1.0;
      xzlarf(9 - i, 8 - i, b_i + 1, tau[knt], V, b_i + 11, work);
      ix0 = b_i + 2;
      c_i = (b_i - i) + 9;
      for (k = ix0; k <= c_i; k++) {
        V[k - 1] *= -tau[knt];
      }
    }

    V[b_i] = 1.0 - tau[knt];
    for (ix0 = 0; ix0 < i; ix0++) {
      V[(b_i - ix0) - 1] = 0.0;
    }

    knt--;
  }

  eml_dlahqr(T, V);
  knt = 4;
  for (ix0 = 0; ix0 < 7; ix0++) {
    if (knt <= 10) {
      std::memset(&T[(ix0 * 10 + knt) + -1], 0, (11 - knt) * sizeof(double));
    }

    knt++;
  }
}

//
// function [n, xs, ys, zs] = solve_3Q3(c, e)
// c -- 3x10 coefficients matrix
// SOLVE_3Q3 Summary of this function goes here
//    Detailed explanation goes here
// Arguments    : const double c[30]
//                double *n
//                double xs_data[]
//                int xs_size[2]
//                double ys_data[]
//                int ys_size[2]
//                double zs_data[]
//                int zs_size[2]
// Return Type  : void
//
static void solve_3Q3(const double c[30], double *n, double xs_data[], int
                      xs_size[2], double ys_data[], int ys_size[2], double
                      zs_data[], int zs_size[2])
{
  double A[9];
  int i;
  double P[27];
  double a21;
  double b_A[9];
  double M[45];
  int r1;
  int r2;
  double mons[5];
  int r3;
  double maxval;
  double P_prime[27];
  double mons_tmp;
  double b_mons_tmp;
  double c_mons_tmp;
  double d_mons_tmp;
  double e_mons_tmp;
  int rtemp;
  double f_mons_tmp;
  int k;
  double g_mons_tmp;
  int b_k;
  double h_mons_tmp;
  double i_mons_tmp;
  double j_mons_tmp;
  double k_mons_tmp;
  double l_mons_tmp;
  double m_mons_tmp;
  double n_mons_tmp;
  double o_mons_tmp;
  double p_mons_tmp;
  double q_mons_tmp;
  double r_mons_tmp;
  double s_mons_tmp;
  double t_mons_tmp;
  double u_mons_tmp;
  double v_mons_tmp;
  double w_mons_tmp;
  double C[13];
  double y[20];
  double b_C[13];
  double c_C[13];
  creal_T xs_complex_data[8];
  int xs_complex_size[1];
  double b_xs_data[8];
  double c_A[9];

  // 'solve_3Q3:4' xs = zeros(1, 0, 'like', c);
  xs_size[0] = 1;
  xs_size[1] = 0;

  // 'solve_3Q3:5' coder.varsize('xs', [1 13], [0 1]);
  // 'solve_3Q3:6' ys = zeros(1, 0, 'like', c);
  ys_size[0] = 1;
  ys_size[1] = 0;

  // 'solve_3Q3:7' coder.varsize('ys', [1 13], [0 1]);
  // 'solve_3Q3:8' zs = zeros(1, 0, 'like', c);
  zs_size[0] = 1;
  zs_size[1] = 0;

  // 'solve_3Q3:9' coder.varsize('zs', [1 13], [0 1]);
  // 'solve_3Q3:11' A = find_A(c);
  // 'solve_3Q3:47' A = [c(1, 2), c(1, 3), c(1, 6);
  // 'solve_3Q3:48'          c(2, 2), c(2, 3), c(2, 6);
  // 'solve_3Q3:49'          c(3, 2), c(3, 3), c(3, 6)];
  A[0] = c[3];
  A[3] = c[6];
  A[6] = c[15];
  A[1] = c[4];
  A[4] = c[7];
  A[7] = c[16];
  A[2] = c[5];
  A[5] = c[8];
  A[8] = c[17];

  // 'solve_3Q3:12' if rcond(A) < eps
 if(c_rcond(A) < 2.2204460492503131E-16) { //todo:change back
  //if(false){
    // 'solve_3Q3:13' n = 0;
    *n = 0.0;
  } else {
    // 'solve_3Q3:16' P = find_P(c);
    // 'solve_3Q3:53' P = zeros(3, 3, 3, 'like', c);
    // [x^2, x, 1]
    // 'solve_3Q3:54' for i = 1 : 3
    for (i = 0; i < 3; i++) {
      // 'solve_3Q3:55' P(i, 1, :) = [0, -c(i, 4), -c(i, 8)];
      P[i] = 0.0;
      P[i + 9] = -c[i + 9];
      P[i + 18] = -c[i + 21];

      // y
      // 'solve_3Q3:56' P(i, 2, :) = [0, -c(i, 5), -c(i, 9)];
      P[i + 3] = 0.0;
      P[i + 12] = -c[i + 12];
      P[i + 21] = -c[i + 24];

      // z
      // 'solve_3Q3:57' P(i, 3, :) = [-c(i, 1), -c(i, 7), -c(i, 10)];
      P[i + 6] = -c[i];
      P[i + 15] = -c[i + 18];
      P[i + 24] = -c[i + 27];

      // 1
    }

    // 'solve_3Q3:17' P_prime = zeros(3, 3, 3, 'like', c);
    // 'solve_3Q3:18' for i = 1 : 3
    a21 = std::abs(c[4]);
    for (i = 0; i < 3; i++) {
      // 'solve_3Q3:19' P_prime(:, :, i) = A\P(:, :, i);
      std::memcpy(&b_A[0], &A[0], 9U * sizeof(double));
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
      b_A[r2 + 3] -= b_A[r2] * b_A[r1 + 3];
      b_A[r3 + 3] -= b_A[r3] * b_A[r1 + 3];
      b_A[r2 + 6] -= b_A[r2] * b_A[r1 + 6];
      b_A[r3 + 6] -= b_A[r3] * b_A[r1 + 6];
      if (std::abs(b_A[r3 + 3]) > std::abs(b_A[r2 + 3])) {
        rtemp = r2;
        r2 = r3;
        r3 = rtemp;
      }

      b_A[r3 + 3] /= b_A[r2 + 3];
      b_A[r3 + 6] -= b_A[r3 + 3] * b_A[r2 + 6];
      rtemp = r1 + 9 * i;
      k = r2 + 9 * i;
      maxval = P[k] - P[rtemp] * b_A[r2];
      b_k = r3 + 9 * i;
      mons_tmp = b_A[r3 + 3];
      b_mons_tmp = b_A[r3 + 6];
      c_mons_tmp = ((P[b_k] - P[rtemp] * b_A[r3]) - maxval * mons_tmp) /
        b_mons_tmp;
      P_prime[9 * i + 2] = c_mons_tmp;
      d_mons_tmp = b_A[r1 + 6];
      e_mons_tmp = b_A[r2 + 6];
      maxval -= c_mons_tmp * e_mons_tmp;
      f_mons_tmp = b_A[r2 + 3];
      maxval /= f_mons_tmp;
      P_prime[9 * i + 1] = maxval;
      g_mons_tmp = b_A[r1 + 3];
      P_prime[9 * i] = ((P[rtemp] - c_mons_tmp * d_mons_tmp) - maxval *
                        g_mons_tmp) / b_A[r1];
      h_mons_tmp = P[rtemp + 3];
      maxval = P[k + 3] - h_mons_tmp * b_A[r2];
      c_mons_tmp = ((P[b_k + 3] - h_mons_tmp * b_A[r3]) - maxval * mons_tmp) /
        b_mons_tmp;
      P_prime[9 * i + 5] = c_mons_tmp;
      h_mons_tmp -= c_mons_tmp * d_mons_tmp;
      maxval -= c_mons_tmp * e_mons_tmp;
      maxval /= f_mons_tmp;
      P_prime[9 * i + 4] = maxval;
      h_mons_tmp -= maxval * g_mons_tmp;
      h_mons_tmp /= b_A[r1];
      P_prime[9 * i + 3] = h_mons_tmp;
      h_mons_tmp = P[rtemp + 6];
      maxval = P[k + 6] - h_mons_tmp * b_A[r2];
      c_mons_tmp = ((P[b_k + 6] - h_mons_tmp * b_A[r3]) - maxval * mons_tmp) /
        b_mons_tmp;
      P_prime[9 * i + 8] = c_mons_tmp;
      h_mons_tmp -= c_mons_tmp * d_mons_tmp;
      maxval -= c_mons_tmp * e_mons_tmp;
      maxval /= f_mons_tmp;
      P_prime[9 * i + 7] = maxval;
      h_mons_tmp -= maxval * g_mons_tmp;
      h_mons_tmp /= b_A[r1];
      P_prime[9 * i + 6] = h_mons_tmp;
    }

    // 'solve_3Q3:21' M = find_M(P_prime);
    // 'find_M:2' M = zeros(5, 3, 3, 'like', P);
    std::memset(&M[0], 0, 45U * sizeof(double));

    // 'find_M:3' M(:, 1, 1) = [0, 0, P(1, 2, 2)*P(2, 1, 2) - P(3, 3, 1) - P(1, 1, 2)*P(3, 1, 2) + P(3, 1, 2)*(P(1, 1, 2) - P(3, 2, 2)), P(1, 2, 2)*P(2, 1, 3) - P(3, 3, 2) + P(1, 2, 3)*P(2, 1, 2) - P(1, 1, 2)*P(3, 1, 3) - P(1, 1, 3)*P(3, 1, 2) + P(3, 1, 3)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 1, 2)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 2, 3)*P(2, 1, 3) - P(3, 3, 3) - P(1, 1, 3)*P(3, 1, 3) + P(3, 1, 3)*(P(1, 1, 3) - P(3, 2, 3))]; 
    mons[0] = 0.0;
    mons[1] = 0.0;
    maxval = P_prime[9] - P_prime[14];
    mons_tmp = P_prime[12] * P_prime[10];
    mons[2] = ((mons_tmp - P_prime[8]) - P_prime[9] * P_prime[11]) + P_prime[11]
      * maxval;
    b_mons_tmp = P_prime[18] - P_prime[23];
    c_mons_tmp = P_prime[12] * P_prime[19];
    d_mons_tmp = P_prime[21] * P_prime[10];
    mons[3] = (((((c_mons_tmp - P_prime[17]) + d_mons_tmp) - P_prime[9] *
                 P_prime[20]) - P_prime[18] * P_prime[11]) + P_prime[20] *
               maxval) + P_prime[11] * b_mons_tmp;
    e_mons_tmp = P_prime[21] * P_prime[19];
    mons[4] = ((e_mons_tmp - P_prime[26]) - P_prime[18] * P_prime[20]) +
      P_prime[20] * b_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp] = mons[rtemp];
    }

    //  degree = 2
    // 'find_M:5' M(:, 1, 2) = [0, 0, P(1, 3, 1) + P(1, 2, 2)*P(2, 2, 2) - P(1, 2, 2)*P(3, 1, 2) + P(3, 2, 2)*(P(1, 1, 2) - P(3, 2, 2)), P(1, 3, 2) + P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - P(1, 2, 2)*P(3, 1, 3) - P(1, 2, 3)*P(3, 1, 2) + P(3, 2, 3)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 2, 2)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 3, 3) + P(1, 2, 3)*P(2, 2, 3) - P(1, 2, 3)*P(3, 1, 3) + P(3, 2, 3)*(P(1, 1, 3) - P(3, 2, 3))]; 
    mons[0] = 0.0;
    mons[1] = 0.0;
    f_mons_tmp = P_prime[12] * P_prime[13];
    mons[2] = ((P_prime[6] + f_mons_tmp) - P_prime[12] * P_prime[11]) + P_prime
      [14] * maxval;
    g_mons_tmp = P_prime[12] * P_prime[22];
    h_mons_tmp = P_prime[21] * P_prime[13];
    mons[3] = (((((P_prime[15] + g_mons_tmp) + h_mons_tmp) - P_prime[12] *
                 P_prime[20]) - P_prime[21] * P_prime[11]) + P_prime[23] *
               maxval) + P_prime[14] * b_mons_tmp;
    a21 = P_prime[21] * P_prime[22];
    mons[4] = ((P_prime[24] + a21) - P_prime[21] * P_prime[20]) + P_prime[23] *
      b_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 15] = mons[rtemp];
    }

    //  degree = 2
    // 'find_M:7' M(:, 1, 3) = [0, P(1, 2, 2)*P(2, 3, 1) - P(1, 3, 1)*P(3, 1, 2) + P(3, 3, 1)*(P(1, 1, 2) - P(3, 2, 2)), P(1, 2, 2)*P(2, 3, 2) + P(1, 2, 3)*P(2, 3, 1) - P(1, 3, 1)*P(3, 1, 3) - P(1, 3, 2)*P(3, 1, 2) + P(3, 3, 2)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 3, 1)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 2, 2)*P(2, 3, 3) + P(1, 2, 3)*P(2, 3, 2) - P(1, 3, 2)*P(3, 1, 3) - P(1, 3, 3)*P(3, 1, 2) + P(3, 3, 3)*(P(1, 1, 2) - P(3, 2, 2)) + P(3, 3, 2)*(P(1, 1, 3) - P(3, 2, 3)), P(1, 2, 3)*P(2, 3, 3) - P(1, 3, 3)*P(3, 1, 3) + P(3, 3, 3)*(P(1, 1, 3) - P(3, 2, 3))]; 
    mons[0] = 0.0;
    i_mons_tmp = P_prime[12] * P_prime[7];
    mons[1] = (i_mons_tmp - P_prime[6] * P_prime[11]) + P_prime[8] * maxval;
    j_mons_tmp = P_prime[12] * P_prime[16];
    k_mons_tmp = P_prime[21] * P_prime[7];
    mons[2] = ((((j_mons_tmp + k_mons_tmp) - P_prime[6] * P_prime[20]) -
                P_prime[15] * P_prime[11]) + P_prime[17] * maxval) + P_prime[8] *
      b_mons_tmp;
    l_mons_tmp = P_prime[12] * P_prime[25];
    m_mons_tmp = P_prime[21] * P_prime[16];
    mons[3] = ((((l_mons_tmp + m_mons_tmp) - P_prime[15] * P_prime[20]) -
                P_prime[24] * P_prime[11]) + P_prime[26] * maxval) + P_prime[17]
      * b_mons_tmp;
    maxval = P_prime[21] * P_prime[25];
    mons[4] = (maxval - P_prime[24] * P_prime[20]) + P_prime[26] * b_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 30] = mons[rtemp];
    }

    //  degree = 3
    // 'find_M:9' M(:, 2, 1) = [0, 0, P(2, 1, 2)*P(3, 2, 2) - P(1, 1, 2)*P(2, 1, 2) - P(2, 3, 1) - P(3, 1, 2)*(P(2, 2, 2) - P(3, 1, 2)), P(2, 1, 2)*P(3, 2, 3) - P(1, 1, 2)*P(2, 1, 3) - P(1, 1, 3)*P(2, 1, 2) - P(2, 3, 2) + P(2, 1, 3)*P(3, 2, 2) - P(3, 1, 3)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 1, 2)*(P(2, 2, 3) - P(3, 1, 3)), P(2, 1, 3)*P(3, 2, 3) - P(1, 1, 3)*P(2, 1, 3) - P(2, 3, 3) - P(3, 1, 3)*(P(2, 2, 3) - P(3, 1, 3))]; 
    mons[0] = 0.0;
    mons[1] = 0.0;
    b_mons_tmp = P_prime[13] - P_prime[11];
    n_mons_tmp = P_prime[9] * P_prime[10];
    mons[2] = ((P_prime[10] * P_prime[14] - n_mons_tmp) - P_prime[7]) - P_prime
      [11] * b_mons_tmp;
    o_mons_tmp = P_prime[22] - P_prime[20];
    p_mons_tmp = P_prime[9] * P_prime[19];
    q_mons_tmp = P_prime[18] * P_prime[10];
    mons[3] = (((((P_prime[10] * P_prime[23] - p_mons_tmp) - q_mons_tmp) -
                 P_prime[16]) + P_prime[19] * P_prime[14]) - P_prime[20] *
               b_mons_tmp) - P_prime[11] * o_mons_tmp;
    r_mons_tmp = P_prime[18] * P_prime[19];
    mons[4] = ((P_prime[19] * P_prime[23] - r_mons_tmp) - P_prime[25]) -
      P_prime[20] * o_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 5] = mons[rtemp];
    }

    //  degree = 2
    // 'find_M:11' M(:, 2, 2) = [0, 0, P(3, 3, 1) - P(1, 2, 2)*P(2, 1, 2) + P(2, 2, 2)*P(3, 2, 2) - P(3, 2, 2)*(P(2, 2, 2) - P(3, 1, 2)), P(3, 3, 2) - P(1, 2, 2)*P(2, 1, 3) - P(1, 2, 3)*P(2, 1, 2) + P(2, 2, 2)*P(3, 2, 3) + P(2, 2, 3)*P(3, 2, 2) - P(3, 2, 3)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 2, 2)*(P(2, 2, 3) - P(3, 1, 3)), P(3, 3, 3) - P(1, 2, 3)*P(2, 1, 3) + P(2, 2, 3)*P(3, 2, 3) - P(3, 2, 3)*(P(2, 2, 3) - P(3, 1, 3))]; 
    mons[0] = 0.0;
    mons[1] = 0.0;
    mons[2] = ((P_prime[8] - mons_tmp) + P_prime[13] * P_prime[14]) - P_prime[14]
      * b_mons_tmp;
    mons[3] = (((((P_prime[17] - c_mons_tmp) - d_mons_tmp) + P_prime[13] *
                 P_prime[23]) + P_prime[22] * P_prime[14]) - P_prime[23] *
               b_mons_tmp) - P_prime[14] * o_mons_tmp;
    mons[4] = ((P_prime[26] - e_mons_tmp) + P_prime[22] * P_prime[23]) -
      P_prime[23] * o_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 20] = mons[rtemp];
    }

    //  degree = 2
    // 'find_M:13' M(:, 2, 3) = [0, P(2, 3, 1)*P(3, 2, 2) - P(1, 3, 1)*P(2, 1, 2) - P(3, 3, 1)*(P(2, 2, 2) - P(3, 1, 2)), P(2, 3, 1)*P(3, 2, 3) - P(1, 3, 2)*P(2, 1, 2) - P(1, 3, 1)*P(2, 1, 3) + P(2, 3, 2)*P(3, 2, 2) - P(3, 3, 2)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 3, 1)*(P(2, 2, 3) - P(3, 1, 3)), P(2, 3, 2)*P(3, 2, 3) - P(1, 3, 3)*P(2, 1, 2) - P(1, 3, 2)*P(2, 1, 3) + P(2, 3, 3)*P(3, 2, 2) - P(3, 3, 3)*(P(2, 2, 2) - P(3, 1, 2)) - P(3, 3, 2)*(P(2, 2, 3) - P(3, 1, 3)), P(2, 3, 3)*P(3, 2, 3) - P(1, 3, 3)*P(2, 1, 3) - P(3, 3, 3)*(P(2, 2, 3) - P(3, 1, 3))]; 
    mons[0] = 0.0;
    s_mons_tmp = P_prime[6] * P_prime[10];
    mons[1] = (P_prime[7] * P_prime[14] - s_mons_tmp) - P_prime[8] * b_mons_tmp;
    t_mons_tmp = P_prime[6] * P_prime[19];
    u_mons_tmp = P_prime[15] * P_prime[10];
    mons[2] = ((((P_prime[7] * P_prime[23] - u_mons_tmp) - t_mons_tmp) +
                P_prime[16] * P_prime[14]) - P_prime[17] * b_mons_tmp) -
      P_prime[8] * o_mons_tmp;
    v_mons_tmp = P_prime[15] * P_prime[19];
    w_mons_tmp = P_prime[24] * P_prime[10];
    mons[3] = ((((P_prime[16] * P_prime[23] - w_mons_tmp) - v_mons_tmp) +
                P_prime[25] * P_prime[14]) - P_prime[26] * b_mons_tmp) -
      P_prime[17] * o_mons_tmp;
    b_mons_tmp = P_prime[24] * P_prime[19];
    mons[4] = (P_prime[25] * P_prime[23] - b_mons_tmp) - P_prime[26] *
      o_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 35] = mons[rtemp];
    }

    //  degree = 3
    // 'find_M:15' M(:, 3, 1) = [0, 2*P(3, 1, 2)*P(3, 3, 1) - P(1, 3, 1)*P(2, 1, 2) - P(1, 1, 2)*P(2, 3, 1) - P(1, 1, 2)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 1, 2)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(3, 1, 2)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 1, 2)*P(3, 3, 2) - P(1, 1, 2)*P(2, 3, 2) - P(1, 1, 3)*P(2, 3, 1) - P(1, 3, 1)*P(2, 1, 3) - P(1, 3, 2)*P(2, 1, 2) - P(3, 1, 2)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 1, 3)*P(3, 3, 1) - P(1, 1, 3)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 1, 3)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(1, 1, 2)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 1, 2)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 1, 3)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 1, 2)*P(3, 3, 3) - P(1, 1, 2)*P(2, 3, 3) - P(1, 1, 3)*P(2, 3, 2) - P(1, 3, 2)*P(2, 1, 3) - P(1, 3, 3)*P(2, 1, 2) - P(3, 1, 3)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 1, 3)*P(3, 3, 2) - P(1, 1, 2)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 1, 2)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 1, 3)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 1, 3)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 1, 2)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)), 2*P(3, 1, 3)*P(3, 3, 3) - P(1, 3, 3)*P(2, 1, 3) - P(1, 1, 3)*P(2, 3, 3) - P(1, 1, 3)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 1, 3)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(3, 1, 3)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3))]; 
    mons[0] = 0.0;
    n_mons_tmp -= P_prime[11] * P_prime[11];
    f_mons_tmp -= P_prime[14] * P_prime[14];
    mons_tmp = (P_prime[9] * P_prime[13] + mons_tmp) - 2.0 * P_prime[11] *
      P_prime[14];
    mons[1] = ((((2.0 * P_prime[11] * P_prime[8] - s_mons_tmp) - P_prime[9] *
                 P_prime[7]) - P_prime[9] * n_mons_tmp) - P_prime[10] *
               f_mons_tmp) - P_prime[11] * mons_tmp;
    c_mons_tmp = ((((P_prime[9] * P_prime[22] + P_prime[18] * P_prime[13]) +
                    c_mons_tmp) + d_mons_tmp) - 2.0 * P_prime[11] * P_prime[23])
      - 2.0 * P_prime[20] * P_prime[14];
    d_mons_tmp = (p_mons_tmp + q_mons_tmp) - 2.0 * P_prime[11] * P_prime[20];
    g_mons_tmp = (g_mons_tmp + h_mons_tmp) - 2.0 * P_prime[14] * P_prime[23];
    mons[2] = ((((((((((2.0 * P_prime[11] * P_prime[17] - P_prime[9] * P_prime
                        [16]) - P_prime[18] * P_prime[7]) - t_mons_tmp) -
                     u_mons_tmp) - P_prime[11] * c_mons_tmp) + 2.0 * P_prime[20]
                   * P_prime[8]) - P_prime[18] * n_mons_tmp) - P_prime[19] *
                 f_mons_tmp) - P_prime[9] * d_mons_tmp) - P_prime[10] *
               g_mons_tmp) - P_prime[20] * mons_tmp;
    h_mons_tmp = r_mons_tmp - P_prime[20] * P_prime[20];
    a21 -= P_prime[23] * P_prime[23];
    e_mons_tmp = (P_prime[18] * P_prime[22] + e_mons_tmp) - 2.0 * P_prime[20] *
      P_prime[23];
    mons[3] = ((((((((((2.0 * P_prime[11] * P_prime[26] - P_prime[9] * P_prime
                        [25]) - P_prime[18] * P_prime[16]) - v_mons_tmp) -
                     w_mons_tmp) - P_prime[20] * c_mons_tmp) + 2.0 * P_prime[20]
                   * P_prime[17]) - P_prime[9] * h_mons_tmp) - P_prime[10] * a21)
                - P_prime[18] * d_mons_tmp) - P_prime[19] * g_mons_tmp) -
      P_prime[11] * e_mons_tmp;
    mons[4] = ((((2.0 * P_prime[20] * P_prime[26] - b_mons_tmp) - P_prime[18] *
                 P_prime[25]) - P_prime[18] * h_mons_tmp) - P_prime[19] * a21) -
      P_prime[20] * e_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 10] = mons[rtemp];
    }

    //  degree = 3
    // 'find_M:17' M(:, 3, 2) = [0, 2*P(3, 2, 2)*P(3, 3, 1) - P(1, 3, 1)*P(2, 2, 2) - P(1, 2, 2)*P(2, 3, 1) - P(1, 2, 2)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 2, 2)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(3, 2, 2)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 2, 2)*P(3, 3, 2) - P(1, 2, 2)*P(2, 3, 2) - P(1, 2, 3)*P(2, 3, 1) - P(1, 3, 1)*P(2, 2, 3) - P(1, 3, 2)*P(2, 2, 2) - P(3, 2, 2)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 2, 3)*P(3, 3, 1) - P(1, 2, 3)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 2, 3)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(1, 2, 2)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 2, 2)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 2, 3)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 2, 2)*P(3, 3, 3) - P(1, 2, 2)*P(2, 3, 3) - P(1, 2, 3)*P(2, 3, 2) - P(1, 3, 2)*P(2, 2, 3) - P(1, 3, 3)*P(2, 2, 2) - P(3, 2, 3)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) + 2*P(3, 2, 3)*P(3, 3, 2) - P(1, 2, 2)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 2, 2)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 2, 3)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 2, 3)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 2, 2)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)), 2*P(3, 2, 3)*P(3, 3, 3) - P(1, 3, 3)*P(2, 2, 3) - P(1, 2, 3)*P(2, 3, 3) - P(1, 2, 3)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 2, 3)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(3, 2, 3)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3))]; 
    mons[0] = 0.0;
    mons[1] = ((((2.0 * P_prime[14] * P_prime[8] - P_prime[6] * P_prime[13]) -
                 i_mons_tmp) - P_prime[12] * n_mons_tmp) - P_prime[13] *
               f_mons_tmp) - P_prime[14] * mons_tmp;
    mons[2] = ((((((((((2.0 * P_prime[14] * P_prime[17] - j_mons_tmp) -
                       k_mons_tmp) - P_prime[6] * P_prime[22]) - P_prime[15] *
                     P_prime[13]) - P_prime[14] * c_mons_tmp) + 2.0 * P_prime[23]
                   * P_prime[8]) - P_prime[21] * n_mons_tmp) - P_prime[22] *
                 f_mons_tmp) - P_prime[12] * d_mons_tmp) - P_prime[13] *
               g_mons_tmp) - P_prime[23] * mons_tmp;
    mons[3] = ((((((((((2.0 * P_prime[14] * P_prime[26] - l_mons_tmp) -
                       m_mons_tmp) - P_prime[15] * P_prime[22]) - P_prime[24] *
                     P_prime[13]) - P_prime[23] * c_mons_tmp) + 2.0 * P_prime[23]
                   * P_prime[17]) - P_prime[12] * h_mons_tmp) - P_prime[13] *
                 a21) - P_prime[21] * d_mons_tmp) - P_prime[22] * g_mons_tmp) -
      P_prime[14] * e_mons_tmp;
    mons[4] = ((((2.0 * P_prime[23] * P_prime[26] - P_prime[24] * P_prime[22]) -
                 maxval) - P_prime[21] * h_mons_tmp) - P_prime[22] * a21) -
      P_prime[23] * e_mons_tmp;
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 25] = mons[rtemp];
    }

    //  degree = 3
    // 'find_M:19' M(:, 3, 3) = [P(3, 3, 1)^2 - P(1, 3, 1)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 3, 1)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(3, 3, 1)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)) - P(1, 3, 1)*P(2, 3, 1), 2*P(3, 3, 1)*P(3, 3, 2) - P(1, 3, 1)*P(2, 3, 2) - P(1, 3, 2)*P(2, 3, 1) - P(3, 3, 1)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) - P(1, 3, 2)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(2, 3, 2)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(1, 3, 1)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 3, 1)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 3, 2)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)), 2*P(3, 3, 1)*P(3, 3, 3) - P(1, 3, 1)*P(2, 3, 3) - P(1, 3, 2)*P(2, 3, 2) - P(1, 3, 3)*P(2, 3, 1) - P(3, 3, 2)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) - P(1, 3, 3)*(P(1, 1, 2)*P(2, 1, 2) - P(3, 1, 2)^2) - P(1, 3, 1)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 3, 3)*(P(1, 2, 2)*P(2, 2, 2) - P(3, 2, 2)^2) - P(2, 3, 1)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 3, 2)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 3, 2)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 3, 3)*(P(1, 1, 2)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 2)) - P(3, 3, 1)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)) + P(3, 3, 2)^2, 2*P(3, 3, 2)*P(3, 3, 3) - P(1, 3, 2)*P(2, 3, 3) - P(1, 3, 3)*P(2, 3, 2) - P(3, 3, 3)*(P(1, 1, 2)*P(2, 2, 3) + P(1, 1, 3)*P(2, 2, 2) + P(1, 2, 2)*P(2, 1, 3) + P(1, 2, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 2, 3) - 2*P(3, 1, 3)*P(3, 2, 2)) - P(1, 3, 2)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 3, 2)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(1, 3, 3)*(P(1, 1, 2)*P(2, 1, 3) + P(1, 1, 3)*P(2, 1, 2) - 2*P(3, 1, 2)*P(3, 1, 3)) - P(2, 3, 3)*(P(1, 2, 2)*P(2, 2, 3) + P(1, 2, 3)*P(2, 2, 2) - 2*P(3, 2, 2)*P(3, 2, 3)) - P(3, 3, 2)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)), P(3, 3, 3)^2 - P(1, 3, 3)*(P(1, 1, 3)*P(2, 1, 3) - P(3, 1, 3)^2) - P(2, 3, 3)*(P(1, 2, 3)*P(2, 2, 3) - P(3, 2, 3)^2) - P(3, 3, 3)*(P(1, 1, 3)*P(2, 2, 3) + P(1, 2, 3)*P(2, 1, 3) - 2*P(3, 1, 3)*P(3, 2, 3)) - P(1, 3, 3)*P(2, 3, 3)]; 
    mons[0] = (((P_prime[8] * P_prime[8] - P_prime[6] * n_mons_tmp) - P_prime[7]
                * f_mons_tmp) - P_prime[8] * mons_tmp) - P_prime[6] * P_prime[7];
    mons[1] = (((((((2.0 * P_prime[8] * P_prime[17] - P_prime[6] * P_prime[16])
                    - P_prime[15] * P_prime[7]) - P_prime[8] * c_mons_tmp) -
                  P_prime[15] * n_mons_tmp) - P_prime[16] * f_mons_tmp) -
                P_prime[6] * d_mons_tmp) - P_prime[7] * g_mons_tmp) - P_prime[17]
      * mons_tmp;
    mons[2] = ((((((((((((2.0 * P_prime[8] * P_prime[26] - P_prime[6] * P_prime
                          [25]) - P_prime[15] * P_prime[16]) - P_prime[24] *
                        P_prime[7]) - P_prime[17] * c_mons_tmp) - P_prime[24] *
                      n_mons_tmp) - P_prime[6] * h_mons_tmp) - P_prime[25] *
                    f_mons_tmp) - P_prime[7] * a21) - P_prime[15] * d_mons_tmp)
                 - P_prime[16] * g_mons_tmp) - P_prime[26] * mons_tmp) -
               P_prime[8] * e_mons_tmp) + P_prime[17] * P_prime[17];
    mons[3] = (((((((2.0 * P_prime[17] * P_prime[26] - P_prime[15] * P_prime[25])
                    - P_prime[24] * P_prime[16]) - P_prime[26] * c_mons_tmp) -
                  P_prime[15] * h_mons_tmp) - P_prime[16] * a21) - P_prime[24] *
                d_mons_tmp) - P_prime[25] * g_mons_tmp) - P_prime[17] *
      e_mons_tmp;
    mons[4] = (((P_prime[26] * P_prime[26] - P_prime[24] * h_mons_tmp) -
                P_prime[25] * a21) - P_prime[26] * e_mons_tmp) - P_prime[24] *
      P_prime[25];
    for (rtemp = 0; rtemp < 5; rtemp++) {
      M[rtemp + 40] = mons[rtemp];
    }

    //  degree = 4
    // 'solve_3Q3:22' pol = find_det_M(M);
    // 'solve_3Q3:62' d_ = conv(M(:, 1, 1), find_det2(M(:, 2:3, 2:3))) - ...
    // 'solve_3Q3:63'          conv(M(:, 1, 2), find_det2(cat(3, M(:, 2:3, 1), M(:, 2:3, 3)))) + ... 
    // 'solve_3Q3:64'          conv(M(:, 1, 3), find_det2(M(:, 2:3, 1:2)));
    // 'solve_3Q3:70' d = conv(M(:, 1, 1), M(:, 2, 2)) - conv(M(:, 1, 2), M(:, 2, 1)); 
    std::memset(&A[0], 0, 9U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 5; b_k++) {
        rtemp = k + b_k;
        A[rtemp] += M[k + 40] * M[b_k + 20];
      }
    }

    std::memset(&b_A[0], 0, 9U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 5; b_k++) {
        rtemp = k + b_k;
        b_A[rtemp] += M[k + 25] * M[b_k + 35];
      }
    }

    for (rtemp = 0; rtemp < 9; rtemp++) {
      A[rtemp] -= b_A[rtemp];
    }

    std::memset(&C[0], 0, 13U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 9; b_k++) {
        rtemp = k + b_k;
        C[rtemp] += M[k] * A[b_k];
      }
    }

    rtemp = -1;
    for (r1 = 0; r1 < 10; r1++) {
      rtemp++;
      y[rtemp] = M[r1 % 5 + 5 * (r1 / 5 + 1)];
    }

    for (r1 = 0; r1 < 10; r1++) {
      rtemp++;
      y[rtemp] = M[(r1 % 5 + 5 * (r1 / 5 + 1)) + 30];
    }

    // 'solve_3Q3:70' d = conv(M(:, 1, 1), M(:, 2, 2)) - conv(M(:, 1, 2), M(:, 2, 1)); 
    std::memset(&A[0], 0, 9U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 5; b_k++) {
        rtemp = k + b_k;
        A[rtemp] += y[k + 15] * y[b_k];
      }
    }

    std::memset(&b_A[0], 0, 9U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 5; b_k++) {
        rtemp = k + b_k;
        b_A[rtemp] += y[k + 5] * y[b_k + 10];
      }
    }

    for (rtemp = 0; rtemp < 9; rtemp++) {
      A[rtemp] -= b_A[rtemp];
    }

    std::memset(&b_C[0], 0, 13U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 9; b_k++) {
        rtemp = k + b_k;
        b_C[rtemp] += M[k + 15] * A[b_k];
      }
    }

    // 'solve_3Q3:70' d = conv(M(:, 1, 1), M(:, 2, 2)) - conv(M(:, 1, 2), M(:, 2, 1)); 
    std::memset(&A[0], 0, 9U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 5; b_k++) {
        rtemp = k + b_k;
        A[rtemp] += M[k + 25] * M[b_k + 5];
      }
    }

    std::memset(&b_A[0], 0, 9U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 5; b_k++) {
        rtemp = k + b_k;
        b_A[rtemp] += M[k + 10] * M[b_k + 20];
      }
    }

    for (rtemp = 0; rtemp < 9; rtemp++) {
      A[rtemp] -= b_A[rtemp];
    }

    std::memset(&c_C[0], 0, 13U * sizeof(double));
    for (k = 0; k < 5; k++) {
      for (b_k = 0; b_k < 9; b_k++) {
        rtemp = k + b_k;
        c_C[rtemp] += M[k + 30] * A[b_k];
      }
    }

    //  Due to the structure of M, d is 8-degree polynomial
    // 'solve_3Q3:66' d = d_(end-8:end);
    // 'solve_3Q3:23' assert(numel(pol)==9);
    // 'solve_3Q3:24' if ~isfinite(pol)
    // 'solve_3Q3:29' xs_complex = roots(pol');
    for (rtemp = 0; rtemp < 13; rtemp++) {
      C[rtemp] = (C[rtemp] - b_C[rtemp]) + c_C[rtemp];
    }

    roots(*(double (*)[9])&C[4], xs_complex_data, xs_complex_size);

    // 'solve_3Q3:30' xs = zeros(1, length(xs_complex), 'like', c);
    k = xs_complex_size[0];
    if (0 <= k - 1) {
      std::memset(&xs_data[0], 0, k * sizeof(double));
    }

    // 'solve_3Q3:31' n = 0;
    *n = 0.0;

    // 'solve_3Q3:33' for i = 1 : length(xs_complex)
    rtemp = xs_complex_size[0];
    for (i = 0; i < rtemp; i++) {
      // 'solve_3Q3:34' n = n + 1;
      (*n)++;

      // 'solve_3Q3:35' xs(n) = real(xs_complex(i));
      xs_data[static_cast<int>(*n) - 1] = xs_complex_data[i].re;
    }

    // 'solve_3Q3:38' xs = xs(1:n);
    if (1.0 > *n) {
      k = 0;
    } else {
      k = static_cast<int>(*n);
    }

    if (0 <= k - 1) {
      std::memcpy(&b_xs_data[0], &xs_data[0], k * sizeof(double));
      std::memcpy(&C[0], &xs_data[0], k * sizeof(double));
    }

    xs_size[0] = 1;
    xs_size[1] = k;
    if (0 <= k - 1) {
      std::memcpy(&xs_data[0], &C[0], k * sizeof(double));
    }

    // 'solve_3Q3:39' ys = zeros(1, n, 'like', c);
    ys_size[0] = 1;
    k = static_cast<int>(*n);
    ys_size[1] = k;

    // 'solve_3Q3:40' zs = zeros(1, n, 'like', c);
    zs_size[0] = 1;
    zs_size[1] = k;
    if (0 <= k - 1) {
      std::memset(&ys_data[0], 0, k * sizeof(double));
      std::memset(&zs_data[0], 0, k * sizeof(double));
    }

    // 'solve_3Q3:41' for i = 1 : n
    if (0 <= k - 1) {
      mons[4] = 1.0;
    }

    for (i = 0; i < k; i++) {
      // 'solve_3Q3:42' [ys(i), zs(i)] = find_yz(M, xs(i));
      // 'find_yz:2' M_subs = zeros(3, 'like', x);
      // 'find_yz:3' mons = [x^4, x^3, x^2, x, 1];
      mons[0] = pow(b_xs_data[i], 4.0);
      mons[1] = pow(b_xs_data[i], 3.0);
      mons[2] = b_xs_data[i] * b_xs_data[i];
      mons[3] = xs_data[i];

      // 'find_yz:4' for i = 1 : 3
      // 'find_yz:9' [Q, ~] = qr(M_subs');
      for (b_k = 0; b_k < 3; b_k++) {
        // 'find_yz:5' for j = 1 : 3
        for (r1 = 0; r1 < 3; r1++) {
          // 'find_yz:6' M_subs(i, j) = mons * M(:, i, j);
          maxval = 0.0;
          for (rtemp = 0; rtemp < 5; rtemp++) {
            maxval += mons[rtemp] * M[(rtemp + 5 * b_k) + 15 * r1];
          }

          c_A[r1 + 3 * b_k] = maxval;
        }
      }

      d_qr(c_A, A, b_A);

      // 'find_yz:9' ~
      // 'find_yz:10' y = Q(1, 3) / Q(3, 3);
      ys_data[i] = A[6] / A[8];

      // 'find_yz:11' z = Q(2, 3) / Q(3, 3);
      zs_data[i] = A[7] / A[8];
    }
  }
}

//
// Arguments    : int n
//                double a
//                const double x[2]
//                double y[8]
//                int iy0
// Return Type  : void
//
static void xaxpy(int n, double a, const double x[2], double y[8], int iy0)
{
  int iy;
  if ((n >= 1) && (a != 0.0)) {
    iy = iy0 - 1;
    y[iy] += a * x[1];
  }
}

//
// Arguments    : double *a
//                double *b
//                double *c
//                double *d
//                double *rt1r
//                double *rt1i
//                double *rt2r
//                double *rt2i
//                double *cs
//                double *sn
// Return Type  : void
//
static void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn)
{
  double tau;
  double p;
  double bcmax;
  double bcmis;
  double scale;
  int b_b;
  int b_c;
  double z;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    bcmax = *d;
    *d = *a;
    *a = bcmax;
    *b = -*c;
    *c = 0.0;
  } else {
    tau = *a - *d;
    if ((tau == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      p = 0.5 * tau;
      bcmis = std::abs(*b);
      scale = std::abs(*c);
      if (bcmis > scale) {
        bcmax = bcmis;
      } else {
        bcmax = scale;
      }

      if (bcmis < scale) {
        scale = bcmis;
      }

      if (*b >= 0.0) {
        b_b = 1;
      } else {
        b_b = -1;
      }

      if (*c >= 0.0) {
        b_c = 1;
      } else {
        b_c = -1;
      }

      bcmis = scale * static_cast<double>(b_b) * static_cast<double>(b_c);
      scale = std::abs(p);
      if (scale <= bcmax) {
        scale = bcmax;
      }

      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817841970012523E-16) {
        *a = std::sqrt(scale) * std::sqrt(z);
        if (p < 0.0) {
          *a = -*a;
        }

        z = p + *a;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        tau = rt_hypotd(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0;
      } else {
        bcmis = *b + *c;
        tau = rt_hypotd(bcmis, tau);
        *cs = std::sqrt(0.5 * (std::abs(bcmis) / tau + 1.0));
        if (bcmis >= 0.0) {
          b_b = 1;
        } else {
          b_b = -1;
        }

        *sn = -(p / (tau * *cs)) * static_cast<double>(b_b);
        bcmax = *a * *cs + *b * *sn;
        scale = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        bcmis = -*c * *sn + *d * *cs;
        *b = scale * *cs + bcmis * *sn;
        *c = -bcmax * *sn + z * *cs;
        bcmax = 0.5 * ((bcmax * *cs + z * *sn) + (-scale * *sn + bcmis * *cs));
        *a = bcmax;
        *d = bcmax;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              bcmis = std::sqrt(std::abs(*b));
              z = std::sqrt(std::abs(*c));
              *a = bcmis * z;
              if (*c >= 0.0) {
                p = *a;
              } else {
                p = -*a;
              }

              tau = 1.0 / std::sqrt(std::abs(*b + *c));
              *a = bcmax + p;
              *d = bcmax - p;
              *b -= *c;
              *c = 0.0;
              scale = bcmis * tau;
              bcmis = z * tau;
              bcmax = *cs * scale - *sn * bcmis;
              *sn = *cs * bcmis + *sn * scale;
              *cs = bcmax;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            bcmax = *cs;
            *cs = -*sn;
            *sn = bcmax;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = std::sqrt(std::abs(*b)) * std::sqrt(std::abs(*c));
    *rt2i = -*rt1i;
  }
}

//
// Arguments    : double x1
//                double x2
//                double x3
// Return Type  : double
//
static double xdlapy3(double x1, double x2, double x3)
{
  double y;
  double a;
  double b;
  double c;
  a = std::abs(x1);
  b = std::abs(x2);
  c = std::abs(x3);
  if (a > b) {
    y = a;
  } else {
    y = b;
  }

  if (c > y) {
    y = c;
  }

  if (y > 0.0) {
    a /= y;
    b /= y;
    c /= y;
    y *= std::sqrt((a * a + c * c) + b * b);
  } else {
    y = (a + b) + c;
  }

  return y;
}

//
// Arguments    : creal_T a_data[]
//                const int a_size[2]
// Return Type  : void
//
static void xgehrd(creal_T a_data[], const int a_size[2])
{
  int n;
  int i;
  creal_T work_data[8];
  int b_i;
  int c_i;
  int in;
  int alpha1_tmp;
  creal_T alpha1;
  int n_tmp_tmp;
  int n_tmp;
  creal_T tau_data[7];
  double c_re;
  double beta1;
  int iv0_tmp;
  int jA;
  int lastv;
  int knt;
  int lastc;
  double c_im;
  int i1;
  int rowleft;
  boolean_T exitg2;
  int ix;
  int ia;
  int exitg1;
  int i2;
  double temp_im;
  creal_T b_alpha1;
  double a_im;
  n = a_size[0];
  i = static_cast<signed char>(a_size[0]);
  if (0 <= i - 1) {
    std::memset(&work_data[0], 0, i * sizeof(creal_T));
  }

  b_i = a_size[0];
  for (c_i = 0; c_i <= b_i - 2; c_i++) {
    in = (c_i + 1) * n;
    alpha1_tmp = (c_i + a_size[0] * c_i) + 1;
    alpha1 = a_data[alpha1_tmp];
    i = c_i + 3;
    if (i >= n) {
      i = n;
    }

    i += c_i * n;
    n_tmp_tmp = n - c_i;
    n_tmp = n_tmp_tmp - 3;
    tau_data[c_i].re = 0.0;
    tau_data[c_i].im = 0.0;
    if (n_tmp + 2 > 0) {
      c_re = h_xnrm2(n_tmp + 1, a_data, i);
      if ((c_re != 0.0) || (a_data[alpha1_tmp].im != 0.0)) {
        beta1 = xdlapy3(a_data[alpha1_tmp].re, a_data[alpha1_tmp].im, c_re);
        if (a_data[alpha1_tmp].re >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i1 = i + n_tmp;
          do {
            knt++;
            for (rowleft = i; rowleft <= i1; rowleft++) {
              c_re = 9.9792015476736E+291 * a_data[rowleft - 1].im;
              a_data[rowleft - 1].re *= 9.9792015476736E+291;
              a_data[rowleft - 1].im = c_re;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1.re *= 9.9792015476736E+291;
            alpha1.im *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = xdlapy3(alpha1.re, alpha1.im, h_xnrm2(n_tmp + 1, a_data, i));
          if (alpha1.re >= 0.0) {
            beta1 = -beta1;
          }

          c_re = beta1 - alpha1.re;
          if (0.0 - alpha1.im == 0.0) {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = 0.0;
          } else if (c_re == 0.0) {
            tau_data[c_i].re = 0.0;
            tau_data[c_i].im = (0.0 - alpha1.im) / beta1;
          } else {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = (0.0 - alpha1.im) / beta1;
          }

          b_alpha1.re = alpha1.re - beta1;
          b_alpha1.im = alpha1.im;
          alpha1 = recip(b_alpha1);
          for (rowleft = i; rowleft <= i1; rowleft++) {
            c_re = alpha1.re * a_data[rowleft - 1].im + alpha1.im *
              a_data[rowleft - 1].re;
            a_data[rowleft - 1].re = alpha1.re * a_data[rowleft - 1].re -
              alpha1.im * a_data[rowleft - 1].im;
            a_data[rowleft - 1].im = c_re;
          }

          for (rowleft = 0; rowleft <= knt; rowleft++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        } else {
          c_re = beta1 - a_data[alpha1_tmp].re;
          c_im = 0.0 - a_data[alpha1_tmp].im;
          if (c_im == 0.0) {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = 0.0;
          } else if (c_re == 0.0) {
            tau_data[c_i].re = 0.0;
            tau_data[c_i].im = c_im / beta1;
          } else {
            tau_data[c_i].re = c_re / beta1;
            tau_data[c_i].im = c_im / beta1;
          }

          alpha1.re = a_data[alpha1_tmp].re - beta1;
          alpha1.im = a_data[alpha1_tmp].im;
          alpha1 = recip(alpha1);
          i1 = i + n_tmp;
          for (rowleft = i; rowleft <= i1; rowleft++) {
            c_re = alpha1.re * a_data[rowleft - 1].im + alpha1.im *
              a_data[rowleft - 1].re;
            a_data[rowleft - 1].re = alpha1.re * a_data[rowleft - 1].re -
              alpha1.im * a_data[rowleft - 1].im;
            a_data[rowleft - 1].im = c_re;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        }
      }
    }

    a_data[alpha1_tmp].re = 1.0;
    a_data[alpha1_tmp].im = 0.0;
    iv0_tmp = (c_i + c_i * n) + 1;
    jA = in + 1;
    if ((tau_data[c_i].re != 0.0) || (tau_data[c_i].im != 0.0)) {
      lastv = n_tmp + 1;
      i = iv0_tmp + n_tmp;
      while ((lastv + 1 > 0) && ((a_data[i + 1].re == 0.0) && (a_data[i + 1].im ==
               0.0))) {
        lastv--;
        i--;
      }

      lastc = n;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        rowleft = in + lastc;
        ia = rowleft;
        do {
          exitg1 = 0;
          if ((n > 0) && (ia <= rowleft + lastv * n)) {
            if ((a_data[ia - 1].re != 0.0) || (a_data[ia - 1].im != 0.0)) {
              exitg1 = 1;
            } else {
              ia += n;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = -1;
      lastc = 0;
    }

    if (lastv + 1 > 0) {
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          std::memset(&work_data[0], 0, lastc * sizeof(creal_T));
        }

        ix = iv0_tmp;
        i1 = (in + n * lastv) + 1;
        for (knt = jA; n < 0 ? knt >= i1 : knt <= i1; knt += n) {
          c_re = a_data[ix].re;
          c_im = a_data[ix].im;
          i = 0;
          i2 = (knt + lastc) - 1;
          for (ia = knt; ia <= i2; ia++) {
            work_data[i].re += a_data[ia - 1].re * c_re - a_data[ia - 1].im *
              c_im;
            work_data[i].im += a_data[ia - 1].re * c_im + a_data[ia - 1].im *
              c_re;
            i++;
          }

          ix++;
        }
      }

      c_re = tau_data[c_i].re;
      c_im = tau_data[c_i].im;
      if ((-tau_data[c_i].re != 0.0) || (-tau_data[c_i].im != 0.0)) {
        jA = in;
        knt = iv0_tmp;
        for (i = 0; i <= lastv; i++) {
          if ((a_data[knt].re != 0.0) || (a_data[knt].im != 0.0)) {
            beta1 = a_data[knt].re * -c_re + a_data[knt].im * -c_im;
            temp_im = a_data[knt].re * -c_im - a_data[knt].im * -c_re;
            ix = 0;
            i1 = jA + 1;
            i2 = lastc + jA;
            for (rowleft = i1; rowleft <= i2; rowleft++) {
              a_data[rowleft - 1].re += work_data[ix].re * beta1 - work_data[ix]
                .im * temp_im;
              a_data[rowleft - 1].im += work_data[ix].re * temp_im +
                work_data[ix].im * beta1;
              ix++;
            }
          }

          knt++;
          jA += n;
        }
      }
    }

    jA = (c_i + in) + 2;
    if ((tau_data[c_i].re != 0.0) || (-tau_data[c_i].im != 0.0)) {
      lastv = n_tmp + 2;
      i = iv0_tmp + n_tmp;
      while ((lastv > 0) && ((a_data[i + 1].re == 0.0) && (a_data[i + 1].im ==
               0.0))) {
        lastv--;
        i--;
      }

      lastc = n_tmp_tmp - 2;
      exitg2 = false;
      while ((!exitg2) && (lastc + 1 > 0)) {
        i = jA + lastc * n;
        ia = i;
        do {
          exitg1 = 0;
          if (ia <= (i + lastv) - 1) {
            if ((a_data[ia - 1].re != 0.0) || (a_data[ia - 1].im != 0.0)) {
              exitg1 = 1;
            } else {
              ia++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = -1;
    }

    if (lastv > 0) {
      if (lastc + 1 != 0) {
        if (0 <= lastc) {
          std::memset(&work_data[0], 0, (lastc + 1) * sizeof(creal_T));
        }

        i = 0;
        i1 = jA + n * lastc;
        for (knt = jA; n < 0 ? knt >= i1 : knt <= i1; knt += n) {
          ix = iv0_tmp;
          c_re = 0.0;
          c_im = 0.0;
          i2 = (knt + lastv) - 1;
          for (ia = knt; ia <= i2; ia++) {
            c_re += a_data[ia - 1].re * a_data[ix].re + a_data[ia - 1].im *
              a_data[ix].im;
            c_im += a_data[ia - 1].re * a_data[ix].im - a_data[ia - 1].im *
              a_data[ix].re;
            ix++;
          }

          work_data[i].re += c_re;
          work_data[i].im += c_im;
          i++;
        }
      }

      c_re = tau_data[c_i].re;
      c_im = tau_data[c_i].im;
      if ((-tau_data[c_i].re != 0.0) || (tau_data[c_i].im != 0.0)) {
        knt = 0;
        for (i = 0; i <= lastc; i++) {
          if ((work_data[knt].re != 0.0) || (work_data[knt].im != 0.0)) {
            beta1 = work_data[knt].re * -c_re + work_data[knt].im * c_im;
            temp_im = work_data[knt].re * c_im - work_data[knt].im * -c_re;
            ix = iv0_tmp;
            i1 = lastv + jA;
            for (rowleft = jA; rowleft < i1; rowleft++) {
              a_im = a_data[ix].re * temp_im + a_data[ix].im * beta1;
              a_data[rowleft - 1].re += a_data[ix].re * beta1 - a_data[ix].im *
                temp_im;
              a_data[rowleft - 1].im += a_im;
              ix++;
            }
          }

          knt++;
          jA += n;
        }
      }
    }

    a_data[alpha1_tmp] = alpha1;
  }
}

//
// Arguments    : int m
//                int n
//                double alpha1
//                int ix0
//                const double y[4]
//                double A[12]
//                int ia0
// Return Type  : void
//
static void xgerc(int m, int n, double alpha1, int ix0, const double y[4],
                  double A[12], int ia0)
{
  int jA;
  int jy;
  int j;
  double temp;
  int ix;
  int i;
  int ijA;
  if (alpha1 != 0.0) {
    jA = ia0;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

//
// Arguments    : int n
//                const double x[100]
//                int ix0
// Return Type  : double
//
static double xnrm2(int n, const double x[100], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

//
// Arguments    : int n
//                double x[100]
//                int ix0
//                int iy0
//                double c
//                double s
// Return Type  : void
//
static void xrot(int n, double x[100], int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  if (n >= 1) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      temp = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = temp;
      iy++;
      ix++;
    }
  }
}

//
// Arguments    : double *a
//                double *b
//                double *c
//                double *s
// Return Type  : void
//
static void xrotg(double *a, double *b, double *c, double *s)
{
  double roe;
  double absa;
  double absb;
  double scale;
  double ads;
  double bds;
  roe = *b;
  absa = std::abs(*a);
  absb = std::abs(*b);
  if (absa > absb) {
    roe = *a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    ads = absa / scale;
    bds = absb / scale;
    scale *= std::sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }

    *c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }

    *a = scale;
  }
}

//
// Arguments    : const creal_T A_data[]
//                const int A_size[2]
//                int *info
//                creal_T alpha1_data[]
//                int alpha1_size[1]
//                creal_T beta1_data[]
//                int beta1_size[1]
// Return Type  : void
//
static void xzgeev(const creal_T A_data[], const int A_size[2], int *info,
                   creal_T alpha1_data[], int alpha1_size[1], creal_T
                   beta1_data[], int beta1_size[1])
{
  int At_size[2];
  int ii;
  creal_T At_data[64];
  double anrm;
  int jcolp1;
  boolean_T ilascl;
  double absxk;
  double anrmto;
  boolean_T guard1 = false;
  int ilo;
  double ctoc;
  int ihi;
  boolean_T notdone;
  int exitg2;
  double stemp_im;
  int i;
  int n;
  double cto1;
  int j;
  int jcol;
  double a;
  boolean_T exitg3;
  int jrow;
  int b_i;
  int nzcount;
  creal_T atmp;
  boolean_T exitg4;
  int exitg1;
  At_size[0] = A_size[0];
  At_size[1] = A_size[1];
  ii = A_size[0] * A_size[1];
  if (0 <= ii - 1) {
    std::memcpy(&At_data[0], &A_data[0], ii * sizeof(creal_T));
  }

  anrm = 0.0;
  for (jcolp1 = 0; jcolp1 < ii; jcolp1++) {
    absxk = rt_hypotd(A_data[jcolp1].re, A_data[jcolp1].im);
    if (absxk > anrm) {
      anrm = absxk;
    }
  }

  ilascl = false;
  anrmto = anrm;
  guard1 = false;
  if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
    anrmto = 6.7178761075670888E-139;
    ilascl = true;
    guard1 = true;
  } else {
    if (anrm > 1.4885657073574029E+138) {
      anrmto = 1.4885657073574029E+138;
      ilascl = true;
      guard1 = true;
    }
  }

  if (guard1) {
    absxk = anrm;
    ctoc = anrmto;
    notdone = true;
    while (notdone) {
      stemp_im = absxk * 2.0041683600089728E-292;
      cto1 = ctoc / 4.9896007738368E+291;
      if ((stemp_im > ctoc) && (ctoc != 0.0)) {
        a = 2.0041683600089728E-292;
        absxk = stemp_im;
      } else if (cto1 > absxk) {
        a = 4.9896007738368E+291;
        ctoc = cto1;
      } else {
        a = ctoc / absxk;
        notdone = false;
      }

      ii = At_size[0] * At_size[1];
      for (b_i = 0; b_i < ii; b_i++) {
        At_data[b_i].re *= a;
        At_data[b_i].im *= a;
      }
    }
  }

  ilo = 1;
  ihi = A_size[0];
  if (A_size[0] <= 1) {
    ihi = 1;
  } else {
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      notdone = false;
      ii = ihi;
      exitg3 = false;
      while ((!exitg3) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        jcol = 0;
        exitg4 = false;
        while ((!exitg4) && (jcol <= ihi - 1)) {
          b_i = (ii + At_size[0] * jcol) - 1;
          if ((At_data[b_i].re != 0.0) || (At_data[b_i].im != 0.0) || (ii ==
               jcol + 1)) {
            if (nzcount == 0) {
              j = jcol + 1;
              nzcount = 1;
              jcol++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            jcol++;
          }
        }

        if (nzcount < 2) {
          notdone = true;
          exitg3 = true;
        } else {
          ii--;
        }
      }

      if (!notdone) {
        exitg2 = 2;
      } else {
        n = At_size[0];
        if (i != ihi) {
          for (jcolp1 = 1; jcolp1 <= n; jcolp1++) {
            ii = At_size[0] * (jcolp1 - 1);
            nzcount = (i + ii) - 1;
            atmp = At_data[nzcount];
            jcol = (ihi + ii) - 1;
            At_data[nzcount] = At_data[jcol];
            At_data[jcol] = atmp;
          }
        }

        if (j != ihi) {
          for (jcolp1 = 0; jcolp1 < ihi; jcolp1++) {
            ii = jcolp1 + At_size[0] * (j - 1);
            atmp = At_data[ii];
            jcol = jcolp1 + At_size[0] * (ihi - 1);
            At_data[ii] = At_data[jcol];
            At_data[jcol] = atmp;
          }
        }

        ihi--;
        if (ihi == 1) {
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0);

    if (exitg2 != 1) {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        notdone = false;
        jcol = ilo;
        exitg3 = false;
        while ((!exitg3) && (jcol <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = jcol;
          ii = ilo;
          exitg4 = false;
          while ((!exitg4) && (ii <= ihi)) {
            b_i = (ii + At_size[0] * (jcol - 1)) - 1;
            if ((At_data[b_i].re != 0.0) || (At_data[b_i].im != 0.0) || (ii ==
                 jcol)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                ii++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              ii++;
            }
          }

          if (nzcount < 2) {
            notdone = true;
            exitg3 = true;
          } else {
            jcol++;
          }
        }

        if (!notdone) {
          exitg1 = 1;
        } else {
          n = At_size[0];
          if (i != ilo) {
            for (jcolp1 = ilo; jcolp1 <= n; jcolp1++) {
              ii = At_size[0] * (jcolp1 - 1);
              nzcount = (i + ii) - 1;
              atmp = At_data[nzcount];
              jcol = (ilo + ii) - 1;
              At_data[nzcount] = At_data[jcol];
              At_data[jcol] = atmp;
            }
          }

          if (j != ilo) {
            for (jcolp1 = 0; jcolp1 < ihi; jcolp1++) {
              ii = jcolp1 + At_size[0] * (j - 1);
              atmp = At_data[ii];
              jcol = jcolp1 + At_size[0] * (ilo - 1);
              At_data[ii] = At_data[jcol];
              At_data[jcol] = atmp;
            }
          }

          ilo++;
          if (ilo == ihi) {
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }
  }

  n = A_size[0];
  if ((A_size[0] > 1) && (ihi >= ilo + 2)) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      jcolp1 = jcol + 2;
      for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
        xzlartg(At_data[(jrow + At_size[0] * jcol) - 1], At_data[jrow + At_size
                [0] * jcol], &absxk, &atmp, &At_data[(jrow + At_size[0] * jcol)
                - 1]);
        b_i = jrow + At_size[0] * jcol;
        At_data[b_i].re = 0.0;
        At_data[b_i].im = 0.0;
        for (j = jcolp1; j <= n; j++) {
          ii = jrow + At_size[0] * (j - 1);
          nzcount = ii - 1;
          ctoc = absxk * At_data[nzcount].re + (atmp.re * At_data[ii].re -
            atmp.im * At_data[ii].im);
          stemp_im = absxk * At_data[nzcount].im + (atmp.re * At_data[ii].im +
            atmp.im * At_data[ii].re);
          cto1 = At_data[nzcount].re;
          At_data[ii].re = absxk * At_data[ii].re - (atmp.re * At_data[nzcount].
            re + atmp.im * At_data[nzcount].im);
          At_data[ii].im = absxk * At_data[ii].im - (atmp.re * At_data[nzcount].
            im - atmp.im * cto1);
          At_data[nzcount].re = ctoc;
          At_data[nzcount].im = stemp_im;
        }

        atmp.re = -atmp.re;
        atmp.im = -atmp.im;
        for (i = 1; i <= ihi; i++) {
          ii = (i + At_size[0] * (jrow - 1)) - 1;
          nzcount = (i + At_size[0] * jrow) - 1;
          ctoc = absxk * At_data[nzcount].re + (atmp.re * At_data[ii].re -
            atmp.im * At_data[ii].im);
          stemp_im = absxk * At_data[nzcount].im + (atmp.re * At_data[ii].im +
            atmp.im * At_data[ii].re);
          cto1 = At_data[nzcount].re;
          At_data[ii].re = absxk * At_data[ii].re - (atmp.re * At_data[nzcount].
            re + atmp.im * At_data[nzcount].im);
          At_data[ii].im = absxk * At_data[ii].im - (atmp.re * At_data[nzcount].
            im - atmp.im * cto1);
          At_data[nzcount].re = ctoc;
          At_data[nzcount].im = stemp_im;
        }
      }
    }
  }

  xzhgeqz(At_data, At_size, ilo, ihi, info, alpha1_data, alpha1_size, beta1_data,
          beta1_size);
  if ((*info == 0) && ilascl) {
    notdone = true;
    while (notdone) {
      stemp_im = anrmto * 2.0041683600089728E-292;
      cto1 = anrm / 4.9896007738368E+291;
      if ((stemp_im > anrm) && (anrm != 0.0)) {
        a = 2.0041683600089728E-292;
        anrmto = stemp_im;
      } else if (cto1 > anrmto) {
        a = 4.9896007738368E+291;
        anrm = cto1;
      } else {
        a = anrm / anrmto;
        notdone = false;
      }

      ii = alpha1_size[0];
      for (b_i = 0; b_i < ii; b_i++) {
        alpha1_data[b_i].re *= a;
        alpha1_data[b_i].im *= a;
      }
    }
  }
}

//
// Arguments    : double A[400]
//                int ipiv[20]
//                int *info
// Return Type  : void
//
static void xzgetrf(double A[400], int ipiv[20], int *info)
{
  int i;
  int j;
  int mmj_tmp;
  int b;
  int jj;
  int jp1j;
  int iy;
  int jA;
  int ix;
  double smax;
  int k;
  double s;
  for (i = 0; i < 20; i++) {
    ipiv[i] = i + 1;
  }

  *info = 0;
  for (j = 0; j < 19; j++) {
    mmj_tmp = 18 - j;
    b = j * 21;
    jj = j * 21;
    jp1j = b + 2;
    iy = 20 - j;
    jA = 0;
    ix = b;
    smax = std::abs(A[jj]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (A[jj + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = iy + 1;
        ix = j;
        for (k = 0; k < 20; k++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 20;
          iy += 20;
        }
      }

      i = (jj - j) + 20;
      for (iy = jp1j; iy <= i; iy++) {
        A[iy - 1] /= A[jj];
      }
    } else {
      *info = j + 1;
    }

    iy = b + 20;
    jA = jj;
    for (k = 0; k <= mmj_tmp; k++) {
      smax = A[iy];
      if (A[iy] != 0.0) {
        ix = jj + 1;
        i = jA + 22;
        b = (jA - j) + 40;
        for (jp1j = i; jp1j <= b; jp1j++) {
          A[jp1j - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 20;
      jA += 20;
    }
  }

  if ((*info == 0) && (A[399] == 0.0)) {
    *info = 20;
  }
}

//
// Arguments    : creal_T V[100]
//                int ilo
//                int ihi
//                const int rscale[10]
// Return Type  : void
//
static void xzggbak(creal_T V[100], int ilo, int ihi, const int rscale[10])
{
  int i;
  int b_i;
  int k;
  int j;
  int tmp_re_tmp;
  double tmp_re;
  double tmp_im;
  int i1;
  if (ilo > 1) {
    for (i = ilo - 2; i + 1 >= 1; i--) {
      k = rscale[i] - 1;
      if (rscale[i] != i + 1) {
        for (j = 0; j < 10; j++) {
          tmp_re_tmp = i + 10 * j;
          tmp_re = V[tmp_re_tmp].re;
          tmp_im = V[tmp_re_tmp].im;
          b_i = k + 10 * j;
          V[tmp_re_tmp] = V[b_i];
          V[b_i].re = tmp_re;
          V[b_i].im = tmp_im;
        }
      }
    }
  }

  if (ihi < 10) {
    b_i = ihi + 1;
    for (i = b_i; i < 11; i++) {
      k = rscale[i - 1];
      if (k != i) {
        for (j = 0; j < 10; j++) {
          tmp_re_tmp = (i + 10 * j) - 1;
          tmp_re = V[tmp_re_tmp].re;
          tmp_im = V[tmp_re_tmp].im;
          i1 = (k + 10 * j) - 1;
          V[tmp_re_tmp] = V[i1];
          V[i1].re = tmp_re;
          V[i1].im = tmp_im;
        }
      }
    }
  }
}

//
// Arguments    : creal_T A[100]
//                int *ilo
//                int *ihi
//                int rscale[10]
// Return Type  : void
//
static void xzggbal(creal_T A[100], int *ilo, int *ihi, int rscale[10])
{
  int i;
  int exitg2;
  int j;
  boolean_T found;
  int ii;
  boolean_T exitg3;
  int nzcount;
  int jj;
  boolean_T exitg4;
  double atmp_re;
  double atmp_im;
  int exitg1;
  int A_tmp;
  for (i = 0; i < 10; i++) {
    rscale[i] = 1;
  }

  *ilo = 1;
  *ihi = 10;
  do {
    exitg2 = 0;
    i = 0;
    j = 0;
    found = false;
    ii = *ihi;
    exitg3 = false;
    while ((!exitg3) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = *ihi;
      jj = 0;
      exitg4 = false;
      while ((!exitg4) && (jj <= *ihi - 1)) {
        A_tmp = (ii + 10 * jj) - 1;
        if ((A[A_tmp].re != 0.0) || (A[A_tmp].im != 0.0) || (ii == jj + 1)) {
          if (nzcount == 0) {
            j = jj + 1;
            nzcount = 1;
            jj++;
          } else {
            nzcount = 2;
            exitg4 = true;
          }
        } else {
          jj++;
        }
      }

      if (nzcount < 2) {
        found = true;
        exitg3 = true;
      } else {
        ii--;
      }
    }

    if (!found) {
      exitg2 = 2;
    } else {
      if (i != *ihi) {
        for (nzcount = 0; nzcount < 10; nzcount++) {
          jj = (i + 10 * nzcount) - 1;
          atmp_re = A[jj].re;
          atmp_im = A[jj].im;
          ii = (*ihi + 10 * nzcount) - 1;
          A[jj] = A[ii];
          A[ii].re = atmp_re;
          A[ii].im = atmp_im;
        }
      }

      if (j != *ihi) {
        for (nzcount = 0; nzcount < *ihi; nzcount++) {
          jj = nzcount + 10 * (j - 1);
          atmp_re = A[jj].re;
          atmp_im = A[jj].im;
          ii = nzcount + 10 * (*ihi - 1);
          A[jj] = A[ii];
          A[ii].re = atmp_re;
          A[ii].im = atmp_im;
        }
      }

      rscale[*ihi - 1] = j;
      (*ihi)--;
      if (*ihi == 1) {
        rscale[0] = 1;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 != 1) {
    do {
      exitg1 = 0;
      i = 0;
      j = 0;
      found = false;
      jj = *ilo;
      exitg3 = false;
      while ((!exitg3) && (jj <= *ihi)) {
        nzcount = 0;
        i = *ihi;
        j = jj;
        ii = *ilo;
        exitg4 = false;
        while ((!exitg4) && (ii <= *ihi)) {
          A_tmp = (ii + 10 * (jj - 1)) - 1;
          if ((A[A_tmp].re != 0.0) || (A[A_tmp].im != 0.0) || (ii == jj)) {
            if (nzcount == 0) {
              i = ii;
              nzcount = 1;
              ii++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            ii++;
          }
        }

        if (nzcount < 2) {
          found = true;
          exitg3 = true;
        } else {
          jj++;
        }
      }

      if (!found) {
        exitg1 = 1;
      } else {
        if (i != *ilo) {
          for (nzcount = *ilo; nzcount < 11; nzcount++) {
            ii = 10 * (nzcount - 1);
            jj = (i + ii) - 1;
            atmp_re = A[jj].re;
            atmp_im = A[jj].im;
            ii = (*ilo + ii) - 1;
            A[jj] = A[ii];
            A[ii].re = atmp_re;
            A[ii].im = atmp_im;
          }
        }

        if (j != *ilo) {
          for (nzcount = 0; nzcount < *ihi; nzcount++) {
            jj = nzcount + 10 * (j - 1);
            atmp_re = A[jj].re;
            atmp_im = A[jj].im;
            ii = nzcount + 10 * (*ilo - 1);
            A[jj] = A[ii];
            A[ii].re = atmp_re;
            A[ii].im = atmp_im;
          }
        }

        rscale[*ilo - 1] = j;
        (*ilo)++;
        if (*ilo == *ihi) {
          rscale[*ilo - 1] = *ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : creal_T A[100]
//                int *info
//                creal_T alpha1[10]
//                creal_T beta1[10]
//                creal_T V[100]
// Return Type  : void
//
static void xzggev(creal_T A[100], int *info, creal_T alpha1[10], creal_T beta1
                   [10], creal_T V[100])
{
  double anrm;
  int k;
  boolean_T ilascl;
  double absxk;
  double anrmto;
  boolean_T guard1 = false;
  int ihi;
  int rscale[10];
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double a;
  int jr;
  anrm = 0.0;
  for (k = 0; k < 100; k++) {
    absxk = rt_hypotd(A[k].re, A[k].im);
    if (absxk > anrm) {
      anrm = absxk;
    }
  }

  ilascl = false;
  anrmto = anrm;
  guard1 = false;
  if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
    anrmto = 6.7178761075670888E-139;
    ilascl = true;
    guard1 = true;
  } else {
    if (anrm > 1.4885657073574029E+138) {
      anrmto = 1.4885657073574029E+138;
      ilascl = true;
      guard1 = true;
    }
  }

  if (guard1) {
    absxk = anrm;
    ctoc = anrmto;
    notdone = true;
    while (notdone) {
      cfrom1 = absxk * 2.0041683600089728E-292;
      cto1 = ctoc / 4.9896007738368E+291;
      if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
        a = 2.0041683600089728E-292;
        absxk = cfrom1;
      } else if (cto1 > absxk) {
        a = 4.9896007738368E+291;
        ctoc = cto1;
      } else {
        a = ctoc / absxk;
        notdone = false;
      }

      for (k = 0; k < 100; k++) {
        A[k].re *= a;
        A[k].im *= a;
      }
    }
  }

  xzggbal(A, &k, &ihi, rscale);
  xzgghrd(k, ihi, A, V);
  c_xzhgeqz(A, k, ihi, V, info, alpha1, beta1);
  if (*info == 0) {
    xztgevc(A, V);
    xzggbak(V, k, ihi, rscale);
    for (ihi = 0; ihi < 10; ihi++) {
      absxk = std::abs(V[10 * ihi].re) + std::abs(V[10 * ihi].im);
      for (jr = 0; jr < 9; jr++) {
        k = (jr + 10 * ihi) + 1;
        ctoc = std::abs(V[k].re) + std::abs(V[k].im);
        if (ctoc > absxk) {
          absxk = ctoc;
        }
      }

      if (absxk >= 6.7178761075670888E-139) {
        absxk = 1.0 / absxk;
        for (jr = 0; jr < 10; jr++) {
          k = jr + 10 * ihi;
          V[k].re *= absxk;
          V[k].im *= absxk;
        }
      }
    }

    if (ilascl) {
      notdone = true;
      while (notdone) {
        cfrom1 = anrmto * 2.0041683600089728E-292;
        cto1 = anrm / 4.9896007738368E+291;
        if ((cfrom1 > anrm) && (anrm != 0.0)) {
          a = 2.0041683600089728E-292;
          anrmto = cfrom1;
        } else if (cto1 > anrmto) {
          a = 4.9896007738368E+291;
          anrm = cto1;
        } else {
          a = anrm / anrmto;
          notdone = false;
        }

        for (k = 0; k < 10; k++) {
          alpha1[k].re *= a;
          alpha1[k].im *= a;
        }
      }
    }
  }
}

//
// Arguments    : int ilo
//                int ihi
//                creal_T A[100]
//                creal_T Z[100]
// Return Type  : void
//
static void xzgghrd(int ilo, int ihi, creal_T A[100], creal_T Z[100])
{
  signed char b_I[100];
  int k;
  int jcol;
  int jcolp1;
  int jrow;
  double c;
  creal_T s;
  int s_re_tmp;
  int stemp_re_tmp;
  double stemp_re;
  double stemp_im;
  double d;
  double d1;
  double d2;
  std::memset(&b_I[0], 0, 100U * sizeof(signed char));
  for (k = 0; k < 10; k++) {
    b_I[k + 10 * k] = 1;
  }

  for (k = 0; k < 100; k++) {
    Z[k].re = b_I[k];
    Z[k].im = 0.0;
  }

  if (ihi >= ilo + 2) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      jcolp1 = jcol + 2;
      for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
        k = jrow + 10 * jcol;
        xzlartg(A[k - 1], A[k], &c, &s, &A[(jrow + 10 * jcol) - 1]);
        A[k].re = 0.0;
        A[k].im = 0.0;
        for (k = jcolp1; k < 11; k++) {
          s_re_tmp = jrow + 10 * (k - 1);
          stemp_re_tmp = s_re_tmp - 1;
          stemp_re = c * A[stemp_re_tmp].re + (s.re * A[s_re_tmp].re - s.im *
            A[s_re_tmp].im);
          stemp_im = c * A[stemp_re_tmp].im + (s.re * A[s_re_tmp].im + s.im *
            A[s_re_tmp].re);
          d = A[stemp_re_tmp].im;
          d1 = A[stemp_re_tmp].re;
          A[s_re_tmp].re = c * A[s_re_tmp].re - (s.re * A[stemp_re_tmp].re +
            s.im * A[stemp_re_tmp].im);
          A[s_re_tmp].im = c * A[s_re_tmp].im - (s.re * d - s.im * d1);
          A[stemp_re_tmp].re = stemp_re;
          A[stemp_re_tmp].im = stemp_im;
        }

        s.re = -s.re;
        s.im = -s.im;
        for (k = 1; k <= ihi; k++) {
          s_re_tmp = (k + 10 * (jrow - 1)) - 1;
          stemp_re_tmp = (k + 10 * jrow) - 1;
          stemp_re = c * A[stemp_re_tmp].re + (s.re * A[s_re_tmp].re - s.im *
            A[s_re_tmp].im);
          stemp_im = c * A[stemp_re_tmp].im + (s.re * A[s_re_tmp].im + s.im *
            A[s_re_tmp].re);
          d = A[stemp_re_tmp].im;
          d1 = A[stemp_re_tmp].re;
          A[s_re_tmp].re = c * A[s_re_tmp].re - (s.re * A[stemp_re_tmp].re +
            s.im * A[stemp_re_tmp].im);
          A[s_re_tmp].im = c * A[s_re_tmp].im - (s.re * d - s.im * d1);
          A[stemp_re_tmp].re = stemp_re;
          A[stemp_re_tmp].im = stemp_im;
        }

        d = s.re;
        d1 = s.im;
        for (k = 0; k < 10; k++) {
          s_re_tmp = k + 10 * (jrow - 1);
          stemp_re_tmp = k + 10 * jrow;
          stemp_re = c * Z[stemp_re_tmp].re + (d * Z[s_re_tmp].re - d1 *
            Z[s_re_tmp].im);
          stemp_im = c * Z[stemp_re_tmp].im + (d * Z[s_re_tmp].im + d1 *
            Z[s_re_tmp].re);
          d2 = Z[stemp_re_tmp].re;
          Z[s_re_tmp].re = c * Z[s_re_tmp].re - (d * Z[stemp_re_tmp].re + d1 *
            Z[stemp_re_tmp].im);
          Z[s_re_tmp].im = c * Z[s_re_tmp].im - (d * Z[stemp_re_tmp].im - d1 *
            d2);
          Z[stemp_re_tmp].re = stemp_re;
          Z[stemp_re_tmp].im = stemp_im;
        }
      }
    }
  }
}

//
// Arguments    : const creal_T A_data[]
//                const int A_size[2]
//                int ilo
//                int ihi
//                int *info
//                creal_T alpha1_data[]
//                int alpha1_size[1]
//                creal_T beta1_data[]
//                int beta1_size[1]
// Return Type  : void
//
static void xzhgeqz(const creal_T A_data[], const int A_size[2], int ilo, int
                    ihi, int *info, creal_T alpha1_data[], int alpha1_size[1],
                    creal_T beta1_data[], int beta1_size[1])
{
  int A_size_idx_0;
  int jm1;
  creal_T b_A_data[64];
  int n;
  int ctemp_tmp;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double b_atol;
  boolean_T firstNonZero;
  int j;
  double ascale;
  int i;
  double bscale;
  double temp1;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  double temp2;
  int ifirst;
  int istart;
  int ilast;
  int ilastm1;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T b_guard1 = false;
  boolean_T guard3 = false;
  boolean_T exitg2;
  double ascale_re;
  creal_T shift;
  double ascale_im;
  double ad22_re;
  double ad22_im;
  double t1_re;
  double t1_im;
  double t1_im_tmp;
  creal_T b_ascale;
  A_size_idx_0 = A_size[0];
  jm1 = A_size[0] * A_size[1];
  if (0 <= jm1 - 1) {
    std::memcpy(&b_A_data[0], &A_data[0], jm1 * sizeof(creal_T));
  }

  *info = 0;
  if ((A_size[0] == 1) && (A_size[1] == 1)) {
    ihi = 1;
  }

  n = A_size[0];
  alpha1_size[0] = A_size[0];
  jm1 = A_size[0];
  if (0 <= jm1 - 1) {
    std::memset(&alpha1_data[0], 0, jm1 * sizeof(creal_T));
  }

  beta1_size[0] = A_size[0];
  jm1 = A_size[0];
  for (ctemp_tmp = 0; ctemp_tmp < jm1; ctemp_tmp++) {
    beta1_data[ctemp_tmp].re = 1.0;
    beta1_data[ctemp_tmp].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (ilo <= ihi) {
    scale = 0.0;
    anorm = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ctemp_tmp = j + 1;
      if (ihi < j + 1) {
        ctemp_tmp = ihi;
      }

      for (i = ilo; i <= ctemp_tmp; i++) {
        jm1 = (i + A_size[0] * (j - 1)) - 1;
        if (A_data[jm1].re != 0.0) {
          temp1 = std::abs(A_data[jm1].re);
          if (firstNonZero) {
            anorm = 1.0;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }

        if (A_data[jm1].im != 0.0) {
          temp1 = std::abs(A_data[jm1].im);
          if (firstNonZero) {
            anorm = 1.0;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * std::sqrt(anorm);
  }

  scale = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (scale > 2.2250738585072014E-308) {
    b_atol = scale;
  }

  scale = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    scale = anorm;
  }

  ascale = 1.0 / scale;
  bscale = 1.0 / std::sqrt(static_cast<double>(A_size[0]));
  firstNonZero = true;
  ctemp_tmp = ihi + 1;
  for (j = ctemp_tmp; j <= n; j++) {
    alpha1_data[j - 1] = A_data[(j + A_size[0] * (j - 1)) - 1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          ctemp_tmp = ilast + A_size_idx_0 * ilastm1;
          if (std::abs(b_A_data[ctemp_tmp].re) + std::abs(b_A_data[ctemp_tmp].im)
              <= b_atol) {
            b_A_data[ctemp_tmp].re = 0.0;
            b_A_data[ctemp_tmp].im = 0.0;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                ctemp_tmp = j + A_size_idx_0 * (j - 1);
                if (std::abs(b_A_data[ctemp_tmp].re) + std::abs
                    (b_A_data[ctemp_tmp].im) <= b_atol) {
                  b_A_data[ctemp_tmp].re = 0.0;
                  b_A_data[ctemp_tmp].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }

            if (goto70) {
              b_guard1 = true;
            } else {
              jm1 = alpha1_size[0];
              if (0 <= jm1 - 1) {
                std::memset(&alpha1_data[0], 0, jm1 * sizeof(creal_T));
              }

              jm1 = beta1_size[0];
              if (0 <= jm1 - 1) {
                std::memset(&beta1_data[0], 0, jm1 * sizeof(creal_T));
              }

              *info = 1;
              exitg1 = 1;
            }
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1_data[ilast] = b_A_data[ilast + A_size_idx_0 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              ilastm = ilast + 1;
              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                jm1 = ilastm1 + A_size_idx_0 * ilastm1;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0) {
                  shift.re = anorm / bscale;
                  shift.im = 0.0;
                } else if (anorm == 0.0) {
                  shift.re = 0.0;
                  shift.im = scale / bscale;
                } else {
                  shift.re = anorm / bscale;
                  shift.im = scale / bscale;
                }

                jm1 = ilast + A_size_idx_0 * ilast;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0) {
                  ad22_re = anorm / bscale;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = scale / bscale;
                } else {
                  ad22_re = anorm / bscale;
                  ad22_im = scale / bscale;
                }

                t1_re = 0.5 * (shift.re + ad22_re);
                t1_im = 0.5 * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                jm1 = ilastm1 + A_size_idx_0 * ilast;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0) {
                  ascale_re = anorm / bscale;
                  ascale_im = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  ascale_im = scale / bscale;
                } else {
                  ascale_re = anorm / bscale;
                  ascale_im = scale / bscale;
                }

                jm1 = ilast + A_size_idx_0 * ilastm1;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0) {
                  temp2 = anorm / bscale;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  temp2 = 0.0;
                  anorm = scale / bscale;
                } else {
                  temp2 = anorm / bscale;
                  anorm = scale / bscale;
                }

                scale = shift.re * ad22_re - shift.im * ad22_im;
                temp1 = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (ascale_re * temp2
                  - ascale_im * anorm)) - scale;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (ascale_re * anorm +
                  ascale_im * temp2)) - temp1;
                b_sqrt(&shift);
                if ((t1_re - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                    0.0) {
                  shift.re += t1_re;
                  shift.im += t1_im;
                } else {
                  shift.re = t1_re - shift.re;
                  shift.im = t1_im - shift.im;
                }
              } else {
                jm1 = ilast + A_size_idx_0 * ilastm1;
                anorm = ascale * b_A_data[jm1].re;
                scale = ascale * b_A_data[jm1].im;
                if (scale == 0.0) {
                  ascale_re = anorm / bscale;
                  ascale_im = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  ascale_im = scale / bscale;
                } else {
                  ascale_re = anorm / bscale;
                  ascale_im = scale / bscale;
                }

                eshift_re += ascale_re;
                eshift_im += ascale_im;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              n = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                jm1 = A_size_idx_0 * j;
                ctemp_tmp = j + jm1;
                ctemp.re = ascale * b_A_data[ctemp_tmp].re - shift.re * bscale;
                ctemp.im = ascale * b_A_data[ctemp_tmp].im - shift.im * bscale;
                anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
                jm1 += n;
                temp2 = ascale * (std::abs(b_A_data[jm1].re) + std::abs
                                  (b_A_data[jm1].im));
                scale = anorm;
                if (temp2 > anorm) {
                  scale = temp2;
                }

                if ((scale < 1.0) && (scale != 0.0)) {
                  anorm /= scale;
                  temp2 /= scale;
                }

                ctemp_tmp = j + A_size_idx_0 * (j - 1);
                if ((std::abs(b_A_data[ctemp_tmp].re) + std::abs
                     (b_A_data[ctemp_tmp].im)) * temp2 <= anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  n = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp_tmp = (ifirst + A_size_idx_0 * (ifirst - 1)) - 1;
                ctemp.re = ascale * b_A_data[ctemp_tmp].re - shift.re * bscale;
                ctemp.im = ascale * b_A_data[ctemp_tmp].im - shift.im * bscale;
              }

              goto90 = false;
              jm1 = istart + A_size_idx_0 * (istart - 1);
              b_ascale.re = ascale * b_A_data[jm1].re;
              b_ascale.im = ascale * b_A_data[jm1].im;
              b_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              jm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  xzlartg(b_A_data[(j + A_size_idx_0 * jm1) - 1], b_A_data[j +
                          A_size_idx_0 * jm1], &anorm, &shift, &b_A_data[(j +
                           A_size_idx_0 * jm1) - 1]);
                  ctemp_tmp = j + A_size_idx_0 * jm1;
                  b_A_data[ctemp_tmp].re = 0.0;
                  b_A_data[ctemp_tmp].im = 0.0;
                }

                for (n = j; n <= ilastm; n++) {
                  jm1 = j + A_size_idx_0 * (n - 1);
                  ctemp_tmp = jm1 - 1;
                  ad22_re = anorm * b_A_data[ctemp_tmp].re + (shift.re *
                    b_A_data[jm1].re - shift.im * b_A_data[jm1].im);
                  ad22_im = anorm * b_A_data[ctemp_tmp].im + (shift.re *
                    b_A_data[jm1].im + shift.im * b_A_data[jm1].re);
                  scale = b_A_data[ctemp_tmp].re;
                  b_A_data[jm1].re = anorm * b_A_data[jm1].re - (shift.re *
                    b_A_data[ctemp_tmp].re + shift.im * b_A_data[ctemp_tmp].im);
                  b_A_data[jm1].im = anorm * b_A_data[jm1].im - (shift.re *
                    b_A_data[ctemp_tmp].im - shift.im * scale);
                  b_A_data[ctemp_tmp].re = ad22_re;
                  b_A_data[ctemp_tmp].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                n = j;
                if (ilast + 1 < j + 2) {
                  n = ilast - 1;
                }

                for (i = ifirst; i <= n + 2; i++) {
                  jm1 = (i + A_size_idx_0 * (j - 1)) - 1;
                  ctemp_tmp = (i + A_size_idx_0 * j) - 1;
                  ad22_re = anorm * b_A_data[ctemp_tmp].re + (shift.re *
                    b_A_data[jm1].re - shift.im * b_A_data[jm1].im);
                  ad22_im = anorm * b_A_data[ctemp_tmp].im + (shift.re *
                    b_A_data[jm1].im + shift.im * b_A_data[jm1].re);
                  scale = b_A_data[ctemp_tmp].re;
                  b_A_data[jm1].re = anorm * b_A_data[jm1].re - (shift.re *
                    b_A_data[ctemp_tmp].re + shift.im * b_A_data[ctemp_tmp].im);
                  b_A_data[jm1].im = anorm * b_A_data[jm1].im - (shift.re *
                    b_A_data[ctemp_tmp].im - shift.im * scale);
                  b_A_data[ctemp_tmp].re = ad22_re;
                  b_A_data[ctemp_tmp].im = ad22_im;
                }

                jm1 = j - 1;
                j++;
              }
            }

            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (firstNonZero) {
      *info = ilast + 1;
      if (0 <= ilast) {
        std::memset(&alpha1_data[0], 0, (ilast + 1) * sizeof(creal_T));
        std::memset(&beta1_data[0], 0, (ilast + 1) * sizeof(creal_T));
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1_data[j] = b_A_data[j + A_size_idx_0 * j];
    }
  }
}

//
// Arguments    : int m
//                int n
//                int iv0
//                double tau
//                double C[100]
//                int ic0
//                double work[10]
// Return Type  : void
//
static void xzlarf(int m, int n, int iv0, double tau, double C[100], int ic0,
                   double work[10])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int jy;
  int b_i;
  int j;
  int ia;
  int ix;
  int exitg1;
  double c;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      i = ic0 + lastc * 10;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    if (lastc + 1 != 0) {
      if (0 <= lastc) {
        std::memset(&work[0], 0, (lastc + 1) * sizeof(double));
      }

      i = 0;
      b_i = ic0 + 10 * lastc;
      for (jy = ic0; jy <= b_i; jy += 10) {
        ix = iv0;
        c = 0.0;
        j = (jy + lastv) - 1;
        for (ia = jy; ia <= j; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[i] += c;
        i++;
      }
    }

    if (-tau != 0.0) {
      i = ic0;
      jy = 0;
      for (j = 0; j <= lastc; j++) {
        if (work[jy] != 0.0) {
          c = work[jy] * -tau;
          ix = iv0;
          b_i = lastv + i;
          for (ia = i; ia < b_i; ia++) {
            C[ia - 1] += C[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += 10;
      }
    }
  }
}

//
// Arguments    : const creal_T f
//                const creal_T g
//                double *cs
//                creal_T *sn
//                creal_T *r
// Return Type  : void
//
static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r)
{
  double scale_tmp;
  double f2s;
  double scale;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  double f2;
  double g2;
  scale_tmp = std::abs(f.re);
  f2s = std::abs(f.im);
  if (f2s > scale_tmp) {
    scale_tmp = f2s;
  }

  f2s = std::abs(g.re);
  scale = std::abs(g.im);
  if (scale > f2s) {
    f2s = scale;
  }

  scale = scale_tmp;
  if (f2s > scale_tmp) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = -1;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      count++;
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    scale = g2;
    if (1.0 > g2) {
      scale = 1.0;
    }

    if (f2 <= scale * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd(g.re, g.im);
        r->im = 0.0;
        f2 = rt_hypotd(gs_re, gs_im);
        sn->re = gs_re / f2;
        sn->im = -gs_im / f2;
      } else {
        g2 = std::sqrt(g2);
        *cs = rt_hypotd(fs_re, fs_im) / g2;
        if (scale_tmp > 1.0) {
          f2 = rt_hypotd(f.re, f.im);
          fs_re = f.re / f2;
          fs_im = f.im / f2;
        } else {
          scale = 7.4428285367870146E+137 * f.re;
          f2s = 7.4428285367870146E+137 * f.im;
          f2 = rt_hypotd(scale, f2s);
          fs_re = scale / f2;
          fs_im = f2s / f2;
        }

        gs_re /= g2;
        gs_im = -gs_im / g2;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = std::sqrt(g2 / f2 + 1.0);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      f2 += g2;
      f2s = r->re / f2;
      scale = r->im / f2;
      sn->re = f2s * gs_re - scale * -gs_im;
      sn->im = f2s * -gs_im + scale * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 0; rescaledir <= count; rescaledir++) {
            r->re *= 1.3435752215134178E-138;
            r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

//
// Arguments    : const creal_T A[100]
//                creal_T V[100]
// Return Type  : void
//
static void xztgevc(const creal_T A[100], creal_T V[100])
{
  double rworka[10];
  double anorm;
  int j;
  double xmx;
  int i;
  int re_tmp;
  double ascale;
  double d_re;
  int je;
  int x_tmp_tmp;
  double temp;
  double salpha_re;
  double salpha_im;
  double acoeff;
  boolean_T lscalea;
  double z;
  boolean_T lscaleb;
  double scale;
  creal_T work1[10];
  double dmin;
  int b_i;
  int jr;
  creal_T work2[10];
  double d_im;
  double brm;
  std::memset(&rworka[0], 0, 10U * sizeof(double));
  anorm = std::abs(A[0].re) + std::abs(A[0].im);
  for (j = 0; j < 9; j++) {
    for (i = 0; i <= j; i++) {
      re_tmp = i + 10 * (j + 1);
      rworka[j + 1] += std::abs(A[re_tmp].re) + std::abs(A[re_tmp].im);
    }

    i = (j + 10 * (j + 1)) + 1;
    d_re = rworka[j + 1] + (std::abs(A[i].re) + std::abs(A[i].im));
    if (d_re > anorm) {
      anorm = d_re;
    }
  }

  xmx = anorm;
  if (2.2250738585072014E-308 > anorm) {
    xmx = 2.2250738585072014E-308;
  }

  ascale = 1.0 / xmx;
  for (je = 0; je < 10; je++) {
    x_tmp_tmp = 10 * (9 - je);
    i = (x_tmp_tmp - je) + 9;
    xmx = (std::abs(A[i].re) + std::abs(A[i].im)) * ascale;
    if (1.0 > xmx) {
      xmx = 1.0;
    }

    temp = 1.0 / xmx;
    salpha_re = ascale * (temp * A[i].re);
    salpha_im = ascale * (temp * A[i].im);
    acoeff = temp * ascale;
    if ((temp >= 2.2250738585072014E-308) && (acoeff < 1.0020841800044864E-291))
    {
      lscalea = true;
    } else {
      lscalea = false;
    }

    z = std::abs(salpha_re) + std::abs(salpha_im);
    if ((z >= 2.2250738585072014E-308) && (z < 1.0020841800044864E-291)) {
      lscaleb = true;
    } else {
      lscaleb = false;
    }

    scale = 1.0;
    if (lscalea) {
      xmx = anorm;
      if (9.9792015476736E+290 < anorm) {
        xmx = 9.9792015476736E+290;
      }

      scale = 1.0020841800044864E-291 / temp * xmx;
    }

    if (lscaleb) {
      d_re = 1.0020841800044864E-291 / z;
      if (d_re > scale) {
        scale = d_re;
      }
    }

    if (lscalea || lscaleb) {
      xmx = acoeff;
      if (1.0 > acoeff) {
        xmx = 1.0;
      }

      if (z > xmx) {
        xmx = z;
      }

      d_re = 1.0 / (2.2250738585072014E-308 * xmx);
      if (d_re < scale) {
        scale = d_re;
      }

      if (lscalea) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      salpha_re *= scale;
      salpha_im *= scale;
    }

    std::memset(&work1[0], 0, 10U * sizeof(creal_T));
    work1[9 - je].re = 1.0;
    work1[9 - je].im = 0.0;
    dmin = 2.2204460492503131E-16 * acoeff * anorm;
    d_re = 2.2204460492503131E-16 * (std::abs(salpha_re) + std::abs(salpha_im));
    if (d_re > dmin) {
      dmin = d_re;
    }

    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }

    b_i = 8 - je;
    for (jr = 0; jr <= b_i; jr++) {
      i = jr + x_tmp_tmp;
      work1[jr].re = acoeff * A[i].re;
      work1[jr].im = acoeff * A[i].im;
    }

    work1[9 - je].re = 1.0;
    work1[9 - je].im = 0.0;
    b_i = static_cast<int>((((-1.0 - ((-static_cast<double>(je) + 10.0) - 1.0))
      + 1.0) / -1.0));
    for (j = 0; j < b_i; j++) {
      re_tmp = 8 - (je + j);
      i = re_tmp + 10 * re_tmp;
      d_re = acoeff * A[i].re - salpha_re;
      d_im = acoeff * A[i].im - salpha_im;
      if (std::abs(d_re) + std::abs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0;
      }

      brm = std::abs(d_re);
      scale = std::abs(d_im);
      xmx = brm + scale;
      if (xmx < 1.0) {
        z = std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im);
        if (z >= 4.49423283715579E+306 * xmx) {
          temp = 1.0 / z;
          i = 9 - je;
          for (jr = 0; jr <= i; jr++) {
            work1[jr].re *= temp;
            work1[jr].im *= temp;
          }
        }
      }

      if (d_im == 0.0) {
        if (-work1[re_tmp].im == 0.0) {
          scale = -work1[re_tmp].re / d_re;
          xmx = 0.0;
        } else if (-work1[re_tmp].re == 0.0) {
          scale = 0.0;
          xmx = -work1[re_tmp].im / d_re;
        } else {
          scale = -work1[re_tmp].re / d_re;
          xmx = -work1[re_tmp].im / d_re;
        }
      } else if (d_re == 0.0) {
        if (-work1[re_tmp].re == 0.0) {
          scale = -work1[re_tmp].im / d_im;
          xmx = 0.0;
        } else if (-work1[re_tmp].im == 0.0) {
          scale = 0.0;
          xmx = -(-work1[re_tmp].re / d_im);
        } else {
          scale = -work1[re_tmp].im / d_im;
          xmx = -(-work1[re_tmp].re / d_im);
        }
      } else if (brm > scale) {
        z = d_im / d_re;
        xmx = d_re + z * d_im;
        scale = (-work1[re_tmp].re + z * -work1[re_tmp].im) / xmx;
        xmx = (-work1[re_tmp].im - z * -work1[re_tmp].re) / xmx;
      } else if (scale == brm) {
        if (d_re > 0.0) {
          z = 0.5;
        } else {
          z = -0.5;
        }

        if (d_im > 0.0) {
          xmx = 0.5;
        } else {
          xmx = -0.5;
        }

        scale = (-work1[re_tmp].re * z + -work1[re_tmp].im * xmx) / brm;
        xmx = (-work1[re_tmp].im * z - -work1[re_tmp].re * xmx) / brm;
      } else {
        z = d_re / d_im;
        xmx = d_im + z * d_re;
        scale = (z * -work1[re_tmp].re + -work1[re_tmp].im) / xmx;
        xmx = (z * -work1[re_tmp].im - (-work1[re_tmp].re)) / xmx;
      }

      work1[re_tmp].re = scale;
      work1[re_tmp].im = xmx;
      if (re_tmp + 1 > 1) {
        if (std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im) > 1.0) {
          temp = 1.0 / (std::abs(work1[re_tmp].re) + std::abs(work1[re_tmp].im));
          if (acoeff * rworka[re_tmp] >= 4.49423283715579E+306 * temp) {
            i = 9 - je;
            for (jr = 0; jr <= i; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        d_re = acoeff * work1[re_tmp].re;
        d_im = acoeff * work1[re_tmp].im;
        for (jr = 0; jr < re_tmp; jr++) {
          i = jr + 10 * re_tmp;
          work1[jr].re += d_re * A[i].re - d_im * A[i].im;
          work1[jr].im += d_re * A[i].im + d_im * A[i].re;
        }
      }
    }

    std::memset(&work2[0], 0, 10U * sizeof(creal_T));
    b_i = 9 - je;
    for (i = 0; i <= b_i; i++) {
      xmx = work1[i].re;
      z = work1[i].im;
      for (jr = 0; jr < 10; jr++) {
        re_tmp = jr + 10 * i;
        work2[jr].re += V[re_tmp].re * xmx - V[re_tmp].im * z;
        work2[jr].im += V[re_tmp].re * z + V[re_tmp].im * xmx;
      }
    }

    xmx = std::abs(work2[0].re) + std::abs(work2[0].im);
    for (jr = 0; jr < 9; jr++) {
      d_re = std::abs(work2[jr + 1].re) + std::abs(work2[jr + 1].im);
      if (d_re > xmx) {
        xmx = d_re;
      }
    }

    if (xmx > 2.2250738585072014E-308) {
      temp = 1.0 / xmx;
      for (jr = 0; jr < 10; jr++) {
        b_i = jr + x_tmp_tmp;
        V[b_i].re = temp * work2[jr].re;
        V[b_i].im = temp * work2[jr].im;
      }
    } else {
      std::memset(&V[x_tmp_tmp], 0, 10U * sizeof(creal_T));
    }
  }
}

//
// function [n,f,r,t] = p35p_double(X, x, y, e)
// Arguments    : const double X[12]
//                const double x[4]
//                const double y[4]
//                double e
//                int *n
//                double f_data[]
//                int f_size[2]
//                double r_data[]
//                int r_size[3]
//                double t_data[]
//                int t_size[2]
// Return Type  : void
//
void p35p_double(const double X[12], const double x[4], const double y[4],
                 double e, int *n, double f_data[], int f_size[2], double
                 r_data[], int r_size[3], double t_data[], int t_size[2])
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
  double A[400];
  double C[200];
  double M[100];
  int b_i;
  double b_M[100];
  creal_T W[100];
  creal_T D[100];
  int i1;
  double qx;
  int ar_tmp;
  int br_tmp;
  double f;
  double W_re;
  double bim;
  double s;
  double fc_set_data[2];
  int fc_set_size[2];
  double fs_set_data[2];
  int fs_set_size[2];
  int j;
  double s_tmp;
  double R_xy[9];
  double R_xy_tmp;
  double b_x;
  double T[3];
  double fc_set[9];
  double b_fc_set[12];
  double c_fc_set[9];
  double p4[3];
  double b_X[4];
  double R_curr[9];
  double tmp_data[99];
  int tmp_size[3];
  double b_t_data[30];
  int t_data_tmp;

  // 'p35p_double:2' assert(isa(X, 'double'));
  // 'p35p_double:3' assert(isa(x, 'double'));
  // 'p35p_double:4' assert(isa(y, 'double'));
  // 'p35p_double:5' assert(isa(e, 'double'));
  // 'p35p_double:7' [n, f, r, t] = p35p_solver(X, x, y, e);
  // 'p35p_solver:2' type = class(X);
  // 'p35p_solver:3' assert(strcmp(type, class(y)));
  // 'p35p_solver:4' assert(strcmp(type, class(x)));
  // 'p35p_solver:5' assert(strcmp(type, class(e)));
  // 'p35p_solver:6' solution_num = int32(0);
  *n = 0;

  // 'p35p_solver:7' f_sol = zeros(1, 0, 'like', X);
  f_size[0] = 1;
  f_size[1] = 0;

  // 'p35p_solver:8' coder.varsize('f_sol', [1 10], [0 1]);
  // 'p35p_solver:9' R_sol = zeros(3, 3, 0, 'like', X);
  r_size[0] = 3;
  r_size[1] = 3;
  r_size[2] = 0;

  // 'p35p_solver:10' coder.varsize('R_sol', [3 3 10], [0 0 1]);
  // 'p35p_solver:11' T_sol = zeros(3, 0, 'like', X);
  t_size[0] = 3;
  t_size[1] = 0;

  // 'p35p_solver:12' coder.varsize('T_sol', [3 10], [0 1]);
  // 'p35p_solver:14' R = init_R(type);
  // 'p35p_solver:16' F = init_F(x, y, X, R);
  // 'init_F:2' F = zeros(4, 3, 6, 'like', x);
  std::memset(&F[0], 0, 72U * sizeof(double));

  // 'init_F:3' F(1, :, :) = quadruple_constraint(1, 2, 3, x, y, X, R);
  quadruple_constraint(1.0, 2.0, 3.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i] = dv[3 * i];
    F[12 * i + 4] = dv[3 * i + 1];
    F[12 * i + 8] = dv[3 * i + 2];
  }

  // 'init_F:4' F(2, :, :) = quadruple_constraint(1, 3, 2, x, y, X, R);
  quadruple_constraint(1.0, 3.0, 2.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 1] = dv[3 * i];
    F[12 * i + 5] = dv[3 * i + 1];
    F[12 * i + 9] = dv[3 * i + 2];
  }

  // 'init_F:5' F(3, :, :) = quadruple_constraint(2, 4, 3, x, y, X, R);
  quadruple_constraint(2.0, 4.0, 3.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 2] = dv[3 * i];
    F[12 * i + 6] = dv[3 * i + 1];
    F[12 * i + 10] = dv[3 * i + 2];
  }

  // 'init_F:6' F(4, :, :) = quadruple_constraint(3, 4, 2, x, y, X, R);
  quadruple_constraint(3.0, 4.0, 2.0, x, y, X, R, dv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 3] = dv[3 * i];
    F[12 * i + 7] = dv[3 * i + 1];
    F[12 * i + 11] = dv[3 * i + 2];
  }

  // 4x3x6
  // 'p35p_solver:17' G4 = equations_for_groebner(F);
  // 4x28
  // 'p35p_solver:19' G20 = mult_for_groebner(G4);
  equations_for_groebner(F, dv1);
  mult_for_groebner(dv1, G20);

  // 20x45
  // 'p35p_solver:20' A = [G20(:, 1:2), G20(:, 10:13), G20(:, 18:22), G20(:, 25:28), G20(:, 31:35)]; 
  for (i = 0; i < 2; i++) {
    std::memcpy(&A[i * 20], &G20[i * 20], 20U * sizeof(double));
  }

  for (i = 0; i < 4; i++) {
    std::memcpy(&A[i * 20 + 40], &G20[i * 20 + 180], 20U * sizeof(double));
  }

  for (i = 0; i < 5; i++) {
    std::memcpy(&A[i * 20 + 120], &G20[i * 20 + 340], 20U * sizeof(double));
  }

  for (i = 0; i < 4; i++) {
    std::memcpy(&A[i * 20 + 220], &G20[i * 20 + 480], 20U * sizeof(double));
  }

  for (i = 0; i < 5; i++) {
    std::memcpy(&A[i * 20 + 300], &G20[i * 20 + 600], 20U * sizeof(double));
  }

  // 'p35p_solver:21' B = G20(:, 36:45);
  // 'p35p_solver:22' if rcond(A) < e
  if (rcond(A) >= e) { //todo: change back
  //if(true){
    // 'p35p_solver:26' C = A\B;
    for (i = 0; i < 10; i++) {
      std::memcpy(&C[i * 20], &G20[i * 20 + 700], 20U * sizeof(double));
    }

    b_mldivide(A, C);

    // 'p35p_solver:27' M = make_mult_matrix(C);
    // this function creates a matrix for multiplication by x in the monomial
    // basis B
    // monomial basis B = {x^3, ...., 1} -- monomials up to the 3d degree, #B = 10 
    // x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1
    // 'make_mult_matrix:6' M = zeros(10, 10, 'like', C);
    std::memset(&M[0], 0, 100U * sizeof(double));

    // 'make_mult_matrix:7' for i = 1 : 4
    for (b_i = 0; b_i < 4; b_i++) {
      // 'make_mult_matrix:8' M(:, i) = -C(15 + i, :)';
      for (i = 0; i < 10; i++) {
        M[i + 10 * b_i] = -C[(b_i + 20 * i) + 15];
      }
    }

    // 'make_mult_matrix:10' M(1, 5) = 1;
    M[40] = 1.0;

    // 'make_mult_matrix:11' M(2, 6) = 1;
    M[51] = 1.0;

    // 'make_mult_matrix:12' M(3, 7) = 1;
    M[62] = 1.0;

    // 'make_mult_matrix:13' M(5, 8) = 1;
    M[74] = 1.0;

    // 'make_mult_matrix:14' M(6, 9) = 1;
    M[85] = 1.0;

    // 'make_mult_matrix:15' M(8, 10) = 1;
    M[97] = 1.0;

    // 'p35p_solver:29' [W,D] = eig(M');
    for (i = 0; i < 10; i++) {
      for (i1 = 0; i1 < 10; i1++) {
        b_M[i1 + 10 * i] = M[i + 10 * i1];
      }
    }

    eig(b_M, W, D);

    // 'p35p_solver:31' for i = 1 : 10
    for (b_i = 0; b_i < 10; b_i++) {
      // 'p35p_solver:32' qx = D(i, i);
      // 'p35p_solver:33' qy = W(9, i) / W(10, i);
      // 'p35p_solver:35' qx = real(qx);
      qx = D[b_i + 10 * b_i].re;

      // 'p35p_solver:36' qy = real(qy);
      ar_tmp = 10 * b_i + 8;
      br_tmp = 10 * b_i + 9;
      if (W[br_tmp].im == 0.0) {
        if (W[ar_tmp].im == 0.0) {
          W_re = W[ar_tmp].re / W[br_tmp].re;
        } else if (W[ar_tmp].re == 0.0) {
          W_re = 0.0;
        } else {
          W_re = W[ar_tmp].re / W[br_tmp].re;
        }
      } else if (W[br_tmp].re == 0.0) {
        if (W[ar_tmp].re == 0.0) {
          W_re = W[ar_tmp].im / W[br_tmp].im;
        } else if (W[ar_tmp].im == 0.0) {
          W_re = 0.0;
        } else {
          W_re = W[ar_tmp].im / W[br_tmp].im;
        }
      } else {
        f = std::abs(W[br_tmp].re);
        bim = std::abs(W[br_tmp].im);
        if (f > bim) {
          s = W[br_tmp].im / W[br_tmp].re;
          W_re = (W[ar_tmp].re + s * W[ar_tmp].im) / (W[br_tmp].re + s *
            W[br_tmp].im);
        } else if (bim == f) {
          if (W[br_tmp].re > 0.0) {
            W_re = 0.5;
          } else {
            W_re = -0.5;
          }

          if (W[br_tmp].im > 0.0) {
            bim = 0.5;
          } else {
            bim = -0.5;
          }

          W_re = (W[ar_tmp].re * W_re + W[ar_tmp].im * bim) / f;
        } else {
          s = W[br_tmp].re / W[br_tmp].im;
          W_re = (s * W[ar_tmp].re + W[ar_tmp].im) / (W[br_tmp].im + s *
            W[br_tmp].re);
        }
      }

      // 'p35p_solver:37' [fc_set, fs_set, f_num] = find_f(F, qx, qy, e);
      find_f(F, qx, W_re, e, fc_set_data, fc_set_size, fs_set_data, fs_set_size,
             &f);

      // 'p35p_solver:39' for j = 1 : f_num
      i = static_cast<int>(f);
      if (0 <= i - 1) {
        bim = W_re * W_re;
        s_tmp = qx * qx;
        s = 1.0 / ((s_tmp + 1.0) + bim);
        R_xy[0] = 1.0 - 2.0 * s * bim;
        f = 2.0 * s * qx;
        R_xy_tmp = f * W_re;
        R_xy[3] = R_xy_tmp;
        R_xy[6] = 2.0 * s * W_re;
        R_xy[1] = R_xy_tmp;
        R_xy[4] = 1.0 - 2.0 * s * s_tmp;
        R_xy[7] = -2.0 * s * qx;
        R_xy[2] = -2.0 * s * W_re;
        R_xy[5] = f;
        R_xy[8] = 1.0 - 2.0 * s * (s_tmp + bim);
        b_x = x[1];
      }

      for (j = 0; j < i; j++) {
        // 'p35p_solver:40' fc = fc_set(j);
        // 'p35p_solver:41' fs = fs_set(j);
        // 'p35p_solver:43' s = 1/(1 + qx^2 + qy^2);
        // 'p35p_solver:44' R_xy = [1 - 2*s*qy^2, 2*s*qx*qy,     2*s*qy;
        // 'p35p_solver:45'                       2*s*qx*qy,    1 - 2*s*qx^2, -2*s*qx; 
        // 'p35p_solver:46'                      -2*s*qy,       2*s*qx,        1 - 2*s*(qx^2 + qy^2)]; 
        // li = find_lamda(R, fc, fs, Xi, Xj, xi, xj) -- signature
        // 'p35p_solver:49' lambda1 = find_lamda(R_xy, fc, fs, X(:, 1), X(:, 2), x(1), x(2)); 
        // 'find_lamda:2' li = -(fc*R(1, :) - fs*R(2, :) - xj*R(3, :))*(Xj - Xi)/(xi - xj); 
        R_xy_tmp = fc_set_data[j] * R_xy[0] - fs_set_data[j] * R_xy[1];
        f = fc_set_data[j] * R_xy[3] - fs_set_data[j] * R_xy[4];
        bim = fc_set_data[j] * R_xy[6] - fs_set_data[j] * R_xy[7];
        s_tmp = ((-(R_xy_tmp - b_x * R_xy[2]) * (X[3] - X[0]) + -(f - b_x *
                   R_xy[5]) * (X[4] - X[1])) + -(bim - b_x * R_xy[8]) * (X[5] -
                  X[2])) / (x[0] - x[1]);

        // T = find_translation(R, fc, fs, li, Xi, xi, yi)
        // 'p35p_solver:51' T = find_translation(R_xy, fc, fs, lambda1, X(:, 1), x(1), y(1)); 
        // 'find_translation:2' T = [li*xi - (fc*R(1, :) - fs*R(2, :))*Xi;
        // 'find_translation:3'          li*yi - (fs*R(1, :) + fc*R(2, :))*Xi;
        // 'find_translation:4'          li - R(3, :)*Xi];
        T[0] = s_tmp * x[0] - ((R_xy_tmp * X[0] + f * X[1]) + bim * X[2]);
        T[1] = s_tmp * y[0] - (((fs_set_data[j] * R_xy[0] + fc_set_data[j] *
          R_xy[1]) * X[0] + (fs_set_data[j] * R_xy[3] + fc_set_data[j] * R_xy[4])
          * X[1]) + (fs_set_data[j] * R_xy[6] + fc_set_data[j] * R_xy[7]) * X[2]);
        T[2] = s_tmp - ((R_xy[2] * X[0] + R_xy[5] * X[1]) + R_xy[8] * X[2]);

        // 'p35p_solver:53' f = hypot(fs, fc);
        f = rt_hypotd(fs_set_data[j], fc_set_data[j]);

        // 'p35p_solver:54' K = [fc, -fs, 0;
        // 'p35p_solver:55'                  fs,  fc, 0;
        // 'p35p_solver:56'                    0,  0, 1];
        // 'p35p_solver:57' P = [K*R_xy, T];
        // 'p35p_solver:59' p4 = P*[X(:, 4); 1];
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
          bim = fc_set[i1 + 3];
          ar_tmp = static_cast<int>(fc_set[i1 + 6]);
          for (br_tmp = 0; br_tmp < 3; br_tmp++) {
            c_fc_set[i1 + 3 * br_tmp] = (fc_set[i1] * R_xy[3 * br_tmp] + bim *
              R_xy[3 * br_tmp + 1]) + static_cast<double>(ar_tmp) * R_xy[3 *
              br_tmp + 2];
          }
        }

        for (i1 = 0; i1 < 3; i1++) {
          b_fc_set[3 * i1] = c_fc_set[3 * i1];
          ar_tmp = 3 * i1 + 1;
          b_fc_set[ar_tmp] = c_fc_set[ar_tmp];
          ar_tmp = 3 * i1 + 2;
          b_fc_set[ar_tmp] = c_fc_set[ar_tmp];
          b_fc_set[i1 + 9] = T[i1];
          b_X[i1] = X[i1 + 9];
        }

        for (i1 = 0; i1 < 3; i1++) {
          p4[i1] = ((b_fc_set[i1] * b_X[0] + b_fc_set[i1 + 3] * b_X[1]) +
                    b_fc_set[i1 + 6] * b_X[2]) + b_fc_set[i1 + 9];
        }

        // 'p35p_solver:60' y4 = p4(2)/p4(3);
        // 'p35p_solver:62' if abs(y4 - y(4)) < 0.01*f
        if (std::abs(p4[1] / p4[2] - y[3]) < 0.01 * f) {
          // 'p35p_solver:63' solution_num = solution_num + 1;
          (*n)++;

          // 'p35p_solver:64' R_z = [fc/f, -fs/f, 0;
          // 'p35p_solver:65'                        fs/f,  fc/f, 0;
          // 'p35p_solver:66'                        0,        0, 1];
          // 'p35p_solver:67' R_curr = R_z*R_xy;
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

          // 'p35p_solver:68' f_sol = [f_sol, f];
          i1 = f_size[1];
          f_size[1]++;
          f_data[i1] = f;

          // 'p35p_solver:69' T = -R_curr'*([1/f; 1/f; 1].*T);
          for (i1 = 0; i1 < 3; i1++) {
            ar_tmp = static_cast<int>(fc_set[i1 + 6]);
            bim = fc_set[i1 + 3];
            for (br_tmp = 0; br_tmp < 3; br_tmp++) {
              s_tmp = (fc_set[i1] * R_xy[3 * br_tmp] + bim * R_xy[3 * br_tmp + 1])
                + static_cast<double>(ar_tmp) * R_xy[3 * br_tmp + 2];
              R_curr[i1 + 3 * br_tmp] = s_tmp;
              c_fc_set[br_tmp + 3 * i1] = -s_tmp;
            }
          }

          bim = 1.0 / f * T[0];
          s_tmp = 1.0 / f * T[1];
          f = T[2];
          for (i1 = 0; i1 < 3; i1++) {
            T[i1] = (c_fc_set[i1] * bim + c_fc_set[i1 + 3] * s_tmp) +
              c_fc_set[i1 + 6] * f;
          }

          // 'p35p_solver:70' if solution_num == 1
          if (*n == 1) {
            // 'p35p_solver:71' R_sol = R_curr;
            r_size[0] = 3;
            r_size[1] = 3;
            r_size[2] = 1;
            std::memcpy(&r_data[0], &R_curr[0], 9U * sizeof(double));
          } else {
            // 'p35p_solver:72' else
            // 'p35p_solver:73' R_sol = cat(3, R_sol, R_curr);
            cat(r_data, r_size, R_curr, tmp_data, tmp_size);
            r_size[0] = 3;
            r_size[1] = 3;
            r_size[2] = tmp_size[2];
            ar_tmp = tmp_size[0] * tmp_size[1] * tmp_size[2];
            if (0 <= ar_tmp - 1) {
              std::memcpy(&r_data[0], &tmp_data[0], ar_tmp * sizeof(double));
            }
          }

          // 'p35p_solver:75' T_sol = [T_sol, T];
          ar_tmp = t_size[1];
          br_tmp = t_size[1] + 1;
          for (i1 = 0; i1 < ar_tmp; i1++) {
            b_t_data[3 * i1] = t_data[3 * i1];
            t_data_tmp = 3 * i1 + 1;
            b_t_data[t_data_tmp] = t_data[t_data_tmp];
            t_data_tmp = 3 * i1 + 2;
            b_t_data[t_data_tmp] = t_data[t_data_tmp];
          }

          b_t_data[3 * t_size[1]] = T[0];
          b_t_data[3 * t_size[1] + 1] = T[1];
          b_t_data[3 * t_size[1] + 2] = T[2];
          t_size[0] = 3;
          t_size[1] = br_tmp;
          ar_tmp = 3 * br_tmp;
          if (0 <= ar_tmp - 1) {
            std::memcpy(&t_data[0], &b_t_data[0], ar_tmp * sizeof(double));
          }
        }
      }
    }
  }

  // 'p35p_double:8' assert(isa(n, 'int32'));
  // 'p35p_double:9' assert(isa(f, 'double'));
  // 'p35p_double:10' assert(isa(r, 'double'));
  // 'p35p_double:11' assert(isa(t, 'double'));
}

//
// function [n,f,r,t] = p35p_single(X, x, y, e)
// Arguments    : const float X[12]
//                const float x[4]
//                const float y[4]
//                float e
//                int *n
//                float f_data[]
//                int f_size[2]
//                float r_data[]
//                int r_size[3]
//                float t_data[]
//                int t_size[2]
// Return Type  : void
//
void p35p_single(const float X[12], const float x[4], const float y[4], float e,
                 int *n, float f_data[], int f_size[2], float r_data[], int
                 r_size[3], float t_data[], int t_size[2])
{
  float F[72];
  static const float R[54] = { 1.0F, 0.0F, 0.0F, 0.0F, -1.0F, 0.0F, 0.0F, 0.0F,
    -1.0F, 0.0F, 2.0F, 0.0F, 2.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, -1.0F, 0.0F,
    0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, -1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
    2.0F, 0.0F, -2.0F, 0.0F, 0.0F, 0.0F, -2.0F, 0.0F, 0.0F, 0.0F, 2.0F, 0.0F,
    0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 1.0F };

  float fv[18];
  int i;
  float fv1[112];
  float G20[900];
  float A[400];
  float C[200];
  float M[100];
  int b_i;
  float b_M[100];
  creal32_T W[100];
  creal32_T D[100];
  int i1;
  float qx;
  int ar_tmp;
  int br_tmp;
  float brm;
  float W_re;
  float f;
  float s;
  float fc_set_data[2];
  int fc_set_size[2];
  float fs_set_data[2];
  int fs_set_size[2];
  double f_num;
  int j;
  float R_xy_tmp;
  float b_R_xy_tmp;
  float R_xy[9];
  float b_x;
  float T[3];
  float fc_set[9];
  float b_fc_set[12];
  float R_curr[9];
  float p4[3];
  float b_X[4];
  float tmp_data[99];
  int tmp_size[3];
  float b_t_data[30];
  int t_data_tmp;

  // 'p35p_single:2' assert(isa(X, 'single'));
  // 'p35p_single:3' assert(isa(x, 'single'));
  // 'p35p_single:4' assert(isa(y, 'single'));
  // 'p35p_single:5' assert(isa(e, 'single'));
  // 'p35p_single:7' [n, f, r, t] = p35p_solver(X, x, y, e);
  // 'p35p_solver:2' type = class(X);
  // 'p35p_solver:3' assert(strcmp(type, class(y)));
  // 'p35p_solver:4' assert(strcmp(type, class(x)));
  // 'p35p_solver:5' assert(strcmp(type, class(e)));
  // 'p35p_solver:6' solution_num = int32(0);
  *n = 0;

  // 'p35p_solver:7' f_sol = zeros(1, 0, 'like', X);
  f_size[0] = 1;
  f_size[1] = 0;

  // 'p35p_solver:8' coder.varsize('f_sol', [1 10], [0 1]);
  // 'p35p_solver:9' R_sol = zeros(3, 3, 0, 'like', X);
  r_size[0] = 3;
  r_size[1] = 3;
  r_size[2] = 0;

  // 'p35p_solver:10' coder.varsize('R_sol', [3 3 10], [0 0 1]);
  // 'p35p_solver:11' T_sol = zeros(3, 0, 'like', X);
  t_size[0] = 3;
  t_size[1] = 0;

  // 'p35p_solver:12' coder.varsize('T_sol', [3 10], [0 1]);
  // 'p35p_solver:14' R = init_R(type);
  // 'p35p_solver:16' F = init_F(x, y, X, R);
  // 'init_F:2' F = zeros(4, 3, 6, 'like', x);
  std::memset(&F[0], 0, 72U * sizeof(float));

  // 'init_F:3' F(1, :, :) = quadruple_constraint(1, 2, 3, x, y, X, R);
  b_quadruple_constraint(1.0, 2.0, 3.0, x, y, X, R, fv);
  for (i = 0; i < 6; i++) {
    F[12 * i] = fv[3 * i];
    F[12 * i + 4] = fv[3 * i + 1];
    F[12 * i + 8] = fv[3 * i + 2];
  }

  // 'init_F:4' F(2, :, :) = quadruple_constraint(1, 3, 2, x, y, X, R);
  b_quadruple_constraint(1.0, 3.0, 2.0, x, y, X, R, fv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 1] = fv[3 * i];
    F[12 * i + 5] = fv[3 * i + 1];
    F[12 * i + 9] = fv[3 * i + 2];
  }

  // 'init_F:5' F(3, :, :) = quadruple_constraint(2, 4, 3, x, y, X, R);
  b_quadruple_constraint(2.0, 4.0, 3.0, x, y, X, R, fv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 2] = fv[3 * i];
    F[12 * i + 6] = fv[3 * i + 1];
    F[12 * i + 10] = fv[3 * i + 2];
  }

  // 'init_F:6' F(4, :, :) = quadruple_constraint(3, 4, 2, x, y, X, R);
  b_quadruple_constraint(3.0, 4.0, 2.0, x, y, X, R, fv);
  for (i = 0; i < 6; i++) {
    F[12 * i + 3] = fv[3 * i];
    F[12 * i + 7] = fv[3 * i + 1];
    F[12 * i + 11] = fv[3 * i + 2];
  }

  // 4x3x6
  // 'p35p_solver:17' G4 = equations_for_groebner(F);
  // 4x28
  // 'p35p_solver:19' G20 = mult_for_groebner(G4);
  b_equations_for_groebner(F, fv1);
  b_mult_for_groebner(fv1, G20);

  // 20x45
  // 'p35p_solver:20' A = [G20(:, 1:2), G20(:, 10:13), G20(:, 18:22), G20(:, 25:28), G20(:, 31:35)]; 
  for (i = 0; i < 2; i++) {
    std::memcpy(&A[i * 20], &G20[i * 20], 20U * sizeof(float));
  }

  for (i = 0; i < 4; i++) {
    std::memcpy(&A[i * 20 + 40], &G20[i * 20 + 180], 20U * sizeof(float));
  }

  for (i = 0; i < 5; i++) {
    std::memcpy(&A[i * 20 + 120], &G20[i * 20 + 340], 20U * sizeof(float));
  }

  for (i = 0; i < 4; i++) {
    std::memcpy(&A[i * 20 + 220], &G20[i * 20 + 480], 20U * sizeof(float));
  }

  for (i = 0; i < 5; i++) {
    std::memcpy(&A[i * 20 + 300], &G20[i * 20 + 600], 20U * sizeof(float));
  }

  // 'p35p_solver:21' B = G20(:, 36:45);
  // 'p35p_solver:22' if rcond(A) < e
  if (b_rcond(A) >= e) { //todo: change back
  //if (true) {
    // 'p35p_solver:26' C = A\B;
    for (i = 0; i < 10; i++) {
      std::memcpy(&C[i * 20], &G20[i * 20 + 700], 20U * sizeof(float));
    }

    c_mldivide(A, C);

    // 'p35p_solver:27' M = make_mult_matrix(C);
    // this function creates a matrix for multiplication by x in the monomial
    // basis B
    // monomial basis B = {x^3, ...., 1} -- monomials up to the 3d degree, #B = 10 
    // x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1
    // 'make_mult_matrix:6' M = zeros(10, 10, 'like', C);
    std::memset(&M[0], 0, 100U * sizeof(float));

    // 'make_mult_matrix:7' for i = 1 : 4
    for (b_i = 0; b_i < 4; b_i++) {
      // 'make_mult_matrix:8' M(:, i) = -C(15 + i, :)';
      for (i = 0; i < 10; i++) {
        M[i + 10 * b_i] = -C[(b_i + 20 * i) + 15];
      }
    }

    // 'make_mult_matrix:10' M(1, 5) = 1;
    M[40] = 1.0F;

    // 'make_mult_matrix:11' M(2, 6) = 1;
    M[51] = 1.0F;

    // 'make_mult_matrix:12' M(3, 7) = 1;
    M[62] = 1.0F;

    // 'make_mult_matrix:13' M(5, 8) = 1;
    M[74] = 1.0F;

    // 'make_mult_matrix:14' M(6, 9) = 1;
    M[85] = 1.0F;

    // 'make_mult_matrix:15' M(8, 10) = 1;
    M[97] = 1.0F;

    // 'p35p_solver:29' [W,D] = eig(M');
    for (i = 0; i < 10; i++) {
      for (i1 = 0; i1 < 10; i1++) {
        b_M[i1 + 10 * i] = M[i + 10 * i1];
      }
    }

    b_eig(b_M, W, D);

    // 'p35p_solver:31' for i = 1 : 10
    for (b_i = 0; b_i < 10; b_i++) {
      // 'p35p_solver:32' qx = D(i, i);
      // 'p35p_solver:33' qy = W(9, i) / W(10, i);
      // 'p35p_solver:35' qx = real(qx);
      qx = D[b_i + 10 * b_i].re;

      // 'p35p_solver:36' qy = real(qy);
      ar_tmp = 10 * b_i + 8;
      br_tmp = 10 * b_i + 9;
      if (W[br_tmp].im == 0.0F) {
        if (W[ar_tmp].im == 0.0F) {
          W_re = W[ar_tmp].re / W[br_tmp].re;
        } else if (W[ar_tmp].re == 0.0F) {
          W_re = 0.0F;
        } else {
          W_re = W[ar_tmp].re / W[br_tmp].re;
        }
      } else if (W[br_tmp].re == 0.0F) {
        if (W[ar_tmp].re == 0.0F) {
          W_re = W[ar_tmp].im / W[br_tmp].im;
        } else if (W[ar_tmp].im == 0.0F) {
          W_re = 0.0F;
        } else {
          W_re = W[ar_tmp].im / W[br_tmp].im;
        }
      } else {
        brm = std::abs(W[br_tmp].re);
        f = std::abs(W[br_tmp].im);
        if (brm > f) {
          s = W[br_tmp].im / W[br_tmp].re;
          W_re = (W[ar_tmp].re + s * W[ar_tmp].im) / (W[br_tmp].re + s *
            W[br_tmp].im);
        } else if (f == brm) {
          if (W[br_tmp].re > 0.0F) {
            W_re = 0.5F;
          } else {
            W_re = -0.5F;
          }

          if (W[br_tmp].im > 0.0F) {
            f = 0.5F;
          } else {
            f = -0.5F;
          }

          W_re = (W[ar_tmp].re * W_re + W[ar_tmp].im * f) / brm;
        } else {
          s = W[br_tmp].re / W[br_tmp].im;
          W_re = (s * W[ar_tmp].re + W[ar_tmp].im) / (W[br_tmp].im + s *
            W[br_tmp].re);
        }
      }

      // 'p35p_solver:37' [fc_set, fs_set, f_num] = find_f(F, qx, qy, e);
      b_find_f(F, qx, W_re, e, fc_set_data, fc_set_size, fs_set_data,
               fs_set_size, &f_num);

      // 'p35p_solver:39' for j = 1 : f_num
      i = static_cast<int>(f_num);
      if (0 <= i - 1) {
        brm = W_re * W_re;
        f = qx * qx;
        s = 1.0F / ((f + 1.0F) + brm);
        R_xy_tmp = 2.0F * s;
        b_R_xy_tmp = -2.0F * s;
        R_xy[0] = 1.0F - R_xy_tmp * brm;
        R_xy[3] = R_xy_tmp * qx * W_re;
        R_xy[6] = R_xy_tmp * W_re;
        R_xy[1] = 2.0F * s * qx * W_re;
        R_xy[4] = 1.0F - R_xy_tmp * f;
        R_xy[7] = b_R_xy_tmp * qx;
        R_xy[2] = b_R_xy_tmp * W_re;
        R_xy[5] = 2.0F * s * qx;
        R_xy[8] = 1.0F - R_xy_tmp * (f + brm);
        b_x = x[1];
      }

      for (j = 0; j < i; j++) {
        // 'p35p_solver:40' fc = fc_set(j);
        // 'p35p_solver:41' fs = fs_set(j);
        // 'p35p_solver:43' s = 1/(1 + qx^2 + qy^2);
        // 'p35p_solver:44' R_xy = [1 - 2*s*qy^2, 2*s*qx*qy,     2*s*qy;
        // 'p35p_solver:45'                       2*s*qx*qy,    1 - 2*s*qx^2, -2*s*qx; 
        // 'p35p_solver:46'                      -2*s*qy,       2*s*qx,        1 - 2*s*(qx^2 + qy^2)]; 
        // li = find_lamda(R, fc, fs, Xi, Xj, xi, xj) -- signature
        // 'p35p_solver:49' lambda1 = find_lamda(R_xy, fc, fs, X(:, 1), X(:, 2), x(1), x(2)); 
        // 'find_lamda:2' li = -(fc*R(1, :) - fs*R(2, :) - xj*R(3, :))*(Xj - Xi)/(xi - xj); 
        brm = fc_set_data[j] * R_xy[0] - fs_set_data[j] * R_xy[1];
        f = fc_set_data[j] * R_xy[3] - fs_set_data[j] * R_xy[4];
        R_xy_tmp = fc_set_data[j] * R_xy[6] - fs_set_data[j] * R_xy[7];
        b_R_xy_tmp = ((-(brm - b_x * R_xy[2]) * (X[3] - X[0]) + -(f - b_x *
          R_xy[5]) * (X[4] - X[1])) + -(R_xy_tmp - b_x * R_xy[8]) * (X[5] - X[2]))
          / (x[0] - x[1]);

        // T = find_translation(R, fc, fs, li, Xi, xi, yi)
        // 'p35p_solver:51' T = find_translation(R_xy, fc, fs, lambda1, X(:, 1), x(1), y(1)); 
        // 'find_translation:2' T = [li*xi - (fc*R(1, :) - fs*R(2, :))*Xi;
        // 'find_translation:3'          li*yi - (fs*R(1, :) + fc*R(2, :))*Xi;
        // 'find_translation:4'          li - R(3, :)*Xi];
        T[0] = b_R_xy_tmp * x[0] - ((brm * X[0] + f * X[1]) + R_xy_tmp * X[2]);
        T[1] = b_R_xy_tmp * y[0] - (((fs_set_data[j] * R_xy[0] + fc_set_data[j] *
          R_xy[1]) * X[0] + (fs_set_data[j] * R_xy[3] + fc_set_data[j] * R_xy[4])
          * X[1]) + (fs_set_data[j] * R_xy[6] + fc_set_data[j] * R_xy[7]) * X[2]);
        T[2] = b_R_xy_tmp - ((R_xy[2] * X[0] + R_xy[5] * X[1]) + R_xy[8] * X[2]);

        // 'p35p_solver:53' f = hypot(fs, fc);
        f = rt_hypotf(fs_set_data[j], fc_set_data[j]);

        // 'p35p_solver:54' K = [fc, -fs, 0;
        // 'p35p_solver:55'                  fs,  fc, 0;
        // 'p35p_solver:56'                    0,  0, 1];
        // 'p35p_solver:57' P = [K*R_xy, T];
        // 'p35p_solver:59' p4 = P*[X(:, 4); 1];
        fc_set[0] = fc_set_data[j];
        fc_set[3] = -fs_set_data[j];
        fc_set[6] = 0.0F;
        fc_set[1] = fs_set_data[j];
        fc_set[4] = fc_set_data[j];
        fc_set[7] = 0.0F;
        fc_set[2] = 0.0F;
        fc_set[5] = 0.0F;
        fc_set[8] = 1.0F;
        for (i1 = 0; i1 < 3; i1++) {
          R_xy_tmp = fc_set[i1 + 3];
          ar_tmp = static_cast<int>(fc_set[i1 + 6]);
          for (br_tmp = 0; br_tmp < 3; br_tmp++) {
            R_curr[i1 + 3 * br_tmp] = (fc_set[i1] * R_xy[3 * br_tmp] + R_xy_tmp *
              R_xy[3 * br_tmp + 1]) + static_cast<float>(ar_tmp) * R_xy[3 *
              br_tmp + 2];
          }
        }

        for (i1 = 0; i1 < 3; i1++) {
          b_fc_set[3 * i1] = R_curr[3 * i1];
          ar_tmp = 3 * i1 + 1;
          b_fc_set[ar_tmp] = R_curr[ar_tmp];
          ar_tmp = 3 * i1 + 2;
          b_fc_set[ar_tmp] = R_curr[ar_tmp];
          b_fc_set[i1 + 9] = T[i1];
          b_X[i1] = X[i1 + 9];
        }

        for (i1 = 0; i1 < 3; i1++) {
          p4[i1] = ((b_fc_set[i1] * b_X[0] + b_fc_set[i1 + 3] * b_X[1]) +
                    b_fc_set[i1 + 6] * b_X[2]) + b_fc_set[i1 + 9];
        }

        // 'p35p_solver:60' y4 = p4(2)/p4(3);
        // 'p35p_solver:62' if abs(y4 - y(4)) < 0.01*f
        if (std::abs(p4[1] / p4[2] - y[3]) < 0.01F * f) {
          // 'p35p_solver:63' solution_num = solution_num + 1;
          (*n)++;

          // 'p35p_solver:64' R_z = [fc/f, -fs/f, 0;
          // 'p35p_solver:65'                        fs/f,  fc/f, 0;
          // 'p35p_solver:66'                        0,        0, 1];
          brm = fc_set_data[j] / f;

          // 'p35p_solver:67' R_curr = R_z*R_xy;
          fc_set[0] = brm;
          fc_set[3] = -fs_set_data[j] / f;
          fc_set[6] = 0.0F;
          fc_set[1] = fs_set_data[j] / f;
          fc_set[4] = brm;
          fc_set[7] = 0.0F;
          fc_set[2] = 0.0F;
          fc_set[5] = 0.0F;
          fc_set[8] = 1.0F;
          for (i1 = 0; i1 < 3; i1++) {
            R_xy_tmp = fc_set[i1 + 3];
            ar_tmp = static_cast<int>(fc_set[i1 + 6]);
            for (br_tmp = 0; br_tmp < 3; br_tmp++) {
              R_curr[i1 + 3 * br_tmp] = (fc_set[i1] * R_xy[3 * br_tmp] +
                R_xy_tmp * R_xy[3 * br_tmp + 1]) + static_cast<float>(ar_tmp) *
                R_xy[3 * br_tmp + 2];
            }
          }

          // 'p35p_solver:68' f_sol = [f_sol, f];
          i1 = f_size[1];
          f_size[1]++;
          f_data[i1] = f;

          // 'p35p_solver:69' T = -R_curr'*([1/f; 1/f; 1].*T);
          brm = 1.0F / f;
          for (i1 = 0; i1 < 3; i1++) {
            fc_set[3 * i1] = -R_curr[i1];
            fc_set[3 * i1 + 1] = -R_curr[i1 + 3];
            fc_set[3 * i1 + 2] = -R_curr[i1 + 6];
          }

          R_xy_tmp = brm * T[0];
          brm *= T[1];
          f = T[2];
          for (i1 = 0; i1 < 3; i1++) {
            T[i1] = (fc_set[i1] * R_xy_tmp + fc_set[i1 + 3] * brm) + fc_set[i1 +
              6] * f;
          }

          // 'p35p_solver:70' if solution_num == 1
          if (*n == 1) {
            // 'p35p_solver:71' R_sol = R_curr;
            r_size[0] = 3;
            r_size[1] = 3;
            r_size[2] = 1;
            for (i1 = 0; i1 < 9; i1++) {
              r_data[i1] = R_curr[i1];
            }
          } else {
            // 'p35p_solver:72' else
            // 'p35p_solver:73' R_sol = cat(3, R_sol, R_curr);
            b_cat(r_data, r_size, R_curr, tmp_data, tmp_size);
            r_size[0] = 3;
            r_size[1] = 3;
            r_size[2] = tmp_size[2];
            ar_tmp = tmp_size[0] * tmp_size[1] * tmp_size[2];
            if (0 <= ar_tmp - 1) {
              std::memcpy(&r_data[0], &tmp_data[0], ar_tmp * sizeof(float));
            }
          }

          // 'p35p_solver:75' T_sol = [T_sol, T];
          ar_tmp = t_size[1];
          br_tmp = t_size[1] + 1;
          for (i1 = 0; i1 < ar_tmp; i1++) {
            b_t_data[3 * i1] = t_data[3 * i1];
            t_data_tmp = 3 * i1 + 1;
            b_t_data[t_data_tmp] = t_data[t_data_tmp];
            t_data_tmp = 3 * i1 + 2;
            b_t_data[t_data_tmp] = t_data[t_data_tmp];
          }

          b_t_data[3 * t_size[1]] = T[0];
          b_t_data[3 * t_size[1] + 1] = T[1];
          b_t_data[3 * t_size[1] + 2] = T[2];
          t_size[0] = 3;
          t_size[1] = br_tmp;
          ar_tmp = 3 * br_tmp;
          if (0 <= ar_tmp - 1) {
            std::memcpy(&t_data[0], &b_t_data[0], ar_tmp * sizeof(float));
          }
        }
      }
    }
  }

  // 'p35p_single:8' assert(isa(n, 'int32'));
  // 'p35p_single:9' assert(isa(t, 'single'));
  // 'p35p_single:10' assert(isa(r, 'single'));
  // 'p35p_single:11' assert(isa(f, 'single'));
}

//
// function [n,f,r,t] = p4pf_double(X, x, y, e)
// Arguments    : const double X[12]
//                const double x[4]
//                const double y[4]
//                double e
//                int *n
//                double f_data[]
//                int f_size[2]
//                double r_data[]
//                int r_size[3]
//                double t_data[]
//                int t_size[2]
// Return Type  : void
//
void p4pf_double(const double X[12], const double x[4], const double y[4],
                 double e, int *n, double f_data[], int f_size[2], double
                 r_data[], int r_size[3], double t_data[], int t_size[2])
{
  int i;
  double A[32];
  double Q[64];
  double b_A[32];
  int X_tmp;
  double b_x;
  double b_X[16];
  double d;
  double d1;
  double B[16];
  double C[16];
  double scale;
  int offset;
  int b_i;
  double ND[48];
  int B_tmp;
  int j;
  double eqs[40];
  int b_B_tmp;
  double b_ND[10];
  double c_ND[10];
  double d_ND[10];
  double e_ND[10];
  double f_ND[10];
  double g_ND[10];
  double b_eqs[30];
  double xs_data[13];
  int xs_size[2];
  double ys_data[13];
  int ys_size[2];
  double zs_data[13];
  int zs_size[2];
  double P1[3];
  double P3[3];
  double b_Q;
  double w;
  double absxk;
  double t;
  double alpha;
  double R[9];
  double c_x[9];
  int ipiv[3];
  boolean_T isodd;
  double c_A[36];
  double b[12];
  double c_Q[3];
  double d_x[12];
  double b_b[8];
  double y_data[99];

  // 'p4pf_double:2' assert(isa(X, 'double'));
  // 'p4pf_double:3' assert(isa(x, 'double'));
  // 'p4pf_double:4' assert(isa(y, 'double'));
  // 'p4pf_double:5' assert(isa(e, 'double'));
  // 'p4pf_double:7' [n, f, r, t] = solve_P4Pf(X, x, y, e);
  // SOLVE_P4Pf Summary of this function goes here
  //        X = [p1, p2, p3, p4], pi = [4, 1]; X(:, i) <-> (u(i), v(i))
  //        if f is a correct foal length, then [R, T] = [R, T] / sign(d)*abs(d)^(1/3); 
  //        where d = det(R)
  // 'solve_P4Pf:6' X = [X; ones(1, 4, 'like', X)];
  // 'solve_P4Pf:7' A = find_A(X, u, v);
  // 'solve_P4Pf:8' [Q, ~] = qr(A');
  // [p11, p12, p13, p14, p21, p22, p23, p24]'
  // 'solve_P4Pf:70' A = zeros(4, 8, 'like', X);
  // 'solve_P4Pf:71' for i = 1 : 4
  for (i = 0; i < 4; i++) {
    X_tmp = i << 2;
    b_x = X[3 * i];
    b_X[X_tmp] = b_x;
    d = X[3 * i + 1];
    b_X[X_tmp + 1] = d;
    d1 = X[3 * i + 2];
    b_X[X_tmp + 2] = d1;
    b_X[X_tmp + 3] = 1.0;

    // 'solve_P4Pf:72' A(i,  :) = [-v(i)*X(:, i)', u(i)*X(:, i)'];
    b_A[i] = -y[i] * b_x;
    b_A[i + 16] = x[i] * b_x;
    b_A[i + 4] = -y[i] * d;
    b_A[i + 20] = x[i] * d;
    b_A[i + 8] = -y[i] * d1;
    b_A[i + 24] = x[i] * d1;
    b_A[i + 12] = -y[i];
    b_A[i + 28] = x[i];
    for (b_i = 0; b_i < 8; b_i++) {
      A[b_i + (i << 3)] = b_A[i + (b_i << 2)];
    }
  }

  c_qr(A, Q, b_A);

  // 'solve_P4Pf:8' ~
  // 'solve_P4Pf:9' N = Q(:, 5:end);
  // nullspace
  // 'solve_P4Pf:10' D = find_D(X, u, v, N, e);
  // 'solve_P4Pf:77' B = zeros(4, 'double');
  // 'solve_P4Pf:78' C = zeros(4, 'double');
  // 'solve_P4Pf:79' for i = 1 : 4
  for (i = 0; i < 4; i++) {
    // 'solve_P4Pf:80' if abs(u(i)) < e
    if (std::abs(x[i]) < e) {
      // 'solve_P4Pf:81' fctr = v(i);
      scale = y[i];

      // 'solve_P4Pf:82' offset = 4;
      offset = 4;
    } else {
      // 'solve_P4Pf:83' else
      // 'solve_P4Pf:84' fctr = u(i);
      scale = x[i];

      // 'solve_P4Pf:85' offset = 0;
      offset = 0;
    }

    // 'solve_P4Pf:87' B(i, :) = fctr*X(:, i)';
    // 'solve_P4Pf:88' for j = 1 : 4
    B_tmp = i << 2;
    for (j = 0; j < 4; j++) {
      b_B_tmp = i + (j << 2);
      B[b_B_tmp] = scale * b_X[j + B_tmp];

      // 'solve_P4Pf:89' C(i, j) = X(:, i)'*ns((1:4)+offset, j);
      X_tmp = offset + ((j + 4) << 3);
      C[b_B_tmp] = ((b_X[B_tmp] * Q[X_tmp] + b_X[B_tmp + 1] * Q[X_tmp + 1]) +
                    b_X[B_tmp + 2] * Q[X_tmp + 2]) + b_X[B_tmp + 3] * Q[X_tmp +
        3];
    }
  }

  // 'solve_P4Pf:92' D_ = B\C;
  d_mldivide(B, C);

  // 'solve_P4Pf:93' if isa(X, 'single')
  // 'solve_P4Pf:95' else
  // 'solve_P4Pf:96' D = D_;
  // 'solve_P4Pf:11' eqs = find_eqs([N; D]);
  // 'solve_P4Pf:12' [n, xs, ys, zs] = solve_3Q3(eqs(1:3, :), e);
  for (b_i = 0; b_i < 4; b_i++) {
    std::memcpy(&ND[b_i * 12], &Q[b_i * 8 + 32], 8U * sizeof(double));
    X_tmp = b_i << 2;
    ND[12 * b_i + 8] = C[X_tmp];
    ND[12 * b_i + 9] = C[X_tmp + 1];
    ND[12 * b_i + 10] = C[X_tmp + 2];
    ND[12 * b_i + 11] = C[X_tmp + 3];
  }

  // ND = [N; D]
  // 'find_eqs:2' eqs = zeros(4, 10, 'like', ND);
  std::memset(&eqs[0], 0, 40U * sizeof(double));

  // 'find_eqs:3' eqs(1, :) = mult_pp(1, 1, 2, 1, ND) + mult_pp(1, 2, 2, 2, ND) + mult_pp(1, 3, 2, 3, ND); 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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
  for (b_i = 0; b_i < 10; b_i++) {
    eqs[b_i << 2] = (b_ND[b_i] + c_ND[b_i]) + d_ND[b_i];
  }

  // 'find_eqs:4' eqs(2, :) = mult_pp(3, 1, 1, 1, ND) + mult_pp(3, 2, 1, 2, ND) + mult_pp(3, 3, 1, 3, ND); 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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
  for (b_i = 0; b_i < 10; b_i++) {
    eqs[(b_i << 2) + 1] = (b_ND[b_i] + c_ND[b_i]) + d_ND[b_i];
  }

  // 'find_eqs:5' eqs(3, :) = mult_pp(3, 1, 2, 1, ND) + mult_pp(3, 2, 2, 2, ND) + mult_pp(3, 3, 2, 3, ND); 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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
  for (b_i = 0; b_i < 10; b_i++) {
    eqs[(b_i << 2) + 2] = (b_ND[b_i] + c_ND[b_i]) + d_ND[b_i];
  }

  // 'find_eqs:6' eqs(4, :) = sq_pp(1, 1, ND) + sq_pp(1, 2, ND) + sq_pp(1, 3, ND) - ... 
  // 'find_eqs:7'         sq_pp(2, 1, ND) - sq_pp(2, 2, ND) - sq_pp(2, 3, ND);
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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
  for (b_i = 0; b_i < 10; b_i++) {
    X_tmp = b_i << 2;
    eqs[X_tmp + 3] = ((((b_ND[b_i] + c_ND[b_i]) + d_ND[b_i]) - e_ND[b_i]) -
                      f_ND[b_i]) - g_ND[b_i];
    b_eqs[3 * b_i] = eqs[X_tmp];
    b_eqs[3 * b_i + 1] = eqs[X_tmp + 1];
    b_eqs[3 * b_i + 2] = eqs[X_tmp + 2];
  }

  solve_3Q3(b_eqs, &scale, xs_data, xs_size, ys_data, ys_size, zs_data, zs_size);

  // we may choes another 3 out of 4 eq-s
  // 'solve_P4Pf:13' fs = zeros(1, 0, 'like', X);
  f_size[0] = 1;
  f_size[1] = 0;

  // 'solve_P4Pf:14' coder.varsize('fs', [1 10], [0 1]);
  // 'solve_P4Pf:15' Rs = zeros(3, 3, 0, 'like', X);
  r_size[0] = 3;
  r_size[1] = 3;
  r_size[2] = 0;

  // 'solve_P4Pf:16' coder.varsize('Rs', [3 3 10], [0 0 1]);
  // 'solve_P4Pf:17' Ts = zeros(3, 0, 'like', X);
  t_size[0] = 3;
  t_size[1] = 0;

  // 'solve_P4Pf:18' coder.varsize('Ts', [3 10], [0 1]);
  // 'solve_P4Pf:19' solution_num = int32(0);
  *n = 0;

  // 'solve_P4Pf:20' for i = 1 : n
  b_i = static_cast<int>(scale);
  for (i = 0; i < b_i; i++) {
    // 'solve_P4Pf:21' P1 = (N(1:3, :)*[xs(i); ys(i); zs(i); 1])';
    // 'solve_P4Pf:22' P3 = (D(1:3, :)*[xs(i); ys(i); zs(i); 1])';
    b_x = xs_data[i];
    d = ys_data[i];
    d1 = zs_data[i];
    for (b_B_tmp = 0; b_B_tmp < 3; b_B_tmp++) {
      P1[b_B_tmp] = ((Q[b_B_tmp + 32] * b_x + Q[b_B_tmp + 40] * d) + Q[b_B_tmp +
                     48] * d1) + Q[b_B_tmp + 56];
      P3[b_B_tmp] = ((C[b_B_tmp] * b_x + C[b_B_tmp + 4] * d) + C[b_B_tmp + 8] *
                     d1) + C[b_B_tmp + 12];
    }

    // 'solve_P4Pf:23' P21 = N(5, :)*[xs(i); ys(i); zs(i); 1];
    d = ys_data[i];
    d1 = zs_data[i];
    b_Q = ((Q[36] * xs_data[i] + Q[44] * ys_data[i]) + Q[52] * zs_data[i]) + Q
      [60];

    // 'solve_P4Pf:25' w = sqrt(sum(P3(1:3).^2)/sum(P1(1:3).^2));
    w = std::sqrt(((P3[0] * P3[0] + P3[1] * P3[1]) + P3[2] * P3[2]) / ((P1[0] *
      P1[0] + P1[1] * P1[1]) + P1[2] * P1[2]));

    // 'solve_P4Pf:27' alpha = norm(P1);
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
      alpha = alpha * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      alpha += t * t;
    }

    absxk = std::abs(P1[2]);
    if (absxk > scale) {
      t = scale / absxk;
      alpha = alpha * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      alpha += t * t;
    }

    alpha = scale * std::sqrt(alpha);

    // 'solve_P4Pf:28' R1 = P1/alpha;
    // 'solve_P4Pf:29' R3 = P3/(w*alpha);
    scale = w * alpha;
    for (b_B_tmp = 0; b_B_tmp < 3; b_B_tmp++) {
      P1[b_B_tmp] = (((Q[b_B_tmp + 32] * b_x + Q[b_B_tmp + 40] * d) + Q[b_B_tmp
                      + 48] * d1) + Q[b_B_tmp + 56]) / alpha;
      P3[b_B_tmp] = (((C[b_B_tmp] * b_x + C[b_B_tmp + 4] * d) + C[b_B_tmp + 8] *
                      d1) + C[b_B_tmp + 12]) / scale;
    }

    // 'solve_P4Pf:30' R2 = cross(R1, R3);
    scale = P1[1] * P3[2] - P1[2] * P3[1];
    absxk = P1[2] * P3[0] - P1[0] * P3[2];
    t = P1[0] * P3[1] - P1[1] * P3[0];

    // R2 = -R2/norm(R2);
    // 'solve_P4Pf:32' if sign(R2(1)) == sign(P21)
    b_x = scale;
    if (scale < 0.0) {
      b_x = -1.0;
    } else {
      if (scale > 0.0) {
        b_x = 1.0;
      }
    }

    if (b_Q < 0.0) {
      b_Q = -1.0;
    } else {
      if (b_Q > 0.0) {
        b_Q = 1.0;
      }
    }

    if (b_x == b_Q) {
      // 'solve_P4Pf:33' R = [R1; R2; R3];
      R[0] = P1[0];
      R[1] = scale;
      R[2] = P3[0];
      R[3] = P1[1];
      R[4] = absxk;
      R[5] = P3[1];
      R[6] = P1[2];
      R[7] = t;
      R[8] = P3[2];
    } else {
      // 'solve_P4Pf:34' else
      // 'solve_P4Pf:35' R = [R1; -R2; R3];
      R[0] = P1[0];
      R[1] = -scale;
      R[2] = P3[0];
      R[3] = P1[1];
      R[4] = -absxk;
      R[5] = P3[1];
      R[6] = P1[2];
      R[7] = -t;
      R[8] = P3[2];
    }

    // 'solve_P4Pf:37' if det(R) < 0
    std::memcpy(&c_x[0], &R[0], 9U * sizeof(double));
    c_xzgetrf(c_x, ipiv, &X_tmp);
    isodd = (ipiv[0] > 1);
    scale = c_x[0] * c_x[4] * c_x[8];
    if (ipiv[1] > 2) {
      isodd = !isodd;
    }

    if (isodd) {
      scale = -scale;
    }

    if (scale < 0.0) {
      // 'solve_P4Pf:38' R = -R;
      for (b_B_tmp = 0; b_B_tmp < 9; b_B_tmp++) {
        R[b_B_tmp] = -R[b_B_tmp];
      }
    }

    // 'solve_P4Pf:40' T_ = find_T(X(1:3, :), u, v, R, w);
    // 'find_T:2' A = zeros(12, 3, 'double');
    std::memset(&c_A[0], 0, 36U * sizeof(double));

    // 'find_T:3' b = zeros(12, 1, 'double');
    std::memset(&b[0], 0, 12U * sizeof(double));

    // 'find_T:5' for i = 1 : 4
    for (B_tmp = 0; B_tmp < 4; B_tmp++) {
      // 'find_T:6' A(3*i - 2:3*i, :) = [     0,   -1,  w*v(i);
      // 'find_T:7'                                   1,    0, -w*u(i);
      // 'find_T:8'                               -v(i), u(i),      0];
      X_tmp = 3 * (B_tmp + 1);
      P3[0] = static_cast<double>(X_tmp) + -2.0;
      P3[1] = static_cast<double>(X_tmp) + -1.0;
      P3[2] = X_tmp;
      offset = X_tmp + -2;
      c_A[offset - 1] = 0.0;
      c_A[offset + 11] = -1.0;
      c_A[offset + 23] = w * y[B_tmp];
      offset = X_tmp + -1;
      c_A[offset - 1] = 1.0;
      c_A[offset + 11] = 0.0;
      c_A[offset + 23] = -w * x[B_tmp];
      c_A[X_tmp - 1] = -y[B_tmp];
      c_A[X_tmp + 11] = x[B_tmp];
      c_A[X_tmp + 23] = 0.0;

      // 'find_T:9' b(3*i - 2:3*i) = -A(3*i - 2:3*i, :)*R*X(:, i);
      b_B_tmp = X_tmp - 3;
      for (offset = 0; offset < 3; offset++) {
        X_tmp = b_B_tmp + 12 * offset;
        c_x[3 * offset] = -c_A[X_tmp];
        c_x[3 * offset + 1] = -c_A[X_tmp + 1];
        c_x[3 * offset + 2] = -c_A[X_tmp + 2];
      }

      for (b_B_tmp = 0; b_B_tmp < 3; b_B_tmp++) {
        b_x = 0.0;
        d = c_x[b_B_tmp + 3];
        d1 = c_x[b_B_tmp + 6];
        for (offset = 0; offset < 3; offset++) {
          b_x += ((c_x[b_B_tmp] * R[3 * offset] + d * R[3 * offset + 1]) + d1 *
                  R[3 * offset + 2]) * b_X[offset + (B_tmp << 2)];
        }

        b[static_cast<int>(P3[b_B_tmp]) - 1] = b_x;
      }
    }

    // 'find_T:11' T = A\b;
    mldivide(c_A, b, c_Q);
    P3[0] = c_Q[0];
    P3[1] = c_Q[1];
    P3[2] = c_Q[2];

    // 'solve_P4Pf:41' if isa(X,'single')
    // 'solve_P4Pf:43' else
    // 'solve_P4Pf:44' T = T_;
    // 'solve_P4Pf:47' P = diag([1, 1, w])*[R, T];
    P1[0] = 1.0;
    P1[1] = 1.0;
    P1[2] = w;
    std::memset(&c_x[0], 0, 9U * sizeof(double));

    // 'solve_P4Pf:48' U_eval = P*X;
    for (j = 0; j < 3; j++) {
      c_x[j + 3 * j] = P1[j];
      b[3 * j] = R[3 * j];
      X_tmp = 3 * j + 1;
      b[X_tmp] = R[X_tmp];
      X_tmp = 3 * j + 2;
      b[X_tmp] = R[X_tmp];
      b[j + 9] = P3[j];
    }

    for (b_B_tmp = 0; b_B_tmp < 3; b_B_tmp++) {
      b_x = c_x[b_B_tmp + 3];
      d = c_x[b_B_tmp + 6];
      for (offset = 0; offset < 4; offset++) {
        d_x[b_B_tmp + 3 * offset] = (c_x[b_B_tmp] * b[3 * offset] + b_x * b[3 *
          offset + 1]) + d * b[3 * offset + 2];
      }
    }

    for (b_B_tmp = 0; b_B_tmp < 3; b_B_tmp++) {
      b_x = d_x[b_B_tmp + 3];
      d = d_x[b_B_tmp + 6];
      d1 = d_x[b_B_tmp + 9];
      for (offset = 0; offset < 4; offset++) {
        X_tmp = offset << 2;
        b[b_B_tmp + 3 * offset] = ((d_x[b_B_tmp] * b_X[X_tmp] + b_x * b_X[X_tmp
          + 1]) + d * b_X[X_tmp + 2]) + d1 * b_X[X_tmp + 3];
      }
    }

    // 'solve_P4Pf:49' for j = 1 : 4
    for (j = 0; j < 4; j++) {
      // 'solve_P4Pf:50' U_eval(:, j) = U_eval(:, j) ./U_eval(3, j);
      b_B_tmp = 3 * j + 2;
      X_tmp = 3 * j + 1;
      c_Q[1] = b[X_tmp] / b[b_B_tmp];
      c_Q[2] = b[b_B_tmp] / b[b_B_tmp];
      b[3 * j] /= b[b_B_tmp];
      b[X_tmp] = c_Q[1];
      b[b_B_tmp] = c_Q[2];
    }

    // U_eval = U_eval ./ U_eval(3, :);
    // U_eval = bsxfun(@rdivide, U_eval, U_eval(3, :));
    // reprojection error check
    // 'solve_P4Pf:55' if norm(U_eval(1:2, :) - [u;v])*w < 0.01
    b_b[0] = b[0] - x[0];
    b_b[1] = b[1] - y[0];
    b_b[2] = b[3] - x[1];
    b_b[3] = b[4] - y[1];
    b_b[4] = b[6] - x[2];
    b_b[5] = b[7] - y[2];
    b_b[6] = b[9] - x[3];
    b_b[7] = b[10] - y[3];
    if (b_norm(b_b) * w < 0.01) {
      // 'solve_P4Pf:56' solution_num = solution_num + 1;
      (*n)++;

      // 'solve_P4Pf:57' fs = [fs, 1/w];
      b_B_tmp = f_size[1];
      f_size[1]++;
      f_data[b_B_tmp] = 1.0 / w;

      // 'solve_P4Pf:58' if solution_num == 1
      if (*n == 1) {
        // 'solve_P4Pf:59' Rs = R;
        r_size[0] = 3;
        r_size[1] = 3;
        r_size[2] = 1;
        std::memcpy(&r_data[0], &R[0], 9U * sizeof(double));
      } else {
        // 'solve_P4Pf:60' else
        // 'solve_P4Pf:61' Rs = cat(3, Rs, R);
        X_tmp = static_cast<signed char>((r_size[2] + 1));
        offset = -1;
        b_B_tmp = 9 * r_size[2];
        for (j = 0; j < b_B_tmp; j++) {
          offset++;
          y_data[offset] = r_data[j];
        }

        for (j = 0; j < 9; j++) {
          offset++;
          y_data[offset] = R[j];
        }

        r_size[0] = 3;
        r_size[1] = 3;
        r_size[2] = X_tmp;
        X_tmp *= 9;
        if (0 <= X_tmp - 1) {
          std::memcpy(&r_data[0], &y_data[0], X_tmp * sizeof(double));
        }
      }

      // 'solve_P4Pf:63' Ts = [Ts, -R'*T];
      for (b_B_tmp = 0; b_B_tmp < 3; b_B_tmp++) {
        c_x[3 * b_B_tmp] = -R[b_B_tmp];
        c_x[3 * b_B_tmp + 1] = -R[b_B_tmp + 3];
        c_x[3 * b_B_tmp + 2] = -R[b_B_tmp + 6];
      }

      for (b_B_tmp = 0; b_B_tmp < 3; b_B_tmp++) {
        c_Q[b_B_tmp] = (c_x[b_B_tmp] * P3[0] + c_x[b_B_tmp + 3] * P3[1]) +
          c_x[b_B_tmp + 6] * P3[2];
      }

      X_tmp = t_size[1];
      offset = t_size[1] + 1;
      for (b_B_tmp = 0; b_B_tmp < X_tmp; b_B_tmp++) {
        b_eqs[3 * b_B_tmp] = t_data[3 * b_B_tmp];
        B_tmp = 3 * b_B_tmp + 1;
        b_eqs[B_tmp] = t_data[B_tmp];
        B_tmp = 3 * b_B_tmp + 2;
        b_eqs[B_tmp] = t_data[B_tmp];
      }

      b_eqs[3 * t_size[1]] = c_Q[0];
      b_eqs[3 * t_size[1] + 1] = c_Q[1];
      b_eqs[3 * t_size[1] + 2] = c_Q[2];
      t_size[0] = 3;
      t_size[1] = offset;
      X_tmp = 3 * offset;
      if (0 <= X_tmp - 1) {
        std::memcpy(&t_data[0], &b_eqs[0], X_tmp * sizeof(double));
      }
    }
  }

  // 'p4pf_double:8' assert(isa(n, 'int32'));
  // 'p4pf_double:9' assert(isa(f, 'double'));
  // 'p4pf_double:10' assert(isa(r, 'double'));
  // 'p4pf_double:11' assert(isa(t, 'double'));
}

//
// function [n,f,r,t] = p4pf_single(X, x, y, e)
// Arguments    : const float X[12]
//                const float x[4]
//                const float y[4]
//                float e
//                int *n
//                float f_data[]
//                int f_size[2]
//                float r_data[]
//                int r_size[3]
//                float t_data[]
//                int t_size[2]
// Return Type  : void
//
void p4pf_single(const float X[12], const float x[4], const float y[4], float e,
                 int *n, float f_data[], int f_size[2], float r_data[], int
                 r_size[3], float t_data[], int t_size[2])
{
  int i;
  float A[32];
  float Q[64];
  float b_A[32];
  int X_tmp;
  float b_x;
  float b_X[16];
  float f;
  float f1;
  double B[16];
  double C[16];
  float scale;
  int b_i;
  int offset;
  float D[16];
  int B_tmp;
  int b_B_tmp;
  int j;
  float eqs[40];
  float ND[48];
  float b_ND[10];
  float c_ND[10];
  float d_ND[10];
  float e_ND[10];
  float f_ND[10];
  float g_ND[10];
  float b_eqs[30];
  double b_n;
  float xs_data[13];
  int xs_size[2];
  float ys_data[13];
  int ys_size[2];
  float zs_data[13];
  int zs_size[2];
  float P1[3];
  float P3[3];
  float b_Q;
  float w;
  float absxk;
  float t;
  float alpha;
  float R[9];
  float c_x[9];
  int ipiv[3];
  boolean_T isodd;
  double c_A[36];
  double b[12];
  double A_tmp[3];
  float U_eval[12];
  float d_x[12];
  float b_U_eval[8];
  float y_data[99];

  // 'p4pf_single:2' assert(isa(X, 'single'));
  // 'p4pf_single:3' assert(isa(x, 'single'));
  // 'p4pf_single:4' assert(isa(y, 'single'));
  // 'p4pf_single:5' assert(isa(e, 'single'));
  //  TODO: rename solve_P4Pf to match p3.5p
  // 'p4pf_single:8' [n, f, r, t] = solve_P4Pf(X, x, y, e);
  // SOLVE_P4Pf Summary of this function goes here
  //        X = [p1, p2, p3, p4], pi = [4, 1]; X(:, i) <-> (u(i), v(i))
  //        if f is a correct foal length, then [R, T] = [R, T] / sign(d)*abs(d)^(1/3); 
  //        where d = det(R)
  // 'solve_P4Pf:6' X = [X; ones(1, 4, 'like', X)];
  // 'solve_P4Pf:7' A = find_A(X, u, v);
  // [p11, p12, p13, p14, p21, p22, p23, p24]'
  // 'solve_P4Pf:70' A = zeros(4, 8, 'like', X);
  // 'solve_P4Pf:71' for i = 1 : 4
  // 'solve_P4Pf:8' [Q, ~] = qr(A');
  for (i = 0; i < 4; i++) {
    X_tmp = i << 2;
    b_x = X[3 * i];
    b_X[X_tmp] = b_x;
    f = X[3 * i + 1];
    b_X[X_tmp + 1] = f;
    f1 = X[3 * i + 2];
    b_X[X_tmp + 2] = f1;
    b_X[X_tmp + 3] = 1.0F;

    // 'solve_P4Pf:72' A(i,  :) = [-v(i)*X(:, i)', u(i)*X(:, i)'];
    b_A[i] = -y[i] * b_x;
    b_A[i + 16] = x[i] * b_x;
    b_A[i + 4] = -y[i] * f;
    b_A[i + 20] = x[i] * f;
    b_A[i + 8] = -y[i] * f1;
    b_A[i + 24] = x[i] * f1;
    b_A[i + 12] = -y[i];
    b_A[i + 28] = x[i];
    for (b_i = 0; b_i < 8; b_i++) {
      A[b_i + (i << 3)] = b_A[i + (b_i << 2)];
    }
  }

  e_qr(A, Q, b_A);

  // 'solve_P4Pf:8' ~
  // 'solve_P4Pf:9' N = Q(:, 5:end);
  // nullspace
  // 'solve_P4Pf:10' D = find_D(X, u, v, N, e);
  // 'solve_P4Pf:77' B = zeros(4, 'double');
  // 'solve_P4Pf:78' C = zeros(4, 'double');
  // 'solve_P4Pf:79' for i = 1 : 4
  for (i = 0; i < 4; i++) {
    // 'solve_P4Pf:80' if abs(u(i)) < e
    if (std::abs(x[i]) < e) {
      // 'solve_P4Pf:81' fctr = v(i);
      scale = y[i];

      // 'solve_P4Pf:82' offset = 4;
      offset = 4;
    } else {
      // 'solve_P4Pf:83' else
      // 'solve_P4Pf:84' fctr = u(i);
      scale = x[i];

      // 'solve_P4Pf:85' offset = 0;
      offset = 0;
    }

    // 'solve_P4Pf:87' B(i, :) = fctr*X(:, i)';
    // 'solve_P4Pf:88' for j = 1 : 4
    b_B_tmp = i << 2;
    for (j = 0; j < 4; j++) {
      B_tmp = i + (j << 2);
      B[B_tmp] = scale * b_X[j + b_B_tmp];

      // 'solve_P4Pf:89' C(i, j) = X(:, i)'*ns((1:4)+offset, j);
      X_tmp = offset + ((j + 4) << 3);
      C[B_tmp] = ((b_X[b_B_tmp] * Q[X_tmp] + b_X[b_B_tmp + 1] * Q[X_tmp + 1]) +
                  b_X[b_B_tmp + 2] * Q[X_tmp + 2]) + b_X[b_B_tmp + 3] * Q[X_tmp
        + 3];
    }
  }

  // 'solve_P4Pf:92' D_ = B\C;
  // 'solve_P4Pf:93' if isa(X, 'single')
  // 'solve_P4Pf:94' D = single(D_);
  d_mldivide(B, C);
  for (b_i = 0; b_i < 16; b_i++) {
    D[b_i] = static_cast<float>(C[b_i]);
  }

  // 'solve_P4Pf:11' eqs = find_eqs([N; D]);
  for (b_i = 0; b_i < 4; b_i++) {
    for (B_tmp = 0; B_tmp < 8; B_tmp++) {
      ND[B_tmp + 12 * b_i] = Q[B_tmp + ((b_i + 4) << 3)];
    }

    X_tmp = b_i << 2;
    ND[12 * b_i + 8] = D[X_tmp];
    ND[12 * b_i + 9] = D[X_tmp + 1];
    ND[12 * b_i + 10] = D[X_tmp + 2];
    ND[12 * b_i + 11] = D[X_tmp + 3];
  }

  // ND = [N; D]
  // 'find_eqs:2' eqs = zeros(4, 10, 'like', ND);
  std::memset(&eqs[0], 0, 40U * sizeof(float));

  // 'find_eqs:3' eqs(1, :) = mult_pp(1, 1, 2, 1, ND) + mult_pp(1, 2, 2, 2, ND) + mult_pp(1, 3, 2, 3, ND); 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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
  for (b_i = 0; b_i < 10; b_i++) {
    eqs[b_i << 2] = (b_ND[b_i] + c_ND[b_i]) + d_ND[b_i];
  }

  // 'find_eqs:4' eqs(2, :) = mult_pp(3, 1, 1, 1, ND) + mult_pp(3, 2, 1, 2, ND) + mult_pp(3, 3, 1, 3, ND); 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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
  for (b_i = 0; b_i < 10; b_i++) {
    eqs[(b_i << 2) + 1] = (b_ND[b_i] + c_ND[b_i]) + d_ND[b_i];
  }

  // 'find_eqs:5' eqs(3, :) = mult_pp(3, 1, 2, 1, ND) + mult_pp(3, 2, 2, 2, ND) + mult_pp(3, 3, 2, 3, ND); 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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
  for (b_i = 0; b_i < 10; b_i++) {
    eqs[(b_i << 2) + 2] = (b_ND[b_i] + c_ND[b_i]) + d_ND[b_i];
  }

  // 'find_eqs:6' eqs(4, :) = sq_pp(1, 1, ND) + sq_pp(1, 2, ND) + sq_pp(1, 3, ND) - ... 
  // 'find_eqs:7'         sq_pp(2, 1, ND) - sq_pp(2, 2, ND) - sq_pp(2, 3, ND);
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
  // 'find_eqs:19' prod = mult_pp(i, j, i, j, ND);
  //  we need coeffs for [x^2, y^2, z^2, xy, xz, yz, x, y, z, 1)]
  // 'find_eqs:12' a = ND(4*(i - 1) + j, :);
  // 'find_eqs:13' b = ND(4*(m - 1) + k, :);
  // 'find_eqs:14' prod = [a(1)*b(1), a(2)*b(2), a(3)*b(3), a(1)*b(2) + a(2)*b(1), a(1)*b(3) + a(3)*b(1), a(2)*b(3) + a(3)*b(2), a(1)*b(4) + a(4)*b(1), a(2)*b(4) + a(4)*b(2), a(3)*b(4) + a(4)*b(3), a(4)*b(4)]; 
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

  // 'solve_P4Pf:12' [n, xs, ys, zs] = solve_3Q3(eqs(1:3, :), e);
  for (b_i = 0; b_i < 10; b_i++) {
    X_tmp = b_i << 2;
    eqs[X_tmp + 3] = ((((b_ND[b_i] + c_ND[b_i]) + d_ND[b_i]) - e_ND[b_i]) -
                      f_ND[b_i]) - g_ND[b_i];
    b_eqs[3 * b_i] = eqs[X_tmp];
    b_eqs[3 * b_i + 1] = eqs[X_tmp + 1];
    b_eqs[3 * b_i + 2] = eqs[X_tmp + 2];
  }

  b_solve_3Q3(b_eqs, &b_n, xs_data, xs_size, ys_data, ys_size, zs_data, zs_size);

  // we may choes another 3 out of 4 eq-s
  // 'solve_P4Pf:13' fs = zeros(1, 0, 'like', X);
  f_size[0] = 1;
  f_size[1] = 0;

  // 'solve_P4Pf:14' coder.varsize('fs', [1 10], [0 1]);
  // 'solve_P4Pf:15' Rs = zeros(3, 3, 0, 'like', X);
  r_size[0] = 3;
  r_size[1] = 3;
  r_size[2] = 0;

  // 'solve_P4Pf:16' coder.varsize('Rs', [3 3 10], [0 0 1]);
  // 'solve_P4Pf:17' Ts = zeros(3, 0, 'like', X);
  t_size[0] = 3;
  t_size[1] = 0;

  // 'solve_P4Pf:18' coder.varsize('Ts', [3 10], [0 1]);
  // 'solve_P4Pf:19' solution_num = int32(0);
  *n = 0;

  // 'solve_P4Pf:20' for i = 1 : n
  b_i = static_cast<int>(b_n);
  for (i = 0; i < b_i; i++) {
    // 'solve_P4Pf:21' P1 = (N(1:3, :)*[xs(i); ys(i); zs(i); 1])';
    // 'solve_P4Pf:22' P3 = (D(1:3, :)*[xs(i); ys(i); zs(i); 1])';
    b_x = xs_data[i];
    f = ys_data[i];
    f1 = zs_data[i];
    for (B_tmp = 0; B_tmp < 3; B_tmp++) {
      P1[B_tmp] = ((Q[B_tmp + 32] * b_x + Q[B_tmp + 40] * f) + Q[B_tmp + 48] *
                   f1) + Q[B_tmp + 56];
      P3[B_tmp] = ((D[B_tmp] * b_x + D[B_tmp + 4] * f) + D[B_tmp + 8] * f1) +
        D[B_tmp + 12];
    }

    // 'solve_P4Pf:23' P21 = N(5, :)*[xs(i); ys(i); zs(i); 1];
    f = ys_data[i];
    f1 = zs_data[i];
    b_Q = ((Q[36] * xs_data[i] + Q[44] * ys_data[i]) + Q[52] * zs_data[i]) + Q
      [60];

    // 'solve_P4Pf:25' w = sqrt(sum(P3(1:3).^2)/sum(P1(1:3).^2));
    w = std::sqrt(((P3[0] * P3[0] + P3[1] * P3[1]) + P3[2] * P3[2]) / ((P1[0] *
      P1[0] + P1[1] * P1[1]) + P1[2] * P1[2]));

    // 'solve_P4Pf:27' alpha = norm(P1);
    scale = 1.29246971E-26F;
    absxk = std::abs(P1[0]);
    if (absxk > 1.29246971E-26F) {
      alpha = 1.0F;
      scale = absxk;
    } else {
      t = absxk / 1.29246971E-26F;
      alpha = t * t;
    }

    absxk = std::abs(P1[1]);
    if (absxk > scale) {
      t = scale / absxk;
      alpha = alpha * t * t + 1.0F;
      scale = absxk;
    } else {
      t = absxk / scale;
      alpha += t * t;
    }

    absxk = std::abs(P1[2]);
    if (absxk > scale) {
      t = scale / absxk;
      alpha = alpha * t * t + 1.0F;
      scale = absxk;
    } else {
      t = absxk / scale;
      alpha += t * t;
    }

    alpha = scale * std::sqrt(alpha);

    // 'solve_P4Pf:28' R1 = P1/alpha;
    // 'solve_P4Pf:29' R3 = P3/(w*alpha);
    scale = w * alpha;
    for (B_tmp = 0; B_tmp < 3; B_tmp++) {
      P1[B_tmp] = (((Q[B_tmp + 32] * b_x + Q[B_tmp + 40] * f) + Q[B_tmp + 48] *
                    f1) + Q[B_tmp + 56]) / alpha;
      P3[B_tmp] = (((D[B_tmp] * b_x + D[B_tmp + 4] * f) + D[B_tmp + 8] * f1) +
                   D[B_tmp + 12]) / scale;
    }

    // 'solve_P4Pf:30' R2 = cross(R1, R3);
    scale = P1[1] * P3[2] - P1[2] * P3[1];
    absxk = P1[2] * P3[0] - P1[0] * P3[2];
    t = P1[0] * P3[1] - P1[1] * P3[0];

    // R2 = -R2/norm(R2);
    // 'solve_P4Pf:32' if sign(R2(1)) == sign(P21)
    b_x = scale;
    if (scale < 0.0F) {
      b_x = -1.0F;
    } else {
      if (scale > 0.0F) {
        b_x = 1.0F;
      }
    }

    if (b_Q < 0.0F) {
      b_Q = -1.0F;
    } else {
      if (b_Q > 0.0F) {
        b_Q = 1.0F;
      }
    }

    if (b_x == b_Q) {
      // 'solve_P4Pf:33' R = [R1; R2; R3];
      R[0] = P1[0];
      R[1] = scale;
      R[2] = P3[0];
      R[3] = P1[1];
      R[4] = absxk;
      R[5] = P3[1];
      R[6] = P1[2];
      R[7] = t;
      R[8] = P3[2];
    } else {
      // 'solve_P4Pf:34' else
      // 'solve_P4Pf:35' R = [R1; -R2; R3];
      R[0] = P1[0];
      R[1] = -scale;
      R[2] = P3[0];
      R[3] = P1[1];
      R[4] = -absxk;
      R[5] = P3[1];
      R[6] = P1[2];
      R[7] = -t;
      R[8] = P3[2];
    }

    // 'solve_P4Pf:37' if det(R) < 0
    for (B_tmp = 0; B_tmp < 9; B_tmp++) {
      c_x[B_tmp] = R[B_tmp];
    }

    d_xzgetrf(c_x, ipiv, &X_tmp);
    isodd = (ipiv[0] > 1);
    scale = c_x[0] * c_x[4] * c_x[8];
    if (ipiv[1] > 2) {
      isodd = !isodd;
    }

    if (isodd) {
      scale = -scale;
    }

    if (scale < 0.0F) {
      // 'solve_P4Pf:38' R = -R;
      for (B_tmp = 0; B_tmp < 9; B_tmp++) {
        R[B_tmp] = -R[B_tmp];
      }
    }

    // 'solve_P4Pf:40' T_ = find_T(X(1:3, :), u, v, R, w);
    // 'find_T:2' A = zeros(12, 3, 'double');
    std::memset(&c_A[0], 0, 36U * sizeof(double));

    // 'find_T:3' b = zeros(12, 1, 'double');
    std::memset(&b[0], 0, 12U * sizeof(double));

    // 'find_T:5' for i = 1 : 4
    for (b_B_tmp = 0; b_B_tmp < 4; b_B_tmp++) {
      // 'find_T:6' A(3*i - 2:3*i, :) = [     0,   -1,  w*v(i);
      // 'find_T:7'                                   1,    0, -w*u(i);
      // 'find_T:8'                               -v(i), u(i),      0];
      X_tmp = 3 * (b_B_tmp + 1);
      A_tmp[0] = static_cast<double>(X_tmp) + -2.0;
      A_tmp[1] = static_cast<double>(X_tmp) + -1.0;
      A_tmp[2] = X_tmp;
      offset = X_tmp + -2;
      c_A[offset - 1] = 0.0;
      c_A[offset + 11] = -1.0;
      c_A[offset + 23] = w * y[b_B_tmp];
      offset = X_tmp + -1;
      c_A[offset - 1] = 1.0;
      c_A[offset + 11] = 0.0;
      c_A[offset + 23] = -w * x[b_B_tmp];
      c_A[X_tmp - 1] = -y[b_B_tmp];
      c_A[X_tmp + 11] = x[b_B_tmp];
      c_A[X_tmp + 23] = 0.0;

      // 'find_T:9' b(3*i - 2:3*i) = -A(3*i - 2:3*i, :)*R*X(:, i);
      B_tmp = X_tmp - 3;
      for (offset = 0; offset < 3; offset++) {
        X_tmp = B_tmp + 12 * offset;
        c_x[3 * offset] = static_cast<float>(-c_A[X_tmp]);
        c_x[3 * offset + 1] = static_cast<float>(-c_A[X_tmp + 1]);
        c_x[3 * offset + 2] = static_cast<float>(-c_A[X_tmp + 2]);
      }

      for (B_tmp = 0; B_tmp < 3; B_tmp++) {
        b_x = 0.0F;
        f = c_x[B_tmp + 3];
        f1 = c_x[B_tmp + 6];
        for (offset = 0; offset < 3; offset++) {
          b_x += ((c_x[B_tmp] * R[3 * offset] + f * R[3 * offset + 1]) + f1 * R
                  [3 * offset + 2]) * b_X[offset + (b_B_tmp << 2)];
        }

        b[static_cast<int>(A_tmp[B_tmp]) - 1] = b_x;
      }
    }

    // 'find_T:11' T = A\b;
    // 'solve_P4Pf:41' if isa(X,'single')
    // 'solve_P4Pf:42' T = single(T_);
    mldivide(c_A, b, A_tmp);
    P1[0] = static_cast<float>(A_tmp[0]);
    P1[1] = static_cast<float>(A_tmp[1]);
    P1[2] = static_cast<float>(A_tmp[2]);

    // 'solve_P4Pf:47' P = diag([1, 1, w])*[R, T];
    P3[0] = 1.0F;
    P3[1] = 1.0F;
    P3[2] = w;
    for (B_tmp = 0; B_tmp < 9; B_tmp++) {
      c_x[B_tmp] = 0.0F;
    }

    // 'solve_P4Pf:48' U_eval = P*X;
    for (j = 0; j < 3; j++) {
      c_x[j + 3 * j] = P3[j];
      U_eval[3 * j] = R[3 * j];
      X_tmp = 3 * j + 1;
      U_eval[X_tmp] = R[X_tmp];
      X_tmp = 3 * j + 2;
      U_eval[X_tmp] = R[X_tmp];
      U_eval[j + 9] = P1[j];
    }

    for (B_tmp = 0; B_tmp < 3; B_tmp++) {
      b_x = c_x[B_tmp + 3];
      f = c_x[B_tmp + 6];
      for (offset = 0; offset < 4; offset++) {
        d_x[B_tmp + 3 * offset] = (c_x[B_tmp] * U_eval[3 * offset] + b_x *
          U_eval[3 * offset + 1]) + f * U_eval[3 * offset + 2];
      }
    }

    for (B_tmp = 0; B_tmp < 3; B_tmp++) {
      b_x = d_x[B_tmp + 3];
      f = d_x[B_tmp + 6];
      f1 = d_x[B_tmp + 9];
      for (offset = 0; offset < 4; offset++) {
        X_tmp = offset << 2;
        U_eval[B_tmp + 3 * offset] = ((d_x[B_tmp] * b_X[X_tmp] + b_x * b_X[X_tmp
          + 1]) + f * b_X[X_tmp + 2]) + f1 * b_X[X_tmp + 3];
      }
    }

    // 'solve_P4Pf:49' for j = 1 : 4
    // U_eval = U_eval ./ U_eval(3, :);
    // U_eval = bsxfun(@rdivide, U_eval, U_eval(3, :));
    // reprojection error check
    // 'solve_P4Pf:55' if norm(U_eval(1:2, :) - [u;v])*w < 0.01
    for (j = 0; j < 4; j++) {
      // 'solve_P4Pf:50' U_eval(:, j) = U_eval(:, j) ./U_eval(3, j);
      B_tmp = 3 * j + 2;
      b_x = U_eval[B_tmp];
      f = U_eval[3 * j] / U_eval[B_tmp];
      U_eval[3 * j] = f;
      X_tmp = 3 * j + 1;
      f1 = U_eval[X_tmp] / b_x;
      U_eval[X_tmp] = f1;
      b_x /= b_x;
      U_eval[B_tmp] = b_x;
      X_tmp = j << 1;
      b_U_eval[X_tmp] = f - x[j];
      b_U_eval[X_tmp + 1] = f1 - y[j];
    }

    if (c_norm(b_U_eval) * w < 0.01) {
      // 'solve_P4Pf:56' solution_num = solution_num + 1;
      (*n)++;

      // 'solve_P4Pf:57' fs = [fs, 1/w];
      B_tmp = f_size[1];
      f_size[1]++;
      f_data[B_tmp] = 1.0F / w;

      // 'solve_P4Pf:58' if solution_num == 1
      if (*n == 1) {
        // 'solve_P4Pf:59' Rs = R;
        r_size[0] = 3;
        r_size[1] = 3;
        r_size[2] = 1;
        for (B_tmp = 0; B_tmp < 9; B_tmp++) {
          r_data[B_tmp] = R[B_tmp];
        }
      } else {
        // 'solve_P4Pf:60' else
        // 'solve_P4Pf:61' Rs = cat(3, Rs, R);
        X_tmp = static_cast<signed char>((r_size[2] + 1));
        offset = -1;
        B_tmp = 9 * r_size[2];
        for (j = 0; j < B_tmp; j++) {
          offset++;
          y_data[offset] = r_data[j];
        }

        for (j = 0; j < 9; j++) {
          offset++;
          y_data[offset] = R[j];
        }

        r_size[0] = 3;
        r_size[1] = 3;
        r_size[2] = X_tmp;
        X_tmp *= 9;
        if (0 <= X_tmp - 1) {
          std::memcpy(&r_data[0], &y_data[0], X_tmp * sizeof(float));
        }
      }

      // 'solve_P4Pf:63' Ts = [Ts, -R'*T];
      for (B_tmp = 0; B_tmp < 3; B_tmp++) {
        c_x[3 * B_tmp] = -R[B_tmp];
        c_x[3 * B_tmp + 1] = -R[B_tmp + 3];
        c_x[3 * B_tmp + 2] = -R[B_tmp + 6];
      }

      for (B_tmp = 0; B_tmp < 3; B_tmp++) {
        P3[B_tmp] = (c_x[B_tmp] * P1[0] + c_x[B_tmp + 3] * P1[1]) + c_x[B_tmp +
          6] * P1[2];
      }

      X_tmp = t_size[1];
      offset = t_size[1] + 1;
      for (B_tmp = 0; B_tmp < X_tmp; B_tmp++) {
        b_eqs[3 * B_tmp] = t_data[3 * B_tmp];
        b_B_tmp = 3 * B_tmp + 1;
        b_eqs[b_B_tmp] = t_data[b_B_tmp];
        b_B_tmp = 3 * B_tmp + 2;
        b_eqs[b_B_tmp] = t_data[b_B_tmp];
      }

      b_eqs[3 * t_size[1]] = P3[0];
      b_eqs[3 * t_size[1] + 1] = P3[1];
      b_eqs[3 * t_size[1] + 2] = P3[2];
      t_size[0] = 3;
      t_size[1] = offset;
      X_tmp = 3 * offset;
      if (0 <= X_tmp - 1) {
        std::memcpy(&t_data[0], &b_eqs[0], X_tmp * sizeof(float));
      }
    }
  }

  // 'p4pf_single:9' assert(isa(n, 'int32'));
  // 'p4pf_single:10' assert(isa(t, 'single'));
  // 'p4pf_single:11' assert(isa(r, 'single'));
  // 'p4pf_single:12' assert(isa(f, 'single'));
}

//
// Arguments    : void
// Return Type  : void
//
void pnpf_initialize()
{
}

//
// Arguments    : void
// Return Type  : void
//
void pnpf_terminate()
{
  // (no terminate code required)
}

//
// File trailer for pnpf.cpp
//
// [EOF]
//
