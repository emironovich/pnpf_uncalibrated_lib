//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: p35p_solver.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Oct-2019 04:11:40
//

// Include Files
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Variable Definitions
static boolean_T isInitialized_p35p_solver_double = false;

// Function Declarations
static void b_sqrt(creal_T *x);
static double b_xnrm2(int n, const double x[3]);
static void b_xrot(double x[100], int ix0, int iy0, double c, double s);
static void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn);
static double c_xnrm2(int n, const double x[12], int ix0);
static void cat(const double varargin_1_data[], const int varargin_1_size[3],
                const double varargin_2[9], double y_data[], int y_size[3]);
static void eig(const double A[100], creal_T V[100], creal_T D[100]);
static int eml_dlahqr(double h[100], double z[100]);
static void equations_for_groebner(const double F[72], double G[112]);
static void find_f(const double F[72], double x, double y, double e, double
                   fc_data[], int fc_size[2], double fs_data[], int fs_size[2],
                   double *n);
static void mldivide(const double A[400], double B[200]);
static void mult_for_groebner(const double G4[112], double G20[900]);
static void mult_poly22(const double a[6], const double b[6], double c[15]);
static void mult_poly42(const double d[15], const double a[6], double p[28]);
static void qr(const double A[12], double Q[9], double R[12]);
static void quadruple_constraint(double i, double j, double k, const double x[4],
  const double y[4], const double X[12], const double R[54], double F_row[18]);
static double rt_hypotd_snf(double u0, double u1);
static void schur(const double A[100], double V[100], double T[100]);
static void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn);
static void xgerc(int m, int n, double alpha1, int ix0, const double y[4],
                  double A[12], int ia0);
static double xnrm2(int n, const double x[100], int ix0);
static void xrot(int n, double x[100], int ix0, int iy0, double c, double s);
static void xzggbak(creal_T V[100], int ilo, int ihi, const int rscale[10]);
static void xzggbal(creal_T A[100], int *ilo, int *ihi, int rscale[10]);
static void xzggev(creal_T A[100], int *info, creal_T alpha1[10], creal_T beta1
                   [10], creal_T V[100]);
static void xzgghrd(int ilo, int ihi, creal_T A[100], creal_T Z[100]);
static void xzhgeqz(creal_T A[100], int ilo, int ihi, creal_T Z[100], int *info,
                    creal_T alpha1[10], creal_T beta1[10]);
static void xzlarf(int m, int n, int iv0, double tau, double C[100], int ic0,
                   double work[10]);
static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r);
static void xztgevc(const creal_T A[100], creal_T V[100]);

// Function Definitions

//
// Arguments    : creal_T *x
// Return Type  : void
//
static void b_sqrt(creal_T *x)
{
  double xr;
  double xi;
  double absxi;
  double absxr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = std::sqrt(-xr);
    } else {
      absxi = std::sqrt(xr);
      xr = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxi = std::sqrt(-xi / 2.0);
      xr = -absxi;
    } else {
      absxi = std::sqrt(xi / 2.0);
      xr = absxi;
    }
  } else if (rtIsNaN(xr)) {
    absxi = xr;
  } else if (rtIsNaN(xi)) {
    absxi = xi;
    xr = xi;
  } else if (rtIsInf(xi)) {
    absxi = std::abs(xi);
    xr = xi;
  } else if (rtIsInf(xr)) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = xi * -xr;
    } else {
      absxi = xr;
      xr = 0.0;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi * 0.5);
      if (absxi > absxr) {
        absxi = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0);
      } else {
        absxi = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxi = std::sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (xr > 0.0) {
      xr = 0.5 * (xi / absxi);
    } else {
      if (xi < 0.0) {
        xr = -absxi;
      } else {
        xr = absxi;
      }

      absxi = 0.5 * (xi / xr);
    }
  }

  x->re = absxi;
  x->im = xr;
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
        g2 = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / g2;
        sn->im = -gs_im / g2;
      } else {
        g2s = std::sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        if (scale_tmp > 1.0) {
          g2 = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / g2;
          fs_im = f.im / g2;
        } else {
          f2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          g2 = rt_hypotd_snf(f2, scale);
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
// Arguments    : const double A[100]
//                creal_T V[100]
//                creal_T D[100]
// Return Type  : void
//
static void eig(const double A[100], creal_T V[100], creal_T D[100])
{
  boolean_T p;
  int k;
  int info;
  boolean_T exitg2;
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
  for (k = 0; k < 100; k++) {
    if ((!p) || (rtIsInf(A[k]) || rtIsNaN(A[k]))) {
      p = false;
    }
  }

  if (!p) {
    for (info = 0; info < 100; info++) {
      V[info].re = rtNaN;
      V[info].im = 0.0;
      D[info].re = 0.0;
      D[info].im = 0.0;
    }

    for (k = 0; k < 10; k++) {
      info = k + 10 * k;
      D[info].re = rtNaN;
      D[info].im = 0.0;
    }
  } else {
    p = true;
    k = 0;
    exitg2 = false;
    while ((!exitg2) && (k < 10)) {
      info = 0;
      do {
        exitg1 = 0;
        if (info <= k) {
          if (!(A[info + 10 * k] == A[k + 10 * info])) {
            p = false;
            exitg1 = 1;
          } else {
            info++;
          }
        } else {
          k++;
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

      for (k = 0; k < 9; k++) {
        b_D[(k + 10 * k) + 1] = 0.0;
        for (info = 0; info <= k; info++) {
          b_D[info + 10 * (k + 1)] = 0.0;
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
        for (k = coltop + 1; k <= info; k++) {
          absxk = std::abs(V[k - 1].re);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }

          absxk = std::abs(V[k - 1].im);
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
        for (k = coltop + 1; k <= info; k++) {
          absxk = V[k - 1].re;
          scale = V[k - 1].im;
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

          V[k - 1].re = absxk;
          V[k - 1].im = scale;
        }
      }

      std::memset(&D[0], 0, 100U * sizeof(creal_T));
      for (k = 0; k < 10; k++) {
        if (beta1[k].im == 0.0) {
          if (alpha1[k].im == 0.0) {
            info = k + 10 * k;
            D[info].re = alpha1[k].re / beta1[k].re;
            D[info].im = 0.0;
          } else if (alpha1[k].re == 0.0) {
            info = k + 10 * k;
            D[info].re = 0.0;
            D[info].im = alpha1[k].im / beta1[k].re;
          } else {
            info = k + 10 * k;
            D[info].re = alpha1[k].re / beta1[k].re;
            D[info].im = alpha1[k].im / beta1[k].re;
          }
        } else if (beta1[k].re == 0.0) {
          if (alpha1[k].re == 0.0) {
            info = k + 10 * k;
            D[info].re = alpha1[k].im / beta1[k].im;
            D[info].im = 0.0;
          } else if (alpha1[k].im == 0.0) {
            info = k + 10 * k;
            D[info].re = 0.0;
            D[info].im = -(alpha1[k].re / beta1[k].im);
          } else {
            info = k + 10 * k;
            D[info].re = alpha1[k].im / beta1[k].im;
            D[info].im = -(alpha1[k].re / beta1[k].im);
          }
        } else {
          t = std::abs(beta1[k].re);
          scale = std::abs(beta1[k].im);
          if (t > scale) {
            scale = beta1[k].im / beta1[k].re;
            absxk = beta1[k].re + scale * beta1[k].im;
            info = k + 10 * k;
            D[info].re = (alpha1[k].re + scale * alpha1[k].im) / absxk;
            D[info].im = (alpha1[k].im - scale * alpha1[k].re) / absxk;
          } else if (scale == t) {
            if (beta1[k].re > 0.0) {
              scale = 0.5;
            } else {
              scale = -0.5;
            }

            if (beta1[k].im > 0.0) {
              absxk = 0.5;
            } else {
              absxk = -0.5;
            }

            info = k + 10 * k;
            D[info].re = (alpha1[k].re * scale + alpha1[k].im * absxk) / t;
            D[info].im = (alpha1[k].im * scale - alpha1[k].re * absxk) / t;
          } else {
            scale = beta1[k].re / beta1[k].im;
            absxk = beta1[k].im + scale * beta1[k].re;
            info = k + 10 * k;
            D[info].re = (scale * alpha1[k].re + alpha1[k].im) / absxk;
            D[info].im = (scale * alpha1[k].im - alpha1[k].re) / absxk;
          }
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
            if ((1.0020841800044864E-291 > tst) || rtIsNaN(tst)) {
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
              ab = rt_hypotd_snf(v[0], tst);
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

                ab = rt_hypotd_snf(aa, b_xnrm2(nr - 1, v));
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
  std::memset(&F_eval[0], 0, 12U * sizeof(double));
  mons[0] = x * x;
  mons[1] = x * y;
  mons[2] = y * y;
  mons[3] = x;
  mons[4] = y;
  mons[5] = 1.0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      b_i = i + (j << 2);
      fc_tmp = F_eval[b_i];
      for (k = 0; k < 6; k++) {
        fc_tmp += F[b_i + 12 * k] * mons[k];
      }

      F_eval[b_i] = fc_tmp;
      b_F_eval[j + 3 * i] = fc_tmp;
    }
  }

  qr(b_F_eval, Q, F_eval);
  if (std::abs(F_eval[4]) < e) {
    // rankF = 1
    *n = 2.0;
    fc_size[0] = 1;
    fc_size[1] = 2;
    fc_tmp = Q[3] / Q[5];
    fc_data[0] = fc_tmp;
    fc_data[1] = fc_tmp;
    fs_size[0] = 1;
    fs_size[1] = 2;
    fs_data[0] = Q[4] / Q[5];
    fs_data[1] = Q[7] / Q[8];
  } else {
    *n = 1.0;
    fc_size[0] = 1;
    fc_size[1] = 1;
    fc_data[0] = Q[6] / Q[8];
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
static void mldivide(const double A[400], double B[200])
{
  double b_A[400];
  int i;
  int j;
  signed char ipiv[20];
  int mmj_tmp;
  int b;
  int jA;
  int jj;
  int k;
  int jp1j;
  int iy;
  int ix;
  double smax;
  int i1;
  double s;
  std::memcpy(&b_A[0], &A[0], 400U * sizeof(double));
  for (i = 0; i < 20; i++) {
    ipiv[i] = static_cast<signed char>((i + 1));
  }

  for (j = 0; j < 19; j++) {
    mmj_tmp = 18 - j;
    b = j * 21;
    jj = j * 21;
    jp1j = b + 2;
    iy = 20 - j;
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
        ix = j;
        for (k = 0; k < 20; k++) {
          smax = b_A[ix];
          b_A[ix] = b_A[iy];
          b_A[iy] = smax;
          ix += 20;
          iy += 20;
        }
      }

      i = (jj - j) + 20;
      for (ix = jp1j; ix <= i; ix++) {
        b_A[ix - 1] /= b_A[jj];
      }
    }

    iy = b + 20;
    jA = jj;
    for (b = 0; b <= mmj_tmp; b++) {
      smax = b_A[iy];
      if (b_A[iy] != 0.0) {
        ix = jj + 1;
        i = jA + 22;
        i1 = (jA - j) + 40;
        for (jp1j = i; jp1j <= i1; jp1j++) {
          b_A[jp1j - 1] += b_A[ix] * -smax;
          ix++;
        }
      }

      iy += 20;
      jA += 20;
    }

    if (ipiv[j] != j + 1) {
      for (iy = 0; iy < 10; iy++) {
        jA = j + 20 * iy;
        smax = B[jA];
        i = (ipiv[j] + 20 * iy) - 1;
        B[jA] = B[i];
        B[i] = smax;
      }
    }
  }

  for (j = 0; j < 10; j++) {
    jA = 20 * j;
    for (k = 0; k < 20; k++) {
      iy = 20 * k;
      i = k + jA;
      if (B[i] != 0.0) {
        i1 = k + 2;
        for (ix = i1; ix < 21; ix++) {
          b = (ix + jA) - 1;
          B[b] -= B[i] * b_A[(ix + iy) - 1];
        }
      }
    }
  }

  for (j = 0; j < 10; j++) {
    jA = 20 * j;
    for (k = 19; k >= 0; k--) {
      iy = 20 * k;
      i = k + jA;
      if (B[i] != 0.0) {
        B[i] /= b_A[k + iy];
        for (ix = 0; ix < k; ix++) {
          i1 = ix + jA;
          B[i1] -= B[i] * b_A[ix + iy];
        }
      }
    }
  }
}

//
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
  std::memset(&G20[0], 0, 900U * sizeof(double));

  // multiply by qx^2
  // multiply by qxqy
  // multiply by qx
  // multiply by qy
  // multiply by 1
  for (i = 0; i < 4; i++) {
    G20[i] = G4[i];
    G20_tmp = G4[i + 4];
    G20[i + 20] = G20_tmp;
    b_G20_tmp = G4[i + 8];
    G20[i + 40] = b_G20_tmp;
    c_G20_tmp = G4[i + 12];
    G20[i + 60] = c_G20_tmp;
    d_G20_tmp = G4[i + 16];
    G20[i + 80] = d_G20_tmp;
    e_G20_tmp = G4[i + 20];
    G20[i + 100] = e_G20_tmp;
    f_G20_tmp = G4[i + 24];
    G20[i + 120] = f_G20_tmp;
    g_G20_tmp = G4[i + 28];
    G20[i + 180] = g_G20_tmp;
    h_G20_tmp = G4[i + 32];
    G20[i + 200] = h_G20_tmp;
    i_G20_tmp = G4[i + 36];
    G20[i + 220] = i_G20_tmp;
    j_G20_tmp = G4[i + 40];
    G20[i + 240] = j_G20_tmp;
    k_G20_tmp = G4[i + 44];
    G20[i + 260] = k_G20_tmp;
    l_G20_tmp = G4[i + 48];
    G20[i + 280] = l_G20_tmp;
    m_G20_tmp = G4[i + 52];
    G20[i + 340] = m_G20_tmp;
    n_G20_tmp = G4[i + 56];
    G20[i + 360] = n_G20_tmp;
    o_G20_tmp = G4[i + 60];
    G20[i + 380] = o_G20_tmp;
    p_G20_tmp = G4[i + 64];
    G20[i + 400] = p_G20_tmp;
    q_G20_tmp = G4[i + 68];
    G20[i + 420] = q_G20_tmp;
    r_G20_tmp = G4[i + 72];
    G20[i + 480] = r_G20_tmp;
    s_G20_tmp = G4[i + 76];
    G20[i + 500] = s_G20_tmp;
    t_G20_tmp = G4[i + 80];
    G20[i + 520] = t_G20_tmp;
    u_G20_tmp = G4[i + 84];
    G20[i + 540] = u_G20_tmp;
    v_G20_tmp = G4[i + 88];
    G20[i + 600] = v_G20_tmp;
    w_G20_tmp = G4[i + 92];
    G20[i + 620] = w_G20_tmp;
    x_G20_tmp = G4[i + 96];
    G20[i + 640] = x_G20_tmp;
    y_G20_tmp = G4[i + 100];
    G20[i + 700] = y_G20_tmp;
    ab_G20_tmp = G4[i + 104];
    G20[i + 720] = ab_G20_tmp;
    bb_G20_tmp = G4[i + 108];
    G20[i + 780] = bb_G20_tmp;
    G20[i + 24] = G4[i];
    G20[i + 44] = G20_tmp;
    G20[i + 64] = b_G20_tmp;
    G20[i + 84] = c_G20_tmp;
    G20[i + 104] = d_G20_tmp;
    G20[i + 124] = e_G20_tmp;
    G20[i + 144] = f_G20_tmp;
    G20[i + 204] = g_G20_tmp;
    G20[i + 224] = h_G20_tmp;
    G20[i + 244] = i_G20_tmp;
    G20[i + 264] = j_G20_tmp;
    G20[i + 284] = k_G20_tmp;
    G20[i + 304] = l_G20_tmp;
    G20[i + 364] = m_G20_tmp;
    G20[i + 384] = n_G20_tmp;
    G20[i + 404] = o_G20_tmp;
    G20[i + 424] = p_G20_tmp;
    G20[i + 444] = q_G20_tmp;
    G20[i + 504] = r_G20_tmp;
    G20[i + 524] = s_G20_tmp;
    G20[i + 544] = t_G20_tmp;
    G20[i + 564] = u_G20_tmp;
    G20[i + 624] = v_G20_tmp;
    G20[i + 644] = w_G20_tmp;
    G20[i + 664] = x_G20_tmp;
    G20[i + 724] = y_G20_tmp;
    G20[i + 744] = ab_G20_tmp;
    G20[i + 804] = bb_G20_tmp;
    G20[i + 188] = G4[i];
    G20[i + 208] = G20_tmp;
    G20[i + 228] = b_G20_tmp;
    G20[i + 248] = c_G20_tmp;
    G20[i + 268] = d_G20_tmp;
    G20[i + 288] = e_G20_tmp;
    G20[i + 308] = f_G20_tmp;
    G20[i + 348] = g_G20_tmp;
    G20[i + 368] = h_G20_tmp;
    G20[i + 388] = i_G20_tmp;
    G20[i + 408] = j_G20_tmp;
    G20[i + 428] = k_G20_tmp;
    G20[i + 448] = l_G20_tmp;
    G20[i + 488] = m_G20_tmp;
    G20[i + 508] = n_G20_tmp;
    G20[i + 528] = o_G20_tmp;
    G20[i + 548] = p_G20_tmp;
    G20[i + 568] = q_G20_tmp;
    G20[i + 608] = r_G20_tmp;
    G20[i + 628] = s_G20_tmp;
    G20[i + 648] = t_G20_tmp;
    G20[i + 668] = u_G20_tmp;
    G20[i + 708] = v_G20_tmp;
    G20[i + 728] = w_G20_tmp;
    G20[i + 748] = x_G20_tmp;
    G20[i + 788] = y_G20_tmp;
    G20[i + 808] = ab_G20_tmp;
    G20[i + 848] = bb_G20_tmp;
    G20[i + 212] = G4[i];
    G20[i + 232] = G20_tmp;
    G20[i + 252] = b_G20_tmp;
    G20[i + 272] = c_G20_tmp;
    G20[i + 292] = d_G20_tmp;
    G20[i + 312] = e_G20_tmp;
    G20[i + 332] = f_G20_tmp;
    G20[i + 372] = g_G20_tmp;
    G20[i + 392] = h_G20_tmp;
    G20[i + 412] = i_G20_tmp;
    G20[i + 432] = j_G20_tmp;
    G20[i + 452] = k_G20_tmp;
    G20[i + 472] = l_G20_tmp;
    G20[i + 512] = m_G20_tmp;
    G20[i + 532] = n_G20_tmp;
    G20[i + 552] = o_G20_tmp;
    G20[i + 572] = p_G20_tmp;
    G20[i + 592] = q_G20_tmp;
    G20[i + 632] = r_G20_tmp;
    G20[i + 652] = s_G20_tmp;
    G20[i + 672] = t_G20_tmp;
    G20[i + 692] = u_G20_tmp;
    G20[i + 732] = v_G20_tmp;
    G20[i + 752] = w_G20_tmp;
    G20[i + 772] = x_G20_tmp;
    G20[i + 812] = y_G20_tmp;
    G20[i + 832] = ab_G20_tmp;
    G20[i + 872] = bb_G20_tmp;
    for (b_i = 0; b_i < 28; b_i++) {
      G20[(i + 20 * (b_i + 17)) + 16] = G4[i + (b_i << 2)];
    }
  }
}

//
// Arguments    : const double a[6]
//                const double b[6]
//                double c[15]
// Return Type  : void
//
static void mult_poly22(const double a[6], const double b[6], double c[15])
{
  c[0] = a[0] * b[0];

  // x^4
  c[1] = a[0] * b[1] + a[1] * b[0];

  // x^3*y
  c[2] = (a[0] * b[2] + a[1] * b[1]) + a[2] * b[0];

  // x^2*y^2
  c[3] = a[1] * b[2] + a[2] * b[1];

  // x*y^3
  c[4] = a[2] * b[2];

  // y^4
  c[5] = a[0] * b[3] + a[3] * b[0];

  // x^3
  c[6] = ((a[0] * b[4] + a[1] * b[3]) + a[3] * b[1]) + a[4] * b[0];

  // x^2*y
  c[7] = ((a[1] * b[4] + a[2] * b[3]) + a[3] * b[2]) + a[4] * b[1];

  // x*y^2
  c[8] = a[2] * b[4] + a[4] * b[2];

  // y^3
  c[9] = (a[0] * b[5] + a[5] * b[0]) + a[3] * b[3];

  // x^2
  c[10] = ((a[1] * b[5] + a[5] * b[1]) + a[3] * b[4]) + a[4] * b[3];

  // x*y
  c[11] = (a[2] * b[5] + a[5] * b[2]) + a[4] * b[4];

  // y^2
  c[12] = a[3] * b[5] + a[5] * b[3];

  // x
  c[13] = a[4] * b[5] + a[5] * b[4];

  // y
  c[14] = a[5] * b[5];

  // 1
}

//
// Arguments    : const double d[15]
//                const double a[6]
//                double p[28]
// Return Type  : void
//
static void mult_poly42(const double d[15], const double a[6], double p[28])
{
  p[0] = a[0] * d[0];

  // x^6
  p[1] = a[0] * d[1] + a[1] * d[0];

  // x^5*y
  p[2] = (a[0] * d[2] + a[1] * d[1]) + a[2] * d[0];

  // x^4*y^2
  p[3] = (a[0] * d[3] + a[1] * d[2]) + a[2] * d[1];

  // x^3*y^3
  p[4] = (a[0] * d[4] + a[1] * d[3]) + a[2] * d[2];

  // x^2*y^4
  p[5] = a[1] * d[4] + a[2] * d[3];

  // x*y^5
  p[6] = a[2] * d[4];

  // y^6
  p[7] = a[3] * d[0] + a[0] * d[5];

  // x^5
  p[8] = ((a[3] * d[1] + a[4] * d[0]) + a[0] * d[6]) + a[1] * d[5];

  // x^4*y
  p[9] = (((a[3] * d[2] + a[4] * d[1]) + a[0] * d[7]) + a[1] * d[6]) + a[2] * d
    [5];

  // x^3*y^2
  p[10] = (((a[3] * d[3] + a[4] * d[2]) + a[0] * d[8]) + a[1] * d[7]) + a[2] *
    d[6];

  // x^2*y^3
  p[11] = ((a[3] * d[4] + a[4] * d[3]) + a[1] * d[8]) + a[2] * d[7];

  // x*y^4
  p[12] = a[4] * d[4] + a[2] * d[8];

  // y^5
  p[13] = (a[5] * d[0] + a[3] * d[5]) + a[0] * d[9];

  // x^4
  p[14] = (((a[5] * d[1] + a[3] * d[6]) + a[4] * d[5]) + a[0] * d[10]) + a[1] *
    d[9];

  // x^3*y
  p[15] = ((((a[5] * d[2] + a[3] * d[7]) + a[4] * d[6]) + a[0] * d[11]) + a[1] *
           d[10]) + a[2] * d[9];

  // x^2*y^2
  p[16] = (((a[5] * d[3] + a[3] * d[8]) + a[4] * d[7]) + a[1] * d[11]) + a[2] *
    d[10];

  // x*y^3
  p[17] = (a[5] * d[4] + a[4] * d[8]) + a[2] * d[11];

  // y^4
  p[18] = (a[5] * d[5] + a[0] * d[12]) + a[3] * d[9];

  // x^3
  p[19] = (((a[5] * d[6] + a[0] * d[13]) + a[1] * d[12]) + a[3] * d[10]) + a[4] *
    d[9];

  // x^2*y
  p[20] = (((a[5] * d[7] + a[1] * d[13]) + a[2] * d[12]) + a[3] * d[11]) + a[4] *
    d[10];

  // x*y^2
  p[21] = (a[5] * d[8] + a[2] * d[13]) + a[4] * d[11];

  // y^3
  p[22] = (a[0] * d[14] + a[5] * d[9]) + a[3] * d[12];

  // x^2
  p[23] = ((a[1] * d[14] + a[5] * d[10]) + a[3] * d[13]) + a[4] * d[12];

  // x*y
  p[24] = (a[2] * d[14] + a[5] * d[11]) + a[4] * d[13];

  // y^2
  p[25] = a[3] * d[14] + a[5] * d[12];

  // x
  p[26] = a[4] * d[14] + a[5] * d[13];

  // y
  p[27] = a[5] * d[14];

  // 1
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
        beta1 = rt_hypotd_snf(b_A[ii], c);
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

          beta1 = rt_hypotd_snf(atmp, c_xnrm2(2 - i, b_A, ii + 2));
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

  // fc:
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
  for (b_i = 0; b_i < 6; b_i++) {
    b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * F_row_tmp[b_i];
    }
  }

  a_tmp = x[F_row_tmp_tmp_tmp] - x[b_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] *
        b_F_row_tmp[b_i];
    }
  }

  // fs:
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i] = a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 1] * F_row_tmp[b_i];
    }
  }

  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b_b[F_row_tmp_tmp] += R[3 * b_i + 9 * F_row_tmp_tmp] * b_F_row_tmp[b_i];
    }
  }

  // 1:
  a = -a_tmp_tmp * x[b_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    F_row[3 * b_i + 1] = -a_tmp_tmp * b[b_i] - a_tmp * b_b[b_i];
    b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
    for (F_row_tmp_tmp = 0; F_row_tmp_tmp < 6; F_row_tmp_tmp++) {
      b[F_row_tmp_tmp] += R[(3 * b_i + 9 * F_row_tmp_tmp) + 2] * F_row_tmp[b_i];
    }
  }

  a_tmp_tmp = a_tmp * y[c_F_row_tmp_tmp_tmp];
  for (b_i = 0; b_i < 6; b_i++) {
    b_b[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 3; b_i++) {
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
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  a = std::abs(u0);
  y = std::abs(u1);
  if (a < y) {
    a /= y;
    y *= std::sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = a * std::sqrt(y * y + 1.0);
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
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
  boolean_T p;
  int k;
  int i;
  double work[10];
  int b_i;
  int knt;
  int im1n_tmp;
  int in;
  int alpha1_tmp;
  int ix0;
  double alpha1;
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
  p = true;
  for (k = 0; k < 100; k++) {
    if ((!p) || (rtIsInf(A[k]) || rtIsNaN(A[k]))) {
      p = false;
    }
  }

  if (!p) {
    for (i = 0; i < 100; i++) {
      V[i] = rtNaN;
    }

    knt = 2;
    for (k = 0; k < 9; k++) {
      if (knt <= 10) {
        std::memset(&V[(k * 10 + knt) + -1], 0, (11 - knt) * sizeof(double));
      }

      knt++;
    }

    for (i = 0; i < 100; i++) {
      T[i] = rtNaN;
    }
  } else {
    std::memcpy(&T[0], &A[0], 100U * sizeof(double));
    std::memset(&work[0], 0, 10U * sizeof(double));
    for (b_i = 0; b_i < 9; b_i++) {
      im1n_tmp = b_i * 10 + 2;
      in = (b_i + 1) * 10;
      alpha1_tmp = (b_i + 10 * b_i) + 1;
      alpha1 = T[alpha1_tmp];
      if (b_i + 3 < 10) {
        c_i = b_i + 1;
      } else {
        c_i = 8;
      }

      ix0 = c_i + im1n_tmp;
      tau[b_i] = 0.0;
      xnorm = xnrm2(8 - b_i, T, ix0);
      if (xnorm != 0.0) {
        beta1 = rt_hypotd_snf(T[alpha1_tmp], xnorm);
        if (T[alpha1_tmp] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i = (ix0 - b_i) + 7;
          do {
            knt++;
            for (k = ix0; k <= i; k++) {
              T[k - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1 *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = rt_hypotd_snf(alpha1, xnrm2(8 - b_i, T, ix0));
          if (alpha1 >= 0.0) {
            beta1 = -beta1;
          }

          tau[b_i] = (beta1 - alpha1) / beta1;
          xnorm = 1.0 / (alpha1 - beta1);
          i = (ix0 - b_i) + 7;
          for (k = ix0; k <= i; k++) {
            T[k - 1] *= xnorm;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1 = beta1;
        } else {
          tau[b_i] = (beta1 - T[alpha1_tmp]) / beta1;
          xnorm = 1.0 / (T[alpha1_tmp] - beta1);
          i = (ix0 - b_i) + 7;
          for (k = ix0; k <= i; k++) {
            T[k - 1] *= xnorm;
          }

          alpha1 = beta1;
        }
      }

      T[alpha1_tmp] = 1.0;
      jy = (b_i + im1n_tmp) - 1;
      k = in + 1;
      if (tau[b_i] != 0.0) {
        lastv = 8 - b_i;
        c_i = (jy - b_i) + 8;
        while ((lastv + 1 > 0) && (T[c_i] == 0.0)) {
          lastv--;
          c_i--;
        }

        lastc = 10;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = in + lastc;
          ix0 = knt;
          do {
            exitg1 = 0;
            if (ix0 <= knt + lastv * 10) {
              if (T[ix0 - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ix0 += 10;
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
          i = (in + 10 * lastv) + 1;
          for (knt = k; knt <= i; knt += 10) {
            c_i = 0;
            i1 = (knt + lastc) - 1;
            for (ix0 = knt; ix0 <= i1; ix0++) {
              work[c_i] += T[ix0 - 1] * T[ix];
              c_i++;
            }

            ix++;
          }
        }

        if (!(-tau[b_i] == 0.0)) {
          knt = in;
          for (k = 0; k <= lastv; k++) {
            if (T[jy] != 0.0) {
              xnorm = T[jy] * -tau[b_i];
              ix = 0;
              i = knt + 1;
              i1 = lastc + knt;
              for (c_i = i; c_i <= i1; c_i++) {
                T[c_i - 1] += work[ix] * xnorm;
                ix++;
              }
            }

            jy++;
            knt += 10;
          }
        }
      }

      xzlarf(9 - b_i, 9 - b_i, b_i + im1n_tmp, tau[b_i], T, (b_i + in) + 2, work);
      T[alpha1_tmp] = alpha1;
    }

    std::memcpy(&V[0], &T[0], 100U * sizeof(double));
    for (k = 8; k >= 0; k--) {
      ix0 = (k + 1) * 10;
      for (b_i = 0; b_i <= k; b_i++) {
        V[ix0 + b_i] = 0.0;
      }

      i = k + 3;
      for (b_i = i; b_i < 11; b_i++) {
        knt = ix0 + b_i;
        V[knt - 1] = V[knt - 11];
      }
    }

    std::memset(&V[0], 0, 10U * sizeof(double));
    V[0] = 1.0;
    knt = 8;
    std::memset(&work[0], 0, 10U * sizeof(double));
    for (b_i = 8; b_i >= 0; b_i--) {
      c_i = (b_i + b_i * 10) + 11;
      if (b_i + 1 < 9) {
        V[c_i] = 1.0;
        xzlarf(9 - b_i, 8 - b_i, c_i + 1, tau[knt], V, c_i + 11, work);
        ix0 = c_i + 2;
        i = (c_i - b_i) + 9;
        for (k = ix0; k <= i; k++) {
          V[k - 1] *= -tau[knt];
        }
      }

      V[c_i] = 1.0 - tau[knt];
      for (k = 0; k < b_i; k++) {
        V[(c_i - k) - 1] = 0.0;
      }

      knt--;
    }

    eml_dlahqr(T, V);
    knt = 4;
    for (k = 0; k < 7; k++) {
      if (knt <= 10) {
        std::memset(&T[(k * 10 + knt) + -1], 0, (11 - knt) * sizeof(double));
      }

      knt++;
    }
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
  double z;
  double scale;
  double bcmis;
  double bcmax;
  int b_b;
  int b_c;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    z = *d;
    *d = *a;
    *a = z;
    *b = -*c;
    *c = 0.0;
  } else {
    tau = *a - *d;
    if ((tau == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      p = 0.5 * tau;
      scale = std::abs(*b);
      bcmis = std::abs(*c);
      if ((scale > bcmis) || rtIsNaN(bcmis)) {
        bcmax = scale;
      } else {
        bcmax = bcmis;
      }

      if ((scale < bcmis) || rtIsNaN(bcmis)) {
        bcmis = scale;
      }

      if (!(*b < 0.0)) {
        b_b = 1;
      } else {
        b_b = -1;
      }

      if (!(*c < 0.0)) {
        b_c = 1;
      } else {
        b_c = -1;
      }

      bcmis = bcmis * static_cast<double>(b_b) * static_cast<double>(b_c);
      scale = std::abs(p);
      if ((!(scale > bcmax)) && (!rtIsNaN(bcmax))) {
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
        tau = rt_hypotd_snf(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0;
      } else {
        scale = *b + *c;
        tau = rt_hypotd_snf(scale, tau);
        *cs = std::sqrt(0.5 * (std::abs(scale) / tau + 1.0));
        if (!(scale < 0.0)) {
          b_b = 1;
        } else {
          b_b = -1;
        }

        *sn = -(p / (tau * *cs)) * static_cast<double>(b_b);
        bcmax = *a * *cs + *b * *sn;
        bcmis = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        scale = -*c * *sn + *d * *cs;
        *b = bcmis * *cs + scale * *sn;
        *c = -bcmax * *sn + z * *cs;
        z = 0.5 * ((bcmax * *cs + z * *sn) + (-bcmis * *sn + scale * *cs));
        *a = z;
        *d = z;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              scale = std::sqrt(std::abs(*b));
              bcmis = std::sqrt(std::abs(*c));
              *a = scale * bcmis;
              if (!(*c < 0.0)) {
                p = *a;
              } else {
                p = -*a;
              }

              tau = 1.0 / std::sqrt(std::abs(*b + *c));
              *a = z + p;
              *d = z - p;
              *b -= *c;
              *c = 0.0;
              bcmax = scale * tau;
              scale = bcmis * tau;
              z = *cs * bcmax - *sn * scale;
              *sn = *cs * scale + *sn * bcmax;
              *cs = z;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            z = *cs;
            *cs = -*sn;
            *sn = z;
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
  if (!(alpha1 == 0.0)) {
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
  int i;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
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
  *info = 0;
  anrm = 0.0;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i < 100)) {
    absxk = rt_hypotd_snf(A[i].re, A[i].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      i++;
    }
  }

  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    for (i = 0; i < 10; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }

    for (i = 0; i < 100; i++) {
      V[i].re = rtNaN;
      V[i].im = 0.0;
    }
  } else {
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

        for (i = 0; i < 100; i++) {
          A[i].re *= a;
          A[i].im *= a;
        }
      }
    }

    xzggbal(A, &i, &ihi, rscale);
    xzgghrd(i, ihi, A, V);
    xzhgeqz(A, i, ihi, V, info, alpha1, beta1);
    if (*info == 0) {
      xztgevc(A, V);
      xzggbak(V, i, ihi, rscale);
      for (ihi = 0; ihi < 10; ihi++) {
        absxk = std::abs(V[10 * ihi].re) + std::abs(V[10 * ihi].im);
        for (jr = 0; jr < 9; jr++) {
          i = (jr + 10 * ihi) + 1;
          ctoc = std::abs(V[i].re) + std::abs(V[i].im);
          if (ctoc > absxk) {
            absxk = ctoc;
          }
        }

        if (absxk >= 6.7178761075670888E-139) {
          absxk = 1.0 / absxk;
          for (jr = 0; jr < 10; jr++) {
            i = jr + 10 * ihi;
            V[i].re *= absxk;
            V[i].im *= absxk;
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

          for (i = 0; i < 10; i++) {
            alpha1[i].re *= a;
            alpha1[i].im *= a;
          }
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
// Arguments    : creal_T A[100]
//                int ilo
//                int ihi
//                creal_T Z[100]
//                int *info
//                creal_T alpha1[10]
//                creal_T beta1[10]
// Return Type  : void
//
static void xzhgeqz(creal_T A[100], int ilo, int ihi, creal_T Z[100], int *info,
                    creal_T alpha1[10], creal_T beta1[10])
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
              for (i = 0; i < 10; i++) {
                alpha1[i].re = rtNaN;
                alpha1[i].im = 0.0;
                beta1[i].re = rtNaN;
                beta1[i].im = 0.0;
              }

              for (ctemp_tmp = 0; ctemp_tmp < 100; ctemp_tmp++) {
                Z[ctemp_tmp].re = rtNaN;
                Z[ctemp_tmp].im = 0.0;
              }

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
      for (jp1 = 0; jp1 <= ilast; jp1++) {
        alpha1[jp1].re = rtNaN;
        alpha1[jp1].im = 0.0;
        beta1[jp1].re = rtNaN;
        beta1[jp1].im = 0.0;
      }

      for (ctemp_tmp = 0; ctemp_tmp < 100; ctemp_tmp++) {
        Z[ctemp_tmp].re = rtNaN;
        Z[ctemp_tmp].im = 0.0;
      }
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

    if (!(-tau == 0.0)) {
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
        r->re = rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        f2 = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / f2;
        sn->im = -gs_im / f2;
      } else {
        g2 = std::sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2;
        if (scale_tmp > 1.0) {
          f2 = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / f2;
          fs_im = f.im / f2;
        } else {
          scale = 7.4428285367870146E+137 * f.re;
          f2s = 7.4428285367870146E+137 * f.im;
          f2 = rt_hypotd_snf(scale, f2s);
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
  if (isInitialized_p35p_solver_double == false) {
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
// Arguments    : void
// Return Type  : void
//
void p35p_solver_initialize()
{
  rt_InitInfAndNaN();
  isInitialized_p35p_solver_double = true;
}

//
// Arguments    : void
// Return Type  : void
//
void p35p_solver_terminate()
{
  // (no terminate code required)
  isInitialized_p35p_solver_double = false;
}

//
// File trailer for p35p_solver.cpp
//
// [EOF]
//
