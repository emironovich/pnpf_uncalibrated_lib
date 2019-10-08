//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: p35p_solver.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Oct-2019 05:30:59
//

// Include Files
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Variable Definitions
static boolean_T isInitialized_p35p_solver_single = false;

// Function Declarations
static void b_sqrt(creal32_T *x);
static float b_xnrm2(int n, const float x[3]);
static void b_xrot(float x[100], int ix0, int iy0, float c, float s);
static void b_xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
                      *sn);
static float c_xnrm2(int n, const float x[12], int ix0);
static void cat(const float varargin_1_data[], const int varargin_1_size[3],
                const float varargin_2[9], float y_data[], int y_size[3]);
static void eig(const float A[100], creal32_T V[100], creal32_T D[100]);
static int eml_dlahqr(float h[100], float z[100]);
static void find_f(const float F[72], float x, float y, float e, float fc_data[],
                   int fc_size[2], float fs_data[], int fs_size[2], float *n);
static void mldivide(const float A[400], float B[200]);
static void qr(const float A[12], float Q[9], float R[12]);
static float rt_hypotf_snf(float u0, float u1);
static void schur(float A[100], float V[100]);
static void xdlanv2(float *a, float *b, float *c, float *d, float *rt1r, float
                    *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);
static void xgerc(int m, int n, float alpha1, int ix0, const float y[4], float
                  A[12], int ia0);
static float xnrm2(int n, const float x[100], int ix0);
static void xrot(int n, float x[100], int ix0, int iy0, float c, float s);
static void xzggev(creal32_T A[100], int *info, creal32_T alpha1[10], creal32_T
                   beta1[10], creal32_T V[100]);
static void xzhgeqz(creal32_T A[100], int ilo, int ihi, creal32_T Z[100], int
                    *info, creal32_T alpha1[10], creal32_T beta1[10]);
static void xzlarf(int m, int n, int iv0, float tau, float C[100], int ic0,
                   float work[10]);
static void xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
                    *sn, creal32_T *r);
static void xztgevc(const creal32_T A[100], creal32_T V[100]);

// Function Definitions

//
// Arguments    : creal32_T *x
// Return Type  : void
//
static void b_sqrt(creal32_T *x)
{
  float xr;
  float xi;
  float absxi;
  float absxr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0F) {
    if (xr < 0.0F) {
      absxi = 0.0F;
      xr = std::sqrt(-xr);
    } else {
      absxi = std::sqrt(xr);
      xr = 0.0F;
    }
  } else if (xr == 0.0F) {
    if (xi < 0.0F) {
      absxi = std::sqrt(-xi / 2.0F);
      xr = -absxi;
    } else {
      absxi = std::sqrt(xi / 2.0F);
      xr = absxi;
    }
  } else if (rtIsNaNF(xr)) {
    absxi = xr;
  } else if (rtIsNaNF(xi)) {
    absxi = xi;
    xr = xi;
  } else if (rtIsInfF(xi)) {
    absxi = std::abs(xi);
    xr = xi;
  } else if (rtIsInfF(xr)) {
    if (xr < 0.0F) {
      absxi = 0.0F;
      xr = xi * -xr;
    } else {
      absxi = xr;
      xr = 0.0F;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 8.50705867E+37F) || (absxi > 8.50705867E+37F)) {
      absxr *= 0.5F;
      absxi = rt_hypotf_snf(absxr, absxi * 0.5F);
      if (absxi > absxr) {
        absxi = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0F);
      } else {
        absxi = std::sqrt(absxi) * 1.41421354F;
      }
    } else {
      absxi = std::sqrt((rt_hypotf_snf(absxr, absxi) + absxr) * 0.5F);
    }

    if (xr > 0.0F) {
      xr = 0.5F * (xi / absxi);
    } else {
      if (xi < 0.0F) {
        xr = -absxi;
      } else {
        xr = absxi;
      }

      absxi = 0.5F * (xi / xr);
    }
  }

  x->re = absxi;
  x->im = xr;
}

//
// Arguments    : int n
//                const float x[3]
// Return Type  : float
//
static float b_xnrm2(int n, const float x[3])
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
// Arguments    : float x[100]
//                int ix0
//                int iy0
//                float c
//                float s
// Return Type  : void
//
static void b_xrot(float x[100], int ix0, int iy0, float c, float s)
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
// Arguments    : const creal32_T f
//                const creal32_T g
//                float *cs
//                creal32_T *sn
// Return Type  : void
//
static void b_xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
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
        g2 = rt_hypotf_snf(gs_re, gs_im);
        sn->re = gs_re / g2;
        sn->im = -gs_im / g2;
      } else {
        g2s = std::sqrt(g2);
        *cs = rt_hypotf_snf(fs_re, fs_im) / g2s;
        if (scale_tmp > 1.0F) {
          g2 = rt_hypotf_snf(f.re, f.im);
          fs_re = f.re / g2;
          fs_im = f.im / g2;
        } else {
          f2 = 5.49755814E+11F * f.re;
          scale = 5.49755814E+11F * f.im;
          g2 = rt_hypotf_snf(f2, scale);
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
// Arguments    : int n
//                const float x[12]
//                int ix0
// Return Type  : float
//
static float c_xnrm2(int n, const float x[12], int ix0)
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
// Arguments    : const float varargin_1_data[]
//                const int varargin_1_size[3]
//                const float varargin_2[9]
//                float y_data[]
//                int y_size[3]
// Return Type  : void
//
static void cat(const float varargin_1_data[], const int varargin_1_size[3],
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
static void eig(const float A[100], creal32_T V[100], creal32_T D[100])
{
  boolean_T p;
  int k;
  int info;
  boolean_T exitg2;
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
  for (k = 0; k < 100; k++) {
    if ((!p) || (rtIsInfF(A[k]) || rtIsNaNF(A[k]))) {
      p = false;
    }
  }

  if (!p) {
    for (info = 0; info < 100; info++) {
      V[info].re = rtNaNF;
      V[info].im = 0.0F;
      D[info].re = 0.0F;
      D[info].im = 0.0F;
    }

    for (k = 0; k < 10; k++) {
      info = k + 10 * k;
      D[info].re = rtNaNF;
      D[info].im = 0.0F;
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
      std::memcpy(&b_D[0], &A[0], 100U * sizeof(float));
      schur(b_D, b_V);
      for (info = 0; info < 100; info++) {
        V[info].re = b_V[info];
        V[info].im = 0.0F;
      }

      for (k = 0; k < 9; k++) {
        b_D[(k + 10 * k) + 1] = 0.0F;
        for (info = 0; info <= k; info++) {
          b_D[info + 10 * (k + 1)] = 0.0F;
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

      xzggev(At, &info, alpha1, beta1, V);
      for (coltop = 0; coltop <= 90; coltop += 10) {
        colnorm = 0.0F;
        scale = 1.29246971E-26F;
        info = coltop + 10;
        for (k = coltop + 1; k <= info; k++) {
          absxk = std::abs(V[k - 1].re);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0F;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }

          absxk = std::abs(V[k - 1].im);
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
        for (k = coltop + 1; k <= info; k++) {
          absxk = V[k - 1].re;
          scale = V[k - 1].im;
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

          V[k - 1].re = absxk;
          V[k - 1].im = scale;
        }
      }

      std::memset(&D[0], 0, 100U * sizeof(creal32_T));
      for (k = 0; k < 10; k++) {
        if (beta1[k].im == 0.0F) {
          if (alpha1[k].im == 0.0F) {
            info = k + 10 * k;
            D[info].re = alpha1[k].re / beta1[k].re;
            D[info].im = 0.0F;
          } else if (alpha1[k].re == 0.0F) {
            info = k + 10 * k;
            D[info].re = 0.0F;
            D[info].im = alpha1[k].im / beta1[k].re;
          } else {
            info = k + 10 * k;
            D[info].re = alpha1[k].re / beta1[k].re;
            D[info].im = alpha1[k].im / beta1[k].re;
          }
        } else if (beta1[k].re == 0.0F) {
          if (alpha1[k].re == 0.0F) {
            info = k + 10 * k;
            D[info].re = alpha1[k].im / beta1[k].im;
            D[info].im = 0.0F;
          } else if (alpha1[k].im == 0.0F) {
            info = k + 10 * k;
            D[info].re = 0.0F;
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
            if (beta1[k].re > 0.0F) {
              scale = 0.5F;
            } else {
              scale = -0.5F;
            }

            if (beta1[k].im > 0.0F) {
              absxk = 0.5F;
            } else {
              absxk = -0.5F;
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
// Arguments    : float h[100]
//                float z[100]
// Return Type  : int
//
static int eml_dlahqr(float h[100], float z[100])
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
            if ((9.86076132E-31F > tst) || rtIsNaNF(tst)) {
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
            tst = b_xnrm2(nr - 1, v);
            if (tst != 0.0F) {
              ab = rt_hypotf_snf(v[0], tst);
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

                ab = rt_hypotf_snf(aa, b_xnrm2(nr - 1, v));
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
        xdlanv2(&h[(b_i + 10 * (b_i - 1)) - 1], &s, &f, &tst, &aa, &ab, &bb, &ba,
                &h22, &rt1r);
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
        xrot(b_i - 1, h, ix, hoffset, h22, rt1r);
        b_xrot(z, ix, hoffset, h22, rt1r);
      }

      b_i = L - 2;
    }
  }

  return info;
}

//
// Arguments    : const float F[72]
//                float x
//                float y
//                float e
//                float fc_data[]
//                int fc_size[2]
//                float fs_data[]
//                int fs_size[2]
//                float *n
// Return Type  : void
//
static void find_f(const float F[72], float x, float y, float e, float fc_data[],
                   int fc_size[2], float fs_data[], int fs_size[2], float *n)
{
  int i;
  float mons[6];
  float F_eval[12];
  int b_i;
  float b_F_eval[12];
  float Q[9];
  int j;
  float fc_tmp;
  int k;
  for (i = 0; i < 12; i++) {
    F_eval[i] = 0.0F;
  }

  mons[0] = x * x;
  mons[1] = x * y;
  mons[2] = y * y;
  mons[3] = x;
  mons[4] = y;
  mons[5] = 1.0F;
  for (b_i = 0; b_i < 4; b_i++) {
    for (j = 0; j < 3; j++) {
      i = b_i + (j << 2);
      fc_tmp = F_eval[i];
      for (k = 0; k < 6; k++) {
        fc_tmp += F[i + 12 * k] * mons[k];
      }

      F_eval[i] = fc_tmp;
      b_F_eval[j + 3 * b_i] = fc_tmp;
    }
  }

  qr(b_F_eval, Q, F_eval);
  if (std::abs(F_eval[4]) < e) {
    // rankF = 1
    *n = 2.0F;
    fc_tmp = Q[3] / Q[5];
    fc_size[0] = 1;
    fc_size[1] = 2;
    fc_data[0] = fc_tmp;
    fc_data[1] = fc_tmp;
    fs_size[0] = 1;
    fs_size[1] = 2;
    fs_data[0] = Q[4] / Q[5];
    fs_data[1] = Q[7] / Q[8];
  } else {
    *n = 1.0F;
    fc_size[0] = 1;
    fc_size[1] = 1;
    fc_data[0] = Q[6] / Q[8];
    fs_size[0] = 1;
    fs_size[1] = 1;
    fs_data[0] = Q[7] / Q[8];
  }
}

//
// Arguments    : const float A[400]
//                float B[200]
// Return Type  : void
//
static void mldivide(const float A[400], float B[200])
{
  float b_A[400];
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
  float smax;
  int i1;
  float s;
  std::memcpy(&b_A[0], &A[0], 400U * sizeof(float));
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

    if (b_A[jj + jA] != 0.0F) {
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
      if (b_A[iy] != 0.0F) {
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
      if (B[i] != 0.0F) {
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
      if (B[i] != 0.0F) {
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
// Arguments    : const float A[12]
//                float Q[9]
//                float R[12]
// Return Type  : void
//
static void qr(const float A[12], float Q[9], float R[12])
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
      c = c_xnrm2(2 - b_i, b_A, ii + 2);
      if (c != 0.0F) {
        beta1 = rt_hypotf_snf(b_A[ii], c);
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

          beta1 = rt_hypotf_snf(atmp, c_xnrm2(2 - b_i, b_A, ii + 2));
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

      xgerc(lastv, lastc, -tau[b_i], ii + 1, work, b_A, ii + 4);
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

        xgerc(lastv, lastc, -tau[ii], iaii - 3, work, b_A, iaii);
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
// Arguments    : float u0
//                float u1
// Return Type  : float
//
static float rt_hypotf_snf(float u0, float u1)
{
  float y;
  float a;
  a = std::abs(u0);
  y = std::abs(u1);
  if (a < y) {
    a /= y;
    y *= std::sqrt(a * a + 1.0F);
  } else if (a > y) {
    y /= a;
    y = a * std::sqrt(y * y + 1.0F);
  } else {
    if (!rtIsNaNF(y)) {
      y = a * 1.41421354F;
    }
  }

  return y;
}

//
// Arguments    : float A[100]
//                float V[100]
// Return Type  : void
//
static void schur(float A[100], float V[100])
{
  boolean_T p;
  int k;
  int i;
  int b_i;
  float xnorm;
  float work[10];
  int knt;
  int im1n_tmp;
  int in;
  int alpha1_tmp;
  int ix0;
  float alpha1;
  int c_i;
  float tau[9];
  float beta1;
  int jy;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int ix;
  int exitg1;
  int i1;
  p = true;
  for (k = 0; k < 100; k++) {
    if (p) {
      xnorm = A[k];
      if (rtIsInfF(xnorm) || rtIsNaNF(xnorm)) {
        p = false;
      }
    } else {
      p = false;
    }
  }

  if (!p) {
    for (b_i = 0; b_i < 100; b_i++) {
      V[b_i] = rtNaNF;
    }

    knt = 2;
    for (k = 0; k < 9; k++) {
      if (knt <= 10) {
        std::memset(&V[(k * 10 + knt) + -1], 0, (11 - knt) * sizeof(float));
      }

      knt++;
    }

    for (b_i = 0; b_i < 100; b_i++) {
      A[b_i] = rtNaNF;
    }
  } else {
    for (i = 0; i < 10; i++) {
      work[i] = 0.0F;
    }

    for (i = 0; i < 9; i++) {
      im1n_tmp = i * 10 + 2;
      in = (i + 1) * 10;
      alpha1_tmp = (i + 10 * i) + 1;
      alpha1 = A[alpha1_tmp];
      if (i + 3 < 10) {
        c_i = i + 1;
      } else {
        c_i = 8;
      }

      ix0 = c_i + im1n_tmp;
      tau[i] = 0.0F;
      xnorm = xnrm2(8 - i, A, ix0);
      if (xnorm != 0.0F) {
        beta1 = rt_hypotf_snf(A[alpha1_tmp], xnorm);
        if (A[alpha1_tmp] >= 0.0F) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 9.86076132E-32F) {
          knt = -1;
          b_i = (ix0 - i) + 7;
          do {
            knt++;
            for (k = ix0; k <= b_i; k++) {
              A[k - 1] *= 1.01412048E+31F;
            }

            beta1 *= 1.01412048E+31F;
            alpha1 *= 1.01412048E+31F;
          } while (!(std::abs(beta1) >= 9.86076132E-32F));

          beta1 = rt_hypotf_snf(alpha1, xnrm2(8 - i, A, ix0));
          if (alpha1 >= 0.0F) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - alpha1) / beta1;
          xnorm = 1.0F / (alpha1 - beta1);
          b_i = (ix0 - i) + 7;
          for (k = ix0; k <= b_i; k++) {
            A[k - 1] *= xnorm;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 9.86076132E-32F;
          }

          alpha1 = beta1;
        } else {
          tau[i] = (beta1 - A[alpha1_tmp]) / beta1;
          xnorm = 1.0F / (A[alpha1_tmp] - beta1);
          b_i = (ix0 - i) + 7;
          for (k = ix0; k <= b_i; k++) {
            A[k - 1] *= xnorm;
          }

          alpha1 = beta1;
        }
      }

      A[alpha1_tmp] = 1.0F;
      jy = (i + im1n_tmp) - 1;
      k = in + 1;
      if (tau[i] != 0.0F) {
        lastv = 8 - i;
        c_i = (jy - i) + 8;
        while ((lastv + 1 > 0) && (A[c_i] == 0.0F)) {
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
              if (A[ix0 - 1] != 0.0F) {
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
            std::memset(&work[0], 0, lastc * sizeof(float));
          }

          ix = jy;
          b_i = (in + 10 * lastv) + 1;
          for (knt = k; knt <= b_i; knt += 10) {
            c_i = 0;
            i1 = (knt + lastc) - 1;
            for (ix0 = knt; ix0 <= i1; ix0++) {
              work[c_i] += A[ix0 - 1] * A[ix];
              c_i++;
            }

            ix++;
          }
        }

        if (!(-tau[i] == 0.0F)) {
          knt = in;
          for (k = 0; k <= lastv; k++) {
            if (A[jy] != 0.0F) {
              xnorm = A[jy] * -tau[i];
              ix = 0;
              b_i = knt + 1;
              i1 = lastc + knt;
              for (c_i = b_i; c_i <= i1; c_i++) {
                A[c_i - 1] += work[ix] * xnorm;
                ix++;
              }
            }

            jy++;
            knt += 10;
          }
        }
      }

      xzlarf(9 - i, 9 - i, i + im1n_tmp, tau[i], A, (i + in) + 2, work);
      A[alpha1_tmp] = alpha1;
    }

    std::memcpy(&V[0], &A[0], 100U * sizeof(float));
    for (k = 8; k >= 0; k--) {
      ix0 = (k + 1) * 10;
      for (i = 0; i <= k; i++) {
        V[ix0 + i] = 0.0F;
      }

      b_i = k + 3;
      for (i = b_i; i < 11; i++) {
        knt = ix0 + i;
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
      c_i = (i + i * 10) + 11;
      if (i + 1 < 9) {
        V[c_i] = 1.0F;
        xzlarf(9 - i, 8 - i, c_i + 1, tau[knt], V, c_i + 11, work);
        ix0 = c_i + 2;
        b_i = (c_i - i) + 9;
        for (k = ix0; k <= b_i; k++) {
          V[k - 1] *= -tau[knt];
        }
      }

      V[c_i] = 1.0F - tau[knt];
      for (k = 0; k < i; k++) {
        V[(c_i - k) - 1] = 0.0F;
      }

      knt--;
    }

    eml_dlahqr(A, V);
    knt = 4;
    for (k = 0; k < 7; k++) {
      if (knt <= 10) {
        std::memset(&A[(k * 10 + knt) + -1], 0, (11 - knt) * sizeof(float));
      }

      knt++;
    }
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
static void xdlanv2(float *a, float *b, float *c, float *d, float *rt1r, float
                    *rt1i, float *rt2r, float *rt2i, float *cs, float *sn)
{
  float tau;
  float p;
  float z;
  float scale;
  float bcmis;
  float bcmax;
  int b_b;
  int b_c;
  if (*c == 0.0F) {
    *cs = 1.0F;
    *sn = 0.0F;
  } else if (*b == 0.0F) {
    *cs = 0.0F;
    *sn = 1.0F;
    z = *d;
    *d = *a;
    *a = z;
    *b = -*c;
    *c = 0.0F;
  } else {
    tau = *a - *d;
    if ((tau == 0.0F) && ((*b < 0.0F) != (*c < 0.0F))) {
      *cs = 1.0F;
      *sn = 0.0F;
    } else {
      p = 0.5F * tau;
      scale = std::abs(*b);
      bcmis = std::abs(*c);
      if ((scale > bcmis) || rtIsNaNF(bcmis)) {
        bcmax = scale;
      } else {
        bcmax = bcmis;
      }

      if ((scale < bcmis) || rtIsNaNF(bcmis)) {
        bcmis = scale;
      }

      if (!(*b < 0.0F)) {
        b_b = 1;
      } else {
        b_b = -1;
      }

      if (!(*c < 0.0F)) {
        b_c = 1;
      } else {
        b_c = -1;
      }

      bcmis = bcmis * static_cast<float>(b_b) * static_cast<float>(b_c);
      scale = std::abs(p);
      if ((!(scale > bcmax)) && (!rtIsNaNF(bcmax))) {
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
        tau = rt_hypotf_snf(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0F;
      } else {
        scale = *b + *c;
        tau = rt_hypotf_snf(scale, tau);
        *cs = std::sqrt(0.5F * (std::abs(scale) / tau + 1.0F));
        if (!(scale < 0.0F)) {
          b_b = 1;
        } else {
          b_b = -1;
        }

        *sn = -(p / (tau * *cs)) * static_cast<float>(b_b);
        bcmax = *a * *cs + *b * *sn;
        bcmis = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        scale = -*c * *sn + *d * *cs;
        *b = bcmis * *cs + scale * *sn;
        *c = -bcmax * *sn + z * *cs;
        z = 0.5F * ((bcmax * *cs + z * *sn) + (-bcmis * *sn + scale * *cs));
        *a = z;
        *d = z;
        if (*c != 0.0F) {
          if (*b != 0.0F) {
            if ((*b < 0.0F) == (*c < 0.0F)) {
              scale = std::sqrt(std::abs(*b));
              bcmis = std::sqrt(std::abs(*c));
              *a = scale * bcmis;
              if (!(*c < 0.0F)) {
                p = *a;
              } else {
                p = -*a;
              }

              tau = 1.0F / std::sqrt(std::abs(*b + *c));
              *a = z + p;
              *d = z - p;
              *b -= *c;
              *c = 0.0F;
              bcmax = scale * tau;
              scale = bcmis * tau;
              z = *cs * bcmax - *sn * scale;
              *sn = *cs * scale + *sn * bcmax;
              *cs = z;
            }
          } else {
            *b = -*c;
            *c = 0.0F;
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
  if (*c == 0.0F) {
    *rt1i = 0.0F;
    *rt2i = 0.0F;
  } else {
    *rt1i = std::sqrt(std::abs(*b)) * std::sqrt(std::abs(*c));
    *rt2i = -*rt1i;
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
static void xgerc(int m, int n, float alpha1, int ix0, const float y[4], float
                  A[12], int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i;
  int ijA;
  if (!(alpha1 == 0.0F)) {
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
//                const float x[100]
//                int ix0
// Return Type  : float
//
static float xnrm2(int n, const float x[100], int ix0)
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
// Arguments    : int n
//                float x[100]
//                int ix0
//                int iy0
//                float c
//                float s
// Return Type  : void
//
static void xrot(int n, float x[100], int ix0, int iy0, float c, float s)
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
// Arguments    : creal32_T A[100]
//                int *info
//                creal32_T alpha1[10]
//                creal32_T beta1[10]
//                creal32_T V[100]
// Return Type  : void
//
static void xzggev(creal32_T A[100], int *info, creal32_T alpha1[10], creal32_T
                   beta1[10], creal32_T V[100])
{
  float anrm;
  int k;
  boolean_T exitg1;
  float absxk;
  boolean_T ilascl;
  int i;
  float anrmto;
  boolean_T guard1 = false;
  int jcolp1;
  float ctoc;
  boolean_T notdone;
  int ilo;
  int rscale[10];
  int ihi;
  float stemp_im;
  int exitg3;
  float cto1;
  int j;
  float a;
  int ii;
  int nzcount;
  int jcol;
  boolean_T exitg4;
  creal32_T atmp;
  signed char b_I[100];
  int exitg2;
  int jrow;
  float f;
  *info = 0;
  anrm = 0.0F;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 100)) {
    absxk = rt_hypotf_snf(A[k].re, A[k].im);
    if (rtIsNaNF(absxk)) {
      anrm = rtNaNF;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      k++;
    }
  }

  if (rtIsInfF(anrm) || rtIsNaNF(anrm)) {
    for (i = 0; i < 10; i++) {
      alpha1[i].re = rtNaNF;
      alpha1[i].im = 0.0F;
      beta1[i].re = rtNaNF;
      beta1[i].im = 0.0F;
    }

    for (jcolp1 = 0; jcolp1 < 100; jcolp1++) {
      V[jcolp1].re = rtNaNF;
      V[jcolp1].im = 0.0F;
    }
  } else {
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
      exitg3 = 0;
      i = 0;
      j = 0;
      notdone = false;
      ii = ihi;
      exitg1 = false;
      while ((!exitg1) && (ii > 0)) {
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
          exitg1 = true;
        } else {
          ii--;
        }
      }

      if (!notdone) {
        exitg3 = 2;
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
          exitg3 = 1;
        }
      }
    } while (exitg3 == 0);

    if (exitg3 != 1) {
      do {
        exitg2 = 0;
        i = 0;
        j = 0;
        notdone = false;
        k = ilo;
        exitg1 = false;
        while ((!exitg1) && (k <= ihi)) {
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
            exitg1 = true;
          } else {
            k++;
          }
        }

        if (!notdone) {
          exitg2 = 1;
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
            exitg2 = 1;
          }
        }
      } while (exitg2 == 0);
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
          xzlartg(A[k - 1], A[k], &absxk, &atmp, &A[(jrow + 10 * jcol) - 1]);
          A[k].re = 0.0F;
          A[k].im = 0.0F;
          for (j = jcolp1; j < 11; j++) {
            ii = jrow + 10 * (j - 1);
            nzcount = ii - 1;
            ctoc = absxk * A[nzcount].re + (atmp.re * A[ii].re - atmp.im * A[ii]
              .im);
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
            ctoc = absxk * A[nzcount].re + (atmp.re * A[ii].re - atmp.im * A[ii]
              .im);
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
            V[ii].re = absxk * V[ii].re - (cto1 * V[nzcount].re + a * V[nzcount]
              .im);
            V[ii].im = absxk * V[ii].im - (cto1 * V[nzcount].im - a * f);
            V[nzcount].re = ctoc;
            V[nzcount].im = stemp_im;
          }
        }
      }
    }

    xzhgeqz(A, ilo, ihi, V, info, alpha1, beta1);
    if (*info == 0) {
      xztgevc(A, V);
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
static void xzhgeqz(creal32_T A[100], int ilo, int ihi, creal32_T Z[100], int
                    *info, creal32_T alpha1[10], creal32_T beta1[10])
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
              for (i = 0; i < 10; i++) {
                alpha1[i].re = rtNaNF;
                alpha1[i].im = 0.0F;
                beta1[i].re = rtNaNF;
                beta1[i].im = 0.0F;
              }

              for (ctemp_tmp = 0; ctemp_tmp < 100; ctemp_tmp++) {
                Z[ctemp_tmp].re = rtNaNF;
                Z[ctemp_tmp].im = 0.0F;
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
                b_sqrt(&shift);
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
              b_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              jp1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  ctemp_tmp = j + 10 * jp1;
                  xzlartg(A[ctemp_tmp - 1], A[ctemp_tmp], &anorm, &shift, &A[(j
                           + 10 * jp1) - 1]);
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
      for (jp1 = 0; jp1 <= ilast; jp1++) {
        alpha1[jp1].re = rtNaNF;
        alpha1[jp1].im = 0.0F;
        beta1[jp1].re = rtNaNF;
        beta1[jp1].im = 0.0F;
      }

      for (ctemp_tmp = 0; ctemp_tmp < 100; ctemp_tmp++) {
        Z[ctemp_tmp].re = rtNaNF;
        Z[ctemp_tmp].im = 0.0F;
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
//                float tau
//                float C[100]
//                int ic0
//                float work[10]
// Return Type  : void
//
static void xzlarf(int m, int n, int iv0, float tau, float C[100], int ic0,
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

    if (!(-tau == 0.0F)) {
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
// Arguments    : const creal32_T f
//                const creal32_T g
//                float *cs
//                creal32_T *sn
//                creal32_T *r
// Return Type  : void
//
static void xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T
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
        r->re = rt_hypotf_snf(g.re, g.im);
        r->im = 0.0F;
        f2 = rt_hypotf_snf(gs_re, gs_im);
        sn->re = gs_re / f2;
        sn->im = -gs_im / f2;
      } else {
        g2 = std::sqrt(g2);
        *cs = rt_hypotf_snf(fs_re, fs_im) / g2;
        if (scale_tmp > 1.0F) {
          f2 = rt_hypotf_snf(f.re, f.im);
          fs_re = f.re / f2;
          fs_im = f.im / f2;
        } else {
          scale = 5.49755814E+11F * f.re;
          f2s = 5.49755814E+11F * f.im;
          f2 = rt_hypotf_snf(scale, f2s);
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
// Arguments    : const creal32_T A[100]
//                creal32_T V[100]
// Return Type  : void
//
static void xztgevc(const creal32_T A[100], creal32_T V[100])
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
// Arguments    : const float X[12]
//                const float x[4]
//                const float y[4]
//                float e
//                float *solution_num
//                float f_sol_data[]
//                int f_sol_size[2]
//                float R_sol_data[]
//                int R_sol_size[3]
//                float T_sol_data[]
//                int T_sol_size[2]
// Return Type  : void
//
void p35p_solver(const float X[12], const float x[4], const float y[4], float e,
                 float *solution_num, float f_sol_data[], int f_sol_size[2],
                 float R_sol_data[], int R_sol_size[3], float T_sol_data[], int
                 T_sol_size[2])
{
  float F[72];
  float F_row[18];
  float f;
  float T_tmp;
  float T[3];
  float p4_tmp;
  float p4[3];
  float b_T_tmp;
  float b_p4_tmp;
  float c_T_tmp;
  float c_p4_tmp;
  float F_row_tmp_tmp;
  int i;
  int j;
  float b[6];
  int b_i;
  float b_b[6];
  static const signed char R[54] = { 1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 2, 0, 2, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0,
    -2, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  float a;
  float qy_im;
  float G20_tmp;
  float G[112];
  float b_a[15];
  float c_b[15];
  float c_a[28];
  float M[24];
  int M_tmp;
  float d_b[28];
  float e_b[28];
  float b_M[54];
  int i1;
  int b_M_tmp;
  float G20[900];
  float C[200];
  float b_G20[400];
  float qy_re;
  float b_G20_tmp;
  float c_M[100];
  float c_G20_tmp;
  float d_G20_tmp;
  float e_G20_tmp;
  float f_G20_tmp;
  float g_G20_tmp;
  float d_M[100];
  creal32_T W[100];
  creal32_T D[100];
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
  float fc_set_data[2];
  int fc_set_size[2];
  float fs_set_data[2];
  int fs_set_size[2];
  float s_G20_tmp;
  float t_G20_tmp;
  int fc_set_tmp;
  float R_xy[9];
  float b_x;
  float fc_set[9];
  float b_fc_set[12];
  float R_curr[9];
  float b_X[4];
  float tmp_data[99];
  int tmp_size[3];
  float b_T_sol_data[30];
  if (isInitialized_p35p_solver_single == false) {
    p35p_solver_initialize();
  }

  std::memset(&F[0], 0, 72U * sizeof(float));
  std::memset(&F_row[0], 0, 18U * sizeof(float));

  // fc:
  f = y[0] - y[2];
  T_tmp = X[3] - X[0];
  T[0] = T_tmp;
  p4_tmp = X[6] - X[0];
  p4[0] = p4_tmp;
  b_T_tmp = X[4] - X[1];
  T[1] = b_T_tmp;
  b_p4_tmp = X[7] - X[1];
  p4[1] = b_p4_tmp;
  c_T_tmp = X[5] - X[2];
  T[2] = c_T_tmp;
  c_p4_tmp = X[8] - X[2];
  p4[2] = c_p4_tmp;
  F_row_tmp_tmp = x[0] - x[1];
  for (i = 0; i < 6; i++) {
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * p4[j];
    }
  }

  // fs:
  for (i = 0; i < 6; i++) {
    F_row[3 * i] = f * b[i] - F_row_tmp_tmp * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * p4[j];
    }
  }

  // 1:
  a = -f * x[1];
  for (i = 0; i < 6; i++) {
    F_row[3 * i + 1] = -f * b[i] - F_row_tmp_tmp * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * T[j];
    }
  }

  qy_im = F_row_tmp_tmp * y[2];
  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * p4[j];
    }
  }

  for (b_i = 0; b_i < 6; b_i++) {
    G20_tmp = a * b[b_i] + qy_im * b_b[b_i];
    b[b_i] = G20_tmp;
    F_row[3 * b_i + 2] = G20_tmp;
    F[12 * b_i] = F_row[3 * b_i];
    F[12 * b_i + 4] = F_row[3 * b_i + 1];
    F[12 * b_i + 8] = G20_tmp;
  }

  std::memset(&F_row[0], 0, 18U * sizeof(float));

  // fc:
  f = y[0] - y[1];
  T[0] = p4_tmp;
  p4[0] = T_tmp;
  T[1] = b_p4_tmp;
  p4[1] = b_T_tmp;
  T[2] = c_p4_tmp;
  p4[2] = c_T_tmp;
  qy_im = x[0] - x[2];
  for (i = 0; i < 6; i++) {
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * p4[j];
    }
  }

  // fs:
  for (i = 0; i < 6; i++) {
    F_row[3 * i] = f * b[i] - qy_im * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * p4[j];
    }
  }

  // 1:
  a = -f * x[2];
  for (i = 0; i < 6; i++) {
    F_row[3 * i + 1] = -f * b[i] - qy_im * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * T[j];
    }
  }

  qy_im *= y[1];
  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * p4[j];
    }
  }

  for (b_i = 0; b_i < 6; b_i++) {
    G20_tmp = a * b[b_i] + qy_im * b_b[b_i];
    b[b_i] = G20_tmp;
    F_row[3 * b_i + 2] = G20_tmp;
    F[12 * b_i + 1] = F_row[3 * b_i];
    F[12 * b_i + 5] = F_row[3 * b_i + 1];
    F[12 * b_i + 9] = G20_tmp;
  }

  std::memset(&F_row[0], 0, 18U * sizeof(float));

  // fc:
  f = y[1] - y[2];
  T[0] = X[9] - X[3];
  p4[0] = X[6] - X[3];
  T[1] = X[10] - X[4];
  p4[1] = X[7] - X[4];
  T[2] = X[11] - X[5];
  p4[2] = X[8] - X[5];
  qy_im = x[1] - x[3];
  for (i = 0; i < 6; i++) {
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * p4[j];
    }
  }

  // fs:
  for (i = 0; i < 6; i++) {
    F_row[3 * i] = f * b[i] - qy_im * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * p4[j];
    }
  }

  // 1:
  a = -f * x[3];
  for (i = 0; i < 6; i++) {
    F_row[3 * i + 1] = -f * b[i] - qy_im * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * T[j];
    }
  }

  qy_im *= y[2];
  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * p4[j];
    }
  }

  for (b_i = 0; b_i < 6; b_i++) {
    G20_tmp = a * b[b_i] + qy_im * b_b[b_i];
    b[b_i] = G20_tmp;
    F_row[3 * b_i + 2] = G20_tmp;
    F[12 * b_i + 2] = F_row[3 * b_i];
    F[12 * b_i + 6] = F_row[3 * b_i + 1];
    F[12 * b_i + 10] = G20_tmp;
  }

  std::memset(&F_row[0], 0, 18U * sizeof(float));

  // fc:
  f = y[2] - y[1];
  T[0] = X[9] - X[6];
  p4[0] = X[3] - X[6];
  T[1] = X[10] - X[7];
  p4[1] = X[4] - X[7];
  T[2] = X[11] - X[8];
  p4[2] = X[5] - X[8];
  qy_im = x[2] - x[3];
  for (i = 0; i < 6; i++) {
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * p4[j];
    }
  }

  // fs:
  for (i = 0; i < 6; i++) {
    F_row[3 * i] = f * b[i] - qy_im * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 1]) * T[j];
    }
  }

  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[3 * j + 9 * b_i]) * p4[j];
    }
  }

  // 1:
  a = -f * x[3];
  for (i = 0; i < 6; i++) {
    F_row[3 * i + 1] = -f * b[i] - qy_im * b_b[i];
    b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * T[j];
    }
  }

  qy_im *= y[1];
  for (i = 0; i < 6; i++) {
    b_b[i] = 0.0F;
  }

  for (j = 0; j < 3; j++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_b[b_i] += static_cast<float>(R[(3 * j + 9 * b_i) + 2]) * p4[j];
    }
  }

  for (b_i = 0; b_i < 6; b_i++) {
    G20_tmp = a * b[b_i] + qy_im * b_b[b_i];
    b[b_i] = G20_tmp;
    F_row[3 * b_i + 2] = G20_tmp;
    F[12 * b_i + 3] = F_row[3 * b_i];
    F[12 * b_i + 7] = F_row[3 * b_i + 1];
    F[12 * b_i + 11] = G20_tmp;
  }

  // 4x3x6
  std::memset(&G[0], 0, 112U * sizeof(float));
  b_a[0] = F[6] * F[11];

  // x^4
  b_a[1] = F[6] * F[23] + F[18] * F[11];

  // x^3*y
  b_a[2] = (F[6] * F[35] + F[18] * F[23]) + F[30] * F[11];

  // x^2*y^2
  b_a[3] = F[18] * F[35] + F[30] * F[23];

  // x*y^3
  b_a[4] = F[30] * F[35];

  // y^4
  b_a[5] = F[6] * F[47] + F[42] * F[11];

  // x^3
  b_a[6] = ((F[6] * F[59] + F[18] * F[47]) + F[42] * F[23]) + F[54] * F[11];

  // x^2*y
  b_a[7] = ((F[18] * F[59] + F[30] * F[47]) + F[42] * F[35]) + F[54] * F[23];

  // x*y^2
  b_a[8] = F[30] * F[59] + F[54] * F[35];

  // y^3
  b_a[9] = (F[6] * F[71] + F[66] * F[11]) + F[42] * F[47];

  // x^2
  b_a[10] = ((F[18] * F[71] + F[66] * F[23]) + F[42] * F[59]) + F[54] * F[47];

  // x*y
  b_a[11] = (F[30] * F[71] + F[66] * F[35]) + F[54] * F[59];

  // y^2
  b_a[12] = F[42] * F[71] + F[66] * F[47];

  // x
  b_a[13] = F[54] * F[71] + F[66] * F[59];

  // y
  b_a[14] = F[66] * F[71];

  // 1
  c_b[0] = F[10] * F[7];

  // x^4
  c_b[1] = F[10] * F[19] + F[22] * F[7];

  // x^3*y
  c_b[2] = (F[10] * F[31] + F[22] * F[19]) + F[34] * F[7];

  // x^2*y^2
  c_b[3] = F[22] * F[31] + F[34] * F[19];

  // x*y^3
  c_b[4] = F[34] * F[31];

  // y^4
  c_b[5] = F[10] * F[43] + F[46] * F[7];

  // x^3
  c_b[6] = ((F[10] * F[55] + F[22] * F[43]) + F[46] * F[19]) + F[58] * F[7];

  // x^2*y
  c_b[7] = ((F[22] * F[55] + F[34] * F[43]) + F[46] * F[31]) + F[58] * F[19];

  // x*y^2
  c_b[8] = F[34] * F[55] + F[58] * F[31];

  // y^3
  c_b[9] = (F[10] * F[67] + F[70] * F[7]) + F[46] * F[43];

  // x^2
  c_b[10] = ((F[22] * F[67] + F[70] * F[19]) + F[46] * F[55]) + F[58] * F[43];

  // x*y
  c_b[11] = (F[34] * F[67] + F[70] * F[31]) + F[58] * F[55];

  // y^2
  c_b[12] = F[46] * F[67] + F[70] * F[43];

  // x
  c_b[13] = F[58] * F[67] + F[70] * F[55];

  // y
  c_b[14] = F[70] * F[67];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  c_a[0] = F[1] * b_a[0];

  // x^6
  c_a[1] = F[1] * b_a[1] + F[13] * b_a[0];

  // x^5*y
  c_a[2] = (F[1] * b_a[2] + F[13] * b_a[1]) + F[25] * b_a[0];

  // x^4*y^2
  c_a[3] = (F[1] * b_a[3] + F[13] * b_a[2]) + F[25] * b_a[1];

  // x^3*y^3
  c_a[4] = (F[1] * b_a[4] + F[13] * b_a[3]) + F[25] * b_a[2];

  // x^2*y^4
  c_a[5] = F[13] * b_a[4] + F[25] * b_a[3];

  // x*y^5
  c_a[6] = F[25] * b_a[4];

  // y^6
  c_a[7] = F[37] * b_a[0] + F[1] * b_a[5];

  // x^5
  c_a[8] = ((F[37] * b_a[1] + F[49] * b_a[0]) + F[1] * b_a[6]) + F[13] * b_a[5];

  // x^4*y
  c_a[9] = (((F[37] * b_a[2] + F[49] * b_a[1]) + F[1] * b_a[7]) + F[13] * b_a[6])
    + F[25] * b_a[5];

  // x^3*y^2
  c_a[10] = (((F[37] * b_a[3] + F[49] * b_a[2]) + F[1] * b_a[8]) + F[13] * b_a[7])
    + F[25] * b_a[6];

  // x^2*y^3
  c_a[11] = ((F[37] * b_a[4] + F[49] * b_a[3]) + F[13] * b_a[8]) + F[25] * b_a[7];

  // x*y^4
  c_a[12] = F[49] * b_a[4] + F[25] * b_a[8];

  // y^5
  c_a[13] = (F[61] * b_a[0] + F[37] * b_a[5]) + F[1] * b_a[9];

  // x^4
  c_a[14] = (((F[61] * b_a[1] + F[37] * b_a[6]) + F[49] * b_a[5]) + F[1] * b_a
             [10]) + F[13] * b_a[9];

  // x^3*y
  c_a[15] = ((((F[61] * b_a[2] + F[37] * b_a[7]) + F[49] * b_a[6]) + F[1] * b_a
              [11]) + F[13] * b_a[10]) + F[25] * b_a[9];

  // x^2*y^2
  c_a[16] = (((F[61] * b_a[3] + F[37] * b_a[8]) + F[49] * b_a[7]) + F[13] * b_a
             [11]) + F[25] * b_a[10];

  // x*y^3
  c_a[17] = (F[61] * b_a[4] + F[49] * b_a[8]) + F[25] * b_a[11];

  // y^4
  c_a[18] = (F[61] * b_a[5] + F[1] * b_a[12]) + F[37] * b_a[9];

  // x^3
  c_a[19] = (((F[61] * b_a[6] + F[1] * b_a[13]) + F[13] * b_a[12]) + F[37] *
             b_a[10]) + F[49] * b_a[9];

  // x^2*y
  c_a[20] = (((F[61] * b_a[7] + F[13] * b_a[13]) + F[25] * b_a[12]) + F[37] *
             b_a[11]) + F[49] * b_a[10];

  // x*y^2
  c_a[21] = (F[61] * b_a[8] + F[25] * b_a[13]) + F[49] * b_a[11];

  // y^3
  c_a[22] = (F[1] * b_a[14] + F[61] * b_a[9]) + F[37] * b_a[12];

  // x^2
  c_a[23] = ((F[13] * b_a[14] + F[61] * b_a[10]) + F[37] * b_a[13]) + F[49] *
    b_a[12];

  // x*y
  c_a[24] = (F[25] * b_a[14] + F[61] * b_a[11]) + F[49] * b_a[13];

  // y^2
  c_a[25] = F[37] * b_a[14] + F[61] * b_a[12];

  // x
  c_a[26] = F[49] * b_a[14] + F[61] * b_a[13];

  // y
  c_a[27] = F[61] * b_a[14];

  // 1
  for (b_i = 0; b_i < 6; b_i++) {
    M_tmp = b_i << 2;
    M[M_tmp] = F[12 * b_i + 2];
    M[M_tmp + 2] = F[12 * b_i + 10];
    M[M_tmp + 1] = F[12 * b_i + 3];
    M[M_tmp + 3] = F[12 * b_i + 11];
  }

  b_a[0] = M[0] * M[3];

  // x^4
  b_a[1] = M[0] * M[7] + M[4] * M[3];

  // x^3*y
  b_a[2] = (M[0] * M[11] + M[4] * M[7]) + M[8] * M[3];

  // x^2*y^2
  b_a[3] = M[4] * M[11] + M[8] * M[7];

  // x*y^3
  b_a[4] = M[8] * M[11];

  // y^4
  b_a[5] = M[0] * M[15] + M[12] * M[3];

  // x^3
  b_a[6] = ((M[0] * M[19] + M[4] * M[15]) + M[12] * M[7]) + M[16] * M[3];

  // x^2*y
  b_a[7] = ((M[4] * M[19] + M[8] * M[15]) + M[12] * M[11]) + M[16] * M[7];

  // x*y^2
  b_a[8] = M[8] * M[19] + M[16] * M[11];

  // y^3
  b_a[9] = (M[0] * M[23] + M[20] * M[3]) + M[12] * M[15];

  // x^2
  b_a[10] = ((M[4] * M[23] + M[20] * M[7]) + M[12] * M[19]) + M[16] * M[15];

  // x*y
  b_a[11] = (M[8] * M[23] + M[20] * M[11]) + M[16] * M[19];

  // y^2
  b_a[12] = M[12] * M[23] + M[20] * M[15];

  // x
  b_a[13] = M[16] * M[23] + M[20] * M[19];

  // y
  b_a[14] = M[20] * M[23];

  // 1
  c_b[0] = M[2] * M[1];

  // x^4
  c_b[1] = M[2] * M[5] + M[6] * M[1];

  // x^3*y
  c_b[2] = (M[2] * M[9] + M[6] * M[5]) + M[10] * M[1];

  // x^2*y^2
  c_b[3] = M[6] * M[9] + M[10] * M[5];

  // x*y^3
  c_b[4] = M[10] * M[9];

  // y^4
  c_b[5] = M[2] * M[13] + M[14] * M[1];

  // x^3
  c_b[6] = ((M[2] * M[17] + M[6] * M[13]) + M[14] * M[5]) + M[18] * M[1];

  // x^2*y
  c_b[7] = ((M[6] * M[17] + M[10] * M[13]) + M[14] * M[9]) + M[18] * M[5];

  // x*y^2
  c_b[8] = M[10] * M[17] + M[18] * M[9];

  // y^3
  c_b[9] = (M[2] * M[21] + M[22] * M[1]) + M[14] * M[13];

  // x^2
  c_b[10] = ((M[6] * M[21] + M[22] * M[5]) + M[14] * M[17]) + M[18] * M[13];

  // x*y
  c_b[11] = (M[10] * M[21] + M[22] * M[9]) + M[18] * M[17];

  // y^2
  c_b[12] = M[14] * M[21] + M[22] * M[13];

  // x
  c_b[13] = M[18] * M[21] + M[22] * M[17];

  // y
  c_b[14] = M[22] * M[21];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  d_b[0] = F[5] * b_a[0];

  // x^6
  d_b[1] = F[5] * b_a[1] + F[17] * b_a[0];

  // x^5*y
  d_b[2] = (F[5] * b_a[2] + F[17] * b_a[1]) + F[29] * b_a[0];

  // x^4*y^2
  d_b[3] = (F[5] * b_a[3] + F[17] * b_a[2]) + F[29] * b_a[1];

  // x^3*y^3
  d_b[4] = (F[5] * b_a[4] + F[17] * b_a[3]) + F[29] * b_a[2];

  // x^2*y^4
  d_b[5] = F[17] * b_a[4] + F[29] * b_a[3];

  // x*y^5
  d_b[6] = F[29] * b_a[4];

  // y^6
  d_b[7] = F[41] * b_a[0] + F[5] * b_a[5];

  // x^5
  d_b[8] = ((F[41] * b_a[1] + F[53] * b_a[0]) + F[5] * b_a[6]) + F[17] * b_a[5];

  // x^4*y
  d_b[9] = (((F[41] * b_a[2] + F[53] * b_a[1]) + F[5] * b_a[7]) + F[17] * b_a[6])
    + F[29] * b_a[5];

  // x^3*y^2
  d_b[10] = (((F[41] * b_a[3] + F[53] * b_a[2]) + F[5] * b_a[8]) + F[17] * b_a[7])
    + F[29] * b_a[6];

  // x^2*y^3
  d_b[11] = ((F[41] * b_a[4] + F[53] * b_a[3]) + F[17] * b_a[8]) + F[29] * b_a[7];

  // x*y^4
  d_b[12] = F[53] * b_a[4] + F[29] * b_a[8];

  // y^5
  d_b[13] = (F[65] * b_a[0] + F[41] * b_a[5]) + F[5] * b_a[9];

  // x^4
  d_b[14] = (((F[65] * b_a[1] + F[41] * b_a[6]) + F[53] * b_a[5]) + F[5] * b_a
             [10]) + F[17] * b_a[9];

  // x^3*y
  d_b[15] = ((((F[65] * b_a[2] + F[41] * b_a[7]) + F[53] * b_a[6]) + F[5] * b_a
              [11]) + F[17] * b_a[10]) + F[29] * b_a[9];

  // x^2*y^2
  d_b[16] = (((F[65] * b_a[3] + F[41] * b_a[8]) + F[53] * b_a[7]) + F[17] * b_a
             [11]) + F[29] * b_a[10];

  // x*y^3
  d_b[17] = (F[65] * b_a[4] + F[53] * b_a[8]) + F[29] * b_a[11];

  // y^4
  d_b[18] = (F[65] * b_a[5] + F[5] * b_a[12]) + F[41] * b_a[9];

  // x^3
  d_b[19] = (((F[65] * b_a[6] + F[5] * b_a[13]) + F[17] * b_a[12]) + F[41] *
             b_a[10]) + F[53] * b_a[9];

  // x^2*y
  d_b[20] = (((F[65] * b_a[7] + F[17] * b_a[13]) + F[29] * b_a[12]) + F[41] *
             b_a[11]) + F[53] * b_a[10];

  // x*y^2
  d_b[21] = (F[65] * b_a[8] + F[29] * b_a[13]) + F[53] * b_a[11];

  // y^3
  d_b[22] = (F[5] * b_a[14] + F[65] * b_a[9]) + F[41] * b_a[12];

  // x^2
  d_b[23] = ((F[17] * b_a[14] + F[65] * b_a[10]) + F[41] * b_a[13]) + F[53] *
    b_a[12];

  // x*y
  d_b[24] = (F[29] * b_a[14] + F[65] * b_a[11]) + F[53] * b_a[13];

  // y^2
  d_b[25] = F[41] * b_a[14] + F[65] * b_a[12];

  // x
  d_b[26] = F[53] * b_a[14] + F[65] * b_a[13];

  // y
  d_b[27] = F[65] * b_a[14];

  // 1
  b_a[0] = F[2] * F[7];

  // x^4
  b_a[1] = F[2] * F[19] + F[14] * F[7];

  // x^3*y
  b_a[2] = (F[2] * F[31] + F[14] * F[19]) + F[26] * F[7];

  // x^2*y^2
  b_a[3] = F[14] * F[31] + F[26] * F[19];

  // x*y^3
  b_a[4] = F[26] * F[31];

  // y^4
  b_a[5] = F[2] * F[43] + F[38] * F[7];

  // x^3
  b_a[6] = ((F[2] * F[55] + F[14] * F[43]) + F[38] * F[19]) + F[50] * F[7];

  // x^2*y
  b_a[7] = ((F[14] * F[55] + F[26] * F[43]) + F[38] * F[31]) + F[50] * F[19];

  // x*y^2
  b_a[8] = F[26] * F[55] + F[50] * F[31];

  // y^3
  b_a[9] = (F[2] * F[67] + F[62] * F[7]) + F[38] * F[43];

  // x^2
  b_a[10] = ((F[14] * F[67] + F[62] * F[19]) + F[38] * F[55]) + F[50] * F[43];

  // x*y
  b_a[11] = (F[26] * F[67] + F[62] * F[31]) + F[50] * F[55];

  // y^2
  b_a[12] = F[38] * F[67] + F[62] * F[43];

  // x
  b_a[13] = F[50] * F[67] + F[62] * F[55];

  // y
  b_a[14] = F[62] * F[67];

  // 1
  c_b[0] = F[6] * F[3];

  // x^4
  c_b[1] = F[6] * F[15] + F[18] * F[3];

  // x^3*y
  c_b[2] = (F[6] * F[27] + F[18] * F[15]) + F[30] * F[3];

  // x^2*y^2
  c_b[3] = F[18] * F[27] + F[30] * F[15];

  // x*y^3
  c_b[4] = F[30] * F[27];

  // y^4
  c_b[5] = F[6] * F[39] + F[42] * F[3];

  // x^3
  c_b[6] = ((F[6] * F[51] + F[18] * F[39]) + F[42] * F[15]) + F[54] * F[3];

  // x^2*y
  c_b[7] = ((F[18] * F[51] + F[30] * F[39]) + F[42] * F[27]) + F[54] * F[15];

  // x*y^2
  c_b[8] = F[30] * F[51] + F[54] * F[27];

  // y^3
  c_b[9] = (F[6] * F[63] + F[66] * F[3]) + F[42] * F[39];

  // x^2
  c_b[10] = ((F[18] * F[63] + F[66] * F[15]) + F[42] * F[51]) + F[54] * F[39];

  // x*y
  c_b[11] = (F[30] * F[63] + F[66] * F[27]) + F[54] * F[51];

  // y^2
  c_b[12] = F[42] * F[63] + F[66] * F[39];

  // x
  c_b[13] = F[54] * F[63] + F[66] * F[51];

  // y
  c_b[14] = F[66] * F[63];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  e_b[0] = F[9] * b_a[0];

  // x^6
  e_b[1] = F[9] * b_a[1] + F[21] * b_a[0];

  // x^5*y
  e_b[2] = (F[9] * b_a[2] + F[21] * b_a[1]) + F[33] * b_a[0];

  // x^4*y^2
  e_b[3] = (F[9] * b_a[3] + F[21] * b_a[2]) + F[33] * b_a[1];

  // x^3*y^3
  e_b[4] = (F[9] * b_a[4] + F[21] * b_a[3]) + F[33] * b_a[2];

  // x^2*y^4
  e_b[5] = F[21] * b_a[4] + F[33] * b_a[3];

  // x*y^5
  e_b[6] = F[33] * b_a[4];

  // y^6
  e_b[7] = F[45] * b_a[0] + F[9] * b_a[5];

  // x^5
  e_b[8] = ((F[45] * b_a[1] + F[57] * b_a[0]) + F[9] * b_a[6]) + F[21] * b_a[5];

  // x^4*y
  e_b[9] = (((F[45] * b_a[2] + F[57] * b_a[1]) + F[9] * b_a[7]) + F[21] * b_a[6])
    + F[33] * b_a[5];

  // x^3*y^2
  e_b[10] = (((F[45] * b_a[3] + F[57] * b_a[2]) + F[9] * b_a[8]) + F[21] * b_a[7])
    + F[33] * b_a[6];

  // x^2*y^3
  e_b[11] = ((F[45] * b_a[4] + F[57] * b_a[3]) + F[21] * b_a[8]) + F[33] * b_a[7];

  // x*y^4
  e_b[12] = F[57] * b_a[4] + F[33] * b_a[8];

  // y^5
  e_b[13] = (F[69] * b_a[0] + F[45] * b_a[5]) + F[9] * b_a[9];

  // x^4
  e_b[14] = (((F[69] * b_a[1] + F[45] * b_a[6]) + F[57] * b_a[5]) + F[9] * b_a
             [10]) + F[21] * b_a[9];

  // x^3*y
  e_b[15] = ((((F[69] * b_a[2] + F[45] * b_a[7]) + F[57] * b_a[6]) + F[9] * b_a
              [11]) + F[21] * b_a[10]) + F[33] * b_a[9];

  // x^2*y^2
  e_b[16] = (((F[69] * b_a[3] + F[45] * b_a[8]) + F[57] * b_a[7]) + F[21] * b_a
             [11]) + F[33] * b_a[10];

  // x*y^3
  e_b[17] = (F[69] * b_a[4] + F[57] * b_a[8]) + F[33] * b_a[11];

  // y^4
  e_b[18] = (F[69] * b_a[5] + F[9] * b_a[12]) + F[45] * b_a[9];

  // x^3
  e_b[19] = (((F[69] * b_a[6] + F[9] * b_a[13]) + F[21] * b_a[12]) + F[45] *
             b_a[10]) + F[57] * b_a[9];

  // x^2*y
  e_b[20] = (((F[69] * b_a[7] + F[21] * b_a[13]) + F[33] * b_a[12]) + F[45] *
             b_a[11]) + F[57] * b_a[10];

  // x*y^2
  e_b[21] = (F[69] * b_a[8] + F[33] * b_a[13]) + F[57] * b_a[11];

  // y^3
  e_b[22] = (F[9] * b_a[14] + F[69] * b_a[9]) + F[45] * b_a[12];

  // x^2
  e_b[23] = ((F[21] * b_a[14] + F[69] * b_a[10]) + F[45] * b_a[13]) + F[57] *
    b_a[12];

  // x*y
  e_b[24] = (F[33] * b_a[14] + F[69] * b_a[11]) + F[57] * b_a[13];

  // y^2
  e_b[25] = F[45] * b_a[14] + F[69] * b_a[12];

  // x
  e_b[26] = F[57] * b_a[14] + F[69] * b_a[13];

  // y
  e_b[27] = F[69] * b_a[14];

  // 1
  for (b_i = 0; b_i < 28; b_i++) {
    G[b_i << 2] = (c_a[b_i] - d_b[b_i]) + e_b[b_i];
  }

  for (b_i = 0; b_i < 6; b_i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * b_i;
      b_M_tmp = 3 * i1 + 9 * b_i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 2];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  b_a[0] = b_M[4] * b_M[8];

  // x^4
  b_a[1] = b_M[4] * b_M[17] + b_M[13] * b_M[8];

  // x^3*y
  b_a[2] = (b_M[4] * b_M[26] + b_M[13] * b_M[17]) + b_M[22] * b_M[8];

  // x^2*y^2
  b_a[3] = b_M[13] * b_M[26] + b_M[22] * b_M[17];

  // x*y^3
  b_a[4] = b_M[22] * b_M[26];

  // y^4
  b_a[5] = b_M[4] * b_M[35] + b_M[31] * b_M[8];

  // x^3
  b_a[6] = ((b_M[4] * b_M[44] + b_M[13] * b_M[35]) + b_M[31] * b_M[17]) + b_M[40]
    * b_M[8];

  // x^2*y
  b_a[7] = ((b_M[13] * b_M[44] + b_M[22] * b_M[35]) + b_M[31] * b_M[26]) + b_M
    [40] * b_M[17];

  // x*y^2
  b_a[8] = b_M[22] * b_M[44] + b_M[40] * b_M[26];

  // y^3
  b_a[9] = (b_M[4] * b_M[53] + b_M[49] * b_M[8]) + b_M[31] * b_M[35];

  // x^2
  b_a[10] = ((b_M[13] * b_M[53] + b_M[49] * b_M[17]) + b_M[31] * b_M[44]) + b_M
    [40] * b_M[35];

  // x*y
  b_a[11] = (b_M[22] * b_M[53] + b_M[49] * b_M[26]) + b_M[40] * b_M[44];

  // y^2
  b_a[12] = b_M[31] * b_M[53] + b_M[49] * b_M[35];

  // x
  b_a[13] = b_M[40] * b_M[53] + b_M[49] * b_M[44];

  // y
  b_a[14] = b_M[49] * b_M[53];

  // 1
  c_b[0] = b_M[7] * b_M[5];

  // x^4
  c_b[1] = b_M[7] * b_M[14] + b_M[16] * b_M[5];

  // x^3*y
  c_b[2] = (b_M[7] * b_M[23] + b_M[16] * b_M[14]) + b_M[25] * b_M[5];

  // x^2*y^2
  c_b[3] = b_M[16] * b_M[23] + b_M[25] * b_M[14];

  // x*y^3
  c_b[4] = b_M[25] * b_M[23];

  // y^4
  c_b[5] = b_M[7] * b_M[32] + b_M[34] * b_M[5];

  // x^3
  c_b[6] = ((b_M[7] * b_M[41] + b_M[16] * b_M[32]) + b_M[34] * b_M[14]) + b_M[43]
    * b_M[5];

  // x^2*y
  c_b[7] = ((b_M[16] * b_M[41] + b_M[25] * b_M[32]) + b_M[34] * b_M[23]) + b_M
    [43] * b_M[14];

  // x*y^2
  c_b[8] = b_M[25] * b_M[41] + b_M[43] * b_M[23];

  // y^3
  c_b[9] = (b_M[7] * b_M[50] + b_M[52] * b_M[5]) + b_M[34] * b_M[32];

  // x^2
  c_b[10] = ((b_M[16] * b_M[50] + b_M[52] * b_M[14]) + b_M[34] * b_M[41]) + b_M
    [43] * b_M[32];

  // x*y
  c_b[11] = (b_M[25] * b_M[50] + b_M[52] * b_M[23]) + b_M[43] * b_M[41];

  // y^2
  c_b[12] = b_M[34] * b_M[50] + b_M[52] * b_M[32];

  // x
  c_b[13] = b_M[43] * b_M[50] + b_M[52] * b_M[41];

  // y
  c_b[14] = b_M[52] * b_M[50];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  c_a[0] = b_M[0] * b_a[0];

  // x^6
  c_a[1] = b_M[0] * b_a[1] + b_M[9] * b_a[0];

  // x^5*y
  c_a[2] = (b_M[0] * b_a[2] + b_M[9] * b_a[1]) + b_M[18] * b_a[0];

  // x^4*y^2
  c_a[3] = (b_M[0] * b_a[3] + b_M[9] * b_a[2]) + b_M[18] * b_a[1];

  // x^3*y^3
  c_a[4] = (b_M[0] * b_a[4] + b_M[9] * b_a[3]) + b_M[18] * b_a[2];

  // x^2*y^4
  c_a[5] = b_M[9] * b_a[4] + b_M[18] * b_a[3];

  // x*y^5
  c_a[6] = b_M[18] * b_a[4];

  // y^6
  c_a[7] = b_M[27] * b_a[0] + b_M[0] * b_a[5];

  // x^5
  c_a[8] = ((b_M[27] * b_a[1] + b_M[36] * b_a[0]) + b_M[0] * b_a[6]) + b_M[9] *
    b_a[5];

  // x^4*y
  c_a[9] = (((b_M[27] * b_a[2] + b_M[36] * b_a[1]) + b_M[0] * b_a[7]) + b_M[9] *
            b_a[6]) + b_M[18] * b_a[5];

  // x^3*y^2
  c_a[10] = (((b_M[27] * b_a[3] + b_M[36] * b_a[2]) + b_M[0] * b_a[8]) + b_M[9] *
             b_a[7]) + b_M[18] * b_a[6];

  // x^2*y^3
  c_a[11] = ((b_M[27] * b_a[4] + b_M[36] * b_a[3]) + b_M[9] * b_a[8]) + b_M[18] *
    b_a[7];

  // x*y^4
  c_a[12] = b_M[36] * b_a[4] + b_M[18] * b_a[8];

  // y^5
  c_a[13] = (b_M[45] * b_a[0] + b_M[27] * b_a[5]) + b_M[0] * b_a[9];

  // x^4
  c_a[14] = (((b_M[45] * b_a[1] + b_M[27] * b_a[6]) + b_M[36] * b_a[5]) + b_M[0]
             * b_a[10]) + b_M[9] * b_a[9];

  // x^3*y
  c_a[15] = ((((b_M[45] * b_a[2] + b_M[27] * b_a[7]) + b_M[36] * b_a[6]) + b_M[0]
              * b_a[11]) + b_M[9] * b_a[10]) + b_M[18] * b_a[9];

  // x^2*y^2
  c_a[16] = (((b_M[45] * b_a[3] + b_M[27] * b_a[8]) + b_M[36] * b_a[7]) + b_M[9]
             * b_a[11]) + b_M[18] * b_a[10];

  // x*y^3
  c_a[17] = (b_M[45] * b_a[4] + b_M[36] * b_a[8]) + b_M[18] * b_a[11];

  // y^4
  c_a[18] = (b_M[45] * b_a[5] + b_M[0] * b_a[12]) + b_M[27] * b_a[9];

  // x^3
  c_a[19] = (((b_M[45] * b_a[6] + b_M[0] * b_a[13]) + b_M[9] * b_a[12]) + b_M[27]
             * b_a[10]) + b_M[36] * b_a[9];

  // x^2*y
  c_a[20] = (((b_M[45] * b_a[7] + b_M[9] * b_a[13]) + b_M[18] * b_a[12]) + b_M
             [27] * b_a[11]) + b_M[36] * b_a[10];

  // x*y^2
  c_a[21] = (b_M[45] * b_a[8] + b_M[18] * b_a[13]) + b_M[36] * b_a[11];

  // y^3
  c_a[22] = (b_M[0] * b_a[14] + b_M[45] * b_a[9]) + b_M[27] * b_a[12];

  // x^2
  c_a[23] = ((b_M[9] * b_a[14] + b_M[45] * b_a[10]) + b_M[27] * b_a[13]) + b_M
    [36] * b_a[12];

  // x*y
  c_a[24] = (b_M[18] * b_a[14] + b_M[45] * b_a[11]) + b_M[36] * b_a[13];

  // y^2
  c_a[25] = b_M[27] * b_a[14] + b_M[45] * b_a[12];

  // x
  c_a[26] = b_M[36] * b_a[14] + b_M[45] * b_a[13];

  // y
  c_a[27] = b_M[45] * b_a[14];

  // 1
  for (b_i = 0; b_i < 6; b_i++) {
    M_tmp = b_i << 2;
    M[M_tmp] = b_M[9 * b_i + 1];
    M[M_tmp + 2] = b_M[9 * b_i + 7];
    M[M_tmp + 1] = b_M[9 * b_i + 2];
    M[M_tmp + 3] = b_M[9 * b_i + 8];
  }

  b_a[0] = M[0] * M[3];

  // x^4
  b_a[1] = M[0] * M[7] + M[4] * M[3];

  // x^3*y
  b_a[2] = (M[0] * M[11] + M[4] * M[7]) + M[8] * M[3];

  // x^2*y^2
  b_a[3] = M[4] * M[11] + M[8] * M[7];

  // x*y^3
  b_a[4] = M[8] * M[11];

  // y^4
  b_a[5] = M[0] * M[15] + M[12] * M[3];

  // x^3
  b_a[6] = ((M[0] * M[19] + M[4] * M[15]) + M[12] * M[7]) + M[16] * M[3];

  // x^2*y
  b_a[7] = ((M[4] * M[19] + M[8] * M[15]) + M[12] * M[11]) + M[16] * M[7];

  // x*y^2
  b_a[8] = M[8] * M[19] + M[16] * M[11];

  // y^3
  b_a[9] = (M[0] * M[23] + M[20] * M[3]) + M[12] * M[15];

  // x^2
  b_a[10] = ((M[4] * M[23] + M[20] * M[7]) + M[12] * M[19]) + M[16] * M[15];

  // x*y
  b_a[11] = (M[8] * M[23] + M[20] * M[11]) + M[16] * M[19];

  // y^2
  b_a[12] = M[12] * M[23] + M[20] * M[15];

  // x
  b_a[13] = M[16] * M[23] + M[20] * M[19];

  // y
  b_a[14] = M[20] * M[23];

  // 1
  c_b[0] = M[2] * M[1];

  // x^4
  c_b[1] = M[2] * M[5] + M[6] * M[1];

  // x^3*y
  c_b[2] = (M[2] * M[9] + M[6] * M[5]) + M[10] * M[1];

  // x^2*y^2
  c_b[3] = M[6] * M[9] + M[10] * M[5];

  // x*y^3
  c_b[4] = M[10] * M[9];

  // y^4
  c_b[5] = M[2] * M[13] + M[14] * M[1];

  // x^3
  c_b[6] = ((M[2] * M[17] + M[6] * M[13]) + M[14] * M[5]) + M[18] * M[1];

  // x^2*y
  c_b[7] = ((M[6] * M[17] + M[10] * M[13]) + M[14] * M[9]) + M[18] * M[5];

  // x*y^2
  c_b[8] = M[10] * M[17] + M[18] * M[9];

  // y^3
  c_b[9] = (M[2] * M[21] + M[22] * M[1]) + M[14] * M[13];

  // x^2
  c_b[10] = ((M[6] * M[21] + M[22] * M[5]) + M[14] * M[17]) + M[18] * M[13];

  // x*y
  c_b[11] = (M[10] * M[21] + M[22] * M[9]) + M[18] * M[17];

  // y^2
  c_b[12] = M[14] * M[21] + M[22] * M[13];

  // x
  c_b[13] = M[18] * M[21] + M[22] * M[17];

  // y
  c_b[14] = M[22] * M[21];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  d_b[0] = b_M[3] * b_a[0];

  // x^6
  d_b[1] = b_M[3] * b_a[1] + b_M[12] * b_a[0];

  // x^5*y
  d_b[2] = (b_M[3] * b_a[2] + b_M[12] * b_a[1]) + b_M[21] * b_a[0];

  // x^4*y^2
  d_b[3] = (b_M[3] * b_a[3] + b_M[12] * b_a[2]) + b_M[21] * b_a[1];

  // x^3*y^3
  d_b[4] = (b_M[3] * b_a[4] + b_M[12] * b_a[3]) + b_M[21] * b_a[2];

  // x^2*y^4
  d_b[5] = b_M[12] * b_a[4] + b_M[21] * b_a[3];

  // x*y^5
  d_b[6] = b_M[21] * b_a[4];

  // y^6
  d_b[7] = b_M[30] * b_a[0] + b_M[3] * b_a[5];

  // x^5
  d_b[8] = ((b_M[30] * b_a[1] + b_M[39] * b_a[0]) + b_M[3] * b_a[6]) + b_M[12] *
    b_a[5];

  // x^4*y
  d_b[9] = (((b_M[30] * b_a[2] + b_M[39] * b_a[1]) + b_M[3] * b_a[7]) + b_M[12] *
            b_a[6]) + b_M[21] * b_a[5];

  // x^3*y^2
  d_b[10] = (((b_M[30] * b_a[3] + b_M[39] * b_a[2]) + b_M[3] * b_a[8]) + b_M[12]
             * b_a[7]) + b_M[21] * b_a[6];

  // x^2*y^3
  d_b[11] = ((b_M[30] * b_a[4] + b_M[39] * b_a[3]) + b_M[12] * b_a[8]) + b_M[21]
    * b_a[7];

  // x*y^4
  d_b[12] = b_M[39] * b_a[4] + b_M[21] * b_a[8];

  // y^5
  d_b[13] = (b_M[48] * b_a[0] + b_M[30] * b_a[5]) + b_M[3] * b_a[9];

  // x^4
  d_b[14] = (((b_M[48] * b_a[1] + b_M[30] * b_a[6]) + b_M[39] * b_a[5]) + b_M[3]
             * b_a[10]) + b_M[12] * b_a[9];

  // x^3*y
  d_b[15] = ((((b_M[48] * b_a[2] + b_M[30] * b_a[7]) + b_M[39] * b_a[6]) + b_M[3]
              * b_a[11]) + b_M[12] * b_a[10]) + b_M[21] * b_a[9];

  // x^2*y^2
  d_b[16] = (((b_M[48] * b_a[3] + b_M[30] * b_a[8]) + b_M[39] * b_a[7]) + b_M[12]
             * b_a[11]) + b_M[21] * b_a[10];

  // x*y^3
  d_b[17] = (b_M[48] * b_a[4] + b_M[39] * b_a[8]) + b_M[21] * b_a[11];

  // y^4
  d_b[18] = (b_M[48] * b_a[5] + b_M[3] * b_a[12]) + b_M[30] * b_a[9];

  // x^3
  d_b[19] = (((b_M[48] * b_a[6] + b_M[3] * b_a[13]) + b_M[12] * b_a[12]) + b_M
             [30] * b_a[10]) + b_M[39] * b_a[9];

  // x^2*y
  d_b[20] = (((b_M[48] * b_a[7] + b_M[12] * b_a[13]) + b_M[21] * b_a[12]) + b_M
             [30] * b_a[11]) + b_M[39] * b_a[10];

  // x*y^2
  d_b[21] = (b_M[48] * b_a[8] + b_M[21] * b_a[13]) + b_M[39] * b_a[11];

  // y^3
  d_b[22] = (b_M[3] * b_a[14] + b_M[48] * b_a[9]) + b_M[30] * b_a[12];

  // x^2
  d_b[23] = ((b_M[12] * b_a[14] + b_M[48] * b_a[10]) + b_M[30] * b_a[13]) + b_M
    [39] * b_a[12];

  // x*y
  d_b[24] = (b_M[21] * b_a[14] + b_M[48] * b_a[11]) + b_M[39] * b_a[13];

  // y^2
  d_b[25] = b_M[30] * b_a[14] + b_M[48] * b_a[12];

  // x
  d_b[26] = b_M[39] * b_a[14] + b_M[48] * b_a[13];

  // y
  d_b[27] = b_M[48] * b_a[14];

  // 1
  b_a[0] = b_M[1] * b_M[5];

  // x^4
  b_a[1] = b_M[1] * b_M[14] + b_M[10] * b_M[5];

  // x^3*y
  b_a[2] = (b_M[1] * b_M[23] + b_M[10] * b_M[14]) + b_M[19] * b_M[5];

  // x^2*y^2
  b_a[3] = b_M[10] * b_M[23] + b_M[19] * b_M[14];

  // x*y^3
  b_a[4] = b_M[19] * b_M[23];

  // y^4
  b_a[5] = b_M[1] * b_M[32] + b_M[28] * b_M[5];

  // x^3
  b_a[6] = ((b_M[1] * b_M[41] + b_M[10] * b_M[32]) + b_M[28] * b_M[14]) + b_M[37]
    * b_M[5];

  // x^2*y
  b_a[7] = ((b_M[10] * b_M[41] + b_M[19] * b_M[32]) + b_M[28] * b_M[23]) + b_M
    [37] * b_M[14];

  // x*y^2
  b_a[8] = b_M[19] * b_M[41] + b_M[37] * b_M[23];

  // y^3
  b_a[9] = (b_M[1] * b_M[50] + b_M[46] * b_M[5]) + b_M[28] * b_M[32];

  // x^2
  b_a[10] = ((b_M[10] * b_M[50] + b_M[46] * b_M[14]) + b_M[28] * b_M[41]) + b_M
    [37] * b_M[32];

  // x*y
  b_a[11] = (b_M[19] * b_M[50] + b_M[46] * b_M[23]) + b_M[37] * b_M[41];

  // y^2
  b_a[12] = b_M[28] * b_M[50] + b_M[46] * b_M[32];

  // x
  b_a[13] = b_M[37] * b_M[50] + b_M[46] * b_M[41];

  // y
  b_a[14] = b_M[46] * b_M[50];

  // 1
  c_b[0] = b_M[4] * b_M[2];

  // x^4
  c_b[1] = b_M[4] * b_M[11] + b_M[13] * b_M[2];

  // x^3*y
  c_b[2] = (b_M[4] * b_M[20] + b_M[13] * b_M[11]) + b_M[22] * b_M[2];

  // x^2*y^2
  c_b[3] = b_M[13] * b_M[20] + b_M[22] * b_M[11];

  // x*y^3
  c_b[4] = b_M[22] * b_M[20];

  // y^4
  c_b[5] = b_M[4] * b_M[29] + b_M[31] * b_M[2];

  // x^3
  c_b[6] = ((b_M[4] * b_M[38] + b_M[13] * b_M[29]) + b_M[31] * b_M[11]) + b_M[40]
    * b_M[2];

  // x^2*y
  c_b[7] = ((b_M[13] * b_M[38] + b_M[22] * b_M[29]) + b_M[31] * b_M[20]) + b_M
    [40] * b_M[11];

  // x*y^2
  c_b[8] = b_M[22] * b_M[38] + b_M[40] * b_M[20];

  // y^3
  c_b[9] = (b_M[4] * b_M[47] + b_M[49] * b_M[2]) + b_M[31] * b_M[29];

  // x^2
  c_b[10] = ((b_M[13] * b_M[47] + b_M[49] * b_M[11]) + b_M[31] * b_M[38]) + b_M
    [40] * b_M[29];

  // x*y
  c_b[11] = (b_M[22] * b_M[47] + b_M[49] * b_M[20]) + b_M[40] * b_M[38];

  // y^2
  c_b[12] = b_M[31] * b_M[47] + b_M[49] * b_M[29];

  // x
  c_b[13] = b_M[40] * b_M[47] + b_M[49] * b_M[38];

  // y
  c_b[14] = b_M[49] * b_M[47];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  e_b[0] = b_M[6] * b_a[0];

  // x^6
  e_b[1] = b_M[6] * b_a[1] + b_M[15] * b_a[0];

  // x^5*y
  e_b[2] = (b_M[6] * b_a[2] + b_M[15] * b_a[1]) + b_M[24] * b_a[0];

  // x^4*y^2
  e_b[3] = (b_M[6] * b_a[3] + b_M[15] * b_a[2]) + b_M[24] * b_a[1];

  // x^3*y^3
  e_b[4] = (b_M[6] * b_a[4] + b_M[15] * b_a[3]) + b_M[24] * b_a[2];

  // x^2*y^4
  e_b[5] = b_M[15] * b_a[4] + b_M[24] * b_a[3];

  // x*y^5
  e_b[6] = b_M[24] * b_a[4];

  // y^6
  e_b[7] = b_M[33] * b_a[0] + b_M[6] * b_a[5];

  // x^5
  e_b[8] = ((b_M[33] * b_a[1] + b_M[42] * b_a[0]) + b_M[6] * b_a[6]) + b_M[15] *
    b_a[5];

  // x^4*y
  e_b[9] = (((b_M[33] * b_a[2] + b_M[42] * b_a[1]) + b_M[6] * b_a[7]) + b_M[15] *
            b_a[6]) + b_M[24] * b_a[5];

  // x^3*y^2
  e_b[10] = (((b_M[33] * b_a[3] + b_M[42] * b_a[2]) + b_M[6] * b_a[8]) + b_M[15]
             * b_a[7]) + b_M[24] * b_a[6];

  // x^2*y^3
  e_b[11] = ((b_M[33] * b_a[4] + b_M[42] * b_a[3]) + b_M[15] * b_a[8]) + b_M[24]
    * b_a[7];

  // x*y^4
  e_b[12] = b_M[42] * b_a[4] + b_M[24] * b_a[8];

  // y^5
  e_b[13] = (b_M[51] * b_a[0] + b_M[33] * b_a[5]) + b_M[6] * b_a[9];

  // x^4
  e_b[14] = (((b_M[51] * b_a[1] + b_M[33] * b_a[6]) + b_M[42] * b_a[5]) + b_M[6]
             * b_a[10]) + b_M[15] * b_a[9];

  // x^3*y
  e_b[15] = ((((b_M[51] * b_a[2] + b_M[33] * b_a[7]) + b_M[42] * b_a[6]) + b_M[6]
              * b_a[11]) + b_M[15] * b_a[10]) + b_M[24] * b_a[9];

  // x^2*y^2
  e_b[16] = (((b_M[51] * b_a[3] + b_M[33] * b_a[8]) + b_M[42] * b_a[7]) + b_M[15]
             * b_a[11]) + b_M[24] * b_a[10];

  // x*y^3
  e_b[17] = (b_M[51] * b_a[4] + b_M[42] * b_a[8]) + b_M[24] * b_a[11];

  // y^4
  e_b[18] = (b_M[51] * b_a[5] + b_M[6] * b_a[12]) + b_M[33] * b_a[9];

  // x^3
  e_b[19] = (((b_M[51] * b_a[6] + b_M[6] * b_a[13]) + b_M[15] * b_a[12]) + b_M
             [33] * b_a[10]) + b_M[42] * b_a[9];

  // x^2*y
  e_b[20] = (((b_M[51] * b_a[7] + b_M[15] * b_a[13]) + b_M[24] * b_a[12]) + b_M
             [33] * b_a[11]) + b_M[42] * b_a[10];

  // x*y^2
  e_b[21] = (b_M[51] * b_a[8] + b_M[24] * b_a[13]) + b_M[42] * b_a[11];

  // y^3
  e_b[22] = (b_M[6] * b_a[14] + b_M[51] * b_a[9]) + b_M[33] * b_a[12];

  // x^2
  e_b[23] = ((b_M[15] * b_a[14] + b_M[51] * b_a[10]) + b_M[33] * b_a[13]) + b_M
    [42] * b_a[12];

  // x*y
  e_b[24] = (b_M[24] * b_a[14] + b_M[51] * b_a[11]) + b_M[42] * b_a[13];

  // y^2
  e_b[25] = b_M[33] * b_a[14] + b_M[51] * b_a[12];

  // x
  e_b[26] = b_M[42] * b_a[14] + b_M[51] * b_a[13];

  // y
  e_b[27] = b_M[51] * b_a[14];

  // 1
  for (b_i = 0; b_i < 28; b_i++) {
    G[(b_i << 2) + 1] = (c_a[b_i] - d_b[b_i]) + e_b[b_i];
  }

  for (b_i = 0; b_i < 6; b_i++) {
    for (i1 = 0; i1 < 3; i1++) {
      M_tmp = (i1 << 2) + 12 * b_i;
      b_M_tmp = 3 * i1 + 9 * b_i;
      b_M[b_M_tmp] = F[M_tmp];
      b_M[b_M_tmp + 1] = F[M_tmp + 1];
      b_M[b_M_tmp + 2] = F[M_tmp + 3];
    }
  }

  b_a[0] = b_M[4] * b_M[8];

  // x^4
  b_a[1] = b_M[4] * b_M[17] + b_M[13] * b_M[8];

  // x^3*y
  b_a[2] = (b_M[4] * b_M[26] + b_M[13] * b_M[17]) + b_M[22] * b_M[8];

  // x^2*y^2
  b_a[3] = b_M[13] * b_M[26] + b_M[22] * b_M[17];

  // x*y^3
  b_a[4] = b_M[22] * b_M[26];

  // y^4
  b_a[5] = b_M[4] * b_M[35] + b_M[31] * b_M[8];

  // x^3
  b_a[6] = ((b_M[4] * b_M[44] + b_M[13] * b_M[35]) + b_M[31] * b_M[17]) + b_M[40]
    * b_M[8];

  // x^2*y
  b_a[7] = ((b_M[13] * b_M[44] + b_M[22] * b_M[35]) + b_M[31] * b_M[26]) + b_M
    [40] * b_M[17];

  // x*y^2
  b_a[8] = b_M[22] * b_M[44] + b_M[40] * b_M[26];

  // y^3
  b_a[9] = (b_M[4] * b_M[53] + b_M[49] * b_M[8]) + b_M[31] * b_M[35];

  // x^2
  b_a[10] = ((b_M[13] * b_M[53] + b_M[49] * b_M[17]) + b_M[31] * b_M[44]) + b_M
    [40] * b_M[35];

  // x*y
  b_a[11] = (b_M[22] * b_M[53] + b_M[49] * b_M[26]) + b_M[40] * b_M[44];

  // y^2
  b_a[12] = b_M[31] * b_M[53] + b_M[49] * b_M[35];

  // x
  b_a[13] = b_M[40] * b_M[53] + b_M[49] * b_M[44];

  // y
  b_a[14] = b_M[49] * b_M[53];

  // 1
  c_b[0] = b_M[7] * b_M[5];

  // x^4
  c_b[1] = b_M[7] * b_M[14] + b_M[16] * b_M[5];

  // x^3*y
  c_b[2] = (b_M[7] * b_M[23] + b_M[16] * b_M[14]) + b_M[25] * b_M[5];

  // x^2*y^2
  c_b[3] = b_M[16] * b_M[23] + b_M[25] * b_M[14];

  // x*y^3
  c_b[4] = b_M[25] * b_M[23];

  // y^4
  c_b[5] = b_M[7] * b_M[32] + b_M[34] * b_M[5];

  // x^3
  c_b[6] = ((b_M[7] * b_M[41] + b_M[16] * b_M[32]) + b_M[34] * b_M[14]) + b_M[43]
    * b_M[5];

  // x^2*y
  c_b[7] = ((b_M[16] * b_M[41] + b_M[25] * b_M[32]) + b_M[34] * b_M[23]) + b_M
    [43] * b_M[14];

  // x*y^2
  c_b[8] = b_M[25] * b_M[41] + b_M[43] * b_M[23];

  // y^3
  c_b[9] = (b_M[7] * b_M[50] + b_M[52] * b_M[5]) + b_M[34] * b_M[32];

  // x^2
  c_b[10] = ((b_M[16] * b_M[50] + b_M[52] * b_M[14]) + b_M[34] * b_M[41]) + b_M
    [43] * b_M[32];

  // x*y
  c_b[11] = (b_M[25] * b_M[50] + b_M[52] * b_M[23]) + b_M[43] * b_M[41];

  // y^2
  c_b[12] = b_M[34] * b_M[50] + b_M[52] * b_M[32];

  // x
  c_b[13] = b_M[43] * b_M[50] + b_M[52] * b_M[41];

  // y
  c_b[14] = b_M[52] * b_M[50];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  c_a[0] = b_M[0] * b_a[0];

  // x^6
  c_a[1] = b_M[0] * b_a[1] + b_M[9] * b_a[0];

  // x^5*y
  c_a[2] = (b_M[0] * b_a[2] + b_M[9] * b_a[1]) + b_M[18] * b_a[0];

  // x^4*y^2
  c_a[3] = (b_M[0] * b_a[3] + b_M[9] * b_a[2]) + b_M[18] * b_a[1];

  // x^3*y^3
  c_a[4] = (b_M[0] * b_a[4] + b_M[9] * b_a[3]) + b_M[18] * b_a[2];

  // x^2*y^4
  c_a[5] = b_M[9] * b_a[4] + b_M[18] * b_a[3];

  // x*y^5
  c_a[6] = b_M[18] * b_a[4];

  // y^6
  c_a[7] = b_M[27] * b_a[0] + b_M[0] * b_a[5];

  // x^5
  c_a[8] = ((b_M[27] * b_a[1] + b_M[36] * b_a[0]) + b_M[0] * b_a[6]) + b_M[9] *
    b_a[5];

  // x^4*y
  c_a[9] = (((b_M[27] * b_a[2] + b_M[36] * b_a[1]) + b_M[0] * b_a[7]) + b_M[9] *
            b_a[6]) + b_M[18] * b_a[5];

  // x^3*y^2
  c_a[10] = (((b_M[27] * b_a[3] + b_M[36] * b_a[2]) + b_M[0] * b_a[8]) + b_M[9] *
             b_a[7]) + b_M[18] * b_a[6];

  // x^2*y^3
  c_a[11] = ((b_M[27] * b_a[4] + b_M[36] * b_a[3]) + b_M[9] * b_a[8]) + b_M[18] *
    b_a[7];

  // x*y^4
  c_a[12] = b_M[36] * b_a[4] + b_M[18] * b_a[8];

  // y^5
  c_a[13] = (b_M[45] * b_a[0] + b_M[27] * b_a[5]) + b_M[0] * b_a[9];

  // x^4
  c_a[14] = (((b_M[45] * b_a[1] + b_M[27] * b_a[6]) + b_M[36] * b_a[5]) + b_M[0]
             * b_a[10]) + b_M[9] * b_a[9];

  // x^3*y
  c_a[15] = ((((b_M[45] * b_a[2] + b_M[27] * b_a[7]) + b_M[36] * b_a[6]) + b_M[0]
              * b_a[11]) + b_M[9] * b_a[10]) + b_M[18] * b_a[9];

  // x^2*y^2
  c_a[16] = (((b_M[45] * b_a[3] + b_M[27] * b_a[8]) + b_M[36] * b_a[7]) + b_M[9]
             * b_a[11]) + b_M[18] * b_a[10];

  // x*y^3
  c_a[17] = (b_M[45] * b_a[4] + b_M[36] * b_a[8]) + b_M[18] * b_a[11];

  // y^4
  c_a[18] = (b_M[45] * b_a[5] + b_M[0] * b_a[12]) + b_M[27] * b_a[9];

  // x^3
  c_a[19] = (((b_M[45] * b_a[6] + b_M[0] * b_a[13]) + b_M[9] * b_a[12]) + b_M[27]
             * b_a[10]) + b_M[36] * b_a[9];

  // x^2*y
  c_a[20] = (((b_M[45] * b_a[7] + b_M[9] * b_a[13]) + b_M[18] * b_a[12]) + b_M
             [27] * b_a[11]) + b_M[36] * b_a[10];

  // x*y^2
  c_a[21] = (b_M[45] * b_a[8] + b_M[18] * b_a[13]) + b_M[36] * b_a[11];

  // y^3
  c_a[22] = (b_M[0] * b_a[14] + b_M[45] * b_a[9]) + b_M[27] * b_a[12];

  // x^2
  c_a[23] = ((b_M[9] * b_a[14] + b_M[45] * b_a[10]) + b_M[27] * b_a[13]) + b_M
    [36] * b_a[12];

  // x*y
  c_a[24] = (b_M[18] * b_a[14] + b_M[45] * b_a[11]) + b_M[36] * b_a[13];

  // y^2
  c_a[25] = b_M[27] * b_a[14] + b_M[45] * b_a[12];

  // x
  c_a[26] = b_M[36] * b_a[14] + b_M[45] * b_a[13];

  // y
  c_a[27] = b_M[45] * b_a[14];

  // 1
  for (b_i = 0; b_i < 6; b_i++) {
    M_tmp = b_i << 2;
    M[M_tmp] = b_M[9 * b_i + 1];
    M[M_tmp + 2] = b_M[9 * b_i + 7];
    M[M_tmp + 1] = b_M[9 * b_i + 2];
    M[M_tmp + 3] = b_M[9 * b_i + 8];
  }

  b_a[0] = M[0] * M[3];

  // x^4
  b_a[1] = M[0] * M[7] + M[4] * M[3];

  // x^3*y
  b_a[2] = (M[0] * M[11] + M[4] * M[7]) + M[8] * M[3];

  // x^2*y^2
  b_a[3] = M[4] * M[11] + M[8] * M[7];

  // x*y^3
  b_a[4] = M[8] * M[11];

  // y^4
  b_a[5] = M[0] * M[15] + M[12] * M[3];

  // x^3
  b_a[6] = ((M[0] * M[19] + M[4] * M[15]) + M[12] * M[7]) + M[16] * M[3];

  // x^2*y
  b_a[7] = ((M[4] * M[19] + M[8] * M[15]) + M[12] * M[11]) + M[16] * M[7];

  // x*y^2
  b_a[8] = M[8] * M[19] + M[16] * M[11];

  // y^3
  b_a[9] = (M[0] * M[23] + M[20] * M[3]) + M[12] * M[15];

  // x^2
  b_a[10] = ((M[4] * M[23] + M[20] * M[7]) + M[12] * M[19]) + M[16] * M[15];

  // x*y
  b_a[11] = (M[8] * M[23] + M[20] * M[11]) + M[16] * M[19];

  // y^2
  b_a[12] = M[12] * M[23] + M[20] * M[15];

  // x
  b_a[13] = M[16] * M[23] + M[20] * M[19];

  // y
  b_a[14] = M[20] * M[23];

  // 1
  c_b[0] = M[2] * M[1];

  // x^4
  c_b[1] = M[2] * M[5] + M[6] * M[1];

  // x^3*y
  c_b[2] = (M[2] * M[9] + M[6] * M[5]) + M[10] * M[1];

  // x^2*y^2
  c_b[3] = M[6] * M[9] + M[10] * M[5];

  // x*y^3
  c_b[4] = M[10] * M[9];

  // y^4
  c_b[5] = M[2] * M[13] + M[14] * M[1];

  // x^3
  c_b[6] = ((M[2] * M[17] + M[6] * M[13]) + M[14] * M[5]) + M[18] * M[1];

  // x^2*y
  c_b[7] = ((M[6] * M[17] + M[10] * M[13]) + M[14] * M[9]) + M[18] * M[5];

  // x*y^2
  c_b[8] = M[10] * M[17] + M[18] * M[9];

  // y^3
  c_b[9] = (M[2] * M[21] + M[22] * M[1]) + M[14] * M[13];

  // x^2
  c_b[10] = ((M[6] * M[21] + M[22] * M[5]) + M[14] * M[17]) + M[18] * M[13];

  // x*y
  c_b[11] = (M[10] * M[21] + M[22] * M[9]) + M[18] * M[17];

  // y^2
  c_b[12] = M[14] * M[21] + M[22] * M[13];

  // x
  c_b[13] = M[18] * M[21] + M[22] * M[17];

  // y
  c_b[14] = M[22] * M[21];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  d_b[0] = b_M[3] * b_a[0];

  // x^6
  d_b[1] = b_M[3] * b_a[1] + b_M[12] * b_a[0];

  // x^5*y
  d_b[2] = (b_M[3] * b_a[2] + b_M[12] * b_a[1]) + b_M[21] * b_a[0];

  // x^4*y^2
  d_b[3] = (b_M[3] * b_a[3] + b_M[12] * b_a[2]) + b_M[21] * b_a[1];

  // x^3*y^3
  d_b[4] = (b_M[3] * b_a[4] + b_M[12] * b_a[3]) + b_M[21] * b_a[2];

  // x^2*y^4
  d_b[5] = b_M[12] * b_a[4] + b_M[21] * b_a[3];

  // x*y^5
  d_b[6] = b_M[21] * b_a[4];

  // y^6
  d_b[7] = b_M[30] * b_a[0] + b_M[3] * b_a[5];

  // x^5
  d_b[8] = ((b_M[30] * b_a[1] + b_M[39] * b_a[0]) + b_M[3] * b_a[6]) + b_M[12] *
    b_a[5];

  // x^4*y
  d_b[9] = (((b_M[30] * b_a[2] + b_M[39] * b_a[1]) + b_M[3] * b_a[7]) + b_M[12] *
            b_a[6]) + b_M[21] * b_a[5];

  // x^3*y^2
  d_b[10] = (((b_M[30] * b_a[3] + b_M[39] * b_a[2]) + b_M[3] * b_a[8]) + b_M[12]
             * b_a[7]) + b_M[21] * b_a[6];

  // x^2*y^3
  d_b[11] = ((b_M[30] * b_a[4] + b_M[39] * b_a[3]) + b_M[12] * b_a[8]) + b_M[21]
    * b_a[7];

  // x*y^4
  d_b[12] = b_M[39] * b_a[4] + b_M[21] * b_a[8];

  // y^5
  d_b[13] = (b_M[48] * b_a[0] + b_M[30] * b_a[5]) + b_M[3] * b_a[9];

  // x^4
  d_b[14] = (((b_M[48] * b_a[1] + b_M[30] * b_a[6]) + b_M[39] * b_a[5]) + b_M[3]
             * b_a[10]) + b_M[12] * b_a[9];

  // x^3*y
  d_b[15] = ((((b_M[48] * b_a[2] + b_M[30] * b_a[7]) + b_M[39] * b_a[6]) + b_M[3]
              * b_a[11]) + b_M[12] * b_a[10]) + b_M[21] * b_a[9];

  // x^2*y^2
  d_b[16] = (((b_M[48] * b_a[3] + b_M[30] * b_a[8]) + b_M[39] * b_a[7]) + b_M[12]
             * b_a[11]) + b_M[21] * b_a[10];

  // x*y^3
  d_b[17] = (b_M[48] * b_a[4] + b_M[39] * b_a[8]) + b_M[21] * b_a[11];

  // y^4
  d_b[18] = (b_M[48] * b_a[5] + b_M[3] * b_a[12]) + b_M[30] * b_a[9];

  // x^3
  d_b[19] = (((b_M[48] * b_a[6] + b_M[3] * b_a[13]) + b_M[12] * b_a[12]) + b_M
             [30] * b_a[10]) + b_M[39] * b_a[9];

  // x^2*y
  d_b[20] = (((b_M[48] * b_a[7] + b_M[12] * b_a[13]) + b_M[21] * b_a[12]) + b_M
             [30] * b_a[11]) + b_M[39] * b_a[10];

  // x*y^2
  d_b[21] = (b_M[48] * b_a[8] + b_M[21] * b_a[13]) + b_M[39] * b_a[11];

  // y^3
  d_b[22] = (b_M[3] * b_a[14] + b_M[48] * b_a[9]) + b_M[30] * b_a[12];

  // x^2
  d_b[23] = ((b_M[12] * b_a[14] + b_M[48] * b_a[10]) + b_M[30] * b_a[13]) + b_M
    [39] * b_a[12];

  // x*y
  d_b[24] = (b_M[21] * b_a[14] + b_M[48] * b_a[11]) + b_M[39] * b_a[13];

  // y^2
  d_b[25] = b_M[30] * b_a[14] + b_M[48] * b_a[12];

  // x
  d_b[26] = b_M[39] * b_a[14] + b_M[48] * b_a[13];

  // y
  d_b[27] = b_M[48] * b_a[14];

  // 1
  b_a[0] = b_M[1] * b_M[5];

  // x^4
  b_a[1] = b_M[1] * b_M[14] + b_M[10] * b_M[5];

  // x^3*y
  b_a[2] = (b_M[1] * b_M[23] + b_M[10] * b_M[14]) + b_M[19] * b_M[5];

  // x^2*y^2
  b_a[3] = b_M[10] * b_M[23] + b_M[19] * b_M[14];

  // x*y^3
  b_a[4] = b_M[19] * b_M[23];

  // y^4
  b_a[5] = b_M[1] * b_M[32] + b_M[28] * b_M[5];

  // x^3
  b_a[6] = ((b_M[1] * b_M[41] + b_M[10] * b_M[32]) + b_M[28] * b_M[14]) + b_M[37]
    * b_M[5];

  // x^2*y
  b_a[7] = ((b_M[10] * b_M[41] + b_M[19] * b_M[32]) + b_M[28] * b_M[23]) + b_M
    [37] * b_M[14];

  // x*y^2
  b_a[8] = b_M[19] * b_M[41] + b_M[37] * b_M[23];

  // y^3
  b_a[9] = (b_M[1] * b_M[50] + b_M[46] * b_M[5]) + b_M[28] * b_M[32];

  // x^2
  b_a[10] = ((b_M[10] * b_M[50] + b_M[46] * b_M[14]) + b_M[28] * b_M[41]) + b_M
    [37] * b_M[32];

  // x*y
  b_a[11] = (b_M[19] * b_M[50] + b_M[46] * b_M[23]) + b_M[37] * b_M[41];

  // y^2
  b_a[12] = b_M[28] * b_M[50] + b_M[46] * b_M[32];

  // x
  b_a[13] = b_M[37] * b_M[50] + b_M[46] * b_M[41];

  // y
  b_a[14] = b_M[46] * b_M[50];

  // 1
  c_b[0] = b_M[4] * b_M[2];

  // x^4
  c_b[1] = b_M[4] * b_M[11] + b_M[13] * b_M[2];

  // x^3*y
  c_b[2] = (b_M[4] * b_M[20] + b_M[13] * b_M[11]) + b_M[22] * b_M[2];

  // x^2*y^2
  c_b[3] = b_M[13] * b_M[20] + b_M[22] * b_M[11];

  // x*y^3
  c_b[4] = b_M[22] * b_M[20];

  // y^4
  c_b[5] = b_M[4] * b_M[29] + b_M[31] * b_M[2];

  // x^3
  c_b[6] = ((b_M[4] * b_M[38] + b_M[13] * b_M[29]) + b_M[31] * b_M[11]) + b_M[40]
    * b_M[2];

  // x^2*y
  c_b[7] = ((b_M[13] * b_M[38] + b_M[22] * b_M[29]) + b_M[31] * b_M[20]) + b_M
    [40] * b_M[11];

  // x*y^2
  c_b[8] = b_M[22] * b_M[38] + b_M[40] * b_M[20];

  // y^3
  c_b[9] = (b_M[4] * b_M[47] + b_M[49] * b_M[2]) + b_M[31] * b_M[29];

  // x^2
  c_b[10] = ((b_M[13] * b_M[47] + b_M[49] * b_M[11]) + b_M[31] * b_M[38]) + b_M
    [40] * b_M[29];

  // x*y
  c_b[11] = (b_M[22] * b_M[47] + b_M[49] * b_M[20]) + b_M[40] * b_M[38];

  // y^2
  c_b[12] = b_M[31] * b_M[47] + b_M[49] * b_M[29];

  // x
  c_b[13] = b_M[40] * b_M[47] + b_M[49] * b_M[38];

  // y
  c_b[14] = b_M[49] * b_M[47];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  e_b[0] = b_M[6] * b_a[0];

  // x^6
  e_b[1] = b_M[6] * b_a[1] + b_M[15] * b_a[0];

  // x^5*y
  e_b[2] = (b_M[6] * b_a[2] + b_M[15] * b_a[1]) + b_M[24] * b_a[0];

  // x^4*y^2
  e_b[3] = (b_M[6] * b_a[3] + b_M[15] * b_a[2]) + b_M[24] * b_a[1];

  // x^3*y^3
  e_b[4] = (b_M[6] * b_a[4] + b_M[15] * b_a[3]) + b_M[24] * b_a[2];

  // x^2*y^4
  e_b[5] = b_M[15] * b_a[4] + b_M[24] * b_a[3];

  // x*y^5
  e_b[6] = b_M[24] * b_a[4];

  // y^6
  e_b[7] = b_M[33] * b_a[0] + b_M[6] * b_a[5];

  // x^5
  e_b[8] = ((b_M[33] * b_a[1] + b_M[42] * b_a[0]) + b_M[6] * b_a[6]) + b_M[15] *
    b_a[5];

  // x^4*y
  e_b[9] = (((b_M[33] * b_a[2] + b_M[42] * b_a[1]) + b_M[6] * b_a[7]) + b_M[15] *
            b_a[6]) + b_M[24] * b_a[5];

  // x^3*y^2
  e_b[10] = (((b_M[33] * b_a[3] + b_M[42] * b_a[2]) + b_M[6] * b_a[8]) + b_M[15]
             * b_a[7]) + b_M[24] * b_a[6];

  // x^2*y^3
  e_b[11] = ((b_M[33] * b_a[4] + b_M[42] * b_a[3]) + b_M[15] * b_a[8]) + b_M[24]
    * b_a[7];

  // x*y^4
  e_b[12] = b_M[42] * b_a[4] + b_M[24] * b_a[8];

  // y^5
  e_b[13] = (b_M[51] * b_a[0] + b_M[33] * b_a[5]) + b_M[6] * b_a[9];

  // x^4
  e_b[14] = (((b_M[51] * b_a[1] + b_M[33] * b_a[6]) + b_M[42] * b_a[5]) + b_M[6]
             * b_a[10]) + b_M[15] * b_a[9];

  // x^3*y
  e_b[15] = ((((b_M[51] * b_a[2] + b_M[33] * b_a[7]) + b_M[42] * b_a[6]) + b_M[6]
              * b_a[11]) + b_M[15] * b_a[10]) + b_M[24] * b_a[9];

  // x^2*y^2
  e_b[16] = (((b_M[51] * b_a[3] + b_M[33] * b_a[8]) + b_M[42] * b_a[7]) + b_M[15]
             * b_a[11]) + b_M[24] * b_a[10];

  // x*y^3
  e_b[17] = (b_M[51] * b_a[4] + b_M[42] * b_a[8]) + b_M[24] * b_a[11];

  // y^4
  e_b[18] = (b_M[51] * b_a[5] + b_M[6] * b_a[12]) + b_M[33] * b_a[9];

  // x^3
  e_b[19] = (((b_M[51] * b_a[6] + b_M[6] * b_a[13]) + b_M[15] * b_a[12]) + b_M
             [33] * b_a[10]) + b_M[42] * b_a[9];

  // x^2*y
  e_b[20] = (((b_M[51] * b_a[7] + b_M[15] * b_a[13]) + b_M[24] * b_a[12]) + b_M
             [33] * b_a[11]) + b_M[42] * b_a[10];

  // x*y^2
  e_b[21] = (b_M[51] * b_a[8] + b_M[24] * b_a[13]) + b_M[42] * b_a[11];

  // y^3
  e_b[22] = (b_M[6] * b_a[14] + b_M[51] * b_a[9]) + b_M[33] * b_a[12];

  // x^2
  e_b[23] = ((b_M[15] * b_a[14] + b_M[51] * b_a[10]) + b_M[33] * b_a[13]) + b_M
    [42] * b_a[12];

  // x*y
  e_b[24] = (b_M[24] * b_a[14] + b_M[51] * b_a[11]) + b_M[42] * b_a[13];

  // y^2
  e_b[25] = b_M[33] * b_a[14] + b_M[51] * b_a[12];

  // x
  e_b[26] = b_M[42] * b_a[14] + b_M[51] * b_a[13];

  // y
  e_b[27] = b_M[51] * b_a[14];

  // 1
  for (b_i = 0; b_i < 28; b_i++) {
    G[(b_i << 2) + 2] = (c_a[b_i] - d_b[b_i]) + e_b[b_i];
  }

  b_a[0] = F[5] * F[10];

  // x^4
  b_a[1] = F[5] * F[22] + F[17] * F[10];

  // x^3*y
  b_a[2] = (F[5] * F[34] + F[17] * F[22]) + F[29] * F[10];

  // x^2*y^2
  b_a[3] = F[17] * F[34] + F[29] * F[22];

  // x*y^3
  b_a[4] = F[29] * F[34];

  // y^4
  b_a[5] = F[5] * F[46] + F[41] * F[10];

  // x^3
  b_a[6] = ((F[5] * F[58] + F[17] * F[46]) + F[41] * F[22]) + F[53] * F[10];

  // x^2*y
  b_a[7] = ((F[17] * F[58] + F[29] * F[46]) + F[41] * F[34]) + F[53] * F[22];

  // x*y^2
  b_a[8] = F[29] * F[58] + F[53] * F[34];

  // y^3
  b_a[9] = (F[5] * F[70] + F[65] * F[10]) + F[41] * F[46];

  // x^2
  b_a[10] = ((F[17] * F[70] + F[65] * F[22]) + F[41] * F[58]) + F[53] * F[46];

  // x*y
  b_a[11] = (F[29] * F[70] + F[65] * F[34]) + F[53] * F[58];

  // y^2
  b_a[12] = F[41] * F[70] + F[65] * F[46];

  // x
  b_a[13] = F[53] * F[70] + F[65] * F[58];

  // y
  b_a[14] = F[65] * F[70];

  // 1
  c_b[0] = F[9] * F[6];

  // x^4
  c_b[1] = F[9] * F[18] + F[21] * F[6];

  // x^3*y
  c_b[2] = (F[9] * F[30] + F[21] * F[18]) + F[33] * F[6];

  // x^2*y^2
  c_b[3] = F[21] * F[30] + F[33] * F[18];

  // x*y^3
  c_b[4] = F[33] * F[30];

  // y^4
  c_b[5] = F[9] * F[42] + F[45] * F[6];

  // x^3
  c_b[6] = ((F[9] * F[54] + F[21] * F[42]) + F[45] * F[18]) + F[57] * F[6];

  // x^2*y
  c_b[7] = ((F[21] * F[54] + F[33] * F[42]) + F[45] * F[30]) + F[57] * F[18];

  // x*y^2
  c_b[8] = F[33] * F[54] + F[57] * F[30];

  // y^3
  c_b[9] = (F[9] * F[66] + F[69] * F[6]) + F[45] * F[42];

  // x^2
  c_b[10] = ((F[21] * F[66] + F[69] * F[18]) + F[45] * F[54]) + F[57] * F[42];

  // x*y
  c_b[11] = (F[33] * F[66] + F[69] * F[30]) + F[57] * F[54];

  // y^2
  c_b[12] = F[45] * F[66] + F[69] * F[42];

  // x
  c_b[13] = F[57] * F[66] + F[69] * F[54];

  // y
  c_b[14] = F[69] * F[66];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  c_a[0] = F[0] * b_a[0];

  // x^6
  c_a[1] = F[0] * b_a[1] + F[12] * b_a[0];

  // x^5*y
  c_a[2] = (F[0] * b_a[2] + F[12] * b_a[1]) + F[24] * b_a[0];

  // x^4*y^2
  c_a[3] = (F[0] * b_a[3] + F[12] * b_a[2]) + F[24] * b_a[1];

  // x^3*y^3
  c_a[4] = (F[0] * b_a[4] + F[12] * b_a[3]) + F[24] * b_a[2];

  // x^2*y^4
  c_a[5] = F[12] * b_a[4] + F[24] * b_a[3];

  // x*y^5
  c_a[6] = F[24] * b_a[4];

  // y^6
  c_a[7] = F[36] * b_a[0] + F[0] * b_a[5];

  // x^5
  c_a[8] = ((F[36] * b_a[1] + F[48] * b_a[0]) + F[0] * b_a[6]) + F[12] * b_a[5];

  // x^4*y
  c_a[9] = (((F[36] * b_a[2] + F[48] * b_a[1]) + F[0] * b_a[7]) + F[12] * b_a[6])
    + F[24] * b_a[5];

  // x^3*y^2
  c_a[10] = (((F[36] * b_a[3] + F[48] * b_a[2]) + F[0] * b_a[8]) + F[12] * b_a[7])
    + F[24] * b_a[6];

  // x^2*y^3
  c_a[11] = ((F[36] * b_a[4] + F[48] * b_a[3]) + F[12] * b_a[8]) + F[24] * b_a[7];

  // x*y^4
  c_a[12] = F[48] * b_a[4] + F[24] * b_a[8];

  // y^5
  c_a[13] = (F[60] * b_a[0] + F[36] * b_a[5]) + F[0] * b_a[9];

  // x^4
  c_a[14] = (((F[60] * b_a[1] + F[36] * b_a[6]) + F[48] * b_a[5]) + F[0] * b_a
             [10]) + F[12] * b_a[9];

  // x^3*y
  c_a[15] = ((((F[60] * b_a[2] + F[36] * b_a[7]) + F[48] * b_a[6]) + F[0] * b_a
              [11]) + F[12] * b_a[10]) + F[24] * b_a[9];

  // x^2*y^2
  c_a[16] = (((F[60] * b_a[3] + F[36] * b_a[8]) + F[48] * b_a[7]) + F[12] * b_a
             [11]) + F[24] * b_a[10];

  // x*y^3
  c_a[17] = (F[60] * b_a[4] + F[48] * b_a[8]) + F[24] * b_a[11];

  // y^4
  c_a[18] = (F[60] * b_a[5] + F[0] * b_a[12]) + F[36] * b_a[9];

  // x^3
  c_a[19] = (((F[60] * b_a[6] + F[0] * b_a[13]) + F[12] * b_a[12]) + F[36] *
             b_a[10]) + F[48] * b_a[9];

  // x^2*y
  c_a[20] = (((F[60] * b_a[7] + F[12] * b_a[13]) + F[24] * b_a[12]) + F[36] *
             b_a[11]) + F[48] * b_a[10];

  // x*y^2
  c_a[21] = (F[60] * b_a[8] + F[24] * b_a[13]) + F[48] * b_a[11];

  // y^3
  c_a[22] = (F[0] * b_a[14] + F[60] * b_a[9]) + F[36] * b_a[12];

  // x^2
  c_a[23] = ((F[12] * b_a[14] + F[60] * b_a[10]) + F[36] * b_a[13]) + F[48] *
    b_a[12];

  // x*y
  c_a[24] = (F[24] * b_a[14] + F[60] * b_a[11]) + F[48] * b_a[13];

  // y^2
  c_a[25] = F[36] * b_a[14] + F[60] * b_a[12];

  // x
  c_a[26] = F[48] * b_a[14] + F[60] * b_a[13];

  // y
  c_a[27] = F[60] * b_a[14];

  // 1
  for (b_i = 0; b_i < 6; b_i++) {
    M_tmp = b_i << 2;
    M[M_tmp] = F[12 * b_i + 1];
    M[M_tmp + 2] = F[12 * b_i + 9];
    M[M_tmp + 1] = F[12 * b_i + 2];
    M[M_tmp + 3] = F[12 * b_i + 10];
  }

  b_a[0] = M[0] * M[3];

  // x^4
  b_a[1] = M[0] * M[7] + M[4] * M[3];

  // x^3*y
  b_a[2] = (M[0] * M[11] + M[4] * M[7]) + M[8] * M[3];

  // x^2*y^2
  b_a[3] = M[4] * M[11] + M[8] * M[7];

  // x*y^3
  b_a[4] = M[8] * M[11];

  // y^4
  b_a[5] = M[0] * M[15] + M[12] * M[3];

  // x^3
  b_a[6] = ((M[0] * M[19] + M[4] * M[15]) + M[12] * M[7]) + M[16] * M[3];

  // x^2*y
  b_a[7] = ((M[4] * M[19] + M[8] * M[15]) + M[12] * M[11]) + M[16] * M[7];

  // x*y^2
  b_a[8] = M[8] * M[19] + M[16] * M[11];

  // y^3
  b_a[9] = (M[0] * M[23] + M[20] * M[3]) + M[12] * M[15];

  // x^2
  b_a[10] = ((M[4] * M[23] + M[20] * M[7]) + M[12] * M[19]) + M[16] * M[15];

  // x*y
  b_a[11] = (M[8] * M[23] + M[20] * M[11]) + M[16] * M[19];

  // y^2
  b_a[12] = M[12] * M[23] + M[20] * M[15];

  // x
  b_a[13] = M[16] * M[23] + M[20] * M[19];

  // y
  b_a[14] = M[20] * M[23];

  // 1
  c_b[0] = M[2] * M[1];

  // x^4
  c_b[1] = M[2] * M[5] + M[6] * M[1];

  // x^3*y
  c_b[2] = (M[2] * M[9] + M[6] * M[5]) + M[10] * M[1];

  // x^2*y^2
  c_b[3] = M[6] * M[9] + M[10] * M[5];

  // x*y^3
  c_b[4] = M[10] * M[9];

  // y^4
  c_b[5] = M[2] * M[13] + M[14] * M[1];

  // x^3
  c_b[6] = ((M[2] * M[17] + M[6] * M[13]) + M[14] * M[5]) + M[18] * M[1];

  // x^2*y
  c_b[7] = ((M[6] * M[17] + M[10] * M[13]) + M[14] * M[9]) + M[18] * M[5];

  // x*y^2
  c_b[8] = M[10] * M[17] + M[18] * M[9];

  // y^3
  c_b[9] = (M[2] * M[21] + M[22] * M[1]) + M[14] * M[13];

  // x^2
  c_b[10] = ((M[6] * M[21] + M[22] * M[5]) + M[14] * M[17]) + M[18] * M[13];

  // x*y
  c_b[11] = (M[10] * M[21] + M[22] * M[9]) + M[18] * M[17];

  // y^2
  c_b[12] = M[14] * M[21] + M[22] * M[13];

  // x
  c_b[13] = M[18] * M[21] + M[22] * M[17];

  // y
  c_b[14] = M[22] * M[21];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  d_b[0] = F[4] * b_a[0];

  // x^6
  d_b[1] = F[4] * b_a[1] + F[16] * b_a[0];

  // x^5*y
  d_b[2] = (F[4] * b_a[2] + F[16] * b_a[1]) + F[28] * b_a[0];

  // x^4*y^2
  d_b[3] = (F[4] * b_a[3] + F[16] * b_a[2]) + F[28] * b_a[1];

  // x^3*y^3
  d_b[4] = (F[4] * b_a[4] + F[16] * b_a[3]) + F[28] * b_a[2];

  // x^2*y^4
  d_b[5] = F[16] * b_a[4] + F[28] * b_a[3];

  // x*y^5
  d_b[6] = F[28] * b_a[4];

  // y^6
  d_b[7] = F[40] * b_a[0] + F[4] * b_a[5];

  // x^5
  d_b[8] = ((F[40] * b_a[1] + F[52] * b_a[0]) + F[4] * b_a[6]) + F[16] * b_a[5];

  // x^4*y
  d_b[9] = (((F[40] * b_a[2] + F[52] * b_a[1]) + F[4] * b_a[7]) + F[16] * b_a[6])
    + F[28] * b_a[5];

  // x^3*y^2
  d_b[10] = (((F[40] * b_a[3] + F[52] * b_a[2]) + F[4] * b_a[8]) + F[16] * b_a[7])
    + F[28] * b_a[6];

  // x^2*y^3
  d_b[11] = ((F[40] * b_a[4] + F[52] * b_a[3]) + F[16] * b_a[8]) + F[28] * b_a[7];

  // x*y^4
  d_b[12] = F[52] * b_a[4] + F[28] * b_a[8];

  // y^5
  d_b[13] = (F[64] * b_a[0] + F[40] * b_a[5]) + F[4] * b_a[9];

  // x^4
  d_b[14] = (((F[64] * b_a[1] + F[40] * b_a[6]) + F[52] * b_a[5]) + F[4] * b_a
             [10]) + F[16] * b_a[9];

  // x^3*y
  d_b[15] = ((((F[64] * b_a[2] + F[40] * b_a[7]) + F[52] * b_a[6]) + F[4] * b_a
              [11]) + F[16] * b_a[10]) + F[28] * b_a[9];

  // x^2*y^2
  d_b[16] = (((F[64] * b_a[3] + F[40] * b_a[8]) + F[52] * b_a[7]) + F[16] * b_a
             [11]) + F[28] * b_a[10];

  // x*y^3
  d_b[17] = (F[64] * b_a[4] + F[52] * b_a[8]) + F[28] * b_a[11];

  // y^4
  d_b[18] = (F[64] * b_a[5] + F[4] * b_a[12]) + F[40] * b_a[9];

  // x^3
  d_b[19] = (((F[64] * b_a[6] + F[4] * b_a[13]) + F[16] * b_a[12]) + F[40] *
             b_a[10]) + F[52] * b_a[9];

  // x^2*y
  d_b[20] = (((F[64] * b_a[7] + F[16] * b_a[13]) + F[28] * b_a[12]) + F[40] *
             b_a[11]) + F[52] * b_a[10];

  // x*y^2
  d_b[21] = (F[64] * b_a[8] + F[28] * b_a[13]) + F[52] * b_a[11];

  // y^3
  d_b[22] = (F[4] * b_a[14] + F[64] * b_a[9]) + F[40] * b_a[12];

  // x^2
  d_b[23] = ((F[16] * b_a[14] + F[64] * b_a[10]) + F[40] * b_a[13]) + F[52] *
    b_a[12];

  // x*y
  d_b[24] = (F[28] * b_a[14] + F[64] * b_a[11]) + F[52] * b_a[13];

  // y^2
  d_b[25] = F[40] * b_a[14] + F[64] * b_a[12];

  // x
  d_b[26] = F[52] * b_a[14] + F[64] * b_a[13];

  // y
  d_b[27] = F[64] * b_a[14];

  // 1
  b_a[0] = F[1] * F[6];

  // x^4
  b_a[1] = F[1] * F[18] + F[13] * F[6];

  // x^3*y
  b_a[2] = (F[1] * F[30] + F[13] * F[18]) + F[25] * F[6];

  // x^2*y^2
  b_a[3] = F[13] * F[30] + F[25] * F[18];

  // x*y^3
  b_a[4] = F[25] * F[30];

  // y^4
  b_a[5] = F[1] * F[42] + F[37] * F[6];

  // x^3
  b_a[6] = ((F[1] * F[54] + F[13] * F[42]) + F[37] * F[18]) + F[49] * F[6];

  // x^2*y
  b_a[7] = ((F[13] * F[54] + F[25] * F[42]) + F[37] * F[30]) + F[49] * F[18];

  // x*y^2
  b_a[8] = F[25] * F[54] + F[49] * F[30];

  // y^3
  b_a[9] = (F[1] * F[66] + F[61] * F[6]) + F[37] * F[42];

  // x^2
  b_a[10] = ((F[13] * F[66] + F[61] * F[18]) + F[37] * F[54]) + F[49] * F[42];

  // x*y
  b_a[11] = (F[25] * F[66] + F[61] * F[30]) + F[49] * F[54];

  // y^2
  b_a[12] = F[37] * F[66] + F[61] * F[42];

  // x
  b_a[13] = F[49] * F[66] + F[61] * F[54];

  // y
  b_a[14] = F[61] * F[66];

  // 1
  c_b[0] = F[5] * F[2];

  // x^4
  c_b[1] = F[5] * F[14] + F[17] * F[2];

  // x^3*y
  c_b[2] = (F[5] * F[26] + F[17] * F[14]) + F[29] * F[2];

  // x^2*y^2
  c_b[3] = F[17] * F[26] + F[29] * F[14];

  // x*y^3
  c_b[4] = F[29] * F[26];

  // y^4
  c_b[5] = F[5] * F[38] + F[41] * F[2];

  // x^3
  c_b[6] = ((F[5] * F[50] + F[17] * F[38]) + F[41] * F[14]) + F[53] * F[2];

  // x^2*y
  c_b[7] = ((F[17] * F[50] + F[29] * F[38]) + F[41] * F[26]) + F[53] * F[14];

  // x*y^2
  c_b[8] = F[29] * F[50] + F[53] * F[26];

  // y^3
  c_b[9] = (F[5] * F[62] + F[65] * F[2]) + F[41] * F[38];

  // x^2
  c_b[10] = ((F[17] * F[62] + F[65] * F[14]) + F[41] * F[50]) + F[53] * F[38];

  // x*y
  c_b[11] = (F[29] * F[62] + F[65] * F[26]) + F[53] * F[50];

  // y^2
  c_b[12] = F[41] * F[62] + F[65] * F[38];

  // x
  c_b[13] = F[53] * F[62] + F[65] * F[50];

  // y
  c_b[14] = F[65] * F[62];

  // 1
  for (b_i = 0; b_i < 15; b_i++) {
    b_a[b_i] -= c_b[b_i];
  }

  e_b[0] = F[8] * b_a[0];

  // x^6
  e_b[1] = F[8] * b_a[1] + F[20] * b_a[0];

  // x^5*y
  e_b[2] = (F[8] * b_a[2] + F[20] * b_a[1]) + F[32] * b_a[0];

  // x^4*y^2
  e_b[3] = (F[8] * b_a[3] + F[20] * b_a[2]) + F[32] * b_a[1];

  // x^3*y^3
  e_b[4] = (F[8] * b_a[4] + F[20] * b_a[3]) + F[32] * b_a[2];

  // x^2*y^4
  e_b[5] = F[20] * b_a[4] + F[32] * b_a[3];

  // x*y^5
  e_b[6] = F[32] * b_a[4];

  // y^6
  e_b[7] = F[44] * b_a[0] + F[8] * b_a[5];

  // x^5
  e_b[8] = ((F[44] * b_a[1] + F[56] * b_a[0]) + F[8] * b_a[6]) + F[20] * b_a[5];

  // x^4*y
  e_b[9] = (((F[44] * b_a[2] + F[56] * b_a[1]) + F[8] * b_a[7]) + F[20] * b_a[6])
    + F[32] * b_a[5];

  // x^3*y^2
  e_b[10] = (((F[44] * b_a[3] + F[56] * b_a[2]) + F[8] * b_a[8]) + F[20] * b_a[7])
    + F[32] * b_a[6];

  // x^2*y^3
  e_b[11] = ((F[44] * b_a[4] + F[56] * b_a[3]) + F[20] * b_a[8]) + F[32] * b_a[7];

  // x*y^4
  e_b[12] = F[56] * b_a[4] + F[32] * b_a[8];

  // y^5
  e_b[13] = (F[68] * b_a[0] + F[44] * b_a[5]) + F[8] * b_a[9];

  // x^4
  e_b[14] = (((F[68] * b_a[1] + F[44] * b_a[6]) + F[56] * b_a[5]) + F[8] * b_a
             [10]) + F[20] * b_a[9];

  // x^3*y
  e_b[15] = ((((F[68] * b_a[2] + F[44] * b_a[7]) + F[56] * b_a[6]) + F[8] * b_a
              [11]) + F[20] * b_a[10]) + F[32] * b_a[9];

  // x^2*y^2
  e_b[16] = (((F[68] * b_a[3] + F[44] * b_a[8]) + F[56] * b_a[7]) + F[20] * b_a
             [11]) + F[32] * b_a[10];

  // x*y^3
  e_b[17] = (F[68] * b_a[4] + F[56] * b_a[8]) + F[32] * b_a[11];

  // y^4
  e_b[18] = (F[68] * b_a[5] + F[8] * b_a[12]) + F[44] * b_a[9];

  // x^3
  e_b[19] = (((F[68] * b_a[6] + F[8] * b_a[13]) + F[20] * b_a[12]) + F[44] *
             b_a[10]) + F[56] * b_a[9];

  // x^2*y
  e_b[20] = (((F[68] * b_a[7] + F[20] * b_a[13]) + F[32] * b_a[12]) + F[44] *
             b_a[11]) + F[56] * b_a[10];

  // x*y^2
  e_b[21] = (F[68] * b_a[8] + F[32] * b_a[13]) + F[56] * b_a[11];

  // y^3
  e_b[22] = (F[8] * b_a[14] + F[68] * b_a[9]) + F[44] * b_a[12];

  // x^2
  e_b[23] = ((F[20] * b_a[14] + F[68] * b_a[10]) + F[44] * b_a[13]) + F[56] *
    b_a[12];

  // x*y
  e_b[24] = (F[32] * b_a[14] + F[68] * b_a[11]) + F[56] * b_a[13];

  // y^2
  e_b[25] = F[44] * b_a[14] + F[68] * b_a[12];

  // x
  e_b[26] = F[56] * b_a[14] + F[68] * b_a[13];

  // y
  e_b[27] = F[68] * b_a[14];

  // 1
  for (b_i = 0; b_i < 28; b_i++) {
    G[(b_i << 2) + 3] = (c_a[b_i] - d_b[b_i]) + e_b[b_i];
  }

  // 4x28
  std::memset(&G20[0], 0, 900U * sizeof(float));

  // multiply by qx^2
  // multiply by qxqy
  // multiply by qx
  // multiply by qy
  // multiply by 1
  for (i = 0; i < 4; i++) {
    G20[i] = G[i];
    qy_im = G[i + 4];
    G20[i + 20] = qy_im;
    f = G[i + 8];
    G20[i + 40] = f;
    p4_tmp = G[i + 12];
    G20[i + 60] = p4_tmp;
    b_p4_tmp = G[i + 16];
    G20[i + 80] = b_p4_tmp;
    c_p4_tmp = G[i + 20];
    G20[i + 100] = c_p4_tmp;
    a = G[i + 24];
    G20[i + 120] = a;
    G20_tmp = G[i + 28];
    G20[i + 180] = G20_tmp;
    qy_re = G[i + 32];
    G20[i + 200] = qy_re;
    b_G20_tmp = G[i + 36];
    G20[i + 220] = b_G20_tmp;
    c_G20_tmp = G[i + 40];
    G20[i + 240] = c_G20_tmp;
    d_G20_tmp = G[i + 44];
    G20[i + 260] = d_G20_tmp;
    e_G20_tmp = G[i + 48];
    G20[i + 280] = e_G20_tmp;
    f_G20_tmp = G[i + 52];
    G20[i + 340] = f_G20_tmp;
    g_G20_tmp = G[i + 56];
    G20[i + 360] = g_G20_tmp;
    h_G20_tmp = G[i + 60];
    G20[i + 380] = h_G20_tmp;
    i_G20_tmp = G[i + 64];
    G20[i + 400] = i_G20_tmp;
    j_G20_tmp = G[i + 68];
    G20[i + 420] = j_G20_tmp;
    k_G20_tmp = G[i + 72];
    G20[i + 480] = k_G20_tmp;
    l_G20_tmp = G[i + 76];
    G20[i + 500] = l_G20_tmp;
    m_G20_tmp = G[i + 80];
    G20[i + 520] = m_G20_tmp;
    n_G20_tmp = G[i + 84];
    G20[i + 540] = n_G20_tmp;
    o_G20_tmp = G[i + 88];
    G20[i + 600] = o_G20_tmp;
    p_G20_tmp = G[i + 92];
    G20[i + 620] = p_G20_tmp;
    q_G20_tmp = G[i + 96];
    G20[i + 640] = q_G20_tmp;
    r_G20_tmp = G[i + 100];
    G20[i + 700] = r_G20_tmp;
    s_G20_tmp = G[i + 104];
    G20[i + 720] = s_G20_tmp;
    t_G20_tmp = G[i + 108];
    G20[i + 780] = t_G20_tmp;
    G20[i + 24] = G[i];
    G20[i + 44] = qy_im;
    G20[i + 64] = f;
    G20[i + 84] = p4_tmp;
    G20[i + 104] = b_p4_tmp;
    G20[i + 124] = c_p4_tmp;
    G20[i + 144] = a;
    G20[i + 204] = G20_tmp;
    G20[i + 224] = qy_re;
    G20[i + 244] = b_G20_tmp;
    G20[i + 264] = c_G20_tmp;
    G20[i + 284] = d_G20_tmp;
    G20[i + 304] = e_G20_tmp;
    G20[i + 364] = f_G20_tmp;
    G20[i + 384] = g_G20_tmp;
    G20[i + 404] = h_G20_tmp;
    G20[i + 424] = i_G20_tmp;
    G20[i + 444] = j_G20_tmp;
    G20[i + 504] = k_G20_tmp;
    G20[i + 524] = l_G20_tmp;
    G20[i + 544] = m_G20_tmp;
    G20[i + 564] = n_G20_tmp;
    G20[i + 624] = o_G20_tmp;
    G20[i + 644] = p_G20_tmp;
    G20[i + 664] = q_G20_tmp;
    G20[i + 724] = r_G20_tmp;
    G20[i + 744] = s_G20_tmp;
    G20[i + 804] = t_G20_tmp;
    G20[i + 188] = G[i];
    G20[i + 208] = qy_im;
    G20[i + 228] = f;
    G20[i + 248] = p4_tmp;
    G20[i + 268] = b_p4_tmp;
    G20[i + 288] = c_p4_tmp;
    G20[i + 308] = a;
    G20[i + 348] = G20_tmp;
    G20[i + 368] = qy_re;
    G20[i + 388] = b_G20_tmp;
    G20[i + 408] = c_G20_tmp;
    G20[i + 428] = d_G20_tmp;
    G20[i + 448] = e_G20_tmp;
    G20[i + 488] = f_G20_tmp;
    G20[i + 508] = g_G20_tmp;
    G20[i + 528] = h_G20_tmp;
    G20[i + 548] = i_G20_tmp;
    G20[i + 568] = j_G20_tmp;
    G20[i + 608] = k_G20_tmp;
    G20[i + 628] = l_G20_tmp;
    G20[i + 648] = m_G20_tmp;
    G20[i + 668] = n_G20_tmp;
    G20[i + 708] = o_G20_tmp;
    G20[i + 728] = p_G20_tmp;
    G20[i + 748] = q_G20_tmp;
    G20[i + 788] = r_G20_tmp;
    G20[i + 808] = s_G20_tmp;
    G20[i + 848] = t_G20_tmp;
    G20[i + 212] = G[i];
    G20[i + 232] = qy_im;
    G20[i + 252] = f;
    G20[i + 272] = p4_tmp;
    G20[i + 292] = b_p4_tmp;
    G20[i + 312] = c_p4_tmp;
    G20[i + 332] = a;
    G20[i + 372] = G20_tmp;
    G20[i + 392] = qy_re;
    G20[i + 412] = b_G20_tmp;
    G20[i + 432] = c_G20_tmp;
    G20[i + 452] = d_G20_tmp;
    G20[i + 472] = e_G20_tmp;
    G20[i + 512] = f_G20_tmp;
    G20[i + 532] = g_G20_tmp;
    G20[i + 552] = h_G20_tmp;
    G20[i + 572] = i_G20_tmp;
    G20[i + 592] = j_G20_tmp;
    G20[i + 632] = k_G20_tmp;
    G20[i + 652] = l_G20_tmp;
    G20[i + 672] = m_G20_tmp;
    G20[i + 692] = n_G20_tmp;
    G20[i + 732] = o_G20_tmp;
    G20[i + 752] = p_G20_tmp;
    G20[i + 772] = q_G20_tmp;
    G20[i + 812] = r_G20_tmp;
    G20[i + 832] = s_G20_tmp;
    G20[i + 872] = t_G20_tmp;
    for (b_i = 0; b_i < 28; b_i++) {
      G20[(i + 20 * (b_i + 17)) + 16] = G[i + (b_i << 2)];
    }
  }

  // 20x45
  for (b_i = 0; b_i < 10; b_i++) {
    std::memcpy(&C[b_i * 20], &G20[b_i * 20 + 700], 20U * sizeof(float));
  }

  for (b_i = 0; b_i < 2; b_i++) {
    std::memcpy(&b_G20[b_i * 20], &G20[b_i * 20], 20U * sizeof(float));
  }

  for (b_i = 0; b_i < 4; b_i++) {
    std::memcpy(&b_G20[b_i * 20 + 40], &G20[b_i * 20 + 180], 20U * sizeof(float));
  }

  for (b_i = 0; b_i < 5; b_i++) {
    std::memcpy(&b_G20[b_i * 20 + 120], &G20[b_i * 20 + 340], 20U * sizeof(float));
  }

  for (b_i = 0; b_i < 4; b_i++) {
    std::memcpy(&b_G20[b_i * 20 + 220], &G20[b_i * 20 + 480], 20U * sizeof(float));
  }

  for (b_i = 0; b_i < 5; b_i++) {
    std::memcpy(&b_G20[b_i * 20 + 300], &G20[b_i * 20 + 600], 20U * sizeof(float));
  }

  mldivide(b_G20, C);

  // this function creates a matrix for multiplication by x in the monomial
  // basis B
  // monomial basis B = {x^3, ...., 1} -- monomials up to the 3d degree, #B = 10 
  // x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1
  std::memset(&c_M[0], 0, 100U * sizeof(float));
  for (i = 0; i < 4; i++) {
    for (b_i = 0; b_i < 10; b_i++) {
      c_M[b_i + 10 * i] = -C[(i + 20 * b_i) + 15];
    }
  }

  c_M[40] = 1.0F;
  c_M[51] = 1.0F;
  c_M[62] = 1.0F;
  c_M[74] = 1.0F;
  c_M[85] = 1.0F;
  c_M[97] = 1.0F;
  for (b_i = 0; b_i < 10; b_i++) {
    for (i1 = 0; i1 < 10; i1++) {
      d_M[i1 + 10 * b_i] = c_M[b_i + 10 * i1];
    }
  }

  eig(d_M, W, D);
  *solution_num = 0.0F;
  f_sol_size[0] = 1;
  f_sol_size[1] = 0;
  R_sol_size[0] = 3;
  R_sol_size[1] = 3;
  R_sol_size[2] = 0;
  T_sol_size[0] = 3;
  T_sol_size[1] = 0;
  for (i = 0; i < 10; i++) {
    M_tmp = 10 * i + 8;
    b_M_tmp = 10 * i + 9;
    if (W[b_M_tmp].im == 0.0F) {
      if (W[M_tmp].im == 0.0F) {
        qy_re = W[M_tmp].re / W[b_M_tmp].re;
        qy_im = 0.0F;
      } else if (W[M_tmp].re == 0.0F) {
        qy_re = 0.0F;
        qy_im = W[M_tmp].im / W[b_M_tmp].re;
      } else {
        qy_re = W[M_tmp].re / W[b_M_tmp].re;
        qy_im = W[M_tmp].im / W[b_M_tmp].re;
      }
    } else if (W[b_M_tmp].re == 0.0F) {
      if (W[M_tmp].re == 0.0F) {
        qy_re = W[M_tmp].im / W[b_M_tmp].im;
        qy_im = 0.0F;
      } else if (W[M_tmp].im == 0.0F) {
        qy_re = 0.0F;
        qy_im = -(W[M_tmp].re / W[b_M_tmp].im);
      } else {
        qy_re = W[M_tmp].im / W[b_M_tmp].im;
        qy_im = -(W[M_tmp].re / W[b_M_tmp].im);
      }
    } else {
      b_p4_tmp = std::abs(W[b_M_tmp].re);
      qy_im = std::abs(W[b_M_tmp].im);
      if (b_p4_tmp > qy_im) {
        G20_tmp = W[b_M_tmp].im / W[b_M_tmp].re;
        qy_im = W[b_M_tmp].re + G20_tmp * W[b_M_tmp].im;
        qy_re = (W[M_tmp].re + G20_tmp * W[M_tmp].im) / qy_im;
        qy_im = (W[M_tmp].im - G20_tmp * W[M_tmp].re) / qy_im;
      } else if (qy_im == b_p4_tmp) {
        if (W[b_M_tmp].re > 0.0F) {
          qy_im = 0.5F;
        } else {
          qy_im = -0.5F;
        }

        if (W[b_M_tmp].im > 0.0F) {
          f = 0.5F;
        } else {
          f = -0.5F;
        }

        qy_re = (W[M_tmp].re * qy_im + W[M_tmp].im * f) / b_p4_tmp;
        qy_im = (W[M_tmp].im * qy_im - W[M_tmp].re * f) / b_p4_tmp;
      } else {
        G20_tmp = W[b_M_tmp].re / W[b_M_tmp].im;
        qy_im = W[b_M_tmp].im + G20_tmp * W[b_M_tmp].re;
        qy_re = (G20_tmp * W[M_tmp].re + W[M_tmp].im) / qy_im;
        qy_im = (G20_tmp * W[M_tmp].im - W[M_tmp].re) / qy_im;
      }
    }

    b_i = i + 10 * i;
    if ((!(std::abs(D[b_i].im) > e)) && (!(std::abs(qy_im) > e))) {
      find_f(F, D[b_i].re, qy_re, e, fc_set_data, fc_set_size, fs_set_data,
             fs_set_size, &qy_im);
      i1 = static_cast<int>(qy_im);
      if (0 <= i1 - 1) {
        f = qy_re * qy_re;
        p4_tmp = D[b_i].re * D[b_i].re;
        G20_tmp = 1.0F / ((p4_tmp + 1.0F) + f);
        b_p4_tmp = 2.0F * G20_tmp;
        c_p4_tmp = b_p4_tmp * D[b_i].re;
        a = c_p4_tmp * qy_re;
        qy_im = -2.0F * G20_tmp;
        R_xy[0] = 1.0F - b_p4_tmp * f;
        R_xy[3] = a;
        R_xy[6] = b_p4_tmp * qy_re;
        R_xy[1] = a;
        R_xy[4] = 1.0F - b_p4_tmp * p4_tmp;
        R_xy[7] = qy_im * D[b_i].re;
        R_xy[2] = qy_im * qy_re;
        R_xy[5] = c_p4_tmp;
        R_xy[8] = 1.0F - b_p4_tmp * (p4_tmp + f);
        b_x = x[1];
      }

      for (j = 0; j < i1; j++) {
        // li = find_lamda(R, fc, fs, Xi, Xj, xi, xj) -- signature
        fc_set_tmp = static_cast<int>((static_cast<float>(j) + 1.0F)) - 1;
        qy_im = fc_set_data[fc_set_tmp] * R_xy[0] - fs_set_data[fc_set_tmp] *
          R_xy[1];
        f = fc_set_data[fc_set_tmp] * R_xy[3] - fs_set_data[fc_set_tmp] * R_xy[4];
        b_p4_tmp = fc_set_data[fc_set_tmp] * R_xy[6] - fs_set_data[fc_set_tmp] *
          R_xy[7];
        p4_tmp = ((-(qy_im - b_x * R_xy[2]) * T_tmp + -(f - b_x * R_xy[5]) *
                   b_T_tmp) + -(b_p4_tmp - b_x * R_xy[8]) * c_T_tmp) /
          F_row_tmp_tmp;

        // T = find_translation(R, fc, fs, li, Xi, xi, yi)
        T[0] = p4_tmp * x[0] - ((qy_im * X[0] + f * X[1]) + b_p4_tmp * X[2]);
        T[1] = p4_tmp * y[0] - (((fs_set_data[fc_set_tmp] * R_xy[0] +
          fc_set_data[fc_set_tmp] * R_xy[1]) * X[0] + (fs_set_data[fc_set_tmp] *
          R_xy[3] + fc_set_data[fc_set_tmp] * R_xy[4]) * X[1]) +
          (fs_set_data[fc_set_tmp] * R_xy[6] + fc_set_data[fc_set_tmp] * R_xy[7])
          * X[2]);
        T[2] = p4_tmp - ((R_xy[2] * X[0] + R_xy[5] * X[1]) + R_xy[8] * X[2]);
        f = rt_hypotf_snf(fs_set_data[fc_set_tmp], fc_set_data[fc_set_tmp]);
        fc_set[0] = fc_set_data[fc_set_tmp];
        fc_set[3] = -fs_set_data[fc_set_tmp];
        fc_set[6] = 0.0F;
        fc_set[1] = fs_set_data[fc_set_tmp];
        fc_set[4] = fc_set_data[fc_set_tmp];
        fc_set[7] = 0.0F;
        fc_set[2] = 0.0F;
        fc_set[5] = 0.0F;
        fc_set[8] = 1.0F;
        for (b_i = 0; b_i < 3; b_i++) {
          G20_tmp = fc_set[b_i + 3];
          M_tmp = static_cast<int>(fc_set[b_i + 6]);
          for (b_M_tmp = 0; b_M_tmp < 3; b_M_tmp++) {
            R_curr[b_i + 3 * b_M_tmp] = (fc_set[b_i] * R_xy[3 * b_M_tmp] +
              G20_tmp * R_xy[3 * b_M_tmp + 1]) + static_cast<float>(M_tmp) *
              R_xy[3 * b_M_tmp + 2];
          }
        }

        for (b_i = 0; b_i < 3; b_i++) {
          b_fc_set[3 * b_i] = R_curr[3 * b_i];
          M_tmp = 3 * b_i + 1;
          b_fc_set[M_tmp] = R_curr[M_tmp];
          M_tmp = 3 * b_i + 2;
          b_fc_set[M_tmp] = R_curr[M_tmp];
          b_fc_set[b_i + 9] = T[b_i];
          b_X[b_i] = X[b_i + 9];
        }

        for (b_i = 0; b_i < 3; b_i++) {
          p4[b_i] = ((b_fc_set[b_i] * b_X[0] + b_fc_set[b_i + 3] * b_X[1]) +
                     b_fc_set[b_i + 6] * b_X[2]) + b_fc_set[b_i + 9];
        }

        if (std::abs(p4[1] / p4[2] - y[3]) < 0.01F * f) {
          (*solution_num)++;
          qy_im = fc_set_data[fc_set_tmp] / f;
          fc_set[0] = qy_im;
          fc_set[3] = -fs_set_data[fc_set_tmp] / f;
          fc_set[6] = 0.0F;
          fc_set[1] = fs_set_data[fc_set_tmp] / f;
          fc_set[4] = qy_im;
          fc_set[7] = 0.0F;
          fc_set[2] = 0.0F;
          fc_set[5] = 0.0F;
          fc_set[8] = 1.0F;
          for (b_i = 0; b_i < 3; b_i++) {
            G20_tmp = fc_set[b_i + 3];
            M_tmp = static_cast<int>(fc_set[b_i + 6]);
            for (b_M_tmp = 0; b_M_tmp < 3; b_M_tmp++) {
              R_curr[b_i + 3 * b_M_tmp] = (fc_set[b_i] * R_xy[3 * b_M_tmp] +
                G20_tmp * R_xy[3 * b_M_tmp + 1]) + static_cast<float>(M_tmp) *
                R_xy[3 * b_M_tmp + 2];
            }
          }

          b_i = f_sol_size[1];
          f_sol_size[1]++;
          f_sol_data[b_i] = f;
          if (*solution_num == 1.0F) {
            R_sol_size[0] = 3;
            R_sol_size[1] = 3;
            R_sol_size[2] = 1;
            for (b_i = 0; b_i < 9; b_i++) {
              R_sol_data[b_i] = R_curr[b_i];
            }
          } else {
            cat(R_sol_data, R_sol_size, R_curr, tmp_data, tmp_size);
            R_sol_size[0] = 3;
            R_sol_size[1] = 3;
            R_sol_size[2] = tmp_size[2];
            M_tmp = tmp_size[0] * tmp_size[1] * tmp_size[2];
            if (0 <= M_tmp - 1) {
              std::memcpy(&R_sol_data[0], &tmp_data[0], M_tmp * sizeof(float));
            }
          }

          M_tmp = T_sol_size[1];
          b_M_tmp = T_sol_size[1] + 1;
          for (b_i = 0; b_i < M_tmp; b_i++) {
            b_T_sol_data[3 * b_i] = T_sol_data[3 * b_i];
            fc_set_tmp = 3 * b_i + 1;
            b_T_sol_data[fc_set_tmp] = T_sol_data[fc_set_tmp];
            fc_set_tmp = 3 * b_i + 2;
            b_T_sol_data[fc_set_tmp] = T_sol_data[fc_set_tmp];
          }

          b_T_sol_data[3 * T_sol_size[1]] = T[0];
          b_T_sol_data[3 * T_sol_size[1] + 1] = T[1];
          b_T_sol_data[3 * T_sol_size[1] + 2] = T[2];
          T_sol_size[0] = 3;
          T_sol_size[1] = b_M_tmp;
          M_tmp = 3 * b_M_tmp;
          if (0 <= M_tmp - 1) {
            std::memcpy(&T_sol_data[0], &b_T_sol_data[0], M_tmp * sizeof(float));
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
  isInitialized_p35p_solver_single = true;
}

//
// Arguments    : void
// Return Type  : void
//
void p35p_solver_terminate()
{
  // (no terminate code required)
  isInitialized_p35p_solver_single = false;
}

//
// File trailer for p35p_solver.cpp
//
// [EOF]
//
