/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xdhseqr.cpp
 *
 * Code generation for function 'xdhseqr'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xdhseqr.h"
#include "xdlanv2.h"
#include "xzlarfg.h"

/* Function Definitions */
int eml_dlahqr(float h[100], float z[100])
{
  int info;
  float v[3];
  int ix;
  int i19;
  int i;
  boolean_T exitg1;
  int L;
  boolean_T goto150;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  int nr;
  float f0;
  int b_ix;
  float s;
  float tst;
  float htmp1;
  float aa;
  float ab;
  float ba;
  float rt2r;
  float rt1r;
  int iy;
  int m;
  int b_k;
  float f1;
  info = 0;
  v[0] = 0.0F;
  v[1] = 0.0F;
  v[2] = 0.0F;
  for (ix = 0; ix < 7; ix++) {
    i19 = ix + 10 * ix;
    h[i19 + 2] = 0.0F;
    h[i19 + 3] = 0.0F;
  }

  h[79] = 0.0F;
  i = 9;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        i19 = k + 10 * (k - 1);
        f0 = std::abs(h[i19]);
        if (f0 <= 9.86076132E-31F) {
          exitg3 = true;
        } else {
          nr = k + 10 * k;
          b_ix = i19 - 1;
          tst = std::abs(h[b_ix]) + std::abs(h[nr]);
          if (tst == 0.0F) {
            if (k - 1 >= 1) {
              tst = std::abs(h[(k + 10 * (k - 2)) - 1]);
            }

            if (k + 2 <= 10) {
              tst += std::abs(h[(k + 10 * k) + 1]);
            }
          }

          if (f0 <= 1.1920929E-7F * tst) {
            htmp1 = std::abs(h[i19]);
            tst = std::abs(h[nr - 1]);
            if (htmp1 > tst) {
              ab = htmp1;
              ba = tst;
            } else {
              ab = tst;
              ba = htmp1;
            }

            htmp1 = std::abs(h[k + 10 * k]);
            tst = std::abs(h[b_ix] - h[k + 10 * k]);
            if (htmp1 > tst) {
              aa = htmp1;
              htmp1 = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            tst = 1.1920929E-7F * (htmp1 * (aa / s));
            if ((9.86076132E-31F > tst) || rtIsNaNF(tst)) {
              f1 = 9.86076132E-31F;
            } else {
              f1 = tst;
            }

            if (ba * (ab / s) <= f1) {
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

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          s = std::abs(h[(k + 10 * k) + 1]) + std::abs(h[(k + 10 * (k + 1)) + 2]);
          tst = 0.75F * s + h[k + 10 * k];
          aa = -0.4375F * s;
          htmp1 = s;
          ba = tst;
        } else if (its == 20) {
          s = std::abs(h[i + 10 * (i - 1)]) + std::abs(h[(i + 10 * (i - 2)) - 1]);
          tst = 0.75F * s + h[i + 10 * i];
          aa = -0.4375F * s;
          htmp1 = s;
          ba = tst;
        } else {
          tst = h[(i + 10 * (i - 1)) - 1];
          htmp1 = h[i + 10 * (i - 1)];
          aa = h[(i + 10 * i) - 1];
          ba = h[i + 10 * i];
        }

        s = ((std::abs(tst) + std::abs(aa)) + std::abs(htmp1)) + std::abs(ba);
        if (s == 0.0F) {
          rt1r = 0.0F;
          tst = 0.0F;
          rt2r = 0.0F;
          htmp1 = 0.0F;
        } else {
          tst /= s;
          htmp1 /= s;
          aa /= s;
          ba /= s;
          ab = (tst + ba) / 2.0F;
          tst = (tst - ab) * (ba - ab) - aa * htmp1;
          htmp1 = std::sqrt(std::abs(tst));
          if (tst >= 0.0F) {
            rt1r = ab * s;
            rt2r = rt1r;
            tst = htmp1 * s;
            htmp1 = -tst;
          } else {
            rt1r = ab + htmp1;
            rt2r = ab - htmp1;
            if (std::abs(rt1r - ba) <= std::abs(rt2r - ba)) {
              rt1r *= s;
              rt2r = rt1r;
            } else {
              rt2r *= s;
              rt1r = rt2r;
            }

            tst = 0.0F;
            htmp1 = 0.0F;
          }
        }

        m = i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          nr = m + 10 * (m - 1);
          aa = h[nr - 1];
          ab = aa - rt2r;
          s = (std::abs(ab) + std::abs(htmp1)) + std::abs(h[nr]);
          ba = h[m + 10 * (m - 1)] / s;
          nr = m + 10 * m;
          v[0] = (ba * h[nr - 1] + (aa - rt1r) * (ab / s)) - tst * (htmp1 / s);
          v[1] = ba * (((aa + h[nr]) - rt1r) - rt2r);
          v[2] = ba * h[nr + 1];
          s = (std::abs(v[0]) + std::abs(v[1])) + std::abs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if (m == k + 1) {
            exitg3 = true;
          } else {
            i19 = m + 10 * (m - 2);
            if (std::abs(h[i19 - 1]) * (std::abs(v[1]) + std::abs(v[2])) <=
                1.1920929E-7F * std::abs(v[0]) * ((std::abs(h[i19 - 2]) + std::
                  abs(aa)) + std::abs(h[m + 10 * m]))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (b_k = m; b_k <= i; b_k++) {
          nr = (i - b_k) + 2;
          if (3 < nr) {
            nr = 3;
          }

          if (b_k > m) {
            b_ix = (b_k + 10 * (b_k - 2)) - 1;
            for (ix = 0; ix < nr; ix++) {
              v[ix] = h[ix + b_ix];
            }
          }

          tst = v[0];
          ab = xzlarfg(nr, &tst, v);
          v[0] = tst;
          if (b_k > m) {
            h[(b_k + 10 * (b_k - 2)) - 1] = tst;
            h[b_k + 10 * (b_k - 2)] = 0.0F;
            if (b_k < i) {
              h[(b_k + 10 * (b_k - 2)) + 1] = 0.0F;
            }
          } else {
            if (m > k + 1) {
              h[(b_k + 10 * (b_k - 2)) - 1] *= 1.0F - ab;
            }
          }

          f0 = v[1];
          htmp1 = ab * v[1];
          if (nr == 3) {
            s = v[2];
            tst = ab * v[2];
            for (ix = b_k; ix < 11; ix++) {
              nr = b_k + 10 * (ix - 1);
              b_ix = nr - 1;
              iy = nr + 1;
              aa = (h[b_ix] + f0 * h[nr]) + s * h[iy];
              h[b_ix] -= aa * ab;
              h[nr] -= aa * htmp1;
              h[iy] -= aa * tst;
            }

            if (b_k + 3 < i + 1) {
              i19 = b_k + 2;
            } else {
              i19 = i;
            }

            for (ix = 0; ix <= i19; ix++) {
              nr = ix + 10 * (b_k - 1);
              b_ix = ix + 10 * b_k;
              iy = ix + 10 * (b_k + 1);
              aa = (h[nr] + f0 * h[b_ix]) + s * h[iy];
              h[nr] -= aa * ab;
              h[b_ix] -= aa * htmp1;
              h[iy] -= aa * tst;
            }

            for (ix = 0; ix < 10; ix++) {
              nr = ix + 10 * (b_k - 1);
              b_ix = ix + 10 * b_k;
              iy = ix + 10 * (b_k + 1);
              aa = (z[nr] + f0 * z[b_ix]) + s * z[iy];
              z[nr] -= aa * ab;
              z[b_ix] -= aa * htmp1;
              z[iy] -= aa * tst;
            }
          } else {
            if (nr == 2) {
              for (ix = b_k; ix < 11; ix++) {
                nr = b_k + 10 * (ix - 1);
                b_ix = nr - 1;
                tst = h[b_ix];
                aa = tst + f0 * h[nr];
                h[b_ix] = tst - aa * ab;
                h[nr] -= aa * htmp1;
              }

              for (ix = 0; ix <= i; ix++) {
                nr = ix + 10 * (b_k - 1);
                b_ix = ix + 10 * b_k;
                aa = h[nr] + f0 * h[b_ix];
                h[nr] -= aa * ab;
                h[b_ix] -= aa * htmp1;
              }

              for (ix = 0; ix < 10; ix++) {
                b_ix = ix + 10 * (b_k - 1);
                tst = z[b_ix];
                nr = ix + 10 * b_k;
                aa = tst + f0 * z[nr];
                z[b_ix] = tst - aa * ab;
                z[nr] -= aa * htmp1;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = i + 1;
      exitg1 = true;
    } else {
      if ((L != i + 1) && (L == i)) {
        i19 = i + 10 * i;
        nr = i19 - 1;
        f0 = h[nr];
        ix = 10 * (i - 1);
        b_ix = i + ix;
        s = h[b_ix];
        tst = h[i19];
        xdlanv2(&h[(i + 10 * (i - 1)) - 1], &f0, &s, &tst, &htmp1, &aa, &ab, &ba,
                &rt2r, &rt1r);
        h[nr] = f0;
        h[b_ix] = s;
        h[i19] = tst;
        if (10 > i + 1) {
          nr = 8 - i;
          iy = i + (i + 1) * 10;
          b_ix = iy - 1;
          for (k = 0; k <= nr; k++) {
            tst = rt2r * h[b_ix] + rt1r * h[iy];
            h[iy] = rt2r * h[iy] - rt1r * h[b_ix];
            h[b_ix] = tst;
            iy += 10;
            b_ix += 10;
          }
        }

        if (i - 1 >= 1) {
          b_ix = ix;
          iy = i * 10;
          for (k = 0; k <= i - 2; k++) {
            tst = rt2r * h[b_ix] + rt1r * h[iy];
            h[iy] = rt2r * h[iy] - rt1r * h[b_ix];
            h[b_ix] = tst;
            iy++;
            b_ix++;
          }
        }

        iy = i * 10;
        for (k = 0; k < 10; k++) {
          tst = rt2r * z[ix] + rt1r * z[iy];
          z[iy] = rt2r * z[iy] - rt1r * z[ix];
          z[ix] = tst;
          iy++;
          ix++;
        }
      }

      i = L - 2;
    }
  }

  return info;
}

/* End of code generation (xdhseqr.cpp) */
