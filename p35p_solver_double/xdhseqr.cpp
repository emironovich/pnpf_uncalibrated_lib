//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xdhseqr.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xdhseqr.h"
#include "p35p_solver.h"
#include "p35p_solver_rtwutil.h"
#include "rt_nonfinite.h"
#include "xdlanv2.h"
#include "xnrm2.h"
#include "xrot.h"
#include <cmath>

// Function Definitions

//
// Arguments    : double h[100]
//                double z[100]
// Return Type  : int
//
int eml_dlahqr(double h[100], double z[100])
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
// File trailer for xdhseqr.cpp
//
// [EOF]
//
