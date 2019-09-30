/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhseqr.cpp
 *
 * Code generation for function 'xzhseqr'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzhseqr.h"
#include "xgeqp3.h"
#include "xzlarfg.h"
#include "sqrt.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
int eml_zlahqr(creal32_T h_data[], int h_size[2])
{
  int info;
  int n;
  int u1;
  int itmax;
  int ldh;
  int i17;
  int j;
  int ix0_tmp;
  int i;
  float SMLNUM;
  boolean_T exitg1;
  float tst;
  float htmp2;
  float br;
  int L;
  boolean_T goto140;
  creal32_T sc;
  int its;
  boolean_T exitg2;
  float ab;
  float bb;
  int k;
  boolean_T exitg3;
  float ba;
  float t_re;
  float t_im;
  int ix0;
  boolean_T goto70;
  creal32_T x2;
  int m;
  float u_re;
  float u_im;
  float s;
  int b_k;
  creal32_T v[2];
  float b_SMLNUM;
  n = h_size[0];
  u1 = h_size[0];
  if (10 > u1) {
    u1 = 10;
  }

  itmax = 30 * u1;
  ldh = h_size[0];
  info = 0;
  if (1 != h_size[0]) {
    i17 = h_size[0];
    for (j = 0; j <= i17 - 4; j++) {
      u1 = j + h_size[0] * j;
      ix0_tmp = u1 + 2;
      h_data[ix0_tmp].re = 0.0F;
      h_data[ix0_tmp].im = 0.0F;
      u1 += 3;
      h_data[u1].re = 0.0F;
      h_data[u1].im = 0.0F;
    }

    if (1 <= n - 2) {
      i17 = (n + h_size[0] * (n - 3)) - 1;
      h_data[i17].re = 0.0F;
      h_data[i17].im = 0.0F;
    }

    for (i = 2; i <= n; i++) {
      i17 = (i + h_size[0] * (i - 2)) - 1;
      if (h_data[i17].im != 0.0F) {
        tst = h_data[(i + h_size[0] * (i - 2)) - 1].re;
        htmp2 = h_data[(i + h_size[0] * (i - 2)) - 1].im;
        br = std::abs(h_data[(i + h_size[0] * (i - 2)) - 1].re) + std::abs
          (h_data[(i + h_size[0] * (i - 2)) - 1].im);
        if (htmp2 == 0.0F) {
          sc.re = tst / br;
          sc.im = 0.0F;
        } else if (tst == 0.0F) {
          sc.re = 0.0F;
          sc.im = htmp2 / br;
        } else {
          sc.re = tst / br;
          sc.im = htmp2 / br;
        }

        ab = sc.re;
        bb = -sc.im;
        tst = sc.re;
        htmp2 = sc.im;
        br = rt_hypotf_snf(tst, htmp2);
        if (bb == 0.0F) {
          sc.re = ab / br;
          sc.im = 0.0F;
        } else if (ab == 0.0F) {
          sc.re = 0.0F;
          sc.im = bb / br;
        } else {
          sc.re = ab / br;
          sc.im = bb / br;
        }

        h_data[i17].re = rt_hypotf_snf(h_data[i17].re, h_data[(i + h_size[0] *
          (i - 2)) - 1].im);
        h_data[i17].im = 0.0F;
        ix0_tmp = (i - 1) * ldh;
        ix0 = i + ix0_tmp;
        i17 = ix0 + ldh * (n - i);
        for (k = ix0; ldh < 0 ? k >= i17 : k <= i17; k += ldh) {
          htmp2 = h_data[k - 1].re;
          tst = h_data[k - 1].im;
          h_data[k - 1].re = sc.re * htmp2 - sc.im * tst;
          h_data[k - 1].im = sc.re * tst + sc.im * htmp2;
        }

        ix0 = ix0_tmp + 1;
        sc.im = -sc.im;
        u1 = i + 1;
        if (n < u1) {
          u1 = n;
        }

        i17 = ix0_tmp + u1;
        for (k = ix0; k <= i17; k++) {
          htmp2 = h_data[k - 1].re;
          tst = h_data[k - 1].im;
          h_data[k - 1].re = sc.re * htmp2 - sc.im * tst;
          h_data[k - 1].im = sc.re * tst + sc.im * htmp2;
        }
      }
    }

    SMLNUM = 1.17549435E-38F * (static_cast<float>(n) / 1.1920929E-7F);
    i = n - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 >= 1)) {
      L = -1;
      goto140 = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its <= itmax)) {
        k = i;
        exitg3 = false;
        while ((!exitg3) && (k + 1 > L + 2)) {
          i17 = k + h_size[0] * (k - 1);
          htmp2 = std::abs(h_data[i17].re);
          ba = htmp2 + std::abs(h_data[k + h_size[0] * (k - 1)].im);
          if (ba <= SMLNUM) {
            exitg3 = true;
          } else {
            u1 = k + h_size[0] * k;
            bb = std::abs(h_data[u1].re) + std::abs(h_data[k + h_size[0] * k].im);
            tst = (std::abs(h_data[i17 - 1].re) + std::abs(h_data[(k + h_size[0]
                     * (k - 1)) - 1].im)) + bb;
            if (tst == 0.0F) {
              if (k - 1 >= 1) {
                tst = std::abs(h_data[(k + h_size[0] * (k - 2)) - 1].re);
              }

              if (k + 2 <= n) {
                tst += std::abs(h_data[u1 + 1].re);
              }
            }

            if (htmp2 <= 1.1920929E-7F * tst) {
              htmp2 = std::abs(h_data[u1 - 1].re) + std::abs(h_data[(k + h_size
                [0] * k) - 1].im);
              if (ba > htmp2) {
                ab = ba;
                ba = htmp2;
              } else {
                ab = htmp2;
              }

              htmp2 = std::abs(h_data[(k + h_size[0] * (k - 1)) - 1].re -
                               h_data[k + h_size[0] * k].re) + std::abs(h_data
                [(k + h_size[0] * (k - 1)) - 1].im - h_data[k + h_size[0] * k].
                im);
              if (bb > htmp2) {
                tst = bb;
                bb = htmp2;
              } else {
                tst = htmp2;
              }

              s = tst + ab;
              htmp2 = 1.1920929E-7F * (bb * (tst / s));
              if ((SMLNUM > htmp2) || rtIsNaNF(htmp2)) {
                b_SMLNUM = SMLNUM;
              } else {
                b_SMLNUM = htmp2;
              }

              if (ba * (ab / s) <= b_SMLNUM) {
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
          h_data[k + h_size[0] * (k - 1)].re = 0.0F;
          h_data[k + h_size[0] * (k - 1)].im = 0.0F;
        }

        if (k + 1 >= i + 1) {
          goto140 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            t_re = 0.75F * std::abs(h_data[(k + h_size[0] * k) + 1].re) +
              h_data[k + h_size[0] * k].re;
            t_im = h_data[k + h_size[0] * k].im;
          } else if (its == 20) {
            t_re = 0.75F * std::abs(h_data[i + h_size[0] * (i - 1)].re) +
              h_data[i + h_size[0] * i].re;
            t_im = h_data[i + h_size[0] * i].im;
          } else {
            u1 = i + h_size[0] * i;
            t_re = h_data[u1].re;
            t_im = h_data[i + h_size[0] * i].im;
            x2 = h_data[u1 - 1];
            c_sqrt(&x2);
            sc = h_data[i + h_size[0] * (i - 1)];
            c_sqrt(&sc);
            u_re = x2.re * sc.re - x2.im * sc.im;
            u_im = x2.re * sc.im + x2.im * sc.re;
            s = std::abs(u_re) + std::abs(u_im);
            if (s != 0.0F) {
              t_re = 0.5F * (h_data[(i + h_size[0] * (i - 1)) - 1].re -
                             h_data[u1].re);
              t_im = 0.5F * (h_data[(i + h_size[0] * (i - 1)) - 1].im - h_data[i
                             + h_size[0] * i].im);
              ba = std::abs(t_re) + std::abs(t_im);
              if ((!(s > ba)) && (!rtIsNaNF(ba))) {
                s = ba;
              }

              if (t_im == 0.0F) {
                x2.re = t_re / s;
                x2.im = 0.0F;
              } else if (t_re == 0.0F) {
                x2.re = 0.0F;
                x2.im = t_im / s;
              } else {
                x2.re = t_re / s;
                x2.im = t_im / s;
              }

              tst = x2.re;
              htmp2 = x2.im;
              ab = x2.re;
              bb = x2.im;
              x2.re = tst * ab - htmp2 * bb;
              x2.im = tst * bb + htmp2 * ab;
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

              ab = sc.re;
              bb = sc.im;
              tst = sc.re;
              htmp2 = sc.im;
              x2.re += ab * tst - bb * htmp2;
              x2.im += ab * htmp2 + bb * tst;
              c_sqrt(&x2);
              sc.re = s * x2.re;
              sc.im = s * x2.im;
              if (ba > 0.0F) {
                if (t_im == 0.0F) {
                  x2.re = t_re / ba;
                  x2.im = 0.0F;
                } else if (t_re == 0.0F) {
                  x2.re = 0.0F;
                  x2.im = t_im / ba;
                } else {
                  x2.re = t_re / ba;
                  x2.im = t_im / ba;
                }

                if (x2.re * sc.re + x2.im * sc.im < 0.0F) {
                  sc.re = -sc.re;
                  sc.im = -sc.im;
                }
              }

              br = t_re + sc.re;
              ab = t_im + sc.im;
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
                    htmp2 = 0.5F;
                  } else {
                    htmp2 = -0.5F;
                  }

                  if (ab > 0.0F) {
                    tst = 0.5F;
                  } else {
                    tst = -0.5F;
                  }

                  ba = (u_re * htmp2 + u_im * tst) / bb;
                  tst = (u_im * htmp2 - u_re * tst) / bb;
                } else {
                  s = br / ab;
                  tst = ab + s * br;
                  ba = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h_data[i + h_size[0] * i].re - (u_re * ba - u_im * tst);
              t_im = h_data[i + h_size[0] * i].im - (u_re * tst + u_im * ba);
            }
          }

          goto70 = false;
          m = i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            u1 = m + h_size[0] * (m - 1);
            ix0_tmp = u1 - 1;
            sc.re = h_data[ix0_tmp].re - t_re;
            sc.im = h_data[(m + h_size[0] * (m - 1)) - 1].im - t_im;
            tst = h_data[u1].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            ab = sc.re;
            bb = sc.im;
            if (bb == 0.0F) {
              sc.re = ab / s;
              sc.im = 0.0F;
            } else if (ab == 0.0F) {
              sc.re = 0.0F;
              sc.im = bb / s;
            } else {
              sc.re = ab / s;
              sc.im = bb / s;
            }

            tst /= s;
            v[0] = sc;
            v[1].re = tst;
            v[1].im = 0.0F;
            if (std::abs(h_data[(m + h_size[0] * (m - 2)) - 1].re) * std::abs
                (tst) <= 1.1920929E-7F * ((std::abs(sc.re) + std::abs(sc.im)) *
                 ((std::abs(h_data[ix0_tmp].re) + std::abs(h_data[(m + h_size[0]
                     * (m - 1)) - 1].im)) + (std::abs(h_data[m + h_size[0] * m].
                    re) + std::abs(h_data[m + h_size[0] * m].im))))) {
              goto70 = true;
              exitg3 = true;
            } else {
              m--;
            }
          }

          if (!goto70) {
            sc.re = h_data[k + h_size[0] * k].re - t_re;
            sc.im = h_data[k + h_size[0] * k].im - t_im;
            tst = h_data[(k + h_size[0] * k) + 1].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            if (sc.im == 0.0F) {
              v[0].re = sc.re / s;
              v[0].im = 0.0F;
            } else if (sc.re == 0.0F) {
              v[0].re = 0.0F;
              v[0].im = sc.im / s;
            } else {
              v[0].re = sc.re / s;
              v[0].im = sc.im / s;
            }

            tst /= s;
            v[1].re = tst;
            v[1].im = 0.0F;
          }

          for (b_k = m; b_k <= i; b_k++) {
            if (b_k > m) {
              u1 = b_k + h_size[0] * (b_k - 2);
              v[0] = h_data[u1 - 1];
              v[1] = h_data[u1];
            }

            sc = xzlarfg(&v[0], &v[1]);
            if (b_k > m) {
              h_data[(b_k + h_size[0] * (b_k - 2)) - 1] = v[0];
              h_data[b_k + h_size[0] * (b_k - 2)].re = 0.0F;
              h_data[b_k + h_size[0] * (b_k - 2)].im = 0.0F;
            }

            t_re = v[1].re;
            t_im = v[1].im;
            ab = sc.re * v[1].re - sc.im * v[1].im;
            for (j = b_k; j <= n; j++) {
              ix0_tmp = b_k + h_size[0] * (j - 1);
              u1 = ix0_tmp - 1;
              x2.re = (sc.re * h_data[u1].re - -sc.im * h_data[(b_k + h_size[0] *
                        (j - 1)) - 1].im) + ab * h_data[ix0_tmp].re;
              x2.im = (sc.re * h_data[(b_k + h_size[0] * (j - 1)) - 1].im +
                       -sc.im * h_data[(b_k + h_size[0] * (j - 1)) - 1].re) + ab
                * h_data[b_k + h_size[0] * (j - 1)].im;
              h_data[u1].re = h_data[(b_k + h_size[0] * (j - 1)) - 1].re - x2.re;
              h_data[u1].im = h_data[(b_k + h_size[0] * (j - 1)) - 1].im - x2.im;
              h_data[ix0_tmp].re = h_data[b_k + h_size[0] * (j - 1)].re - (x2.re
                * t_re - x2.im * t_im);
              h_data[ix0_tmp].im = h_data[b_k + h_size[0] * (j - 1)].im - (x2.re
                * t_im + x2.im * t_re);
            }

            if (b_k + 2 < i + 1) {
              i17 = b_k + 1;
            } else {
              i17 = i;
            }

            for (j = 0; j <= i17; j++) {
              ix0_tmp = j + h_size[0] * (b_k - 1);
              u1 = j + h_size[0] * b_k;
              x2.re = (sc.re * h_data[ix0_tmp].re - sc.im * h_data[j + h_size[0]
                       * (b_k - 1)].im) + ab * h_data[u1].re;
              x2.im = (sc.re * h_data[j + h_size[0] * (b_k - 1)].im + sc.im *
                       h_data[j + h_size[0] * (b_k - 1)].re) + ab * h_data[j +
                h_size[0] * b_k].im;
              h_data[ix0_tmp].re = h_data[j + h_size[0] * (b_k - 1)].re - x2.re;
              h_data[ix0_tmp].im = h_data[j + h_size[0] * (b_k - 1)].im - x2.im;
              h_data[u1].re = h_data[j + h_size[0] * b_k].re - (x2.re * t_re -
                x2.im * -t_im);
              h_data[u1].im = h_data[j + h_size[0] * b_k].im - (x2.re * -t_im +
                x2.im * t_re);
            }

            if ((b_k == m) && (m > k + 1)) {
              br = rt_hypotf_snf(1.0F - sc.re, 0.0F - sc.im);
              if (0.0F - sc.im == 0.0F) {
                t_re = (1.0F - sc.re) / br;
                t_im = 0.0F;
              } else if (1.0F - sc.re == 0.0F) {
                t_re = 0.0F;
                t_im = (0.0F - sc.im) / br;
              } else {
                t_re = (1.0F - sc.re) / br;
                t_im = (0.0F - sc.im) / br;
              }

              htmp2 = h_data[m + h_size[0] * (m - 1)].re;
              tst = h_data[m + h_size[0] * (m - 1)].im;
              h_data[m + h_size[0] * (m - 1)].re = htmp2 * t_re - tst * -t_im;
              h_data[m + h_size[0] * (m - 1)].im = htmp2 * -t_im + tst * t_re;
              if (m + 2 <= i + 1) {
                u1 = (m + h_size[0] * m) + 1;
                htmp2 = h_data[u1].re;
                tst = h_data[(m + h_size[0] * m) + 1].im;
                h_data[u1].re = htmp2 * t_re - tst * t_im;
                h_data[u1].im = htmp2 * t_im + tst * t_re;
              }

              for (j = m; j <= i + 1; j++) {
                if (j != m + 1) {
                  if (n > j) {
                    ix0 = j + j * ldh;
                    i17 = ix0 + ldh * ((n - j) - 1);
                    for (u1 = ix0; ldh < 0 ? u1 >= i17 : u1 <= i17; u1 += ldh) {
                      htmp2 = h_data[u1 - 1].re;
                      tst = h_data[u1 - 1].im;
                      h_data[u1 - 1].re = t_re * htmp2 - t_im * tst;
                      h_data[u1 - 1].im = t_re * tst + t_im * htmp2;
                    }
                  }

                  ix0_tmp = (j - 1) * ldh;
                  ix0 = ix0_tmp + 1;
                  i17 = (ix0_tmp + j) - 1;
                  for (u1 = ix0; u1 <= i17; u1++) {
                    htmp2 = h_data[u1 - 1].re;
                    tst = h_data[u1 - 1].im;
                    h_data[u1 - 1].re = t_re * htmp2 - -t_im * tst;
                    h_data[u1 - 1].im = t_re * tst + -t_im * htmp2;
                  }
                }
              }
            }
          }

          t_re = h_data[i + h_size[0] * (i - 1)].re;
          t_im = h_data[i + h_size[0] * (i - 1)].im;
          if (h_data[i + h_size[0] * (i - 1)].im != 0.0F) {
            tst = rt_hypotf_snf(h_data[i + h_size[0] * (i - 1)].re, h_data[i +
                                h_size[0] * (i - 1)].im);
            h_data[i + h_size[0] * (i - 1)].re = tst;
            h_data[i + h_size[0] * (i - 1)].im = 0.0F;
            if (t_im == 0.0F) {
              t_re /= tst;
              t_im = 0.0F;
            } else if (t_re == 0.0F) {
              t_re = 0.0F;
              t_im /= tst;
            } else {
              t_re /= tst;
              t_im /= tst;
            }

            if (n > i + 1) {
              ix0 = (i + (i + 1) * ldh) + 1;
              i17 = ix0 + ldh * ((n - i) - 2);
              for (k = ix0; ldh < 0 ? k >= i17 : k <= i17; k += ldh) {
                htmp2 = h_data[k - 1].re;
                tst = h_data[k - 1].im;
                h_data[k - 1].re = t_re * htmp2 - -t_im * tst;
                h_data[k - 1].im = t_re * tst + -t_im * htmp2;
              }
            }

            ix0_tmp = i * ldh;
            ix0 = ix0_tmp + 1;
            i17 = ix0_tmp + i;
            for (k = ix0; k <= i17; k++) {
              htmp2 = h_data[k - 1].re;
              tst = h_data[k - 1].im;
              h_data[k - 1].re = t_re * htmp2 - t_im * tst;
              h_data[k - 1].im = t_re * tst + t_im * htmp2;
            }
          }

          its++;
        }
      }

      if (!goto140) {
        info = i + 1;
        exitg1 = true;
      } else {
        i = L;
      }
    }
  }

  return info;
}

/* End of code generation (xzhseqr.cpp) */
