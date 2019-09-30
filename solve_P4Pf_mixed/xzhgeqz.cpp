/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.cpp
 *
 * Code generation for function 'xzhgeqz'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzhgeqz.h"
#include "xzlartg.h"
#include "sqrt.h"

/* Function Definitions */
void xzhgeqz(const creal32_T A_data[], const int A_size[2], int ilo, int ihi,
             int *info, creal32_T alpha1_data[], int alpha1_size[1], creal32_T
             beta1_data[], int beta1_size[1])
{
  int A_size_idx_0;
  int jm1;
  creal32_T b_A_data[144];
  int n;
  int ctemp_tmp;
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
  float ascale;
  int i;
  float bscale;
  float imAij;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
  int ilast;
  int ilastm1;
  int ifrstm;
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
  creal32_T b_ascale;
  creal32_T shift;
  float ascale_im;
  float ad22_re;
  float ad22_im;
  float t1_re;
  float t1_im;
  float t1_im_tmp;
  A_size_idx_0 = A_size[0];
  jm1 = A_size[0] * A_size[1];
  if (0 <= jm1 - 1) {
    memcpy(&b_A_data[0], &A_data[0], (unsigned int)(jm1 * static_cast<int>
            (sizeof(creal32_T))));
  }

  *info = 0;
  if ((A_size[0] == 1) && (A_size[1] == 1)) {
    ihi = 1;
  }

  n = A_size[0];
  alpha1_size[0] = A_size[0];
  if (0 <= A_size[0] - 1) {
    memset(&alpha1_data[0], 0, (unsigned int)(A_size[0] * static_cast<int>
            (sizeof(creal32_T))));
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
    sumsq = 0.0F;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ctemp_tmp = j + 1;
      if (ihi < j + 1) {
        ctemp_tmp = ihi;
      }

      for (i = ilo; i <= ctemp_tmp; i++) {
        reAij = A_data[(i + A_size[0] * (j - 1)) - 1].re;
        imAij = A_data[(i + A_size[0] * (j - 1)) - 1].im;
        if (reAij != 0.0F) {
          reAij = std::abs(reAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            scale = reAij;
            firstNonZero = false;
          } else if (scale < reAij) {
            anorm = scale / reAij;
            sumsq = 1.0F + sumsq * anorm * anorm;
            scale = reAij;
          } else {
            anorm = reAij / scale;
            sumsq += anorm * anorm;
          }
        }

        if (imAij != 0.0F) {
          reAij = std::abs(imAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            scale = reAij;
            firstNonZero = false;
          } else if (scale < reAij) {
            anorm = scale / reAij;
            sumsq = 1.0F + sumsq * anorm * anorm;
            scale = reAij;
          } else {
            anorm = reAij / scale;
            sumsq += anorm * anorm;
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
  bscale = 1.0F / std::sqrt((float)A_size[0]);
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
    ifrstm = ilo;
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
          if (std::abs(b_A_data[ctemp_tmp].re) + std::abs(b_A_data[ilast +
               A_size_idx_0 * ilastm1].im) <= b_atol) {
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
                if (std::abs(b_A_data[ctemp_tmp].re) + std::abs(b_A_data[j +
                     A_size_idx_0 * (j - 1)].im) <= b_atol) {
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
              for (ctemp_tmp = 0; ctemp_tmp < jm1; ctemp_tmp++) {
                alpha1_data[ctemp_tmp].re = rtNaNF;
                alpha1_data[ctemp_tmp].im = 0.0F;
              }

              jm1 = beta1_size[0];
              for (ctemp_tmp = 0; ctemp_tmp < jm1; ctemp_tmp++) {
                beta1_data[ctemp_tmp].re = rtNaNF;
                beta1_data[ctemp_tmp].im = 0.0F;
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
              if (ifrstm > ilast + 1) {
                ifrstm = ilo;
              }

              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              ifrstm = ifirst;
              if (iiter - iiter / 10 * 10 != 0) {
                anorm = ascale * b_A_data[ilastm1 + A_size_idx_0 * ilastm1].re;
                reAij = ascale * b_A_data[ilastm1 + A_size_idx_0 * ilastm1].im;
                if (reAij == 0.0F) {
                  shift.re = anorm / bscale;
                  shift.im = 0.0F;
                } else if (anorm == 0.0F) {
                  shift.re = 0.0F;
                  shift.im = reAij / bscale;
                } else {
                  shift.re = anorm / bscale;
                  shift.im = reAij / bscale;
                }

                anorm = ascale * b_A_data[ilast + A_size_idx_0 * ilast].re;
                reAij = ascale * b_A_data[ilast + A_size_idx_0 * ilast].im;
                if (reAij == 0.0F) {
                  ad22_re = anorm / bscale;
                  ad22_im = 0.0F;
                } else if (anorm == 0.0F) {
                  ad22_re = 0.0F;
                  ad22_im = reAij / bscale;
                } else {
                  ad22_re = anorm / bscale;
                  ad22_im = reAij / bscale;
                }

                t1_re = 0.5F * (shift.re + ad22_re);
                t1_im = 0.5F * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                jm1 = ilastm1 + A_size_idx_0 * ilast;
                anorm = ascale * b_A_data[jm1].re;
                reAij = ascale * b_A_data[jm1].im;
                if (reAij == 0.0F) {
                  sumsq = anorm / bscale;
                  ascale_im = 0.0F;
                } else if (anorm == 0.0F) {
                  sumsq = 0.0F;
                  ascale_im = reAij / bscale;
                } else {
                  sumsq = anorm / bscale;
                  ascale_im = reAij / bscale;
                }

                anorm = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].re;
                reAij = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].im;
                if (reAij == 0.0F) {
                  scale = anorm / bscale;
                  anorm = 0.0F;
                } else if (anorm == 0.0F) {
                  scale = 0.0F;
                  anorm = reAij / bscale;
                } else {
                  scale = anorm / bscale;
                  anorm = reAij / bscale;
                }

                reAij = shift.re * ad22_re - shift.im * ad22_im;
                imAij = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (sumsq * scale -
                  ascale_im * anorm)) - reAij;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (sumsq * anorm + ascale_im
                  * scale)) - imAij;
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
                anorm = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].re;
                reAij = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].im;
                if (reAij == 0.0F) {
                  sumsq = anorm / bscale;
                  ascale_im = 0.0F;
                } else if (anorm == 0.0F) {
                  sumsq = 0.0F;
                  ascale_im = reAij / bscale;
                } else {
                  sumsq = anorm / bscale;
                  ascale_im = reAij / bscale;
                }

                eshift_re += sumsq;
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
                reAij = ascale * (std::abs(b_A_data[n + jm1].re) + std::abs
                                  (b_A_data[n + A_size_idx_0 * j].im));
                imAij = anorm;
                if (reAij > anorm) {
                  imAij = reAij;
                }

                if ((imAij < 1.0F) && (imAij != 0.0F)) {
                  anorm /= imAij;
                  reAij /= imAij;
                }

                if ((std::abs(b_A_data[j + A_size_idx_0 * (j - 1)].re) + std::
                     abs(b_A_data[j + A_size_idx_0 * (j - 1)].im)) * reAij <=
                    anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  n = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp.re = ascale * b_A_data[(ifirst + A_size_idx_0 * (ifirst -
                  1)) - 1].re - shift.re * bscale;
                ctemp.im = ascale * b_A_data[(ifirst + A_size_idx_0 * (ifirst -
                  1)) - 1].im - shift.im * bscale;
                goto90 = true;
              }
            }

            if (goto90) {
              goto90 = false;
              jm1 = istart + A_size_idx_0 * (istart - 1);
              b_ascale.re = ascale * b_A_data[jm1].re;
              b_ascale.im = ascale * b_A_data[jm1].im;
              b_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              jm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  n = j + A_size_idx_0 * jm1;
                  xzlartg(b_A_data[n - 1], b_A_data[n], &anorm, &shift,
                          &b_A_data[(j + A_size_idx_0 * jm1) - 1]);
                  b_A_data[n].re = 0.0F;
                  b_A_data[n].im = 0.0F;
                }

                for (n = j; n <= ilastm; n++) {
                  jm1 = j + A_size_idx_0 * (n - 1);
                  ctemp_tmp = jm1 - 1;
                  ad22_re = anorm * b_A_data[ctemp_tmp].re + (shift.re *
                    b_A_data[jm1].re - shift.im * b_A_data[j + A_size_idx_0 * (n
                    - 1)].im);
                  ad22_im = anorm * b_A_data[(j + A_size_idx_0 * (n - 1)) - 1].
                    im + (shift.re * b_A_data[j + A_size_idx_0 * (n - 1)].im +
                          shift.im * b_A_data[j + A_size_idx_0 * (n - 1)].re);
                  reAij = b_A_data[(j + A_size_idx_0 * (n - 1)) - 1].re;
                  b_A_data[jm1].re = anorm * b_A_data[jm1].re - (shift.re *
                    b_A_data[ctemp_tmp].re + shift.im * b_A_data[(j +
                    A_size_idx_0 * (n - 1)) - 1].im);
                  b_A_data[jm1].im = anorm * b_A_data[jm1].im - (shift.re *
                    b_A_data[(j + A_size_idx_0 * (n - 1)) - 1].im - shift.im *
                    reAij);
                  b_A_data[ctemp_tmp].re = ad22_re;
                  b_A_data[ctemp_tmp].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                n = j;
                if (ilast + 1 < j + 2) {
                  n = ilast - 1;
                }

                for (i = ifrstm; i <= n + 2; i++) {
                  jm1 = (i + A_size_idx_0 * (j - 1)) - 1;
                  ctemp_tmp = (i + A_size_idx_0 * j) - 1;
                  ad22_re = anorm * b_A_data[ctemp_tmp].re + (shift.re *
                    b_A_data[jm1].re - shift.im * b_A_data[(i + A_size_idx_0 *
                    (j - 1)) - 1].im);
                  ad22_im = anorm * b_A_data[(i + A_size_idx_0 * j) - 1].im +
                    (shift.re * b_A_data[(i + A_size_idx_0 * (j - 1)) - 1].im +
                     shift.im * b_A_data[(i + A_size_idx_0 * (j - 1)) - 1].re);
                  reAij = b_A_data[ctemp_tmp].re;
                  b_A_data[jm1].re = anorm * b_A_data[jm1].re - (shift.re *
                    b_A_data[(i + A_size_idx_0 * j) - 1].re + shift.im *
                    b_A_data[(i + A_size_idx_0 * j) - 1].im);
                  b_A_data[jm1].im = anorm * b_A_data[jm1].im - (shift.re *
                    b_A_data[ctemp_tmp].im - shift.im * reAij);
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
      for (jm1 = 0; jm1 <= ilast; jm1++) {
        alpha1_data[jm1].re = rtNaNF;
        alpha1_data[jm1].im = 0.0F;
        beta1_data[jm1].re = rtNaNF;
        beta1_data[jm1].im = 0.0F;
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

/* End of code generation (xzhgeqz.cpp) */
