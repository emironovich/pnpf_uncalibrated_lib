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
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzhgeqz.h"
#include "xzlartg.h"
#include "sqrt.h"

/* Function Definitions */
void xzhgeqz(creal_T A[100], int ilo, int ihi, creal_T Z[100], int *info,
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
  int jp1;
  double ascale;
  double imAij;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
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
  creal_T b_ascale;
  creal_T shift;
  double ascale_re;
  double ad22_re;
  double ad22_im;
  int shift_re_tmp;
  int ad22_re_tmp;
  double t1_re;
  double t1_im;
  double t1_im_tmp;
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
      jp1 = j + 1;
      if (ihi < j + 1) {
        jp1 = ihi;
      }

      for (i = ilo; i <= jp1; i++) {
        reAij = A[(i + 10 * (j - 1)) - 1].re;
        imAij = A[(i + 10 * (j - 1)) - 1].im;
        if (reAij != 0.0) {
          reAij = std::abs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = reAij;
            firstNonZero = false;
          } else if (scale < reAij) {
            anorm = scale / reAij;
            sumsq = 1.0 + sumsq * anorm * anorm;
            scale = reAij;
          } else {
            anorm = reAij / scale;
            sumsq += anorm * anorm;
          }
        }

        if (imAij != 0.0) {
          reAij = std::abs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = reAij;
            firstNonZero = false;
          } else if (scale < reAij) {
            anorm = scale / reAij;
            sumsq = 1.0 + sumsq * anorm * anorm;
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
  jp1 = ihi + 1;
  for (j = jp1; j < 11; j++) {
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
          jp1 = ilast + 10 * ilastm1;
          if (std::abs(A[jp1].re) + std::abs(A[ilast + 10 * ilastm1].im) <=
              b_atol) {
            A[jp1].re = 0.0;
            A[jp1].im = 0.0;
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
                jp1 = j + 10 * (j - 1);
                if (std::abs(A[jp1].re) + std::abs(A[j + 10 * (j - 1)].im) <=
                    b_atol) {
                  A[jp1].re = 0.0;
                  A[jp1].im = 0.0;
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

              for (jp1 = 0; jp1 < 100; jp1++) {
                Z[jp1].re = rtNaN;
                Z[jp1].im = 0.0;
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
                anorm = ascale * A[ilastm1 + 10 * ilastm1].re;
                reAij = ascale * A[ilastm1 + 10 * ilastm1].im;
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

                anorm = ascale * A[ilast + 10 * ilast].re;
                reAij = ascale * A[ilast + 10 * ilast].im;
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
                anorm = ascale * A[ilastm1 + 10 * ilast].re;
                reAij = ascale * A[ilastm1 + 10 * ilast].im;
                if (reAij == 0.0) {
                  ascale_re = anorm / 0.31622776601683794;
                  sumsq = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  sumsq = reAij / 0.31622776601683794;
                } else {
                  ascale_re = anorm / 0.31622776601683794;
                  sumsq = reAij / 0.31622776601683794;
                }

                anorm = ascale * A[ilast + 10 * ilastm1].re;
                reAij = ascale * A[ilast + 10 * ilastm1].im;
                if (reAij == 0.0) {
                  scale = anorm / 0.31622776601683794;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  scale = 0.0;
                  anorm = reAij / 0.31622776601683794;
                } else {
                  scale = anorm / 0.31622776601683794;
                  anorm = reAij / 0.31622776601683794;
                }

                reAij = shift.re * ad22_re - shift.im * ad22_im;
                imAij = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (ascale_re * scale
                  - sumsq * anorm)) - reAij;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (ascale_re * anorm + sumsq
                  * scale)) - imAij;
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
                anorm = ascale * A[ilast + 10 * ilastm1].re;
                reAij = ascale * A[ilast + 10 * ilastm1].im;
                if (reAij == 0.0) {
                  ascale_re = anorm / 0.31622776601683794;
                  sumsq = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  sumsq = reAij / 0.31622776601683794;
                } else {
                  ascale_re = anorm / 0.31622776601683794;
                  sumsq = reAij / 0.31622776601683794;
                }

                eshift_re += ascale_re;
                eshift_im += sumsq;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp.re = ascale * A[j + 10 * j].re - shift.re *
                  0.31622776601683794;
                ctemp.im = ascale * A[j + 10 * j].im - shift.im *
                  0.31622776601683794;
                anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
                reAij = ascale * (std::abs(A[jp1 + 10 * j].re) + std::abs(A[jp1
                  + 10 * j].im));
                imAij = anorm;
                if (reAij > anorm) {
                  imAij = reAij;
                }

                if ((imAij < 1.0) && (imAij != 0.0)) {
                  anorm /= imAij;
                  reAij /= imAij;
                }

                if ((std::abs(A[j + 10 * (j - 1)].re) + std::abs(A[j + 10 * (j -
                       1)].im)) * reAij <= anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp.re = ascale * A[(ifirst + 10 * (ifirst - 1)) - 1].re -
                  shift.re * 0.31622776601683794;
                ctemp.im = ascale * A[(ifirst + 10 * (ifirst - 1)) - 1].im -
                  shift.im * 0.31622776601683794;
                goto90 = true;
              }
            }

            if (goto90) {
              goto90 = false;
              b_ascale.re = ascale * A[istart + 10 * (istart - 1)].re;
              b_ascale.im = ascale * A[istart + 10 * (istart - 1)].im;
              b_xzlartg(ctemp, b_ascale, &imAij, &shift);
              j = istart;
              jp1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  i = j + 10 * jp1;
                  xzlartg(A[i - 1], A[i], &imAij, &shift, &A[(j + 10 * jp1) - 1]);
                  A[i].re = 0.0;
                  A[i].im = 0.0;
                }

                for (jp1 = j; jp1 < 11; jp1++) {
                  shift_re_tmp = j + 10 * (jp1 - 1);
                  ad22_re_tmp = shift_re_tmp - 1;
                  ad22_re = imAij * A[ad22_re_tmp].re + (shift.re *
                    A[shift_re_tmp].re - shift.im * A[j + 10 * (jp1 - 1)].im);
                  ad22_im = imAij * A[(j + 10 * (jp1 - 1)) - 1].im + (shift.re *
                    A[j + 10 * (jp1 - 1)].im + shift.im * A[j + 10 * (jp1 - 1)].
                    re);
                  anorm = A[(j + 10 * (jp1 - 1)) - 1].im;
                  reAij = A[(j + 10 * (jp1 - 1)) - 1].re;
                  A[shift_re_tmp].re = imAij * A[j + 10 * (jp1 - 1)].re -
                    (shift.re * A[(j + 10 * (jp1 - 1)) - 1].re + shift.im * A[(j
                      + 10 * (jp1 - 1)) - 1].im);
                  A[shift_re_tmp].im = imAij * A[j + 10 * (jp1 - 1)].im -
                    (shift.re * anorm - shift.im * reAij);
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
                  shift_re_tmp = (i + 10 * (j - 1)) - 1;
                  ad22_re_tmp = (i + 10 * j) - 1;
                  ad22_re = imAij * A[ad22_re_tmp].re + (shift.re *
                    A[shift_re_tmp].re - shift.im * A[(i + 10 * (j - 1)) - 1].im);
                  ad22_im = imAij * A[(i + 10 * j) - 1].im + (shift.re * A[(i +
                    10 * (j - 1)) - 1].im + shift.im * A[(i + 10 * (j - 1)) - 1]
                    .re);
                  anorm = A[(i + 10 * j) - 1].im;
                  reAij = A[(i + 10 * j) - 1].re;
                  A[shift_re_tmp].re = imAij * A[(i + 10 * (j - 1)) - 1].re -
                    (shift.re * A[(i + 10 * j) - 1].re + shift.im * A[(i + 10 *
                      j) - 1].im);
                  A[shift_re_tmp].im = imAij * A[(i + 10 * (j - 1)) - 1].im -
                    (shift.re * anorm - shift.im * reAij);
                  A[ad22_re_tmp].re = ad22_re;
                  A[ad22_re_tmp].im = ad22_im;
                }

                for (i = 0; i < 10; i++) {
                  shift_re_tmp = i + 10 * (j - 1);
                  ad22_re_tmp = i + 10 * j;
                  ad22_re = imAij * Z[ad22_re_tmp].re + (shift.re *
                    Z[shift_re_tmp].re - shift.im * Z[i + 10 * (j - 1)].im);
                  ad22_im = imAij * Z[i + 10 * j].im + (shift.re * Z[i + 10 * (j
                    - 1)].im + shift.im * Z[i + 10 * (j - 1)].re);
                  anorm = Z[i + 10 * j].im;
                  reAij = Z[i + 10 * j].re;
                  Z[shift_re_tmp].re = imAij * Z[i + 10 * (j - 1)].re -
                    (shift.re * Z[i + 10 * j].re + shift.im * Z[i + 10 * j].im);
                  Z[shift_re_tmp].im = imAij * Z[i + 10 * (j - 1)].im -
                    (shift.re * anorm - shift.im * reAij);
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

      for (jp1 = 0; jp1 < 100; jp1++) {
        Z[jp1].re = rtNaN;
        Z[jp1].im = 0.0;
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

/* End of code generation (xzhgeqz.cpp) */
