/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * roots.cpp
 *
 * Code generation for function 'roots'
 *
 */

/* Include files */
#include <cmath>
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "roots.h"
#include "schur.h"
#include "xzgeev.h"
#include "anyNonFinite.h"

/* Function Definitions */
void roots(const float c[13], creal32_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int nTrailingZeros;
  int companDim;
  boolean_T exitg1;
  int j;
  int a_size[2];
  boolean_T exitg2;
  creal32_T a_data[144];
  float ctmp[13];
  int i5;
  boolean_T p;
  creal32_T eiga_data[12];
  int eiga_size[1];
  creal32_T beta1_data[12];
  int beta1_size[1];
  creal32_T T_data[144];
  int T_size[2];
  int exitg3;
  float eiga_data_re;
  float brm;
  float bim;
  float s;
  memset(&r_data[0], 0, 12U * sizeof(creal32_T));
  k1 = 1;
  while ((k1 <= 13) && (!(c[k1 - 1] != 0.0F))) {
    k1++;
  }

  k2 = 13;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0F))) {
    k2--;
  }

  nTrailingZeros = 12 - k2;
  if (k1 < k2) {
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp[j] = c[k1 + j] / c[k1 - 1];
        if (rtIsInfF(std::abs(ctmp[j]))) {
          exitg2 = true;
        } else {
          j++;
        }
      }

      if (j + 1 > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }

    if (companDim < 1) {
      if (1 > 13 - k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = 13 - k2;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      memset(&a_data[0], 0, (unsigned int)(companDim * companDim * static_cast<
              int>(sizeof(creal32_T))));
      for (j = 0; j <= companDim - 2; j++) {
        i5 = companDim * j;
        a_data[i5].re = -ctmp[j];
        a_data[i5].im = 0.0F;
        i5 = (j + i5) + 1;
        a_data[i5].re = 1.0F;
        a_data[i5].im = 0.0F;
      }

      i5 = companDim * (companDim - 1);
      a_data[i5].re = -ctmp[companDim - 1];
      a_data[i5].im = 0.0F;
      for (j = 0; j <= nTrailingZeros; j++) {
        r_data[j].re = 0.0F;
        r_data[j].im = 0.0F;
      }

      if (anyNonFinite(a_data, a_size)) {
        if (companDim == 1) {
          eiga_data[0].re = rtNaNF;
          eiga_data[0].im = 0.0F;
        } else {
          for (i5 = 0; i5 < companDim; i5++) {
            eiga_data[i5].re = rtNaNF;
            eiga_data[i5].im = 0.0F;
          }
        }
      } else if (companDim == 1) {
        eiga_data[0] = a_data[0];
      } else {
        p = true;
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j <= companDim - 1)) {
          k1 = 0;
          do {
            exitg3 = 0;
            if (k1 <= j) {
              if ((!(a_data[k1 + companDim * j].re == a_data[j + companDim * k1]
                     .re)) || (!(a_data[k1 + companDim * j].im == -a_data[j +
                                 companDim * k1].im))) {
                p = false;
                exitg3 = 1;
              } else {
                k1++;
              }
            } else {
              j++;
              exitg3 = 2;
            }
          } while (exitg3 == 0);

          if (exitg3 == 1) {
            exitg1 = true;
          }
        }

        if (p) {
          schur(a_data, a_size, T_data, T_size);
          k1 = T_size[0];
          for (j = 0; j < k1; j++) {
            eiga_data[j] = T_data[j + T_size[0] * j];
          }
        } else {
          xzgeev(a_data, a_size, &k1, eiga_data, eiga_size, beta1_data,
                 beta1_size);
          k1 = eiga_size[0];
          for (i5 = 0; i5 < k1; i5++) {
            eiga_data_re = eiga_data[i5].re;
            if (beta1_data[i5].im == 0.0F) {
              if (eiga_data[i5].im == 0.0F) {
                eiga_data[i5].re /= beta1_data[i5].re;
                eiga_data[i5].im = 0.0F;
              } else if (eiga_data[i5].re == 0.0F) {
                eiga_data[i5].re = 0.0F;
                eiga_data[i5].im /= beta1_data[i5].re;
              } else {
                eiga_data[i5].re /= beta1_data[i5].re;
                eiga_data[i5].im /= beta1_data[i5].re;
              }
            } else if (beta1_data[i5].re == 0.0F) {
              if (eiga_data[i5].re == 0.0F) {
                eiga_data[i5].re = eiga_data[i5].im / beta1_data[i5].im;
                eiga_data[i5].im = 0.0F;
              } else if (eiga_data[i5].im == 0.0F) {
                eiga_data[i5].re = 0.0F;
                eiga_data[i5].im = -(eiga_data_re / beta1_data[i5].im);
              } else {
                eiga_data[i5].re = eiga_data[i5].im / beta1_data[i5].im;
                eiga_data[i5].im = -(eiga_data_re / beta1_data[i5].im);
              }
            } else {
              brm = std::abs(beta1_data[i5].re);
              bim = std::abs(beta1_data[i5].im);
              if (brm > bim) {
                s = beta1_data[i5].im / beta1_data[i5].re;
                brm = beta1_data[i5].re + s * beta1_data[i5].im;
                eiga_data[i5].re = (eiga_data[i5].re + s * eiga_data[i5].im) /
                  brm;
                eiga_data[i5].im = (eiga_data[i5].im - s * eiga_data_re) / brm;
              } else if (bim == brm) {
                if (beta1_data[i5].re > 0.0F) {
                  bim = 0.5F;
                } else {
                  bim = -0.5F;
                }

                if (beta1_data[i5].im > 0.0F) {
                  s = 0.5F;
                } else {
                  s = -0.5F;
                }

                eiga_data[i5].re = (eiga_data[i5].re * bim + eiga_data[i5].im *
                                    s) / brm;
                eiga_data[i5].im = (eiga_data[i5].im * bim - eiga_data_re * s) /
                  brm;
              } else {
                s = beta1_data[i5].re / beta1_data[i5].im;
                brm = beta1_data[i5].im + s * beta1_data[i5].re;
                bim = s * eiga_data[i5].im - eiga_data[i5].re;
                eiga_data[i5].re = (s * eiga_data[i5].re + eiga_data[i5].im) /
                  brm;
                eiga_data[i5].im = bim / brm;
              }
            }
          }
        }
      }

      for (j = 0; j < companDim; j++) {
        r_data[(j - k2) + 13] = eiga_data[j];
      }

      r_size[0] = (companDim - k2) + 13;
    }
  } else if (1 > 13 - k2) {
    r_size[0] = 0;
  } else {
    r_size[0] = 13 - k2;
  }
}

/* End of code generation (roots.cpp) */
