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
#include "eig.h"

/* Function Definitions */
void roots(const double c[13], creal_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int nTrailingZeros;
  int companDim;
  boolean_T exitg1;
  int j;
  int a_size[2];
  boolean_T exitg2;
  creal_T a_data[144];
  double ctmp[13];
  creal_T eiga_data[12];
  int eiga_size[1];
  memset(&r_data[0], 0, 12U * sizeof(creal_T));
  k1 = 1;
  while ((k1 <= 13) && (!(c[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = 13;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0))) {
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
        if (rtIsInf(std::abs(ctmp[j]))) {
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
              int>(sizeof(creal_T))));
      for (k1 = 0; k1 <= companDim - 2; k1++) {
        j = companDim * k1;
        a_data[j].re = -ctmp[k1];
        a_data[j].im = 0.0;
        j = (k1 + j) + 1;
        a_data[j].re = 1.0;
        a_data[j].im = 0.0;
      }

      j = companDim * (companDim - 1);
      a_data[j].re = -ctmp[companDim - 1];
      a_data[j].im = 0.0;
      for (k1 = 0; k1 <= nTrailingZeros; k1++) {
        r_data[k1].re = 0.0;
        r_data[k1].im = 0.0;
      }

      eig(a_data, a_size, eiga_data, eiga_size);
      for (k1 = 0; k1 < companDim; k1++) {
        r_data[(k1 - k2) + 13] = eiga_data[k1];
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
