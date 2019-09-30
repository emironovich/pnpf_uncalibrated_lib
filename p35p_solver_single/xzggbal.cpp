/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzggbal.cpp
 *
 * Code generation for function 'xzggbal'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzggbal.h"

/* Function Definitions */
void xzggbal(creal32_T A[100], int *ilo, int *ihi, int rscale[10])
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
  float atmp_re;
  float atmp_im;
  int exitg1;
  int A_re_tmp;
  for (i = 0; i < 10; i++) {
    rscale[i] = 1;
  }

  *ilo = 1;
  *ihi = 10;
  do {
    exitg2 = 0;
    i = 0;
    j = -1;
    found = false;
    ii = *ihi;
    exitg3 = false;
    while ((!exitg3) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = *ihi - 1;
      jj = 0;
      exitg4 = false;
      while ((!exitg4) && (jj <= *ihi - 1)) {
        A_re_tmp = (ii + 10 * jj) - 1;
        if ((A[A_re_tmp].re != 0.0F) || (A[A_re_tmp].im != 0.0F) || (ii == jj +
             1)) {
          if (nzcount == 0) {
            j = jj;
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
        for (ii = 0; ii < 10; ii++) {
          nzcount = (i + 10 * ii) - 1;
          atmp_re = A[nzcount].re;
          atmp_im = A[(i + 10 * ii) - 1].im;
          jj = (*ihi + 10 * ii) - 1;
          A[nzcount] = A[jj];
          A[jj].re = atmp_re;
          A[jj].im = atmp_im;
        }
      }

      if (j + 1 != *ihi) {
        for (ii = 0; ii < *ihi; ii++) {
          nzcount = ii + 10 * j;
          atmp_re = A[nzcount].re;
          atmp_im = A[nzcount].im;
          jj = ii + 10 * (*ihi - 1);
          A[nzcount] = A[jj];
          A[jj].re = atmp_re;
          A[jj].im = atmp_im;
        }
      }

      rscale[*ihi - 1] = j + 1;
      (*ihi)--;
      if (*ihi == 1) {
        rscale[0] = 1;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
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
          A_re_tmp = (ii + 10 * (jj - 1)) - 1;
          if ((A[A_re_tmp].re != 0.0F) || (A[A_re_tmp].im != 0.0F) || (ii == jj))
          {
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
          for (ii = *ilo; ii < 11; ii++) {
            nzcount = 10 * (ii - 1);
            A_re_tmp = (i + nzcount) - 1;
            atmp_re = A[A_re_tmp].re;
            atmp_im = A[(i + 10 * (ii - 1)) - 1].im;
            jj = (*ilo + nzcount) - 1;
            A[A_re_tmp] = A[jj];
            A[jj].re = atmp_re;
            A[jj].im = atmp_im;
          }
        }

        if (j != *ilo) {
          for (ii = 0; ii < *ihi; ii++) {
            nzcount = ii + 10 * (j - 1);
            atmp_re = A[nzcount].re;
            atmp_im = A[nzcount].im;
            jj = ii + 10 * (*ilo - 1);
            A[nzcount] = A[jj];
            A[jj].re = atmp_re;
            A[jj].im = atmp_im;
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

/* End of code generation (xzggbal.cpp) */
