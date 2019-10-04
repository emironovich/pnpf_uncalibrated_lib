//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzggbal.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xzggbal.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : creal_T A[100]
//                int *ilo
//                int *ihi
//                int rscale[10]
// Return Type  : void
//
void xzggbal(creal_T A[100], int *ilo, int *ihi, int rscale[10])
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
// File trailer for xzggbal.cpp
//
// [EOF]
//
