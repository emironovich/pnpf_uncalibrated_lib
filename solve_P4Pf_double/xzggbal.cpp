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
#include "solve_P4Pf.h"
#include "xzggbal.h"

/* Function Definitions */
void xzggbal(creal_T A_data[], int A_size[2], int *ilo, int *ihi, int
             rscale_data[], int rscale_size[1])
{
  int ii;
  int i28;
  int exitg2;
  int i;
  int j;
  boolean_T found;
  boolean_T exitg3;
  int nzcount;
  int jj;
  boolean_T exitg4;
  int atmp_re_tmp;
  int exitg1;
  double atmp_re;
  double atmp_im;
  rscale_size[0] = A_size[0];
  ii = A_size[0];
  for (i28 = 0; i28 < ii; i28++) {
    rscale_data[i28] = 1;
  }

  *ilo = 1;
  *ihi = A_size[0];
  if (A_size[0] <= 1) {
    *ihi = 1;
  } else {
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
          if ((A_data[(ii + A_size[0] * jj) - 1].re != 0.0) || (A_data[(ii +
                A_size[0] * jj) - 1].im != 0.0) || (ii == jj + 1)) {
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
        ii = A_size[0];
        if (i != *ihi) {
          for (nzcount = 1; nzcount <= ii; nzcount++) {
            jj = A_size[0] * (nzcount - 1);
            atmp_re_tmp = (i + jj) - 1;
            atmp_re = A_data[atmp_re_tmp].re;
            atmp_im = A_data[(i + A_size[0] * (nzcount - 1)) - 1].im;
            i28 = (*ihi + jj) - 1;
            A_data[atmp_re_tmp] = A_data[i28];
            A_data[i28].re = atmp_re;
            A_data[i28].im = atmp_im;
          }
        }

        if (j + 1 != *ihi) {
          for (nzcount = 0; nzcount < *ihi; nzcount++) {
            jj = nzcount + A_size[0] * j;
            atmp_re = A_data[jj].re;
            atmp_im = A_data[jj].im;
            i28 = nzcount + A_size[0] * (*ihi - 1);
            A_data[jj] = A_data[i28];
            A_data[i28].re = atmp_re;
            A_data[i28].im = atmp_im;
          }
        }

        rscale_data[*ihi - 1] = j + 1;
        (*ihi)--;
        if (*ihi == 1) {
          rscale_data[0] = 1;
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
            if ((A_data[(ii + A_size[0] * (jj - 1)) - 1].re != 0.0) || (A_data
                 [(ii + A_size[0] * (jj - 1)) - 1].im != 0.0) || (ii == jj)) {
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
          ii = A_size[0];
          if (i != *ilo) {
            for (nzcount = *ilo; nzcount <= ii; nzcount++) {
              jj = A_size[0] * (nzcount - 1);
              atmp_re_tmp = (i + jj) - 1;
              atmp_re = A_data[atmp_re_tmp].re;
              atmp_im = A_data[(i + A_size[0] * (nzcount - 1)) - 1].im;
              i28 = (*ilo + jj) - 1;
              A_data[atmp_re_tmp] = A_data[i28];
              A_data[i28].re = atmp_re;
              A_data[i28].im = atmp_im;
            }
          }

          if (j != *ilo) {
            for (nzcount = 0; nzcount < *ihi; nzcount++) {
              jj = nzcount + A_size[0] * (j - 1);
              atmp_re = A_data[jj].re;
              atmp_im = A_data[jj].im;
              i28 = nzcount + A_size[0] * (*ilo - 1);
              A_data[jj] = A_data[i28];
              A_data[i28].re = atmp_re;
              A_data[i28].im = atmp_im;
            }
          }

          rscale_data[*ilo - 1] = j;
          (*ilo)++;
          if (*ilo == *ihi) {
            rscale_data[*ilo - 1] = *ilo;
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }
  }
}

/* End of code generation (xzggbal.cpp) */
