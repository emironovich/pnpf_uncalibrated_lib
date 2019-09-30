/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzggev.cpp
 *
 * Code generation for function 'xzggev'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzggev.h"
#include "xzlascl.h"
#include "xzggbak.h"
#include "xztgevc.h"
#include "xzhgeqz.h"
#include "xzgghrd.h"
#include "xzggbal.h"
#include "isfinite.h"
#include "xzlangeM.h"

/* Function Definitions */
void xzggev(creal32_T A[100], int *info, creal32_T alpha1[10], creal32_T beta1
            [10], creal32_T V[100])
{
  float anrm;
  boolean_T ilascl;
  int i;
  float anrmto;
  int i20;
  int ihi;
  int rscale[10];
  float vtemp;
  float y;
  *info = 0;
  anrm = xzlangeM(A);
  if (!b_isfinite(anrm)) {
    for (i = 0; i < 10; i++) {
      alpha1[i].re = rtNaNF;
      alpha1[i].im = 0.0F;
      beta1[i].re = rtNaNF;
      beta1[i].im = 0.0F;
    }

    for (i20 = 0; i20 < 100; i20++) {
      V[i20].re = rtNaNF;
      V[i20].im = 0.0F;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0F) && (anrm < 9.09494702E-13F)) {
      anrmto = 9.09494702E-13F;
      ilascl = true;
    } else {
      if (anrm > 1.09951163E+12F) {
        anrmto = 1.09951163E+12F;
        ilascl = true;
      }
    }

    if (ilascl) {
      xzlascl(anrm, anrmto, A);
    }

    xzggbal(A, &i, &ihi, rscale);
    xzgghrd(i, ihi, A, V);
    xzhgeqz(A, i, ihi, V, info, alpha1, beta1);
    if (*info == 0) {
      xztgevc(A, V);
      xzggbak(V, i, ihi, rscale);
      for (i = 0; i < 10; i++) {
        vtemp = std::abs(V[10 * i].re) + std::abs(V[10 * i].im);
        for (ihi = 0; ihi < 9; ihi++) {
          y = std::abs(V[(ihi + 10 * i) + 1].re) + std::abs(V[(ihi + 10 * i) + 1]
            .im);
          if (y > vtemp) {
            vtemp = y;
          }
        }

        if (vtemp >= 9.09494702E-13F) {
          vtemp = 1.0F / vtemp;
          for (ihi = 0; ihi < 10; ihi++) {
            i20 = ihi + 10 * i;
            V[i20].re *= vtemp;
            V[i20].im *= vtemp;
          }
        }
      }

      if (ilascl) {
        b_xzlascl(anrmto, anrm, alpha1);
      }
    }
  }
}

/* End of code generation (xzggev.cpp) */
