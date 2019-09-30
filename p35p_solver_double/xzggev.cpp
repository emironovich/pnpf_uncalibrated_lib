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
void xzggev(creal_T A[100], int *info, creal_T alpha1[10], creal_T beta1[10],
            creal_T V[100])
{
  double anrm;
  boolean_T ilascl;
  int i;
  double anrmto;
  int i26;
  int ihi;
  int rscale[10];
  double vtemp;
  double y;
  *info = 0;
  anrm = xzlangeM(A);
  if (!b_isfinite(anrm)) {
    for (i = 0; i < 10; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }

    for (i26 = 0; i26 < 100; i26++) {
      V[i26].re = rtNaN;
      V[i26].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
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

        if (vtemp >= 6.7178761075670888E-139) {
          vtemp = 1.0 / vtemp;
          for (ihi = 0; ihi < 10; ihi++) {
            i26 = ihi + 10 * i;
            V[i26].re *= vtemp;
            V[i26].im *= vtemp;
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
