//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzggev.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "xzggev.h"
#include "p35p_solver.h"
#include "p35p_solver_rtwutil.h"
#include "rt_nonfinite.h"
#include "xzggbak.h"
#include "xzggbal.h"
#include "xzgghrd.h"
#include "xzhgeqz.h"
#include "xztgevc.h"
#include <cmath>

// Function Definitions

//
// Arguments    : creal_T A[100]
//                int *info
//                creal_T alpha1[10]
//                creal_T beta1[10]
//                creal_T V[100]
// Return Type  : void
//
void xzggev(creal_T A[100], int *info, creal_T alpha1[10], creal_T beta1[10],
            creal_T V[100])
{
  double anrm;
  int i;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  boolean_T guard1 = false;
  int ihi;
  int rscale[10];
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double a;
  int jr;
  *info = 0;
  anrm = 0.0;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i < 100)) {
    absxk = rt_hypotd_snf(A[i].re, A[i].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      i++;
    }
  }

  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    for (i = 0; i < 10; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }

    for (i = 0; i < 100; i++) {
      V[i].re = rtNaN;
      V[i].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    guard1 = false;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
      guard1 = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
        guard1 = true;
      }
    }

    if (guard1) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          a = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          a = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          a = ctoc / absxk;
          notdone = false;
        }

        for (i = 0; i < 100; i++) {
          A[i].re *= a;
          A[i].im *= a;
        }
      }
    }

    xzggbal(A, &i, &ihi, rscale);
    xzgghrd(i, ihi, A, V);
    xzhgeqz(A, i, ihi, V, info, alpha1, beta1);
    if (*info == 0) {
      xztgevc(A, V);
      xzggbak(V, i, ihi, rscale);
      for (ihi = 0; ihi < 10; ihi++) {
        absxk = std::abs(V[10 * ihi].re) + std::abs(V[10 * ihi].im);
        for (jr = 0; jr < 9; jr++) {
          i = (jr + 10 * ihi) + 1;
          ctoc = std::abs(V[i].re) + std::abs(V[i].im);
          if (ctoc > absxk) {
            absxk = ctoc;
          }
        }

        if (absxk >= 6.7178761075670888E-139) {
          absxk = 1.0 / absxk;
          for (jr = 0; jr < 10; jr++) {
            i = jr + 10 * ihi;
            V[i].re *= absxk;
            V[i].im *= absxk;
          }
        }
      }

      if (ilascl) {
        notdone = true;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            a = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            a = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            a = anrm / anrmto;
            notdone = false;
          }

          for (i = 0; i < 10; i++) {
            alpha1[i].re *= a;
            alpha1[i].im *= a;
          }
        }
      }
    }
  }
}

//
// File trailer for xzggev.cpp
//
// [EOF]
//
