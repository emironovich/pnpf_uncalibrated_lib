/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlascl.cpp
 *
 * Code generation for function 'xzlascl'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzlascl.h"

/* Function Definitions */
void b_xzlascl(float cfrom, float cto, creal32_T A[10])
{
  float cfromc;
  float ctoc;
  boolean_T notdone;
  float cfrom1;
  float cto1;
  float a;
  int i26;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    cfrom1 = cfromc * 1.97215226E-31F;
    cto1 = ctoc / 5.0706024E+30F;
    if ((std::abs(cfrom1) > std::abs(ctoc)) && (ctoc != 0.0F)) {
      a = 1.97215226E-31F;
      cfromc = cfrom1;
    } else if (std::abs(cto1) > std::abs(cfromc)) {
      a = 5.0706024E+30F;
      ctoc = cto1;
    } else {
      a = ctoc / cfromc;
      notdone = false;
    }

    for (i26 = 0; i26 < 10; i26++) {
      A[i26].re *= a;
      A[i26].im *= a;
    }
  }
}

void xzlascl(float cfrom, float cto, creal32_T A[100])
{
  float cfromc;
  float ctoc;
  boolean_T notdone;
  float cfrom1;
  float cto1;
  float a;
  int i21;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    cfrom1 = cfromc * 1.97215226E-31F;
    cto1 = ctoc / 5.0706024E+30F;
    if ((std::abs(cfrom1) > std::abs(ctoc)) && (ctoc != 0.0F)) {
      a = 1.97215226E-31F;
      cfromc = cfrom1;
    } else if (std::abs(cto1) > std::abs(cfromc)) {
      a = 5.0706024E+30F;
      ctoc = cto1;
    } else {
      a = ctoc / cfromc;
      notdone = false;
    }

    for (i21 = 0; i21 < 100; i21++) {
      A[i21].re *= a;
      A[i21].im *= a;
    }
  }
}

/* End of code generation (xzlascl.cpp) */
