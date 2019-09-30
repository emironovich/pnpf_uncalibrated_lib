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
void b_xzlascl(double cfrom, double cto, creal_T A[10])
{
  double cfromc;
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double a;
  int i32;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((std::abs(cfrom1) > std::abs(ctoc)) && (ctoc != 0.0)) {
      a = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (std::abs(cto1) > std::abs(cfromc)) {
      a = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      a = ctoc / cfromc;
      notdone = false;
    }

    for (i32 = 0; i32 < 10; i32++) {
      A[i32].re *= a;
      A[i32].im *= a;
    }
  }
}

void xzlascl(double cfrom, double cto, creal_T A[100])
{
  double cfromc;
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double a;
  int i27;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((std::abs(cfrom1) > std::abs(ctoc)) && (ctoc != 0.0)) {
      a = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (std::abs(cto1) > std::abs(cfromc)) {
      a = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      a = ctoc / cfromc;
      notdone = false;
    }

    for (i27 = 0; i27 < 100; i27++) {
      A[i27].re *= a;
      A[i27].im *= a;
    }
  }
}

/* End of code generation (xzlascl.cpp) */
