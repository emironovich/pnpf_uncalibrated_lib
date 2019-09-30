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
#include "solve_P4Pf.h"
#include "xzlascl.h"

/* Function Definitions */
void b_xzlascl(double cfrom, double cto, creal_T A_data[], int A_size[1])
{
  double cfromc;
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double a;
  int loop_ub;
  int i29;
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

    loop_ub = A_size[0];
    for (i29 = 0; i29 < loop_ub; i29++) {
      A_data[i29].re *= a;
      A_data[i29].im *= a;
    }
  }
}

void xzlascl(double cfrom, double cto, creal_T A_data[], int A_size[2])
{
  double cfromc;
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double a;
  int loop_ub;
  int i25;
  int b_loop_ub;
  int i26;
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

    loop_ub = A_size[1];
    for (i25 = 0; i25 < loop_ub; i25++) {
      b_loop_ub = A_size[0];
      for (i26 = 0; i26 < b_loop_ub; i26++) {
        i27 = i26 + A_size[0] * i25;
        A_data[i27].re *= a;
        A_data[i27].im = a * A_data[i26 + A_size[0] * i25].im;
      }
    }
  }
}

/* End of code generation (xzlascl.cpp) */
