/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlarfg.cpp
 *
 * Code generation for function 'xzlarfg'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzlarfg.h"
#include "xnrm2.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
float xzlarfg(int n, float *alpha1, float x[3])
{
  float tau;
  float xnorm;
  float beta1;
  int knt;
  int k;
  tau = 0.0F;
  if (n > 0) {
    xnorm = b_xnrm2(n - 1, x);
    if (xnorm != 0.0F) {
      beta1 = rt_hypotf_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0F) {
        beta1 = -beta1;
      }

      if (std::abs(beta1) < 9.86076132E-32F) {
        knt = -1;
        do {
          knt++;
          for (k = 2; k <= n; k++) {
            x[k - 1] *= 1.01412048E+31F;
          }

          beta1 *= 1.01412048E+31F;
          *alpha1 *= 1.01412048E+31F;
        } while (!(std::abs(beta1) >= 9.86076132E-32F));

        beta1 = rt_hypotf_snf(*alpha1, b_xnrm2(n - 1, x));
        if (*alpha1 >= 0.0F) {
          beta1 = -beta1;
        }

        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0F / (*alpha1 - beta1);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= xnorm;
        }

        for (k = 0; k <= knt; k++) {
          beta1 *= 9.86076132E-32F;
        }

        *alpha1 = beta1;
      } else {
        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0F / (*alpha1 - beta1);
        for (k = 2; k <= n; k++) {
          x[k - 1] *= xnorm;
        }

        *alpha1 = beta1;
      }
    }
  }

  return tau;
}

/* End of code generation (xzlarfg.cpp) */
