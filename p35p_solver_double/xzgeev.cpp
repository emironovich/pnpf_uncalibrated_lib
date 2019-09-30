/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeev.cpp
 *
 * Code generation for function 'xzgeev'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xzgeev.h"
#include "xnrm2.h"
#include "xzggev.h"

/* Function Definitions */
void xzgeev(const double A[100], int *info, creal_T alpha1[10], creal_T beta1[10],
            creal_T V[100])
{
  int i9;
  creal_T At[100];
  int coltop;
  double colnorm;
  int j;
  double V_re;
  double V_im;
  for (i9 = 0; i9 < 100; i9++) {
    At[i9].re = A[i9];
    At[i9].im = 0.0;
  }

  xzggev(At, info, alpha1, beta1, V);
  for (coltop = 0; coltop <= 90; coltop += 10) {
    colnorm = c_xnrm2(V, coltop + 1);
    i9 = coltop + 10;
    for (j = coltop + 1; j <= i9; j++) {
      V_re = V[j - 1].re;
      V_im = V[j - 1].im;
      if (V_im == 0.0) {
        V[j - 1].re = V_re / colnorm;
        V[j - 1].im = 0.0;
      } else if (V_re == 0.0) {
        V[j - 1].re = 0.0;
        V[j - 1].im = V_im / colnorm;
      } else {
        V[j - 1].re = V_re / colnorm;
        V[j - 1].im = V_im / colnorm;
      }
    }
  }
}

/* End of code generation (xzgeev.cpp) */
