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
#include <string.h>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzgeev.h"
#include "xzlascl.h"
#include "xzhgeqz.h"
#include "xzgghrd.h"
#include "xzggbal.h"
#include "isfinite.h"
#include "xzlangeM.h"

/* Function Definitions */
void xzgeev(const creal_T A_data[], const int A_size[2], int *info, creal_T
            alpha1_data[], int alpha1_size[1], creal_T beta1_data[], int
            beta1_size[1])
{
  int At_size[2];
  int loop_ub;
  creal_T At_data[144];
  double anrm;
  boolean_T ilascl;
  double anrmto;
  int ihi;
  int rscale_data[12];
  int rscale_size[1];
  At_size[0] = A_size[0];
  At_size[1] = A_size[1];
  loop_ub = A_size[0] * A_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&At_data[0], &A_data[0], (unsigned int)(loop_ub * static_cast<int>
            (sizeof(creal_T))));
  }

  *info = 0;
  anrm = xzlangeM(A_data, A_size);
  if (!b_isfinite(anrm)) {
    alpha1_size[0] = A_size[0];
    loop_ub = A_size[0];
    for (ihi = 0; ihi < loop_ub; ihi++) {
      alpha1_data[ihi].re = rtNaN;
      alpha1_data[ihi].im = 0.0;
    }

    beta1_size[0] = A_size[0];
    loop_ub = A_size[0];
    for (ihi = 0; ihi < loop_ub; ihi++) {
      beta1_data[ihi].re = rtNaN;
      beta1_data[ihi].im = 0.0;
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
      xzlascl(anrm, anrmto, At_data, At_size);
    }

    xzggbal(At_data, At_size, &loop_ub, &ihi, rscale_data, rscale_size);
    xzgghrd(loop_ub, ihi, At_data, At_size);
    xzhgeqz(At_data, At_size, loop_ub, ihi, info, alpha1_data, alpha1_size,
            beta1_data, beta1_size);
    if ((*info == 0) && ilascl) {
      b_xzlascl(anrmto, anrm, alpha1_data, alpha1_size);
    }
  }
}

/* End of code generation (xzgeev.cpp) */
