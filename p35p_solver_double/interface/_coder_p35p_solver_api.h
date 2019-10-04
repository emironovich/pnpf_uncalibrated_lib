/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_p35p_solver_api.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 04-Oct-2019 01:44:03
 */

#ifndef _CODER_P35P_SOLVER_API_H
#define _CODER_P35P_SOLVER_API_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void p35p_solver(real_T X[12], real_T x[4], real_T y[4], real_T e, real_T
  *solution_num, real_T f_sol_data[], int32_T f_sol_size[2], real_T R_sol_data[],
  int32_T R_sol_size[3], real_T T_sol_data[], int32_T T_sol_size[2]);
extern void p35p_solver_api(const mxArray * const prhs[4], int32_T nlhs, const
  mxArray *plhs[4]);
extern void p35p_solver_atexit(void);
extern void p35p_solver_initialize(void);
extern void p35p_solver_terminate(void);
extern void p35p_solver_xil_shutdown(void);
extern void p35p_solver_xil_terminate(void);

#endif

/*
 * File trailer for _coder_p35p_solver_api.h
 *
 * [EOF]
 */
