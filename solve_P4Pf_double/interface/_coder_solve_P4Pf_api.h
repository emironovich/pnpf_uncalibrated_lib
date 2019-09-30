/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_solve_P4Pf_api.h
 *
 * Code generation for function '_coder_solve_P4Pf_api'
 *
 */

#ifndef _CODER_SOLVE_P4PF_API_H
#define _CODER_SOLVE_P4PF_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_solve_P4Pf_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void solve_P4Pf(real_T X[12], real_T u[4], real_T v[4], real_T e, real_T *
  solution_num, real_T fs_data[], int32_T fs_size[2], real_T Rs_data[], int32_T
  Rs_size[3], real_T Ts_data[], int32_T Ts_size[2]);
extern void solve_P4Pf_api(const mxArray * const prhs[4], int32_T nlhs, const
  mxArray *plhs[4]);
extern void solve_P4Pf_atexit(void);
extern void solve_P4Pf_initialize(void);
extern void solve_P4Pf_terminate(void);
extern void solve_P4Pf_xil_shutdown(void);
extern void solve_P4Pf_xil_terminate(void);

#endif

/* End of code generation (_coder_solve_P4Pf_api.h) */
