/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_p35p_solver_api.c
 *
 * Code generation for function '_coder_p35p_solver_api'
 *
 */

/* Include files */
#include "tmwtypes.h"
#include "_coder_p35p_solver_api.h"
#include "_coder_p35p_solver_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131482U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "p35p_solver",                       /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[12];
static const mxArray *b_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2]);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *x, const
  char_T *identifier))[4];
static const mxArray *c_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[3]);
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4];
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *e, const
  char_T *identifier);
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *X, const
  char_T *identifier))[12];
static const mxArray *emlrt_marshallOut(const real_T u);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[12];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4];
static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[12]
{
  real_T (*y)[12];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *b_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2])
{
  const mxArray *y;
  const mxArray *m1;
  static const int32_T iv0[2] = { 0, 0 };

  y = NULL;
  m1 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m1, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)m1, *(int32_T (*)[2])&u_size[0], 2);
  emlrtAssign(&y, m1);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *x, const
  char_T *identifier))[4]
{
  real_T (*y)[4];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(x), &thisId);
  emlrtDestroyArray(&x);
  return y;
}
  static const mxArray *c_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[3])
{
  const mxArray *y;
  const mxArray *m2;
  static const int32_T iv1[3] = { 0, 0, 0 };

  y = NULL;
  m2 = emlrtCreateNumericArray(3, iv1, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m2, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)m2, *(int32_T (*)[3])&u_size[0], 3);
  emlrtAssign(&y, m2);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4]
{
  real_T (*y)[4];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *e, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(e), &thisId);
  emlrtDestroyArray(&e);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *X, const
  char_T *identifier))[12]
{
  real_T (*y)[12];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(X), &thisId);
  emlrtDestroyArray(&X);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[12]
{
  real_T (*ret)[12];
  static const int32_T dims[2] = { 3, 4 };

  emlrtCheckBuiltInR2012b(sp, (const emlrtMsgIdentifier *)msgId, src, "double",
    false, 2U, *(int32_T (*)[2])&dims[0]);
  ret = (real_T (*)[12])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4]
{
  real_T (*ret)[4];
  static const int32_T dims[2] = { 1, 4 };

  emlrtCheckBuiltInR2012b(sp, (const emlrtMsgIdentifier *)msgId, src, "double",
    false, 2U, *(int32_T (*)[2])&dims[0]);
  ret = (real_T (*)[4])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, (const emlrtMsgIdentifier *)msgId, src, "double",
    false, 0U, (int32_T *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void p35p_solver_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray *
                     plhs[4])
{
  real_T (*f_sol_data)[10];
  real_T (*R_sol_data)[90];
  real_T (*T_sol_data)[30];
  real_T (*X)[12];
  real_T (*x)[4];
  real_T (*y)[4];
  real_T e;
  real_T solution_num;
  int32_T f_sol_size[2];
  int32_T R_sol_size[3];
  int32_T T_sol_size[2];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  f_sol_data = (real_T (*)[10])mxMalloc(sizeof(real_T [10]));
  R_sol_data = (real_T (*)[90])mxMalloc(sizeof(real_T [90]));
  T_sol_data = (real_T (*)[30])mxMalloc(sizeof(real_T [30]));

  /* Marshall function inputs */
  X = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "X");
  x = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "x");
  y = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "y");
  e = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "e");

  /* Invoke the target function */
  p35p_solver(*X, *x, *y, e, &solution_num, *f_sol_data, f_sol_size, *R_sol_data,
              R_sol_size, *T_sol_data, T_sol_size);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(solution_num);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(*f_sol_data, f_sol_size);
  }

  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(*R_sol_data, R_sol_size);
  }

  if (nlhs > 3) {
    plhs[3] = b_emlrt_marshallOut(*T_sol_data, T_sol_size);
  }
}

void p35p_solver_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  p35p_solver_xil_terminate();
  p35p_solver_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void p35p_solver_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void p35p_solver_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_p35p_solver_api.c) */
