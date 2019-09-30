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
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[12]);
static const mxArray *b_emlrt_marshallOut(const real32_T u_data[], const int32_T
  u_size[2]);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *x, const
  char_T *identifier, real32_T y[4]);
static const mxArray *c_emlrt_marshallOut(const real32_T u_data[], const int32_T
  u_size[3]);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[4]);
static const mxArray *d_emlrt_marshallOut(const real32_T u_data[], const int32_T
  u_size[2]);
static real32_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *e, const
  char_T *identifier);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *X, const
  char_T *identifier, real32_T y[12]);
static const mxArray *emlrt_marshallOut(const real32_T u);
static real32_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[12]);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[4]);
static real32_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[12])
{
  g_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const real32_T u_data[], const int32_T
  u_size[2])
{
  const mxArray *y;
  const mxArray *m1;
  real32_T *pData;
  int32_T i0;
  int32_T i;
  y = NULL;
  m1 = emlrtCreateNumericArray(2, u_size, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m1);
  i0 = 0;
  for (i = 0; i < u_size[1]; i++) {
    pData[i0] = u_data[i];
    i0++;
  }

  emlrtAssign(&y, m1);
  return y;
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *x, const
  char_T *identifier, real32_T y[4])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(x), &thisId, y);
  emlrtDestroyArray(&x);
}

static const mxArray *c_emlrt_marshallOut(const real32_T u_data[], const int32_T
  u_size[3])
{
  const mxArray *y;
  const mxArray *m2;
  real32_T *pData;
  int32_T i1;
  int32_T i;
  int32_T b_i;
  int32_T i2;
  y = NULL;
  m2 = emlrtCreateNumericArray(3, u_size, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m2);
  i1 = 0;
  for (i = 0; i < u_size[2]; i++) {
    for (b_i = 0; b_i < 3; b_i++) {
      i2 = 3 * b_i + 9 * i;
      pData[i1] = u_data[i2];
      i1++;
      pData[i1] = u_data[i2 + 1];
      i1++;
      pData[i1] = u_data[i2 + 2];
      i1++;
    }
  }

  emlrtAssign(&y, m2);
  return y;
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real32_T y[4])
{
  h_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *d_emlrt_marshallOut(const real32_T u_data[], const int32_T
  u_size[2])
{
  const mxArray *y;
  const mxArray *m3;
  real32_T *pData;
  int32_T i3;
  int32_T i;
  y = NULL;
  m3 = emlrtCreateNumericArray(2, u_size, mxSINGLE_CLASS, mxREAL);
  pData = (real32_T *)emlrtMxGetData(m3);
  i3 = 0;
  for (i = 0; i < u_size[1]; i++) {
    pData[i3] = u_data[3 * i];
    i3++;
    pData[i3] = u_data[1 + 3 * i];
    i3++;
    pData[i3] = u_data[2 + 3 * i];
    i3++;
  }

  emlrtAssign(&y, m3);
  return y;
}

static real32_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *e, const
  char_T *identifier)
{
  real32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(e), &thisId);
  emlrtDestroyArray(&e);
  return y;
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *X, const
  char_T *identifier, real32_T y[12])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(X), &thisId, y);
  emlrtDestroyArray(&X);
}

static const mxArray *emlrt_marshallOut(const real32_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
  *(real32_T *)emlrtMxGetData(m0) = u;
  emlrtAssign(&y, m0);
  return y;
}

static real32_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real32_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[12])
{
  static const int32_T dims[2] = { 3, 4 };

  emlrtCheckBuiltInR2012b(sp, (const emlrtMsgIdentifier *)msgId, src,
    "single|double", false, 2U, *(int32_T (*)[2])&dims[0]);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real32_T ret[4])
{
  static const int32_T dims[2] = { 1, 4 };

  emlrtCheckBuiltInR2012b(sp, (const emlrtMsgIdentifier *)msgId, src,
    "single|double", false, 2U, *(int32_T (*)[2])&dims[0]);
  emlrtImportArrayR2015b(sp, src, ret, 4, false);
  emlrtDestroyArray(&src);
}

static real32_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  real32_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, (const emlrtMsgIdentifier *)msgId, src,
    "single|double", false, 0U, (int32_T *)&dims);
  emlrtImportArrayR2015b(sp, src, &ret, 4, false);
  emlrtDestroyArray(&src);
  return ret;
}

void p35p_solver_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray *
                     plhs[4])
{
  real32_T X[12];
  real32_T x[4];
  real32_T y[4];
  real32_T e;
  real32_T solution_num;
  real32_T f_sol_data[10];
  int32_T f_sol_size[2];
  real32_T R_sol_data[90];
  int32_T R_sol_size[3];
  real32_T T_sol_data[30];
  int32_T T_sol_size[2];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "X", X);
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "x", x);
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "y", y);
  e = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "e");

  /* Invoke the target function */
  p35p_solver(X, x, y, e, &solution_num, f_sol_data, f_sol_size, R_sol_data,
              R_sol_size, T_sol_data, T_sol_size);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(solution_num);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(f_sol_data, f_sol_size);
  }

  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(R_sol_data, R_sol_size);
  }

  if (nlhs > 3) {
    plhs[3] = d_emlrt_marshallOut(T_sol_data, T_sol_size);
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
