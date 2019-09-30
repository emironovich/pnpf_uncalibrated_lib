/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgemv.h
 *
 * Code generation for function 'xgemv'
 *
 */

#ifndef XGEMV_H
#define XGEMV_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "p35p_solver_types.h"

/* Function Declarations */
extern void b_xgemv(int m, int n, const float A[12], int ia0, const float x[12],
                    int ix0, float y[4]);
extern void xgemv(int m, int n, const float A[100], int ia0, const float x[100],
                  int ix0, float y[10]);

#endif

/* End of code generation (xgemv.h) */
