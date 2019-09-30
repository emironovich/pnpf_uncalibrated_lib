/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xaxpy.h
 *
 * Code generation for function 'xaxpy'
 *
 */

#ifndef XAXPY_H
#define XAXPY_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "solve_P4Pf_types.h"

/* Function Declarations */
extern void b_xaxpy(int n, float a, const float x[8], int ix0, float y[2]);
extern void c_xaxpy(int n, float a, const float x[2], float y[8], int iy0);
extern void xaxpy(int n, float a, int ix0, float y[8], int iy0);

#endif

/* End of code generation (xaxpy.h) */
