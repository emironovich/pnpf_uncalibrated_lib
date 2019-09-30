/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeev.h
 *
 * Code generation for function 'xzgeev'
 *
 */

#ifndef XZGEEV_H
#define XZGEEV_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "p35p_solver_types.h"

/* Function Declarations */
extern void xzgeev(const float A[100], int *info, creal32_T alpha1[10],
                   creal32_T beta1[10], creal32_T V[100]);

#endif

/* End of code generation (xzgeev.h) */
