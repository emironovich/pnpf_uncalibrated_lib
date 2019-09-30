/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlarf.h
 *
 * Code generation for function 'xzlarf'
 *
 */

#ifndef XZLARF_H
#define XZLARF_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "solve_P4Pf_types.h"

/* Function Declarations */
extern void b_xzlarf(int m, int n, int iv0, const creal32_T tau, creal32_T
                     C_data[], int ic0, int ldc, creal32_T work_data[]);
extern void xzlarf(int m, int n, int iv0, const creal32_T tau, creal32_T C_data[],
                   int ic0, int ldc, creal32_T work_data[]);

#endif

/* End of code generation (xzlarf.h) */
