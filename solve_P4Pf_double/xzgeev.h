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
#include "solve_P4Pf_types.h"

/* Function Declarations */
extern void xzgeev(const creal_T A_data[], const int A_size[2], int *info,
                   creal_T alpha1_data[], int alpha1_size[1], creal_T
                   beta1_data[], int beta1_size[1]);

#endif

/* End of code generation (xzgeev.h) */
