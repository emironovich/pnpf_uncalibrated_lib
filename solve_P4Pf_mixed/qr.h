/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * qr.h
 *
 * Code generation for function 'qr'
 *
 */

#ifndef QR_H
#define QR_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "solve_P4Pf_types.h"

/* Function Declarations */
extern void b_qr(const float A[9], float Q[9], float R[9]);
extern void qr(const float A[32], float Q[64], float R[32]);

#endif

/* End of code generation (qr.h) */
