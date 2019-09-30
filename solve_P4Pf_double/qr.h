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
extern void b_qr(const double A[9], double Q[9], double R[9]);
extern void qr(const double A[32], double Q[64], double R[32]);

#endif

/* End of code generation (qr.h) */
