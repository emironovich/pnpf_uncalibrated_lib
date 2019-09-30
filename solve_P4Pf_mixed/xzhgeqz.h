/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.h
 *
 * Code generation for function 'xzhgeqz'
 *
 */

#ifndef XZHGEQZ_H
#define XZHGEQZ_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "solve_P4Pf_types.h"

/* Function Declarations */
extern void xzhgeqz(const creal32_T A_data[], const int A_size[2], int ilo, int
                    ihi, int *info, creal32_T alpha1_data[], int alpha1_size[1],
                    creal32_T beta1_data[], int beta1_size[1]);

#endif

/* End of code generation (xzhgeqz.h) */
