/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * solve_P4Pf.h
 *
 * Code generation for function 'solve_P4Pf'
 *
 */

#ifndef SOLVE_P4PF_H
#define SOLVE_P4PF_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "solve_P4Pf_types.h"

/* Function Declarations */
extern void solve_P4Pf(const float X[12], const float u[4], const float v[4],
  float e, float *solution_num, float fs_data[], int fs_size[2], float Rs_data[],
  int Rs_size[3], float Ts_data[], int Ts_size[2]);

#endif

/* End of code generation (solve_P4Pf.h) */
