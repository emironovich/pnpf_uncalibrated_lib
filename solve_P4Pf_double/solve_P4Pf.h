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
extern void solve_P4Pf(const double X[12], const double u[4], const double v[4],
  double e, double *solution_num, double fs_data[], int fs_size[2], double
  Rs_data[], int Rs_size[3], double Ts_data[], int Ts_size[2]);

#endif

/* End of code generation (solve_P4Pf.h) */
