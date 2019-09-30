/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadruple_constraint.h
 *
 * Code generation for function 'quadruple_constraint'
 *
 */

#ifndef QUADRUPLE_CONSTRAINT_H
#define QUADRUPLE_CONSTRAINT_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "p35p_solver_types.h"

/* Function Declarations */
extern void quadruple_constraint(double i, double j, double k, const double x[4],
  const double y[4], const double X[12], const double R[54], double F_row[18]);

#endif

/* End of code generation (quadruple_constraint.h) */
