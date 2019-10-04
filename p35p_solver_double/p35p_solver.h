//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: p35p_solver.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//
#ifndef P35P_SOLVER_H
#define P35P_SOLVER_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "p35p_solver_types.h"

// Function Declarations
extern void p35p_solver(const double X[12], const double x[4], const double y[4],
  double e, double *solution_num, double f_sol_data[], int f_sol_size[2], double
  R_sol_data[], int R_sol_size[3], double T_sol_data[], int T_sol_size[2]);

#endif

//
// File trailer for p35p_solver.h
//
// [EOF]
//
