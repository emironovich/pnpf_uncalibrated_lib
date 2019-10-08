//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: p35p_solver.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Oct-2019 05:30:59
//
#ifndef P35P_SOLVER_H
#define P35P_SOLVER_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "p35p_solver_types.h"

// Function Declarations
extern void p35p_solver(const float X[12], const float x[4], const float y[4],
  float e, float *solution_num, float f_sol_data[], int f_sol_size[2], float
  R_sol_data[], int R_sol_size[3], float T_sol_data[], int T_sol_size[2]);
extern void p35p_solver_initialize();
extern void p35p_solver_terminate();

#endif

//
// File trailer for p35p_solver.h
//
// [EOF]
//
