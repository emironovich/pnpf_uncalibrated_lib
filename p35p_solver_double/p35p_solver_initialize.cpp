//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: p35p_solver_initialize.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "p35p_solver_initialize.h"
#include "p35p_solver.h"
#include "p35p_solver_data.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void p35p_solver_initialize()
{
  rt_InitInfAndNaN();
  isInitialized_p35p_solver = true;
}

//
// File trailer for p35p_solver_initialize.cpp
//
// [EOF]
//
