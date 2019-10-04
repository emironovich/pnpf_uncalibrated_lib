//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: cat.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "cat.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : const double varargin_1_data[]
//                const int varargin_1_size[3]
//                const double varargin_2[9]
//                double y_data[]
//                int y_size[3]
// Return Type  : void
//
void cat(const double varargin_1_data[], const int varargin_1_size[3], const
         double varargin_2[9], double y_data[], int y_size[3])
{
  int iy;
  int i;
  int j;
  y_size[0] = 3;
  y_size[1] = 3;
  y_size[2] = static_cast<signed char>((varargin_1_size[2] + 1));
  iy = -1;
  i = 9 * varargin_1_size[2];
  for (j = 0; j < i; j++) {
    iy++;
    y_data[iy] = varargin_1_data[j];
  }

  for (j = 0; j < 9; j++) {
    iy++;
    y_data[iy] = varargin_2[j];
  }
}

//
// File trailer for cat.cpp
//
// [EOF]
//
