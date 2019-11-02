//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: pnpf.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 02-Nov-2019 09:23:45
//
#ifndef PNPF_H
#define PNPF_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "pnpf_types.h"

// Function Declarations
extern void p35pf_double(const double X[12], const double x[4], const double y[4],
  double e, int *n, double f_data[], int f_size[2], double r_data[], int r_size
  [3], double t_data[], int t_size[2]);
extern void p35pf_single(const float X[12], const float x[4], const float y[4],
  float e, int *n, float f_data[], int f_size[2], float r_data[], int r_size[3],
  float t_data[], int t_size[2]);
extern void p4pf_double(const double X[12], const double x[4], const double y[4],
  double e, int *n, double f_data[], int f_size[2], double r_data[], int r_size
  [3], double t_data[], int t_size[2]);
extern void p4pf_single(const float X[12], const float x[4], const float y[4],
  float e, int *n, float f_data[], int f_size[2], float r_data[], int r_size[3],
  float t_data[], int t_size[2]);
extern void pnpf_initialize();
extern void pnpf_terminate();

#endif

//
// File trailer for pnpf.h
//
// [EOF]
//
