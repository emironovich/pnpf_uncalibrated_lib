/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * init_F.cpp
 *
 * Code generation for function 'init_F'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "init_F.h"
#include "quadruple_constraint.h"

/* Function Definitions */
void init_F(const double x[4], const double y[4], const double X[12], const
            double R[54], double F[72])
{
  double dv5[18];
  int i2;
  quadruple_constraint(1.0, 2.0, 3.0, x, y, X, R, dv5);
  for (i2 = 0; i2 < 6; i2++) {
    F[12 * i2] = dv5[3 * i2];
    F[4 + 12 * i2] = dv5[1 + 3 * i2];
    F[8 + 12 * i2] = dv5[2 + 3 * i2];
  }

  quadruple_constraint(1.0, 3.0, 2.0, x, y, X, R, dv5);
  for (i2 = 0; i2 < 6; i2++) {
    F[1 + 12 * i2] = dv5[3 * i2];
    F[12 * i2 + 5] = dv5[1 + 3 * i2];
    F[12 * i2 + 9] = dv5[2 + 3 * i2];
  }

  quadruple_constraint(2.0, 4.0, 3.0, x, y, X, R, dv5);
  for (i2 = 0; i2 < 6; i2++) {
    F[2 + 12 * i2] = dv5[3 * i2];
    F[12 * i2 + 6] = dv5[1 + 3 * i2];
    F[12 * i2 + 10] = dv5[2 + 3 * i2];
  }

  quadruple_constraint(3.0, 4.0, 2.0, x, y, X, R, dv5);
  for (i2 = 0; i2 < 6; i2++) {
    F[3 + 12 * i2] = dv5[3 * i2];
    F[12 * i2 + 7] = dv5[1 + 3 * i2];
    F[12 * i2 + 11] = dv5[2 + 3 * i2];
  }
}

/* End of code generation (init_F.cpp) */
