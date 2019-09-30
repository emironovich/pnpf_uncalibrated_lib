/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * make_mult_matrix.cpp
 *
 * Code generation for function 'make_mult_matrix'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "make_mult_matrix.h"

/* Function Definitions */
void make_mult_matrix(const double C[200], double M[100])
{
  int i;
  int i6;

  /* this function creates a matrix for multiplication by x in the monomial  */
  /* basis B */
  /* monomial basis B = {x^3, ...., 1} -- monomials up to the 3d degree, #B = 10 */
  /* x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1 */
  memset(&M[0], 0, 100U * sizeof(double));
  for (i = 0; i < 4; i++) {
    for (i6 = 0; i6 < 10; i6++) {
      M[i6 + 10 * i] = -C[(i + 20 * i6) + 15];
    }
  }

  M[40] = 1.0;
  M[51] = 1.0;
  M[62] = 1.0;
  M[74] = 1.0;
  M[85] = 1.0;
  M[97] = 1.0;
}

/* End of code generation (make_mult_matrix.cpp) */
