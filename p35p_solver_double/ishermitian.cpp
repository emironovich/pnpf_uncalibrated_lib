/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ishermitian.cpp
 *
 * Code generation for function 'ishermitian'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "ishermitian.h"

/* Function Definitions */
boolean_T ishermitian(const double A[100])
{
  boolean_T p;
  int j;
  boolean_T exitg2;
  int i;
  int exitg1;
  p = true;
  j = 0;
  exitg2 = false;
  while ((!exitg2) && (j < 10)) {
    i = 0;
    do {
      exitg1 = 0;
      if (i <= j) {
        if (!(A[i + 10 * j] == A[j + 10 * i])) {
          p = false;
          exitg1 = 1;
        } else {
          i++;
        }
      } else {
        j++;
        exitg1 = 2;
      }
    } while (exitg1 == 0);

    if (exitg1 == 1) {
      exitg2 = true;
    }
  }

  return p;
}

/* End of code generation (ishermitian.cpp) */
