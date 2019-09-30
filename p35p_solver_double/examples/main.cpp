/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.cpp
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "main.h"
#include "p35p_solver_terminate.h"
#include "p35p_solver_initialize.h"

/* Function Declarations */
static void argInit_1x4_real_T(double result[4]);
static void argInit_3x4_real_T(double result[12]);
static double argInit_real_T();
static void main_p35p_solver();

/* Function Definitions */
static void argInit_1x4_real_T(double result[4])
{
  double result_tmp;

  /* Loop over the array to initialize each element. */
  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result[0] = result_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[1] = result_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[2] = argInit_real_T();

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[3] = argInit_real_T();
}

static void argInit_3x4_real_T(double result[12])
{
  int idx0;
  double result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_real_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 3] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 6] = argInit_real_T();

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 9] = argInit_real_T();
  }
}

static double argInit_real_T()
{
  return 0.0;
}

static void main_p35p_solver()
{
  double X[12];
  double x_tmp[4];
  double solution_num;
  double f_sol_data[10];
  int f_sol_size[2];
  double R_sol_data[90];
  int R_sol_size[3];
  double T_sol_data[30];
  int T_sol_size[2];

  /* Initialize function 'p35p_solver' input arguments. */
  /* Initialize function input argument 'X'. */
  argInit_3x4_real_T(X);

  /* Initialize function input argument 'x'. */
  argInit_1x4_real_T(x_tmp);

  /* Initialize function input argument 'y'. */
  /* Call the entry-point 'p35p_solver'. */
  p35p_solver(X, x_tmp, x_tmp, argInit_real_T(), &solution_num, f_sol_data,
              f_sol_size, R_sol_data, R_sol_size, T_sol_data, T_sol_size);
}

int main(int, const char * const [])
{
  /* Initialize the application.
     You do not need to do this more than one time. */
  p35p_solver_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_p35p_solver();

  /* Terminate the application.
     You do not need to do this more than one time. */
  p35p_solver_terminate();
  return 0;
}

/* End of code generation (main.cpp) */
