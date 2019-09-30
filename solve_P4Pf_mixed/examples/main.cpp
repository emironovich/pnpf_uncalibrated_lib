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
#include "solve_P4Pf.h"
#include "main.h"
#include "solve_P4Pf_terminate.h"
#include "solve_P4Pf_initialize.h"

/* Function Declarations */
static void argInit_1x4_real32_T(float result[4]);
static void argInit_3x4_real32_T(float result[12]);
static float argInit_real32_T();
static void main_solve_P4Pf();

/* Function Definitions */
static void argInit_1x4_real32_T(float result[4])
{
  float result_tmp;

  /* Loop over the array to initialize each element. */
  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result_tmp = argInit_real32_T();
  result[0] = result_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[1] = result_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[2] = argInit_real32_T();

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[3] = argInit_real32_T();
}

static void argInit_3x4_real32_T(float result[12])
{
  int idx0;
  float result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_real32_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 3] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 6] = argInit_real32_T();

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 9] = argInit_real32_T();
  }
}

static float argInit_real32_T()
{
  return 0.0F;
}

static void main_solve_P4Pf()
{
  float X[12];
  float u_tmp[4];
  float solution_num;
  float fs_data[10];
  int fs_size[2];
  float Rs_data[90];
  int Rs_size[3];
  float Ts_data[30];
  int Ts_size[2];

  /* Initialize function 'solve_P4Pf' input arguments. */
  /* Initialize function input argument 'X'. */
  argInit_3x4_real32_T(X);

  /* Initialize function input argument 'u'. */
  argInit_1x4_real32_T(u_tmp);

  /* Initialize function input argument 'v'. */
  /* Call the entry-point 'solve_P4Pf'. */
  solve_P4Pf(X, u_tmp, u_tmp, argInit_real32_T(), &solution_num, fs_data,
             fs_size, Rs_data, Rs_size, Ts_data, Ts_size);
}

int main(int, const char * const [])
{
  /* Initialize the application.
     You do not need to do this more than one time. */
  solve_P4Pf_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_solve_P4Pf();

  /* Terminate the application.
     You do not need to do this more than one time. */
  solve_P4Pf_terminate();
  return 0;
}

/* End of code generation (main.cpp) */
