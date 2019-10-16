#pragma once
#include "pnpf.h"
template<class T> 
  void p35p_solver(const T X[12], const T x[4], const T y[4],
  T e, int *n, T f_data[], int f_size[2], T r_data[], int r_size
  [3], T t_data[], int t_size[2]){}
	
  template<>	
  void p35p_solver<double>(const double X[12], const double x[4], const double y[4],
  double e, int *n, double f_data[], int f_size[2], double r_data[], int r_size
  [3], double t_data[], int t_size[2]) {
  		p35p_double(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);
  }

  template<>
  void p35p_solver<float>(const float X[12], const float x[4], const float y[4],
  float e, int *n, float f_data[], int f_size[2], float r_data[], int r_size[3],
  float t_data[], int t_size[2]) {
	  p35p_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);
  }



template<class T> 
  void p4p_solver(const T X[12], const T x[4], const T y[4],
  T e, int *n, T f_data[], int f_size[2], T r_data[], int r_size
  [3], T t_data[], int t_size[2]){}
	
  template<>	
  void p4p_solver<double>(const double X[12], const double x[4], const double y[4],
  double e, int *n, double f_data[], int f_size[2], double r_data[], int r_size
  [3], double t_data[], int t_size[2]) {
  		p4pf_double(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);
  }

  template<>
  void p4p_solver<float>(const float X[12], const float x[4], const float y[4],
  float e, int *n, float f_data[], int f_size[2], float r_data[], int r_size[3],
  float t_data[], int t_size[2]) {
	  p4pf_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);
  }