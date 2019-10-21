#pragma once
#include "pnpf.h"
#include <Eigen/Dense>
//#include <iostream>

using namespace Eigen;
template <class T>
void p35p_solver(const Matrix<T, 3, 4> &points_3d,
                 const Matrix<T, 2, 4> &points_2d, T e, int *n, T *fs,
                 Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs) {}

template <>
void p35p_solver<double>(const Matrix<double, 3, 4> &points_3d,
                         const Matrix<double, 2, 4> &points_2d, double e,
                         int *n, double *fs, Matrix<double, 3, 3> *Rs,
                         Matrix<double, 3, 1> *Cs) {
  double X[12];
  int ind = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
      X[ind] = points_3d(i, j);
      ind++;
    }
  }
  double x[4], y[4];
  for (int i = 0; i < 4; ++i) {
    x[i] = points_2d(0, i);
    y[i] = points_2d(1, i);
  }

  int f_size[2], r_size[3], t_size[2];
  double f_data[10], r_data[90], t_data[30]; // todo: check sizes

  p35p_double(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

  for (int i = 0; i < *n; ++i) {
    fs[i] = f_data[i];
    for (int j = 0; j < 3; ++j) {
      Cs[i](j) = t_data[3 * i + j];
    }
    ind = 9 * i;
    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
        Rs[i](j, k) = r_data[ind];
        ind++;
      }
    }
  }
}

template <>
void p35p_solver<float>(const Matrix<float, 3, 4> &points_3d,
                        const Matrix<float, 2, 4> &points_2d, float e, int *n,
                        float *fs, Matrix<float, 3, 3> *Rs,
                        Matrix<float, 3, 1> *Cs) {
  float X[12];
  int ind = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
      X[ind] = points_3d(i, j);
      ind++;
    }
  }
  float x[4], y[4];
  for (int i = 0; i < 4; ++i) {
    x[i] = points_2d(0, i);
    y[i] = points_2d(1, i);
  }

  int f_size[2], r_size[3], t_size[2];
  float f_data[10], r_data[90], t_data[30]; // todo: check sizes

  p35p_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

  for (int i = 0; i < *n; ++i) {
    fs[i] = f_data[i];
    for (int j = 0; j < 3; ++j) {
      Cs[i](j) = t_data[3 * i + j];
    }
    ind = 9 * i;
    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
        Rs[i](j, k) = r_data[ind];
        ind++;
      }
    }
  }
}

//  template<>
//  void p35p_solver<float>(const float X[12], const float x[4], const float
//  y[4], float e, int *n, float f_data[], int f_size[2], float r_data[], int
//  r_size[3], float t_data[], int t_size[2]) {
//	  p35p_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data,
// t_size);
//  }

template <class T>
void p4p_solver(const Matrix<T, 3, 4> &points_3d,
                const Matrix<T, 2, 4> &points_2d, T e, int *n, T *fs,
                Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs) {}

template <>
void p4p_solver<double>(const Matrix<double, 3, 4> &points_3d,
                        const Matrix<double, 2, 4> &points_2d, double e, int *n,
                        double *fs, Matrix<double, 3, 3> *Rs,
                        Matrix<double, 3, 1> *Cs) {
  double X[12];
  int ind = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
      X[ind] = points_3d(i, j);
      ind++;
    }
  }
  double x[4], y[4];
  for (int i = 0; i < 4; ++i) {
    x[i] = points_2d(0, i);
    y[i] = points_2d(1, i);
  }

  int f_size[2], r_size[3], t_size[2];
  double f_data[10], r_data[90], t_data[30]; // todo: check sizes

  p4pf_double(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

  for (int i = 0; i < *n; ++i) {
    fs[i] = f_data[i];
    for (int j = 0; j < 3; ++j) {
      Cs[i](j) = t_data[3 * i + j];
    }
    ind = 9 * i;
    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
        Rs[i](j, k) = r_data[ind];
        ind++;
      }
    }
  }
}

template <>
void p4p_solver<float>(const Matrix<float, 3, 4> &points_3d,
                       const Matrix<float, 2, 4> &points_2d, float e, int *n,
                       float *fs, Matrix<float, 3, 3> *Rs,
                       Matrix<float, 3, 1> *Cs) {
  float X[12];
  int ind = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
      X[ind] = points_3d(i, j);
      ind++;
    }
  }
  float x[4], y[4];
  for (int i = 0; i < 4; ++i) {
    x[i] = points_2d(0, i);
    y[i] = points_2d(1, i);
  }

  int f_size[2], r_size[3], t_size[2];
  float f_data[10], r_data[90], t_data[30]; // todo: check sizes

  p4pf_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

  for (int i = 0; i < *n; ++i) {
    fs[i] = f_data[i];
    for (int j = 0; j < 3; ++j) {
      Cs[i](j) = t_data[3 * i + j];
    }
    ind = 9 * i;
    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
        Rs[i](j, k) = r_data[ind];
        ind++;
      }
    }
  }
}

//
// template<class T>
//  void p4p_solver(const T X[12], const T x[4], const T y[4],
//  T e, int *n, T f_data[], int f_size[2], T r_data[], int r_size
//  [3], T t_data[], int t_size[2]){}
//
//  template<>
//  void p4p_solver<double>(const double X[12], const double x[4], const double
//  y[4], double e, int *n, double f_data[], int f_size[2], double r_data[], int
//  r_size [3], double t_data[], int t_size[2]) { 		p4pf_double(X, x,
//  y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);
//  }
//
//  template<>
//  void p4p_solver<float>(const float X[12], const float x[4], const float
//  y[4], float e, int *n, float f_data[], int f_size[2], float r_data[], int
//  r_size[3], float t_data[], int t_size[2]) {
//	  p4pf_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data,
// t_size);
//  }